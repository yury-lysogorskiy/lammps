/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Yury Lysogorskiy (ICAMS)
   grace/runtime/kk, pace/runtime/kk, see pair_grace_runtime_kokkos.h. Per compute(): ev_init, per-atom arrays, atomKK
   views, neighbour list; the engine of the model kind runs the chunked stages (2L / 3L: with the message exchange
   through LAMMPS comm) and returns the global energy and virial, added to eng_vdwl / virial.
------------------------------------------------------------------------- */

#ifndef NO_GRACE_RUNTIME    // GRACE_RUNTIME=OFF: ML-PACE without the native GRACE styles
#include "pair_grace_runtime_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "comm_kokkos.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor_kokkos.h"
#include "update.h"

#include "grace/kk/support.h"

#include <algorithm>
#include <cstring>
#include <string>
#include <type_traits>
#include <vector>

using namespace LAMMPS_NS;

// f(engine) for the engine of the model; nothing before pair_coeff
template<class DeviceType> template<class F> void PairGRACERuntimeKokkos<DeviceType>::with_engine(F &&f)
{
  std::visit(
      [&](auto &e) {
        if constexpr (!std::is_same_v<std::decay_t<decltype(e)>, std::monostate>) f(e);
      },
      eng);
}

template<class DeviceType>
PairGRACERuntimeKokkos<DeviceType>::PairGRACERuntimeKokkos(LAMMPS *lmp, bool pace_) :
    Pair(lmp), pace(pace_), name(pace_ ? "pace/runtime/kk" : "grace/runtime/kk")
{
  respa_enable = 0;
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  kokkosable = 1;
  reverse_comm_device = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
  no_virial_fdotr_compute = 1;
  if (!pace) {    // as the CPU style
    nextra = 1;    // dE/dlambda of pair_modify decouple/4d (compute pair ... -> c_ID[1])
    pvector = dc.pvector;
  }
  uq.gamma_only = pace;
  exchange2.p = this;
  exchange3.p = this;
}

template<class DeviceType> PairGRACERuntimeKokkos<DeviceType>::~PairGRACERuntimeKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(k_eatom, eatom);
  memoryKK->destroy_kokkos(k_vatom, vatom);
  memoryKK->destroy_kokkos(k_cvatom, cvatom);
  if (allocated) {    // map is freed by Pair::~Pair
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

template<class DeviceType> void PairGRACERuntimeKokkos<DeviceType>::settings(int narg, char **arg)
{
  if (strcmp("metal", update->unit_style) != 0) error->all(FLERR, "pair_style {} requires 'metal' units", name);
  // the kernels are written for GPUs; on CPUs grace/runtime and pace/runtime (libgrace CPU drivers) are the supported
  // path, a host Kokkos backend (Serial, OpenMP) is not
  if (execution_space != Device)
    error->all(FLERR, "pair_style {} needs a GPU Kokkos backend; on CPUs use pair_style {}", name,
               name.substr(0, name.find("/kk")));
  GraceRuntimeSettings s;
  s.parse(lmp, name, GraceRuntimeSettings::KK, narg, arg);
  if (coeff_args.empty()) {
    set = s;
    return;
  }
  // pair_style again after pair_coeff: the model of pair_coeff with the new keywords, its kind detected and checked
  // again (the file may have changed); nothing changes before the checks pass
  if (comm->me == 0) utils::logmesg(lmp, "[{}] new pair_style keywords: rebuilding the model of pair_coeff\n", name);
  std::vector<char *> av;
  for (auto &a : coeff_args) av.push_back(a.data());
  const GraceRuntime::ModelFile m = GraceRuntime::model_file(lmp, name, pace, s, (int) av.size(), av.data());
  set = s;
  mf = m;
  build();
}

template<class DeviceType> void *PairGRACERuntimeKokkos<DeviceType>::extract(const char *str, int &dim)
{
  dim = 0;    // every value exposed here is a global scalar (fix pair triggers, flags, decouple/4d lambda)
  if (!pace)    // PACE: no decouple/4d
    if (void *p = dc.extract(str, dim)) return p;
  // UQ trigger flags only when the model carries usable UQ artifacts, so fix pair ... gamma errors cleanly otherwise
  return has_uq ? uq.flag(str) : nullptr;
}

template<class DeviceType> void *PairGRACERuntimeKokkos<DeviceType>::extract_peratom(const char *str, int &ncol)
{
  ncol = 0;    // all UQ fields are per-atom scalars
  return has_uq ? uq.peratom(str) : nullptr;
}

// the engine of kind E: model loaded, engine settings, LAMMPS type -> model element
template<class DeviceType>
template<class E>
void PairGRACERuntimeKokkos<DeviceType>::load(const std::string &path, const std::string &asi)
{
  using Model = typename E::model_type;
  auto model = std::make_shared<Model>();
  try {
    if constexpr (std::is_same_v<E, EP>)
      model->load(path, asi);
    else
      model->load(path);
  } catch (const std::exception &e) {
    error->all(FLERR, "pair_style {}: cannot load {}: {}", name, path, e.what());
  }
  type2mu = GraceRuntime::type2mu(lmp, name, map, elements, model->element_names);
  E &e = eng.template emplace<E>();
  e.log = grace::Logger{GraceRuntime::log, lmp};
  e.chunksize = set.chunksize;
  e.debug_no_energy_only_calc = set.debug_no_energy_only_calc;
  e.spline_n = set.spline;
  if constexpr (std::is_same_v<E, E1> || std::is_same_v<E, E2>) e.no_shape = set.noshape;
  if constexpr (std::is_same_v<E, E2>) e.store_l1 = set.store_l1;
  if constexpr (std::is_same_v<E, E3>) {
    for (int k = 0; k < 3; k++) e.spline_n_layer[k] = set.spline_layer[k];
    e.spline_mb = set.spline_mb;
    e.store_fwd = set.store_fwd;
    e.fc_fold = set.fc_fold;
    e.fuse_radial = set.fuse_radial;
  }
  try {
    e.set_model(model, type2mu);
  } catch (const grace::Error &err) {
    error->all(FLERR, "pair_style {}: {}", name, err.what());
  }
  has_uq = e.uq_available();
}

template<class DeviceType> void PairGRACERuntimeKokkos<DeviceType>::allocate()
{
  allocated = 1;
  const int n = atom->ntypes + 1;
  memory->create(setflag, n, n, "pair:setflag");
  memory->create(cutsq, n, n, "pair:cutsq");
  map = new int[n];
}

template<class DeviceType> void PairGRACERuntimeKokkos<DeviceType>::coeff(int narg, char **arg)
{
  if (!allocated) allocate();
  // pair_coeff * * model El ...; PACE: model.yaml [model.asi] El ..., or its npz
  const GraceRuntime::ModelFile m = GraceRuntime::model_file(lmp, name, pace, set, narg, arg);
  map_element2type(narg - 3 - m.nskip, arg + 3 + m.nskip);
  mf = m;
  coeff_args.assign(arg, arg + narg);
  build();
}

// the engine of the model of pair_coeff (mf) with the keywords of pair_style (set), for the LAMMPS types (map)
template<class DeviceType> void PairGRACERuntimeKokkos<DeviceType>::build()
{
  eng.template emplace<std::monostate>();    // the old engine first: its device memory
  switch (mf.kind) {
    case grace::ModelKind::Grace1L: load<E1>(mf.path, mf.asi); break;
    case grace::ModelKind::Grace2L: load<E2>(mf.path, mf.asi); break;
    case grace::ModelKind::Grace3L: load<E3>(mf.path, mf.asi); break;
    case grace::ModelKind::Pace: load<EP>(mf.path, mf.asi); break;
  }
  if (comm->me == 0)
    for (int i = 1; i <= atom->ntypes; i++)
      if (type2mu[i] >= 0)
        utils::logmesg(lmp, "[{}] LAMMPS type {} -> {} (model index {})\n", name, i, elements[map[i]], type2mu[i]);
}

template<class DeviceType> void PairGRACERuntimeKokkos<DeviceType>::modify_params(int narg, char **arg)
{
  GraceGlue::Decouple::modify_params(this, narg, arg);
}

template<class DeviceType> int PairGRACERuntimeKokkos<DeviceType>::modify_param(int narg, char **arg)
{
  return pace ? 0 : dc.modify_param(lmp, name.c_str(), narg, arg);    // PACE: no decouple/4d
}

template<class DeviceType> void PairGRACERuntimeKokkos<DeviceType>::init_style()
{
  if (std::holds_alternative<std::monostate>(eng))
    error->all(FLERR, "pair_style {}: pair_coeff not set", name);
  if (atom->tag_enable == 0) error->all(FLERR, "pair_style {} requires atom IDs", name);
  if (force->newton_pair == 0) error->all(FLERR, "pair_style {} requires newton pair on", name);
  const int neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->add_request(this, NeighConst::REQ_FULL);
  request->set_kokkos_host(std::is_same_v<DeviceType, LMPHostType> && !std::is_same_v<DeviceType, LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType, LMPDeviceType>);
  if (neighflag == FULL) error->all(FLERR, "pair_style {} needs a half neighbor list (-pk kokkos neigh half)", name);
  if (dc.eq) atomKK->sync(Host, MASK_MASK);    // dc.init reads the host masks
  dc.init(lmp, name.c_str());
  if (ranks_per_device == 0) {
    // chunksize auto: the ranks on this rank's GPU share its free memory. Counted once per pair style: the ranks of
    // this node (shared-memory communicator) whose GPU has the same grace::kk::device_key() (device UUID, else PCI
    // bus id), compared exactly. With one visible GPU per rank (CUDA_VISIBLE_DEVICES / ROCR_VISIBLE_DEVICES per
    // rank) every rank's device index is 0, whichever GPU it is. Host execution space, or no key: the ranks of the
    // node without one. Serial build (STUBS): one rank.
    ranks_per_device = GraceRuntime::ranks_sharing(
        world, std::is_same_v<DeviceType, LMPDeviceType> ? grace::kk::device_key() : std::string());
  }
  with_engine([&](auto &e) {
    e.ranks_per_device = ranks_per_device;
    try {
      e.setup();
    } catch (const grace::Error &err) {    // per GPU (e.g. a GEMM self-test): may fail on some ranks only
      error->one(FLERR, "pair_style {}: {}", name, err.what());
    }
  });
  // MPI communication sizes (known after setup): one message row per atom (2L), a row of either layer (3L: eq1_norm /
  // eq2_norm, parity-doubled) or its adjoint; 1L, PACE: no layer messages
  comm_forward = comm_reverse = 0;
  if (auto *e = std::get_if<E2>(&eng)) comm_forward = comm_reverse = e->comm_width();
  if (auto *e = std::get_if<E3>(&eng)) comm_forward = comm_reverse = std::max(e->comm_width(0), e->comm_width(1));
}

template<class DeviceType> double PairGRACERuntimeKokkos<DeviceType>::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");
  double rc = 0.0;
  with_engine([&](auto &e) { rc = e.cutoff(type2mu[i], type2mu[j]); });
  return rc;
}

template<class DeviceType> double PairGRACERuntimeKokkos<DeviceType>::memory_usage()
{
  double b = 0.0;
  with_engine([&](auto &e) { b = e.memory_usage(); });
  return b;
}

template<class DeviceType> void PairGRACERuntimeKokkos<DeviceType>::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag, 0);
  with_engine([&](auto &e) { run(e); });
}

// per-atom arrays: reallocated (zeroed) only when maxeatom etc. (ev_init) outgrow them, otherwise zeroed on the device;
// the engine accumulates into them
template<class DeviceType> void PairGRACERuntimeKokkos<DeviceType>::prepare_peratom()
{
  if (eflag_atom) {
    if ((int) k_eatom.extent(0) < maxeatom) {
      memoryKK->destroy_kokkos(k_eatom, eatom);
      memoryKK->create_kokkos(k_eatom, eatom, maxeatom, "pair:eatom");
    } else
      Kokkos::deep_copy(k_eatom.template view<DeviceType>(), 0.0);
  }
  if (vflag_atom) {
    if ((int) k_vatom.extent(0) < maxvatom) {
      memoryKK->destroy_kokkos(k_vatom, vatom);
      memoryKK->create_kokkos(k_vatom, vatom, maxvatom, "pair:vatom");
    } else
      Kokkos::deep_copy(k_vatom.template view<DeviceType>(), 0.0);
  }
  if (cvflag_atom) {
    if ((int) k_cvatom.extent(0) < maxcvatom) {
      memoryKK->destroy_kokkos(k_cvatom, cvatom);
      memoryKK->create_kokkos(k_cvatom, cvatom, maxcvatom, "pair:cvatom");
    } else
      Kokkos::deep_copy(k_cvatom.template view<DeviceType>(), 0.0);
  }
}

// the per-atom arrays of the engine to the host
template<class DeviceType> void PairGRACERuntimeKokkos<DeviceType>::finish_peratom()
{
  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }
  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }
  if (cvflag_atom) {
    k_cvatom.template modify<DeviceType>();
    k_cvatom.sync_host();
  }
}

template<class DeviceType> template<class E> void PairGRACERuntimeKokkos<DeviceType>::run(E &e)
{
  prepare_peratom();

  atomKK->sync(execution_space, X_MASK | F_MASK | TYPE_MASK);
  auto *k_list = static_cast<NeighListKokkos<DeviceType> *>(list);
  typename E::RunInput in;
  in.x = atomKK->k_x.template view<DeviceType>();
  in.f = atomKK->k_f.template view<DeviceType>();
  in.type = atomKK->k_type.template view<DeviceType>();
  in.ilist = k_list->d_ilist;
  in.numneigh = k_list->d_numneigh;
  in.neighbors = k_list->d_neighbors;
  in.inum = list->inum;
  in.nlocal = atom->nlocal;
  in.ntypes = atom->ntypes;
  in.eflag_global = eflag_global;
  in.eflag_atom = eflag_atom;
  in.vflag_global = vflag_global;
  in.vflag_atom = vflag_atom;
  in.cvflag_atom = cvflag_atom;
  in.energy_only = eflag_only;
  if (dc.eq) {    // group masks of local and ghost atoms (border comm), on the device
    atomKK->sync(execution_space, MASK_MASK);
    in.couple_eq = dc.eq;
    in.couple_lambda = dc.lambda;
    in.couple_bit1 = dc.bit1;
    in.couple_bit2 = dc.bit2;
    in.mask = atomKK->k_mask.template view<DeviceType>();
  }
  if (eflag_atom) in.eatom = k_eatom.template view<DeviceType>();
  if (vflag_atom) in.vatom = k_vatom.template view<DeviceType>();
  if (cvflag_atom) in.cvatom = k_cvatom.template view<DeviceType>();
  const bool any_uq = uq.any(has_uq);
  in.uq = any_uq;
  if (any_uq) uq.grow(atom->nlocal);
  if constexpr (std::is_same_v<E, E2>) {    // the layer messages of the ghost atoms
    in.nall = atom->nlocal + atom->nghost;
    in.exchange = &exchange2;
  } else if constexpr (std::is_same_v<E, E3>) {
    in.nall = atom->nlocal + atom->nghost;
    in.ex = &exchange3;
  }

  typename E::RunOutput out;
  try {
    out = e.run(in);
  } catch (const grace::Error &err) {
    error->one(FLERR, "pair_style {}: {}", name, err.what());
  }
  eng_vdwl += out.eng_vdwl;
  for (int k = 0; k < 6; k++) virial[k] += out.virial[k];
  dc.pvector[0] = out.dedl;

  // per-atom UQ results device -> host arrays (indexed by local atom)
  if (any_uq) {
    const int nlocal = atom->nlocal;
    auto copy_out = [&](double *host, const Kokkos::View<double *, DeviceType> &dview) {
      if (dview.extent(0) == 0) return;    // an output the engine does not have (PACE: gamma only)
      Kokkos::View<double *, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> hwrap(host, nlocal);
      Kokkos::deep_copy(hwrap, Kokkos::subview(dview, Kokkos::make_pair(0, nlocal)));
    };
    copy_out(uq.gamma, e.uq_gamma());
    copy_out(uq.sigma, e.uq_sigma());
    copy_out(uq.cluster, e.uq_cluster());
  }

  finish_peratom();
  atomKK->modified(execution_space, F_MASK);
}

/* ----------------------------------------------------------------------
   layer messages: 2L (one layer) and 3L (layer comm_layer, in the order forward 0, forward 1, reverse 1, reverse 0)
------------------------------------------------------------------------- */

template<class DeviceType>
template<class ToHost, class FromHost>
void PairGRACERuntimeKokkos<DeviceType>::pair_exchange(bool reverse, ToHost &&to_host, FromHost &&from_host)
{
  // CommKokkos packs on the device: in legacy pair-comm mode (GPU-aware MPI off, several ranks) through
  // *_comm_pair_device (only the packed rows of the swaps with another rank go through the host), otherwise through
  // the usual comm call. Any other comm (comm_style tiled: CommTiledKokkos is not a CommKokkos) calls the host
  // pack / unpack, which read and write the host copy of the rows.
  auto *ck = dynamic_cast<CommKokkos *>(comm);
  if (ck && (reverse ? lmp->kokkos->reverse_pair_comm_legacy : lmp->kokkos->forward_pair_comm_legacy)) {
    if (reverse) ck->reverse_comm_pair_device(this);
    else ck->forward_comm_pair_device(this);
    return;
  }
  if (!ck) to_host();
  if (reverse) comm->reverse_comm(this);
  else comm->forward_comm(this);
  if (!ck) from_host();
}

template<class DeviceType> void PairGRACERuntimeKokkos<DeviceType>::Exchange2::forward()
{
  auto &e = std::get<E2>(p->eng);
  p->pair_exchange(false, [&] { e.copy_rows_to_host(false); }, [&] { e.copy_rows_from_host(false); });
}

template<class DeviceType> void PairGRACERuntimeKokkos<DeviceType>::Exchange2::reverse()
{
  auto &e = std::get<E2>(p->eng);
  p->pair_exchange(true, [&] { e.copy_rows_to_host(true); }, [&] { e.copy_rows_from_host(true); });
}

template<class DeviceType> void PairGRACERuntimeKokkos<DeviceType>::Exchange3::forward(int layer)
{
  auto &e = std::get<E3>(p->eng);
  p->comm_layer = layer;
  p->pair_exchange(false, [&] { e.copy_rows_to_host(layer, false); }, [&] { e.copy_rows_from_host(layer, false); });
}

template<class DeviceType> void PairGRACERuntimeKokkos<DeviceType>::Exchange3::reverse(int layer)
{
  auto &e = std::get<E3>(p->eng);
  p->comm_layer = layer;
  p->pair_exchange(true, [&] { e.copy_rows_to_host(layer, true); }, [&] { e.copy_rows_from_host(layer, true); });
}

// the comm callbacks run only during an exchange of a 2L or 3L engine
template<class DeviceType>
int PairGRACERuntimeKokkos<DeviceType>::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                                          int * /*pbc*/)
{
  if (auto *e = std::get_if<E2>(&eng)) return e->pack_forward_host(n, list, buf);
  return std::get<E3>(eng).pack_forward_host(comm_layer, n, list, buf);
}

template<class DeviceType> void PairGRACERuntimeKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  if (auto *e = std::get_if<E2>(&eng)) e->unpack_forward_host(n, first, buf);
  else std::get<E3>(eng).unpack_forward_host(comm_layer, n, first, buf);
}

template<class DeviceType> int PairGRACERuntimeKokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  if (auto *e = std::get_if<E2>(&eng)) return e->pack_reverse_host(n, first, buf);
  return std::get<E3>(eng).pack_reverse_host(comm_layer, n, first, buf);
}

template<class DeviceType> void PairGRACERuntimeKokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
{
  if (auto *e = std::get_if<E2>(&eng)) e->unpack_reverse_host(n, list, buf);
  else std::get<E3>(eng).unpack_reverse_host(comm_layer, n, list, buf);
}

template<class DeviceType>
int PairGRACERuntimeKokkos<DeviceType>::pack_forward_comm_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                                                                 DAT::tdual_double_1d &buf, int /*pbc_flag*/,
                                                                 int * /*pbc*/)
{
  const auto list = k_sendlist.view<DeviceType>();
  const auto b = buf.view<DeviceType>();
  if (auto *e = std::get_if<E2>(&eng)) return e->pack_forward(n, list, b);
  return std::get<E3>(eng).pack_forward(comm_layer, n, list, b);
}

template<class DeviceType>
void PairGRACERuntimeKokkos<DeviceType>::unpack_forward_comm_kokkos(int n, int first, DAT::tdual_double_1d &buf)
{
  const auto b = buf.view<DeviceType>();
  if (auto *e = std::get_if<E2>(&eng)) e->unpack_forward(n, first, b);
  else std::get<E3>(eng).unpack_forward(comm_layer, n, first, b);
}

template<class DeviceType>
int PairGRACERuntimeKokkos<DeviceType>::pack_reverse_comm_kokkos(int n, int first, DAT::tdual_double_1d &buf)
{
  const auto b = buf.view<DeviceType>();
  if (auto *e = std::get_if<E2>(&eng)) return e->pack_reverse(n, first, b);
  return std::get<E3>(eng).pack_reverse(comm_layer, n, first, b);
}

template<class DeviceType>
void PairGRACERuntimeKokkos<DeviceType>::unpack_reverse_comm_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                                                                    DAT::tdual_double_1d &buf)
{
  const auto list = k_sendlist.view<DeviceType>();
  const auto b = buf.view<DeviceType>();
  if (auto *e = std::get_if<E2>(&eng)) e->unpack_reverse(n, list, b);
  else std::get<E3>(eng).unpack_reverse(comm_layer, n, list, b);
}

namespace LAMMPS_NS {
template class PairGRACERuntimeKokkos<LMPDeviceType>;
template class PairPACERuntimeKokkos<LMPDeviceType>;
}    // namespace LAMMPS_NS
#endif    // NO_GRACE_RUNTIME
