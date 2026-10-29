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

   grace/runtime, pace/runtime (CPU): glue over the libgrace CPU drivers (grace-runtime), one per model kind:
     GRACE-1L  grace::Grace1LFused   type-homogeneous batches per OpenMP thread, no layer messages
     GRACE-2L  grace::Grace2LFused   one layer message I: forward comm to the ghosts, dE/dI reverse comm to the owners
     GRACE-3L  grace::Grace3LFused   two layer messages (I1, I2): forward 0, forward 1, reverse 1, reverse 0
     PACE      grace::PaceCpu        no layer messages; extrapolation grade gamma with an active set
   The kind comes from the model file at pair_coeff (grace::model_kind); the drivers hold the kernels (fp32 for GRACE
   with fp64 bond geometry and accumulation, fp64 for PACE), this file moves LAMMPS data in and out.
------------------------------------------------------------------------- */

#ifndef NO_GRACE_RUNTIME    // GRACE_RUNTIME=OFF: ML-PACE without the native GRACE styles
#include "pair_grace_runtime.h"

#include "grace/cpu_simd.h"     // isa_hint
#include "grace/fused_1l.h"
#include "grace/fused_2l.h"
#include "grace/fused_3l.h"
#include "grace/fused_cpu.h"    // profile counters (G2LF_PROF)
#include "grace/pace_cpu.h"
#include "grace/pace_model.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "exceptions.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"
#include "utils.h"

#include <algorithm>
#include <cstring>
#include <type_traits>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace LAMMPS_NS;

namespace {
using D1 = std::unique_ptr<grace::Grace1LFused>;
using D2 = std::unique_ptr<grace::Grace2LFused>;
using D3 = std::unique_ptr<grace::Grace3LFused>;
using DP = std::unique_ptr<grace::PaceCpu>;

// the layer messages over LAMMPS comm: forward = I_k to the ghosts, reverse = dE/dI_k of the ghosts summed into the
// owners (2L: k = 0 only)
struct LammpsExchange final : grace::Exchange {
  PairGRACERuntime *pair;
  Comm *comm;
  void forward(int k) override
  {
    pair->comm_layer = k;
    comm->forward_comm(pair);
  }
  void reverse(int k) override
  {
    pair->comm_layer = k;
    comm->reverse_comm(pair);
  }
};

// the 2L message: nI rows of PR floats per atom, CR of them used; I and dE/dI rows travel as raw fp32 pairs inside the
// double buffer, ceil(nI CR / 2) doubles per atom
int msg2_doubles(const grace::Grace2LFused &d)
{
  return (d.message_rows() * d.message_channels() + 1) / 2;
}
}    // namespace

// f(driver) for the driver of the model; nothing before pair_coeff
template<class F> void PairGRACERuntime::with_driver(F &&f)
{
  std::visit(
      [&](auto &d) {
        if constexpr (!std::is_same_v<std::decay_t<decltype(d)>, std::monostate>) f(d);
      },
      drv);
}

/* ---------------------------------------------------------------------- */

PairGRACERuntime::PairGRACERuntime(LAMMPS *lmp, bool pace_) :
    Pair(lmp), pace(pace_), name(pace_ ? "pace/runtime" : "grace/runtime")
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  uq.gamma_only = pace;
  if (!pace) {
    nextra = 1;    // dE/dlambda of pair_modify decouple/4d (compute pair ... -> c_ID[1])
    pvector = dc.pvector;
  }
}

PairGRACERuntime::~PairGRACERuntime()
{
  using namespace grace::fused;
  // the 1L / 2L profile (G2LF_PROF) goes to the log here, the 3L driver logs its own when it is destroyed
  const bool l1 = std::holds_alternative<D1>(drv), l2 = std::holds_alternative<D2>(drv);
  if (g_prof_on && comm->me == 0 && (l1 || l2)) {
    const double tot = l1 ? g_prof[T_P3]
                          : g_prof[T_P1] + g_prof[T_P2] + g_prof[T_P3] + g_prof[T_P4] + g_prof[T_P5] + g_prof[T_P6];
    utils::logmesg(lmp, "[{}] profile (thread 0 categories; {}), s:\n", name, l1 ? "total wall" : "phases wall");
    for (int q = 0; q < T_N; q++)
      if (l2 || g_prof[q] > 0.0)
        utils::logmesg(lmp, "  {:16s} {:10.4f}  {:5.1f} %{}\n", T_NAME[q], g_prof[q],
                       100.0 * g_prof[q] / std::max(tot, 1e-30),
                       g_flops[q] > 0 ? fmt::format("  {:6.1f} GF/s (mv)", g_flops[q] / g_prof[q] * 1e-9) : "");
  }
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

void PairGRACERuntime::allocate()
{
  allocated = 1;
  const int n = atom->ntypes + 1;
  memory->create(setflag, n, n, "pair:setflag");
  memory->create(cutsq, n, n, "pair:cutsq");
  map = new int[n];
}

/* ---------------------------------------------------------------------- */

void PairGRACERuntime::settings(int narg, char **arg)
{
  if (strcmp("metal", update->unit_style) != 0) error->all(FLERR, "pair_style {} requires 'metal' units", name);
  GraceRuntimeSettings s;
  s.parse(lmp, name, GraceRuntimeSettings::CPU, narg, arg);
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

/* ---------------------------------------------------------------------- */

void PairGRACERuntime::coeff(int narg, char **arg)
{
  if (!allocated) allocate();
  // pair_coeff * * model El ...; PACE: model.yaml [model.asi] El ... (as pair_style pace / pace/extrapolation), or
  // its npz
  const GraceRuntime::ModelFile m = GraceRuntime::model_file(lmp, name, pace, set, narg, arg);
  map_element2type(narg - 3 - m.nskip, arg + 3 + m.nskip);
  mf = m;
  coeff_args.assign(arg, arg + narg);
  build();
}

// the driver of the model of pair_coeff (mf) with the keywords of pair_style (set), for the LAMMPS types (map)
void PairGRACERuntime::build()
{
  drv = std::monostate();    // the old driver first: its memory, and the memory budgets of the new one (auto)
  const std::string &path = mf.path, &asi = mf.asi;

  // model load and driver setup throw grace::Error / std::exception; the element mapping between them reports through
  // error->all itself, so it stays outside the try blocks
  const int ntypes = atom->ntypes;
  auto load = [&](auto model) {
    try {
      if constexpr (std::is_same_v<std::decay_t<decltype(*model)>, grace::PaceModel>)
        model->load(path, asi);
      else
        model->load(path);
    } catch (const std::exception &e) {
      error->all(FLERR, "pair_style {}: cannot load {}: {}", name, path, e.what());
    }
    return model;
  };
  // the driver of the model, set up for the LAMMPS types; GRACE: its UQ head
  auto build = [&](auto model, auto make) {
    const std::vector<int> type2mu = GraceRuntime::type2mu(lmp, name, map, elements, model->element_names);
    has_uq = false;
    try {
      auto d = make(model);
      d->setup_types(ntypes, type2mu);
      if constexpr (!std::is_same_v<std::decay_t<decltype(d)>, DP>) {
        if (model->has_uq) {
          try {
            d->enable_uq();
            has_uq = true;
          } catch (const grace::Error &e) {
            if (comm->me == 0)
              utils::logmesg(lmp, "[{}] UQ head not usable, running without UQ: {}\n", name, e.what());
          }
        }
      } else {    // the extrapolation grade with an active set (.asi, or in the npz)
        has_uq = model->has_uq;
        if (comm->me == 0 && has_uq)
          utils::logmesg(lmp, "[{}] active set found: extrapolation grade available\n", name);
      }
      drv = std::move(d);
    } catch (const std::exception &e) {    // no LAMMPS error is raised in here
      error->all(FLERR, "pair_style {}: model setup failed: {}", name, e.what());
    }
  };
  const grace::Logger logger{GraceRuntime::log, lmp};
  grace::FusedSettings fs;
  fs.spline_n = set.spline;
  fs.ranks_per_node = GraceRuntime::ranks_sharing(world, "");    // the auto memory budgets share the host memory
  switch (mf.kind) {
    case grace::ModelKind::Grace1L:
      build(load(std::make_shared<grace::GRACE1LModel>()),
            [&](auto m) { return std::make_unique<grace::Grace1LFused>(m, fs, logger); });
      break;
    case grace::ModelKind::Grace2L:
      // store budget: libgrace's auto (half the available host memory, divided by the ranks of the node)
      fs.store_l1 = set.store_l1;
      build(load(std::make_shared<grace::GRACE2LModel>()),
            [&](auto m) { return std::make_unique<grace::Grace2LFused>(m, fs, logger); });
      break;
    case grace::ModelKind::Grace3L:
      for (int k = 0; k < 3; k++) fs.spline_n_layer[k] = set.spline_layer[k];
      fs.block = set.block;
      fs.spline_mb = set.spline_mb;
      fs.store_fwd = set.store_fwd != 0;    // auto: yes
      fs.fc_fold = set.fc_fold;
      build(load(std::make_shared<grace::GRACE3LModel>()),
            [&](auto m) { return std::make_unique<grace::Grace3LFused>(m, fs, logger); });
      break;
    case grace::ModelKind::Pace:
      build(load(std::make_shared<grace::PaceModel>()),
            [&](auto m) { return std::make_unique<grace::PaceCpu>(m, set.spline, logger); });
      break;
  }
  if (comm->me == 0 && !isa_warned)
    if (const char *isa = grace::simd::isa_hint())
      error->warning(FLERR,
                     "pair_style {}: the CPU supports {}, but LAMMPS was compiled without it; recompile with "
                     "-march=native (or the -march of this CPU) for full speed",
                     name, isa);
  isa_warned = true;
}

double PairGRACERuntime::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");
  double rc = 0.0;
  with_driver([&](auto &d) { rc = d->cutoff(i, j); });
  cutsq[i][j] = cutsq[j][i] = rc * rc;
  return rc;
}

/* ---------------------------------------------------------------------- */

void PairGRACERuntime::init_style()
{
  if (std::holds_alternative<std::monostate>(drv))
    error->all(FLERR, "pair_style {}: pair_coeff not set", name);
  if (atom->tag_enable == 0) error->all(FLERR, "pair_style {} requires atom IDs", name);
  if (force->newton_pair == 0) error->all(FLERR, "pair_style {} requires newton pair on", name);
  neighbor->add_request(this, NeighConst::REQ_FULL);
  dc.init(lmp, name.c_str());
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#else
  nthreads = 1;
#endif
  comm_forward = comm_reverse = 0;    // 1L, PACE: no layer messages
  if (auto *d = std::get_if<D2>(&drv)) comm_forward = comm_reverse = msg2_doubles(**d);
  if (auto *d = std::get_if<D3>(&drv)) {
    comm_forward = comm_reverse = std::max((*d)->message_width(0), (*d)->message_width(1)) / 2;
    if (comm->me == 0)
      utils::logmesg(lmp, "[{}] GEMM: {}, block {}, {} OpenMP threads\n", name, (*d)->gemm_backend(), set.block,
                     nthreads);
  }
}

/* ---------------------------------------------------------------------- */

void PairGRACERuntime::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);
  const int nlocal = atom->nlocal, nall = nlocal + atom->nghost;
  grace::LammpsBonds src{atom->x, atom->type, list->numneigh, list->firstneigh, NEIGHMASK};
  grace::Outputs out;
  out.f = atom->f[0];
  out.eatom = eflag_atom ? eatom : nullptr;
  out.vatom = vflag_atom ? vatom : nullptr;
  grace::ComputeFlags fl;
  if (dc.eq) {    // the masks of ghost atoms come with the border comm
    fl.couple.eq = dc.eq;
    fl.couple.lambda = dc.lambda;
    fl.couple.mask = atom->mask;
    fl.couple.bit1 = dc.bit1;
    fl.couple.bit2 = dc.bit2;
  }
  fl.uq = uq.any(has_uq);
  if (fl.uq) {
    uq.grow(nlocal);
    // the driver writes the atoms of its list only (pair_style hybrid: other atoms are not in it) -> 0 elsewhere
    std::fill(uq.gamma, uq.gamma + nlocal, 0.0);
    out.gamma = uq.gamma;
    if (!pace) {    // PACE: the extrapolation grade only
      std::fill(uq.sigma, uq.sigma + nlocal, 0.0);
      std::fill(uq.cluster, uq.cluster + nlocal, 0.0);
      out.atomic_sigma = uq.sigma;
      out.gmm_cluster = uq.cluster;
    }
  }
  fl.eflag_atom = eflag_atom;
  fl.vflag_atom = vflag_atom;
  fl.energy_only = eflag_only;    // LAMMPS ENERGY_ONLY (e.g. Monte Carlo trial moves): no forces, no virials
  LammpsExchange ex;
  ex.pair = this;
  ex.comm = comm;
  try {
    with_driver([&](auto &d) {
      using D = std::decay_t<decltype(d)>;
      if constexpr (std::is_same_v<D, D2> || std::is_same_v<D, D3>)
        d->compute(list->inum, list->ilist, nall, atom->nmax, src, ex, out, fl, nthreads);
      else
        d->compute(list->inum, list->ilist, nall, atom->nmax, src, out, fl, nthreads);
    });
  } catch (const LAMMPSException &) {    // from the layer-message comm
    throw;
  } catch (const std::exception &e) {
    error->one(FLERR, "pair_style {}: {}", name, e.what());
  }
  dc.pvector[0] = out.dedl;
  if (eflag_global) eng_vdwl += out.energy;

  if (vflag_fdotr && !fl.energy_only) virial_fdotr_compute();
  else if (vflag_global && vflag_atom && !fl.energy_only) {    // fdotr is off when per-atom virials are requested
    for (int i = 0; i < nall; i++)
      for (int k = 0; k < 6; k++) virial[k] += vatom[i][k];
  }
}

/* ----------------------------------------------------------------------
   layer messages: 2L rows without their channel padding, 3L rows of message_width(comm_layer) floats; both as raw
   fp32 inside the double buffer
------------------------------------------------------------------------- */

int PairGRACERuntime::pack_forward_comm(int n, int *list_, double *buf, int, int *)
{
  if (auto *dp = std::get_if<D2>(&drv)) {
    auto &d = **dp;
    const int nI = d.message_rows(), CR = d.message_channels(), PR = d.message_stride(), per = msg2_doubles(d);
    const size_t row = (size_t) CR * sizeof(float);
    const float *I_glob = d.message();
    for (int q = 0; q < n; q++) {
      const float *src = I_glob + (size_t) list_[q] * nI * PR;
      char *dst = reinterpret_cast<char *>(buf + (size_t) q * per);
      for (int fn = 0; fn < nI; fn++) std::memcpy(dst + fn * row, src + fn * PR, row);
    }
    return n * per;
  }
  auto &d = *std::get<D3>(drv);
  const int W = d.message_width(comm_layer);
  const float *src = d.message(comm_layer);
  float *fb = reinterpret_cast<float *>(buf);
#pragma omp parallel for schedule(static)
  for (int q = 0; q < n; q++) std::memcpy(fb + (size_t) q * W, src + (size_t) list_[q] * W, W * sizeof(float));
  return n * (W / 2);
}

void PairGRACERuntime::unpack_forward_comm(int n, int first, double *buf)
{
  if (auto *dp = std::get_if<D2>(&drv)) {
    auto &d = **dp;
    const int nI = d.message_rows(), CR = d.message_channels(), PR = d.message_stride(), per = msg2_doubles(d);
    const size_t row = (size_t) CR * sizeof(float);
    float *I_glob = d.message();
    for (int q = 0; q < n; q++) {
      float *dst = I_glob + (size_t) (first + q) * nI * PR;
      const char *src = reinterpret_cast<const char *>(buf + (size_t) q * per);
      for (int fn = 0; fn < nI; fn++) {
        std::memcpy(dst + fn * PR, src + fn * row, row);
        for (int c = CR; c < PR; c++) dst[fn * PR + c] = 0.0f;
      }
    }
    return;
  }
  auto &d = *std::get<D3>(drv);
  const int W = d.message_width(comm_layer);
  float *dst = d.message(comm_layer);
  const float *fb = reinterpret_cast<const float *>(buf);
#pragma omp parallel for schedule(static)
  for (int q = 0; q < n; q++) std::memcpy(dst + (size_t) (first + q) * W, fb + (size_t) q * W, W * sizeof(float));
}

int PairGRACERuntime::pack_reverse_comm(int n, int first, double *buf)
{
  if (auto *dp = std::get_if<D2>(&drv)) {
    auto &d = **dp;
    const int nI = d.message_rows(), CR = d.message_channels(), PR = d.message_stride(), per = msg2_doubles(d);
    const size_t row = (size_t) CR * sizeof(float);
    const float *gI_sum = d.message_grad();
    for (int q = 0; q < n; q++) {
      const float *src = gI_sum + (size_t) (first + q) * nI * PR;
      char *dst = reinterpret_cast<char *>(buf + (size_t) q * per);
      for (int fn = 0; fn < nI; fn++) std::memcpy(dst + fn * row, src + fn * PR, row);
    }
    return n * per;
  }
  auto &d = *std::get<D3>(drv);
  const int W = d.message_width(comm_layer);
  const float *src = d.message_grad(comm_layer);
  float *fb = reinterpret_cast<float *>(buf);
#pragma omp parallel for schedule(static)
  for (int q = 0; q < n; q++) std::memcpy(fb + (size_t) q * W, src + (size_t) (first + q) * W, W * sizeof(float));
  return n * (W / 2);
}

void PairGRACERuntime::unpack_reverse_comm(int n, int *list_, double *buf)
{
  // an owner can appear more than once in the list (small periodic cells): serial
  if (auto *dp = std::get_if<D2>(&drv)) {
    auto &d = **dp;
    const int nI = d.message_rows(), CR = d.message_channels(), PR = d.message_stride(), per = msg2_doubles(d);
    float *gI_sum = d.message_grad();
    for (int q = 0; q < n; q++) {
      float *dst = gI_sum + (size_t) list_[q] * nI * PR;
      const float *src = reinterpret_cast<const float *>(buf + (size_t) q * per);    // rows of CR floats, as 3L
      for (int fn = 0; fn < nI; fn++)
        for (int c = 0; c < CR; c++) dst[fn * PR + c] += src[fn * CR + c];
    }
    return;
  }
  auto &d = *std::get<D3>(drv);
  const int W = d.message_width(comm_layer);
  float *dst = d.message_grad(comm_layer);
  const float *fb = reinterpret_cast<const float *>(buf);
  for (int q = 0; q < n; q++) {
    float *dq = dst + (size_t) list_[q] * W;
    const float *s = fb + (size_t) q * W;
#pragma omp simd
    for (int c = 0; c < W; c++) dq[c] += s[c];
  }
}

/* ---------------------------------------------------------------------- */

double PairGRACERuntime::memory_usage()
{
  double b = Pair::memory_usage();
  if (auto *d = std::get_if<D2>(&drv)) b += (*d)->memory_bytes();
  if (auto *d = std::get_if<D3>(&drv)) b += (*d)->memory_bytes();
  return b;
}

void *PairGRACERuntime::extract(const char *str, int &dim)
{
  dim = 0;
  if (!pace)    // PACE: no decouple/4d
    if (void *p = dc.extract(str, dim)) return p;
  return has_uq ? uq.flag(str) : nullptr;
}

void *PairGRACERuntime::extract_peratom(const char *str, int &ncol)
{
  ncol = 0;
  return has_uq ? uq.peratom(str) : nullptr;
}

void PairGRACERuntime::modify_params(int narg, char **arg)
{
  GraceGlue::Decouple::modify_params(this, narg, arg);
}

int PairGRACERuntime::modify_param(int narg, char **arg)
{
  return pace ? 0 : dc.modify_param(lmp, name.c_str(), narg, arg);    // PACE: no decouple/4d
}
#endif    // NO_GRACE_RUNTIME
