/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef NO_GRACE_TF

#include "pair_grace_1layer_chunk.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include "yaml-cpp/yaml.h"
#include <cstring>
#include <numeric>

#include "utils_grace.h"
#include "utils_pace.h"

// CppFlow headers
#include <cppflow/model.h>
#include <cppflow/ops.h>
#include <cppflow/tensor.h>
#include <string>
#include <tensorflow/c/c_api.h>
#include <unistd.h>

using namespace LAMMPS_NS;
using namespace MathConst;

namespace LAMMPS_NS {

struct GRACE1LayerChunkImpl {
  GRACE1LayerChunkImpl() : model(nullptr) {}
  ~GRACE1LayerChunkImpl() { delete model; }

  cppflow::model *model;
  GRACE::GracePaddingDimension atom_padding;
  GRACE::GracePaddingDimension neighbor_padding;

  // Translation maps
  std::vector<int> global_to_chunk_map;
  std::vector<int> chunk_to_global_map;

  // Persistent vectors to avoid re-allocation
  std::vector<int32_t> ind_i;
  std::vector<int32_t> ind_j;
  std::vector<double> bond_vector;
  std::vector<int32_t> mu_i;
  std::vector<int32_t> mu_j;
  std::vector<int32_t> atomic_mu_i;
  std::vector<int32_t> atomic_mu_i_local;
};

}    // namespace LAMMPS_NS

/* ---------------------------------------------------------------------- */

PairGRACE1LayerChunk::PairGRACE1LayerChunk(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  graceimpl = new GRACE1LayerChunkImpl;
  scale = nullptr;
  chunksize = 4096;

  total_timer.init();
  data_timer.init();
  tp_timer.init();
}

/* ---------------------------------------------------------------------- */

PairGRACE1LayerChunk::~PairGRACE1LayerChunk()
{
  if (copymode) return;

  delete graceimpl;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(scale);
    delete[] map;
  }
}

/* ---------------------------------------------------------------------- */

void PairGRACE1LayerChunk::compute(int eflag, int vflag)
{
  total_timer.start();
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  if (comm->me == 0)
    utils::logmesg(lmp, "[GRACE/DEBUG] Compute started: nlocal={}, nall={}, eflag={}, vflag={}\n",
                   nlocal, nall, eflag, vflag);

  // Ensure maps are correctly sized
  if (graceimpl->global_to_chunk_map.size() < (size_t) nall) {
    if (comm->me == 0)
      utils::logmesg(lmp, "[GRACE/DEBUG] Resizing global_to_chunk_map to {}\n", nall);
    graceimpl->global_to_chunk_map.assign(nall, -1);
  }

  const std::string forward_name = compute_function_name;
  if (graceimpl->model->signatures.count(forward_name) == 0) {
    error->all(FLERR, "[GRACE-ERROR] Signature '{}' not found in model", forward_name);
  }
  auto sig = graceimpl->model->signatures.at(forward_name);

  // Chunking Loop
  for (int chunk_offset = 0; chunk_offset < nlocal; chunk_offset += chunksize) {
    data_timer.start();

    int current_chunk_size = std::min(chunksize, nlocal - chunk_offset);
    if (comm->me == 0)
      utils::logmesg(lmp, "[GRACE/DEBUG] Processing chunk: offset={}, size={}\n", chunk_offset,
                     current_chunk_size);
    graceimpl->chunk_to_global_map.clear();

    // Mapping Step A
    // Phase 1: Real atoms in chunk
    for (int k = 0; k < current_chunk_size; ++k) {
      int i = chunk_offset + k;
      graceimpl->global_to_chunk_map[i] = k;
      graceimpl->chunk_to_global_map.push_back(i);
    }

    int n_real_chunk = current_chunk_size;
    int n_total_chunk = n_real_chunk;

    // Phase 2: Neighbors
    graceimpl->ind_i.clear();
    graceimpl->ind_j.clear();
    graceimpl->bond_vector.clear();
    graceimpl->mu_i.clear();
    graceimpl->mu_j.clear();

    double cutoff_sq = cutoff * cutoff;

    for (int k = 0; k < current_chunk_size; ++k) {
      int i = chunk_offset + k;
      int i_chunk = k;
      int type_i = type[i];
      double xtmp = x[i][0];
      double ytmp = x[i][1];
      double ztmp = x[i][2];

      int *jlist = firstneigh[i];
      int jnum = numneigh[i];

      for (int jj = 0; jj < jnum; ++jj) {
        int j = jlist[jj];
        j &= NEIGHMASK;

        int type_j = type[j];
        if (is_custom_cutoffs) {
          double cur_cutoff = cutoff_matrix_per_lammps_type[type_i][type_j];
          cutoff_sq = cur_cutoff * cur_cutoff;
        }

        double delx = xtmp - x[j][0];
        double dely = ytmp - x[j][1];
        double delz = ztmp - x[j][2];
        double rsq = delx * delx + dely * dely + delz * delz;

        if (rsq < cutoff_sq) {
          int j_chunk = graceimpl->global_to_chunk_map[j];
          if (j_chunk == -1) {
            j_chunk = n_total_chunk++;
            graceimpl->global_to_chunk_map[j] = j_chunk;
            graceimpl->chunk_to_global_map.push_back(j);
          }

          graceimpl->ind_i.push_back(i_chunk);
          graceimpl->ind_j.push_back(j_chunk);
          graceimpl->mu_i.push_back(element_type_mapping[type_i]);
          graceimpl->mu_j.push_back(element_type_mapping[type_j]);
          graceimpl->bond_vector.push_back(x[j][0] - x[i][0]);
          graceimpl->bond_vector.push_back(x[j][1] - x[i][1]);
          graceimpl->bond_vector.push_back(x[j][2] - x[i][2]);
        }
      }
    }

    if (comm->me == 0)
      utils::logmesg(lmp, "[GRACE/DEBUG] Mapping done: n_total_chunk={}, n_bonds_real={}\n",
                     n_total_chunk, (int) graceimpl->ind_i.size());

    int n_bonds_real = graceimpl->ind_i.size();

    // Step B: Padding
    int n_atoms_padded = graceimpl->atom_padding.update(n_total_chunk);
    int n_bonds_padded = graceimpl->neighbor_padding.update(n_bonds_real);

    // Padding fake atoms and bonds
    graceimpl->atomic_mu_i.assign(n_atoms_padded, 0);
    for (int k = 0; k < n_total_chunk; ++k) {
      graceimpl->atomic_mu_i[k] = element_type_mapping[type[graceimpl->chunk_to_global_map[k]]];
    }

    graceimpl->atomic_mu_i_local.assign(n_atoms_padded, 0);
    for (int k = 0; k < n_real_chunk; ++k) {
      graceimpl->atomic_mu_i_local[k] =
          element_type_mapping[type[graceimpl->chunk_to_global_map[k]]];
    }

    for (int k = n_bonds_real; k < n_bonds_padded; ++k) {
      graceimpl->ind_i.push_back(n_atoms_padded - 1);
      graceimpl->ind_j.push_back(n_atoms_padded - 1);
      graceimpl->mu_i.push_back(0);
      graceimpl->mu_j.push_back(0);
      graceimpl->bond_vector.push_back(1e6);
      graceimpl->bond_vector.push_back(0.0);
      graceimpl->bond_vector.push_back(0.0);
    }

    if (comm->me == 0)
      utils::logmesg(lmp, "[GRACE/DEBUG] Padding done: n_atoms_padded={}, n_bonds_padded={}\n",
                     n_atoms_padded, n_bonds_padded);

    // Step C: Tensor Construction
    std::vector<std::tuple<std::string, cppflow::tensor>> inputs;
    auto add_input = [&](const std::string &key, const cppflow::tensor &t) {
      if (sig.inputs.count(key)) inputs.emplace_back(sig.inputs.at(key).name, t);
    };

    add_input("atomic_mu_i", cppflow::tensor(graceimpl->atomic_mu_i, {n_atoms_padded}));
    add_input("atomic_mu_i_local", cppflow::tensor(graceimpl->atomic_mu_i_local, {n_atoms_padded}));
    add_input("ind_i", cppflow::tensor(graceimpl->ind_i, {n_bonds_padded}));
    add_input("ind_j", cppflow::tensor(graceimpl->ind_j, {n_bonds_padded}));
    add_input("mu_i", cppflow::tensor(graceimpl->mu_i, {n_bonds_padded}));
    add_input("mu_j", cppflow::tensor(graceimpl->mu_j, {n_bonds_padded}));
    add_input("bond_vector", cppflow::tensor(graceimpl->bond_vector, {n_bonds_padded, 3}));
    add_input("batch_tot_nat_real", cppflow::tensor(std::vector<int32_t>{n_real_chunk}, {}));

    if (comm->me == 0) utils::logmesg(lmp, "[GRACE/DEBUG] Tensors built, calling model...\n");
    data_timer.stop();
    tp_timer.start();

    // Step D: Execution
    std::vector<std::string> out_names = {sig.outputs.at("atomic_energy").name,
                                          sig.outputs.at("z_pair_f").name};
    bool do_virial = vflag_global && (sig.outputs.count("virial") > 0);
    if (do_virial) out_names.push_back(sig.outputs.at("virial").name);

    auto outputs = graceimpl->model->operator()(inputs, out_names);
    if (comm->me == 0) utils::logmesg(lmp, "[GRACE/DEBUG] Model call successful\n");
    tp_timer.stop();
    data_timer.start();

    // Output Processing
    auto e_tens = outputs[0].get_tensor();
    auto f_tens = outputs[1].get_tensor();
    const double *e_data = static_cast<double *>(TF_TensorData(e_tens.get()));
    const double *f_pair_data = static_cast<double *>(TF_TensorData(f_tens.get()));

    // Energy tally
    if (eflag_either) {
      for (int k = 0; k < n_real_chunk; k++) {
        int i_glob = graceimpl->chunk_to_global_map[k];
        ev_tally_full(i_glob, 2.0 * e_data[k], 0.0, 0.0, 0.0, 0.0, 0.0);
      }
    }

    // Forces via Newton's 3rd Law
    for (int k = 0; k < n_bonds_real; k++) {
      int i_glob = graceimpl->chunk_to_global_map[graceimpl->ind_i[k]];
      int j_glob = graceimpl->chunk_to_global_map[graceimpl->ind_j[k]];
      double fx = f_pair_data[k * 3 + 0];
      double fy = f_pair_data[k * 3 + 1];
      double fz = f_pair_data[k * 3 + 2];

      f[i_glob][0] += fx;
      f[i_glob][1] += fy;
      f[i_glob][2] += fz;
      f[j_glob][0] -= fx;
      f[j_glob][1] -= fy;
      f[j_glob][2] -= fz;

      if (vflag_either) {
        double delx = x[i_glob][0] - x[j_glob][0];
        double dely = x[i_glob][1] - x[j_glob][1];
        double delz = x[i_glob][2] - x[j_glob][2];
        ev_tally_xyz(i_glob, j_glob, nlocal, 1, 0.0, 0.0, fx, fy, fz, delx, dely, delz);
      }
    }

    // Virial global
    if (do_virial) {
      auto v_tens = outputs[2].get_tensor();
      const double *v_data = static_cast<double *>(TF_TensorData(v_tens.get()));
      if (comm->me == 0)
        utils::logmesg(lmp, "[GRACE/DEBUG] Virial from model: {} {} {} {} {} {}\n", v_data[0],
                       v_data[1], v_data[2], v_data[3], v_data[4], v_data[5]);
      virial[0] += v_data[0];
      virial[1] += v_data[1];
      virial[2] += v_data[2];
      virial[3] += v_data[3];
      virial[4] += v_data[4];
      virial[5] += v_data[5];
    }

    // Step E: Cleanup Maps
    for (int glob_id : graceimpl->chunk_to_global_map) {
      graceimpl->global_to_chunk_map[glob_id] = -1;
    }
    data_timer.stop();
  }

  total_timer.stop();
}

/* ---------------------------------------------------------------------- */

void PairGRACE1LayerChunk::settings(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR, "Illegal pair_style command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "chunksize") == 0) {
      chunksize = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Unknown pair_style grace/1layer/chunk keyword: {}", arg[iarg]);
    }
  }

  if (strcmp("metal", update->unit_style) != 0)
    error->all(FLERR, "GRACE potentials require 'metal' units");
}

/* ---------------------------------------------------------------------- */

void PairGRACE1LayerChunk::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  int ntypes = atom->ntypes;
  map_element2type(narg - 3, arg + 3);

  std::string potential_path = arg[2];
  delete graceimpl->model;

  if (comm->me == 0) utils::logmesg(lmp, "[GRACE/CHUNK] Loading {}\n", potential_path);
  graceimpl->model = new cppflow::model(potential_path);

  // Read metadata.yaml
  YAML_PACE::Node metadata_yaml = YAML_PACE::LoadFile(potential_path + "/metadata.yaml");
  elements_name = metadata_yaml["chemical_symbols"].as<std::vector<std::string>>();
  nelements = (int) elements_name.size();
  for (int mu = 0; mu < nelements; mu++) { elements_to_index_map[elements_name.at(mu)] = mu; }
  cutoff = metadata_yaml["cutoff"].as<double>();

  if (metadata_yaml["cutoff_matrix"]) {
    cutoff_matrix = metadata_yaml["cutoff_matrix"].as<vector<vector<double>>>();
    is_custom_cutoffs = true;
  }

  element_type_mapping.resize(ntypes + 1);
  for (int i = 1; i <= ntypes; i++) {
    char *elemname = arg[2 + i];
    if (strcmp(elemname, "NULL") == 0) {
      element_type_mapping[i] = -1;
    } else {
      int mu = elements_to_index_map.at(elemname);
      element_type_mapping[i] = mu;
    }
  }

  for (int i = 1; i <= ntypes; i++) {
    for (int j = i; j <= ntypes; j++) scale[i][j] = 1.0;
  }

  if (is_custom_cutoffs) {
    cutoff_matrix_per_lammps_type.resize(ntypes + 1, vector<double>(ntypes + 1));
    for (int i = 1; i <= ntypes; i++) {
      for (int j = 1; j <= ntypes; j++) {
        cutoff_matrix_per_lammps_type[i][j] =
            cutoff_matrix[element_type_mapping[i]][element_type_mapping[j]];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairGRACE1LayerChunk::init_style()
{
  if (atom->tag_enable == 0) error->all(FLERR, "Pair style grace requires atom IDs");
  if (force->newton_pair == 0) error->all(FLERR, "Pair style grace requires newton pair on");
  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

double PairGRACE1LayerChunk::init_one(int i, int j)
{
  scale[j][i] = scale[i][j];
  if (is_custom_cutoffs) return cutoff_matrix_per_lammps_type[i][j];
  return cutoff;
}

/* ---------------------------------------------------------------------- */

void PairGRACE1LayerChunk::allocate()
{
  allocated = 1;
  int n = atom->ntypes + 1;
  memory->create(setflag, n, n, "pair:setflag");
  memory->create(cutsq, n, n, "pair:cutsq");
  memory->create(scale, n, n, "pair:scale");
  map = new int[n];
}

#endif
