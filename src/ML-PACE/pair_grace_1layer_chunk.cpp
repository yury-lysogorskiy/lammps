/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/*
Copyright 2025 Yury Lysogorskiy^1,  Anton Bochkarev^1, Ralf Drautz^1

^1: Ruhr-University Bochum, Bochum, Germany
*/

#ifndef NO_GRACE_TF
// #define GRACE_CHUNK_DEBUG
// #define GRACE_PROFILE

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

namespace LAMMPS_NS {

struct GRACE1LayerChunkImpl {
  cppflow::model *model = nullptr;
  GRACE::GracePaddingDimension atom_padding;
  GRACE::GracePaddingDimension real_atom_padding;
  GRACE::GracePaddingDimension neighbor_padding;

  std::map<std::string, cppflow::TensorInfo> compute_inputs_sig;

  // Persistent buffers for efficient indexing and marshalling
  std::vector<int> global_to_chunk_map;   // size nall, init to -1
  std::vector<int> chunk_to_global_map;   // size chunksize + buffer

  std::vector<int32_t> atomic_mu_i_local; // size n_real_padded
  std::vector<int32_t> ind_i;
  std::vector<int32_t> ind_j;
  std::vector<int32_t> mu_i;
  std::vector<int32_t> mu_j;
  std::vector<double> bond_vector;
  std::vector<int32_t> map_atoms_to_structure;

  GRACE1LayerChunkImpl() = default;
  ~GRACE1LayerChunkImpl() { delete model; }
};

}    // namespace LAMMPS_NS

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairGRACE1LayerChunk::PairGRACE1LayerChunk(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  impl = new GRACE1LayerChunkImpl;

  scale = nullptr;
  chunksize = 4096;

  total_timer.init();
  data_timer.init();
  model_timer.init();
  tp_timer.init();

  no_virial_fdotr_compute = 1;
  flag_compute_energy_only = 0;
}

/* ---------------------------------------------------------------------- */

PairGRACE1LayerChunk::~PairGRACE1LayerChunk()
{
  if (copymode) return;

  delete impl;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(scale);
  }
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

/* ---------------------------------------------------------------------- */

void PairGRACE1LayerChunk::settings(int narg, char **arg)
{
  if (narg < 0) utils::missing_cmd_args(FLERR, "pair_style grace/1layer/chunk", error);

  if (strcmp("metal", update->unit_style) != 0)
    error->all(FLERR, "GRACE potentials require 'metal' units");

  if (comm->me == 0) utils::logmesg(lmp, "[GRACE] TF version: {}\n", TF_Version());

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "chunksize") == 0) {
      chunksize = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "padding") == 0) {
      neigh_padding_fraction = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "pad_verbose") == 0) {
      pad_verbose = true;
      iarg += 1;
    } else if (strcmp(arg[iarg], "max_number_of_reduction") == 0) {
      max_number_of_reduction = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "reduce_padding") == 0) {
      reducing_neigh_padding_fraction = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "deny_energy_only_calc") == 0) {
      deny_energy_only_calc = true;
      iarg += 1;
    } else
      error->all(FLERR, "[GRACE] Unknown pair_style grace/1layer/chunk keyword: {}", arg[iarg]);
  }

  do_padding = (neigh_padding_fraction > 0);
  
  // Configure padding helpers
  impl->atom_padding.enabled = do_padding;
  impl->atom_padding.padding_fraction = neigh_padding_fraction;
  impl->atom_padding.reduction_threshold_fraction = reducing_neigh_padding_fraction;
  impl->atom_padding.max_reductions = max_number_of_reduction;
  impl->atom_padding.verbose = pad_verbose;

  impl->real_atom_padding.enabled = do_padding;
  impl->real_atom_padding.padding_fraction = neigh_padding_fraction;
  impl->real_atom_padding.reduction_threshold_fraction = reducing_neigh_padding_fraction;
  impl->real_atom_padding.max_reductions = max_number_of_reduction;
  impl->real_atom_padding.verbose = pad_verbose;

  impl->neighbor_padding.enabled = do_padding;
  impl->neighbor_padding.padding_fraction = neigh_padding_fraction;
  impl->neighbor_padding.reduction_threshold_fraction = reducing_neigh_padding_fraction;
  impl->neighbor_padding.max_reductions = max_number_of_reduction;
  impl->neighbor_padding.verbose = pad_verbose;

  centroidstressflag = CENTROID_AVAIL;
}

/* ---------------------------------------------------------------------- */

void PairGRACE1LayerChunk::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  map_element2type(narg - 3, arg + 3);
  
  auto potential_path = std::string(arg[2]);

  if (impl->model) {
      delete impl->model;
      impl->model = nullptr;
  }
  
  if (comm->me == 0) utils::logmesg(lmp, "[GRACE] Loading {}\n", potential_path);

  const std::vector<uint8_t> config_bytes = {
      0x32, 0x05, 0x82, 0x01, 0x02, 0x18, 0x00
  };

  impl->model = new cppflow::model(potential_path, config_bytes);
  if (comm->me == 0) std::cerr << "[GRACE] model loaded" << std::endl;

  YAML_PACE::Node metadata_yaml = YAML_PACE::LoadFile(potential_path + "/metadata.yaml");
  elements_name = metadata_yaml["chemical_symbols"].as<std::vector<std::string>>();
  nelements = (int) elements_name.size();
  for (int mu = 0; mu < nelements; mu++) { elements_to_index_map[elements_name.at(mu)] = mu; }
  cutoff = metadata_yaml["cutoff"].as<double>();

  if (metadata_yaml["cutoff_matrix"]) {
    cutoff_matrix = metadata_yaml["cutoff_matrix"].as<std::vector<std::vector<double>>>();
    is_custom_cutoffs = true;
  }

  const int ntypes = atom->ntypes;
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
    cutoff_matrix_per_lammps_type.resize(ntypes + 1, std::vector<double>(ntypes + 1));
    for (int i = 1; i <= ntypes; i++) {
      for (int j = 1; j <= ntypes; j++) {
        cutoff_matrix_per_lammps_type[i][j] = cutoff_matrix[element_type_mapping[i]][element_type_mapping[j]];
      }
    }
  }

  if (impl->model->has_signature("compute")) {
    compute_function_name = "compute";
  } else if (impl->model->has_signature("serving_default")) {
    compute_function_name = "serving_default";
  }
  impl->compute_inputs_sig = impl->model->signatures.at(compute_function_name).inputs;

  if (impl->model->has_signature(COMPUTE_ENERGY_ONLY_KEY)) {
    compute_energy_only_function_name = COMPUTE_ENERGY_ONLY_KEY;
    has_compute_energy_only = true;
  }

  this->DEFAULT_INPUT_PREFIX = this->compute_function_name + "_";
  has_map_atoms_to_structure_op = impl->model->has_graph_input(DEFAULT_INPUT_PREFIX + "map_atoms_to_structure");
  has_nstruct_total_op = impl->model->has_graph_input(DEFAULT_INPUT_PREFIX + "n_struct_total");
  has_mu_i_op = impl->model->has_graph_input(DEFAULT_INPUT_PREFIX + "mu_i");
  has_batch_tot_nat = impl->model->has_graph_input(DEFAULT_INPUT_PREFIX + "batch_tot_nat");
  has_atomic_mu_i_local = impl->model->has_graph_input(DEFAULT_INPUT_PREFIX + "atomic_mu_i_local");
}

/* ---------------------------------------------------------------------- */

void PairGRACE1LayerChunk::init_style()
{
  if (atom->tag_enable == 0) error->all(FLERR, "Pair style grace/1layer/chunk requires atom IDs");
  if (force->newton_pair == 0) error->all(FLERR, "Pair style grace/1layer/chunk requires newton pair on");

  neighbor->add_request(this, NeighConst::REQ_FULL);

  if (atom->map_style == Atom::MAP_NONE) {
    atom->map_init();
    atom->map_set();
  }
}

/* ---------------------------------------------------------------------- */

double PairGRACE1LayerChunk::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");
  scale[j][i] = scale[i][j];
  if (is_custom_cutoffs) return cutoff_matrix_per_lammps_type[i][j];
  return cutoff;
}

/* ---------------------------------------------------------------------- */

void *PairGRACE1LayerChunk::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str, "compute_energy_only") == 0) return (void *) &flag_compute_energy_only;
  dim = 2;
  if (strcmp(str, "scale") == 0) return (void *) scale;
  return nullptr;
}

/* ---------------------------------------------------------------------- */

void PairGRACE1LayerChunk::compute(int eflag, int vflag)
{
  total_timer.init(); data_timer.init(); model_timer.init();
  total_timer.start();
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  int *ilist = list->ilist;
  int inum = list->inum;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  bool do_energy_only = flag_compute_energy_only & !deny_energy_only_calc;
  auto compute_inputs_sig = impl->compute_inputs_sig;

  // Initialize/resize global-to-chunk map
  if (impl->global_to_chunk_map.size() < (size_t)nall) {
    impl->global_to_chunk_map.assign(nall, -1);
  }
  if (comm->me == 0) std::cerr << "[GRACE-CHUNK] compute() starts, inum=" << inum << ", nall=" << nall << std::endl;

  int chunk_offset = 0;
  int chunk_idx = 0;

  while (chunk_offset < inum) {
    int current_chunk_size = std::min(chunksize, inum - chunk_offset);
    data_timer.start();

    // -- Phase 1: Topology Discovery --
    int node_counter = 0;
    impl->chunk_to_global_map.clear();

    // Map real atoms in chunk
    for (int ii = 0; ii < current_chunk_size; ++ii) {
      int i = ilist[chunk_offset + ii];
      impl->global_to_chunk_map[i] = node_counter;
      impl->chunk_to_global_map.push_back(i);
      node_counter++;
    }
    int n_real = node_counter;

    // Scan neighbors to discover pseudo-ghosts and bonds
    int n_bonds_real = 0;
    for (int ii = 0; ii < current_chunk_size; ++ii) {
      int i = ilist[chunk_offset + ii];
      int i_type = type[i];
      int i_neigh_num = numneigh[i];
      int *i_neigh_list = firstneigh[i];
      double xtmp = x[i][0];
      double ytmp = x[i][1];
      double ztmp = x[i][2];

      for (int jj = 0; jj < i_neigh_num; ++jj) {
        int j = i_neigh_list[jj] & NEIGHMASK;
        int j_type = type[j];
        
        double delx = xtmp - x[j][0];
        double dely = ytmp - x[j][1];
        double delz = ztmp - x[j][2];
        double rsq = delx*delx + dely*dely + delz*delz;
        double cutsq_ij = is_custom_cutoffs ? cutoff_matrix_per_lammps_type[i_type][j_type] * cutoff_matrix_per_lammps_type[i_type][j_type] : cutoff * cutoff;

        if (rsq < cutsq_ij) {
          if (impl->global_to_chunk_map[j] == -1) {
            impl->global_to_chunk_map[j] = node_counter;
            impl->chunk_to_global_map.push_back(j);
            node_counter++;
          }
          n_bonds_real++;
        }
      }
    }
    int n_nodes_in_chunk = node_counter;

    // Padding
    int n_nodes_padded = impl->atom_padding.update(n_nodes_in_chunk);
    int n_real_padded = impl->real_atom_padding.update(n_real);
    int n_bonds_padded = impl->neighbor_padding.update(n_bonds_real);

#ifdef GRACE_CHUNK_DEBUG
    {
      utils::logmesg(lmp, "[CHUNK-DBG] chunk_idx={}, offset={}, size={}\n", chunk_idx, chunk_offset, current_chunk_size);
      utils::logmesg(lmp, "[CHUNK-DBG] n_nodes={} (n_real={}), n_padded={}, n_real_padded={}\n", n_nodes_in_chunk, n_real, n_nodes_padded, n_real_padded);
      utils::logmesg(lmp, "[CHUNK-DBG] n_bonds_real={}, n_bonds_padded={}\n", n_bonds_real, n_bonds_padded);
    }
#endif

    // -- Phase 2: Build Tensors --
    impl->atomic_mu_i_local.assign(n_real_padded, element_type_mapping[type[0]]);
    impl->ind_i.assign(n_bonds_padded, n_nodes_padded - 1);
    impl->ind_j.assign(n_bonds_padded, n_nodes_padded - 1);
    impl->mu_i.assign(n_bonds_padded, 0);
    impl->mu_j.assign(n_bonds_padded, 0);
    impl->bond_vector.assign(3 * n_bonds_padded, 1e6);
    impl->map_atoms_to_structure.assign(n_nodes_padded, 0);

    // Fill real nodes data
    for (int k = 0; k < n_real; ++k) {
      int g_idx = impl->chunk_to_global_map[k];
      impl->atomic_mu_i_local[k] = element_type_mapping[type[g_idx]];
    }

    // Fill real bonds data
    int bond_idx = 0;
    for (int ii = 0; ii < current_chunk_size; ++ii) {
      int i = ilist[chunk_offset + ii];
      int i_type = type[i];
      int i_neigh_num = numneigh[i];
      int *i_neigh_list = firstneigh[i];
      int i_chunk = impl->global_to_chunk_map[i];

      for (int jj = 0; jj < i_neigh_num; ++jj) {
        int j = i_neigh_list[jj] & NEIGHMASK;
        int j_type = type[j];
        double dx = x[j][0] - x[i][0];
        double dy = x[j][1] - x[i][1];
        double dz = x[j][2] - x[i][2];
        double rsq = dx*dx + dy*dy + dz*dz;
        double cutsq_ij = is_custom_cutoffs ? cutoff_matrix_per_lammps_type[i_type][j_type] * cutoff_matrix_per_lammps_type[i_type][j_type] : cutoff * cutoff;

        if (rsq < cutsq_ij) {
          int j_chunk = impl->global_to_chunk_map[j];
          impl->ind_i[bond_idx] = i_chunk;
          impl->ind_j[bond_idx] = j_chunk;
          impl->mu_i[bond_idx] = element_type_mapping[type[i]];
          impl->mu_j[bond_idx] = element_type_mapping[type[j]];
          impl->bond_vector[3 * bond_idx + 0] = dx;
          impl->bond_vector[3 * bond_idx + 1] = dy;
          impl->bond_vector[3 * bond_idx + 2] = dz;
          bond_idx++;
        }
      }
    }

    std::vector<std::tuple<std::string, cppflow::tensor>> inputs;
    inputs.emplace_back(compute_inputs_sig.at("atomic_mu_i").name, cppflow::tensor(impl->atomic_mu_i_local, {n_real_padded}));
    if (has_atomic_mu_i_local) 
      inputs.emplace_back(compute_inputs_sig.at("atomic_mu_i_local").name, cppflow::tensor(impl->atomic_mu_i_local, {n_real_padded}));
    inputs.emplace_back(compute_inputs_sig.at("ind_i").name, cppflow::tensor(impl->ind_i, {n_bonds_padded}));
    inputs.emplace_back(compute_inputs_sig.at("ind_j").name, cppflow::tensor(impl->ind_j, {n_bonds_padded}));
    if (has_mu_i_op)
      inputs.emplace_back(compute_inputs_sig.at("mu_i").name, cppflow::tensor(impl->mu_i, {n_bonds_padded}));
    inputs.emplace_back(compute_inputs_sig.at("mu_j").name, cppflow::tensor(impl->mu_j, {n_bonds_padded}));
    inputs.emplace_back(compute_inputs_sig.at("bond_vector").name, cppflow::tensor(impl->bond_vector, {n_bonds_padded, 3}));
    inputs.emplace_back(compute_inputs_sig.at("batch_tot_nat_real").name, cppflow::tensor(std::vector<int32_t>{n_real}, {}));
    if (has_batch_tot_nat) inputs.emplace_back(compute_inputs_sig.at("batch_tot_nat").name, cppflow::tensor(std::vector<int32_t>{n_nodes_padded}, {}));
    if (has_nstruct_total_op) inputs.emplace_back(compute_inputs_sig.at("n_struct_total").name, cppflow::tensor(std::vector<int32_t>{1}, {}));
    if (has_map_atoms_to_structure_op) inputs.emplace_back(compute_inputs_sig.at("map_atoms_to_structure").name, cppflow::tensor(impl->map_atoms_to_structure, {n_nodes_padded}));

#ifdef GRACE_CHUNK_DEBUG
    GRACE::print_tf_inputs(inputs, comm->me, lmp, false);
#endif
    data_timer.stop();
    model_timer.start();

    std::vector<std::string> output_names;
    auto sig_outputs = impl->model->signatures.at(compute_function_name).outputs;
    output_names.push_back(sig_outputs.at("atomic_energy").name);
    output_names.push_back(sig_outputs.at("z_pair_f").name);

    auto outputs = impl->model->operator()(inputs, output_names);
    model_timer.stop();
    data_timer.start();

    const double *e_data = static_cast<const double *>(TF_TensorData(outputs[0].get_tensor().get()));
    const double *f_data = static_cast<const double *>(TF_TensorData(outputs[1].get_tensor().get()));

    // -- Phase 4: Scattering --
    bond_idx = 0;
    for (int ii = 0; ii < current_chunk_size; ++ii) {
      int i = ilist[chunk_offset + ii];
      int i_type = type[i];
      double i_scale = scale[i_type][i_type];
      
      // Energy
      if (eflag_either) {
        double evdwl = i_scale * e_data[ii];
        ev_tally_full(i, 2.0 * evdwl, 0.0, 0.0, 0.0, 0.0, 0.0);
      }

      // Forces and Virial
      int i_neigh_num = numneigh[i];
      int *i_neigh_list = firstneigh[i];
      for (int jj = 0; jj < i_neigh_num; ++jj) {
        int j = i_neigh_list[jj] & NEIGHMASK;
        int j_type = type[j];
        double dx = x[j][0] - x[i][0];
        double dy = x[j][1] - x[i][1];
        double dz = x[j][2] - x[i][2];
        double rsq = dx*dx + dy*dy + dz*dz;
        double cutsq_ij = is_custom_cutoffs ? cutoff_matrix_per_lammps_type[i_type][j_type] * cutoff_matrix_per_lammps_type[i_type][j_type] : cutoff * cutoff;

        if (rsq < cutsq_ij) {
          double fij[3];
          fij[0] = -i_scale * f_data[3 * bond_idx + 0];
          fij[1] = -i_scale * f_data[3 * bond_idx + 1];
          fij[2] = -i_scale * f_data[3 * bond_idx + 2];

          f[i][0] += fij[0]; f[i][1] += fij[1]; f[i][2] += fij[2];
          f[j][0] -= fij[0]; f[j][1] -= fij[1]; f[j][2] -= fij[2];

          if (vflag_either) {
            ev_tally_xyz(i, j, nlocal, 1, 0.0, 0.0, fij[0], fij[1], fij[2], -dx, -dy, -dz);
            if (cvflag_atom) {
              const double dx_val = -dx;
              const double dy_val = -dy;
              const double dz_val = -dz;
              cvatom[i][0] += 0.5 * dx_val * fij[0];
              cvatom[i][1] += 0.5 * dy_val * fij[1];
              cvatom[i][2] += 0.5 * dz_val * fij[2];
              cvatom[i][3] += 0.5 * dx_val * fij[1];
              cvatom[i][4] += 0.5 * dx_val * fij[2];
              cvatom[i][5] += 0.5 * dy_val * fij[2];
              cvatom[i][6] += 0.5 * dy_val * fij[0];
              cvatom[i][7] += 0.5 * dz_val * fij[0];
              cvatom[i][8] += 0.5 * dz_val * fij[1];

              cvatom[j][0] += 0.5 * dx_val * fij[0];
              cvatom[j][1] += 0.5 * dy_val * fij[1];
              cvatom[j][2] += 0.5 * dz_val * fij[2];
              cvatom[j][3] += 0.5 * dx_val * fij[1];
              cvatom[j][4] += 0.5 * dx_val * fij[2];
              cvatom[j][5] += 0.5 * dy_val * fij[2];
              cvatom[j][6] += 0.5 * dy_val * fij[0];
              cvatom[j][7] += 0.5 * dz_val * fij[0];
              cvatom[j][8] += 0.5 * dz_val * fij[1];
            }
          }
          bond_idx++;
        }
      }
    }

    // -- Phase 5: Cleanup --
    for (int k = 0; k < n_nodes_in_chunk; ++k) {
      impl->global_to_chunk_map[impl->chunk_to_global_map[k]] = -1;
    }

    data_timer.stop();
    chunk_offset += current_chunk_size;
    chunk_idx++;
  }

  total_timer.stop();

#ifdef GRACE_PROFILE
  if (inum > 0) {
      double d_t = data_timer.as_microseconds();
      double m_t = model_timer.as_microseconds();
      double total_t = total_timer.as_microseconds();
      auto pct = [&](double t) { return (total_t > 0) ? (t / total_t * 100.0) : 0.0; };

      fprintf(stderr, "[GRACE-PROFILE] [Rank %d] Timings (mcs): Data: %.1f (%.1f%%) | Model: %.1f (%.1f%%) | Total: %.1f\n",
              comm->me, d_t, pct(d_t), m_t, pct(m_t), total_t);
  }
#endif
}

#endif
