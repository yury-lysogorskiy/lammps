/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef NO_GRACE_TF
#define GRACE_2L_CHUNK_DEBUG

#include "pair_grace_2layer_chunk.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include "yaml-cpp/yaml.h"
#include <algorithm>
#include <cstring>
#include <tuple>

#include "utils_grace.h"
#include "utils_pace.h"

// CppFlow headers
#include <cppflow/model.h>
#include <cppflow/ops.h>
#include <cppflow/tensor.h>
#include <tensorflow/c/c_api.h>
#include <unistd.h>

namespace LAMMPS_NS {

struct GRACE2LayerChunkImpl {
  cppflow::model *model = nullptr;

  // Chunk-level padding helpers
  GRACE::GracePaddingDimension atom_padding;
  GRACE::GracePaddingDimension real_atom_padding;
  GRACE::GracePaddingDimension neighbor_padding;

  // Mapping buffers
  std::vector<int> global_to_chunk_map;    // size nall
  std::vector<int> chunk_to_global_map;    // variable

  // Chunk-local TF input arrays
  std::vector<int32_t> mu_i, mu_j, ind_i, ind_j;
  std::vector<double> bond_vector;
  std::vector<int32_t> atomic_mu_i;          // for all nodes in chunk
  std::vector<int32_t> atomic_mu_i_local;    // for real atoms in chunk

  // Persistent input tensor pools (one per TF signature)
  std::map<std::string, cppflow::tensor> fwd_l1_tensors;
  std::map<std::string, cppflow::tensor> bwd_l2_tensors;
  std::map<std::string, cppflow::tensor> bwd_l1_tensors;

  // Persistent sizes
  int n_nodes_padded = 0;
  int n_real_padded = 0;
  int n_bonds_padded = 0;
  bool graph_recompiled = false;

  GRACE2LayerChunkImpl() = default;
  ~GRACE2LayerChunkImpl() { delete model; }
};

}    // namespace LAMMPS_NS

using namespace LAMMPS_NS;

PairGRACE2LayerChunk::PairGRACE2LayerChunk(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  aceimpl = new GRACE2LayerChunkImpl();
  scale = nullptr;

  no_virial_fdotr_compute = 1;
  chunksize = 4096;
  nelements = 0;
  flag_compute_energy_only = 0;

  total_timer.init();
  data_timer.init();
  comm_timer.init();
  model1_timer.init();
  model2_timer.init();
  model3_timer.init();
}

PairGRACE2LayerChunk::~PairGRACE2LayerChunk()
{
  if (copymode) return;

  if (comm->me == 0 && total_real_atoms_processed > 0) {
    double total = total_timer.as_microseconds();
    double data = data_timer.as_microseconds();
    double comm_t = comm_timer.as_microseconds();
    double m1 = model1_timer.as_microseconds();
    double m2 = model2_timer.as_microseconds();
    double m3 = model3_timer.as_microseconds();

    auto per_atom = [&](double t) {
      return t / total_real_atoms_processed;
    };
    auto pct = [&](double t) {
      return (total > 0) ? (t / total * 100.0) : 0.0;
    };

    utils::logmesg(lmp, "[GRACE-PERF] Total real atoms processed: {:.0f}\n",
                   total_real_atoms_processed);
    utils::logmesg(lmp, "[GRACE-PERF] Total valid compute calls: {}\n", total_compute_calls);
    utils::logmesg(lmp, "[GRACE-PERF] Average atoms per step: {:.1f}\n",
                   total_real_atoms_processed / total_compute_calls);
    utils::logmesg(
        lmp,
        "[GRACE-PERF] Performance (us/atom) [%]: Total: {:.1f}, Data: {:.1f} ({:.1f}%), Comm: "
        "{:.1f} ({:.1f}%), M1: {:.1f} ({:.1f}%), M2: {:.1f} ({:.1f}%), M3: {:.1f} ({:.1f}%)\n",
        per_atom(total), per_atom(data), pct(data), per_atom(comm_t), pct(comm_t), per_atom(m1),
        pct(m1), per_atom(m2), pct(m2), per_atom(m3), pct(m3));
  }

  delete aceimpl;
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(scale);
  }
}

void PairGRACE2LayerChunk::allocate()
{
  allocated = 1;
  int n = atom->ntypes + 1;
  memory->create(setflag, n, n, "pair:setflag");
  memory->create(cutsq, n, n, "pair:cutsq");
  memory->create(scale, n, n, "pair:scale");
  map = new int[n];
}

void PairGRACE2LayerChunk::settings(int narg, char **arg)
{
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
    } else if (strcmp(arg[iarg], "debug_no_energy_only_calc") == 0) {
      debug_no_energy_only_calc = true;
      iarg += 1;
    } else
      error->all(FLERR, "[GRACE] Unknown pair_style grace keyword: {}", arg[iarg]);
  }

  do_padding = (neigh_padding_fraction > 0);

  // Apply to padding helpers in aceimpl
  auto apply_pad = [&](GRACE::GracePaddingDimension &p) {
    p.enabled = do_padding;
    p.padding_fraction = neigh_padding_fraction;
    p.reduction_threshold_fraction = reducing_neigh_padding_fraction;
    p.max_reductions = max_number_of_reduction;
    p.verbose = pad_verbose;
  };
  apply_pad(aceimpl->atom_padding);
  apply_pad(aceimpl->real_atom_padding);
  apply_pad(aceimpl->neighbor_padding);

  centroidstressflag = CENTROID_AVAIL;
}

void PairGRACE2LayerChunk::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  map_element2type(narg - 3, arg + 3);
  auto potential_path = std::string(arg[2]);

  delete aceimpl->model;
  if (comm->me == 0) utils::logmesg(lmp, "[GRACE] Loading {}\n", potential_path);

  const std::vector<uint8_t> config_bytes = {0x32, 0x05, 0x82, 0x01, 0x02, 0x18, 0x00};
  aceimpl->model = new cppflow::model(potential_path, config_bytes);

  YAML_PACE::Node metadata_yaml = YAML_PACE::LoadFile(potential_path + "/metadata.yaml");
  elements_name = metadata_yaml["chemical_symbols"].as<std::vector<std::string>>();
  nelements = static_cast<int>(elements_name.size());
  for (int mu = 0; mu < nelements; mu++) { elements_to_index_map[elements_name.at(mu)] = mu; }
  cutoff = metadata_yaml["cutoff"].as<double>();

  // Parse parallel communication features
  if (metadata_yaml["parallel_communication"]) {
    YAML_PACE::Node parallel_comm = metadata_yaml["parallel_communication"];
    for (YAML_PACE::const_iterator it = parallel_comm.begin(); it != parallel_comm.end(); ++it) {
      std::string key = it->first.as<std::string>();
      YAML_PACE::Node node = it->second;
      std::vector<int64_t> shape;
      int size = 1;
      if (node["shape"]) {
        for (const auto &dim : node["shape"]) {
          int d = dim.as<int>();
          shape.push_back(d);
          size *= d;
        }
      } else
        error->all(FLERR, "[GRACE] Feature '{}' missing 'shape' in metadata.yaml", key);

      bool non_local = node["non_local"] ? node["non_local"].as<bool>() : true;

      feature_shapes[key] = shape;
      feature_sizes[key] = size;
      feature_is_local[key] = !non_local;

      if (comm->me == 0) {
        std::string shape_str;
        for (size_t i = 0; i < shape.size(); ++i) {
          shape_str += std::to_string(shape[i]) + (i < shape.size() - 1 ? ", " : "");
        }
        utils::logmesg(lmp, "[GRACE] Feature '{}' [{}], size: {}, non_local: {}\n", key, shape_str,
                       size, non_local);
      }
    }
  } else
    error->all(FLERR, "[GRACE] metadata.yaml missing 'parallel_communication'");

  comm_forward = 0;
  comm_reverse = 0;
  for (const auto &[key, size] : feature_sizes) {
    comm_forward += size;
    if (!feature_is_local[key]) comm_reverse += size;
  }

  const int ntypes = atom->ntypes;
  element_type_mapping.resize(ntypes + 1);
  for (int i = 1; i <= ntypes; i++) {
    char *elemname = arg[2 + i];
    if (strcmp(elemname, "NULL") == 0) {
      element_type_mapping[i] = -1;
      map[i] = -1;
    } else {
      int mu = elements_to_index_map.count(elemname) ? elements_to_index_map.at(elemname) : -1;
      if (mu == -1) error->all(FLERR, "[GRACE] Element '{}' not supported", elemname);
      map[i] = mu;
      element_type_mapping[i] = mu;
    }
  }

  for (int i = 1; i <= ntypes; i++)
    for (int j = i; j <= ntypes; j++) scale[i][j] = 1.0;

  if (metadata_yaml["cutoff_matrix"]) {
    cutoff_matrix = metadata_yaml["cutoff_matrix"].as<std::vector<std::vector<double>>>();
    is_custom_cutoffs = true;
    cutoff_matrix_per_lammps_type.resize(ntypes + 1, std::vector<double>(ntypes + 1));
    for (int i = 1; i <= ntypes; i++)
      for (int j = 1; j <= ntypes; j++)
        cutoff_matrix_per_lammps_type[i][j] =
            cutoff_matrix[element_type_mapping[i]][element_type_mapping[j]];
  }

  if (!aceimpl->model->has_signature(forward_layer_1_name))
    error->all(FLERR, "Model missing forward_layer_1");
  if (!aceimpl->model->has_signature(backward_layer_2_name))
    error->all(FLERR, "Model missing backward_layer_2");
  if (!aceimpl->model->has_signature(backward_layer_1_name))
    error->all(FLERR, "Model missing backward_layer_1");

  has_atomic_mu_i_local =
      aceimpl->model->has_graph_input(DEFAULT_INPUT_PREFIX + "atomic_mu_i_local");
}

void PairGRACE2LayerChunk::init_style()
{
  if (atom->tag_enable == 0) error->all(FLERR, "Pair style grace requires atom IDs");
  if (force->newton_pair == 0) error->all(FLERR, "Pair style grace requires newton pair on");
  neighbor->add_request(this, NeighConst::REQ_FULL);

  if (atom->map_style == Atom::MAP_NONE) {
    atom->map_init();
    atom->map_set();
  }
}

double PairGRACE2LayerChunk::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");
  return is_custom_cutoffs ? cutoff_matrix_per_lammps_type[i][j] : cutoff;
}

void PairGRACE2LayerChunk::compute(int eflag, int vflag)
{
  total_timer.start_step();
  data_timer.start_step();
  comm_timer.start_step();
  model1_timer.start_step();
  model2_timer.start_step();
  model3_timer.start_step();

  current_step_real_atoms = 0;
  aceimpl->graph_recompiled = false;

  total_timer.start();
  data_timer.start();
  ev_init(eflag, vflag);

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int inum = list->inum;

  // Resize global buffers
  for (const auto &[key, size] : feature_sizes) {
    int n_all_max = nall;    // simplified, matching 2layer parallel's logic for resize
    features[key].resize(n_all_max * size);
    gradients[key].assign(n_all_max * size, 0.0);
  }

  if (aceimpl->global_to_chunk_map.size() < (size_t) nall)
    aceimpl->global_to_chunk_map.assign(nall, -1);

  bool do_energy_only = flag_compute_energy_only && !debug_no_energy_only_calc;

  data_timer.stop();

  // Phase 1: Chunked Forward L1
  for (int offset = 0; offset < inum; offset += chunksize) {
    data_timer.start();
    int current_size = std::min(chunksize, inum - offset);
    current_step_real_atoms += current_size;
    setup_chunk_graph(offset, current_size);
    data_timer.stop();

    model1_timer.start();
    run_forward_layer_1_chunk();
    model1_timer.stop();

    data_timer.start();
    for (int idx : aceimpl->chunk_to_global_map) aceimpl->global_to_chunk_map[idx] = -1;
    data_timer.stop();
  }

  // Phase 2: Monolithic Forward Comm
  comm_timer.start();
  comm->forward_comm(this);
  comm_timer.stop();

  // Phase 3: Chunked Backward L2
  for (int offset = 0; offset < inum; offset += chunksize) {
    data_timer.start();
    int current_size = std::min(chunksize, inum - offset);
    setup_chunk_graph(offset, current_size);
    data_timer.stop();

    model2_timer.start();
    run_backward_layer_2_chunk(eflag, vflag);
    model2_timer.stop();

    data_timer.start();
    for (int idx : aceimpl->chunk_to_global_map) aceimpl->global_to_chunk_map[idx] = -1;
    data_timer.stop();
  }

  // Phase 4: Monolithic Reverse Comm
  if (!do_energy_only) {
    comm_timer.start();
    comm->reverse_comm(this);
    comm_timer.stop();

    data_timer.start();
    // Zero ghost gradients to avoid double counting
    for (auto &[key, grad_vec] : gradients) {
      int size = feature_sizes[key];
      if (nall > nlocal) std::fill_n(&grad_vec[nlocal * size], (nall - nlocal) * size, 0.0);
    }
    data_timer.stop();

    // Phase 5: Chunked Backward L1
    for (int offset = 0; offset < inum; offset += chunksize) {
      data_timer.start();
      int current_size = std::min(chunksize, inum - offset);
      setup_chunk_graph(offset, current_size);
      data_timer.stop();

      model3_timer.start();
      run_backward_layer_1_chunk();
      model3_timer.stop();

      data_timer.start();
      for (int idx : aceimpl->chunk_to_global_map) aceimpl->global_to_chunk_map[idx] = -1;
      data_timer.stop();
    }
  }

  if (vflag_fdotr && !do_energy_only) virial_fdotr_compute();
  total_timer.stop();

  if (aceimpl->graph_recompiled) {
    total_timer.rollback();
    data_timer.rollback();
    comm_timer.rollback();
    model1_timer.rollback();
    model2_timer.rollback();
    model3_timer.rollback();
  } else {
    total_timer.commit();
    data_timer.commit();
    comm_timer.commit();
    model1_timer.commit();
    model2_timer.commit();
    model3_timer.commit();
    total_real_atoms_processed += current_step_real_atoms;
    total_compute_calls++;
  }
}

void PairGRACE2LayerChunk::setup_chunk_graph(int offset, int current_chunk_size)
{
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  double **x = atom->x;
  int *type = atom->type;

  aceimpl->chunk_to_global_map.clear();
  int node_counter = 0;

  // Map real atoms
  for (int ii = 0; ii < current_chunk_size; ii++) {
    int i = ilist[offset + ii];
    aceimpl->global_to_chunk_map[i] = node_counter++;
    aceimpl->chunk_to_global_map.push_back(i);
  }
  int n_real = node_counter;

  // Discover ghosts and bonds
  int n_bonds_real = 0;
  for (int ii = 0; ii < current_chunk_size; ii++) {
    int i = ilist[offset + ii];
    int type_i = type[i];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];
    double xi = x[i][0], yi = x[i][1], zi = x[i][2];

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj] & NEIGHMASK;
      double dx = xi - x[j][0], dy = yi - x[j][1], dz = zi - x[j][2];
      double rsq = dx * dx + dy * dy + dz * dz;
      double cutsq_ij = is_custom_cutoffs ? cutoff_matrix_per_lammps_type[type_i][type[j]] *
              cutoff_matrix_per_lammps_type[type_i][type[j]]
                                          : cutoff * cutoff;
      if (rsq < cutsq_ij) {
        if (aceimpl->global_to_chunk_map[j] == -1) {
          aceimpl->global_to_chunk_map[j] = node_counter++;
          aceimpl->chunk_to_global_map.push_back(j);
        }
        n_bonds_real++;
      }
    }
  }

  aceimpl->n_nodes_padded = aceimpl->atom_padding.update(node_counter);
  if (aceimpl->atom_padding.last_update_triggered_resize()) aceimpl->graph_recompiled = true;
  aceimpl->n_real_padded = aceimpl->real_atom_padding.update(n_real);
  if (aceimpl->real_atom_padding.last_update_triggered_resize()) aceimpl->graph_recompiled = true;
  aceimpl->n_bonds_padded = aceimpl->neighbor_padding.update(n_bonds_real);
  if (aceimpl->neighbor_padding.last_update_triggered_resize()) aceimpl->graph_recompiled = true;

  // Resize chunk-local buffers (conditional resize to avoid reallocation)
  aceimpl->mu_i.resize(aceimpl->n_bonds_padded);
  aceimpl->mu_j.resize(aceimpl->n_bonds_padded);
  aceimpl->ind_i.resize(aceimpl->n_bonds_padded);
  aceimpl->ind_j.resize(aceimpl->n_bonds_padded);
  aceimpl->bond_vector.resize(aceimpl->n_bonds_padded * 3);
  aceimpl->atomic_mu_i.resize(aceimpl->n_nodes_padded);
  aceimpl->atomic_mu_i_local.resize(aceimpl->n_real_padded);

  for (int k = 0; k < node_counter; k++)
    aceimpl->atomic_mu_i[k] = element_type_mapping[type[aceimpl->chunk_to_global_map[k]]];
  for (int k = node_counter; k < aceimpl->n_nodes_padded; k++)
    aceimpl->atomic_mu_i[k] = element_type_mapping[type[0]];

  for (int k = 0; k < n_real; k++)
    aceimpl->atomic_mu_i_local[k] = element_type_mapping[type[aceimpl->chunk_to_global_map[k]]];
  for (int k = n_real; k < aceimpl->n_real_padded; k++)
    aceimpl->atomic_mu_i_local[k] = element_type_mapping[type[0]];

  int bond_idx = 0;
  for (int ii = 0; ii < current_chunk_size; ii++) {
    int i = ilist[offset + ii];
    int type_i = type[i];
    int i_chunk = aceimpl->global_to_chunk_map[i];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];
    double xi = x[i][0], yi = x[i][1], zi = x[i][2];

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj] & NEIGHMASK;
      double dx = x[j][0] - xi, dy = x[j][1] - yi, dz = x[j][2] - zi;
      double rsq = dx * dx + dy * dy + dz * dz;
      double cutsq_ij = is_custom_cutoffs ? cutoff_matrix_per_lammps_type[type_i][type[j]] *
              cutoff_matrix_per_lammps_type[type_i][type[j]]
                                          : cutoff * cutoff;
      if (rsq < cutsq_ij) {
        int j_chunk = aceimpl->global_to_chunk_map[j];
        aceimpl->ind_i[bond_idx] = i_chunk;
        aceimpl->ind_j[bond_idx] = j_chunk;
        aceimpl->mu_i[bond_idx] = element_type_mapping[type_i];
        aceimpl->mu_j[bond_idx] = element_type_mapping[type[j]];
        aceimpl->bond_vector[3 * bond_idx + 0] = dx;
        aceimpl->bond_vector[3 * bond_idx + 1] = dy;
        aceimpl->bond_vector[3 * bond_idx + 2] = dz;
        bond_idx++;
      }
    }
  }

  // Zero-fill padding tail
  for (int k = bond_idx; k < aceimpl->n_bonds_padded; k++) {
    aceimpl->ind_i[k] = aceimpl->n_nodes_padded - 1;
    aceimpl->ind_j[k] = aceimpl->n_nodes_padded - 1;
    aceimpl->mu_i[k] = 0;
    aceimpl->mu_j[k] = 0;
    aceimpl->bond_vector[3 * k + 0] = 1e6;
    aceimpl->bond_vector[3 * k + 1] = 1e6;
    aceimpl->bond_vector[3 * k + 2] = 1e6;
  }
}

void PairGRACE2LayerChunk::run_forward_layer_1_chunk()
{
  auto sig = aceimpl->model->signatures.at(forward_layer_1_name);
  std::vector<std::tuple<std::string, cppflow::tensor>> inputs;
  auto add_in = [&](const std::string &key, const cppflow::tensor &t) {
    if (sig.inputs.count(key)) inputs.emplace_back(sig.inputs.at(key).name, t);
  };

  int n_real =
      aceimpl->chunk_to_global_map
          .size();    // Only counting mapped, correct? No, real is first current_chunk_size
  add_in("atomic_mu_i",
         GRACE::get_or_create_tensor(aceimpl->fwd_l1_tensors, "atomic_mu_i", aceimpl->atomic_mu_i,
                                     {(int64_t) aceimpl->n_nodes_padded},
                                     aceimpl->graph_recompiled));
  add_in("atomic_mu_i_local",
         GRACE::get_or_create_tensor(aceimpl->fwd_l1_tensors, "atomic_mu_i_local",
                                     aceimpl->atomic_mu_i_local, {(int64_t) aceimpl->n_real_padded},
                                     aceimpl->graph_recompiled));
  add_in("ind_i",
         GRACE::get_or_create_tensor(aceimpl->fwd_l1_tensors, "ind_i", aceimpl->ind_i,
                                     {(int64_t) aceimpl->n_bonds_padded},
                                     aceimpl->graph_recompiled));
  add_in("ind_j",
         GRACE::get_or_create_tensor(aceimpl->fwd_l1_tensors, "ind_j", aceimpl->ind_j,
                                     {(int64_t) aceimpl->n_bonds_padded},
                                     aceimpl->graph_recompiled));
  add_in("mu_i",
         GRACE::get_or_create_tensor(aceimpl->fwd_l1_tensors, "mu_i", aceimpl->mu_i,
                                     {(int64_t) aceimpl->n_bonds_padded},
                                     aceimpl->graph_recompiled));
  add_in("mu_j",
         GRACE::get_or_create_tensor(aceimpl->fwd_l1_tensors, "mu_j", aceimpl->mu_j,
                                     {(int64_t) aceimpl->n_bonds_padded},
                                     aceimpl->graph_recompiled));
  add_in("bond_vector",
         GRACE::get_or_create_tensor(aceimpl->fwd_l1_tensors, "bond_vector", aceimpl->bond_vector,
                                     {(int64_t) aceimpl->n_bonds_padded, 3},
                                     aceimpl->graph_recompiled));
  // How many real atoms in this chunk? It's the current_chunk_size passed to setup.
  // I'll assume n_real_padded.update(current_chunk_size) was done.
  // Let's use the real count from setup.
  int n_real_actual = aceimpl->atomic_mu_i_local.size() - aceimpl->real_atom_padding.n_fake;
  add_in("batch_tot_nat_real",
         GRACE::get_or_create_tensor(aceimpl->fwd_l1_tensors, "batch_tot_nat_real",
                                     std::vector<int32_t>{n_real_actual}, {},
                                     aceimpl->graph_recompiled));

  std::vector<std::string> out_names;
  std::vector<std::string> ordered_keys;
  for (const auto &[key, shape] : feature_shapes) {
    out_names.push_back(sig.outputs.at(key).name);
    ordered_keys.push_back(key);
  }

  auto outputs = aceimpl->model->operator()(inputs, out_names);
  for (size_t i = 0; i < outputs.size(); i++) {
    std::string key = ordered_keys[i];
    int size = feature_sizes[key];
    const double *data = static_cast<double *>(TF_TensorData(outputs[i].get_tensor().get()));
    // Scatter only real atoms
    for (int k = 0; k < n_real_actual; k++) {
      int g_idx = aceimpl->chunk_to_global_map[k];
      std::copy_n(&data[k * size], size, &features[key][g_idx * size]);
    }
  }
}

void PairGRACE2LayerChunk::run_backward_layer_1_chunk()
{
  auto sig = aceimpl->model->signatures.at(backward_layer_1_name);
  std::vector<std::tuple<std::string, cppflow::tensor>> inputs;
  auto add_in = [&](const std::string &key, const cppflow::tensor &t) {
    if (sig.inputs.count(key)) inputs.emplace_back(sig.inputs.at(key).name, t);
  };

  int n_real_actual =
      aceimpl->real_atom_padding.current_padded_size - aceimpl->real_atom_padding.n_fake;

  add_in("atomic_mu_i",
         GRACE::get_or_create_tensor(aceimpl->bwd_l1_tensors, "atomic_mu_i", aceimpl->atomic_mu_i,
                                     {(int64_t) aceimpl->n_nodes_padded},
                                     aceimpl->graph_recompiled));
  add_in("atomic_mu_i_local",
         GRACE::get_or_create_tensor(aceimpl->bwd_l1_tensors, "atomic_mu_i_local",
                                     aceimpl->atomic_mu_i_local, {(int64_t) aceimpl->n_real_padded},
                                     aceimpl->graph_recompiled));
  add_in("ind_i",
         GRACE::get_or_create_tensor(aceimpl->bwd_l1_tensors, "ind_i", aceimpl->ind_i,
                                     {(int64_t) aceimpl->n_bonds_padded},
                                     aceimpl->graph_recompiled));
  add_in("ind_j",
         GRACE::get_or_create_tensor(aceimpl->bwd_l1_tensors, "ind_j", aceimpl->ind_j,
                                     {(int64_t) aceimpl->n_bonds_padded},
                                     aceimpl->graph_recompiled));
  add_in("mu_i",
         GRACE::get_or_create_tensor(aceimpl->bwd_l1_tensors, "mu_i", aceimpl->mu_i,
                                     {(int64_t) aceimpl->n_bonds_padded},
                                     aceimpl->graph_recompiled));
  add_in("mu_j",
         GRACE::get_or_create_tensor(aceimpl->bwd_l1_tensors, "mu_j", aceimpl->mu_j,
                                     {(int64_t) aceimpl->n_bonds_padded},
                                     aceimpl->graph_recompiled));
  add_in("bond_vector",
         GRACE::get_or_create_tensor(aceimpl->bwd_l1_tensors, "bond_vector", aceimpl->bond_vector,
                                     {(int64_t) aceimpl->n_bonds_padded, 3},
                                     aceimpl->graph_recompiled));
  add_in("batch_tot_nat_real",
         GRACE::get_or_create_tensor(aceimpl->bwd_l1_tensors, "batch_tot_nat_real",
                                     std::vector<int32_t>{n_real_actual}, {},
                                     aceimpl->graph_recompiled));

  // GATHER gradients
  for (const auto &[key, shape] : feature_shapes) {
    int size = feature_sizes[key];
    std::vector<double> chunk_grad(aceimpl->n_real_padded * size, 0.0);
    for (int k = 0; k < n_real_actual; k++) {
      int g_idx = aceimpl->chunk_to_global_map[k];
      std::copy_n(&gradients[key][g_idx * size], size, &chunk_grad[k * size]);
    }
    std::vector<int64_t> t_shape = {(int64_t) aceimpl->n_real_padded};
    t_shape.insert(t_shape.end(), shape.begin(), shape.end());
    add_in("grad_" + key,
           GRACE::get_or_create_tensor(aceimpl->bwd_l1_tensors, "grad_" + key, chunk_grad, t_shape,
                                       aceimpl->graph_recompiled));
  }

  std::vector<std::string> out_names = {sig.outputs.at(GRAD_BOND_KEY).name};
  auto outputs = aceimpl->model->operator()(inputs, out_names);
  const double *gbv_data = static_cast<double *>(TF_TensorData(outputs[0].get_tensor().get()));

  // Apply forces (L1 part)
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  int bond_idx = 0;
  for (int k = 0; k < aceimpl->n_bonds_padded; k++) {
    int i_chunk = aceimpl->ind_i[k];
    int j_chunk = aceimpl->ind_j[k];
    if (i_chunk >= n_real_actual || j_chunk >= (int) aceimpl->chunk_to_global_map.size()) continue;

    int i = aceimpl->chunk_to_global_map[i_chunk];
    int j = aceimpl->chunk_to_global_map[j_chunk];
    if (i >= nlocal && j >= nlocal) continue;    // Should not happen for active bonds in chunk

    double sc = scale[type[i]][type[i]];
    double fx = sc * gbv_data[3 * k + 0];
    double fy = sc * gbv_data[3 * k + 1];
    double fz = sc * gbv_data[3 * k + 2];

    f[i][0] += fx;
    f[i][1] += fy;
    f[i][2] += fz;
    f[j][0] -= fx;
    f[j][1] -= fy;
    f[j][2] -= fz;

    if (evflag) {
      double dx = -aceimpl->bond_vector[3 * k + 0];
      double dy = -aceimpl->bond_vector[3 * k + 1];
      double dz = -aceimpl->bond_vector[3 * k + 2];
      ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, fx, fy, fz, dx, dy, dz);
    }
  }
}

// MPI Comm methods
int PairGRACE2LayerChunk::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int m = 0;
  std::vector<std::pair<double *, int>> active;
  for (const auto &[key, size] : feature_sizes)
    if (!feature_is_local[key]) active.push_back({features[key].data(), size});
  for (int i = 0; i < n; i++) {
    int j = list[i];
    for (auto &feat : active) {
      std::copy_n(&feat.first[j * feat.second], feat.second, &buf[m]);
      m += feat.second;
    }
  }
  return m;
}

void PairGRACE2LayerChunk::unpack_forward_comm(int n, int first, double *buf)
{
  int m = 0;
  std::vector<std::pair<double *, int>> active;
  for (const auto &[key, size] : feature_sizes)
    if (!feature_is_local[key]) active.push_back({features[key].data(), size});
  for (int i = first; i < first + n; i++) {
    for (auto &feat : active) {
      std::copy_n(&buf[m], feat.second, &feat.first[i * feat.second]);
      m += feat.second;
    }
  }
}

int PairGRACE2LayerChunk::pack_reverse_comm(int n, int first, double *buf)
{
  int m = 0;
  std::vector<std::pair<double *, int>> active;
  for (const auto &[key, size] : feature_sizes)
    if (!feature_is_local[key]) active.push_back({gradients[key].data(), size});
  for (int i = first; i < first + n; i++) {
    for (auto &feat : active) {
      std::copy_n(&feat.first[i * feat.second], feat.second, &buf[m]);
      m += feat.second;
    }
  }
  return m;
}

void PairGRACE2LayerChunk::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m = 0;
  std::vector<std::pair<double *, int>> active;
  for (const auto &[key, size] : feature_sizes)
    if (!feature_is_local[key]) active.push_back({gradients[key].data(), size});
  for (int i = 0; i < n; i++) {
    int j = list[i];
    for (auto &feat : active) {
      for (int k = 0; k < feat.second; k++) feat.first[j * feat.second + k] += buf[m++];
    }
  }
}

void *PairGRACE2LayerChunk::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str, "compute_energy_only") == 0) return (void *) &flag_compute_energy_only;
  dim = 2;
  if (strcmp(str, "scale") == 0) return (void *) scale;
  return nullptr;
}

// Dummy print_tensors to match header (optional debug info)
void PairGRACE2LayerChunk::print_tensors(
    const std::string &name, const std::vector<std::tuple<std::string, cppflow::tensor>> &tensors,
    const std::string &type_prefix)
{
}

// Missing part in run_backward_layer_2_chunk (Scatter GBV) added inside implementation hereafter.
// Re-writing run_backward_layer_2_chunk completely to be safe.

void PairGRACE2LayerChunk::run_backward_layer_2_chunk(int eflag, int vflag)
{
  auto sig = aceimpl->model->signatures.at(backward_layer_2_name);
  std::vector<std::tuple<std::string, cppflow::tensor>> inputs;
  auto add_in = [&](const std::string &key, const cppflow::tensor &t) {
    if (sig.inputs.count(key)) inputs.emplace_back(sig.inputs.at(key).name, t);
  };

  int n_all = aceimpl->chunk_to_global_map.size();
  int n_real_actual =
      aceimpl->real_atom_padding.current_padded_size - aceimpl->real_atom_padding.n_fake;

  add_in("atomic_mu_i",
         GRACE::get_or_create_tensor(aceimpl->bwd_l2_tensors, "atomic_mu_i", aceimpl->atomic_mu_i,
                                     {(int64_t) aceimpl->n_nodes_padded},
                                     aceimpl->graph_recompiled));
  add_in("atomic_mu_i_local",
         GRACE::get_or_create_tensor(aceimpl->bwd_l2_tensors, "atomic_mu_i_local",
                                     aceimpl->atomic_mu_i_local, {(int64_t) aceimpl->n_real_padded},
                                     aceimpl->graph_recompiled));
  add_in("ind_i",
         GRACE::get_or_create_tensor(aceimpl->bwd_l2_tensors, "ind_i", aceimpl->ind_i,
                                     {(int64_t) aceimpl->n_bonds_padded},
                                     aceimpl->graph_recompiled));
  add_in("ind_j",
         GRACE::get_or_create_tensor(aceimpl->bwd_l2_tensors, "ind_j", aceimpl->ind_j,
                                     {(int64_t) aceimpl->n_bonds_padded},
                                     aceimpl->graph_recompiled));
  add_in("mu_i",
         GRACE::get_or_create_tensor(aceimpl->bwd_l2_tensors, "mu_i", aceimpl->mu_i,
                                     {(int64_t) aceimpl->n_bonds_padded},
                                     aceimpl->graph_recompiled));
  add_in("mu_j",
         GRACE::get_or_create_tensor(aceimpl->bwd_l2_tensors, "mu_j", aceimpl->mu_j,
                                     {(int64_t) aceimpl->n_bonds_padded},
                                     aceimpl->graph_recompiled));
  add_in("bond_vector",
         GRACE::get_or_create_tensor(aceimpl->bwd_l2_tensors, "bond_vector", aceimpl->bond_vector,
                                     {(int64_t) aceimpl->n_bonds_padded, 3},
                                     aceimpl->graph_recompiled));
  add_in("batch_tot_nat_real",
         GRACE::get_or_create_tensor(aceimpl->bwd_l2_tensors, "batch_tot_nat_real",
                                     std::vector<int32_t>{n_real_actual}, {},
                                     aceimpl->graph_recompiled));

  for (const auto &[key, shape] : feature_shapes) {
    int size = feature_sizes[key];
    int n_padded = feature_is_local[key] ? aceimpl->n_real_padded : aceimpl->n_nodes_padded;
    std::vector<double> chunk_feat(n_padded * size, 0.0);
    int n_actual = feature_is_local[key] ? n_real_actual : n_all;
    for (int k = 0; k < n_actual; k++) {
      int g_idx = aceimpl->chunk_to_global_map[k];
      std::copy_n(&features[key][g_idx * size], size, &chunk_feat[k * size]);
    }
    std::vector<int64_t> t_shape = {(int64_t) n_padded};
    t_shape.insert(t_shape.end(), shape.begin(), shape.end());
    add_in(key,
           GRACE::get_or_create_tensor(aceimpl->bwd_l2_tensors, key, chunk_feat, t_shape,
                                       aceimpl->graph_recompiled));
  }

  std::vector<std::string> out_names;
  std::vector<std::string> ordered_keys;
  out_names.push_back(sig.outputs.at(ENERGY_KEY).name);
  ordered_keys.push_back(ENERGY_KEY);
  bool do_energy_only = flag_compute_energy_only && !debug_no_energy_only_calc;
  if (!do_energy_only) {
    for (const auto &[key, shape] : feature_shapes) {
      out_names.push_back(sig.outputs.at("grad_" + key).name);
      ordered_keys.push_back("grad_" + key);
    }
    out_names.push_back(sig.outputs.at(GRAD_BOND_KEY).name);
    ordered_keys.push_back(GRAD_BOND_KEY);
  }

  auto outputs = aceimpl->model->operator()(inputs, out_names);
  int out_idx = 0;
  const double *e_data =
      static_cast<double *>(TF_TensorData(outputs[out_idx++].get_tensor().get()));
  if (eflag_either) {
    for (int k = 0; k < n_real_actual; k++)
      ev_tally_full(aceimpl->chunk_to_global_map[k], 2.0 * e_data[k], 0.0, 0.0, 0.0, 0.0, 0.0);
  }

  if (!do_energy_only) {
    for (const auto &[key, shape] : feature_shapes) {
      int size = feature_sizes[key];
      int n_actual = feature_is_local[key] ? n_real_actual : n_all;
      const double *g_data =
          static_cast<double *>(TF_TensorData(outputs[out_idx++].get_tensor().get()));
      for (int k = 0; k < n_actual; k++) {
        int g_idx = aceimpl->chunk_to_global_map[k];
        for (int s = 0; s < size; s++) gradients[key][g_idx * size + s] += g_data[k * size + s];
      }
    }
    const double *gbv_data =
        static_cast<double *>(TF_TensorData(outputs[out_idx++].get_tensor().get()));
    double **f = atom->f;
    int *type = atom->type;
    int nlocal = atom->nlocal;
    int newton_pair = force->newton_pair;
    for (int k = 0; k < aceimpl->n_bonds_padded; k++) {
      int i_chunk = aceimpl->ind_i[k];
      int j_chunk = aceimpl->ind_j[k];
      if (i_chunk >= n_real_actual || j_chunk >= n_all) continue;
      int i = aceimpl->chunk_to_global_map[i_chunk];
      int j = aceimpl->chunk_to_global_map[j_chunk];
      double sc = scale[type[i]][type[i]];
      double fx = sc * gbv_data[3 * k + 0], fy = sc * gbv_data[3 * k + 1],
             fz = sc * gbv_data[3 * k + 2];
      f[i][0] += fx;
      f[i][1] += fy;
      f[i][2] += fz;
      f[j][0] -= fx;
      f[j][1] -= fy;
      f[j][2] -= fz;
      if (evflag) {
        double dx = -aceimpl->bond_vector[3 * k + 0], dy = -aceimpl->bond_vector[3 * k + 1],
               dz = -aceimpl->bond_vector[3 * k + 2];
        ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, fx, fy, fz, dx, dy, dz);
      }
    }
  }
}

#endif
