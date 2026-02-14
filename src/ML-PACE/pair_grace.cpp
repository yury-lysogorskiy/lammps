//
// Created by Yury Lysogorskiy on 01.12.23.
//
#ifndef NO_GRACE_TF
// #define GRACE_PRINT_DEBUG

#include "pair_grace.h"

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

struct GRACEImpl {
  GRACEImpl() : model(nullptr) {}

  ~GRACEImpl() { delete model; }

  GRACE::GracePaddingDimension atom_padding;
  GRACE::GracePaddingDimension neighbor_padding;

  cppflow::model *model;
  bool graph_recompiled = false;
  std::map<std::string, cppflow::TensorInfo> compute_inputs_sig;
  std::map<std::string, cppflow::TensorInfo> compute_energy_only_inputs_sig;
};

}    // namespace LAMMPS_NS

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */
PairGRACE::PairGRACE(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  graceimpl = new GRACEImpl;

  scale = nullptr;

  chunksize = 4096;

  total_timer.init();
  data_timer.init();

  no_virial_fdotr_compute = 1;
  flag_compute_energy_only = 0;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */
PairGRACE::~PairGRACE()
{
  if (copymode) return;

  if (comm->me == 0 && total_real_atoms_processed > 0) {
    double total = total_timer.as_microseconds();
    double data = data_timer.as_microseconds();
    double model = model_timer.as_microseconds();

    auto per_atom = [&](double t) { return t / total_real_atoms_processed; };
    auto pct = [&](double t) { return (total > 0) ? (t / total * 100.0) : 0.0; };

    utils::logmesg(lmp, "[GRACE-PERF] Total real atoms processed: {:.0f}\n", total_real_atoms_processed);
    utils::logmesg(lmp, "[GRACE-PERF] Total valid compute calls: {}\n", total_compute_calls);
    utils::logmesg(lmp, "[GRACE-PERF] Average atoms per step: {:.1f}\n",
                   total_real_atoms_processed / total_compute_calls);
    utils::logmesg(lmp, "[GRACE-PERF] Performance (us/atom) [%]: Total: {:.1f}, Data: {:.1f} ({:.1f}%), Model: {:.1f} ({:.1f}%)\n",
                   per_atom(total), per_atom(data), pct(data), per_atom(model), pct(model));
  }

  delete graceimpl;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(scale);
  }
}

/* ---------------------------------------------------------------------- */
void PairGRACE::allocate()
{
  allocated = 1;
  int n = atom->ntypes + 1;

  memory->create(setflag, n, n, "pair:setflag");
  memory->create(cutsq, n, n, "pair:cutsq");
  memory->create(scale, n, n, "pair:scale");
  map = new int[n];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */
void PairGRACE::settings(int narg, char **arg)
{
  if (narg > 3) utils::missing_cmd_args(FLERR, "pair_style grace", error);

  // ACE potentials are parameterized in metal units
  if (strcmp("metal", update->unit_style) != 0)
    error->all(FLERR, "GRACE potentials require 'metal' units");

  auto tf_version = TF_Version();
  if (comm->me == 0) utils::logmesg(lmp, "[GRACE] TF version: {}\n", tf_version);

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
    } else if (strcmp(arg[iarg], "pair_forces") == 0) {
      pair_forces = true;
      iarg += 1;
      if (comm->me == 0) utils::logmesg(lmp, "[GRACE] Pair forces are ON \n");
    } else if (strcmp(arg[iarg], "max_number_of_reduction") == 0) {
      max_number_of_reduction = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
      if (comm->me == 0)
        utils::logmesg(lmp,
                       "[GRACE] Maximum number of recompilation during padding reduction: {}\n",
                       max_number_of_reduction);
    } else if (strcmp(arg[iarg], "reduce_padding") == 0) {
      reducing_neigh_padding_fraction = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
      if (comm->me == 0)
        utils::logmesg(lmp, "[GRACE] Reducing padding fraction: {}\n",
                       reducing_neigh_padding_fraction);
    } else if (strcmp(arg[iarg], "deny_energy_only_calc") == 0) {
      deny_energy_only_calc = true;
      iarg += 1;
    } else
      error->all(FLERR, "[GRACE] Unknown pair_style grace keyword: {}", arg[iarg]);
  }

  do_padding = (neigh_padding_fraction > 0);
  if (do_padding)
    if (comm->me == 0)
      utils::logmesg(lmp,
                     "[GRACE] Neighbour padding is ON, padding fraction: {}, max padding fraction "
                     "before reduction: {}, max number of reduction(s): {}\n",
                     neigh_padding_fraction, reducing_neigh_padding_fraction,
                     max_number_of_reduction);

  // Apply to padding helpers
  graceimpl->atom_padding.enabled = do_padding;
  graceimpl->atom_padding.padding_fraction = neigh_padding_fraction;
  graceimpl->atom_padding.reduction_threshold_fraction = reducing_neigh_padding_fraction;
  graceimpl->atom_padding.max_reductions = max_number_of_reduction;
  graceimpl->atom_padding.verbose = pad_verbose;

  graceimpl->neighbor_padding.enabled = do_padding;
  graceimpl->neighbor_padding.padding_fraction = neigh_padding_fraction;
  graceimpl->neighbor_padding.reduction_threshold_fraction = reducing_neigh_padding_fraction;
  graceimpl->neighbor_padding.max_reductions = max_number_of_reduction;
  graceimpl->neighbor_padding.verbose = pad_verbose;

  if (!pair_forces && comm->nprocs > 1) {
    pair_forces = true;
    if (comm->me == 0)
      utils::logmesg(lmp,
                     "[GRACE] ENFORCE pair-force mode to ON, because number of processes {} is "
                     "more than one.\n",
                     comm->nprocs);
  }

  if (pair_forces) {
    // mark as available centroid stress flag
    centroidstressflag = CENTROID_AVAIL;
  } else {
    centroidstressflag = CENTROID_NOTAVAIL;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGRACE::coeff(int narg, char **arg)
{

  if (!allocated) allocate();

  map_element2type(narg - 3, arg + 3);
  auto potential_path = std::string(arg[2]);

  //load potential file
  delete graceimpl->model;
  //load potential file
  if (comm->me == 0) utils::logmesg(lmp, "[GRACE] Loading {}\n", potential_path);

  // load cppflow model
  const std::vector<uint8_t> config_bytes = {
      0x32, 0x05,          // Field 6 (GPUOptions), Length 5
      0x82, 0x01, 0x02,    // Field 16 (Experimental), Length 2
      0x18, 0x00           // Field 3 (TF32 Enabled), Value 0 (False)
  };

  graceimpl->model = new cppflow::model(potential_path, config_bytes);
#ifdef GRACE_PRINT_DEBUG
  graceimpl->model->print_signatures();
#endif
  // read elements from metadata.yaml
  YAML_PACE::Node metadata_yaml = YAML_PACE::LoadFile(potential_path + "/metadata.yaml");
  auto elements_yaml = metadata_yaml["chemical_symbols"];
  elements_name = elements_yaml.as<std::vector<std::string>>();
  nelements = (int) elements_name.size();
  for (int mu = 0; mu < nelements; mu++) { elements_to_index_map[elements_name.at(mu)] = mu; }
  cutoff = metadata_yaml["cutoff"].as<double>();

  if (metadata_yaml["cutoff_matrix"]) {
    cutoff_matrix = metadata_yaml["cutoff_matrix"].as<vector<vector<double>>>();
    // assert square size of matrix
    if (cutoff_matrix.size() != nelements)
      error->all(FLERR,
                 "[GRACE] cutoff_matrix is provided, but it's size ({}) is not equal to number of "
                 "elements ({})\n",
                 cutoff_matrix.size(), nelements);

    for (const auto &v : cutoff_matrix)
      if (v.size() != nelements)
        error->all(FLERR,
                   "[GRACE] cutoff_matrix is provided, but it's row size ({}) is not equal to "
                   "number of elements ({})\n",
                   v.size(), nelements);

    if (comm->me == 0) utils::logmesg(lmp, "[GRACE] Custom cutoff matrix is loaded\n");
    is_custom_cutoffs = true;
  }
  if (comm->me == 0) { utils::logmesg(lmp, "[GRACE] Model loaded\n"); }

  // read args that map atom types to PACE elements
  // map[i] = which element the Ith atom type is, -1 if not mapped
  // map[0] is not used

  const int ntypes = atom->ntypes;
  element_type_mapping.resize(ntypes + 1);
  // elements to species-type map
  for (int i = 1; i <= ntypes; i++) {
    char *elemname = arg[2 + i];
    if (strcmp(elemname, "NULL") == 0) {
      // species_type=-1 value will not reach ACE Evaluator::compute_atom,
      // but if it will ,then error will be thrown there
      element_type_mapping[i] = -1;
      map[i] = -1;
      if (comm->me == 0) utils::logmesg(lmp, "[GRACE] Skipping LAMMPS atom type #{}(NULL)\n", i);
    } else {
      int atomic_number = PACE::AtomicNumberByName(elemname);
      if (atomic_number == -1) error->all(FLERR, "[GRACE] '{}' is not a valid element\n", elemname);
      int mu = elements_to_index_map.at(elemname);
      if (mu != -1) {
        if (comm->me == 0)
          utils::logmesg(lmp, "[GRACE] Mapping LAMMPS atom type #{}({}) -> ACE species type #{}\n",
                         i, elemname, mu);
        map[i] = mu;
        // set up LAMMPS atom type to ACE species  mapping for ace evaluator
        element_type_mapping[i] = mu;
      } else {
        error->all(FLERR, "[GRACE] Element {} is not supported by ACE-potential from file {}",
                   elemname, potential_path);
      }
    }
  }

  // initialize scale factor
  for (int i = 1; i <= ntypes; i++) {
    for (int j = i; j <= ntypes; j++) scale[i][j] = 1.0;
  }

  if (is_custom_cutoffs) {
    // matrix of size [ntypes+1][ntypes+1]
    cutoff_matrix_per_lammps_type.resize(ntypes + 1, vector<double>(ntypes + 1));
    double min_cutoff = 1e99, max_cutoff = 0;
    for (int i = 1; i <= ntypes; i++) {
      for (int j = 1; j <= ntypes; j++) {
        auto val = cutoff_matrix[element_type_mapping[i]][element_type_mapping[j]];
        cutoff_matrix_per_lammps_type[i][j] = val;
        if (val < min_cutoff) min_cutoff = val;
        if (val > max_cutoff) max_cutoff = val;
      }
    }

    if (comm->me == 0)
      utils::logmesg(lmp, "[GRACE] Custom cutoffs: min={}, max={}\n", min_cutoff, max_cutoff);
  }

  if (graceimpl->model->has_signature("compute")) {
    compute_function_name = "compute";
  } else if (graceimpl->model->has_signature("serving_default")) {
    compute_function_name = "serving_default";
  }
  graceimpl->compute_inputs_sig = graceimpl->model->signatures.at(compute_function_name).inputs;

  // check for compute_energy_only function
  if (graceimpl->model->has_signature(COMPUTE_ENERGY_ONLY_KEY)) {
    compute_energy_only_function_name = COMPUTE_ENERGY_ONLY_KEY;
    has_compute_energy_only = true;
    if (comm->me == 0) utils::logmesg(lmp, "[GRACE] Compute energy only function is available\n");
    graceimpl->compute_energy_only_inputs_sig =
        graceimpl->model->signatures.at(COMPUTE_ENERGY_ONLY_KEY).inputs;
  } else {
  }

  this->DEFAULT_INPUT_PREFIX = this->compute_function_name + "_";

  //
  has_map_atoms_to_structure_op =
      graceimpl->model->has_graph_input(DEFAULT_INPUT_PREFIX + "map_atoms_to_structure");

  has_nstruct_total_op = graceimpl->model->has_graph_input(DEFAULT_INPUT_PREFIX + "n_struct_total");
  has_mu_i_op = graceimpl->model->has_graph_input(DEFAULT_INPUT_PREFIX + "mu_i");
  has_batch_tot_nat = graceimpl->model->has_graph_input(DEFAULT_INPUT_PREFIX + "batch_tot_nat");

  has_atomic_mu_i_local =
      graceimpl->model->has_graph_input(DEFAULT_INPUT_PREFIX + "atomic_mu_i_local");
  if (has_atomic_mu_i_local) {
    if (comm->me == 0) utils::logmesg(lmp, "[GRACE-DEBUG] atomic_mu_i_local is available\n");
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairGRACE::init_style()
{
  if (atom->tag_enable == 0) error->all(FLERR, "Pair style grace requires atom IDs");
  if (force->newton_pair == 0) error->all(FLERR, "Pair style grace requires newton pair on");

  // request a full neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL);

  // request atom map (maybe?)
  if (atom->map_style == Atom::MAP_NONE) {
    atom->map_init();
    atom->map_set();
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGRACE::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");
  //cutoff from the basis set's radial functions settings
  scale[j][i] = scale[i][j];

  if (is_custom_cutoffs) return cutoff_matrix_per_lammps_type[i][j];
  return cutoff;
}

/* ----------------------------------------------------------------------
    extract method for extracting value of scale variable
 ---------------------------------------------------------------------- */
void *PairGRACE::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str, "compute_energy_only") == 0) return (void *) &flag_compute_energy_only;

  dim = 2;
  if (strcmp(str, "scale") == 0) return (void *) scale;

  return nullptr;
}

/* ---------------------------------------------------------------------- */
/**
 * signature_def['serving_default']:
  The given SavedModel SignatureDef contains the following input(s):
    +inputs['atomic_mu_i'] tensor_info:
        dtype: DT_INT32
        shape: (-1)
        name: serving_default_atomic_mu_i:0
    +inputs['batch_tot_nat'] tensor_info:
        dtype: DT_INT32
        shape: ()
        name: serving_default_batch_tot_nat:0
    +inputs['batch_tot_nat_real'] tensor_info:
        dtype: DT_INT32
        shape: ()
        name: serving_default_batch_tot_nat_real:0
    +inputs['bond_vector'] tensor_info:
        dtype: DT_DOUBLE
        shape: (-1, 3)
        name: serving_default_bond_vector:0
    +inputs['ind_i'] tensor_info:
        dtype: DT_INT32
        shape: (-1)
        name: serving_default_ind_i:0
    +inputs['ind_j'] tensor_info:
        dtype: DT_INT32
        shape: (-1)
        name: serving_default_ind_j:0
    +inputs['map_atoms_to_structure'] tensor_info:
        dtype: DT_INT32
        shape: (-1)
        name: serving_default_map_atoms_to_structure:0
    +inputs['mu_j'] tensor_info:
        dtype: DT_INT32
        shape: (-1)
        name: serving_default_mu_j:0
    +inputs['n_struct_total'] tensor_info:
        dtype: DT_INT32
        shape: ()
        name: serving_default_n_struct_total:0


  The given SavedModel SignatureDef contains the following output(s):
    outputs['atomic_energy'] tensor_info:
        dtype: DT_DOUBLE
        shape: (-1, 1)
        name: StatefulPartitionedCall:0
    outputs['total_energy'] tensor_info:
        dtype: DT_DOUBLE
        shape: (1, 1)
        name: StatefulPartitionedCall:1
    outputs['total_f'] tensor_info:
        dtype: DT_DOUBLE
        shape: (-1, 3)
        name: StatefulPartitionedCall:2
    outputs['virial'] tensor_info:
        dtype: DT_DOUBLE
        shape: (6)
        name: StatefulPartitionedCall:3
    outputs['z_pair_f'] tensor_info:
        dtype: DT_DOUBLE
        shape: (-1, 3)
        name: StatefulPartitionedCall:4

  Method name is: tensorflow/serving/predict
 */
void PairGRACE::compute(int eflag, int vflag)
{
  total_timer.start_step();
  data_timer.start_step();
  model_timer.start_step();

  total_timer.start();
  current_step_real_atoms = 0;
  graceimpl->graph_recompiled = false;
  int i, j, ii, jj, inum, jnum;
  double delx, dely, delz, evdwl;
  double fij[3];
  int *ilist, *jlist, *numneigh, **firstneigh;

  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  //    tagint *tag = atom->tag;

  // number of atoms in cell
  int nlocal = atom->nlocal;
  int n_fake_atoms;    // no fake atoms needed. fake bonds will be 1e6
  int n_real_neighbours;
  int n_fake_neighbours;

  int newton_pair = force->newton_pair;

  // inum: length of the neighborlists list
  inum = list->inum;

  // ilist: list of "i" atoms for which neighbor lists exist
  ilist = list->ilist;

  //numneigh: the length of each these neigbor list
  numneigh = list->numneigh;

  // the pointer to the list of neighbors of "i"
  firstneigh = list->firstneigh;

  bool do_energy_only = flag_compute_energy_only & !deny_energy_only_calc;

  auto compute_inputs_sig = graceimpl->compute_inputs_sig;
  if (do_energy_only) {
    if (has_compute_energy_only) {
      compute_inputs_sig = graceimpl->compute_energy_only_inputs_sig;
    } else {
      do_energy_only = false;
      if (!warning_compute_energy_only_not_avail_shown) {
        utils::logmesg(
            lmp,
            std::string(
                "[GRACE:WARNING] Compute energy only function is not available, but requested. ") +
                "Full compute function will be used. This message is shown only once\n");
        warning_compute_energy_only_not_avail_shown = true;
      }
    }
  }

  // utils::logmesg(lmp, "[GRACE-debug] input_prefix={} \n",input_prefix);

  data_timer.start();
  std::vector<std::tuple<std::string, cppflow::tensor>> inputs;

  tot_atoms = graceimpl->atom_padding.update(nlocal);
  if (graceimpl->atom_padding.last_update_triggered_resize()) graceimpl->graph_recompiled = true;
  if (pad_verbose && graceimpl->atom_padding.last_update_triggered_resize()) {
    utils::logmesg(lmp,
                   "[GRACE] Atoms padding: new num. of atoms = {} (incl. {:.3f}% fake atoms)\n",
                   tot_atoms, 100. * (double) graceimpl->atom_padding.n_fake / tot_atoms);
  }
  current_step_real_atoms = nlocal;

  // atomic_mu_i: per-atom species type + padding with type[0]
  std::vector<int32_t> atomic_mu_i_vector(tot_atoms, element_type_mapping[type[0]]);
  for (i = 0; i < nlocal; ++i) atomic_mu_i_vector[i] = element_type_mapping[type[i]];

  inputs.emplace_back(compute_inputs_sig.at("atomic_mu_i").name,
                      cppflow::tensor(atomic_mu_i_vector, {tot_atoms}));

  if (has_atomic_mu_i_local) {
    inputs.emplace_back(compute_inputs_sig.at("atomic_mu_i_local").name,
                        cppflow::tensor(atomic_mu_i_vector, {tot_atoms}));
  }

  // map_atoms_to_structure
  if (has_map_atoms_to_structure_op) {
    inputs.emplace_back(compute_inputs_sig.at("map_atoms_to_structure").name,
                        cppflow::tensor(std::vector<int32_t>(tot_atoms, 0), {tot_atoms}));
  }

  // batch_nat = number of extened atoms + padding
  if (has_batch_tot_nat) {
    inputs.emplace_back(compute_inputs_sig.at("batch_tot_nat").name,
                        cppflow::tensor(std::vector<int32_t>{tot_atoms}, {}));
  }

  // batch_nreal_atoms_per_structure: number of extened atoms (w/o padding)
  inputs.emplace_back(compute_inputs_sig.at("batch_tot_nat_real").name,
                      cppflow::tensor(std::vector<int32_t>{nlocal}, {}));

  // ind_i, ind_j: bonds
  //determine the maximum number of neighbours (within cutoff)
  std::vector<int> actual_jnum(inum, 0);
  double cutoff_sq = cutoff * cutoff;
  int type_i, type_j;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    type_i = type[i];
    double xtmp = x[i][0];
    double ytmp = x[i][1];
    double ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    int cur_actual_jnum = 0;
    for (jj = 0; jj < jnum; ++jj) {
      j = jlist[jj];
      j &= NEIGHMASK;
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      if (is_custom_cutoffs) {
        type_j = type[j];
        double cur_cutoff = cutoff_matrix_per_lammps_type[type_i][type_j];
        cutoff_sq = cur_cutoff * cur_cutoff;
      }
      const double rsq = delx * delx + dely * dely + delz * delz;
      if (rsq < cutoff_sq) { cur_actual_jnum += 1; }
    }
    actual_jnum[ii] = cur_actual_jnum;
  }
  n_real_neighbours = std::accumulate(actual_jnum.begin(), actual_jnum.end(), 0);

  std::vector<int> actual_jnum_shift(actual_jnum.size(), 0);
  std::partial_sum(actual_jnum.begin(), actual_jnum.end() - 1, actual_jnum_shift.begin() + 1);

  tot_neighbours = graceimpl->neighbor_padding.update(n_real_neighbours);
  if (graceimpl->neighbor_padding.last_update_triggered_resize()) graceimpl->graph_recompiled = true;
  if (pad_verbose && graceimpl->neighbor_padding.last_update_triggered_resize()) {
    utils::logmesg(lmp,
                   "[GRACE] Neighbours padding: extending new num. of neighbours = {} (incl. "
                   "{:.3f}% fake neighbours)\n",
                   tot_neighbours,
                   100. * (double) graceimpl->neighbor_padding.n_fake / tot_neighbours);
  }

  std::vector<int32_t> ind_i_vector(tot_neighbours);
  std::vector<int32_t> ind_j_vector(tot_neighbours);
  std::vector<int32_t> mu_i_vector(tot_neighbours);
  std::vector<int32_t> mu_j_vector(tot_neighbours);
  std::vector<double> bond_vector(3 * tot_neighbours, 1e6);

  for (ii = 0; ii < inum; ++ii) {
    i = ilist[ii];
    type_i = type[i];
    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    int tot_ind = actual_jnum_shift[ii];
    for (jj = 0; jj < jnum; ++jj) {
      j = jlist[jj];
      j &= NEIGHMASK;
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      if (is_custom_cutoffs) {
        type_j = type[j];
        double cur_cutoff = cutoff_matrix_per_lammps_type[type_i][type_j];
        cutoff_sq = cur_cutoff * cur_cutoff;
      }
      const double rsq = delx * delx + dely * dely + delz * delz;
      if (rsq < cutoff_sq) {
        ind_i_vector[tot_ind] = i;
        // remap j to j_local
        int j_local = atom->map(atom->tag[j]);
        ind_j_vector[tot_ind] = j_local;
        mu_i_vector[tot_ind] = element_type_mapping[type[i]];
        mu_j_vector[tot_ind] = element_type_mapping[type[j]];
        double bondx = atom->x[j][0] - atom->x[i][0];
        double bondy = atom->x[j][1] - atom->x[i][1];
        double bondz = atom->x[j][2] - atom->x[i][2];

        bond_vector[3 * tot_ind + 0] = bondx;
        bond_vector[3 * tot_ind + 1] = bondy;
        bond_vector[3 * tot_ind + 2] = bondz;
        ++tot_ind;
      }
    }
  }

  // add fake bonds
  int fake_atom_ind = tot_atoms - 1;
  for (int tot_ind = n_real_neighbours; tot_ind < tot_neighbours; tot_ind++) {
    ind_i_vector[tot_ind] = fake_atom_ind;    // fake atom ind
    ind_j_vector[tot_ind] = fake_atom_ind;    // fake atom ind
    mu_i_vector[tot_ind] = 0;
    mu_j_vector[tot_ind] = 0;
  }

  inputs.emplace_back(
      compute_inputs_sig.at("ind_i").name,    //DEFAULT_INPUT_PREFIX + "ind_i" + ":0",
      cppflow::tensor(ind_i_vector, {tot_neighbours}));
  inputs.emplace_back(
      compute_inputs_sig.at("ind_j").name,    //DEFAULT_INPUT_PREFIX + "ind_j" + ":0",
      cppflow::tensor(ind_j_vector, {tot_neighbours}));

  // mu_i, mu_j: bonds
  if (has_mu_i_op) {
    inputs.emplace_back(
        compute_inputs_sig.at("mu_i").name,    //DEFAULT_INPUT_PREFIX + "mu_i" + ":0",
        cppflow::tensor(mu_i_vector, {tot_neighbours}));
  }
  inputs.emplace_back(compute_inputs_sig.at("mu_j").name,    //DEFAULT_INPUT_PREFIX + "mu_j" + ":0",
                      cppflow::tensor(mu_j_vector, {tot_neighbours}));

  // num_struc: 1
  if (has_nstruct_total_op) {
    inputs.emplace_back(compute_inputs_sig.at("n_struct_total")
                            .name,    //DEFAULT_INPUT_PREFIX + "n_struct_total" + ":0",
                        cppflow::tensor(std::vector<int32_t>{1}, {}));
  }

  // vector_offsets: 1
  inputs.emplace_back(
      compute_inputs_sig.at("bond_vector").name,    //DEFAULT_INPUT_PREFIX + "bond_vector" + ":0",
      cppflow::tensor(bond_vector, {tot_neighbours, 3}));

#ifdef GRACE_PRINT_DEBUG
  print_tf_inputs(inputs, comm->me, lmp, true);
#endif

  vector<string> output_names;
  if (do_energy_only) {
    //TODO: update!
    auto compute_outputs_sig =
        graceimpl->model->signatures.at(this->compute_energy_only_function_name).outputs;
    output_names = {
        compute_outputs_sig.at("atomic_energy")
            .name,    //"StatefulPartitionedCall:0", // atomic_energy [nat,1]
    };
  } else {
    auto compute_outputs_sig = graceimpl->model->signatures.at(this->compute_function_name).outputs;
    output_names = {
        compute_outputs_sig.at("atomic_energy")
            .name,    //"StatefulPartitionedCall:0", // atomic_energy [nat,1]
        compute_outputs_sig.at("total_energy")
            .name,    //"StatefulPartitionedCall:1", // total_energy [-1, 1]
        compute_outputs_sig.at("total_f")
            .name,    //"StatefulPartitionedCall:2", // total_f [n_at, 3]
        compute_outputs_sig.at("virial").name,    //"StatefulPartitionedCall:3", // virial [6]
    };
    // add it optionally
    if (pair_forces)
      output_names.emplace_back(compute_outputs_sig.at("z_pair_f")
                                    .name);    //"StatefulPartitionedCall:4");// pair_f [n_bonds, 3]
  }
  data_timer.stop();

  //CALL MODEL
  model_timer.start();
  std::vector<cppflow::tensor> output = graceimpl->model->operator()(inputs, output_names);
  model_timer.stop();

  data_timer.start();
  auto &e_out = output[0];    // atomic_energy
  auto e_tens = e_out.get_tensor();
  const double *e_data = static_cast<double *>(TF_TensorData(e_tens.get()));

  if (!do_energy_only) {
    //    auto &te_out = output[1]; // total_energy

    auto &total_f_out = output[2];    // total_f
    auto total_f_tens = total_f_out.get_tensor();
    const double *total_f_data = static_cast<double *>(TF_TensorData(total_f_tens.get()));

    if (!pair_forces) {
      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];

        const int itype = type[i];
        double fx = total_f_data[ii * 3 + 0];
        double fy = total_f_data[ii * 3 + 1];
        double fz = total_f_data[ii * 3 + 2];

        f[i][0] += scale[itype][itype] * fx;
        f[i][1] += scale[itype][itype] * fy;
        f[i][2] += scale[itype][itype] * fz;

        // tally energy contribution
        if (eflag_either) {
          // evdwl = energy of atom I
          evdwl = scale[itype][itype] * e_data[i];
          ev_tally_full(i, 2.0 * evdwl, 0.0, 0.0, 0.0, 0.0, 0.0);
        }
      }    // end for(ii)

      // virial order, seems to be OK
      if (vflag_global) {
        auto &v_out = output[3];    // virial
        auto v_tens = v_out.get_tensor();
        const double *v_data = static_cast<double *>(TF_TensorData(v_tens.get()));

        //            ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, fij[0], fij[1], fij[2], -delx, -dely, -delz);
        virial[0] += v_data[0];
        virial[1] += v_data[1];
        virial[2] += v_data[2];
        virial[3] += v_data[3];
        virial[4] += v_data[4];
        virial[5] += v_data[5];
      }
    } else {
      // pair forces
      auto &f_out = output[4];

      auto f_tens = f_out.get_tensor();
      const double *f_data = static_cast<double *>(TF_TensorData(f_tens.get()));
      int tot_ind = 0;
      for (ii = 0; ii < inum; ++ii) {
        i = ilist[ii];
        type_i = type[i];
        const double xtmp = x[i][0];
        const double ytmp = x[i][1];
        const double ztmp = x[i][2];
        jlist = firstneigh[i];
        jnum = numneigh[i];

        for (jj = 0; jj < jnum; ++jj) {
          j = jlist[jj];
          j &= NEIGHMASK;
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          if (is_custom_cutoffs) {
            type_j = type[j];
            double cur_cutoff = cutoff_matrix_per_lammps_type[type_i][type_j];
            cutoff_sq = cur_cutoff * cur_cutoff;
          }
          const double rsq = delx * delx + dely * dely + delz * delz;
          if (rsq < cutoff_sq) {
            // why "-" ?
            fij[0] = -scale[type_i][type_i] * f_data[tot_ind];
            fij[1] = -scale[type_i][type_i] * f_data[tot_ind + 1];
            fij[2] = -scale[type_i][type_i] * f_data[tot_ind + 2];
            tot_ind += 3;

#ifdef GRACE_PRINT_DEBUG
            // Print pair-specific force mapping
            utils::logmesg(lmp,
                           "[GRACE-FORCE] Proc {}: Pair ({}-{}) | Tag ({}-{}) | BondIdx {} | F_tf: "
                           "[{}, {}, {}]\n",
                           comm->me, i, j, atom->tag[i], atom->tag[j], tot_ind / 3,
                           f_data[tot_ind - 3], f_data[tot_ind - 2], f_data[tot_ind - 1]);
#endif
            f[i][0] += fij[0];
            f[i][1] += fij[1];
            f[i][2] += fij[2];
            f[j][0] -= fij[0];
            f[j][1] -= fij[1];
            f[j][2] -= fij[2];

            // tally per-atom virial contribution, OpenMP critical !
            if (vflag_either) {
              ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, fij[0], fij[1], fij[2], delx, dely,
                           delz);

              // Centroid Stress
              if (cvflag_atom) {
                const double fx = fij[0];
                const double fy = fij[1];
                const double fz = fij[2];

                cvatom[i][0] += 0.5 * delx * fx;    // xx
                cvatom[i][1] += 0.5 * dely * fy;    // yy
                cvatom[i][2] += 0.5 * delz * fz;    // zz
                cvatom[i][3] += 0.5 * delx * fy;    // xy
                cvatom[i][4] += 0.5 * delx * fz;    // xz
                cvatom[i][5] += 0.5 * dely * fz;    // yz
                cvatom[i][6] += 0.5 * dely * fx;    // yx
                cvatom[i][7] += 0.5 * delz * fx;    // zx
                cvatom[i][8] += 0.5 * delz * fy;    // zy

                cvatom[j][0] += 0.5 * delx * fx;    // xx
                cvatom[j][1] += 0.5 * dely * fy;    // yy
                cvatom[j][2] += 0.5 * delz * fz;    // zz
                cvatom[j][3] += 0.5 * delx * fy;    // xy
                cvatom[j][4] += 0.5 * delx * fz;    // xz
                cvatom[j][5] += 0.5 * dely * fz;    // yz
                cvatom[j][6] += 0.5 * dely * fx;    // yx
                cvatom[j][7] += 0.5 * delz * fx;    // zx
                cvatom[j][8] += 0.5 * delz * fy;    // zy
              }
            }
          }
        }    // loop over neighbours

        // tally energy contribution
        if (eflag_either) {
          // evdwl = energy of atom I
          evdwl = scale[type_i][type_i] * e_data[i];
          ev_tally_full(i, 2.0 * evdwl, 0.0, 0.0, 0.0, 0.0, 0.0);
        }

      }    // loop over atoms -i

      if (vflag_fdotr) virial_fdotr_compute();
    }
  } else {
    //energy only computes
    for (ii = 0; ii < inum; ++ii) {
      i = ilist[ii];
      type_i = type[i];
      if (eflag_either) {
        // evdwl = energy of atom I
        evdwl = scale[type_i][type_i] * e_data[i];
        ev_tally_full(i, 2.0 * evdwl, 0.0, 0.0, 0.0, 0.0, 0.0);
      }
    }
  }

  data_timer.stop();
  total_timer.stop();

  if (graceimpl->graph_recompiled) {
    total_timer.rollback();
    data_timer.rollback();
    model_timer.rollback();
  } else {
    total_timer.commit();
    data_timer.commit();
    model_timer.commit();
    total_real_atoms_processed += current_step_real_atoms;
    total_compute_calls++;
  }
}

#endif    //#ifndef NO_GRACE_TF
