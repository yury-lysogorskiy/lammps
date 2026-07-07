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
#include "group.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include "yaml-cpp/yaml.h"
#include <cstring>
#include <numeric>
#include <string>
#include <unistd.h>

#include "utils_grace.h"
#include "utils_pace.h"

// CppFlow headers
#include <cppflow/model.h>
#include <cppflow/ops.h>
#include <cppflow/tensor.h>

#ifndef MLPACE_DO_NOT_DISABLE_TFLOAT32
#include "tensorflow/core/platform/tensor_float_32_utils.h"
#endif

#include <tensorflow/c/c_api.h>

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
  std::map<std::string, cppflow::TensorInfo> compute_uq_inputs_sig;
  std::map<std::string, cppflow::TensorInfo> compute_uq_gamma_only_inputs_sig;
  std::map<std::string, cppflow::TensorInfo> compute_dense_inputs_sig;

  // Sticky-max padding for the dense per-atom neighbour-slot count (`width`), so the
  // n_atoms*width tensor shape rarely changes and XLA rarely recompiles.
  GRACE::GracePaddingDimension width_padding;

  // Persistent buffers — only resized when padded counts change
  std::vector<int32_t> atomic_mu_i;
  std::vector<int32_t> ind_i;
  std::vector<int32_t> ind_j;
  std::vector<int32_t> mu_i;
  std::vector<int32_t> mu_j;
  std::vector<double> bond_vector;
  std::vector<int> actual_jnum;
  std::vector<int> actual_jnum_shift;
};

// Slot indices into output[] for UQ-branch tensors. -1 means "not requested
// this step / not present in output_names". Populated alongside the
// emplace_back loop in compute() so there's a single source of truth.
struct UQOutputIdx {
  int gamma = -1;
  int dsigma = -1;
  int gmm = -1;
  int atomic_sigma = -1;
};

}    // namespace LAMMPS_NS

using namespace LAMMPS_NS;
using namespace MathConst;

namespace {
// True if the token is one of the boolean literals utils::logical accepts. Lets
// `bias_virial` accept either a bare-flag form or an explicit yes/no in settings(),
// matching the `pair_modify bias_virial yes|no` spelling.
bool is_logical_token(const char *s)
{
  return !strcmp(s, "yes") || !strcmp(s, "no") || !strcmp(s, "on") || !strcmp(s, "off") ||
      !strcmp(s, "true") || !strcmp(s, "false");
}
}    // namespace

/* ---------------------------------------------------------------------- */
PairGRACE::PairGRACE(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  // Reverse-comm 3 doubles/atom for the per-atom uncertainty_force fold-back
  // (only invoked when uncertainty_force is actually accumulated; see compute()).
  comm_reverse = 3;

  graceimpl = new GRACEImpl;

  scale = nullptr;

  chunksize = 4096;

  total_timer.init();
  data_timer.init();

  no_virial_fdotr_compute = 1;
#ifndef MLPACE_DO_NOT_DISABLE_TFLOAT32
  //disable tensor float 32 execution
  tsl::enable_tensor_float_32_execution(false);
#endif
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */
PairGRACE::~PairGRACE()
{
  if (copymode) return;

  GRACE::log_perf_stats(lmp, "grace", total_real_atoms_processed, total_compute_calls,
                        {{"Total", total_timer.as_microseconds()},
                         {"Data", data_timer.as_microseconds()},
                         {"Model", model_timer.as_microseconds()}});

  delete graceimpl;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(scale);
  }

  if (gamma) memory->destroy(gamma);
  if (gmm_cluster_arr) memory->destroy(gmm_cluster_arr);
  if (uncertainty_force) memory->destroy(uncertainty_force);
  if (atomic_sigma) memory->destroy(atomic_sigma);

  delete[] kappa_group_id;
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
  // ACE potentials are parameterized in metal units
  if (strcmp("metal", update->unit_style) != 0)
    error->all(FLERR, "GRACE potentials require 'metal' units");

  auto tf_version = TF_Version();
  if (comm->me == 0) utils::logmesg(lmp, "[GRACE] TF version: {}\n", tf_version);

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "padding") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "pair_style grace padding", error);
      neigh_padding_fraction = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;

    } else if (strcmp(arg[iarg], "pad_verbose") == 0) {
      pad_verbose = true;
      iarg += 1;
    } else if (strcmp(arg[iarg], "no_pair_forces") == 0) {
      pair_forces = false;
      iarg += 1;
      if (comm->me == 0) utils::logmesg(lmp, "[GRACE] Pair forces are OFF \n");
    } else if (strcmp(arg[iarg], "pair_forces") == 0) {
      pair_forces = true;
      iarg += 1;
      if (comm->me == 0) utils::logmesg(lmp, "[GRACE] Pair forces are ON \n");
    } else if (strcmp(arg[iarg], "max_number_of_reduction") == 0) {
      if (iarg + 1 >= narg)
        utils::missing_cmd_args(FLERR, "pair_style grace max_number_of_reduction", error);
      max_number_of_reduction = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
      if (comm->me == 0)
        utils::logmesg(lmp,
                       "[GRACE] Maximum number of recompilation during padding reduction: {}\n",
                       max_number_of_reduction);
    } else if (strcmp(arg[iarg], "reduce_padding") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "pair_style grace reduce_padding", error);
      reducing_neigh_padding_fraction = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
      if (comm->me == 0)
        utils::logmesg(lmp, "[GRACE] Reducing padding fraction: {}\n",
                       reducing_neigh_padding_fraction);
    } else if (strcmp(arg[iarg], "debug_no_energy_only_calc") == 0) {
      debug_no_energy_only_calc = true;
      iarg += 1;
    } else if (strcmp(arg[iarg], "dense") == 0) {
      dense_enabled = true;
      iarg += 1;
      if (comm->me == 0)
        utils::logmesg(lmp, "[GRACE] Dense (compute_dense) engine requested\n");
    } else if (strcmp(arg[iarg], "kappa") == 0) {
      // UQ / extrapolation knobs (formerly pair_style grace/extrapolation). UQ itself
      // activates implicitly when the model has a UQ head and a `fix pair grace <field>`
      // toggle is set or kappa != 0; these keywords just configure the HAL kappa-rescale.
      if (iarg + 1 >= narg)
        error->all(FLERR,
                   "[GRACE] kappa requires a numeric argument (relative-force uncertainty "
                   "bias coefficient; e.g. 0.1 = 10 % of ||F^phys||)");
      parse_kappa(arg[iarg + 1], "pair_style grace kappa");
      iarg += 2;
    } else if (strcmp(arg[iarg], "bias_virial") == 0) {
      // Accept a bare flag (-> on) or an explicit yes/no, so it matches the
      // `pair_modify bias_virial yes|no` spelling instead of erroring on the value.
      if (iarg + 1 < narg && is_logical_token(arg[iarg + 1])) {
        parse_bias_virial(arg[iarg + 1]);
        iarg += 2;
      } else {
        parse_bias_virial("on");
        iarg += 1;
      }
    } else if (strcmp(arg[iarg], "kappa_norm") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "[GRACE] kappa_norm requires 'max' or 'mean'");
      parse_kappa_norm(arg[iarg + 1], "[GRACE] kappa_norm");
      iarg += 2;
    } else if (strcmp(arg[iarg], "kappa_group") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "[GRACE] kappa_group requires a LAMMPS group name");
      parse_kappa_group(arg[iarg + 1], "[GRACE] kappa_group");
      iarg += 2;
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

  GRACE::warn_padding_reduction_disabled(lmp, neigh_padding_fraction,
                                         reducing_neigh_padding_fraction);

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

  graceimpl->width_padding.enabled = do_padding;
  graceimpl->width_padding.padding_fraction = neigh_padding_fraction;
  graceimpl->width_padding.reduction_threshold_fraction = reducing_neigh_padding_fraction;
  graceimpl->width_padding.max_reductions = max_number_of_reduction;
  graceimpl->width_padding.verbose = pad_verbose;

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

  // check for compute_uq function (UQ / extrapolation head)
  if (graceimpl->model->has_signature(COMPUTE_UQ_KEY)) {
    has_compute_uq = true;
    graceimpl->compute_uq_inputs_sig = graceimpl->model->signatures.at(COMPUTE_UQ_KEY).inputs;
    if (comm->me == 0) utils::logmesg(lmp, "[GRACE] compute_uq signature is available\n");
  }
  // gamma-only UQ head (no dsigma/dr backward pass; used when kappa==0 and
  // uncertainty_force is not requested).
  if (graceimpl->model->has_signature(COMPUTE_UQ_GAMMA_ONLY_KEY)) {
    has_compute_uq_gamma_only = true;
    graceimpl->compute_uq_gamma_only_inputs_sig =
        graceimpl->model->signatures.at(COMPUTE_UQ_GAMMA_ONLY_KEY).inputs;
    if (comm->me == 0)
      utils::logmesg(lmp, "[GRACE] compute_uq_gamma_only signature is available\n");
  }
  // Dense (reshape) engine signature. Same weights / I-O as `compute`; activated by
  // the `dense` keyword (fail-fast at the end of coeff() if requested but absent).
  if (graceimpl->model->has_signature(COMPUTE_DENSE_KEY)) {
    has_compute_dense = true;
    graceimpl->compute_dense_inputs_sig =
        graceimpl->model->signatures.at(COMPUTE_DENSE_KEY).inputs;
    if (comm->me == 0) utils::logmesg(lmp, "[GRACE] compute_dense signature is available\n");
  }
  // NB: UQ-head availability is no longer enforced here. Activation is implicit
  // (parity with grace/kk) -- whether UQ is actually requested is only known at
  // run time (fix-pair flags / kappa), so the "model has no UQ head" error is
  // raised lazily in compute() when a UQ output is first requested.

  if (has_compute_uq) {
    const auto &uq_outs = graceimpl->model->signatures.at(COMPUTE_UQ_KEY).outputs;
    // Raw per-atom sigma -- exported by current UQ-trained models; older ones
    // (e.g. GRACE-2L-OMAT-large_uq_16shards) may omit it. Detected independently.
    has_atomic_sigma = (uq_outs.count("atomic_sigma") > 0);
    if (has_atomic_sigma && comm->me == 0)
      utils::logmesg(lmp, "[GRACE] atomic_sigma head is available\n");
  }

  // Load-time dtype contract: the per-output tensors below are consumed via a
  // hard `static_cast<const double *>` at compute time, so a model exporting any
  // of them as float32 would produce silent garbage. Fail fast at coeff() so the
  // user sees the cause instead of debugging nonsensical forces.
  // (gamma / atomic_sigma have a runtime float<->double dispatch and are
  // intentionally NOT validated here.)
  {
    auto validate_sig = [&](const std::string &sig_name,
                            std::initializer_list<const char *> tensor_names) {
      const auto &outs = graceimpl->model->signatures.at(sig_name).outputs;
      for (const char *tname : tensor_names) {
        auto it = outs.find(tname);
        if (it == outs.end()) continue;    // not exported by this signature; skip
        if (it->second.dtype != TF_DOUBLE)
          error->all(FLERR,
                     "[GRACE] saved-model output '{}' in signature '{}' has dtype id {} "
                     "(expected float64 = {}). This LAMMPS build hard-casts this tensor to "
                     "double*; re-export the model with float64 for this output.",
                     tname, sig_name, static_cast<int>(it->second.dtype),
                     static_cast<int>(TF_DOUBLE));
      }
    };
    validate_sig(compute_function_name, {"atomic_energy", "total_f", "virial", "z_pair_f"});
    if (has_compute_dense)
      validate_sig(COMPUTE_DENSE_KEY, {"atomic_energy", "total_f", "virial", "z_pair_f"});
    if (has_compute_energy_only)
      validate_sig(compute_energy_only_function_name, {"atomic_energy"});
    if (has_compute_uq)
      validate_sig(COMPUTE_UQ_KEY,
                   {"atomic_energy", "total_f", "virial", "pair_f", "dsigma_dr_pair"});
    if (has_compute_uq_gamma_only)
      validate_sig(COMPUTE_UQ_GAMMA_ONLY_KEY, {"atomic_energy", "total_f", "virial", "pair_f"});
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

  if (dense_enabled && !has_compute_dense)
    error->all(FLERR,
               "[GRACE] `dense` was requested but the saved model does not export a "
               "'compute_dense' signature. Re-export the model with the dense engine, or "
               "remove `dense` from the pair_style grace command.");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairGRACE::init_style()
{
  if (atom->tag_enable == 0) error->all(FLERR, "Pair style grace requires atom IDs");
  if (force->newton_pair == 0) error->all(FLERR, "Pair style grace requires newton pair on");

  if (kappa != 0.0 && comm->nprocs > 1)
    error->all(
        FLERR,
        "[GRACE] HAL/kappa-mode (kappa != 0) is not supported with nprocs > 1: "
        "the kappa-rescale uses per-atom force norms that are only complete after the "
        "post-compute reverse_comm, so multi-rank runs would tally biased kappa on "
        "subdomain-boundary atoms. Run single-rank for kappa-mode.");

  if (kappa != 0.0 && !pair_forces)
    error->all(FLERR,
               "[GRACE] kappa != 0 needs per-bond sigma gradients, which require pair-force "
               "mode. Remove 'no_pair_forces' from the pair_style grace command.");

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
  // UQ knobs (only meaningful when the model exports a UQ head; harmless otherwise).
  if (strcmp(str, "gamma_flag") == 0) return (void *) &flag_compute_gamma;
  if (strcmp(str, "gmm_cluster_flag") == 0) return (void *) &flag_compute_gmm_cluster;
  if (strcmp(str, "uncertainty_force_flag") == 0) return (void *) &flag_compute_uncertainty_force;
  if (strcmp(str, "atomic_sigma_flag") == 0) return (void *) &flag_compute_atomic_sigma;
  if (strcmp(str, "kappa") == 0) return (void *) &kappa;

  dim = 2;
  if (strcmp(str, "scale") == 0) return (void *) scale;

  return nullptr;
}

/* ----------------------------------------------------------------------
    runtime modification of UQ knobs via `pair_modify`. Filters out
    `kappa` and `bias_virial`, forwards everything else to Pair::modify_params
    so standard keywords (mix/shift/...) keep working.
 ---------------------------------------------------------------------- */
void PairGRACE::modify_params(int narg, char **arg)
{
  if (narg == 0) utils::missing_cmd_args(FLERR, "pair_modify", error);

  std::vector<char *> forwarded;
  forwarded.reserve(narg);

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "kappa") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "pair_modify kappa", error);
      parse_kappa(arg[iarg + 1], "pair_modify kappa");
      iarg += 2;
    } else if (strcmp(arg[iarg], "bias_virial") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "pair_modify bias_virial", error);
      parse_bias_virial(arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "kappa_norm") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "pair_modify kappa_norm", error);
      parse_kappa_norm(arg[iarg + 1], "pair_modify kappa_norm");
      if (comm->me == 0)
        utils::logmesg(lmp, "[GRACE] kappa_norm = {}\n", kappa_norm_str());
      iarg += 2;
    } else if (strcmp(arg[iarg], "kappa_group") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "pair_modify kappa_group", error);
      parse_kappa_group(arg[iarg + 1], "pair_modify kappa_group");
      if (comm->me == 0) {
        if (kappa_group_id)
          utils::logmesg(lmp,
                         "[GRACE] kappa_group = {} (HAL bias restricted to "
                         "this group; gamma/gmm_cluster/uncertainty_force unaffected)\n",
                         kappa_group_id);
        else
          utils::logmesg(lmp,
                         "[GRACE] kappa_group = all (HAL bias applies to "
                         "every atom; no restriction)\n");
      }
      iarg += 2;
    } else {
      forwarded.push_back(arg[iarg]);
      iarg += 1;
    }
  }

  if (!forwarded.empty()) Pair::modify_params(forwarded.size(), forwarded.data());
}

/* ----------------------------------------------------------------------
    Shared kappa setter used by both PairGRACE::settings and
    PairGRACE::modify_params. Errors if kappa != 0 under nprocs > 1 (the
    kappa-rescale uses per-atom force norms completed only after reverse_comm).
 ---------------------------------------------------------------------- */
void PairGRACE::parse_kappa(const char *val, const char *ctx)
{
  kappa = utils::numeric(FLERR, val, false, lmp);
  if (kappa != 0.0 && comm->nprocs > 1)
    error->all(FLERR,
               "[GRACE] {}: kappa != 0 requires a single MPI rank (the kappa-rescale uses "
               "per-atom force norms completed only after reverse_comm).",
               ctx);
  if (comm->me == 0)
    utils::logmesg(lmp, "[GRACE] kappa = {} (relative-force HAL bias)\n", kappa);
}

/* ----------------------------------------------------------------------
    Shared bias_virial setter used by both PairGRACE::settings and
    PairGRACE::modify_params.
 ---------------------------------------------------------------------- */
void PairGRACE::parse_bias_virial(const char *val)
{
  bias_virial = utils::logical(FLERR, val, false, lmp);
  if (comm->me == 0)
    utils::logmesg(lmp,
                   "[GRACE] bias_virial = {} (kappa contribution {} the global virial; "
                   "per-atom stress unchanged)\n",
                   bias_virial ? "yes" : "no", bias_virial ? "added to" : "excluded from");
}

/* ----------------------------------------------------------------------
    Shared kappa_norm keyword parser used by both PairGRACE::settings and
    PairGRACE::modify_params. Errors (no return) on bad value.
 ---------------------------------------------------------------------- */
void PairGRACE::parse_kappa_norm(const char *mode, const char *ctx)
{
  if (strcmp(mode, "max") == 0)
    kappa_norm_mode = KAPPA_NORM_MAX;
  else if (strcmp(mode, "mean") == 0)
    kappa_norm_mode = KAPPA_NORM_MEAN;
  else
    error->all(FLERR, "{} value must be 'max' or 'mean'", ctx);
}

/* ----------------------------------------------------------------------
    Shared kappa_group keyword parser. The literal name "all" restores the
    default (no mask filter, kappa_groupbit = -1). Any other name is resolved
    to its bitmask via Group::get_bitmask_by_id, which errors out if missing.
 ---------------------------------------------------------------------- */
void PairGRACE::parse_kappa_group(const char *name, const char *ctx)
{
  delete[] kappa_group_id;
  kappa_group_id = nullptr;
  if (strcmp(name, "all") == 0) {
    kappa_groupbit = -1;
    return;
  }
  kappa_groupbit = group->get_bitmask_by_id(FLERR, name, ctx);
  kappa_group_id = utils::strdup(name);
}

/* ----------------------------------------------------------------------
    per-atom UQ extraction (gamma, uncertainty_force)
 ---------------------------------------------------------------------- */
void *PairGRACE::extract_peratom(const char *str, int &ncol)
{
  ncol = 0;
  if (strcmp(str, "gamma") == 0) return (void *) gamma;
  if (strcmp(str, "gmm_cluster") == 0) return (void *) gmm_cluster_arr;
  if (strcmp(str, "atomic_sigma") == 0) return (void *) atomic_sigma;
  if (strcmp(str, "uncertainty_force") == 0) {
    ncol = 3;
    return (void *) uncertainty_force;
  }
  return nullptr;
}

/* ----------------------------------------------------------------------
   reverse-comm of the per-atom uncertainty_force (3 doubles/atom): pack the
   ghost-slot contributions accumulated in the bond loop so comm->reverse_comm()
   can add them onto their local owners. Only invoked from compute() when
   uncertainty_force was actually accumulated (need_dsigma_dr).
------------------------------------------------------------------------- */
int PairGRACE::pack_reverse_comm(int n, int first, double *buf)
{
  // Self-guard: the only call site checks uncertainty_force != null, but comm_reverse=3
  // is advertised unconditionally, so refuse to dereference a null array if ever invoked
  // on a non-UQ run.
  if (!uncertainty_force) return 0;
  int m = 0;
  const int last = first + n;
  for (int i = first; i < last; i++) {
    buf[m++] = uncertainty_force[i][0];
    buf[m++] = uncertainty_force[i][1];
    buf[m++] = uncertainty_force[i][2];
  }
  return m;
}

void PairGRACE::unpack_reverse_comm(int n, int *list, double *buf)
{
  if (!uncertainty_force) return;
  int m = 0;
  for (int i = 0; i < n; i++) {
    const int j = list[i];
    uncertainty_force[j][0] += buf[m++];
    uncertainty_force[j][1] += buf[m++];
    uncertainty_force[j][2] += buf[m++];
  }
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
  int n_real_neighbours;

  int newton_pair = force->newton_pair;

  // inum: length of the neighborlists list
  inum = list->inum;

  // ilist: list of "i" atoms for which neighbor lists exist
  ilist = list->ilist;

  //numneigh: the length of each these neigbor list
  numneigh = list->numneigh;

  // the pointer to the list of neighbors of "i"
  firstneigh = list->firstneigh;

  // Was any UQ output (or kappa) requested? Evaluate from the RAW flags, BEFORE the
  // per-head availability zeroing below, so the "model has no compute_uq head"
  // fail-fast fires consistently for every field -- a gamma request and an
  // atomic_sigma request both error on a non-UQ model, rather than one erroring and
  // the other silently degrading to a confusing downstream "cannot extract" later.
  const bool uq_requested =
      flag_compute_gamma || flag_compute_gmm_cluster || flag_compute_uncertainty_force ||
      flag_compute_atomic_sigma || kappa != 0.0;
  if (uq_requested && !has_compute_uq)
    error->all(FLERR,
               "[GRACE] a UQ output (gamma/uncertainty_force/...) or kappa != 0 was "
               "requested, but the saved model does not export the '{}' signature. Use a GRACE "
               "model trained with a UQ head.",
               COMPUTE_UQ_KEY);

  if (!has_atomic_sigma && flag_compute_atomic_sigma) {
    if (!warning_atomic_sigma_not_avail_shown && comm->me == 0) {
      utils::logmesg(lmp,
                     "[GRACE:WARNING] atomic_sigma head not exported by this saved model; "
                     "atomic_sigma requests are ignored. This message is shown only once.\n");
      warning_atomic_sigma_not_avail_shown = true;
    }
    flag_compute_atomic_sigma = 0;
  }

  // UQ / extrapolation path takes precedence over energy-only when active. do_uq
  // uses the POST-zeroing flags (a request for a head this model lacks was cleared
  // above), so it won't spin up the UQ path for a field that cannot be produced.
  // Activation is implicit (parity with grace/kk): plain `grace` -- and any run with
  // no surviving UQ field requested -- never enters this path and dispatches the
  // regular compute / compute_energy signature below at no extra per-step cost.
  bool do_uq = has_compute_uq &&
      (flag_compute_gamma || flag_compute_gmm_cluster || flag_compute_uncertainty_force ||
       flag_compute_atomic_sigma || kappa != 0.0);
  // dsigma/dr is only needed for the kappa-rescale and for fix-pair uncertainty_force.
  // Otherwise prefer the faster gamma-only signature (skips the σ-grad backward pass)
  // when the model exports it; fall back to compute_uq when it doesn't.
  const bool need_dsigma_dr = do_uq && (kappa != 0.0 || flag_compute_uncertainty_force);
  // kappa-rescale and the per-atom uncertainty_force tally both need per-bond
  // sigma-gradients, which only exist in pair-force mode (the old grace/extrapolation
  // subclass forced this). Plain gamma/sigma/gmm do NOT need forces.
  if (need_dsigma_dr && !pair_forces)
    error->all(FLERR,
               "[GRACE] kappa != 0 or an uncertainty_force request needs per-bond sigma "
               "gradients, which require pair-force mode. Remove 'no_pair_forces' from the "
               "pair_style grace command.");
  const bool use_uq_gamma_only = do_uq && !need_dsigma_dr && has_compute_uq_gamma_only;
  bool do_energy_only = eflag_only && !debug_no_energy_only_calc && !do_uq;

  // Dense engine: only for the standard compute path. UQ / energy-only have no
  // compute_dense signature, so they keep the compact layout + their own signatures.
  const bool use_dense = dense_enabled && has_compute_dense && !do_uq && !do_energy_only;

  if (do_uq) {
    // (re)allocate per-atom UQ arrays sized to atom->nmax (covers ghosts too)
    if (atom->nmax > nmax_uq) {
      // Flag-gated UQ outputs are filled only when their flag is set. Zero on
      // grow so a fix pair / dump that reads them without the trigger sees 0
      // instead of uninitialized memory.
      auto grow_zero = [&](double *&p, const char *tag) {
        memory->grow(p, atom->nmax, tag);
        std::memset(p, 0, sizeof(double) * atom->nmax);
      };
      memory->grow(gamma, atom->nmax, "grace:gamma");
      memory->grow(uncertainty_force, atom->nmax, 3, "grace:uncertainty_force");
      // Zero uncertainty_force on (re)grow: a fix pair / dump can map it without its
      // own trigger set (need_dsigma_dr false below), in which case the per-step
      // memset is skipped and the array is never written -- a reader must see 0, not
      // malloc garbage. (gamma needs no such zeroing: fully overwritten over nlocal.)
      std::memset(&uncertainty_force[0][0], 0, sizeof(double) * 3 * atom->nmax);
      grow_zero(gmm_cluster_arr, "grace:gmm_cluster");
      grow_zero(atomic_sigma, "grace:atomic_sigma");
      nmax_uq = atom->nmax;
    }
    // gamma needs no zeroing: it is fully overwritten over [0, nlocal) post-call,
    // and ghosts are never consumed. uncertainty_force, however, is bond-loop
    // accumulated (locals +=, ghosts -= under newton-on) and only needs to be
    // zeroed when it will actually be populated -- i.e. when need_dsigma_dr is
    // set (kappa-rescale or fix pair uncertainty_force).
    if (need_dsigma_dr && uncertainty_force) {
      const int n_active = atom->nlocal + atom->nghost;
      std::memset(&uncertainty_force[0][0], 0, sizeof(double) * 3 * n_active);
    }
  } else if (do_energy_only && !has_compute_energy_only) {
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

  const auto &compute_inputs_sig = do_uq
      ? (use_uq_gamma_only ? graceimpl->compute_uq_gamma_only_inputs_sig
                           : graceimpl->compute_uq_inputs_sig)
      : (use_dense ? graceimpl->compute_dense_inputs_sig
                   : (do_energy_only && has_compute_energy_only
                          ? graceimpl->compute_energy_only_inputs_sig
                          : graceimpl->compute_inputs_sig));

  // utils::logmesg(lmp, "[GRACE-debug] input_prefix={} \n",input_prefix);

  data_timer.start();
  std::vector<std::tuple<std::string, cppflow::tensor>> inputs;

  tot_atoms = graceimpl->atom_padding.update(nlocal);
  if (graceimpl->atom_padding.last_update_triggered_resize()) graceimpl->graph_recompiled = true;
  if (pad_verbose && graceimpl->atom_padding.last_update_triggered_resize()) {
    utils::logmesg(
        lmp,
        "[GRACE] step {} Atoms padding: {} new num. of atoms = {} (incl. {:.3f}% fake atoms)\n",
        update->ntimestep, graceimpl->atom_padding.action_str(), tot_atoms,
        100. * (double) graceimpl->atom_padding.n_fake / tot_atoms);
  }
  current_step_real_atoms = nlocal;

  // atomic_mu_i: per-atom species type + padding with type[0]
  graceimpl->atomic_mu_i.resize(tot_atoms, element_type_mapping[type[0]]);
  auto &atomic_mu_i_vector = graceimpl->atomic_mu_i;
  std::fill(atomic_mu_i_vector.begin(), atomic_mu_i_vector.end(), element_type_mapping[type[0]]);
  for (i = 0; i < nlocal; ++i) atomic_mu_i_vector[i] = element_type_mapping[type[i]];

  inputs.emplace_back(compute_inputs_sig.at("atomic_mu_i").name,
                      cppflow::tensor(atomic_mu_i_vector, {tot_atoms}));

  if (compute_inputs_sig.count("atomic_mu_i_local")) {
    inputs.emplace_back(compute_inputs_sig.at("atomic_mu_i_local").name,
                        cppflow::tensor(atomic_mu_i_vector, {tot_atoms}));
  }

  // map_atoms_to_structure
  if (compute_inputs_sig.count("map_atoms_to_structure")) {
    inputs.emplace_back(compute_inputs_sig.at("map_atoms_to_structure").name,
                        cppflow::tensor(std::vector<int32_t>(tot_atoms, 0), {tot_atoms}));
  }

  // batch_nat = number of extened atoms + padding
  if (compute_inputs_sig.count("batch_tot_nat")) {
    inputs.emplace_back(compute_inputs_sig.at("batch_tot_nat").name,
                        cppflow::tensor(std::vector<int32_t>{(int32_t) tot_atoms}, {}));
  }

  // batch_nreal_atoms_per_structure: number of extened atoms (w/o padding)
  inputs.emplace_back(compute_inputs_sig.at("batch_tot_nat_real").name,
                      cppflow::tensor(std::vector<int32_t>{nlocal}, {}));

  // ind_i, ind_j: bonds
  //determine the maximum number of neighbours (within cutoff)
  graceimpl->actual_jnum.resize(inum);
  auto &actual_jnum = graceimpl->actual_jnum;
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

  // Bond arrays (persistent buffers). Aliased once here; both engines resize + fill
  // them below. The references stay valid across resize (they alias the vector object).
  auto &ind_i_vector = graceimpl->ind_i;
  auto &ind_j_vector = graceimpl->ind_j;
  auto &mu_i_vector = graceimpl->mu_i;
  auto &mu_j_vector = graceimpl->mu_j;
  auto &bond_vector = graceimpl->bond_vector;

  if (use_dense) {
    // Per-atom-uniform (reshape) layout: n_bonds == tot_atoms * dense_width. width =
    // max neighbours/atom snapped up to DENSE_TIER, then sticky-max padded to bound XLA
    // recompiles. The central atom of each bond is fixed BY POSITION
    // (flat = i*dense_width + slot); real neighbours fill an atom's first slots, the
    // rest stay dummy (bond_vector > cutoff -> zeroed by the radial envelope).
    int true_max = 0;
    for (ii = 0; ii < inum; ii++) true_max = std::max(true_max, actual_jnum[ii]);
    int snapped = ((std::max(true_max, 1) + DENSE_TIER - 1) / DENSE_TIER) * DENSE_TIER;
    dense_width = graceimpl->width_padding.update(snapped);
    if (graceimpl->width_padding.last_update_triggered_resize())
      graceimpl->graph_recompiled = true;
    if (pad_verbose && graceimpl->width_padding.last_update_triggered_resize())
      utils::logmesg(lmp, "[GRACE] step {} Dense width: {} new width = {} (true_max = {})\n",
                     update->ntimestep, graceimpl->width_padding.action_str(), dense_width,
                     true_max);

    tot_neighbours = tot_atoms * dense_width;

    ind_i_vector.resize(tot_neighbours);
    ind_j_vector.resize(tot_neighbours);
    mu_i_vector.resize(tot_neighbours);
    mu_j_vector.resize(tot_neighbours);
    bond_vector.resize(3 * tot_neighbours);

    // Prefill every slot with a dummy: bond_vector > cutoff (radial envelope zeros it),
    // species/indices 0. Real neighbours overwrite the first slots of each atom's block.
    const double dense_dummy = std::max(52.0, cutoff + 1.0);
    std::fill(ind_i_vector.begin(), ind_i_vector.end(), 0);
    std::fill(ind_j_vector.begin(), ind_j_vector.end(), 0);
    std::fill(mu_i_vector.begin(), mu_i_vector.end(), 0);
    std::fill(mu_j_vector.begin(), mu_j_vector.end(), 0);
    std::fill(bond_vector.begin(), bond_vector.end(), dense_dummy);

    for (ii = 0; ii < inum; ++ii) {
      i = ilist[ii];
      type_i = type[i];
      const double xtmp = x[i][0];
      const double ytmp = x[i][1];
      const double ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];
      int slot = 0;
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
          const int dst = i * dense_width + slot;    // atom-major fixed stride
          ind_i_vector[dst] = i;
          ind_j_vector[dst] = atom->map(atom->tag[j]);
          mu_i_vector[dst] = element_type_mapping[type[i]];
          mu_j_vector[dst] = element_type_mapping[type[j]];
          bond_vector[3 * dst + 0] = x[j][0] - x[i][0];
          bond_vector[3 * dst + 1] = x[j][1] - x[i][1];
          bond_vector[3 * dst + 2] = x[j][2] - x[i][2];
          ++slot;
        }
      }
    }
  } else {
    graceimpl->actual_jnum_shift.resize(inum, 0);
    auto &actual_jnum_shift = graceimpl->actual_jnum_shift;
    std::fill(actual_jnum_shift.begin(), actual_jnum_shift.end(), 0);
    std::partial_sum(actual_jnum.begin(), actual_jnum.end() - 1, actual_jnum_shift.begin() + 1);

    tot_neighbours = graceimpl->neighbor_padding.update(n_real_neighbours);
    if (graceimpl->neighbor_padding.last_update_triggered_resize())
      graceimpl->graph_recompiled = true;
    if (pad_verbose && graceimpl->neighbor_padding.last_update_triggered_resize()) {
      utils::logmesg(lmp,
                     "[GRACE] step {} Neighbours padding: {} new num. of neighbours = {} "
                     "(incl. {:.3f}% fake neighbours)\n",
                     update->ntimestep, graceimpl->neighbor_padding.action_str(), tot_neighbours,
                     100. * (double) graceimpl->neighbor_padding.n_fake / tot_neighbours);
    }

    // Resize persistent buffers only when padded neighbour count changes
    ind_i_vector.resize(tot_neighbours);
    ind_j_vector.resize(tot_neighbours);
    mu_i_vector.resize(tot_neighbours);
    mu_j_vector.resize(tot_neighbours);
    bond_vector.resize(3 * tot_neighbours, 1e6);
    // Reset fake-bond sentinel for the padding zone
    std::fill(bond_vector.begin() + 3 * n_real_neighbours, bond_vector.end(), 1e6);

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
  }

  inputs.emplace_back(
      compute_inputs_sig.at("ind_i").name,    //DEFAULT_INPUT_PREFIX + "ind_i" + ":0",
      cppflow::tensor(ind_i_vector, {tot_neighbours}));
  inputs.emplace_back(
      compute_inputs_sig.at("ind_j").name,    //DEFAULT_INPUT_PREFIX + "ind_j" + ":0",
      cppflow::tensor(ind_j_vector, {tot_neighbours}));

  // mu_i, mu_j: bonds
  if (compute_inputs_sig.count("mu_i")) {
    inputs.emplace_back(
        compute_inputs_sig.at("mu_i").name,    //DEFAULT_INPUT_PREFIX + "mu_i" + ":0",
        cppflow::tensor(mu_i_vector, {tot_neighbours}));
  }
  inputs.emplace_back(compute_inputs_sig.at("mu_j").name,    //DEFAULT_INPUT_PREFIX + "mu_j" + ":0",
                      cppflow::tensor(mu_j_vector, {tot_neighbours}));

  // num_struc: 1
  if (compute_inputs_sig.count("n_struct_total")) {
    inputs.emplace_back(compute_inputs_sig.at("n_struct_total")
                            .name,    //DEFAULT_INPUT_PREFIX + "n_struct_total" + ":0",
                        cppflow::tensor(std::vector<int32_t>{1}, {}));
  }

  // vector_offsets: 1
  inputs.emplace_back(
      compute_inputs_sig.at("bond_vector").name,    //DEFAULT_INPUT_PREFIX + "bond_vector" + ":0",
      cppflow::tensor(bond_vector, {tot_neighbours, 3}));

  // map_bonds_to_structure: per-bond -> structure index (single-structure: all zeros)
  if (compute_inputs_sig.count("map_bonds_to_structure")) {
    inputs.emplace_back(compute_inputs_sig.at("map_bonds_to_structure").name,
                        cppflow::tensor(std::vector<int32_t>(tot_neighbours, 0), {tot_neighbours}));
  }

#ifdef GRACE_PRINT_DEBUG
  print_tf_inputs(inputs, comm->me, lmp, true);
#endif

  vector<string> output_names;
  UQOutputIdx uq_idx;
  if (do_energy_only) {
    //TODO: update!
    auto compute_outputs_sig =
        graceimpl->model->signatures.at(this->compute_energy_only_function_name).outputs;
    output_names = {
        compute_outputs_sig.at("atomic_energy").name,    // atomic_energy [nat,1]
    };
  } else if (do_uq) {
    const std::string &uq_key = use_uq_gamma_only ? COMPUTE_UQ_GAMMA_ONLY_KEY : COMPUTE_UQ_KEY;
    auto compute_outputs_sig = graceimpl->model->signatures.at(uq_key).outputs;
    // Append a named tensor to output_names and return its slot index.
    auto append = [&](const char *tname) {
      output_names.emplace_back(compute_outputs_sig.at(tname).name);
      return static_cast<int>(output_names.size()) - 1;
    };
    // Index layout matches the regular branch for [0..4] so shared decoding code stays valid.
    append("atomic_energy");    // [0]
    append("total_energy");     // [1]
    append("total_f");          // [2]
    append("virial");           // [3]
    append("pair_f");           // [4]  (UQ analog of z_pair_f)
    uq_idx.gamma = append("gamma");
    // dsigma_dr_pair only exists in compute_uq (not in compute_uq_gamma_only)
    // and is only consumed when kappa != 0 OR uncertainty_force is requested.
    if (need_dsigma_dr) uq_idx.dsigma = append("dsigma_dr_pair");
    // gmm_cluster is available in both UQ signatures; only fetch on request.
    if (flag_compute_gmm_cluster) uq_idx.gmm = append("gmm_cluster");
    // Raw per-atom sigma -- gated on has_atomic_sigma (set in coeff()).
    if (flag_compute_atomic_sigma) uq_idx.atomic_sigma = append("atomic_sigma");
  } else {
    // Dense reuses the compute_dense signature's outputs (I/O-identical to compute).
    const std::string &out_key = use_dense ? COMPUTE_DENSE_KEY : this->compute_function_name;
    auto compute_outputs_sig = graceimpl->model->signatures.at(out_key).outputs;
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
  const double *e_data = static_cast<const double *>(TF_TensorData(e_tens.get()));

  if (!do_energy_only) {
    //    auto &te_out = output[1]; // total_energy

    auto &total_f_out = output[2];    // total_f
    auto total_f_tens = total_f_out.get_tensor();
    const double *total_f_data = static_cast<const double *>(TF_TensorData(total_f_tens.get()));

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
        const double *v_data = static_cast<const double *>(TF_TensorData(v_tens.get()));

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
      const double *f_data = static_cast<const double *>(TF_TensorData(f_tens.get()));

      // UQ: fetch per-bond dsigma/dr whenever kappa-mixing is on OR fix pair has
      // requested per-atom uncertainty_force (both fold into need_dsigma_dr).
      // The slot index is set by the output_names builder; -1 means absent.
      const double *dsigma_dr_pair_data = nullptr;
      if (need_dsigma_dr) {
        dsigma_dr_pair_data =
            static_cast<const double *>(TF_TensorData(output[uq_idx.dsigma].get_tensor().get()));
      }

      int tot_ind = 0;
      for (ii = 0; ii < inum; ++ii) {
        i = ilist[ii];
        type_i = type[i];
        const double xtmp = x[i][0];
        const double ytmp = x[i][1];
        const double ztmp = x[i][2];
        jlist = firstneigh[i];
        jnum = numneigh[i];
        int slot_i = 0;    // dense per-atom slot counter (unused in compact mode)

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
            // z_pair_f is dE/dr (bond gradient); physical force on i is -dE/dr.
            // The bond loop tallies pure physical forces only -- kappa-mixing is
            // applied per-atom AFTER the loop (it depends on per-atom force
            // norms, which are not bond-local). Virial follows fij.
            // Dense z_pair_f is in per-atom-uniform order (flat = i*width + slot);
            // compact is the running tot_ind stride.
            const int base = use_dense ? 3 * (i * dense_width + slot_i) : tot_ind;
            const double s = scale[type_i][type_i];
            fij[0] = -s * f_data[base + 0];
            fij[1] = -s * f_data[base + 1];
            fij[2] = -s * f_data[base + 2];
            // Per-bond uncertainty force fij_unc = -s*dsigma/dr tallied per-atom
            // (newton-on) for the post-loop kappa-rescale and for fix pair output.
            if (dsigma_dr_pair_data) {
              const double fij_unc_x = -s * dsigma_dr_pair_data[tot_ind];
              const double fij_unc_y = -s * dsigma_dr_pair_data[tot_ind + 1];
              const double fij_unc_z = -s * dsigma_dr_pair_data[tot_ind + 2];
              uncertainty_force[i][0] += fij_unc_x;
              uncertainty_force[i][1] += fij_unc_y;
              uncertainty_force[i][2] += fij_unc_z;
              uncertainty_force[j][0] -= fij_unc_x;
              uncertainty_force[j][1] -= fij_unc_y;
              uncertainty_force[j][2] -= fij_unc_z;
            }
            tot_ind += 3;
            ++slot_i;

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
            // The bond loop's `fij` is the pure physical force in kappa-mode, so
            // the virial and stress are pulled directly from fij with no
            // post-loop fix-up. The kappa contribution affects f[i] only.
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

  // UQ: extract per-atom gamma + optional predicted-error heads + gmm_cluster.
  if (do_uq) {
    // Float32/float64 dispatch -- model dtype varies; a wrong cast silently
    // reinterprets the buffer and returns garbage. Used for every float head.
    auto copy_float_or_double = [&](TF_Tensor *tens, double *dest, const char *name) {
      const TF_DataType dt = TF_TensorType(tens);
      if (dt == TF_FLOAT) {
        const float *src = static_cast<const float *>(TF_TensorData(tens));
        for (i = 0; i < nlocal; ++i) dest[i] = static_cast<double>(src[i]);
      } else if (dt == TF_DOUBLE) {
        const double *src = static_cast<const double *>(TF_TensorData(tens));
        for (i = 0; i < nlocal; ++i) dest[i] = src[i];
      } else {
        error->all(FLERR,
                   "[GRACE] unexpected dtype id {} for output '{}' "
                   "(expected float32 or float64)",
                   static_cast<int>(dt), name);
      }
    };

    copy_float_or_double(output[uq_idx.gamma].get_tensor().get(), gamma, "gamma");

    if (uq_idx.atomic_sigma >= 0)
      copy_float_or_double(output[uq_idx.atomic_sigma].get_tensor().get(), atomic_sigma,
                           "atomic_sigma");

    if (flag_compute_gmm_cluster) {
      // gmm_cluster is int32 in the model and stored as double for fix pair.
      const int32_t *gmm_data =
          static_cast<const int32_t *>(TF_TensorData(output[uq_idx.gmm].get_tensor().get()));
      for (i = 0; i < nlocal; ++i) gmm_cluster_arr[i] = static_cast<double>(gmm_data[i]);
    }
  }

  // Fold ghost uncertainty_force contributions back onto their local owners.
  // The bond loop tallies uncertainty_force[j] -= ... for neighbor j, which on a
  // periodic cell (or multi-rank subdomain) is a ghost image of some local atom;
  // without this those contributions are stranded on ghost slots and the local
  // uncertainty_force is biased (observed Pearson ~0.87 vs the exact non-periodic
  // value). This mirrors how LAMMPS reverse-communicates the physical force f.
  // Run before the kappa-rescale (which reads the per-atom F^sigma norm) and
  // before fix pair / dump read uncertainty_force post-compute.
  if (need_dsigma_dr && uncertainty_force) comm->reverse_comm(this);

  // kappa-rescale (relative-force HAL bias). Per-atom, applied after the bond
  // loop so we can use the assembled per-atom force norms:
  //   F^kappa_i = kappa * (||F^phys_i|| + eps) / (||F^sigma_i|| + eps) * F^sigma_i
  //   f[i] += F^kappa_i
  // do_hal_forces: gate on the actual mutation of f[i].
  // do_hal_virial: opt-in `bias_virial` -- tally sum F^kappa_i*r_i into the GLOBAL
  //   virial via the same one-body F*r convention as
  //   Pair::virial_fdotr_compute. Per-atom stress is intentionally untouched
  //   (the per-atom rescale has no Newton-3 per-bond decomposition).
  // kappa-mode is forbidden under nprocs > 1 (init_style errors out), so the
  // per-atom force norms below are guaranteed complete.
  const bool do_hal_forces = do_uq && kappa != 0.0 && pair_forces;
  const bool do_hal_virial = do_hal_forces && bias_virial && vflag_global;
  if (do_hal_forces) {
    // Reduce per-atom F^phys, F^sigma magnitudes into a single global ratio
    // (see KappaNormMode in pair_grace.h for the WHY of MAX vs MEAN).
    const int *const mask = atom->mask;
    const int kgbit = kappa_groupbit;
    double num = 0.0, den = 0.0;
    if (kappa_norm_mode == KAPPA_NORM_MAX) {
      for (i = 0; i < nlocal; ++i) {
        if (kgbit >= 0 && !(mask[i] & kgbit)) continue;
        num = std::max(num, MathExtra::len3(f[i]));
        den = std::max(den, MathExtra::len3(uncertainty_force[i]));
      }
    } else {    // KAPPA_NORM_MEAN: N's cancel, sum/sum suffices
      for (i = 0; i < nlocal; ++i) {
        if (kgbit >= 0 && !(mask[i] & kgbit)) continue;
        num += MathExtra::len3(f[i]);
        den += MathExtra::len3(uncertainty_force[i]);
      }
    }
    const double s = kappa * (num + KAPPA_EPS) / (den + KAPPA_EPS);
    for (i = 0; i < nlocal; ++i) {
      if (kgbit >= 0 && !(mask[i] & kgbit)) continue;
      const double fkx = s * uncertainty_force[i][0];
      const double fky = s * uncertainty_force[i][1];
      const double fkz = s * uncertainty_force[i][2];
      f[i][0] += fkx;
      f[i][1] += fky;
      f[i][2] += fkz;
      if (do_hal_virial) {
        // F^kappa * r positional virial, mirroring virial_fdotr_compute layout:
        //   xx, yy, zz, xy, xz, yz  <->  fy*x, fz*x, fz*y on the off-diagonals
        virial[0] += fkx * x[i][0];
        virial[1] += fky * x[i][1];
        virial[2] += fkz * x[i][2];
        virial[3] += fky * x[i][0];
        virial[4] += fkz * x[i][0];
        virial[5] += fkz * x[i][1];
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
