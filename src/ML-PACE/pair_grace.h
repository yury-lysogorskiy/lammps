/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/*
Copyright 2024 Yury Lysogorskiy^1,  Anton Bochkarev^1, Ralf Drautz^1

^1: Ruhr-University Bochum, Bochum, Germany
*/

//
// Created by Lysogorskiy Yury on 27.03.24
//

#ifndef NO_GRACE_TF
#ifdef PAIR_CLASS
// clang-format off
PairStyle(grace,PairGRACE);
// clang-format on
#else

#ifndef LMP_PAIR_GRACE_H
#define LMP_PAIR_GRACE_H

#include "pair.h"

#include "utils_pace.h"
#include <map>
#include <set>

namespace LAMMPS_NS {

class PairGRACE : public Pair {
 public:
  PairGRACE(class LAMMPS *);
  ~PairGRACE() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  void *extract(const char *, int &) override;
  void *extract_peratom(const char *, int &) override;

  // Reverse-comm the per-atom uncertainty_force: bond-loop ghost contributions
  // (uncertainty_force[j] -= ... for ghost j, including periodic images) are
  // folded back onto their local owners, exactly as LAMMPS folds the physical
  // force f. Without this the periodic/multi-rank uncertainty force is biased.
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

  // Runtime modification of UQ knobs via `pair_modify`.
  // Recognized: `kappa <val>`, `bias_virial yes|no`.
  // Unknown keywords are forwarded to Pair::modify_params for the standard set.
  void modify_params(int narg, char **arg) override;

 protected:
  struct GRACEImpl *graceimpl;

  std::string DEFAULT_INPUT_PREFIX = "serving_default_";
  std::string compute_function_name = "serving_default";
  const std::string COMPUTE_ENERGY_ONLY_KEY = "compute_energy";
  std::string compute_energy_only_function_name = COMPUTE_ENERGY_ONLY_KEY;
  const std::string COMPUTE_UQ_KEY = "compute_uq";
  const std::string COMPUTE_UQ_GAMMA_ONLY_KEY = "compute_uq_gamma_only";
  bool has_map_atoms_to_structure_op = false;
  bool has_nstruct_total_op = false;
  bool has_mu_i_op = false;
  bool has_batch_tot_nat = false;
  bool has_atomic_mu_i_local = false;
  bool parallel = false;

  // Dense (reshape) engine: opt-in `dense` keyword makes compute() call the model's
  // `compute_dense` signature (per-atom-uniform bond layout, n_bonds == n_atoms*width)
  // instead of `compute` (compact + segment_sum). Same weights, I/O-identical; only the
  // bond-array ordering/padding and z_pair_f indexing differ. UQ / energy-only have no
  // dense signature and fall back to the compact path.
  bool has_compute_dense = false;
  const std::string COMPUTE_DENSE_KEY = "compute_dense";

  bool has_compute_energy_only = false;
  bool warning_compute_energy_only_not_avail_shown = false;
  bool debug_no_energy_only_calc = false;

  // UQ / extrapolation scaffolding. Activation is implicit (parity with grace/kk):
  // UQ runs when the model exports a UQ head (has_compute_uq) AND a per-atom UQ
  // output is requested via `fix pair grace <field>` or kappa != 0. See compute().
  bool has_compute_uq = false;
  // Faster gamma-only signature: skips the dsigma/dr backward pass. When kappa==0
  // and uncertainty_force is not requested, the compute() path prefers this
  // signature over `compute_uq`. Falls back to `compute_uq` when not exported.
  bool has_compute_uq_gamma_only = false;
  // True when the loaded saved model exports the atomic_sigma head. Detected in coeff().
  bool has_atomic_sigma = false;
  bool warning_atomic_sigma_not_avail_shown = false;
  int flag_compute_gamma = 0;             // toggle settable via extract("gamma_flag")
  int flag_compute_gmm_cluster = 0;       // toggle settable via extract("gmm_cluster_flag")
  int flag_compute_uncertainty_force = 0; // toggle settable via extract("uncertainty_force_flag")
  int flag_compute_atomic_sigma = 0;      // toggle settable via extract("atomic_sigma_flag")
  // HAL: f[i] += kappa * (||F_phys_i||+eps) / (||F_sigma_i||+eps) * F_sigma_i (per-atom, post-loop).
  double kappa = 0.0;                     // mutable via extract("kappa")
  // If true, the kappa contribution is also tallied into the GLOBAL virial via a
  // post-loop F^kappa*r positional-virial sum over local atoms; per-atom stress is
  // not affected (the per-atom rescale has no per-bond Newton-3 decomposition).
  bool bias_virial = false;
  // HAL bias gated by (mask[i] & kappa_groupbit). -1 = all atoms, skip the mask
  // check; gamma / gmm_cluster / uncertainty_force are exposed unconditionally.
  int kappa_groupbit = -1;
  char *kappa_group_id = nullptr;
  void parse_kappa_group(const char *name, const char *ctx);
  // Shared setters for kappa / bias_virial, used by both settings() and
  // modify_params() so the two keyword parsers cannot drift.
  void parse_kappa(const char *val, const char *ctx);
  void parse_bias_virial(const char *val);
  const char *kappa_group_str() const {
    return kappa_group_id ? kappa_group_id : "all";
  }
  // Normalization mode for the kappa-rescale denominator. Both forms preserve
  // unit-free kappa; they differ in how the global ratio is computed:
  //   MAX : s = kappa * (max_j||F^phys_j||+eps) / (max_j||F^sigma_j||+eps)
  //   MEAN: s = kappa * (sum_j||F^phys_j||+eps) / (sum_j||F^sigma_j||+eps)
  //         (N's cancel; equivalent to mean/mean).
  // Both apply the resulting scalar uniformly to f[i] += s * F^sigma_i, so
  // per-atom kappa-bias magnitudes stay proportional to ||F^sigma_i||. MEAN
  // amplifies an isolated extrapolating atom more than MAX (because the mean
  // stays small when most atoms are in-distribution), so MEAN is more
  // responsive but also more sensitive to ratio swings.
  enum KappaNormMode { KAPPA_NORM_MAX = 0, KAPPA_NORM_MEAN = 1 };
  KappaNormMode kappa_norm_mode = KAPPA_NORM_MAX;
  // Shared parser used by both PairGRACE::modify_params and the subclass
  // settings(); errors out (no return) on an unknown mode string.
  void parse_kappa_norm(const char *mode, const char *ctx);
  const char *kappa_norm_str() const {
    return kappa_norm_mode == KAPPA_NORM_MAX ? "max" : "mean";
  }
  static constexpr double KAPPA_EPS = 1.0e-8;
  double *gamma = nullptr;                // per-atom extrapolation grade
  double *gmm_cluster_arr = nullptr;      // per-atom GMM cluster index (int32 from TF cast to double for fix pair)
  double **uncertainty_force = nullptr;   // per-atom raw uncertainty force [nmax][3] = -dsigma/dr_i
  double *atomic_sigma = nullptr;         // raw per-atom sigma_i (model output, eV)
  int nmax_uq = 0;                        // current allocation of all per-atom UQ arrays

  virtual void allocate();

  double **scale;
  double cutoff = 6;
  bool is_custom_cutoffs = false;
  vector<vector<double>> cutoff_matrix, cutoff_matrix_per_lammps_type;
  bool pair_forces = true;

  int tot_neighbours = 0;
  int tot_atoms = 0;

  bool dense_enabled = false;              // opt-in `dense` keyword -> use compute_dense
  int dense_width = 0;                     // per-atom neighbour slots this step
  static constexpr int DENSE_TIER = 16;    // width snapped up to a multiple of this

  std::set<int> tot_neighbours_set;
  int max_number_of_reduction = 10, num_of_reductions = 0;

  int chunksize;
  double neigh_padding_fraction = 0.01;
  double reducing_neigh_padding_fraction = 0.2;
  bool do_padding = true;
  bool pad_verbose = false;

  int nelements;
  std::vector<std::string> elements_name;
  std::map<std::string, int> elements_to_index_map;
  std::vector<int> element_type_mapping;    // LAMMPS's type(1,2,3...) to ACE's mu(0,1,2...,89)

  PACE::ACETimer total_timer;
  PACE::ACETimer data_timer;
  PACE::ACETimer model_timer;

  double total_real_atoms_processed = 0.0;
  long long int current_step_real_atoms = 0;
  long long int total_compute_calls = 0;
};
}    // namespace LAMMPS_NS

#endif
#endif
#endif    //#ifndef NO_GRACE_TF
