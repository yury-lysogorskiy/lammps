/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef NO_GRACE_TF
#ifdef PAIR_CLASS
// clang-format off
PairStyle(grace/2layer/parallel,PairGRACE2LayerParallel);
// clang-format on
#else

#ifndef LMP_PAIR_GRACE_2LAYER_PARALLEL_H
#define LMP_PAIR_GRACE_2LAYER_PARALLEL_H

#include "pair.h"
#include <vector>
#include <string>
#include <map>
#include "utils_pace.h"
#include <cppflow/tensor.h>

namespace LAMMPS_NS {

class PairGRACE2LayerParallel : public Pair {
 public:
  PairGRACE2LayerParallel(class LAMMPS *);
  ~PairGRACE2LayerParallel() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

  void *extract(const char *, int &) override;

 protected:
  struct GRACE2LayerImpl *aceimpl;

  std::string DEFAULT_INPUT_PREFIX = "forward_layer_1_";
  std::string forward_layer_1_name = "forward_layer_1";
  std::string backward_layer_2_name = "backward_layer_2";
  std::string backward_layer_1_name = "backward_layer_1";

  const std::string ENERGY_KEY = "atomic_energy";
  const std::string GRAD_BOND_KEY = "grad_bond_vector";

  bool has_mu_i_op = false;
  bool has_batch_tot_nat = false;

  void allocate();
  
  double **scale;
  double cutoff = 6.0;
  bool is_custom_cutoffs = false;
  std::vector<std::vector<double>> cutoff_matrix, cutoff_matrix_per_lammps_type;

  int chunksize;
  double neigh_padding_fraction = 0.01;
  double reducing_neigh_padding_fraction = 0.2;
  int max_number_of_reduction = 10;
  bool do_padding = true;
  bool pad_verbose = false;
  bool pair_forces = false;
  
  int nelements;
  std::vector<std::string> elements_name;
  std::map<std::string, int> elements_to_index_map;
  std::vector<int> element_type_mapping;

  // Per-atom features (forward comm) and gradients (reverse comm)
  std::map<std::string, std::vector<double>> features;
  std::map<std::string, std::vector<double>> gradients;
  std::map<std::string, std::vector<int64_t>> feature_shapes;
  std::map<std::string, int> feature_sizes;
  std::map<std::string, bool> feature_is_local;

  std::vector<double> grad_bv_L2;

// #ifdef GRACE_PROFILE
  PACE::ACETimer total_timer;
  PACE::ACETimer data_timer;
  PACE::ACETimer tp_timer;
  PACE::ACETimer comm_timer;
  PACE::ACETimer model1_timer;
  PACE::ACETimer model2_timer;
  PACE::ACETimer model3_timer;
// #endif

  // Helper methods to match compute() phases
  void run_forward_layer_1();
  void run_backward_layer_2(int eflag, int vflag);
  void run_backward_layer_1();

  void print_tensors(const std::string& name, const std::vector<std::tuple<std::string, cppflow::tensor>>& tensors, const std::string& type_prefix = "Input");

};

} // namespace LAMMPS_NS

#endif
#endif
#endif //#ifndef NO_GRACE_TF
