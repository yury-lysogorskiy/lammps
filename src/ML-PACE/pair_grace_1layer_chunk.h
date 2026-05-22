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
#ifdef PAIR_CLASS
// clang-format off
PairStyle(grace/1layer/chunk,PairGRACE1LayerChunk);
// clang-format on
#else

#ifndef LMP_PAIR_GRACE_1LAYER_CHUNK_H
#define LMP_PAIR_GRACE_1LAYER_CHUNK_H

#include "pair.h"
#include "utils_pace.h"
#include <map>
#include <set>
#include <string>
#include <vector>

namespace LAMMPS_NS {

class PairGRACE1LayerChunk : public Pair {
 public:
  PairGRACE1LayerChunk(class LAMMPS *);
  virtual ~PairGRACE1LayerChunk();

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  void *extract(const char *, int &) override;

 protected:
  struct GRACE1LayerChunkImpl *impl;

  std::string DEFAULT_INPUT_PREFIX = "serving_default_";
  std::string compute_function_name = "serving_default";
  const std::string COMPUTE_ENERGY_ONLY_KEY = "compute_energy";
  std::string compute_energy_only_function_name = COMPUTE_ENERGY_ONLY_KEY;

  bool has_map_atoms_to_structure_op = false;
  bool has_nstruct_total_op = false;
  bool has_mu_i_op = false;
  bool has_batch_tot_nat = false;
  bool has_atomic_mu_i_local = false;

  bool has_compute_energy_only = false;
  bool warning_compute_energy_only_not_avail_shown = false;
  bool debug_no_energy_only_calc = false;

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

  int nelements;
  std::vector<std::string> elements_name;
  std::map<std::string, int> elements_to_index_map;
  std::vector<int> element_type_mapping;    // LAMMPS's type(1,2,3...) to ACE's mu

  PACE::ACETimer total_timer;
  PACE::ACETimer data_timer;
  PACE::ACETimer model_timer;
  PACE::ACETimer tp_timer;

  double total_real_atoms_processed = 0.0;
  long long int current_step_real_atoms = 0;
  long long int total_compute_calls = 0;
};

}    // namespace LAMMPS_NS

#endif
#endif
#endif    //#ifndef NO_GRACE_TF
