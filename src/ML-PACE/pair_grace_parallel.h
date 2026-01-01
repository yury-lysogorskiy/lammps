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
PairStyle(grace/parallel,PairGRACEParallel);
// clang-format on
#else

#ifndef LMP_PAIR_GRACE_PARALLEL_H
#define LMP_PAIR_GRACE_PARALLEL_H

#include "pair.h"

#include <map>
#include <set>
#include "utils_pace.h"

namespace LAMMPS_NS {

class PairGRACEParallel : public Pair {
 public:
  PairGRACEParallel(class LAMMPS *);
  ~PairGRACEParallel() override;


  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  void *extract(const char *, int &) override;

 protected:
  struct GRACEImpl *aceimpl;

  std::string DEFAULT_INPUT_PREFIX = "parallel_compute_";
  bool has_map_atoms_to_structure_op = false;
  bool has_nstruct_total_op = false;
  bool has_mu_i_op = false;
  bool has_batch_tot_nat = false;
  bool parallel = false;

  virtual void allocate();

  double **scale;
  double cutoff = 6;
  bool is_custom_cutoffs = false;
  vector<vector<double>> cutoff_matrix, cutoff_matrix_per_lammps_type;
  bool pair_forces = true;

  //int tot_neighbours = 0;
  //int tot_atoms = 0;


  //std::set<int> tot_neighbours_set;
  //int max_number_of_reduction = 10, num_of_reductions = 0;

  int chunksize;
  //double neigh_padding_fraction = 0.01;
  //double reducing_neigh_padding_fraction = 0.2;
  //bool do_padding = true;
  //bool pad_verbose = false;

  int nelements;
  std::vector<std::string> elements_name;
  std::map<std::string, int> elements_to_index_map;
  std::vector<int> element_type_mapping; // LAMMPS's type(1,2,3...) to ACE's mu(0,1,2...,89)

  PACE::ACETimer data_timer;
  PACE::ACETimer tp_timer;

  std::vector<int> get_atom_shell_mapping() const;
  void print_atomic_neighbours(const std::vector<int>& atom_type_map)  const;
};
}    // namespace LAMMPS_NS

#endif
#endif
#endif //#ifndef NO_GRACE_TF
