/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/*
Copyright 2021 Yury Lysogorskiy^1, Cas van der Oord^2, Anton Bochkarev^1,
 Sarath Menon^1, Matteo Rinaldi^1, Thomas Hammerschmidt^1, Matous Mrovec^1,
 Aidan Thompson^3, Gabor Csanyi^2, Christoph Ortner^4, Ralf Drautz^1

^1: Ruhr-University Bochum, Bochum, Germany
^2: University of Cambridge, Cambridge, United Kingdom
^3: Sandia National Laboratories, Albuquerque, New Mexico, USA
^4: University of British Columbia, Vancouver, BC, Canada
*/

//
// Created by Lysogorskiy Yury on 13.12.23.
//

#ifdef PACE_TP
#ifdef PAIR_CLASS
// clang-format off
PairStyle(pace/tp_cell,PairPACETensorPotentialCell);
// clang-format on
#else

#ifndef LMP_PAIR_PACE_TP_CELL_H
#define LMP_PAIR_PACE_TP_CELL_H

#include "pair.h"
// CppFlow headers
#include <cppflow/ops.h>
#include <cppflow/model.h>
#include <cppflow/tensor.h>
#include <map>

#include "utils_pace.h"

namespace LAMMPS_NS {

class PairPACETensorPotentialCell : public Pair {
 public:
    PairPACETensorPotentialCell(class LAMMPS *);
  ~PairPACETensorPotentialCell() override;


  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  void *extract(const char *, int &) override;

 protected:
  struct ACETPImpl *aceimpl;
  virtual void allocate();

  double **scale;
  double cutoff = 6;

  int tot_neighbours = 0;

  int chunksize;
  double neigh_padding_fraction = 0.01;
  bool do_padding = true;

  int nelements;
  std::vector<std::string> elements_name;
  std::map<std::string, int> elements_to_index_map;
  std::vector<int> element_type_mapping; // LAMMPS's type to ACE's mu

  PACE::ACETimer data_timer;
  PACE::ACETimer tp_timer;
};
}    // namespace LAMMPS_NS

#endif
#endif
#endif //#ifdef