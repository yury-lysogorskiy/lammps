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
PairStyle(grace/1layer/chunk,PairGRACE1LayerChunk);
// clang-format on
#else

#ifndef LMP_PAIR_GRACE_1LAYER_CHUNK_H
#define LMP_PAIR_GRACE_1LAYER_CHUNK_H

#include "pair.h"
#include "utils_pace.h"
#include <map>
#include <vector>

namespace LAMMPS_NS {

class PairGRACE1LayerChunk : public Pair {
 public:
  PairGRACE1LayerChunk(class LAMMPS *);
  ~PairGRACE1LayerChunk() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

 protected:
  struct GRACE1LayerChunkImpl *graceimpl;

  // Model and signature info
  std::string compute_function_name = "serving_default";

  // Cutoffs and scaling
  virtual void allocate();
  double **scale;
  double cutoff = 6.0;
  bool is_custom_cutoffs = false;
  std::vector<std::vector<double>> cutoff_matrix, cutoff_matrix_per_lammps_type;
  bool pair_forces = false;

  // Chunking parameters
  int chunksize = 4096;

  // Element mapping
  std::vector<std::string> elements_name;
  std::map<std::string, int> elements_to_index_map;
  std::vector<int> element_type_mapping;    // LAMMPS type -> ACE species

  // Timers
  PACE::ACETimer total_timer;
  PACE::ACETimer data_timer;
  PACE::ACETimer tp_timer;
};

}    // namespace LAMMPS_NS

#endif
#endif
#endif    //#ifndef NO_GRACE_TF
