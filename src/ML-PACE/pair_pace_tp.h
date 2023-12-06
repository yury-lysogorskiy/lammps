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
// Created by Lysogorskiy Yury on 01.12.23.
//

#ifdef PACE_TP
#ifdef PAIR_CLASS
// clang-format off
PairStyle(pace/tp,PairPACETensorPotential);
// clang-format on
#else

#ifndef LMP_PAIR_PACE_TP_H
#define LMP_PAIR_PACE_TP_H

#include "pair.h"
// CppFlow headers
#include <cppflow/ops.h>
#include <cppflow/model.h>
#include <cppflow/tensor.h>
#include <map>
#include <chrono>

using namespace std::chrono;
using Clock = std::chrono::high_resolution_clock;
using TimePoint = std::chrono::time_point<Clock>;
using Duration = Clock::duration;

//////////////////////////////////////////
/**
 * Helper class for timing the code.
 * The timer should be initialized to reset measured time and
 * then call "start" and "stop" before and after measured code.
 * The measured time is stored in "duration" variable
 */
struct ACETimer {
    Duration duration; ///< measured duration
    TimePoint start_moment; ///< start moment of current measurement

    ACETimer() { init(); };

    /**
     * Reset timer
     */
    void init() { duration = std::chrono::nanoseconds(0); }

    /**
     * Start timer
     */
    void start() { start_moment = Clock::now(); }

    /**
     * Stop timer, update measured "duration"
     */
    void stop() { duration += Clock::now() - start_moment; }

    /**
     * Get duration in microseconds
     */
    long as_microseconds() { return std::chrono::duration_cast<std::chrono::microseconds>(duration).count(); }

    /**
     * Get duration in nanoseconds
     */
    long as_nanoseconds() { return std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count(); }

};

namespace LAMMPS_NS {

class PairPACETensorPotential : public Pair {
 public:
  PairPACETensorPotential(class LAMMPS *);
  ~PairPACETensorPotential() override;


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

  ACETimer data_timer;
  ACETimer tp_timer;
};
}    // namespace LAMMPS_NS

#endif
#endif
#endif //#ifdef PACE_TP