//
// Created by Yury Lysogorskiy on 30.12.25.
//

#ifndef LAMMPS_GRACE_UTILS_H
#define LAMMPS_GRACE_UTILS_H

#include <vector>
#include <cppflow/tensor.h>
#include "lammps.h"

void print_tf_inputs(const std::vector<std::tuple<std::string, cppflow::tensor>>& inputs,
                            int me, LAMMPS_NS::LAMMPS *lmp, bool python_ready = false);

void print_f_data(const double *f_data, int me, LAMMPS_NS::LAMMPS *lmp,
  const std::vector<int>& ind_i_vector, const std::vector<int>& ind_j_vector,
  int * tag);

#endif