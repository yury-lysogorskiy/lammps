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

#endif