#!/bin/bash
source ~/.bashrc
micromamba activate grace
module load mpi/openmpi/gcc/
export TF_cpp_min_log_level=3
LMP=../../build-parallel-comm/lmp

echo "--- Generating Structure ---"
python generate_structure.py

echo "--- Running Serial (grace) ---"
$LMP -in in.lammps.compare -var PAIRSTYLE grace -var SNAME serial > log.serial

echo "--- Running Parallel (1 rank) ---"
mpirun -np 1 $LMP -in in.lammps.compare -var PAIRSTYLE grace/2layer/parallel -var SNAME para1np > log.para1np

echo "--- Running Parallel (2 ranks) ---"
mpirun -np 2 $LMP -in in.lammps.compare -var PAIRSTYLE grace/2layer/parallel -var SNAME para2np > log.para2np

echo "--- Comparing Forces ---"
python compare_forces.py forces.serial.dump forces.para1np.dump 1e-5
echo "Diff Serial vs Para2np:"
python compare_forces.py forces.serial.dump forces.para2np.dump 1e-5

echo "--- Energies ---"
grep -B 1 "Loop time of" log.serial
grep -B 1 "Loop time of" log.para1np
grep -B 1 "Loop time of" log.para2np
