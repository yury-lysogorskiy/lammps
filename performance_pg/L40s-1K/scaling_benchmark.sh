#!/bin/bash
# Strong Scaling Benchmark for Optimized GRACE
source ~/.bashrc
micromamba activate grace
module load mpi/openmpi/gcc/
export TF_cpp_min_log_level=3

LMP="../../build-parallel-comm/lmp"
POT_SERIAL="/home/users/lysogy36/.cache/grace/GRACE-2L-OMAT-large-ft-E"
POT_PARA="/home/users/lysogy36/.cache/grace/GRACE-2L-OMAT-large-ft-E-parallel-comm"

RESULTS="benchmark_results.txt"
echo "=== GRACE Optimization Benchmark - $(date) ===" > $RESULTS
echo "" >> $RESULTS

run_test() {
    local np=$1
    local style=$2
    local replicate=$3
    local label=$4
    
    local atoms=$((2 * replicate * replicate * replicate))
    local pot=$POT_SERIAL
    if [[ "$style" == *"parallel"* ]]; then
        pot=$POT_PARA
    fi
    
    echo "Running: $label (np=$np, style=$style, replicate=$replicate, atoms=$atoms)" | tee -a $RESULTS
    
    # Create temporary input file
    cat > temp_bench.in << EOF
units metal
atom_style atomic
boundary p p p

lattice bcc 3.165
region box block 0 1 0 1 0 1
create_box 1 box
create_atoms 1 box
replicate $replicate $replicate $replicate

pair_style $style
pair_coeff * * $pot W

mass 1 183.84
velocity all create 300.0 12345 mom yes rot yes dist gaussian

timestep 0.001
thermo 100
# Run enough steps to get stable timing
fix 1 all nvt temp 300.0 300.0 0.1
run 20
EOF

    # Run with timing
    cmd=""
    if [ "$np" -eq 1 ]; then
        cmd="$LMP -in temp_bench.in"
    else
        cmd="mpirun -np $np --bind-to none bash -c 'CUDA_VISIBLE_DEVICES=$((OMPI_COMM_WORLD_RANK % 4)) $LMP -in temp_bench.in'"
    fi
    
    # Capture full log
    logname="log.${label// /_}"
    $cmd > $logname 2>&1
    
    # Grep to results
    echo "--- Output tail ---" >> $RESULTS
    tail -20 $logname >> $RESULTS
    grep -E "\[GRACE-PROFILE\]|Loop time|atoms" $logname >> $RESULTS
    echo "" >> $RESULTS
}

echo "========================================" >> $RESULTS
echo "CASE 1: 1024 ATOMS (8x8x8)" >> $RESULTS
echo "========================================" >> $RESULTS

run_test 1 "grace" 8 "Serial Baseline"
run_test 1 "grace/2layer/parallel" 8 "Parallel (Optimized) 1np"
run_test 2 "grace/2layer/parallel" 8 "Parallel (Optimized) 2np"

echo "========================================" >> $RESULTS
echo "CASE 2: 16000 ATOMS (20x20x20)" >> $RESULTS
echo "========================================" >> $RESULTS

run_test 1 "grace" 20 "Serial Baseline"
run_test 1 "grace/2layer/parallel" 20 "Parallel (Optimized) 1np"
run_test 2 "grace/2layer/parallel" 20 "Parallel (Optimized) 2np"
run_test 4 "grace/2layer/parallel" 20 "Parallel (Optimized) 4np"

rm -f temp_bench.in
cat $RESULTS
