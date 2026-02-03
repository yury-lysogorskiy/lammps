#!/bin/bash
# Strong and Weak Scaling Tests for GRACE potentials
# Machine: 4x L40s GPUs

LMP="../../build-parallel-comm/lmp"
POTENTIAL="/home/users/lysogy36/.cache/grace/GRACE-2L-OMAT-large-ft-E-parallel-comm"

# Create result file
RESULTS="scaling_results.txt"
echo "=== GRACE Scaling Tests - $(date) ===" > $RESULTS
echo "" >> $RESULTS

run_test() {
    local np=$1
    local style=$2
    local replicate=$3
    local label=$4
    
    # Calculate expected atoms
    local atoms=$((2 * replicate * replicate * replicate))
    
    echo "Running: $label (np=$np, style=$style, replicate=$replicate, atoms=$atoms)" | tee -a $RESULTS
    
    # Create temporary input file
    cat > temp_scaling.in << EOF
units metal
atom_style atomic
boundary p p p

lattice bcc 3.165
region box block 0 1 0 1 0 1
create_box 1 box
create_atoms 1 box

replicate $replicate $replicate $replicate

pair_style $style
pair_coeff * * $POTENTIAL W

mass 1 183.84
velocity all create 300.0 12345 mom yes rot yes dist gaussian

timestep 0.001
thermo 50

fix 1 all nvt temp 300.0 300.0 0.1
run 3
EOF

    # Run with timing
    if [ "$np" -eq 1 ]; then
        $LMP -in temp_scaling.in 2>&1 | grep -E "\[GRACE-PROFILE\]|Loop time|atoms" | tail -10 >> $RESULTS
    else
        mpirun -np $np $LMP -in temp_scaling.in 2>&1 | grep -E "\[GRACE-PROFILE\]|Loop time|atoms" | tail -10 >> $RESULTS
    fi
    echo "" >> $RESULTS
}

echo "========================================" >> $RESULTS
echo "STRONG SCALING (Fixed 1024 atoms)" >> $RESULTS
echo "========================================" >> $RESULTS

# Strong scaling: 1024 atoms (8x8x8 replicate)
echo "--- pair_grace (serial) ---" >> $RESULTS
run_test 1 "grace" 8 "grace-1np-1024at"

echo "--- pair_grace/2layer/parallel ---" >> $RESULTS
run_test 1 "grace/2layer/parallel" 8 "parallel-1np-1024at"
run_test 2 "grace/2layer/parallel" 8 "parallel-2np-1024at"
run_test 4 "grace/2layer/parallel" 8 "parallel-4np-1024at"

echo "========================================" >> $RESULTS
echo "STRONG SCALING (Fixed 8000 atoms)" >> $RESULTS
echo "========================================" >> $RESULTS

# Strong scaling: 8000 atoms (20x20x20 replicate)
echo "--- pair_grace (serial) ---" >> $RESULTS
run_test 1 "grace" 20 "grace-1np-8000at"

echo "--- pair_grace/2layer/parallel ---" >> $RESULTS
run_test 1 "grace/2layer/parallel" 20 "parallel-1np-8000at"
run_test 2 "grace/2layer/parallel" 20 "parallel-2np-8000at"
run_test 4 "grace/2layer/parallel" 20 "parallel-4np-8000at"

echo "========================================" >> $RESULTS
echo "WEAK SCALING (1024 atoms per rank)" >> $RESULTS
echo "========================================" >> $RESULTS

# Weak scaling: 1024 atoms per rank
echo "--- pair_grace/2layer/parallel ---" >> $RESULTS
run_test 1 "grace/2layer/parallel" 8 "parallel-1np-1024at-weak"   # 1024 atoms
run_test 2 "grace/2layer/parallel" 10 "parallel-2np-2000at-weak"  # 2000 atoms
run_test 4 "grace/2layer/parallel" 13 "parallel-4np-4394at-weak"  # ~4394 atoms

# Cleanup
rm -f temp_scaling.in

echo "========================================" >> $RESULTS
echo "TESTS COMPLETE" >> $RESULTS
echo "========================================" >> $RESULTS

cat $RESULTS
