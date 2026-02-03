# Testing and Benchmarking Guide for GRACE Parallel Implementation

This guide documents the procedures for verifying correctness and measuring performance of the `pair_grace_2layer_parallel` implementation in LAMMPS.

## 1. Environment Requirements

Testing must be performed within the `grace` conda environment to ensure all dependencies (TF, cppflow) are correctly linked.

```bash
source ~/.bashrc
micromamba activate grace
module load mpi/openmpi/gcc/

#export TF_cpp_min_log_level=3
```

## 2. Correctness Verification (Serial vs. Parallel)

This test ensures that the parallel implementation (even on 1 rank) produces results numerically identical to the serial `pair_grace` style.

NOTE: Ensure that lmp is built before running tests.

### Automated Test
Run the `run_compare.sh` script:
```bash
cd performance_pg/L40s-1K
./run_compare.sh
```

### Manual Steps involved:
1. **Generate Rattled Structure**: Uses `generate_structure.py` to create `data.rattled_W_Ta` (ensures non-zero forces).
2. **Run Serial**: `lmp -in in.lammps.compare -var PAIRSTYLE grace -var SNAME serial`
3. **Run Parallel**:
    - 1 Rank: `mpirun -np 1 lmp -in in.lammps.compare -var PAIRSTYLE grace/2layer/parallel -var SNAME para1np`
    - 2 Ranks: `mpirun -np 2 lmp -in in.lammps.compare -var PAIRSTYLE grace/2layer/parallel -var SNAME para2np`
4. **Compare Forces**: Use `compare_forces.py` on the resulting dumps.
   - Threshold for **PASS**: Max Diff < 1e-5 (machine precision).

## 3. Performance Benchmarking

Measures strong scaling and compares parallel overhead against the serial baseline.

### Automated Benchmark
Run the `scaling_benchmark.sh` script:
```bash
cd performance_pg/L40s-1K
./scaling_benchmark.sh
```
NOTE: in order to use more GPUs, use `mpirun -np 4 --bind-to none bash -c 'CUDA_VISIBLE_DEVICES=$((OMPI_COMM_WORLD_RANK % 4)) /path/to/lmp ...`

### Interpreting Results (`benchmark_results.txt`):
- **Serial Baseline**: The target performance using a single fused TF call.
- **Parallel (Optimized) 1np**: Shows the overhead introduced by splitting the model into 3 TF signatures.
- **Scaling (2np, 4np)**: Measures MPI efficiency.

### Known Issues & Performance Expected:
- **Overhead**: Parallel models are expected to be ~3.6x slower on 1np due to multiple TF kernel launches AND dealing with ghost atoms (maybe way to optimization).
- **Memory (OOM)**: Systems larger than ~5000 atoms may hit **GPU RESOURCE_EXHAUSTED** during Layer 1 backward passes due to massive tensor allocations.
