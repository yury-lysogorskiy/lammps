"""
Simulate LAMMPS parallel execution with ghost atoms in pure Python.

This script imitates the LAMMPS `grace/2layer/parallel` pair style using 
the split graph execution. Ghost atoms are represented as additional atoms
with indices beyond the local atom count.

Comparisons:
1. serving_default (no ghosts) vs TPCalculator
2. Parallel graph with ghosts vs TPCalculator
"""

import numpy as np
import tensorflow as tf
from ase.io import read
from tensorpotential.data.databuilder import GeometricalDataBuilder
from tensorpotential.calculator import TPCalculator

# --- Configuration ---
model_path = "/Users/lysogy36/dev/tensorpotential/pgs/GRACE-2L-OMAT-large-ft-E-parallel-comm"
data_file = "/Users/lysogy36/dev/GITHUB/lammps/test/2L-parallel-comm/chain.lammps-data"
rcut = 6.5

# --- Helper Functions ---
def tally_forces(grad_bv, ind_i, ind_j, n_local, ghost_parent_map=None):
    """Tally pair forces from bond gradients to per-atom forces.
    
    If ghost_parent_map is provided, ghost neighbor indices are mapped
    to their parent local atoms for force accumulation.
    """
    f = np.zeros((n_local, 3))
    gbv = grad_bv.numpy() if hasattr(grad_bv, 'numpy') else grad_bv
    ii = ind_i.numpy() if hasattr(ind_i, 'numpy') else ind_i
    jj = ind_j.numpy() if hasattr(ind_j, 'numpy') else ind_j

    for k in range(len(gbv)):
        idx_i = ii[k]
        idx_j = jj[k]
        val = gbv[k]
        
        # Map ghost indices to parent atoms
        if ghost_parent_map is not None:
            if idx_i in ghost_parent_map:
                idx_i = ghost_parent_map[idx_i]
            if idx_j in ghost_parent_map:
                idx_j = ghost_parent_map[idx_j]
        
        if idx_i < n_local:
            f[idx_i] += val
        if idx_j < n_local:
            f[idx_j] -= val
    return -f



def print_tensors_py(key, v, type_prefix="Input"):
    """Print tensor info in atom-wise format (matching C++ output)."""
    shape = v.shape
    shape_str = "[" + ", ".join(map(str, shape)) + "]"
    print(f"  {type_prefix}: {key} | Shape: {shape_str}")
    
    if len(shape) >= 1:
        n_atoms_to_print = min(int(shape[0]), 8)
        elements_per_atom = int(np.prod(shape[1:])) if len(shape) > 1 else 1
        n_elements_to_print = min(elements_per_atom, 8)
        
        flat = tf.reshape(v, [int(shape[0]), -1]).numpy()
        for i in range(n_atoms_to_print):
            print(f"    Atom {i}:", end="")
            for k in range(n_elements_to_print):
                val = flat[i, k]
                if v.dtype == tf.int32:
                    print(f" {int(val)}", end="")
                else:
                    print(f" {val:8.4f}", end="")
            print()
    else:
        # Scalar
        val = v.numpy()
        if v.dtype == tf.int32:
            print(f"    Value: {int(val)}")
        else:
            print(f"    Value: {val:8.4f}")


def print_array(name, arr, max_rows=8):
    """Print numpy array in numpy-like format."""
    print(f"{name}:")
    arr_np = arr.numpy() if hasattr(arr, 'numpy') else arr
    for i, row in enumerate(arr_np[:max_rows]):
        prefix = "[" if i == 0 else " "
        suffix = "]" if i == min(len(arr_np), max_rows) - 1 else ""
        print(f"{prefix}[{row[0]:12.8f} {row[1]:12.8f} {row[2]:12.8f}]{suffix}")


# --- Load Model and Structure ---
print(f"Loading model from {model_path}...")
model = tf.saved_model.load(model_path)
calc = TPCalculator(model_path)

print(f"Loading structure from {data_file}...")
atoms = read(data_file, format='lammps-data', style='atomic')
symbols = ["W"] * len(atoms)
atoms.set_chemical_symbols(symbols)
n_local = len(atoms)
print(f"Atoms: {atoms}")
print(f"Positions:\n{atoms.get_positions()}")

# --- Reference: TPCalculator Forces ---
atoms.calc = calc
f_ref = atoms.get_forces()
print("\n" + "="*60)
print("=== REFERENCE: TPCalculator Forces ===")
print("="*60)
print_array("Forces (TPCalculator)", f_ref)

# --- Prepare Base Data ---
db = GeometricalDataBuilder(
    elements_map={"W": 0, "Ta": 1},
    cutoff=rcut,
    float_dtype="float64"
)
data = db.extract_from_ase_atoms(atoms)

# Patch atomic_mu_i to match LAMMPS C++ mapping (W->83)
correct_mu = np.array([83] * n_local, dtype=np.int32)
data["atomic_mu_i"] = correct_mu

# Re-gather mu_i/mu_j
ind_i = data["ind_i"]
ind_j = data["ind_j"]
data["mu_i"] = tf.gather(data["atomic_mu_i"], ind_i)
data["mu_j"] = tf.gather(data["atomic_mu_i"], ind_j)

inpt_base = {k: tf.convert_to_tensor(v) for k, v in data.items()}


# --- Run Signature with Tensor Printing ---
def run_signature(name, inputs, required_keys):
    """Run a model signature and print input/output tensors."""
    print(f"\n[GRACE-DEBUG] Calling model in {name}")
    
    sig = model.signatures[name]
    feed = {k: inputs[k] for k in required_keys if k in inputs}
    
    # Print inputs
    for k in required_keys:
        if k in inputs:
            print_tensors_py(k, inputs[k], "Input")
    
    # Run
    results = sig(**feed)
    
    # Print outputs
    for k, v in results.items():
        print_tensors_py(k, v, "Output")
    
    return results


# =========================================================================
# PART 1: SERIAL RUN (serving_default, no ghosts)
# =========================================================================
print("\n" + "="*60)
print("=== PART 1: SERIAL RUN (serving_default) ===")
print("="*60)

serial_keys = ["atomic_mu_i", "ind_i", "ind_j", "mu_i", "mu_j", "bond_vector", "batch_tot_nat_real"]
serial_out = run_signature("serving_default", inpt_base, serial_keys)

e_serial = serial_out["atomic_energy"]
f_serial = serial_out["z_pair_f"]
f_serial_per_atom = tally_forces(f_serial, inpt_base["ind_i"], inpt_base["ind_j"], n_local)

print(f"\nSerial Total Energy: {np.sum(e_serial.numpy()):.6f}")
print_array("Serial Pair Forces (z_pair_f)", f_serial)
print_array("Serial Atomic Forces (tallied)", f_serial_per_atom)

print("\n--- Comparison: serving_default vs TPCalculator ---")
diff_serial = f_serial_per_atom - f_ref
print(f"Max Abs Diff: {np.max(np.abs(diff_serial)):.2e}")
if np.max(np.abs(diff_serial)) < 1e-10:
    print("✓ serving_default matches TPCalculator")
else:
    print("✗ serving_default does NOT match TPCalculator")


# =========================================================================
# PART 2: PARALLEL RUN (step-by-step, no ghosts)
# =========================================================================
print("\n" + "="*60)
print("=== PART 2: PARALLEL RUN (step-by-step, no ghosts) ===")
print("="*60)

# Forward Layer 1
l1_keys = ["atomic_mu_i", "ind_i", "ind_j", "mu_i", "mu_j", "bond_vector", "batch_tot_nat_real"]
l1_out = run_signature("forward_layer_1", inpt_base, l1_keys)

# Backward Layer 2
inpt_l2 = {**inpt_base, **l1_out}
l2_keys = ["atomic_mu_i", "ind_i", "ind_j", "mu_i", "mu_j", "bond_vector", "batch_tot_nat_real", "I", "I_nl_LN"]
l2_out = run_signature("backward_layer_2", inpt_l2, l2_keys)

# Backward Layer 1
inpt_l1_back = {**inpt_base, **l2_out}
l1_keys_back = ["atomic_mu_i", "ind_i", "ind_j", "mu_i", "mu_j", "bond_vector", "batch_tot_nat_real", "grad_I", "grad_I_nl_LN"]
l1_back_out = run_signature("backward_layer_1", inpt_l1_back, l1_keys_back)

# Combine Forces
grad_bv_l2 = l2_out["grad_bond_vector"].numpy()
grad_bv_l1 = l1_back_out["grad_bond_vector"].numpy()
grad_bv_total = grad_bv_l2 + grad_bv_l1

print("\n--- Force Comparison (Parallel no-ghost) ---")
print_array("Sum (grad L2 + grad L1)", grad_bv_total)

# Note: grad_bond_vector = -z_pair_f, so negate to match the serving_default output
f_par = tally_forces(-grad_bv_total, inpt_base["ind_i"], inpt_base["ind_j"], n_local)
print_array("Computed Atomic Forces from Parallel Grads", f_par)

print_array("Reference Forces (TPCalculator)", f_ref)

diff_par = f_par - f_ref
print(f"\nMax Abs Diff: {np.max(np.abs(diff_par)):.2e}")
if np.max(np.abs(diff_par)) < 1e-10:
    print("✓ Parallel (no ghost) matches TPCalculator")
else:
    print("✗ Parallel (no ghost) does NOT match TPCalculator")


# =========================================================================
# PART 3: PARALLEL RUN WITH GHOSTS (simulating LAMMPS PBC + comm)
# =========================================================================
print("\n" + "="*60)
print("=== PART 3: PARALLEL RUN WITH GHOSTS ===")
print("="*60)

# From LAMMPS debug output: 4 local + 3 ghosts = 7 atoms
# Ghost mapping (from LAMMPS neighbor list analysis):
#   ghost 4 is image of atom 0 (seen by atom 3 at z+3.1=20)
#   ghost 5 is image of atom 1 (not in neighbor list, but exists)
#   ghost 6 is image of atom 3 (seen by atom 0 at z-3.1=-3.1)
ghost_parent_map = {4: 0, 5: 1, 6: 3}  # ghost_index -> local_parent_index

n_ghost = 3
n_all = n_local + n_ghost

all_mu = np.array([83] * n_all, dtype=np.int32)

# Exact LAMMPS neighbor list (from debug_output)
new_ind_i = np.array([0, 0, 1, 1, 2, 2, 3, 3], dtype=np.int32)
new_ind_j = np.array([6, 1, 0, 2, 1, 3, 2, 4], dtype=np.int32)
new_bond_vectors = np.array([
    [0.0, 0.0, -3.1],
    [0.0, 0.0, 5.95],
    [0.0, 0.0, -5.95],
    [0.0, 0.0, 5.95],
    [0.0, 0.0, -5.95],
    [0.0, 0.0, 5.0],
    [0.0, 0.0, -5.0],
    [0.0, 0.0, 3.1],
], dtype=np.float64)

inpt_ghost = {
    "atomic_mu_i": tf.constant(all_mu, dtype=tf.int32),
    "ind_i": tf.constant(new_ind_i, dtype=tf.int32),
    "ind_j": tf.constant(new_ind_j, dtype=tf.int32),
    "mu_i": tf.gather(all_mu, new_ind_i),
    "mu_j": tf.gather(all_mu, new_ind_j),
    "bond_vector": tf.constant(new_bond_vectors, dtype=tf.float64),
    "batch_tot_nat_real": tf.constant(n_local, dtype=tf.int32),
}

print(f"Local atoms: {n_local}, Ghost atoms: {n_ghost}, Total: {n_all}")
print(f"Ghost parent map: {ghost_parent_map}")
print(f"Bonds: {len(new_ind_i)}")

# Forward Layer 1
l1_out_g = run_signature("forward_layer_1", inpt_ghost, l1_keys)

# ============ SIMULATE FORWARD_COMM ============
# Copy features from local atoms to their ghost copies
print("\n[GRACE-DEBUG] Simulating forward_comm: copying local features to ghosts")
I_np = l1_out_g["I"].numpy().copy()
I_nl_LN_np = l1_out_g["I_nl_LN"].numpy().copy()

print(f"  Before forward_comm:")
print(f"    I[0, 0, :4] = {I_np[0, 0, :4]}")
print(f"    I[4, 0, :4] = {I_np[4, 0, :4]} (ghost of 0)")
print(f"    I[3, 0, :4] = {I_np[3, 0, :4]}")
print(f"    I[6, 0, :4] = {I_np[6, 0, :4]} (ghost of 3)")

for ghost_idx, parent_idx in ghost_parent_map.items():
    I_np[ghost_idx] = I_np[parent_idx]
    I_nl_LN_np[ghost_idx] = I_nl_LN_np[parent_idx]

print(f"  After forward_comm:")
print(f"    I[4, 0, :4] = {I_np[4, 0, :4]} (copied from 0)")
print(f"    I[6, 0, :4] = {I_np[6, 0, :4]} (copied from 3)")

# Update l1_out_g with communicated features
l1_out_g_comm = {
    "I": tf.constant(I_np, dtype=tf.float64),
    "I_nl_LN": tf.constant(I_nl_LN_np, dtype=tf.float64),
}

# Backward Layer 2 (with communicated features)
inpt_l2_g = {**inpt_ghost, **l1_out_g_comm}
l2_out_g = run_signature("backward_layer_2", inpt_l2_g, l2_keys)

# ============ SIMULATE REVERSE_COMM ============
# Correct behavior: 
# - grad_I: accumulate ghost gradients to parent atoms (feature cross-terms)
# - grad_I_nl_LN: just zero ghosts, do NOT accumulate (loss contribution mask)
print("\n[GRACE-DEBUG] Simulating reverse_comm:")
grad_I = l2_out_g["grad_I"].numpy().copy()
grad_I_nl_LN = l2_out_g["grad_I_nl_LN"].numpy().copy()

print(f"  Before reverse_comm:")
print(f"    grad_I[0, 0, :4] = {grad_I[0, 0, :4]}")
print(f"    grad_I[4, 0, :4] = {grad_I[4, 0, :4]} (ghost of 0)")
print(f"    grad_I[3, 0, :4] = {grad_I[3, 0, :4]}")
print(f"    grad_I[6, 0, :4] = {grad_I[6, 0, :4]} (ghost of 3)")
print(f"    grad_I_nl_LN[0:4, 0, 0] = {grad_I_nl_LN[:4, 0, 0]}")
print(f"    grad_I_nl_LN[4:7, 0, 0] = {grad_I_nl_LN[4:, 0, 0]} (ghosts)")

# Accumulate grad_I from ghosts to parents, then zero ghosts
for ghost_idx, parent_idx in ghost_parent_map.items():
    grad_I[parent_idx] += grad_I[ghost_idx]
    grad_I[ghost_idx] = 0.0
    # Just zero ghost grad_I_nl_LN - do NOT accumulate (it's a loss mask, not a gradient sum)
    grad_I_nl_LN[ghost_idx] = 0.0

print(f"  After reverse_comm:")
print(f"    grad_I[0, 0, :4] = {grad_I[0, 0, :4]} (accumulated from ghost 4)")
print(f"    grad_I[4, 0, :4] = {grad_I[4, 0, :4]} (zeroed)")
print(f"    grad_I[3, 0, :4] = {grad_I[3, 0, :4]} (accumulated from ghost 6)")
print(f"    grad_I[6, 0, :4] = {grad_I[6, 0, :4]} (zeroed)")
print(f"    grad_I_nl_LN[0:4, 0, 0] = {grad_I_nl_LN[:4, 0, 0]} (unchanged)")
print(f"    grad_I_nl_LN[4:7, 0, 0] = {grad_I_nl_LN[4:, 0, 0]} (zeroed)")




# Backward Layer 1 (with communicated gradients)
inpt_l1_back_g = {
    **inpt_ghost,
    "grad_I": tf.constant(grad_I, dtype=tf.float64),
    "grad_I_nl_LN": tf.constant(grad_I_nl_LN, dtype=tf.float64),
}
l1_back_out_g = run_signature("backward_layer_1", inpt_l1_back_g, l1_keys_back)

# Combine Forces
grad_bv_l2_g = l2_out_g["grad_bond_vector"].numpy()
grad_bv_l1_g = l1_back_out_g["grad_bond_vector"].numpy()
grad_bv_total_g = grad_bv_l2_g + grad_bv_l1_g

print("\n--- Force Comparison (Parallel with ghosts + comm) ---")
print_array("Sum (grad L2 + grad L1)", grad_bv_total_g)

# Note: grad_bond_vector = -z_pair_f, so negate to match
# Use ghost_parent_map to attribute forces from ghost bonds to parent atoms
f_par_g = tally_forces(-grad_bv_total_g, new_ind_i, new_ind_j, n_local, ghost_parent_map)
print_array("Computed Atomic Forces from Parallel Grads", f_par_g)

print_array("Reference Forces (TPCalculator)", f_ref)

diff_par_g = f_par_g - f_ref
print(f"\nDifference (Parallel with ghosts - Reference):")
print(diff_par_g)
print(f"Max Abs Diff: {np.max(np.abs(diff_par_g)):.2e}")

if np.max(np.abs(diff_par_g)) < 1e-10:
    print("\n✓ SUCCESS: Parallel with ghosts + comm matches TPCalculator!")
else:
    print("\n✗ MISMATCH: Parallel with ghosts does NOT match TPCalculator. Further investigation needed.")

