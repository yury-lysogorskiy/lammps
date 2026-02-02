
import numpy as np
import tensorflow as tf
from ase.io import read
from tensorpotential.data.databuilder import GeometricalDataBuilder
from tensorpotential import constants
import os
from tensorpotential.calculator import TPCalculator

def tally_forces(grad_bv, ind_i, ind_j, n_atoms):
    f = np.zeros((n_atoms, 3))
    gbv = grad_bv.numpy()
    ii = ind_i.numpy()
    jj = ind_j.numpy()
    
    for k in range(len(gbv)):
        idx_i = ii[k]
        idx_j = jj[k]
        val = gbv[k]
        # F_i += val
        # F_j -= val
        if idx_i < n_atoms: f[idx_i] += val
        if idx_j < n_atoms: f[idx_j] -= val
    return -f


# 1. Setup paths
model_path = "/Users/lysogy36/dev/tensorpotential/pgs/GRACE-2L-OMAT-large-ft-E-parallel-comm"
# data_file = "/Users/lysogy36/dev/GITHUB/lammps/test/2L-parallel-comm/dimer.lammps-data"
data_file = "/Users/lysogy36/dev/GITHUB/lammps/test/2L-parallel-comm/chain.lammps-data"

calc = TPCalculator(model_path)

# 2. Load Model
print(f"Loading model from {model_path}...")
model = tf.saved_model.load(model_path)
print("Model loaded.")

# 3. Load Structure
print(f"Loading structure from {data_file}...")
atoms = read(data_file, format='lammps-data', style='atomic')
# ASE reads lammps-data atoms with type 1, 2. We need to map them to species.
# Based on log: Type 1 -> W, Type 2 -> Ta.
# We need to set atomic numbers or symbols for GeometricalDataBuilder to work if it uses them.
# Or if we manage the map manually.
# Let's set symbols for clarity if DataBuilder supports it.
# W (74), Ta (73).
species_map = {1: "W", 2: "Ta"}
symbols = [species_map[t] for t in atoms.get_atomic_numbers()] # atoms.numbers currently holds types 1,2 because style='atomic'
atoms.set_chemical_symbols(symbols)
print("Atoms:", atoms)
print("Positions:", atoms.get_positions())

atoms.calc=calc
f_grace_ase = atoms.get_forces()
print("Computed Atomic Forces from ASE:", f_grace_ase)
# 4. Prepare Data
# We need to know the 'elements' or 'element_map' the model expects.
# We can't easily query the element map from the SavedModel object directly unless stored in assets/metadata.
# However, for a generic GRACE model, we can try to guess or use the builder with known species.
# Assuming standard mapping logic or trying to extract from model if possible?
# Actually, GeometricalDataBuilder needs `elements_map`. W and Ta.
# Let's assume W=0, Ta=1 or similar? 
# Wait, the C++ log said:
# [GRACE] Mapping LAMMPS atom type #1(W) -> ACE species type #83
# [GRACE] Mapping LAMMPS atom type #2(Ta) -> ACE species type #73
# This implies the model expects actual atomic numbers as species types OR it has a map.
# If the inputs to the model are `atomic_mu_i` (atomic species), we need to provide what the model expects.
# Standard GRACE usually uses atomic numbers for `atomic_mu_i`.
# `utils_grace.cpp` maps types to elements.
# Let's try to pass atomic numbers 74 (W) and 73 (Ta) in `atomic_mu_i`.

# Configure DataBuilder
# We set rc to a safe large value to catch neighbors, though exact rc matters for results.
# The log said "master list distance cutoff = 6.5".
rcut = 6.5
# We can manually build the input dictionary if needed, but DataBuilder is safer.
# If we assume the model takes atomic numbers as species (typical for many TPs), we just need to ensure `elements_map` is effectively identity or handles symbols.
# Actually, `GeometricalDataBuilder` usually takes `elements_map` like {"W": 0, "Ta": 1} if it maps to dense integers 0..N.
# BUT, `atomic_mu_i` in some models is the *element type index* (0,1,2...) and in others it's the Z.
# Given the log "ACE species type #83", it sounds like it uses Z (Wait, Bi is 83. W is 74. Ta is 73. So maybe type 1 is Bi? 
# "Mapping LAMMPS atom type #1(W) -> ACE species type #83". 
# Wait, Bi is 83. W is 74. Is there a mismatch in the prompt log vs my periodic table knowledge?
# Or maybe the users potential has specific types.
# User log: "Mapping LAMMPS atom type #1(W) -> ACE species type #83"
# User log: "Mapping LAMMPS atom type #2(Ta) -> ACE species type #73"
# 73 is Ta. 74 is W. 83 is Bi. 
# Maybe type 1 is Bi? The input file said "mass 1 183.84 # W". Bi mass is 208.98. W is 183.84.
# So text says W, but ACE type says 83? That is confusing. 
# Ah, maybe the user explicitly mapped 1 to 83.
# The `pair_coeff` line was: `pair_coeff * * ... W Ta`.
# If the potential file maps "W" to 83 internally, maybe it is a W-Ta potential but trained with funny numbers?
# OR maybe I misread the log snippet.
# Step 190 log:
# [GRACE] Mapping LAMMPS atom type #1(W) -> ACE species type #83
# [GRACE] Mapping LAMMPS atom type #2(Ta) -> ACE species type #73
# OK, I will trust the log. 1->83, 2->73.
# In Python, we need to produce `atomic_mu_i` with these values.
# DataBuilder constructs `type_i` (0,1,...) usually. `atomic_mu_i` is then derived.
# Or maybe `atomic_mu_i` IS proper Z in the model.
# Let's check `GeometricalDataBuilder` defaults. 
# For now, let's manually patch `atomic_mu_i` to be 83/73 based on what we see.

db = GeometricalDataBuilder(
    elements_map={"W": 0, "Ta": 1}, 
    cutoff=rcut, 
    float_dtype="float64"
)
data = db.extract_from_ase_atoms(atoms)

# Patch atomic_mu_i and mu_i (neighs) to match the weird 83/73 if needed.
# Or check what DataBuilder produced.
# DataBuilder usually produces `atomic_mu_i` based on `elements_map` values? No, usually based on Z if `use_atomic_numbers` is set, or 0..N.
# Let's inspect what `atomic_mu_i` holds.
# If the model expects 83/73, we must ensure inputs match.
# In `pair_grace` C++, it constructs `element_type_mapping`. 
# "element_type_mapping[type[i]]" is used.
# If the model works with {83, 73}, we should start with that.
# Let's manually overwrite the species inputs to match the C++ execution exactly:
# Atom 1 (W) -> 83
# Atom 2 (W) -> 83
# Atom 3 (W) -> 83
# Atom 4 (Ta) -> 73
types = atoms.get_atomic_numbers() # These are 74(W) and 73(Ta)
z_map = {74: 83, 73: 73} 
correct_mu = np.array([z_map[t] for t in types], dtype=np.int32)
data["atomic_mu_i"] = correct_mu

# Also need to fix `mu_i` and `mu_j` in the neighbor list!
# `data["mu_i"]` and `data["mu_j"]` come from `ind_i` and `ind_j` indexing into `atomic_mu_i`.
# So if we fix `atomic_mu_i` and re-gather, it should be fine.
ind_i = data["ind_i"]
ind_j = data["ind_j"]
data["mu_i"] = tf.gather(data["atomic_mu_i"], ind_i)
data["mu_j"] = tf.gather(data["atomic_mu_i"], ind_j)

# Prepare input dict
# We wrap everything in Tensors with batch dim if needed, or just single element.
# The SavedModel expects batched inputs? The C++ code uses `cppflow::tensor` with shape {tot_atoms} etc, usually implying flattened or single batch.
# `TPCalculator`/DataBuilder usually adds a batch dim.
# But C++ code passes 1D arrays for `atomic_mu_i`.  Wait, C++ `tensor(..., {tot_atoms})` is rank 1.
# Python `GeometricalDataBuilder` might produce rank 1 or rank 2 (1, N).
# `graph_split_v2.py` line 115: `ds = tf.data.Dataset.from_tensors(data); inpt_tf = ds.get_single_element()`.
# This usually keeps the batch structure if data was batched.
# `extract_from_ase_atoms` returns unbatched dictionaries (rank 1 for arrays).
# `from_tensors` adds a batch dimension! So `inpt_tf` has shape (1, ...).
# However, the C++ code passes UNBATCHED tensors (rank 1 for lists).
# `pair_grace_2layer_parallel.cpp`: `add_input("atomic_mu_i", cppflow::tensor(aceimpl->atomic_mu_i, {aceimpl->tot_atoms}));` -> This is rank 1.
# If I use `tf.data.Dataset.from_tensors(data).get_single_element()`, I get rank 1 inputs inside? 
# No, `from_tensors` treats the whole dictionary as one element. `get_single_element` retrieves it. 
# But `GeometricalDataBuilder` output is unbatched.
# Wait, let's look at `graph_split_v2.py` again.
# `inpt_tf = ds.get_single_element()`
# If `data` is unbatched, `inpt_tf` is unbatched. `from_tensors` creates a dataset containing `t` as an element.
# So yes, it should be fine.

# Convert data to tensors
inpt_tf = {k: tf.convert_to_tensor(v) for k, v in data.items()}

def print_tensors_py(key, v, type_prefix="Input"):
    shape = v.shape
    shape_str = "[" + ", ".join(map(str, shape)) + "]"
    print(f"  {type_prefix}: {key} | Shape: {shape_str}")
    
    if len(shape) >= 1:
        n_atoms_to_print = min(int(shape[0]), 8)
        elements_per_atom = int(np.prod(shape[1:]))
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

# Helper to run signature
def run_signature(name, inputs):
    if "serving_default" not in name:
        # C++ only prints this for layer calls
        pass 
    else:
        print(f"\n--- Running {name} ---")

    sig = model.signatures[name]
    
    # Filter keys matching signature inputs
    # Inspect signature inputs
    # (Since we don't have the structured signature easily available without `structured_input_signature` checks,
    # we can try to pass all valid keys or inspect `sig.inputs`)
    # But `sig` is a ConcreteFunction. `sig.inputs` gives TensorSpecs? No, `sig.structured_input_signature`...
    # Let's just create a filtered dict of *Tensors*.
    # Using the logic from graph_split_v2 lines 226+
    valid_args = {}
    # We can try to match by name if possible, or just pass kwargs and hope TF ignores extras?
    # TF ConcreteFunction usually Errors on extra args.
    # We can rely on `structured_input_signature` if available.
    
    # Simple workaround: Try to match names from `sig.inputs` (tensors) to keys? 
    # The Tensor objects in `sig.inputs` have names like "atomic_mu_i:0".
    # The dictionary keys expected by `func(**kwargs)` match the *argument names* of the function.
    # For a SavedModel, these argument names are usually preserved if saved with signatures.
    
    # Using the helper from graph_split_v2
    func = sig
    if hasattr(func, "structured_input_signature") and func.structured_input_signature:
        spec = func.structured_input_signature[0] # args tuple? No, args is [0], kwargs is [1]?
        # Usually structured_input_signature is (args, kwargs).
        # We need to know what keys are expected.
        # If the model inputs are a flattened list of tensors (common in SavedModel serving), we might need to map them.
        # But `graph_split_v2.py` suggests we can pass a dict.
        pass
        
    # The C++ code passes inputs by name: "atomic_mu_i", etc.
    # So Python should be able to do `func(atomic_mu_i=..., ...)`
    
    # Let's aggressively try to pass strictly the keys we know are needed for each stage,
    # based on C++ implementation knowledge.
    
    required_keys = []
    if name == "forward_layer_1":
        required_keys = ["atomic_mu_i", "ind_i", "ind_j", "mu_i", "mu_j", "bond_vector", "batch_tot_nat_real"]
    elif name == "backward_layer_2":
        required_keys = ["atomic_mu_i", "ind_i", "ind_j", "mu_i", "mu_j", "bond_vector", "batch_tot_nat_real", "I", "I_nl_LN"]
    elif name == "backward_layer_1":
        required_keys = ["atomic_mu_i", "ind_i", "ind_j", "mu_i", "mu_j", "bond_vector", "batch_tot_nat_real", "grad_I", "grad_I_nl_LN"]
    elif name == "serving_default":
        # All basic geometrical inputs
        required_keys = ["atomic_mu_i", "ind_i", "ind_j", "mu_i", "mu_j", "bond_vector", "batch_tot_nat_real"]
    
    # Filter
    feed_dict = {k: inputs[k] for k in required_keys if k in inputs}
    
    if name != "serving_default":
        print(f"[GRACE-DEBUG] Calling model in {name}")
    
    for k, v in feed_dict.items():
        print_tensors_py(k, v, "Input")

    results = func(**feed_dict) if 'func' in locals() else sig(**feed_dict)
    
    # In SavedModel, results is a dict.
    # In C++, we get a vector and map names. 
    # Here we just iterate the dict.
    for k, v in results.items():
        print_tensors_py(k, v, "Output")
        
    return results

# 5. SERIAL RUN
print("\n=== SERIAL MODE (serving_default) ===")
serial_out = run_signature("serving_default", inpt_tf)
e_serial = serial_out["atomic_energy"]
f_serial = serial_out["z_pair_f"]
f_serial_per_atom = tally_forces(f_serial, inpt_tf["ind_i"], inpt_tf["ind_j"], inpt_tf["atomic_mu_i"].shape[0]) 
print(f"\nSerial Total Energy: {np.sum(e_serial)}")
print(f"Serial Pair Forces:\n{f_serial.numpy()}")
print(f"Serial Forces per atom:\n{f_serial_per_atom}")

# --- DEFINITIVE MASKING TEST ---
print(f"\n=== MASKING TEST: Simulate Ghosts in Python ===")
# We have 4 atoms. Let's make 7. 
# Ghosts 4, 5, 6 are images of 0, 1, 3.
# We don't bother with full neighbor consistency, just check if derivatives sum up.
def make_ghosts(d, n_real, n_ghost):
    out = {}
    for k, v in d.items():
        if k == "batch_tot_nat_real": 
            out[k] = tf.constant(n_real, dtype=tf.int32)
        elif len(v.shape) > 0 and v.shape[0] == n_real:
            # Atomic tensor
            ghost_part = v[:n_ghost]
            out[k] = tf.concat([v, ghost_part], axis=0)
        else:
            out[k] = v
    return out

# Prepare inputs for L1
inpt_f1_7 = make_ghosts(inpt_tf, 4, 3)
res_f1_7 = run_signature("forward_layer_1", inpt_f1_7)

# Prepare inputs for L2
inpt_l2_7 = inpt_f1_7.copy()
inpt_l2_7.update(res_f1_7)
res_l2_7 = run_signature("backward_layer_2", inpt_l2_7)

print("\n--- Masking Analysis ---")
gLN_7 = res_l2_7['grad_I_nl_LN'].numpy()
print(f"grad_I_nl_LN for first 4 atoms (7-atom config, real=4):")
print(gLN_7[:4, 0, 0])
print(f"grad_I_nl_LN for ghosts 4,5,6:")
print(gLN_7[4:, 0, 0])

# Reference from serial (already ran with 4 atoms)
res_f1_4 = run_signature("forward_layer_1", inpt_tf)
inpt_l2_4 = inpt_tf.copy()
inpt_l2_4.update(res_f1_4)
res_l2_4 = run_signature("backward_layer_2", inpt_l2_4)
gLN_4 = res_l2_4['grad_I_nl_LN'].numpy()
print(f"\nBaseline grad_I_nl_LN (4-atom config):")
print(gLN_4[:4, 0, 0])

# 6. PARALLEL RUN STEP-BY-STEP
print("\n=== PARALLEL MODE STEP-BY-STEP ===")

# Step A: Forward L1
l1_out = run_signature("forward_layer_1", inpt_tf)

# Prepare inputs for L2
inpt_l2 = inpt_tf.copy()
inpt_l2.update(l1_out)
# Note: "I" and "I_nl_LN" keys might need matching. 
# C++: I_KEY="I", I_LN_KEY="I_nl_LN". 
# Ensure `l1_out` has these keys. The signature outputs should name them correctly?
# C++ `out_names` loop checks signature output names.
# In Python, `results` dict keys are the output names.

# Step B: Backward L2
l2_out = run_signature("backward_layer_2", inpt_l2)

# Step C: Backward L1
inpt_l1_back = inpt_tf.copy()
# Map gradients: grad_I, grad_I_nl_LN
# l2_out should contain them.
inpt_l1_back.update(l2_out)

l1_back_out = run_signature("backward_layer_1", inpt_l1_back)

# Combine Forces
# Look for "grad_bond_vector"
grad_bv_l2 = l2_out["grad_bond_vector"]
grad_bv_l1 = l1_back_out["grad_bond_vector"]

# f_parallel = -(grad_bv_l2 + grad_bv_l1) # Note: Sign?
# In C++ we fixed it to be `+` accumulation into force array (which implies F = -grad).
# Wait, C++ `f[i] += fx`. `fx = scale * grad_bv`.
# If `grad_bv` is dE/dr_ij and F_i = -dE/dr_i.
# r_ij = r_j - r_i.  dE/dr_i = dE/dr_ij * d(r_ij)/dr_i = dE/dr_ij * (-1).
# So F_i = - (dE/dr_ij * -1) = + dE/dr_ij.
# So F_i should be positively proportional to `grad_bond_vector`.
# The standard output `z_pair_f` from `serving_default` is Forces.
# Let's check the sign relation between `z_pair_f` and `grad_bv`.
# Usually `z_pair_f` is computed internally as sum of grads.
# Let's compare `f_parallel` vs `f_serial`.

print("\n--- Force Comparison ---")
# Check if simple sum matches
f_parallel_sum = -(grad_bv_l2 + grad_bv_l1)
print("Sum (grad L2 + grad L1):")
print(f_parallel_sum.numpy()[:4]) # Just vectors matching neighbors?
# Wait, grad_bond_vector is per-bond (neighbor).
# We need to map back to atoms to get F_i.
# `pair_grace.cpp` does `f[i] += fx; f[j] -= fx;`.
# We can replicate this accumulation in Python using `ind_i`, `ind_j`.


# Serial forces are already packed as (N_atoms, 3).
# Parallel outputs are per-bond gradients.
n_atoms = len(atoms)
f_par_atom = tally_forces(f_parallel_sum, inpt_tf['ind_i'], inpt_tf['ind_j'], n_atoms)

print("Computed Atomic Forces from Parallel Grads:")
print(f_par_atom)

print("Serial Atomic Forces:")
print(f_serial_per_atom )

print("Computed Atomic Forces from ASE:", )
print(f_grace_ase )

diff0 = f_serial_per_atom - f_grace_ase
print ("Max diff serial vs ASE:", np.max(np.abs(diff0)))
# Comparison
diff = f_par_atom - f_grace_ase
print("Difference:")
print(diff)
print(f"Max Diff: {np.max(np.abs(diff))}")

