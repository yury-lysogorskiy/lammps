
import numpy as np
from ase.build import bulk
from ase.io import write
import random

# Parameters
structure_file = "data.rattled_W_Ta"
lattice_constant = 3.2  # Approx for W-Ta
rattle_amplitude = 0.1
elements = ['W', 'Ta']

# Create BCC W structure
atoms = bulk('W', 'bcc', a=lattice_constant, cubic=True)
atoms = atoms * (2, 2, 2)

# Randomly assign elements
symbols = []
for i in range(len(atoms)):
    # 50/50 mix roughly
    symbols.append(random.choice(elements))
atoms.set_chemical_symbols(symbols)

# Rattle
atoms.rattle(stdev=rattle_amplitude, seed=42)

# Write to LAMMPS data file
# We need to ensure we map species to types correctly: W->1, Ta->2
# ASE writes types based on order in symbols or spec.
# Let's handle generic writing
write(structure_file, atoms, format='lammps-data', atom_style='atomic')

print(f"Generated {structure_file} with {len(atoms)} atoms.")
print(f"Elements: {atoms.get_chemical_formula()}")
