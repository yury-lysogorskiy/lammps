
import sys
import numpy as np

def read_dump(filename):
    """
    Reads a LAMMPS dump file and returns a dictionary of {id: [fx, fy, fz]}.
    Assumes format: ITEM: ATOMS id type fx fy fz
    """
    atoms = {}
    with open(filename, 'r') as f:
        lines = f.readlines()
        
    start_idx = -1
    for i, line in enumerate(lines):
        if line.startswith("ITEM: ATOMS"):
            start_idx = i + 1
            break
            
    if start_idx == -1:
        raise ValueError(f"Could not find atom data in {filename}")
        
    for line in lines[start_idx:]:
        parts = line.split()
        if len(parts) >= 5:
            atom_id = int(parts[0])
            f_vec = [float(x) for x in parts[2:5]]
            atoms[atom_id] = np.array(f_vec)
            
    return atoms

def compare_forces(file1, file2, tolerance=1e-6):
    atoms1 = read_dump(file1)
    atoms2 = read_dump(file2)
    
    ids1 = set(atoms1.keys())
    ids2 = set(atoms2.keys())
    
    if ids1 != ids2:
        print(f"FAIL: Atom ID mistmatch. Unique in 1: {ids1-ids2}, Unique in 2: {ids2-ids1}")
        return False

    max_diff = 0.0
    sum_sq_diff = 0.0
    count = 0
    
    for aid in ids1:
        f1 = atoms1[aid]
        f2 = atoms2[aid]
        diff = np.linalg.norm(f1 - f2)
        if diff > max_diff:
            max_diff = diff
        sum_sq_diff += diff**2
        count += 1
        
    rmse = np.sqrt(sum_sq_diff / count) if count > 0 else 0.0
    
    print(f"Comparison {file1} vs {file2}:")
    print(f"  Max Diff: {max_diff:.6e}")
    print(f"  RMSE:     {rmse:.6e}")
    
    if max_diff > tolerance:
        print("  RESULT: FAIL (exceeds tolerance)")
        return False
    else:
        print("  RESULT: PASS")
        return True

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python compare_forces.py file1 file2 [tol]")
        sys.exit(1)
        
    f1 = sys.argv[1]
    f2 = sys.argv[2]
    tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-6
    
    compare_forces(f1, f2, tol)
