import numpy as np
import re
import sys
import os
from scipy.spatial import Delaunay, ConvexHull


def parse_qtaim_files(file_paths):
    cps = {'nuclei': {}, 'ccps': {}, 'rho_ccp': {}, 'laplacian_ccp': {}}
    ccp_indices_ordered = []
    rho_ccp_list = []
    pending_laplacian = {}
    laplacian_found = False
    all_lines = []
    full_content = ""

    for fp in file_paths:
        try:
            with open(fp, 'r') as f:
                text = f.read()
                full_content += text + "\n"
                all_lines.extend(text.splitlines())
        except FileNotFoundError:
            print(f"Error: Could not find file '{fp}'. Skipping.")

    rho_match = re.search(r'\(3,\+3\) = \{(.*?)\}', full_content, re.DOTALL)
    if rho_match:
        inner_text = rho_match.group(1).strip()
        if inner_text:
            rho_ccp_list = [float(x) for x in inner_text.split(',') if x.strip()]

    current_cp_idx = None

    for line in all_lines:
        if not line.strip():
            continue

        if "CP" in line and "Type (3," in line:
            match = re.search(r'CP\s+(\d+)', line, re.IGNORECASE)
            if match:
                current_cp_idx = int(match.group(1))
        elif re.match(r'^\s*\d+', line) and "(3," in line:
            try:
                current_cp_idx = int(line.split()[0])
            except Exception:
                pass

        float_matches = re.findall(r'-?\d+\.\d+(?:[eE][+-]?\d+)?', line)
        coords = (np.array([float(float_matches[0]),
                            float(float_matches[1]),
                            float(float_matches[2])])
                  if len(float_matches) >= 3 else np.zeros(3))

        if "(3,-3)" in line:
            atom_match = re.search(r'(\d+)\s*\(\s*([A-Za-z]+)\s*\)', line)
            if atom_match:
                atom_num = int(atom_match.group(1))
                cps['nuclei'][atom_num] = {'coords': coords,
                                           'symbol': atom_match.group(2)}
            elif "Nucleus:" in line:
                nuc_match = re.search(r'Nucleus:\s*([A-Za-z]+)\s*(\d+)', line)
                if nuc_match:
                    atom_num = int(nuc_match.group(2))
                    if atom_num not in cps['nuclei']:
                        cps['nuclei'][atom_num] = {'coords': coords,
                                                   'symbol': nuc_match.group(1)}

        elif "(3,+3)" in line and len(float_matches) >= 3:
            if current_cp_idx is not None:
                cps['ccps'][current_cp_idx] = coords
                if current_cp_idx not in ccp_indices_ordered:
                    ccp_indices_ordered.append(current_cp_idx)

        if "Laplacian of electron density:" in line:
            val_match = re.search(r':\s*(-?\d+\.\d+(?:[eE][+-]?\d+)?)', line)
            if val_match:
                val = float(val_match.group(1))
                if current_cp_idx is not None:
                    pending_laplacian[current_cp_idx] = val
                    laplacian_found = True

    # Map Rho and Laplacian; filter ghost CCPs
    valid_ccps = {}
    for i, cp_idx in enumerate(ccp_indices_ordered):
        if i < len(rho_ccp_list):
            cps['rho_ccp'][cp_idx] = rho_ccp_list[i]
            valid_ccps[cp_idx] = cps['ccps'][cp_idx]
            if cp_idx in pending_laplacian:
                cps['laplacian_ccp'][cp_idx] = pending_laplacian[cp_idx]
            else:
                print(f"  [CRITICAL WARNING] No Laplacian found for CCP {cp_idx}. "
                      f"Sm will be 1.0 — check CPprop.txt.")
                cps['laplacian_ccp'][cp_idx] = 0.0
        else:
            print(f"  [WARNING] Ghost CCP {cp_idx} removed during parsing.")

    cps['ccps'] = valid_ccps

    if not laplacian_found and len(cps['ccps']) > 0:
        print("\n[CRITICAL WARNING] No Laplacian values found in CPprop.txt.\n")

    return cps


# ---------------------------------------------------------------------------
# PARAMETER-FREE CAGE ATOM IDENTIFICATION
# ---------------------------------------------------------------------------

def _is_inside_hull(point, hull):
    """
    True if *point* is on or inside *hull*.
    Uses the half-space representation: all(A·x + b ≤ 0) inside.
    A small tolerance (1e-10) absorbs floating-point rounding on the surface.
    """
    return np.all(hull.equations @ np.append(point, 1) <= 1e-10)


def get_cage_atoms_delaunay(ccp_idx, ccp_coords, nuclei_dict):
    """
    Parameter-free cage atom identification via Delaunay triangulation.

    The CCP is inserted as point 0 into the Delaunay triangulation of
    all heavy nuclei.  Every nucleus that shares at least one Delaunay
    simplex with the CCP is counted as a cage atom.

    No radius, no scaling factor, no empirical threshold — the cage
    geometry is determined entirely by the QTAIM critical point
    topology and the nuclear coordinate set.

    Enclosure guard: the CCP must lie inside the ConvexHull of all
    heavy nuclei.  CCPs outside that hull are surface CPs, not cage
    CPs, and are rejected before the triangulation step.
    """
    heavy = {k: v for k, v in nuclei_dict.items() if v['symbol'] != 'H'}
    if len(heavy) < 4:
        return []

    atom_indices = list(heavy.keys())
    atom_coords  = np.array([heavy[i]['coords'] for i in atom_indices])

    # --- Enclosure check ------------------------------------------------
    try:
        mol_hull = ConvexHull(atom_coords)
    except Exception:
        return []

    if not _is_inside_hull(ccp_coords, mol_hull):
        print(f"  [INFO] CCP {ccp_idx}: outside molecular hull — "
              f"surface CP, not counted as cage.")
        return []

    # --- Delaunay: CCP inserted as vertex 0 -----------------------------
    all_points = np.vstack([ccp_coords.reshape(1, 3), atom_coords])
    try:
        tri = Delaunay(all_points)
    except Exception as e:
        print(f"  [WARNING] CCP {ccp_idx}: Delaunay triangulation failed — {e}")
        return []

    cage_local = set()
    for simplex in tri.simplices:
        if 0 in simplex:
            cage_local.update(simplex)
    cage_local.discard(0)   # remove the CCP itself

    if len(cage_local) < 3:
        print(f"  [WARNING] CCP {ccp_idx}: only {len(cage_local)} Delaunay "
              f"neighbors — cage may be ill-defined.")
        return []

    # Map Delaunay-local indices back to atom numbers
    # (offset by 1 because the CCP occupies index 0)
    return [atom_indices[i - 1] for i in cage_local]


def extract_topological_cages(cps):
    """
    For each valid CCP, identify cage atoms via parameter-free Delaunay
    triangulation and return (ccp_idx, cage_atoms) pairs.
    """
    matched_cages = []
    for ccp_idx, ccp_coords in cps['ccps'].items():
        cage_atoms = get_cage_atoms_delaunay(ccp_idx, ccp_coords, cps['nuclei'])
        if len(cage_atoms) >= 3:
            matched_cages.append((ccp_idx, cage_atoms))
    return matched_cages


def compute_cage_assembly_tax(cage_atoms, rho, laplacian, ccp_idx):
    if not (0 < rho < 1.0):
        raise ValueError(
            f"CRITICAL ERROR: Unphysical CCP density (rho={rho}) for CCP {ccp_idx}.")

    N_cage = len(cage_atoms)
    reduced_laplacian = abs(laplacian) / (rho * N_cage)
    Sm  = 1.0 + reduced_laplacian
    tax = Sm * (-np.log(rho))
    return reduced_laplacian, Sm, tax


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    args = sys.argv[1:]
    file_args = [a for a in args if not a.startswith('--')]

    if not file_args:
        file_args = [f for f in ["summary.txt", "CPprop.txt"] if os.path.exists(f)]
        if not file_args:
            print("Error: No valid QTAIM files found.")
            sys.exit(1)

    try:
        cps = parse_qtaim_files(file_args)

        print("=" * 60)
        print("LAPLACIAN CAGE TAX (Sm) ANALYSIS")
        print("[Parameter-free Delaunay cage identification]")
        print("=" * 60)

        total_cage_tax = 0.0

        if len(cps['ccps']) == 0:
            print("\n--- Zero Cages (CCPs) detected ---")
        else:
            cages = extract_topological_cages(cps)

            print(f"\nTotal valid CCPs found: {len(cages)}")

            for ccp_idx, cage_atoms in cages:
                rho       = cps['rho_ccp'][ccp_idx]
                laplacian = cps['laplacian_ccp'][ccp_idx]

                red_lap, Sm, tax = compute_cage_assembly_tax(
                    cage_atoms, rho, laplacian, ccp_idx)
                total_cage_tax += tax

                print(f"\n--- CCP Index: {ccp_idx} ---")
                print(f"Cage Heavy Atoms (N_cage)      : {len(cage_atoms)}")
                print(f"Cage Atom Indices              : {cage_atoms}")
                print(f"CCP Electron Density (rho)     : {rho:.5f} a.u.")
                print(f"CCP Laplacian (del^2 rho)      : {laplacian:.5f} a.u.")
                print(f"Reduced Laplacian / N_cage     : {red_lap:.4f}")
                print(f"Parameter-Free Sm Multiplier   : {Sm:.5f}")
                print(f"--> Evolution Cage Tax         : {tax:.2f}")

        print("\n" + "-" * 60)
        print(f"TOTAL TOPOLOGICAL CAGE TAX: {total_cage_tax:.2f}")
        print("-" * 60 + "\n")

    except Exception as e:
        import traceback
        print(f"An error occurred: {e}")
        traceback.print_exc()
