import numpy as np
import networkx as nx
import re
import sys
import os

def parse_qtaim_files(file_paths):
    cps = {'nuclei': {}, 'bcps': [], 'bcp_data': {}, 'rcps': {}, 'rho_rcp': {}}
    rcp_indices_ordered = []
    rho_rcp_list = []
    
    bcp_idx_to_edge = {}
    pending_ellipticity = {}
    ellipticity_found_anywhere = False

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

    rho_match = re.search(r'\(3,\+1\) = \{(.*?)\}', full_content, re.DOTALL)
    if rho_match:
        inner_text = rho_match.group(1).strip()
        if inner_text:
            rho_rcp_list = [float(x) for x in inner_text.split(',') if x.strip()]

    current_cp_idx = None

    for line in all_lines:
        if not line.strip(): continue

        if "CP" in line and "Type (3," in line:
            match = re.search(r'CP\s+(\d+)', line, re.IGNORECASE)
            if match:
                current_cp_idx = int(match.group(1))
        elif re.match(r'^\s*\d+', line) and "(3," in line:
            try:
                current_cp_idx = int(line.split()[0])
            except Exception: pass

        float_matches = re.findall(r'-?\d+\.\d+(?:[eE][+-]?\d+)?', line)
        coords = np.array([float(float_matches[0]), float(float_matches[1]), float(float_matches[2])]) if len(float_matches) >= 3 else np.zeros(3)

        if "(3,-3)" in line:
            atom_match = re.search(r'(\d+)\s*\(\s*([A-Za-z]+)\s*\)', line)
            if atom_match:
                atom_num = int(atom_match.group(1))
                atom_symbol = atom_match.group(2)
                cps['nuclei'][atom_num] = {'coords': coords, 'symbol': atom_symbol}
            elif "Nucleus:" in line:
                nuc_match = re.search(r'Nucleus:\s*([A-Za-z]+)\s*(\d+)', line)
                if nuc_match:
                    atom_symbol = nuc_match.group(1)
                    atom_num = int(nuc_match.group(2))
                    if atom_num not in cps['nuclei']:
                        cps['nuclei'][atom_num] = {'coords': coords, 'symbol': atom_symbol}

        elif "(3,-1)" in line and "--" in line:
            atoms = re.findall(r'(\d+)\s*\(', line.split("(3,-1)")[1])
            if len(atoms) >= 2:
                u, v = int(atoms[0]), int(atoms[1])
                cps['bcps'].append((u, v))
                edge_key = tuple(sorted((u, v)))

                if edge_key not in cps['bcp_data']:
                    cps['bcp_data'][edge_key] = {'coords': coords, 'ellipticity': 0.0}
                elif len(float_matches) >= 3:
                    cps['bcp_data'][edge_key]['coords'] = coords

                if current_cp_idx is not None:
                    bcp_idx_to_edge[current_cp_idx] = edge_key
                    if current_cp_idx in pending_ellipticity:
                        cps['bcp_data'][edge_key]['ellipticity'] = pending_ellipticity[current_cp_idx]
                        ellipticity_found_anywhere = True

        elif "(3,+1)" in line and len(float_matches) >= 3:
            if current_cp_idx is not None:
                cps['rcps'][current_cp_idx] = coords
                if current_cp_idx not in rcp_indices_ordered:
                    rcp_indices_ordered.append(current_cp_idx)

        if "Ellipticity of electron density" in line:
            val_match = re.search(r':\s*(-?\d+\.\d+(?:[eE][+-]?\d+)?)', line)
            if val_match:
                val = float(val_match.group(1))
                if current_cp_idx is not None:
                    if current_cp_idx in bcp_idx_to_edge:
                        edge_key = bcp_idx_to_edge[current_cp_idx]
                        cps['bcp_data'][edge_key]['ellipticity'] = val
                        ellipticity_found_anywhere = True
                    else:
                        pending_ellipticity[current_cp_idx] = val

    # Map rho to RCPs and filter ghosts
    valid_rcps = {}
    for i, cp_idx in enumerate(rcp_indices_ordered):
        if i < len(rho_rcp_list):
            cps['rho_rcp'][cp_idx] = rho_rcp_list[i]
            valid_rcps[cp_idx] = cps['rcps'][cp_idx]
        else:
            print(f"  [WARNING] Ghost RCP {cp_idx} removed during parsing.")
    cps['rcps'] = valid_rcps

    if len(cps['bcps']) > 0 and not ellipticity_found_anywhere:
        print("\n[CRITICAL WARNING] No Ellipticity values found in the provided files!")
        print("Ensure you provide the CPprop.txt file.")
        print("All bond ellipticities will default to 0.0 (Sk = 1.0).\n")

    return cps


def find_all_small_cycles(G, max_size=16):
   
    all_cycles = set()

    def dfs(start, current, path, path_set):
        for neighbor in G.neighbors(current):
            if neighbor == start and len(path) >= 3:
                # Found a valid cycle back to start
                all_cycles.add(tuple(sorted(path)))
                continue
            if neighbor in path_set:
                continue
            if len(path) >= max_size:
                continue
            path_set.add(neighbor)
            dfs(start, neighbor, path + [neighbor], path_set)
            path_set.discard(neighbor)

    for node in G.nodes():
        dfs(node, node, [node], {node})

    return [list(c) for c in all_cycles]


def is_valid_cycle(cycle, bcp_set):
    """
    Check that the cycle atoms form a genuine simple ring in the BCP graph.

    Condition: every atom in the cycle must have exactly 2 neighbours
    within the cycle that are connected by real BCPs. This is the
    necessary and sufficient condition for a set of atoms to form a
    simple ring, independent of the order atoms were enumerated.

    This correctly rejects composite path artifacts in fused polycyclic
    systems where the DFS finds closed walks that are not simple rings.
    It also correctly accepts valid rings regardless of atom ordering.
    """
    cycle_set = set(cycle)
    for atom in cycle:
        bcp_neighbors_in_cycle = [
            n for n in cycle_set
            if n != atom and tuple(sorted((atom, n))) in bcp_set
        ]
        if len(bcp_neighbors_in_cycle) != 2:
            return False
    return True


def extract_topological_rings(cps, max_ring_size=12):
    """
    Per-RCP cycle matching using all small cycles up to max_ring_size.

    For each RCP independently, find the VALID cycle (all bonds exist as
    real BCPs) whose geometric centroid is closest to the RCP coordinates.

    Key fix v3.9: cycles where consecutive atoms have no BCP are rejected.
    In fused polycyclic systems the DFS can return composite paths that
    form closed walks in the graph but are not real rings. Bond validation
    filters these out before centroid matching.

    Duplicate mapping detection: if two RCPs map to the same cycle, a
    warning is printed. Increase --maxring if duplicates occur for large
    macrocyclic molecules.
    """
    matched_rings = []
    used_cycles = {}  # cycle_key -> rcp_idx that claimed it first
    G = nx.Graph()
    G.add_edges_from(cps['bcps'])

    # Pre-compute BCP set for O(1) bond lookup during validation
    bcp_set = set(tuple(sorted(e)) for e in cps['bcps'])

    print(f"  [INFO] Generating all cycles up to size {max_ring_size}...")
    all_cycles = find_all_small_cycles(G, max_size=max_ring_size)

    # Filter to only valid cycles where every bond exists as a real BCP
    valid_cycles = [c for c in all_cycles if is_valid_cycle(c, bcp_set)]
    print(f"  [INFO] Found {len(all_cycles)} raw cycles, {len(valid_cycles)} valid (all bonds confirmed) for {len(cps['rcps'])} RCPs.")

    if not valid_cycles:
        print("  [CRITICAL ERROR] No valid cycles found in bond graph.")
        return matched_rings

    for rcp_idx, rcp_coords in cps['rcps'].items():
        best_cycle = None
        min_dist = float('inf')

        for cycle in valid_cycles:
            try:
                coords = np.array([cps['nuclei'][n]['coords'] for n in cycle])
            except KeyError:
                continue

            centroid = np.mean(coords, axis=0)
            dist = np.linalg.norm(centroid - rcp_coords)

            if dist < min_dist:
                min_dist = dist
                best_cycle = cycle

        if best_cycle:
            matched_rings.append((rcp_idx, best_cycle, min_dist))
            cycle_key = tuple(sorted(best_cycle))

            # Duplicate mapping check
            if cycle_key in used_cycles:
                print(f"  [WARNING] RCP {rcp_idx} and RCP {used_cycles[cycle_key]} both map to cycle {list(cycle_key)}.")
                print(f"           Topology ambiguity — Sk for RCP {rcp_idx} is approximate.")
            else:
                used_cycles[cycle_key] = rcp_idx

            # Plausibility check — large distance may mean max_ring_size too small
            if min_dist > 2.0:
                print(f"  [WARNING] RCP {rcp_idx} centroid distance is {min_dist:.4f} Bohr.")
                print(f"           Best matched cycle size: {len(best_cycle)}. If this seems wrong,")
                print(f"           re-run with --maxring {max_ring_size + 4} to search larger cycles.")
        else:
            print(f"  [CRITICAL ERROR] RCP {rcp_idx} could not be matched to ANY cycle!")

    return matched_rings


def sort_ring_atoms(ring_atoms, bcps_list):
    if not ring_atoms: return ring_atoms
    ring_set = set(ring_atoms)
    ring_adj = {a: [] for a in ring_atoms}
    for (u, v) in bcps_list:
        if u in ring_set and v in ring_set:
            ring_adj[u].append(v)
            ring_adj[v].append(u)
    ordered = [ring_atoms[0]]
    visited = {ring_atoms[0]}
    prev = None
    current = ring_atoms[0]
    for _ in range(len(ring_atoms) - 1):
        # Exclude prev AND already-visited atoms to prevent looping
        neighbors = [n for n in ring_adj[current] if n != prev and n not in visited]
        if not neighbors:
            print(f"  [WARNING] Ring walk produced invalid path for atoms {ring_atoms}. Using unsorted.")
            # Validate unsorted version too
            if len(set(ring_atoms)) != len(ring_atoms):
                print(f"  [CRITICAL WARNING] Duplicate atoms in ring, cycle artifact. Sk will be wrong.")
            return ring_atoms
        prev = current
        current = neighbors[0]
        visited.add(current)
        ordered.append(current)

    # Final sanity check
    if len(ordered) != len(ring_atoms) or len(set(ordered)) != len(ring_atoms):
        print(f"  [WARNING] Ring walk produced invalid path for atoms {ring_atoms}. Using unsorted.")
        return ring_atoms

    return ordered


def compute_ring_assembly_tax(ring_atoms, cps_data, rho, rcp_idx, verbose=False):
    if not (0 < rho < 1.0):
        raise ValueError(f"CRITICAL ERROR: Unphysical RCP density (rho={rho}) for RCP {rcp_idx}.")

    # Reject rings with duplicate atoms — indicates cycle artifact not caught upstream
    if len(set(ring_atoms)) != len(ring_atoms):
        print(f"  [CRITICAL WARNING] RCP {rcp_idx}: duplicate atoms in ring {ring_atoms}. Skipping — increase --maxring.")
        return 0.0, 1.0, 0.0

    num_bonds = len(ring_atoms)
    valid_ellipticities = []

    for i in range(num_bonds):
        u = ring_atoms[i]
        v = ring_atoms[(i + 1) % num_bonds]
        edge_key = tuple(sorted((u, v)))

        if edge_key in cps_data['bcp_data']:
            # NO ORGANIC_ELEMENTS FILTER — pure density, no atom identity needed
            ellip = cps_data['bcp_data'][edge_key]['ellipticity']
            valid_ellipticities.append(ellip)
        else:
            print(f"  --> [WARNING] No BCP found for bond {u}-{v}. Sk underestimated!")

    mean_ellipticity = np.mean(valid_ellipticities) if valid_ellipticities else 0.0
    Sk = 1.0 + mean_ellipticity
    ring_assembly_resistance = Sk * (-np.log(rho))

    return mean_ellipticity, Sk, ring_assembly_resistance


if __name__ == "__main__":
    args = sys.argv[1:]
    verbose = '--verbose' in args

    # Parse --maxring N argument (default 12, increase for macrocycles)
    max_ring_size = 12
    for arg in args:
        if arg.startswith('--maxring'):
            try:
                max_ring_size = int(arg.split('=')[1]) if '=' in arg else int(args[args.index(arg) + 1])
            except (IndexError, ValueError):
                print("Warning: Could not parse --maxring value. Using default 12.")

    file_args = [arg for arg in args if not arg.startswith('--')]

    if not file_args:
        file_args = [f for f in ["summary.txt", "CPprop.txt"] if os.path.exists(f)]
        if not file_args:
            print("Error: No valid QTAIM files found or specified.")
            sys.exit(1)

    try:
        cps = parse_qtaim_files(file_args)

        if verbose:
            print("\n--- DIAGNOSTIC CHECK ---")
            print(f"Successfully loaded from {file_args}:")
            print(f"-> {len(cps['nuclei'])} Nuclei, {len(cps['bcps'])} BCPs, {len(cps['rcps'])} RCPs.")
            print("\n--- BCP ELLIPTICITY MAPPING ---")
            for edge, data in cps['bcp_data'].items():
                print(f"BCP {edge} | Ellipticity: {data['ellipticity']:.5f}")
            print("------------------------\n")

        print("=" * 60)
        print("ELLIPTICITY RING TAX (Sk) ANALYSIS")
        print("=" * 60)

        total_molecule_ring_tax = 0.0

        if len(cps['rcps']) == 0:
            print("\n--- Zero Rings (RCPs) detected in this molecule ---")
            print("--> Baseline Acyclic Sk          : 1.00000")
            print("--> Evolution Ring Tax           : 0.00")
        else:
            rings = extract_topological_rings(cps, max_ring_size=max_ring_size)

            for rcp_idx, ring_atoms, min_dist in rings:

                if rcp_idx not in cps['rho_rcp']:
                    print(f"  [WARNING] RCP {rcp_idx} skipped (No Rho value).")
                    continue

                rho = cps['rho_rcp'][rcp_idx]
                ring_atoms_ordered = sort_ring_atoms(ring_atoms, cps['bcps'])

                mean_eps, Sk, tax = compute_ring_assembly_tax(
                    ring_atoms_ordered, cps, rho, rcp_idx, verbose)
                total_molecule_ring_tax += tax

                print(f"\n--- RCP Index: {rcp_idx} ---")
                print(f"Topological Ring Atoms         : {ring_atoms_ordered}")
                print(f"Ring Size (N-membered)         : {len(ring_atoms_ordered)}")
                print(f"RCP Electron Density (rho)     : {rho:.5f} a.u.")
                print(f"Mean Bond Ellipticity <eps>    : {mean_eps:.4f}")
                print(f"Harmonic Sk Multiplier (V3.9)  : {Sk:.5f}")
                print(f"--> Evolution Ring Tax         : {tax:.2f}")

        print("\n" + "-" * 60)
        print(f"TOTAL TOPOLOGICAL RING TAX: {total_molecule_ring_tax:.2f}")
        print("-" * 60 + "\n")
        print(f"-> {len(cps['nuclei'])} Nuclei, {len(cps['bcps'])} BCPs, {len(cps['rcps'])} RCPs.")
    except Exception as e:
        import traceback
        print(f"An execution error occurred: {e}")
        traceback.print_exc()
