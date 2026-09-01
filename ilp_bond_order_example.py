#!/usr/bin/env python3
"""
Minimal ILP example for bond-order assignment from an XYZ input.

Console layout follows Lewis-engine.py: NBO-style electron summary,
$CHOOSE block, octet / valence check, then CBC ligand classification.

XYZ layout: line 1 = atom count, line 2 = comment (required, may be empty),
coordinates always begin on line 3.

Usage:
    python ilp_bond_order_example.py ABATOQ.xyz
"""

from __future__ import annotations

import math
import os
import sys
from collections import defaultdict

LOCAL_VENDOR = os.path.join(os.path.dirname(__file__), "_vendor")
if os.path.isdir(LOCAL_VENDOR) and LOCAL_VENDOR not in sys.path:
    sys.path.insert(0, LOCAL_VENDOR)

try:
    import pulp
except ImportError as exc:  # pragma: no cover - runtime guidance
    raise SystemExit(
        "PuLP is not installed. Install with: pip install pulp"
    ) from exc


# Covalent radii aligned with Lewis-engine.py (Angstrom).
COV_R = {
    "H": 0.31, "He": 0.28,
    "Li": 0.25, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58,
    "Na": 0.25, "Mg": 0.72, "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Ar": 1.06,
    "K": 0.25, "Ca": 0.25,
    "Sc": 1.44, "Ti": 1.36, "V": 1.25, "Cr": 1.22, "Mn": 1.19, "Fe": 1.16, "Co": 1.11,
    "Ni": 1.10, "Cu": 1.12, "Zn": 1.18,
    "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20, "Kr": 1.16,
    "Rb": 0.25, "Sr": 0.25,
    "Y": 1.62, "Zr": 1.48, "Nb": 1.37, "Mo": 1.45, "Tc": 1.56, "Ru": 1.26, "Rh": 1.35,
    "Pd": 1.24, "Ag": 1.45, "Cd": 1.44,
    "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "Xe": 1.40,
    "Cs": 0.25, "Ba": 0.25,
    "La": 1.94, "Ce": 1.83, "Pr": 1.82, "Nd": 1.81, "Pm": 1.80, "Sm": 1.80, "Eu": 1.99,
    "Gd": 1.79, "Tb": 1.76, "Dy": 1.75, "Ho": 1.74, "Er": 1.73, "Tm": 1.72, "Yb": 1.94, "Lu": 1.72,
    "Hf": 1.52, "Ta": 1.46, "W": 1.37, "Re": 1.31, "Os": 1.44, "Ir": 1.41,
    "Pt": 1.36, "Au": 1.36, "Hg": 1.32,
    "Tl": 1.45, "Pb": 1.46, "Bi": 1.48,
}

# Backward-compatible alias used by some existing helpers in this file.
RADII = COV_R

TM_SET = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu",
}

# Common TM oxidation states from tmQM_coordination_results_G_smiles.csv
# Criteria: > 5% of per-metal bracket parses (SMILES_tmQM + SMILES_CSD_fixed).
# Empty frozenset: no parses or no state above 5%. See tmQM_common_oxidation_states.log.
TM_COMMON_OXIDATION_STATES = {
    "Sc": [3],
    "Ti": [3, 4],
    "V": [2, 3, 4, 5],
    "Cr": [2, 3, 4, 6],
    "Mn": [2, 3, 4, 6, 7],
    "Fe": [2, 3],
    "Co": [2, 3],
    "Ni": [2],
    "Cu": [1, 2],
    "Zn": [2],
    "Y": [3],
    "Zr": [4],
    "Nb": [3, 4, 5],
    "Mo": [2, 3, 4, 5, 6],
    "Tc": [2, 3, 4, 5, 6, 7],
    "Ru": [2, 3, 4, 5, 6, 7, 8],
    "Rh": [1, 3],
    "Pd": [2, 4],
    "Ag": [1],
    "Cd": [2],
    "La": [3],
    "Hf": [4],
    "Ta": [3, 4, 5],
    "W": [2, 3, 4, 5, 6],
    "Re": [2, 3, 4, 5, 6, 7],
    "Os": [3, 4, 5, 6, 7, 8],
    "Ir": [1, 3],
    "Pt": [2, 4],
    "Au": [1, 3],
    "Hg": [1, 2],
}

# Heavy alkali (Na, K, …) and alkaline earth (Mg, Ca, …): same ionic limit in the ILP —
# q = ve − bond_sum, no LP variables, no octet slack.  Matches Lewis-engine’s small
# COV_R / cation heuristic for these metals; avoids spurious octet-8 penalties (e.g. K⁺).
# Caveat: Mg and covalent organometallics (Grignard, etc.) are not pure M²⁺; this is a
# deliberate simplification aligned with “ionic s-block when barely bonded”.
ALKALI_HEAVY_IONIC = frozenset({"Na", "K", "Rb", "Cs", "Fr"})
# Be is intentionally omitted: strongly covalent; use full lp/octet model like other non-TM.
ALKALINE_EARTH_IONIC = frozenset({"Mg", "Ca", "Sr", "Ba", "Ra"})


def is_ionlike_s_block_metal(sym: str) -> bool:
    return sym in ALKALI_HEAVY_IONIC or sym in ALKALINE_EARTH_IONIC


# Standard covalent bond capacities — aligned with Lewis-engine.py STD_CAP.
STD_CAP = {
    "H": 1, "He": 0, "Li": 1, "Be": 2, "B": 3, "C": 4, "N": 3, "O": 2, "F": 1, "Ne": 0,
    "Na": 1, "Mg": 2, "Al": 3, "Si": 4, "P": 3, "S": 2, "Cl": 1, "Ar": 0,
    "K": 1, "Ca": 2, "Ga": 3, "Ge": 4, "As": 3, "Se": 2, "Br": 1, "Kr": 0,
    "Rb": 1, "Sr": 2, "In": 3, "Sn": 4, "Sb": 3, "Te": 2, "I": 1, "Xe": 0,
    "Cs": 1, "Ba": 2, "Tl": 3, "Pb": 4, "Bi": 3,
}

# Neutral-atom valence electron totals — aligned with Lewis-engine.py VALENCE.
VALENCE = {
    "H": 1, "He": 2,
    "Li": 1, "Be": 2, "B": 3, "C": 4, "N": 5, "O": 6, "F": 7, "Ne": 8,
    "Na": 1, "Mg": 2, "Al": 3, "Si": 4, "P": 5, "S": 6, "Cl": 7, "Ar": 8,
    "K": 1, "Ca": 2, "Ga": 3, "Ge": 4, "As": 5, "Se": 6, "Br": 7, "Kr": 8,
    "Rb": 1, "Sr": 2, "In": 3, "Sn": 4, "Sb": 5, "Te": 6, "I": 7, "Xe": 8,
    "Cs": 1, "Ba": 2, "Tl": 3, "Pb": 4, "Bi": 5,
    # TM: d+s electrons of neutral ground-state atom (same as Lewis-engine)
    "Sc": 3, "Ti": 4, "V": 5, "Cr": 6, "Mn": 7, "Fe": 8, "Co": 9, "Ni": 10, "Cu": 11, "Zn": 12,
    "Y": 3, "Zr": 4, "Nb": 5, "Mo": 6, "Tc": 7, "Ru": 8, "Rh": 9, "Pd": 10, "Ag": 11, "Cd": 12,
    "Hf": 4, "Ta": 5, "W": 6, "Re": 7, "Os": 8, "Ir": 9, "Pt": 10, "Au": 11, "Hg": 12,
    "La": 3, "Ce": 4, "Pr": 5, "Nd": 6, "Pm": 7, "Sm": 8, "Eu": 9, "Gd": 10, "Tb": 11,
    "Dy": 12, "Ho": 13, "Er": 14, "Tm": 15, "Yb": 16, "Lu": 17,
}

# Same keys as STD_CAP; kept for scripts that still expect this name.
VALENCE_TARGET = STD_CAP

# Viewer bridge expects these globals to exist.
CORE_E = {k: 0 for k in COV_R}
# Pauling electronegativity — aligned with Lewis-engine.py
ENEG = {
    "F": 3.98,
    "O": 3.44,
    "Cl": 3.16,
    "N": 3.04,
    "Br": 2.96,
    "I": 2.66,
    "S": 2.58,
    "Se": 2.55,
    "C": 2.55,
    "P": 2.19,
    "H": 2.20,
    "As": 2.18,
    "Te": 2.10,
    "Si": 1.90,
    "B": 2.04,
    "Ge": 2.01,
    "Sn": 1.96,
    "Sb": 2.05,
    "Pb": 2.33,
    "Al": 1.61,
    "Ga": 1.81,
    "In": 1.78,
    "Tl": 2.04,
}

# Neutral valence electron counts for formal charge (same as Lewis-engine VALENCE).
VALENCE_ELECTRONS = dict(VALENCE)


def _looks_like_xyz_atom_line(line: str) -> bool:
    """True if stripped line has element + three floats (common coord row)."""
    p = line.strip().split()
    if len(p) < 4:
        return False
    try:
        float(p[1])
        float(p[2])
        float(p[3])
    except ValueError:
        return False
    head = p[0]
    if not head:
        return False
    # reject pure numbers as element
    try:
        float(head)
        return False
    except ValueError:
        pass
    return True


def read_xyz(path: str):
    """
    XYZ reader — always uses the same layout (1-based line numbers):

      line 1: integer n (number of atoms)
      line 2: comment / title (ignored for coordinates; may be empty)
      lines 3 .. n+2: element + x y z

    Line 2 must exist. Coordinates always start at line 3, even for minimal files:
    put a placeholder comment on line 2 if you have nothing to say.
    """
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    if not lines:
        raise ValueError(f"empty XYZ file: {path}")

    try:
        n = int(lines[0].strip().split()[0])
    except (ValueError, IndexError) as exc:
        raise ValueError(
            f"{path}: line 1 must be a positive integer atom count, got {lines[0]!r}"
        ) from exc
    if n <= 0:
        raise ValueError(f"{path}: atom count must be positive, got {n}")
    if len(lines) < 2:
        raise ValueError(
            f"{path}: need line 2 (comment/title). Coordinates start at line 3."
        )
    if _looks_like_xyz_atom_line(lines[1]):
        raise ValueError(
            f"{path}: line 2 looks like a coordinate row, but this program always "
            f"treats line 2 as the comment line. Insert a title/comment on line 2 "
            f"and move coordinates to lines 3–{n + 2}."
        )
    need = 2 + n
    if len(lines) < need:
        raise ValueError(
            f"{path}: for n={n}, need at least {need} lines "
            f"(1 count + 1 comment + {n} coords starting at line 3); "
            f"got {len(lines)}"
        )

    atoms = []
    for k in range(n):
        file_line_no = 3 + k
        ln = lines[2 + k]
        p = ln.split()
        if len(p) < 4:
            raise ValueError(
                f"{path}: line {file_line_no}: expected element + x y z, got {ln!r}"
            )
        raw = p[0]
        sym = raw[0].upper() + raw[1:].lower() if len(raw) > 1 else raw.upper()
        try:
            x, y, z = float(p[1]), float(p[2]), float(p[3])
        except ValueError as exc:
            raise ValueError(
                f"{path}: line {file_line_no}: non-numeric coordinates in {ln!r}"
            ) from exc
        atoms.append((k + 1, sym, x, y, z))
    return atoms


def dist(a, b) -> float:
    def _xyz(p):
        if len(p) >= 5:
            return p[2], p[3], p[4]
        return p[0], p[1], p[2]

    ax, ay, az = _xyz(a)
    bx, by, bz = _xyz(b)
    return math.sqrt((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2)


def is_tm(symbol: str) -> bool:
    return symbol in TM_SET


def is_TM(symbol: str) -> bool:
    """Compatibility alias for viewer bridge code."""
    return is_tm(symbol)


def _prune_h_to_closest_nonmetal_neighbor(atoms, edges):
    """
    If H is within cutoff of more than one non-metal, keep only the closest contact.
    TM–H edges are unchanged (not counted toward the one-nonmetal limit).
    """
    by_h: dict[int, list[tuple[tuple[int, int, str, str], float]]] = defaultdict(list)
    for i, j, ei, ej in edges:
        if ei != "H" and ej != "H":
            continue
        if ei == "H":
            h_idx, other_idx, other_sym = i, j, ej
        else:
            h_idx, other_idx, other_sym = j, i, ei
        if is_tm(other_sym):
            continue
        d = dist(atoms[h_idx - 1], atoms[other_idx - 1])
        by_h[h_idx].append(((i, j, ei, ej), d))

    drop = set()
    for h_idx, candidates in by_h.items():
        if len(candidates) <= 1:
            continue
        candidates.sort(key=lambda item: item[1])
        for edge, _d in candidates[1:]:
            drop.add(edge)

    if not drop:
        return edges
    return [e for e in edges if e not in drop]


def connectivity(atoms, factor=1.30):
    """Match Lewis-engine style raw connectivity generation."""
    edges = []
    n = len(atoms)
    for i in range(n):
        ai = atoms[i]
        for j in range(i + 1, n):
            aj = atoms[j]
            ei, ej = ai[1], aj[1]
            if ei not in COV_R or ej not in COV_R:
                continue
            cutoff = COV_R[ei] + COV_R[ej]
            # Match Lewis-engine rule: tighter cutoff for TM-C contacts.
            if (is_tm(ei) and ej == "C") or (is_tm(ej) and ei == "C"):
                local_factor = 1.15
            else:
                local_factor = factor
            cutoff *= local_factor
            if dist(ai, aj) < cutoff:
                edges.append((ai[0], aj[0], ei, ej))
    return _prune_h_to_closest_nonmetal_neighbor(atoms, edges)


def _sigma_tm_c_set(atoms, edges):
    """Mirror Lewis-engine: keep only sigma TM-C; drop haptic/agostic contacts."""
    atom_symbol = {idx: el for idx, el, *_ in atoms}
    tm_c = defaultdict(set)
    adj = defaultdict(set)
    tm_h = defaultdict(set)

    for i, j, _ei, _ej in edges:
        adj[i].add(j)
        adj[j].add(i)
        si, sj = atom_symbol[i], atom_symbol[j]
        if is_tm(si) and sj == "C":
            tm_c[i].add(j)
        if is_tm(sj) and si == "C":
            tm_c[j].add(i)
        if is_tm(si) and sj == "H":
            tm_h[i].add(j)
        if is_tm(sj) and si == "H":
            tm_h[j].add(i)

    sigma = set()
    for tm, c_set in tm_c.items():
        for c in c_set:
            # haptic: C-neighbour also bound to same TM
            haptic = any(
                (nb != tm and (not is_tm(atom_symbol[nb])) and nb in c_set)
                for nb in adj[c]
            )
            if haptic:
                continue

            # agostic: C-H where that H also binds same TM
            h_neighbors = [nb for nb in adj[c] if atom_symbol[nb] == "H" and nb != tm]
            agostic = any(h in tm_h[tm] for h in h_neighbors)
            if agostic:
                continue

            sigma.add((tm, c))
    return sigma


def remove_dative_tm_bonds(atoms, edges):
    """Mirror Lewis-engine TM filter for candidate-edge stage."""
    atom_symbol = {idx: el for idx, el, *_ in atoms}
    sigma = _sigma_tm_c_set(atoms, edges)
    kept = []
    for i, j, ei, ej in edges:
        if is_tm(ei) or is_tm(ej):
            tm, c = (i, j) if is_tm(ei) else (j, i)
            if atom_symbol[c] == "C" and (tm, c) in sigma:
                kept.append((i, j, ei, ej))
        else:
            kept.append((i, j, ei, ej))
    return kept


def fixed_sigma_tm_nonmetal_edges(raw_edges, ilp_edges):
    """
    Compatibility shim.
    TM-X (X = N/O/P/S/halogens) contacts are treated as dative-only and must not
    contribute as fixed Lewis sigma bonds in ILP constraints.
    """
    _ = raw_edges, ilp_edges
    return []


def solve_bond_orders(atoms, edges, mol_charge=0, fixed_sigma_edges=None):
    prob = pulp.LpProblem("BondOrderAssignment", pulp.LpMinimize)
    _ = fixed_sigma_edges  # deprecated: TM-X fixed sigma contributions are disabled

    # Binary edge-choice variables for bond order 0/1/2.
    x = {}
    for i, j, _ei, _ej in edges:
        for k in (0, 1, 2):
            x[(i, j, k)] = pulp.LpVariable(f"x_{i}_{j}_{k}", cat="Binary")
        prob += x[(i, j, 0)] + x[(i, j, 1)] + x[(i, j, 2)] == 1

    inc_edges = defaultdict(list)
    for i, j, *_ in edges:
        inc_edges[i].append((i, j))
        inc_edges[j].append((i, j))

    # Octet / formal-charge variables for non-TM atoms.
    # Ionlike s-block (Na…, Mg, Ca, …): q only (q = ve − bond_sum); TM: skipped;
    # others: lp + octet slack + q.
    lp = {}
    oct_plus, oct_minus = {}, {}
    q, abs_q = {}, {}
    q_neg = {}
    for i, el, *_ in atoms:
        if is_tm(el) or el not in VALENCE_ELECTRONS:
            continue
        q[i] = pulp.LpVariable(f"q_{i}", lowBound=-4, upBound=4, cat="Integer")
        abs_q[i] = pulp.LpVariable(f"absq_{i}", lowBound=0, cat="Integer")
        q_neg[i] = pulp.LpVariable(f"qneg_{i}", lowBound=0, cat="Integer")
        if is_ionlike_s_block_metal(el):
            continue
        lp[i] = pulp.LpVariable(f"lp_{i}", lowBound=0, cat="Integer")
        oct_plus[i] = pulp.LpVariable(f"octp_{i}", lowBound=0, cat="Integer")
        oct_minus[i] = pulp.LpVariable(f"octm_{i}", lowBound=0, cat="Integer")

    for i, el, *_ in atoms:
        if i not in q:
            continue
        bond_sum_terms = []
        for a, b in inc_edges.get(i, []):
            bond_sum_terms.append(1 * x[(a, b, 1)] + 2 * x[(a, b, 2)])
        bond_sum = pulp.lpSum(bond_sum_terms) if bond_sum_terms else 0
        ve = VALENCE_ELECTRONS[el]

        if is_ionlike_s_block_metal(el):
            prob += q[i] == ve - bond_sum
            prob += abs_q[i] >= q[i]
            prob += abs_q[i] >= -q[i]
            prob += q_neg[i] >= -q[i]
            continue

        local_electrons = 2 * lp[i] + 2 * bond_sum
        assigned_electrons = 2 * lp[i] + bond_sum
        # Lewis-engine octet(): 2 for H, He, Li; 8 for other non-TM (TM excluded above).
        oct_target = 2 if el in ("H", "Li") else 8

        prob += local_electrons - oct_target == oct_plus[i] - oct_minus[i]
        prob += q[i] == ve - assigned_electrons
        prob += abs_q[i] >= q[i]
        prob += abs_q[i] >= -q[i]
        prob += q_neg[i] >= -q[i]

    # Avoid high-order bonds to hydrogen and halogen.
    for i, j, ei, ej in edges:
        if "H" in (ei, ej) or "Cl" in (ei, ej):
            prob += x[(i, j, 2)] == 0

    # Molecular-charge consistency and FC-magnitude matching term.
    charge_delta_abs = pulp.LpVariable("charge_delta_abs", lowBound=0, cat="Integer")
    total_q = pulp.lpSum(q.values()) if q else 0
    total_abs_q = pulp.lpSum(abs_q.values()) if abs_q else 0
    prob += charge_delta_abs >= total_q - mol_charge
    prob += charge_delta_abs >= -(total_q - mol_charge)

    absq_match_abs = pulp.LpVariable("absq_match_abs", lowBound=0, cat="Integer")
    target_abs_charge = abs(mol_charge)
    prob += absq_match_abs >= total_abs_q - target_abs_charge
    prob += absq_match_abs >= target_abs_charge - total_abs_q

    # Objective:
    # 1) minimize octet violations
    # 2) minimize sum(|formal charge|)
    # 3) minimize total-charge mismatch and abs-charge mismatch
    # 4) softly bias negative charge toward more electronegative atoms
    octet_penalty = pulp.lpSum(oct_plus.values()) + pulp.lpSum(oct_minus.values())
    formal_charge_penalty = total_abs_q
    double_bond_penalty = pulp.lpSum(x[(i, j, 2)] for i, j, *_ in edges)
    eneg_neg_charge_penalty = pulp.lpSum(
        max(0.0, 4.0 - ENEG.get(el, 2.0)) * q_neg[i]
        for i, el, *_ in atoms
        if i in q_neg
    )
    prob += (
        1000 * octet_penalty
        + 100 * formal_charge_penalty
        + 200 * charge_delta_abs
        + 2 * absq_match_abs
        + 100 * double_bond_penalty
        + 20 * eneg_neg_charge_penalty
    )

    status = prob.solve(pulp.PULP_CBC_CMD(msg=False))
    if pulp.LpStatus[status] not in ("Optimal", "Integer Feasible"):
        raise RuntimeError(f"ILP failed: {pulp.LpStatus[status]}")

    bonds = []
    for i, j, _ei, _ej in edges:
        order = 0
        for k in (1, 2):
            if pulp.value(x[(i, j, k)]) > 0.5:
                order = k
        if order > 0:
            bonds.append((i, j, order))
    bonds.sort()
    lp_out = {}
    fc_out = {}
    for i, var in lp.items():
        lp_out[i] = int(round(pulp.value(var)))
    for i, var in q.items():
        fc_out[i] = int(round(pulp.value(var)))
    return bonds, lp_out, fc_out


def format_choose_style(bonds):
    chunks = []
    for i, j, bo in bonds:
        tag = "S" if bo == 1 else "D"
        chunks.append(f"{tag} {i} {j}")
    return "BOND " + " ".join(chunks) + " END"


def _dative_lp_pairs_per_atom(atoms, bo, metal_adjacency_edges=None, coords=None):
    """
    Per-atom count of lone pairs consumed on dative M–L bonds (ILP bond order 0).
    Indices match atoms/bo/adj (0-based array positions).
    """
    n = len(atoms)
    counts = defaultdict(int)
    if metal_adjacency_edges:
        for i, j, ei, ej in metal_adjacency_edges:
            if not is_TM(ei) ^ is_TM(ej):
                continue
            tm_i, lig_i = (i, j) if is_TM(ei) else (j, i)
            key = (min(tm_i, lig_i), max(tm_i, lig_i))
            if bo.get(key, 0) == 0:
                counts[lig_i] += 1
        return counts
    if coords is None:
        return counts
    for i in range(n):
        if not is_TM(atoms[i]):
            continue
        for j in range(n):
            if i == j:
                continue
            key = (min(i, j), max(i, j))
            if bo.get(key, 0) != 0:
                continue
            ri = COV_R.get(atoms[i], 1.5)
            rj = COV_R.get(atoms[j], 0.77)
            if dist(coords[i], coords[j]) < 1.55 * (ri + rj):
                counts[j] += 1
    return counts


def check_octet_violations(
    atoms, bo, adj, lp=None, fc=None, *, metal_adjacency_edges=None, coords=None
):
    """
    Valence check: bond_e + 2×(free LP + LP on dative M–L) vs VALENCE[sym] − fc.
    Matches aromatic ILP when metal_adjacency_edges is supplied.
    """
    violations = []
    check_syms = {
        "C", "N", "O", "S", "P", "B", "Si", "Se", "Te", "F", "Cl", "Br", "I",
    }
    dative_lp = _dative_lp_pairs_per_atom(
        atoms, bo, metal_adjacency_edges, coords=coords
    )

    for i, sym in enumerate(atoms):
        if sym not in check_syms:
            continue
        bond_e = sum(bo.get((min(i, k), max(i, k)), 0) for k in adj[i])
        lp_free = lp.get(i, 0) if lp else 0
        lp_dative = dative_lp.get(i, 0)
        lp_cnt = lp_free + lp_dative
        fc_val = fc[i] if fc else 0
        total_e = bond_e + 2 * lp_cnt
        exp_val = VALENCE.get(sym, 4) - fc_val

        if total_e == exp_val:
            continue

        neutral_val = VALENCE.get(sym, 4)
        if fc_val != 0 and total_e == neutral_val:
            continue

        if sym in ("B", "Si", "P", "Ge", "Sn", "As", "Sb", "Se", "Te", "Pb", "Bi"):
            if bond_e >= 3:
                continue

        violations.append(
            (i, f"{sym}{i+1}", bond_e, lp_free, lp_dative, lp_cnt, exp_val, fc_val)
        )
    return violations


def print_summary_and_choose_ilp(atom_syms, bonds, lp_out, mol_charge):
    """
    Lewis-engine.py print_summary_and_choose layout for ILP results.
    All valence is reported as Lewis (non-Lewis rows zero); NBO-style table only.
    """
    te = float(sum(VALENCE.get(s, 0) for s in atom_syms) - mol_charge)
    if te <= 0:
        te = 1.0
    le, vnl, ryl = te, 0.0, 0.0
    stats_charge = float(mol_charge)
    print("          -------------------------------")
    print(f"                 Total Lewis{le:10.5f}  ({100 * le / te:8.4f}%)")
    print(f"           Valence non-Lewis{vnl:10.5f}  ({100 * vnl / te:8.4f}%)")
    print(f"           Rydberg non-Lewis{ryl:10.5f}  ({100 * ryl / te:8.4f}%)")
    print("          -------------------------------")
    print(f"               Total unit  1{te:10.5f}  ({100.0:8.4f}%)")
    print(f"              Charge unit  1{stats_charge:10.5f}")
    print()
    print_choose_block(bonds, lp_out)


def print_octet_report(
    atom_syms, bo0, lp_full, fc0, *, metal_adjacency_edges=None, coords=None
):
    """Lewis-engine.py-style octet / valence block after $CHOOSE."""
    adj_check = defaultdict(list)
    for i, j in bo0:
        adj_check[i].append(j)
        adj_check[j].append(i)
    violations = check_octet_violations(
        atom_syms,
        bo0,
        adj_check,
        lp=lp_full,
        fc=fc0,
        metal_adjacency_edges=metal_adjacency_edges,
        coords=coords,
    )
    if violations:
        print()
        print("=" * 72)
        print("  ⚠  Octet / Valence Violations in Best Lewis Structure")
        print("=" * 72)
        print(
            "  Atoms where bond_e + 2×LP ≠ expected valence electrons."
            "  LP = free + on dative M–L bonds."
        )
        print("  These indicate a suboptimal Kekulé pattern or unusual geometry.")
        print()
        print(
            f"  {'Atom':<10} {'Bond-e':>6}  {'LP':>4}  {'Total-e':>7}  {'Exp-e':>6}  {'FC':>4}  Notes"
        )
        print(f"  {'-'*68}")
        for idx, label, be, lp_free, lp_dat, lp_cnt, exp_val, fc_val in violations:
            total_e = be + 2 * lp_cnt
            adj_syms = ", ".join(f"{atom_syms[k]}{k + 1}" for k in adj_check[idx])
            lp_str = f"{lp_cnt}" if lp_dat == 0 else f"{lp_cnt}({lp_free}+{lp_dat}d)"
            print(
                f"  {label:<10} {be:>6}  {lp_str:>4}  {total_e:>7}  {exp_val:>6}  {fc_val:>+4}  bonded to: {adj_syms}"
            )
        print()
        print(f"  → {len(violations)} atom(s) violate expected valence in this Lewis structure.")
    else:
        print()
        print("  ✓  All atoms satisfy their expected valence (octet rule OK).")


def print_choose_block(bonds, lp_out):
    print(" $CHOOSE")
    lone_items = sorted((idx, v) for idx, v in lp_out.items() if v > 0)
    if lone_items:
        print("   LONE " + " ".join(f"{idx} {v}" for idx, v in lone_items) + " END")
    line = "   BOND"
    for i, j, bo in bonds:
        if bo == 1:
            tok = f"S {i} {j}"
        elif bo == 2:
            tok = f"D {i} {j}"
        elif bo == 3:
            tok = f"T {i} {j}"
        else:
            tok = f"S {i} {j}"
        if len(line) + 1 + len(tok) > 72:
            print(line)
            line = "       " + tok
        else:
            line += " " + tok
    print(line + " END")
    print(" $END")


def _pack_atoms_coords(atoms, coords):
    return [
        (idx + 1, atoms[idx], coords[idx][0], coords[idx][1], coords[idx][2])
        for idx in range(len(atoms))
    ]


def find_best_lewis(atoms, coords, charge=0):
    """
    Viewer-compatible API.
    Returns (bo, lp, fc, stats) with 0-based bo/lp indexing.
    """
    packed = _pack_atoms_coords(atoms, coords)
    raw_edges = connectivity(packed)
    edges = remove_dative_tm_bonds(packed, raw_edges)
    bonds, lp_out_1b, fc_out_1b = solve_bond_orders(packed, edges, mol_charge=charge)

    bo = {(i - 1, j - 1): o for i, j, o in bonds}
    lp = {i - 1: v for i, v in lp_out_1b.items() if v > 0}
    fc = [0] * len(atoms)
    for i in range(len(atoms)):
        fc[i] = fc_out_1b.get(i + 1, 0)

    stats = {
        "total_e": float(sum(VALENCE.get(a, 0) for a in atoms) - charge),
        "lewis_e": 0.0,
        "val_nl": 0.0,
        "ryd_nl": 0.0,
        "charge": float(charge),
    }
    return bo, lp, fc, stats


_ORGANIC = {"C", "H", "N", "O", "F", "Cl"}


def is_inorganic(sym):
    return sym not in _ORGANIC


_NORMAL_BONDS = {
    "H": 1,
    "B": 3, "C": 4, "N": 3, "O": 2, "F": 1,
    "Si": 4, "P": 3, "S": 2, "Cl": 1,
    "Ge": 4, "As": 3, "Se": 2, "Br": 1,
    "Sn": 4, "Sb": 3, "Te": 2, "I": 1,
    "Pb": 4, "Bi": 3,
    "Al": 3, "Ga": 3, "In": 3, "Tl": 3,
}


def normal_bonds(sym):
    return _NORMAL_BONDS.get(sym, 4)


_OXSTATE_ENEG = dict(ENEG)
_OXSTATE_ENEG.update({
    "Sc": 1.36, "Ti": 1.54, "V": 1.63, "Cr": 1.66, "Mn": 1.55, "Fe": 1.83, "Co": 1.88,
    "Ni": 1.91, "Cu": 1.90, "Zn": 1.65, "Y": 1.22, "Zr": 1.33, "Nb": 1.60, "Mo": 2.16,
    "Tc": 1.90, "Ru": 2.20, "Rh": 2.28, "Pd": 2.20, "Ag": 1.93, "Cd": 1.69,
    "Hf": 1.30, "Ta": 1.50, "W": 2.36, "Re": 1.90, "Os": 2.20, "Ir": 2.20,
    "Pt": 2.28, "Au": 2.54, "Hg": 2.00,
})


def classify_cbc_ligands(atoms, coords, bo, lp, fc, charge=0, *, metal_adjacency_edges=None):
    n = len(atoms)

    def _full_connectivity_dative(atoms, coords, bo_lewis, lp_lewis, dative_ml_pairs):
        # ILP lp[] is non-donated LPs only; each dative M–L still has one pair on the bond.
        dative_lp_by_lig = defaultdict(int)
        for _tm, lig in dative_ml_pairs:
            dative_lp_by_lig[lig] += 1

        def _effective_lp(lig_idx):
            free = lp_lewis.get(lig_idx, 0) if lp_lewis else 0
            return free + dative_lp_by_lig.get(lig_idx, 0)

        lewis_adj = defaultdict(set)
        for (a, b) in bo_lewis:
            lewis_adj[a].add(b)
            lewis_adj[b].add(a)

        tm_sigma_partners = defaultdict(set)
        for a, b in bo_lewis:
            sa, sb = atoms[a], atoms[b]
            if is_TM(sa) and not is_TM(sb):
                tm_sigma_partners[a].add(b)
            elif is_TM(sb) and not is_TM(sa):
                tm_sigma_partners[b].add(a)

        def is_terminal_co_o(lig_idx, tm_idx):
            sym = atoms[lig_idx]
            if sym not in ("O", "N"):
                return False
            nbrs = list(lewis_adj[lig_idx])
            if len(nbrs) != 1:
                return False
            c_nbr = nbrs[0]
            if atoms[c_nbr] != "C":
                return False
            return tm_idx in lewis_adj[c_nbr]

        def _donor_lp_on_c(c_idx, tm_idx):
            """LP on C: free lone pairs and/or one pair on a dative M–C bond (ILP lp may be 0)."""
            lp_free = lp_lewis.get(c_idx, 0)
            if lp_free > 0:
                return lp_free
            if (tm_idx, c_idx) in dative_ml_pairs:
                return 1
            return 0

        def is_lp_carbene_c_to_tm(c_idx, tm_idx):
            """C coordinating to TM (σ, dative, or lp); substituents on this C are backbone."""
            if atoms[c_idx] != "C":
                return False
            ri = COV_R.get(atoms[tm_idx], 1.5)
            rj = COV_R.get("C", 0.77)
            if dist(coords[c_idx], coords[tm_idx]) >= 1.55 * (ri + rj):
                return False
            if c_idx in tm_sigma_partners[tm_idx]:
                return True
            key = (min(c_idx, tm_idx), max(c_idx, tm_idx))
            if bo_lewis.get(key, 0) > 0:
                return True
            return _donor_lp_on_c(c_idx, tm_idx) > 0

        def is_backbone_atom(lig_idx, tm_idx):
            for nbr in lewis_adj[lig_idx]:
                if nbr in tm_sigma_partners[tm_idx]:
                    return True
                if is_lp_carbene_c_to_tm(nbr, tm_idx):
                    return True
            return False

        def is_saturated_no_lp(lig_idx):
            sym = atoms[lig_idx]
            if sym == "C" or is_TM(sym) or sym == "H":
                return False
            if _effective_lp(lig_idx) > 0:
                return False
            cap = STD_CAP.get(sym, 4)
            bond_e = sum(
                bo_lewis.get((min(lig_idx, k), max(lig_idx, k)), 0)
                for k in lewis_adj[lig_idx]
            )
            return bond_e >= cap

        bonds = []
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                ri = COV_R.get(atoms[i], 0.77)
                rj = COV_R.get(atoms[j], 0.77)
                si, sj = atoms[i], atoms[j]
                tm_i, tm_j = is_TM(si), is_TM(sj)
                d_ij = dist(coords[i], coords[j])

                if (tm_i and not tm_j) or (tm_j and not tm_i):
                    tm = i if tm_i else j
                    lig = j if tm_i else i
                    sym_lig = atoms[lig]
                    lp_count = _effective_lp(lig)

                    if is_terminal_co_o(lig, tm):
                        continue
                    if is_saturated_no_lp(lig):
                        continue

                    threshold_std = 1.3 * (ri + rj)
                    if sym_lig in ("N", "O", "S", "P", "F", "Cl", "Br", "I", "Se", "Te"):
                        if d_ij < threshold_std:
                            bonds.append((i, j))
                        elif lp_count >= 1 and d_ij < 1.6 * (ri + rj):
                            if is_backbone_atom(lig, tm):
                                pass
                            else:
                                if sym_lig in ("Br", "I"):
                                    if d_ij < 1.35 * (ri + rj):
                                        bonds.append((i, j))
                                else:
                                    if d_ij < 1.55 * (ri + rj):
                                        bonds.append((i, j))
                    elif sym_lig == "C":
                        if d_ij < 1.15 * (ri + rj):
                            bonds.append((i, j))
                    else:
                        if d_ij < threshold_std:
                            bonds.append((i, j))
                else:
                    if d_ij < 1.3 * (ri + rj):
                        bonds.append((i, j))

        # Every ILP dative M–L contact (bo=0) must appear in the CBC coordination graph.
        for tm, lig in dative_ml_pairs:
            bonds.append((tm, lig))

        tm_indices = [i for i in range(len(atoms)) if is_TM(atoms[i])]
        bond_set = set(map(tuple, [sorted(b) for b in bonds]))
        adj_tmp = defaultdict(set)
        for i, j in bonds:
            adj_tmp[i].add(j)
            adj_tmp[j].add(i)

        all_C_adj = defaultdict(set)
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                if atoms[i] == "C" and atoms[j] == "C":
                    ri2 = COV_R.get(atoms[i], 0.77)
                    rj2 = COV_R.get(atoms[j], 0.77)
                    if dist(coords[i], coords[j]) < 1.35 * (ri2 + rj2):
                        all_C_adj[i].add(j)
                        all_C_adj[j].add(i)

        for tm in tm_indices:
            tm_C_all = {j for j in adj_tmp[tm] if atoms[j] == "C"}
            tm_C_dative = {c for c in tm_C_all if bo_lewis.get((min(tm, c), max(tm, c)), 0) == 0}
            if not tm_C_dative:
                continue
            visited_c = set(tm_C_all)
            frontier = set(tm_C_dative)
            for _ in range(3):
                new_frontier = set()
                for c in frontier:
                    for nb in all_C_adj[c]:
                        if nb not in visited_c:
                            ri2 = COV_R.get(atoms[tm], 0.77)
                            rj2 = COV_R.get(atoms[nb], 0.77)
                            d_nb = dist(coords[tm], coords[nb])
                            if d_nb < 1.3 * (ri2 + rj2):
                                key = (min(tm, nb), max(tm, nb))
                                if key not in bond_set:
                                    bonds.append((tm, nb))
                                    bond_set.add(key)
                                    adj_tmp[tm].add(nb)
                                    adj_tmp[nb].add(tm)
                                visited_c.add(nb)
                                new_frontier.add(nb)
                frontier = new_frontier

        return bonds

    # M–L contacts with ILP dative bond (order 0): donor LP sits on the coordinate bond.
    dative_ml_pairs = set()
    if metal_adjacency_edges:
        for i, j, ei, ej in metal_adjacency_edges:
            if not is_TM(ei) ^ is_TM(ej):
                continue
            tm_i, lig_i = (i, j) if is_TM(ei) else (j, i)
            key = (min(tm_i, lig_i), max(tm_i, lig_i))
            if bo.get(key, 0) == 0:
                dative_ml_pairs.add((tm_i, lig_i))
    else:
        for i in range(n):
            if not is_TM(atoms[i]):
                continue
            for j in range(n):
                if i == j or is_TM(atoms[j]):
                    continue
                key = (min(i, j), max(i, j))
                if bo.get(key, 0) != 0:
                    continue
                ri = COV_R.get(atoms[i], 1.5)
                rj = COV_R.get(atoms[j], 0.77)
                if dist(coords[i], coords[j]) < 1.55 * (ri + rj):
                    dative_ml_pairs.add((i, j))

    lp = dict(lp) if lp else {}
    for _tm, lig in dative_ml_pairs:
        lp[lig] = lp.get(lig, 0) + 1

    all_bonds_list = _full_connectivity_dative(atoms, coords, bo, lp, dative_ml_pairs)
    adj_full = defaultdict(set)
    for i, j in all_bonds_list:
        adj_full[i].add(j)
        adj_full[j].add(i)

    adj_lewis = defaultdict(set)
    for (i, j) in bo:
        adj_lewis[i].add(j)
        adj_lewis[j].add(i)

    def bo_ij(i, j):
        return bo.get((min(i, j), max(i, j)), 0)

    def _pi_cluster(metal_idx, seed_idx):
        if atoms[seed_idx] not in ("C", "N"):
            return None
        metal_CN = {j for j in adj_full[metal_idx] if atoms[j] in ("C", "N")}
        if seed_idx not in metal_CN:
            return None
        visited = set()
        queue = [seed_idx]
        while queue:
            cur = queue.pop()
            if cur in visited:
                continue
            visited.add(cur)
            for nb in adj_full[cur]:
                if nb in metal_CN and atoms[nb] in ("C", "N") and nb not in visited:
                    queue.append(nb)
        return visited

    def _classify_pi(metal_idx, nbr_idx):
        if atoms[nbr_idx] not in ("C", "N"):
            return None
        cluster = _pi_cluster(metal_idx, nbr_idx)
        if cluster is None or len(cluster) <= 1:
            if bo_ij(nbr_idx, metal_idx) > 0:
                return None
            if atoms[nbr_idx] == "C":
                n_val_c = normal_bonds("C")
                bond_e_c = sum(bo_ij(nbr_idx, k) for k in adj_lewis[nbr_idx])
                lp_c = lp.get(nbr_idx, 0) if lp else 0
                if bond_e_c >= n_val_c and lp_c == 0:
                    return None
                if any(atoms[k] == "O" and bo_ij(nbr_idx, k) >= 2 for k in adj_lewis[nbr_idx]):
                    return ["L"]
                if any(bo_ij(nbr_idx, k) > 1 for k in adj_lewis[nbr_idx]):
                    return ["L"]
                return None
            if atoms[nbr_idx] == "N":
                fc_n = fc[nbr_idx] if fc else 0
                if fc_n < 0:
                    all_nbrs_are_N = all(atoms[k] == "N" for k in adj_lewis[nbr_idx])
                    if all_nbrs_are_N:
                        return None
                if any(bo_ij(nbr_idx, k) > 1 for k in adj_lewis[nbr_idx]):
                    return ["L"]
                return None
            return None

        dative_members = {j for j in cluster if bo_ij(j, metal_idx) == 0}
        if not dative_members:
            return None
        if nbr_idx != min(dative_members):
            return None

        size = len(cluster)

        def _is_x_type(j):
            fc_j = fc[j] if fc else 0
            if fc_j < 0:
                return True
            if bo_ij(j, metal_idx) > 0:
                return True
            sym_j = atoms[j]
            exp_be = 4 if sym_j == "C" else (3 if sym_j == "N" else 2)
            bond_e_j = sum(bo_ij(j, k) for k in adj_lewis[j])
            return bond_e_j < exp_be

        n_X_atoms = sum(1 for j in cluster if _is_x_type(j))
        if size == 2:
            a2, b2 = sorted(cluster)
            pi_order = max(0, bo_ij(a2, b2) - 1)
            if pi_order == 0:
                return None
            if n_X_atoms > 0:
                pi_types = ["L", "X"]
            elif pi_order >= 2:
                pi_types = ["L", "L"]
            else:
                pi_types = ["L"]
        elif size == 3:
            pi_types = ["L", "X"]
        elif size == 4:
            pi_types = ["L", "L", "X"] if n_X_atoms > 0 else ["L", "L"]
        elif size == 5:
            pi_types = ["L", "L", "X"]
        elif size == 6:
            pi_types = ["L", "L", "L"]
        else:
            n_X_forced = size % 2
            n_X_actual = max(n_X_forced, n_X_atoms)
            n_L = (size - n_X_actual) // 2
            pi_types = ["L"] * n_L + ["X"] * n_X_actual

        rep_lp = lp.get(nbr_idx, 0) if lp else 0
        rep_bond_to_M = bo_ij(nbr_idx, metal_idx)
        if rep_lp > 0 and rep_bond_to_M == 0 and "X" not in pi_types:
            pi_types = ["L"] + pi_types
        return pi_types

    def _fc_at(idx):
        if not fc or idx < 0 or idx >= len(fc):
            return 0
        return fc[idx]

    def _ilp_dative_sigma_ml(m_idx, lig_idx):
        """True when the Lewis bond table has no covalent M–L order (ILP dative / omitted M–L)."""
        return bo_ij(m_idx, lig_idx) == 0

    def _dative_ml_instead_of_x(m_idx, lig_idx, lp_lig):
        """ILP-style dative M←L: no M–L bond order but L can donate (lp or anionic)."""
        if not _ilp_dative_sigma_ml(m_idx, lig_idx):
            return False
        return lp_lig > 0 or _fc_at(lig_idx) < 0

    def classify_one(metal_idx, nbr_idx):
        sym_m = atoms[metal_idx]
        sym_n = atoms[nbr_idx]
        if sym_n == "H":
            h_nbrs = list(adj_lewis[nbr_idx])
            if any(is_TM(atoms[k]) for k in h_nbrs):
                return ["X"], None
            return ["X"], None
        if is_TM(sym_n):
            en_m = _OXSTATE_ENEG.get(sym_m, 2.0)
            en_n = _OXSTATE_ENEG.get(sym_n, 2.0)
            return (["L"] if en_n < en_m else ["Z"]), None
        if is_inorganic(sym_n):
            n_val = normal_bonds(sym_n)
            n_sub = sum(bo_ij(nbr_idx, k) for k in adj_lewis[nbr_idx] if k != metal_idx)
            lp_lewis = lp.get(nbr_idx, 0)
            deficit = n_val - n_sub
            if sym_n in ("Br", "I") and n_sub == 0:
                if _dative_ml_instead_of_x(metal_idx, nbr_idx, lp_lewis):
                    return ["L"], None
                return ["X"], None
            if sym_n == "I" and n_sub >= 2 and is_TM(sym_m):
                fc_n = fc[nbr_idx] if fc else 0
                if fc_n > 0:
                    return ["Z"], None
            if sym_n in ("B", "Al", "Ga", "In") and lp_lewis == 0:
                return ["Z"], None
            if deficit <= 0:
                return (["L"] if lp_lewis > 0 else ["Z"]), None
            if deficit == 1:
                if _dative_ml_instead_of_x(metal_idx, nbr_idx, lp_lewis):
                    return ["L"], None
                return ["X"], None
            if lp_lewis > 0:
                return (["L"] + ["X"] * (deficit - 1)), None
            return (["X"] * deficit), None

        pi = _classify_pi(metal_idx, nbr_idx)
        if pi is not None:
            return pi, _pi_cluster(metal_idx, nbr_idx)
        if atoms[nbr_idx] in ("C", "N"):
            cluster_check = _pi_cluster(metal_idx, nbr_idx)
            if cluster_check and len(cluster_check) > 1:
                dative_check = {j for j in cluster_check if bo_ij(j, metal_idx) == 0}
                if dative_check and nbr_idx != min(dative_check):
                    return [], None

        sym_n = atoms[nbr_idx]
        n_val = normal_bonds(sym_n)
        lp_lewis = lp.get(nbr_idx, 0)
        n_sub_organic = sum(bo_ij(nbr_idx, k) for k in adj_lewis[nbr_idx] if not is_inorganic(atoms[k]))
        n_sub_all = sum(bo_ij(nbr_idx, k) for k in adj_lewis[nbr_idx] if k != metal_idx and not is_TM(atoms[k]))
        deficit = n_val - n_sub_organic
        if sym_n == "C":
            if lp_lewis > 0:
                return ["L"], None
            has_CO = any(atoms[k] == "O" and bo_ij(nbr_idx, k) >= 2 for k in adj_lewis[nbr_idx])
            if has_CO:
                other_organic = sum(
                    1 for k in adj_lewis[nbr_idx]
                    if k != metal_idx and not is_inorganic(atoms[k]) and atoms[k] != "O"
                )
                if other_organic == 0:
                    return ["L"], None
            if n_sub_all <= 2:
                return ["L"], None

        if deficit <= 0:
            if lp_lewis > 0:
                return ["L"], None
            return [], None
        if deficit == 1:
            if _dative_ml_instead_of_x(metal_idx, nbr_idx, lp_lewis):
                return ["L"], None
            return ["X"], None

        n_sub_all_total = sum(bo_ij(nbr_idx, k) for k in adj_lewis[nbr_idx] if not is_TM(atoms[k]))
        deficit_real = n_val - n_sub_all_total
        if n_sub_all_total == 0 and sym_n in ("N", "O", "S", "C"):
            return ["X"] * deficit, None
        if deficit_real > 1 and sym_n in ("N", "O", "S", "C"):
            return ["X"] * min(deficit_real, deficit), None
        return ["X"], None

    def _fix_symmetric_chelate_O(metal_idx, nbr_cls, atoms, bo, adj_lewis, adj_full):
        o_L = [
            idx for idx, types in nbr_cls.items()
            if atoms[idx] == "O" and types == ["L"] and bo.get((min(idx, metal_idx), max(idx, metal_idx)), 0) == 0
        ]
        if len(o_L) >= 2:
            for i in range(len(o_L)):
                for j in range(i + 1, len(o_L)):
                    oi, oj = o_L[i], o_L[j]
                    ci = next((k for k in adj_lewis[oi] if atoms[k] == "C" and bo.get((min(oi, k), max(oi, k)), 0) >= 2), None)
                    cj = next((k for k in adj_lewis[oj] if atoms[k] == "C" and bo.get((min(oj, k), max(oj, k)), 0) >= 2), None)
                    if ci is None or cj is None:
                        continue
                    connected = ci == cj or (cj in adj_lewis[ci]) or (ci in adj_lewis[cj])
                    if not connected:
                        for mid in adj_lewis[ci]:
                            if cj in adj_lewis[mid]:
                                connected = True
                                break
                    if connected:
                        nbr_cls[max(oi, oj)] = ["X"]
                        break

    def _expand_cluster(rep_idx, cluster_types, cluster_atoms_sorted, metal_idx):
        records = []
        types_work = list(cluster_types)
        n_atoms = len(cluster_atoms_sorted)
        n_L = types_work.count("L")
        n_L_cluster = n_atoms // 2
        n_L_solo = n_L - n_L_cluster
        for _ in range(n_L_solo):
            records.append(((rep_idx,), "L"))

        cluster_set = set(cluster_atoms_sorted)

        def _is_x_candidate(idx):
            fc_val = fc[idx] if fc else 0
            if fc_val < 0:
                return True
            sym = atoms[idx]
            exp_be = 4 if sym == "C" else (3 if sym == "N" else 2)
            bond_e = sum(bo_ij(idx, k) for k in adj_lewis[idx])
            return bond_e < exp_be

        anionic = {a for a in cluster_set if _is_x_candidate(a)}
        neutral = [a for a in cluster_atoms_sorted if a not in anionic]
        double_pairs = []
        used = set()
        bond_candidates = []
        for a in neutral:
            for b in adj_lewis[a]:
                if b in cluster_set and b not in anionic and b > a:
                    order = bo_ij(a, b)
                    bond_candidates.append((-order, a, b))
        bond_candidates.sort()
        for _, a, b in bond_candidates:
            if a not in used and b not in used and len(double_pairs) < n_L_cluster:
                double_pairs.append((a, b))
                used.add(a)
                used.add(b)
        remaining = [a for a in neutral if a not in used]
        while len(double_pairs) < n_L_cluster and len(remaining) >= 2:
            a = remaining.pop(0)
            partner = next((b for b in remaining if b in adj_lewis[a]), None)
            if partner is None:
                partner = remaining[0]
            double_pairs.append((a, partner))
            remaining.remove(partner)
            used.add(a)
            used.add(partner)
        for a, b in double_pairs:
            records.append(((a, b), "L"))
        x_atoms = sorted(anionic) + [a for a in cluster_atoms_sorted if a not in used and a not in anionic]
        for ax in x_atoms:
            records.append(((ax,), "X"))
        return records

    results = {}
    for metal_idx in range(n):
        sym_m = atoms[metal_idx]
        if not is_inorganic(sym_m):
            continue
        if sym_m in ("Br", "I"):
            n_sub_m = sum(bo_ij(metal_idx, k) for k in adj_lewis[metal_idx])
            tm_nbrs = [j for j in adj_full[metal_idx] if is_TM(atoms[j])]
            if n_sub_m == 0 and tm_nbrs:
                continue
        neighbours = sorted(adj_full[metal_idx])
        if not neighbours:
            continue
        nbr_cls = {}
        nbr_cluster = {}
        for nbr_idx in neighbours:
            cbc, cluster = classify_one(metal_idx, nbr_idx)
            if cbc:
                nbr_cls[nbr_idx] = cbc
                if cluster and len(cluster) > 1:
                    nbr_cluster[nbr_idx] = cluster
        _fix_symmetric_chelate_O(metal_idx, nbr_cls, atoms, bo, adj_lewis, adj_full)
        if not nbr_cls:
            continue
        interaction_records = []
        for nbr_idx in sorted(nbr_cls.keys()):
            types = nbr_cls[nbr_idx]
            cluster = nbr_cluster.get(nbr_idx)
            if cluster and len(cluster) > 1:
                cluster_sorted = sorted(cluster)
                records = _expand_cluster(nbr_idx, types, cluster_sorted, metal_idx)
                interaction_records.extend(records)
                continue
            if len(types) == 1:
                interaction_records.append(((nbr_idx,), types[0]))
            elif len(set(types)) == 1 and len(types) > 1:
                for t in types:
                    interaction_records.append(((nbr_idx,), t))
            else:
                cluster2 = _pi_cluster(metal_idx, nbr_idx)
                if cluster2 and len(cluster2) > 1 and nbr_idx == min(cluster2):
                    cluster_sorted = sorted(cluster2)
                    records = _expand_cluster(nbr_idx, types, cluster_sorted, metal_idx)
                    interaction_records.extend(records)
                else:
                    for t in types:
                        interaction_records.append(((nbr_idx,), t))
        results[metal_idx] = interaction_records

    return results


def print_cbc_report(atoms, coords, bo, lp, fc, charge=0, *, metal_adjacency_edges=None):
    """Print the CBC classification table for every inorganic centre."""
    results = classify_cbc_ligands(
        atoms, coords, bo, lp, fc, charge, metal_adjacency_edges=metal_adjacency_edges
    )

    if not results:
        print("\n  (No inorganic atoms found — no CBC classification.)")
        return

    print()
    print("=" * 72)
    print("  CBC Ligand Classification  (L = dative, X = covalent, Z = Lewis-acid)")
    print("=" * 72)

    adj_lewis = defaultdict(set)
    for (i, j) in bo:
        adj_lewis[i].add(j)
        adj_lewis[j].add(i)

    def bo_ij(i, j):
        return bo.get((min(i, j), max(i, j)), 0)

    def avg_dist(metal_idx, atom_tuple):
        if len(atom_tuple) == 1:
            a = atom_tuple[0]
            return math.sqrt(sum((coords[metal_idx][k] - coords[a][k]) ** 2 for k in range(3)))
        a, b = atom_tuple
        if atoms[b] == "H" and atoms[a] in ("B", "C", "Si", "Al", "Ga"):
            return math.sqrt(sum((coords[metal_idx][k] - coords[b][k]) ** 2 for k in range(3)))
        dists = [
            math.sqrt(sum((coords[metal_idx][k] - coords[x][k]) ** 2 for k in range(3)))
            for x in atom_tuple
        ]
        return sum(dists) / len(dists)

    def row_label(atom_tuple):
        parts = [f"{atoms[a]}{a+1}" for a in atom_tuple]
        return "-".join(parts)

    def row_note(atom_tuple, cbc_char, metal_idx):
        if cbc_char == "Z":
            return "Lewis acid — metal donates electrons to ligand"
        rep = atom_tuple[0]
        sym = atoms[rep]
        n_val = _NORMAL_BONDS.get(sym, 4)
        n_sub = sum(bo_ij(rep, k) for k in adj_lewis[rep] if not is_inorganic(atoms[k]))
        lp_cnt = lp.get(rep, 0) if lp else 0

        if len(atom_tuple) == 2:
            a, b = atom_tuple
            if atoms[a] == "H" and atoms[b] == "H":
                return "η²-H₂ dihydrogen — σ(H-H) dative donor"
            if atoms[b] == "H" and atoms[a] in ("B", "C", "Si", "Al", "Ga"):
                return f"agostic {atoms[a]}-H σ-bond → σ-complex dative donation to M"
            bond_type = "C=C" if (atoms[a] == "C" and atoms[b] == "C") else (
                "N=N" if (atoms[a] == "N" and atoms[b] == "N") else "π"
            )
            return f"η² {bond_type} π-bond → dative π donation to M"

        if sym == "H":
            h_lewis_nbrs = list(adj_lewis[rep])
            if cbc_char == "X":
                if h_lewis_nbrs:
                    return "hydride H — covalent M-H bond"
                return "hydride H — covalent M-H bond"
            if cbc_char == "L":
                if len(h_lewis_nbrs) == 1:
                    parent = h_lewis_nbrs[0]
                    parent_sym = atoms[parent]
                    return f"agostic {parent_sym}-H → σ-complex dative donation to M"
                for other_h in adj_lewis[rep]:
                    if atoms[other_h] == "H":
                        return "η²-H₂ dihydrogen — σ(H-H) dative donor"
                return "H σ-complex — dative donor"

        if cbc_char == "L":
            is_co = (sym == "C" and any(atoms[k] == "O" and bo_ij(rep, k) >= 2 for k in adj_lewis[rep]))
            if is_co:
                return "CO — σ-donor via C lone pair"
            if sym in ("C", "N") and n_sub <= 2:
                return f"{sym} carbene/NHC — lone-pair σ-donor"
            if sym == "O" and n_sub < n_val and lp_cnt > 0:
                fc_rep = fc[rep] if fc and rep < len(fc) else 0
                if fc_rep < 0:
                    return f"O⁻ (oxyanion) — lone-pair dative donor ({lp_cnt} LP)"
                return f"O=N/S (resonance) — lone-pair dative donor ({lp_cnt} LP)"
            if lp_cnt > 0:
                return f"{sym}: {n_sub}/{n_val} bonds → lone-pair donor ({lp_cnt} LP)"
            return f"{sym}: {n_sub}/{n_val} bonds → dative donor"

        if cbc_char == "X":
            n_sub_total = sum(bo_ij(rep, k) for k in adj_lewis[rep] if not is_TM(atoms[k]))
            deficit_real = n_val - n_sub_total
            if deficit_real <= 0:
                return f"{sym}: {n_sub_total}/{n_val} bonds + covalent σ bond to M"
            return f"{sym}: {n_sub_total}/{n_val} bonds → {deficit_real} short → covalent σ to M"
        return ""

    for metal_idx in sorted(results):
        sym_m = atoms[metal_idx]
        interaction_records = results[metal_idx]

        total_L = sum(1 for _, t in interaction_records if t == "L")
        total_X = sum(1 for _, t in interaction_records if t == "X")
        total_Z = sum(1 for _, t in interaction_records if t == "Z")

        parts = []
        if total_L:
            parts.append(f"L{_subscript(total_L)}")
        if total_X:
            parts.append(f"X{_subscript(total_X)}")
        if total_Z:
            parts.append(f"Z{_subscript(total_Z)}")
        designation = "".join(parts) if parts else "—"

        print(f"\n  {sym_m}{metal_idx+1}   [{designation}]")
        print(f"  {'Neighbour':<16} {'Avg.Dist':>9}  {'CBC':<5}  Notes")
        print(f"  {'-'*65}")

        seen_multi = {}
        collapsed = []
        for atom_tuple, cbc_char in interaction_records:
            key = (atom_tuple, cbc_char)
            if key in seen_multi:
                seen_multi[key] += 1
            else:
                seen_multi[key] = 1
                collapsed.append((atom_tuple, cbc_char))

        for atom_tuple, cbc_char in collapsed:
            key = (atom_tuple, cbc_char)
            count = seen_multi[key]
            lbl = row_label(atom_tuple)
            if count > 1:
                lbl = f"{lbl}(×{count})"
            d_str = f"{avg_dist(metal_idx, atom_tuple):.2f} Å"
            note = row_note(atom_tuple, cbc_char, metal_idx)
            if count > 1:
                sym_rep = atoms[atom_tuple[0]]
                bond_names = {2: "double (2×)", 3: "triple (3×)"}
                bond_desc = bond_names.get(count, f"{count}×")
                note = f"{sym_rep}: {bond_desc} covalent M={sym_rep} bond (nitrido/oxo/imido)"
            print(f"  {lbl:<16} {d_str:>9}  {cbc_char:<5}  {note}")

        print(f"\n  Total:  {total_L}L + {total_X}X + {total_Z}Z  →  {designation}")


def _subscript(n):
    subs = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    return str(n).translate(subs) if n > 1 else ""


def main():
    if len(sys.argv) not in (2, 3):
        raise SystemExit(
            "Usage: python ilp_bond_order_example.py <molecule.xyz> [molecular_charge]"
        )

    xyz_path = sys.argv[1]
    mol_charge = int(sys.argv[2]) if len(sys.argv) == 3 else 0
    atoms = read_xyz(xyz_path)
    raw_edges = connectivity(atoms)
    edges = remove_dative_tm_bonds(atoms, raw_edges)
    bonds, lp_out, fc_out = solve_bond_orders(atoms, edges, mol_charge=mol_charge)

    atom_syms = [a[1] for a in atoms]
    print(f"Read {len(atoms)} atoms from {xyz_path}  (charge={mol_charge})")
    print()
    print_summary_and_choose_ilp(atom_syms, bonds, lp_out, mol_charge)

    bo0 = {(i - 1, j - 1): o for i, j, o in bonds}
    lp_full = {i - 1: lp_out.get(i, 0) for i in range(1, len(atom_syms) + 1)}
    fc0 = [0] * len(atom_syms)
    for i in range(len(atom_syms)):
        fc0[i] = fc_out.get(i + 1, 0)
    print_octet_report(atom_syms, bo0, lp_full, fc0)

    lp0 = {i - 1: v for i, v in lp_out.items() if v > 0}
    coords = [[a[2], a[3], a[4]] for a in atoms]
    print_cbc_report(atom_syms, coords, bo0, lp0, fc0, mol_charge)


if __name__ == "__main__":
    main()
