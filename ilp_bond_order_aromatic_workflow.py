#!/usr/bin/env python3
"""
Aromatic-aware ILP workflow (new script, does not modify original backend).

Key ideas:
1) Build raw connectivity from XYZ.
2) Remove metal-nonmetal edges for ring/fused-ring detection.
3) Keep only planar ring/fused-ring systems as aromatic candidates.
4) Add a 4n+2 pi-electron mismatch penalty into ILP objective.
5) For O/N/S in aromatic systems:
   - if no double bond around that atom, one lone pair (2e) can contribute;
   - if double-bonded, pi contribution comes from double bonds (2 per double).
6) For ring carbons not bonded to a transition metal: if unsaturated (any incident
   multiple bond or formal charge <= -1), the same optional lone-pair pi term applies.
7) Penalize lone-pair count on carbons that are not bonded to any transition metal
   (prefer lp on M-bound C); weight W_C_LP_OFF_METAL (tunable).
"""

from __future__ import annotations

import itertools
import math
import sys
from collections import defaultdict

import ilp_bond_order_example as base

pulp = base.pulp

# Objective weight: penalize each lone-pair *unit* on C not adjacent to any TM.
W_C_LP_OFF_METAL = 100


def _remove_metal_nonmetal_edges(raw_edges):
    """For ring detection only: remove edges where exactly one endpoint is TM."""
    keep = []
    for i, j, ei, ej in raw_edges:
        if base.is_TM(ei) ^ base.is_TM(ej):
            continue
        keep.append((i, j, ei, ej))
    return keep


def _build_adj_from_edges(edges):
    adj = defaultdict(set)
    for i, j, *_ in edges:
        adj[i].add(j)
        adj[j].add(i)
    return adj


def _find_simple_rings(adj, atom_symbol, max_size=12):
    """Simple cycle finder (set-unique), adapted from Lewis-engine style."""
    nodes = sorted(adj.keys())
    rings = []
    seen = set()
    for start in nodes:
        if atom_symbol[start] == "H":
            continue
        stack = [(start, [start], {start})]
        while stack:
            node, path, visited = stack.pop()
            for nb in adj[node]:
                if atom_symbol[nb] == "H":
                    continue
                if nb == start and len(path) >= 4:
                    key = frozenset(path)
                    if key not in seen:
                        seen.add(key)
                        rings.append(sorted(path))
                elif nb not in visited and len(path) < max_size:
                    stack.append((nb, path + [nb], visited | {nb}))
    return [frozenset(r) for r in rings]


def _fused_groups(rings):
    """Group rings that share an edge (>=2 common atoms)."""
    rings = list(rings)
    if not rings:
        return []
    parent = list(range(len(rings)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i in range(len(rings)):
        for j in range(i + 1, len(rings)):
            if len(rings[i] & rings[j]) >= 2:
                union(i, j)

    g = defaultdict(list)
    for i in range(len(rings)):
        g[find(i)].append(i)
    groups = []
    for idxs in g.values():
        s = set()
        for k in idxs:
            s |= set(rings[k])
        groups.append(tuple(sorted(s)))
    return groups


def _minimal_rings(rings):
    """
    Keep only smallest rings in an inclusion sense:
    drop any ring that strictly contains another ring.
    """
    rings = [frozenset(r) for r in rings]
    out = []
    for r in rings:
        has_strict_subset = any((s < r) for s in rings if s is not r)
        if not has_strict_subset:
            out.append(r)
    # de-dup preserve stable order by (size, atoms)
    uniq = {}
    for r in sorted(out, key=lambda x: (len(x), tuple(sorted(x)))):
        uniq[r] = None
    return list(uniq.keys())


def _plane_rmsd(points):
    """
    Heuristic planarity check without external deps:
    try all point triplets as candidate planes and take minimal RMS distance.
    """
    if len(points) <= 3:
        return 0.0
    best = float("inf")
    for a, b, c in itertools.combinations(points, 3):
        ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
        vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
        nx, ny, nz = (
            uy * vz - uz * vy,
            uz * vx - ux * vz,
            ux * vy - uy * vx,
        )
        norm = math.sqrt(nx * nx + ny * ny + nz * nz)
        if norm < 1e-8:
            continue
        d0 = -(nx * a[0] + ny * a[1] + nz * a[2])
        sq = 0.0
        for p in points:
            dist = (nx * p[0] + ny * p[1] + nz * p[2] + d0) / norm
            sq += dist * dist
        rms = math.sqrt(sq / len(points))
        if rms < best:
            best = rms
    return best if best < float("inf") else 999.0


def aromatic_candidate_systems(atoms_packed, raw_edges, max_ring_size=12, planarity_rmsd=0.12):
    atom_symbol = {idx: el for idx, el, *_ in atoms_packed}
    atom_xyz = {idx: (x, y, z) for idx, _el, x, y, z in atoms_packed}
    cut_edges = _remove_metal_nonmetal_edges(raw_edges)
    adj = _build_adj_from_edges(cut_edges)
    rings_all = _find_simple_rings(adj, atom_symbol, max_size=max_ring_size)
    rings_min = _minimal_rings(rings_all)

    # Step 1: smallest planar rings
    planar_min_rings = []
    for r in rings_min:
        syms = [atom_symbol[i] for i in r]
        if any(base.is_TM(s) or s == "H" for s in syms):
            continue
        if sum(1 for s in syms if s in ("C", "N", "O", "S")) < 4:
            continue
        pts = [atom_xyz[i] for i in sorted(r)]
        if _plane_rmsd(pts) <= planarity_rmsd:
            planar_min_rings.append(frozenset(r))

    # Step 2: assemble larger fused systems from planar minimal rings
    groups = _fused_groups(planar_min_rings)
    out = []
    for g in groups:
        out.append(tuple(sorted(g)))
    return out


def solve_bond_orders(atoms, edges, aromatic_systems, mol_charge=0):
    prob = pulp.LpProblem("BondOrderAssignmentAromatic", pulp.LpMinimize)

    # Bond orders are defined by connectivity:
    # every candidate edge is at least SINGLE (1). The ILP may upgrade to DOUBLE (2)
    # or TRIPLE (3). No bond order 0 is allowed.
    u2 = {}  # edge is at least double
    u3 = {}  # edge is triple (implies u2)
    for i, j, _ei, _ej in edges:
        key = (i, j) if i < j else (j, i)
        u2[key] = pulp.LpVariable(f"u2_{key[0]}_{key[1]}", cat="Binary")
        u3[key] = pulp.LpVariable(f"u3_{key[0]}_{key[1]}", cat="Binary")
        prob += u3[key] <= u2[key]

    inc_edges = defaultdict(list)
    for i, j, *_ in edges:
        inc_edges[i].append((i, j))
        inc_edges[j].append((i, j))

    # Atoms that share an edge with any transition metal (for aromatic lp_pair on C).
    tm_neighbor_atoms = set()
    for i, j, ei, ej in edges:
        if base.is_TM(ei):
            tm_neighbor_atoms.add(j)
        if base.is_TM(ej):
            tm_neighbor_atoms.add(i)

    lp, oct_plus, oct_minus, q, abs_q, q_neg = {}, {}, {}, {}, {}, {}
    for i, el, *_ in atoms:
        if base.is_tm(el) or el not in base.VALENCE_ELECTRONS:
            continue
        q[i] = pulp.LpVariable(f"q_{i}", lowBound=-4, upBound=4, cat="Integer")
        abs_q[i] = pulp.LpVariable(f"absq_{i}", lowBound=0, cat="Integer")
        q_neg[i] = pulp.LpVariable(f"qneg_{i}", lowBound=0, cat="Integer")
        if base.is_ionlike_s_block_metal(el):
            continue
        lp[i] = pulp.LpVariable(f"lp_{i}", lowBound=0, cat="Integer")
        oct_plus[i] = pulp.LpVariable(f"octp_{i}", lowBound=0, cat="Integer")
        oct_minus[i] = pulp.LpVariable(f"octm_{i}", lowBound=0, cat="Integer")

    for i, el, *_ in atoms:
        if i not in q:
            continue
        bond_sum = pulp.lpSum(
            1 + u2[(min(a, b), max(a, b))] + u3[(min(a, b), max(a, b))]
            for a, b in inc_edges.get(i, [])
        ) if inc_edges.get(i) else 0
        ve = base.VALENCE_ELECTRONS[el]
        if base.is_ionlike_s_block_metal(el):
            prob += q[i] == ve - bond_sum
            prob += abs_q[i] >= q[i]
            prob += abs_q[i] >= -q[i]
            prob += q_neg[i] >= -q[i]
            continue
        oct_target = 2 if el in ("H", "Li") else 8
        local_e = 2 * lp[i] + 2 * bond_sum
        assigned_e = 2 * lp[i] + bond_sum
        prob += local_e - oct_target == oct_plus[i] - oct_minus[i]
        # Enforce octet/duet as a HARD constraint (no violations allowed).
        prob += oct_plus[i] == 0
        prob += oct_minus[i] == 0
        prob += q[i] == ve - assigned_e
        prob += abs_q[i] >= q[i]
        prob += abs_q[i] >= -q[i]
        prob += q_neg[i] >= -q[i]

    for i, j, ei, ej in edges:
        if "H" in (ei, ej) or "Cl" in (ei, ej):
            key = (min(i, j), max(i, j))
            prob += u2[key] == 0
            prob += u3[key] == 0
        # Fix all TM–nonmetal bonds as single (no upgrades).
        if base.is_TM(ei) ^ base.is_TM(ej):
            key = (min(i, j), max(i, j))
            prob += u2[key] == 0
            prob += u3[key] == 0

    # Triple bonds are rare; allow only in small-valence local environments.
    # Relaxed rule: allow u3 only for non-TM, non-H bonds where both endpoints
    # have heavy-degree <= 2 in this ILP edge set.
    heavy_deg = defaultdict(int)
    for i, j, ei, ej in edges:
        if ei != "H":
            heavy_deg[i] += 1
        if ej != "H":
            heavy_deg[j] += 1
    for i, j, ei, ej in edges:
        key = (min(i, j), max(i, j))
        if ei == "H" or ej == "H":
            prob += u3[key] == 0
            continue
        # Never allow triple bonds involving TM (TM–nonmetal upgrades are already fixed off).
        if base.is_TM(ei) or base.is_TM(ej):
            prob += u3[key] == 0
            continue
        # Avoid triple bonds on crowded atoms.
        if heavy_deg[i] > 2 or heavy_deg[j] > 2:
            prob += u3[key] == 0

    # --- Aromatic 4n+2 penalty ---
    atom_el = {i: el for i, el, *_ in atoms}
    aromatic_dev_terms = []
    for sys_idx, sys_atoms in enumerate(aromatic_systems):
        # pi_e = 2 * (#internal double-bond EDGES in system) + 2 * (hetero / ring-C LP contributions)
        # Count double bonds by EDGE (once per bond), not by per-atom incidence.
        sys_atom_set = set(sys_atoms)
        sys_edges = sorted(
            {
                (min(a, b), max(a, b))
                for i in sys_atoms
                for a, b in inc_edges.get(i, [])
                if (a in sys_atom_set and b in sys_atom_set)
            }
        )
        internal_double_edges = pulp.lpSum(
            u2[(a, b)] for a, b in sys_edges
        ) if sys_edges else 0
        internal_triple_edges = pulp.lpSum(
            u3[(a, b)] for a, b in sys_edges
        ) if sys_edges else 0

        double_present = {}
        lp_pair = {}
        for i in sys_atoms:
            dsum = pulp.lpSum(
                u2[(min(a, b), max(a, b))] for a, b in inc_edges.get(i, [])
                if (a in sys_atom_set and b in sys_atom_set)
            ) if inc_edges.get(i) else 0
            dp = pulp.LpVariable(f"sys{sys_idx}_dblp_{i}", cat="Binary")
            double_present[i] = dp
            deg = max(1, len([1 for a, b in inc_edges.get(i, []) if a in sys_atom_set and b in sys_atom_set]))
            prob += dsum >= dp
            prob += dsum <= deg * dp

            if atom_el[i] in ("O", "N", "S") and i in lp:
                lpb = pulp.LpVariable(f"sys{sys_idx}_lppair_{i}", cat="Binary")
                lp_pair[i] = lpb
                prob += lp[i] >= lpb
                prob += lpb <= 1 - dp
            elif atom_el[i] == "C" and i in lp and i not in tm_neighbor_atoms:
                # Like O/N/S: optional 2e from one LP toward pi, only if no internal
                # double from this atom in the system (dp), and only if "unsaturated":
                # any incident multiple bond in the full graph, or q <= -1.
                mult_inc = (
                    pulp.lpSum(
                        u2[(min(a, b), max(a, b))] + u3[(min(a, b), max(a, b))]
                        for a, b in inc_edges.get(i, [])
                    )
                    if inc_edges.get(i)
                    else 0
                )
                mult_any = pulp.LpVariable(f"sys{sys_idx}_multany_{i}", cat="Binary")
                max_mult = 2 * max(1, len(inc_edges.get(i, [])))
                prob += mult_inc >= mult_any
                prob += mult_inc <= max_mult * mult_any
                neg_c = pulp.LpVariable(f"sys{sys_idx}_negc_{i}", cat="Binary")
                prob += q[i] <= -1 + 5 * (1 - neg_c)
                prob += q[i] >= -4 * neg_c
                unsat = pulp.LpVariable(f"sys{sys_idx}_unsat_{i}", cat="Binary")
                prob += mult_any + neg_c <= 2 * unsat
                prob += mult_any + neg_c >= unsat
                lpb = pulp.LpVariable(f"sys{sys_idx}_lppair_{i}", cat="Binary")
                lp_pair[i] = lpb
                prob += lp[i] >= lpb
                prob += lpb <= 1 - dp
                prob += lpb <= unsat
            else:
                lp_pair[i] = 0

        # A double bond contributes 2 pi electrons; a triple contributes 4 (two pi bonds).
        pi_e = (
            2 * (internal_double_edges - internal_triple_edges)
            + 4 * internal_triple_edges
            + pulp.lpSum(2 * lp_pair[i] for i in sys_atoms)
        )
        max_pi = max(2, 4 * len(sys_atoms))
        kmax = max(0, (max_pi - 2) // 4)
        k = pulp.LpVariable(f"sys{sys_idx}_k", lowBound=0, upBound=kmax, cat="Integer")
        dev_p = pulp.LpVariable(f"sys{sys_idx}_devp", lowBound=0, cat="Integer")
        dev_m = pulp.LpVariable(f"sys{sys_idx}_devm", lowBound=0, cat="Integer")
        prob += pi_e - (4 * k + 2) == dev_p - dev_m
        aromatic_dev_terms.append(dev_p + dev_m)

    charge_delta_abs = pulp.LpVariable("charge_delta_abs", lowBound=0, cat="Integer")
    total_q = pulp.lpSum(q.values()) if q else 0
    total_abs_q = pulp.lpSum(abs_q.values()) if abs_q else 0
    prob += charge_delta_abs >= total_q - mol_charge
    prob += charge_delta_abs >= -(total_q - mol_charge)

    absq_match_abs = pulp.LpVariable("absq_match_abs", lowBound=0, cat="Integer")
    prob += absq_match_abs >= total_abs_q - abs(mol_charge)
    prob += absq_match_abs >= abs(mol_charge) - total_abs_q

    octet_penalty = pulp.lpSum(oct_plus.values()) + pulp.lpSum(oct_minus.values())
    formal_charge_penalty = total_abs_q
    double_bond_penalty = pulp.lpSum(
        u2[(min(i, j), max(i, j))] + u3[(min(i, j), max(i, j))] for i, j, *_ in edges
    )
    eneg_neg_charge_penalty = pulp.lpSum(
        max(0.0, 4.0 - base.ENEG.get(el, 2.0)) * q_neg[i]
        for i, el, *_ in atoms if i in q_neg
    )
    aromatic_penalty = pulp.lpSum(aromatic_dev_terms) if aromatic_dev_terms else 0

    c_lp_off_metal_terms = [
        lp[i]
        for i, el, *_ in atoms
        if el == "C" and i in lp and i not in tm_neighbor_atoms
    ]
    c_lp_off_metal_penalty = (
        pulp.lpSum(c_lp_off_metal_terms) if c_lp_off_metal_terms else 0
    )

    prob += (
        100 * formal_charge_penalty
        + 1 * charge_delta_abs
        + 1 * absq_match_abs
        + 100 * double_bond_penalty
        + 10 * eneg_neg_charge_penalty
        + 100 * aromatic_penalty
        + 100 * c_lp_off_metal_penalty
    )

    status = prob.solve(pulp.PULP_CBC_CMD(msg=False))
    if pulp.LpStatus[status] not in ("Optimal", "Integer Feasible"):
        raise RuntimeError(f"ILP failed: {pulp.LpStatus[status]}")

    bonds = []
    for i, j, _ei, _ej in edges:
        key = (min(i, j), max(i, j))
        order = 1 + int(round(pulp.value(u2[key]))) + int(round(pulp.value(u3[key])))
        bonds.append((i, j, order))
    bonds.sort()
    lp_out = {i: int(round(pulp.value(v))) for i, v in lp.items()}
    fc_out = {i: int(round(pulp.value(v))) for i, v in q.items()}
    return bonds, lp_out, fc_out


def main():
    if len(sys.argv) not in (2, 3):
        raise SystemExit(
            "Usage: python ilp_bond_order_aromatic_workflow.py <molecule.xyz> [molecular_charge]"
        )
    xyz = sys.argv[1]
    charge = int(sys.argv[2]) if len(sys.argv) == 3 else 0
    atoms = base.read_xyz(xyz)
    raw = base.connectivity(atoms)
    edges = base.remove_dative_tm_bonds(atoms, raw)
    aromatic_systems = aromatic_candidate_systems(atoms, raw)

    bonds, lp_out, fc_out = solve_bond_orders(atoms, edges, aromatic_systems, mol_charge=charge)
    print(f"Read {len(atoms)} atoms from {xyz} (charge={charge})")
    print(f"Aromatic candidate systems (assembled from minimal planar rings): {len(aromatic_systems)}")
    for i, s in enumerate(aromatic_systems, start=1):
        print(f"  System {i}: atoms {list(s)}")
    print()
    base.print_choose_block(bonds, lp_out)

    atom_syms = [a[1] for a in atoms]
    coords = [[a[2], a[3], a[4]] for a in atoms]
    bo0 = {(i - 1, j - 1): o for i, j, o in bonds}
    lp0 = {i - 1: v for i, v in lp_out.items() if v > 0}
    fc0 = [0] * len(atom_syms)
    for i in range(len(atom_syms)):
        fc0[i] = fc_out.get(i + 1, 0)
    base.print_cbc_report(atom_syms, coords, bo0, lp0, fc0, charge)


if __name__ == "__main__":
    main()

