#!/usr/bin/env python3
"""
Aromatic-aware ILP workflow V2: TM–nonmetal bonds are optimized inside the ILP

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
7) Hard: carbons not bonded to any TM have lp = 0 (lone pairs only on M-bound C).
8) mol_charge matching (hard):
   Σfc(ligand atoms) + Σox(TM) − Σ(covalent M–L bond orders) = mol_charge,
   where each ox(TM) is chosen from TM_COMMON_OXIDATION_STATES in
   ilp_bond_order_example.py (ILP binary choice when several states are listed).
9) Hard: Σox(TM) ≥ Σ(covalent M–L bond orders).

ILP uses full XYZ connectivity (raw): every TM–nonmetal contact is optimized with
dative (d=1, bond order 0) or covalent order 1/2/3 (1+u2+u3), except monatomic
Cl/Br/I/F/H ligands — those are fixed as covalent single bonds (d=0, u2=u3=0).
Each dative M–L bond consumes one LP on the ligand atom (2e in octet / fc accounting;
lp[i] counts only non-donated lone pairs).
Ring detection still strips M–L edges before cycle search.
Non–TM–nonmetal edges: u2/u3 only (order = 1 + u2 + u3).
$CHOOSE omits dative (order 0) M–L edges; covalent M–L prints as S/D/T.
"""

from __future__ import annotations

import itertools
import math
import sys
from collections import defaultdict

import ilp_bond_order_example as base

pulp = base.pulp

# Monatomic ligands forced to covalent M–X single bond in the ILP (MLX X-type).
TM_MONATOMIC_COV_LIGANDS = frozenset({"F", "Cl", "Br", "I", "H"})


def _remove_metal_nonmetal_edges(raw_edges):
    """For ring detection only: remove edges where exactly one endpoint is TM."""
    keep = []
    for i, j, ei, ej in raw_edges:
        if base.is_TM(ei) ^ base.is_TM(ej):
            continue
        keep.append((i, j, ei, ej))
    return keep


def _is_tm_nm_edge(ei, ej):
    """True iff exactly one endpoint is a transition metal."""
    return base.is_TM(ei) ^ base.is_TM(ej)


def _tm_nm_orient(i, j, ei, ej):
    """
    For a TM–nonmetal edge, return (tm_idx, lig_idx).
    Precondition: _is_tm_nm_edge(ei, ej).
    """
    if base.is_TM(ei) and not base.is_TM(ej):
        return i, j
    return j, i


def _atoms_bonded_to_tm(edge_list):
    """All atom indices with at least one TM neighbour in edge_list."""
    out = set()
    for i, j, ei, ej in edge_list:
        if base.is_TM(ei):
            out.add(j)
        if base.is_TM(ej):
            out.add(i)
    return out


def _tm_oxidation_ilp_vars(prob, atoms):
    """
    Per-TM oxidation-state variable constrained to TM_COMMON_OXIDATION_STATES[sym].
    Returns {tm_index: int constant or LpVariable}.
    """
    out = {}
    common = base.TM_COMMON_OXIDATION_STATES
    for i, el, *_ in atoms:
        if not base.is_TM(el):
            continue
        allowed = common.get(el)
        if not allowed:
            out[i] = pulp.LpVariable(f"ox_{i}", lowBound=0, upBound=8, cat="Integer")
            continue
        states = sorted({int(s) for s in allowed})
        if len(states) == 1:
            out[i] = states[0]
            continue
        picks = []
        for s in states:
            y = pulp.LpVariable(f"ox_{i}_{s}", cat="Binary")
            picks.append((y, s))
        prob += pulp.lpSum(y for y, _ in picks) == 1
        ox = pulp.LpVariable(f"ox_{i}", cat="Integer")
        prob += ox == pulp.lpSum(s * y for y, s in picks)
        out[i] = ox
    return out


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


def solve_bond_orders(
    atoms, edges, aromatic_systems, mol_charge=0, *, metal_adjacency_edges=None
):
    prob = pulp.LpProblem("BondOrderAssignmentAromatic", pulp.LpMinimize)

    atom_el = {i: el for i, el, *_ in atoms}

    ilp_edges = [t for t in edges if not _is_tm_nm_edge(t[2], t[3])]
    tm_nm_keys = []
    seen_tm_nm = set()
    for i, j, ei, ej in edges:
        if not _is_tm_nm_edge(ei, ej):
            continue
        tm, lig = _tm_nm_orient(i, j, ei, ej)
        key = (tm, lig)
        if key not in seen_tm_nm:
            seen_tm_nm.add(key)
            tm_nm_keys.append(key)

    u2 = {}
    u3 = {}
    for i, j, _ei, _ej in ilp_edges:
        key = (i, j) if i < j else (j, i)
        u2[key] = pulp.LpVariable(f"u2_{key[0]}_{key[1]}", cat="Binary")
        u3[key] = pulp.LpVariable(f"u3_{key[0]}_{key[1]}", cat="Binary")
        prob += u3[key] <= u2[key]

    # TM–nonmetal: same u2/u3 semantics as other edges (order 1+u2+u3), plus
    # binary d = coordinate/dative (ligand-side bond order 0; forces u2=u3=0).
    tm_dov = {}
    u2_tm = {}
    u3_tm = {}
    b_tm = {}
    for tm, lig in tm_nm_keys:
        d = pulp.LpVariable(f"tmd_{tm}_{lig}", cat="Binary")
        u2t = pulp.LpVariable(f"u2tm_{tm}_{lig}", cat="Binary")
        u3t = pulp.LpVariable(f"u3tm_{tm}_{lig}", cat="Binary")
        prob += u3t <= u2t
        prob += u2t <= 1 - d
        prob += u3t <= 1 - d
        tm_dov[(tm, lig)] = d
        u2_tm[(tm, lig)] = u2t
        u3_tm[(tm, lig)] = u3t
        s = 1 + u2t + u3t
        b = pulp.LpVariable(f"btm_{tm}_{lig}", lowBound=0, upBound=3, cat="Integer")
        prob += b <= s
        prob += b <= 3 * (1 - d)
        prob += b >= s - 3 * d
        prob += b >= 0
        b_tm[(tm, lig)] = b

    # Per ligand: count of dative M–L contacts (each uses one LP = 2e on donor).
    n_dative_ml = {}
    dative_by_lig = defaultdict(list)
    for tm, lig in tm_nm_keys:
        dative_by_lig[lig].append(tm_dov[(tm, lig)])
    for lig, ds in dative_by_lig.items():
        n_dative_ml[lig] = pulp.lpSum(ds)

    inc_edges = defaultdict(list)
    for i, j, *_ in edges:
        inc_edges[i].append((i, j))
        inc_edges[j].append((i, j))

    def _incident_bond_order_sum(idx):
        terms = []
        for a, b in inc_edges.get(idx, []):
            ea, eb = atom_el[a], atom_el[b]
            if _is_tm_nm_edge(ea, eb):
                tm, lig = _tm_nm_orient(a, b, ea, eb)
                if lig != idx:
                    continue
                terms.append(b_tm[(tm, lig)])
            else:
                kk = (min(a, b), max(a, b))
                terms.append(1 + u2[kk] + u3[kk])
        return pulp.lpSum(terms) if terms else 0

    # Atoms bonded to a TM in the full structure (e.g. Cl–Au), for aromatic C rules
    # and for excluding M-bound ligand atoms from mol_charge summation.
    tm_neighbor_atoms = _atoms_bonded_to_tm(
        metal_adjacency_edges if metal_adjacency_edges is not None else edges
    )

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
        bond_sum = _incident_bond_order_sum(i)
        ve = base.VALENCE_ELECTRONS[el]
        if base.is_ionlike_s_block_metal(el):
            prob += q[i] == ve - bond_sum
            prob += abs_q[i] >= q[i]
            prob += abs_q[i] >= -q[i]
            prob += q_neg[i] >= -q[i]
            continue
        oct_target = 2 if el in ("H", "Li") else 8
        dative_pairs = n_dative_ml[i] if i in n_dative_ml else 0
        # Remaining LPs (lp) + Lewis bonds + donated pairs (2e each dative) fill octet.
        local_e = 2 * lp[i] + 2 * bond_sum + 2 * dative_pairs
        assigned_e = 2 * lp[i] + bond_sum + 2 * dative_pairs
        prob += local_e - oct_target == oct_plus[i] - oct_minus[i]
        prob += oct_plus[i] == 0
        prob += oct_minus[i] == 0
        prob += q[i] == ve - assigned_e
        prob += abs_q[i] >= q[i]
        prob += abs_q[i] >= -q[i]
        prob += q_neg[i] >= -q[i]

    for i, j, ei, ej in ilp_edges:
        if "H" in (ei, ej) or "Cl" in (ei, ej):
            key = (min(i, j), max(i, j))
            prob += u2[key] == 0
            prob += u3[key] == 0

    heavy_deg = defaultdict(int)
    for i, j, ei, ej in edges:
        if ei != "H":
            heavy_deg[i] += 1
        if ej != "H":
            heavy_deg[j] += 1

    for i, j, ei, ej in ilp_edges:
        key = (min(i, j), max(i, j))
        if ei == "H" or ej == "H":
            prob += u3[key] == 0
            continue
        if base.is_TM(ei) or base.is_TM(ej):
            prob += u3[key] == 0
            continue
        if heavy_deg[i] > 2 or heavy_deg[j] > 2:
            prob += u3[key] == 0

    # Monatomic F/Cl/Br/I/H: fix covalent single bond (not dative); still in ILP for q/lp.
    for tm, lig in tm_nm_keys:
        if atom_el[lig] not in TM_MONATOMIC_COV_LIGANDS:
            continue
        prob += tm_dov[(tm, lig)] == 0
        prob += u2_tm[(tm, lig)] == 0
        prob += u3_tm[(tm, lig)] == 0

    # --- Aromatic 4n+2 penalty ---
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

    # Molecular charge: Σfc(ligands) + Σox(TM) − Σ(covalent M–L orders) = mol_charge.
    # ox(TM) ∈ TM_COMMON_OXIDATION_STATES; covalent order = b_tm (0 when dative).
    tm_ox_vars = _tm_oxidation_ilp_vars(prob, atoms)
    sigma_cov_ml = (
        pulp.lpSum(b_tm[k] for k in tm_nm_keys) if tm_nm_keys else 0
    )
    ligand_q_sum = pulp.lpSum(q[i] for i in q) if q else 0
    tm_ox_sum = pulp.lpSum(
        v if not isinstance(v, int) else v for v in tm_ox_vars.values()
    )
    prob += ligand_q_sum + tm_ox_sum - sigma_cov_ml == mol_charge

    if tm_ox_vars:
        prob += tm_ox_sum >= sigma_cov_ml

    total_abs_q_all = pulp.lpSum(abs_q.values()) if abs_q else 0
    absq_match_abs = pulp.LpVariable("absq_match_abs", lowBound=0, cat="Integer")
    prob += absq_match_abs >= total_abs_q_all - abs(mol_charge)
    prob += absq_match_abs >= abs(mol_charge) - total_abs_q_all

    octet_penalty = pulp.lpSum(oct_plus.values()) + pulp.lpSum(oct_minus.values())
    formal_charge_penalty = pulp.lpSum(abs_q.values()) if abs_q else 0
    tm_nm_high_order = (
        pulp.lpSum(u2_tm[k] + u3_tm[k] for k in tm_nm_keys) if tm_nm_keys else 0
    )
    double_bond_penalty = pulp.lpSum(
        u2[(min(i, j), max(i, j))] + u3[(min(i, j), max(i, j))] for i, j, *_ in ilp_edges
    ) + tm_nm_high_order
    eneg_neg_charge_penalty = pulp.lpSum(
        max(0.0, 4.0 - base.ENEG.get(el, 2.0)) * q_neg[i]
        for i, el, *_ in atoms if i in q_neg
    )
    aromatic_penalty = pulp.lpSum(aromatic_dev_terms) if aromatic_dev_terms else 0

    for i, el, *_ in atoms:
        if el == "C" and i in lp and i not in tm_neighbor_atoms:
            prob += lp[i] == 0

    prob += (
        100 * formal_charge_penalty
        + 100 * absq_match_abs
        # + 100 * double_bond_penalty
        + 10 * eneg_neg_charge_penalty
        + 100 * aromatic_penalty
    )

    status = prob.solve(pulp.PULP_CBC_CMD(msg=False))
    if pulp.LpStatus[status] not in ("Optimal", "Integer Feasible"):
        raise RuntimeError(
            f"ILP failed: {pulp.LpStatus[status]} "
            "(check mol_charge, connectivity, and dative/LP balance on TM donors)"
        )

    bonds = []
    for i, j, _ei, _ej in ilp_edges:
        key = (min(i, j), max(i, j))
        order = 1 + int(round(pulp.value(u2[key]))) + int(round(pulp.value(u3[key])))
        bonds.append((key[0], key[1], order))
    for tm, lig in tm_nm_keys:
        d = int(round(pulp.value(tm_dov[(tm, lig)])))
        if d == 1:
            continue
        u2t = int(round(pulp.value(u2_tm[(tm, lig)])))
        u3t = int(round(pulp.value(u3_tm[(tm, lig)])))
        order = 1 + u2t + u3t
        a, b = (tm, lig) if tm < lig else (lig, tm)
        bonds.append((a, b, order))
    bonds.sort()
    lp_out = {i: int(round(pulp.value(v))) for i, v in lp.items()}
    fc_out = {i: int(round(pulp.value(v))) for i, v in q.items()}
    for tm_i, ox_v in tm_ox_vars.items():
        if isinstance(ox_v, int):
            fc_out[tm_i] = ox_v
        else:
            fc_out[tm_i] = int(round(pulp.value(ox_v)))
    return bonds, lp_out, fc_out


def tm_oxidation_from_charge_balance(
    atom_syms, fc, mol_charge: int = 0, bo=None
) -> dict[int, tuple[int, int, int]]:
    """
    For each TM (post-ILP check, same identity enforced in the ILP):
      ox(M) = mol_charge − Σfc(ligand atoms) + Σ(covalent M–L bond orders in *bo*).

    Returns {tm_idx: (ox_state, sum_fc_ligands, sum_sigma)}.
    """
    n = len(atom_syms)
    if len(fc) != n:
        raise ValueError(f"fc length {len(fc)} != atom count {n}")
    out: dict[int, tuple[int, int, int]] = {}
    for m, sym in enumerate(atom_syms):
        if not base.is_TM(sym):
            continue
        ligands_sum = sum(
            fc[j] for j in range(n) if j != m and not base.is_TM(atom_syms[j])
        )
        sigma_sum = 0
        if bo is not None:
            for lig in range(n):
                if lig == m or base.is_TM(atom_syms[lig]):
                    continue
                sigma_sum += bo.get((min(m, lig), max(m, lig)), 0)
        ox = int(mol_charge) - int(ligands_sum) + int(sigma_sum)
        out[m] = (ox, int(ligands_sum), int(sigma_sum))
    return out


def print_tm_oxidation_sigma_report(
    atom_syms, fc, mol_charge: int = 0, bo=None
) -> None:
    """Report TM oxidation from charge / σ balance vs TM_COMMON_OXIDATION_STATES."""
    common = base.TM_COMMON_OXIDATION_STATES
    results = tm_oxidation_from_charge_balance(atom_syms, fc, mol_charge, bo=bo)
    if not results:
        return

    print()
    print("=" * 72)
    print("  TM Oxidation State (from formal-charge balance)")
    print("=" * 72)
    print(
        "  Rule (ILP + check): mol_charge = Σfc(ligands) + ox(M)"
        " − Σ(covalent M–L bond orders); ox(M) from TM_COMMON_OXIDATION_STATES."
    )
    print()

    for m in sorted(results):
        sym = atom_syms[m]
        ox_state, ligands_sum, sigma_sum = results[m]
        label = f"{sym}{m + 1}"

        ref = common.get(sym)
        if ref is None:
            status = "—  (no TM_COMMON_OXIDATION_STATES entry)"
        elif ox_state in ref:
            status = f"√  (common for {sym}: {ref})"
        else:
            status = f"⚠ WARNING: oxidation state {ox_state} not in common oxidation states {ref}"

        print(
            f"  {label}:  oxidation state = {ox_state}   "
            f"[{mol_charge} − ({ligands_sum}) + σ={sigma_sum}]"
        )
        print(f"           {status}")
        print()


def main():
    if len(sys.argv) not in (2, 3):
        raise SystemExit(
            "Usage: python ilp_bond_order_aromatic_workflow_V2.py <molecule.xyz> [molecular_charge]"
        )
    xyz = sys.argv[1]
    charge = int(sys.argv[2]) if len(sys.argv) == 3 else 0
    atoms = base.read_xyz(xyz)
    raw = base.connectivity(atoms)
    aromatic_systems = aromatic_candidate_systems(atoms, raw)

    bonds, lp_out, fc_out = solve_bond_orders(
        atoms, raw, aromatic_systems, mol_charge=charge, metal_adjacency_edges=raw
    )
    atom_syms = [a[1] for a in atoms]
    print(f"Read {len(atoms)} atoms from {xyz} (charge={charge})")
    print(f"Aromatic candidate systems (assembled from minimal planar rings): {len(aromatic_systems)}")
    for i, s in enumerate(aromatic_systems, start=1):
        print(f"  System {i}: atoms {list(s)}")
    print()
    base.print_summary_and_choose_ilp(atom_syms, bonds, lp_out, charge)

    coords = [[a[2], a[3], a[4]] for a in atoms]
    bo0 = {(i - 1, j - 1): o for i, j, o in bonds}
    lp_full = {i - 1: lp_out.get(i, 0) for i in range(1, len(atom_syms) + 1)}
    fc0 = [fc_out.get(a[0], 0) for a in atoms]
    idx_to_pos = {a[0]: k for k, a in enumerate(atoms)}
    metal_adj_0 = [
        (idx_to_pos[tm], idx_to_pos[lig], ei, ej)
        for tm, lig, ei, ej in raw
        if base.is_TM(ei) ^ base.is_TM(ej)
    ]
    base.print_octet_report(
        atom_syms, bo0, lp_full, fc0,
        metal_adjacency_edges=metal_adj_0, coords=coords,
    )

    lp_by_arr = {k: lp_out.get(a[0], 0) for k, a in enumerate(atoms)}
    base.print_cbc_report(
        atom_syms, coords, bo0, lp_by_arr, fc0, charge, metal_adjacency_edges=metal_adj_0
    )
    print_tm_oxidation_sigma_report(atom_syms, fc0, charge, bo=bo0)


if __name__ == "__main__":
    main()

