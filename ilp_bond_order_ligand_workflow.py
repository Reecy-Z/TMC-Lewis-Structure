#!/usr/bin/env python3
"""
Ligand-wise ILP workflow based on ilp_bond_order_example.py.

Pipeline:
1) Read XYZ and build raw connectivity.
2) Split metal-neighbor graph into ligand components (cut all metal-related edges).
3) Ask user to assign charge for each ligand component.
4) Run solve_bond_orders() independently for each ligand.
5) Merge ligand Lewis results and analyze bonding around metals (CBC report).

Usage:
    python ilp_bond_order_ligand_workflow.py <molecule.xyz>
"""

from __future__ import annotations

import argparse
from collections import defaultdict

import ilp_bond_order_example as base


def _build_ligand_components(atoms_packed, raw_edges):
    """
    Return ligand components as sorted tuples of 1-based atom indices.

    Rule:
    - Metals are identified by base.is_TM(symbol).
    - Remove every edge touching any metal.
    - Keep only non-metal atoms that are neighbors of at least one metal in raw graph.
    - Connected components in this cut graph are ligands.
    """
    atom_symbol = {idx: el for idx, el, *_ in atoms_packed}
    metal_set = {idx for idx, el, *_ in atoms_packed if base.is_TM(el)}

    # Non-metal atoms directly contacting metals in raw connectivity.
    ligand_seed = set()
    adj = defaultdict(set)
    for i, j, _ei, _ej in raw_edges:
        i_m = i in metal_set
        j_m = j in metal_set
        if i_m and not j_m:
            ligand_seed.add(j)
            continue
        if j_m and not i_m:
            ligand_seed.add(i)
            continue
        if i_m or j_m:
            continue
        if base.is_TM(atom_symbol[i]) or base.is_TM(atom_symbol[j]):
            continue
        adj[i].add(j)
        adj[j].add(i)

    visited = set()
    comps = []
    for s in sorted(ligand_seed):
        if s in visited:
            continue
        stack = [s]
        visited.add(s)
        comp = []
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj.get(u, ()):
                if v not in visited:
                    visited.add(v)
                    stack.append(v)
        comps.append(tuple(sorted(comp)))
    return comps


def _subproblem_for_component(atoms_packed, raw_edges, comp):
    """Build ligand-only atoms/edges and index maps for one component."""
    local_to_global = list(comp)
    global_to_local = {g: i + 1 for i, g in enumerate(local_to_global)}

    sub_atoms = []
    for li, g in enumerate(local_to_global, start=1):
        _idx, el, x, y, z = atoms_packed[g - 1]
        sub_atoms.append((li, el, x, y, z))

    comp_set = set(comp)
    sub_edges = []
    for i, j, ei, ej in raw_edges:
        if i in comp_set and j in comp_set:
            li = global_to_local[i]
            lj = global_to_local[j]
            if li < lj:
                sub_edges.append((li, lj, ei, ej))
            else:
                sub_edges.append((lj, li, ej, ei))
    # remove duplicated edges if any
    sub_edges = sorted(set(sub_edges))
    return sub_atoms, sub_edges, local_to_global, global_to_local


def _prompt_ligand_charges(ligand_components):
    charges = []
    print()
    print("Detected ligand components after cutting all metal-related edges:")
    for idx, comp in enumerate(ligand_components, start=1):
        span = f"{min(comp)}-{max(comp)}"
        members = " ".join(str(a) for a in comp)
        print(f"  Ligand {idx}: atoms [{span}]  members: {members}")
    print()
    for idx, comp in enumerate(ligand_components, start=1):
        while True:
            raw = input(f"Assign charge for ligand {idx} (atoms {list(comp)}): ").strip()
            try:
                q = int(raw)
                charges.append(q)
                break
            except ValueError:
                print("  Please enter an integer charge (e.g. -1, 0, +1).")
    return charges


def run_ligand_workflow(xyz_path):
    atoms_packed = base.read_xyz(xyz_path)
    atom_syms = [a[1] for a in atoms_packed]
    coords = [[a[2], a[3], a[4]] for a in atoms_packed]

    raw_edges = base.connectivity(atoms_packed)
    ligand_components = _build_ligand_components(atoms_packed, raw_edges)

    if not ligand_components:
        raise RuntimeError(
            "No ligand components found after cutting metal-related edges. "
            "Check connectivity or input structure."
        )

    charges = _prompt_ligand_charges(ligand_components)

    # Global merged outputs
    bo_global = {}          # 0-based (i,j) -> order
    lp_global = {}          # 0-based idx -> lp
    fc_global = [0] * len(atom_syms)

    print()
    print("Per-ligand ILP results:")
    for li, (comp, q_lig) in enumerate(zip(ligand_components, charges), start=1):
        sub_atoms, sub_edges, local_to_global, _ = _subproblem_for_component(
            atoms_packed, raw_edges, comp
        )
        bonds_sub, lp_sub, fc_sub = base.solve_bond_orders(
            sub_atoms, sub_edges, mol_charge=q_lig
        )

        print(f"  Ligand {li}: charge={q_lig}, atoms={list(comp)}")
        print(f"    bonds: {len(bonds_sub)}, lone-pair atoms: {len([k for k,v in lp_sub.items() if v > 0])}")

        # Merge local -> global indices
        for a, b, order in bonds_sub:
            ga = local_to_global[a - 1] - 1
            gb = local_to_global[b - 1] - 1
            bo_global[(min(ga, gb), max(ga, gb))] = order

        for a, lpv in lp_sub.items():
            ga = local_to_global[a - 1] - 1
            if lpv > 0:
                lp_global[ga] = lpv

        for a, qv in fc_sub.items():
            ga = local_to_global[a - 1] - 1
            fc_global[ga] = qv

    print()
    print("Merged ligand-only Lewis result:")
    merged_choose = base.format_choose_style(
        [(i + 1, j + 1, o) for (i, j), o in sorted(bo_global.items())]
    )
    print(" ", merged_choose)
    lone_items = sorted((i + 1, v) for i, v in lp_global.items() if v > 0)
    if lone_items:
        print("  LONE", " ".join(f"{i} {v}" for i, v in lone_items), "END")
    else:
        print("  LONE (none)")

    print()
    print("Metal-neighbour bonding analysis (CBC):")
    base.print_cbc_report(atom_syms, coords, bo_global, lp_global, fc_global, charge=0)


def main():
    parser = argparse.ArgumentParser(
        description="Ligand-wise ILP bond-order assignment workflow"
    )
    parser.add_argument("xyz_path", help="Input XYZ file")
    args = parser.parse_args()
    run_ligand_workflow(args.xyz_path)


if __name__ == "__main__":
    main()

