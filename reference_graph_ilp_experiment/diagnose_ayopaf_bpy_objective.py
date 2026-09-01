#!/usr/bin/env python3
"""Compare AYOPAF ILP objective: free solve vs bpy fixed to reference Kekulé."""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

import pulp
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from batch_tmqmg_ilp_benchmark import load_lewis_engine
from run_reference_graph_ilp import decode_edges, pairs_to_engine_edges

CSD_ID = "AYOPAF"
REF_BPY_SMILES = "c1ccc(-c2ccccn2)nc1"
XYZ_PATH = Path("/data/jingyuan_data/tmqmg") / f"{CSD_ID}.xyz"
EDGE_CSV = ROOT / "reference_graph_ilp_experiment" / "reference_edges.csv"
ENGINE_PATH = ROOT / "Lewis-engine-ILP.py"


def load_ayopaf():
    engine = load_lewis_engine(ENGINE_PATH)
    with EDGE_CSV.open(newline="", encoding="utf-8") as handle:
        row = next(r for r in csv.DictReader(handle) if r["IDs"] == CSD_ID)
    charge = int(row["charge"])
    atoms = engine.read_xyz(str(XYZ_PATH))
    edges = pairs_to_engine_edges(atoms, decode_edges(row["edges_0"]))
    aromatic = engine.aromatic_candidate_systems(atoms, edges)
    return engine, atoms, edges, aromatic, charge


def xyz_mol(atoms, edges):
    rw = Chem.RWMol()
    id_to_idx = {}
    for atom_id, el, *_ in atoms:
        idx = rw.AddAtom(Chem.Atom(el))
        id_to_idx[atom_id] = idx
    idx_to_id = {idx: atom_id for atom_id, idx in id_to_idx.items()}
    seen = set()
    for i, j, _ei, _ej in edges:
        key = (min(i, j), max(i, j))
        if key in seen:
            continue
        seen.add(key)
        rw.AddBond(id_to_idx[i], id_to_idx[j], Chem.BondType.SINGLE)
    mol = rw.GetMol()
    mol.UpdatePropertyCache(strict=False)
    return mol, id_to_idx, idx_to_id


def _element_graph(atoms, edges, keep):
    el = {atom_id: sym for atom_id, sym, *_ in atoms}
    adj = defaultdict(set)
    for i, j, _ei, _ej in edges:
        if i in keep and j in keep:
            adj[i].add(j)
            adj[j].add(i)
    return el, adj


def _isomorphism(pat_ids, pat_el, pat_adj, tgt_ids, tgt_el, tgt_adj):
    pat = list(pat_ids)
    tgt = list(tgt_ids)

    def ok(mapping, p, t):
        if pat_el[p] != tgt_el[t]:
            return False
        if len(pat_adj[p]) != len(tgt_adj[t]):
            return False
        return all(
            mapping[nbr] in tgt_adj[t]
            for nbr in pat_adj[p]
            if nbr in mapping
        )

    def rec(i, mapping, used):
        if i == len(pat):
            return dict(mapping)
        p = pat[i]
        for t in tgt:
            if t in used or not ok(mapping, p, t):
                continue
            mapping[p] = t
            used.add(t)
            found = rec(i + 1, mapping, used)
            if found:
                return found
            used.remove(t)
            del mapping[p]
        return None

    return rec(0, {}, set())


def reference_bpy_assignment(atoms, edges):
    pattern = Chem.MolFromSmiles(REF_BPY_SMILES)
    Chem.Kekulize(pattern, clearAromaticFlags=True)
    pat_ids = [atom.GetIdx() for atom in pattern.GetAtoms()]
    pat_el = {atom.GetIdx(): atom.GetSymbol() for atom in pattern.GetAtoms()}
    pat_adj = defaultdict(set)
    pat_bonds = []
    for bond in pattern.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        pat_adj[a].add(b)
        pat_adj[b].add(a)
        pat_bonds.append((a, b, int(bond.GetBondTypeAsDouble())))

    el = {atom_id: sym for atom_id, sym, *_ in atoms}
    heavy = {atom_id for atom_id, sym, *_ in atoms if sym != "H" and not engine_is_tm(sym)}
    tgt_el, tgt_adj = _element_graph(atoms, edges, heavy)
    # Restrict to the unique C10N2 component (bpy), not the aminopyridine.
    components = []
    seen = set()
    for start in heavy:
        if start in seen:
            continue
        stack = [start]
        comp = set()
        while stack:
            cur = stack.pop()
            if cur in seen:
                continue
            seen.add(cur)
            comp.add(cur)
            stack.extend(tgt_adj[cur] - seen)
        components.append(comp)
    bpy_comp = None
    for comp in components:
        formula = "".join(sorted(el[i] for i in comp))
        if formula.count("C") == 10 and formula.count("N") == 2 and len(comp) == 12:
            bpy_comp = comp
            break
    if bpy_comp is None:
        raise RuntimeError(f"No C10N2 bpy component; got {[len(c) for c in components]}")

    mapping = _isomorphism(pat_ids, pat_el, pat_adj, bpy_comp, tgt_el, tgt_adj)
    if mapping is None:
        raise RuntimeError("Could not map reference bpy onto AYOPAF connectivity")
    bond_orders = {
        (min(mapping[a], mapping[b]), max(mapping[a], mapping[b])): order
        for a, b, order in pat_bonds
    }
    bpy_ids = set(mapping.values())
    n_ids = [mapping[i] for i, sym in pat_el.items() if sym == "N"]
    tm_ids = [atom_id for atom_id, sym, *_ in atoms if engine_is_tm(sym)]
    return bond_orders, bpy_ids, n_ids, tm_ids


def engine_is_tm(el: str) -> bool:
    return el in {
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    }


def apply_bpy_reference(prob, ctx, bond_orders, n_ids, tm_ids, *, fix_mn_dative: bool):
    u2, u3, b_tm = ctx["u2"], ctx["u3"], ctx["b_tm"]
    for key, order in bond_orders.items():
        if key not in u2:
            continue
        if order == 1:
            prob += u2[key] == 0
            prob += u3[key] == 0
        elif order == 2:
            prob += u2[key] == 1
            prob += u3[key] == 0
        elif order == 3:
            prob += u2[key] == 0
            prob += u3[key] == 1
        else:
            raise RuntimeError(f"Unexpected bpy bond order {order} on {key}")
    if fix_mn_dative:
        for tm in tm_ids:
            for n_id in n_ids:
                if (tm, n_id) in b_tm:
                    prob += b_tm[(tm, n_id)] == 0


def objective_breakdown(prob, ctx, engine):
    q = {i: int(round(pulp.value(v) or 0)) for i, v in ctx["q"].items()}
    lp = {i: int(round(pulp.value(v) or 0)) for i, v in ctx["lp"].items()}
    abs_q = sum(abs(v) for v in q.values())
    atom_el = ctx["atom_el"]
    eneg = 0.0
    for i, val in q.items():
        if val < 0:
            eneg += max(0.0, 4.0 - engine.ENEG.get(atom_el[i], 2.0)) * (-val)
    arom = 0
    for var in prob.variables():
        name = var.name
        if name.endswith("_devp") or name.endswith("_devm"):
            arom += int(round(pulp.value(var) or 0))
    remote = 0
    for i, var in ctx["remote_c_lp_violation"].items():
        remote += int(round(pulp.value(var) or 0))
    ox = {
        tm: int(round(pulp.value(v) or 0)) if not isinstance(v, int) else v
        for tm, v in ctx["tm_ox_vars"].items()
    }
    btm = {
        key: int(round(pulp.value(var) or 0))
        for key, var in ctx["b_tm"].items()
    }
    return {
        "obj": float(pulp.value(prob.objective) or 0.0),
        "abs_q": abs_q,
        "fc_term": engine.ILP_WEIGHT_FORMAL_CHARGE * abs_q,
        "eneg": eneg,
        "eneg_term": engine.ILP_WEIGHT_ENEG_NEGATIVE_FC * eneg,
        "arom_dev": arom,
        "arom_term": engine.ILP_WEIGHT_AROMATIC_DEVIATION * arom,
        "remote_c_lp": remote,
        "remote_term": engine.ILP_WEIGHT_REMOTE_C_LP_VIOLATION * remote,
        "q": q,
        "lp": lp,
        "ox": ox,
        "b_tm": btm,
    }


def decode_bonds(ctx):
    bonds = []
    for i, j, _ei, _ej in ctx["ilp_edges"]:
        key = (min(i, j), max(i, j))
        order = 1 + int(round(pulp.value(ctx["u2"][key]) or 0)) + int(
            round(pulp.value(ctx["u3"][key]) or 0)
        )
        bonds.append((key[0], key[1], order))
    for (tm, lig), var in ctx["b_tm"].items():
        b_val = int(round(pulp.value(var) or 0))
        if b_val > 0:
            a, b = (tm, lig) if tm < lig else (lig, tm)
            bonds.append((a, b, b_val))
    bonds.sort()
    return bonds


def smiles_from_solution(engine, atoms, edges, ctx, charge, aromatic):
    bonds = decode_bonds(ctx)
    lp_out = {i: int(round(pulp.value(v) or 0)) for i, v in ctx["lp"].items()}
    fc_out = {i: int(round(pulp.value(v) or 0)) for i, v in ctx["q"].items()}
    for tm_i, ox_v in ctx["tm_ox_vars"].items():
        fc_out[tm_i] = int(ox_v) if isinstance(ox_v, int) else int(round(pulp.value(ox_v) or 0))
    bonds, lp_out, fc_out, _labels = engine.apply_heterocyclic_carbene_corrections(
        atoms, bonds, lp_out, fc_out, aromatic, edges, mol_charge=charge
    )
    bonds, lp_out, fc_out = engine.apply_eta_covalent_pi_corrections(
        atoms, bonds, lp_out, fc_out, mol_charge=charge, metal_adjacency_edges=edges
    )
    coords = [[atom[2], atom[3], atom[4]] for atom in atoms]
    atom_symbols = [atom[1] for atom in atoms]
    formal_charges = [fc_out.get(atom[0], 0) for atom in atoms]
    index_to_position = {atom[0]: idx for idx, atom in enumerate(atoms)}
    metal_edges_0 = [
        (index_to_position[begin], index_to_position[end], symbol_i, symbol_j)
        for begin, end, symbol_i, symbol_j in edges
        if engine.is_TM(symbol_i) ^ engine.is_TM(symbol_j)
    ]
    dative_pairs = engine.infer_dative_ml_pairs_cbc(
        atom_symbols, coords, bonds, lp_out, formal_charges, metal_adjacency_edges=metal_edges_0
    )
    return engine.ilp_to_smiles(atoms, bonds, fc_out, dative_ml_pairs=dative_pairs, edges=edges)


def solve_case(engine, atoms, edges, aromatic, charge, label, constrain=None):
    prob, ctx = engine.build_bond_order_ilp(
        atoms, edges, aromatic, mol_charge=charge, metal_adjacency_edges=edges
    )
    if constrain is not None:
        constrain(prob, ctx)
    status = prob.solve(pulp.PULP_CBC_CMD(msg=False, timeLimit=10))
    st = pulp.LpStatus[status]
    print(f"\n=== {label} ===")
    print(f"status: {st}")
    if st not in ("Optimal", "Integer Feasible"):
        return
    br = objective_breakdown(prob, ctx, engine)
    print(
        f"objective={br['obj']:.1f}  "
        f"FC={br['fc_term']:.1f}(|q|={br['abs_q']})  "
        f"eneg={br['eneg_term']:.1f}  "
        f"arom={br['arom_term']:.1f}(dev={br['arom_dev']})  "
        f"remoteC={br['remote_term']:.1f}"
    )
    print(f"ox={br['ox']}")
    nonzero_q = {i: v for i, v in sorted(br["q"].items()) if v}
    print(f"nonzero q: {nonzero_q}")
    smiles = smiles_from_solution(engine, atoms, edges, ctx, charge, aromatic)
    print(f"smiles: {smiles}")
    return br


def fmt_system(atoms, sys_atoms):
    el = {atom_id: sym for atom_id, sym, *_ in atoms}
    return " ".join(f"{el[i]}{i}" for i in sys_atoms)


def main():
    engine, atoms, edges, aromatic, charge = load_ayopaf()
    print(f"{CSD_ID} charge={charge} n_atoms={len(atoms)}")
    print(f"aromatic systems: {len(aromatic)}")
    for idx, sys_atoms in enumerate(aromatic, 1):
        print(f"  {idx}: n={len(sys_atoms)} {fmt_system(atoms, sys_atoms)}")

    bond_orders, bpy_ids, n_ids, tm_ids = reference_bpy_assignment(atoms, edges)
    el = {atom_id: sym for atom_id, sym, *_ in atoms}
    print(
        "mapped bpy "
        + " ".join(f"{el[i]}{i}" for i in sorted(bpy_ids))
        + f"  N={n_ids}  doubles="
        + str(sorted(k for k, o in bond_orders.items() if o >= 2))
    )

    free = solve_case(engine, atoms, edges, aromatic, charge, "free ILP")

    def fix_bonds_only(prob, ctx):
        apply_bpy_reference(prob, ctx, bond_orders, n_ids, tm_ids, fix_mn_dative=False)

    def fix_bonds_and_dative(prob, ctx):
        apply_bpy_reference(prob, ctx, bond_orders, n_ids, tm_ids, fix_mn_dative=True)

    kekule = solve_case(
        engine, atoms, edges, aromatic, charge,
        "bpy bonds = reference Kekulé",
        constrain=fix_bonds_only,
    )
    ref_like = solve_case(
        engine, atoms, edges, aromatic, charge,
        "bpy bonds = reference Kekulé + M–N dative",
        constrain=fix_bonds_and_dative,
    )

    print("\n=== delta vs free ===")
    for name, other in (("Kekulé bonds only", kekule), ("Kekulé + M–N dative", ref_like)):
        if free and other:
            dobj = other["obj"] - free["obj"]
            print(
                f"{name}: Δobj={dobj:+.1f}  "
                f"ΔFC={other['fc_term']-free['fc_term']:+.1f}  "
                f"Δarom={other['arom_term']-free['arom_term']:+.1f}  "
                f"Δeneg={other['eneg_term']-free['eneg_term']:+.1f}"
            )


if __name__ == "__main__":
    main()
