#!/usr/bin/env python3
"""Split xyz2mol fails into Hückel-charge errors vs downstream AC2mol fails."""

from __future__ import annotations

import csv
import logging
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, TimeoutError as FutTimeout
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")
logging.disable(logging.CRITICAL)

sys.path.insert(0, "/home/zhujingyuan/xyz2mol_tm")
from xyz2mol_tm.huckel_to_smiles.xyz2mol_local import AC2mol  # noqa: E402
from xyz2mol_tm.huckel_to_smiles.xyz2mol_tmc import (  # noqa: E402
    MetalNon_Hg,
    TRANSITION_METALS_NUM,
    get_basic_mol,
    get_proposed_ligand_charge,
    params,
)

ROOT = Path("/home/zhujingyuan/TMC_Lewis_Structure")
XYZ_DIR = Path("/data/jingyuan_data/tmqmg")
TM = set(TRANSITION_METALS_NUM)
_BOND_ORDER = {
    Chem.BondType.SINGLE: 1,
    Chem.BondType.DOUBLE: 2,
    Chem.BondType.TRIPLE: 3,
}


def heavy_formula(mol) -> tuple:
    c = Counter()
    for atom in mol.GetAtoms():
        z = atom.GetAtomicNum()
        if z > 1:
            c[z] += 1
    return tuple(sorted(c.items()))


def formula_str(key: tuple) -> str:
    pt = Chem.GetPeriodicTable()
    return "".join(f"{pt.GetElementSymbol(z)}{n if n > 1 else ''}" for z, n in key)


def ilp_ligand_charges(smiles: str):
    """Ionic ligand charges after removing the metal from an ILP SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = Chem.RemoveHs(mol)
    except Exception:
        pass
    tms = [a for a in mol.GetAtoms() if a.GetAtomicNum() in TM]
    if len(tms) != 1:
        return None
    tm = tms[0]
    ox = int(tm.GetFormalCharge())
    x_delta = {}
    for nbr in tm.GetNeighbors():
        bond = mol.GetBondBetweenAtoms(tm.GetIdx(), nbr.GetIdx())
        if bond is None or bond.GetBondType() == Chem.BondType.DATIVE:
            continue
        order = _BOND_ORDER.get(
            bond.GetBondType(), int(round(bond.GetBondTypeAsDouble() or 1))
        )
        if order < 1:
            order = 1
        x_delta[nbr.GetIdx()] = x_delta.get(nbr.GetIdx(), 0) - order
    rw = Chem.RWMol(mol)
    tm_idx = tm.GetIdx()
    for nbr in list(rw.GetAtomWithIdx(tm_idx).GetNeighbors()):
        rw.RemoveBond(tm_idx, nbr.GetIdx())
    rw.RemoveAtom(tm_idx)
    leftover = rw.GetMol()
    frags = Chem.GetMolFrags(leftover, asMols=True, sanitizeFrags=False)
    atom_frags = Chem.GetMolFrags(leftover, asMols=False)
    out = []
    for frag, idxs in zip(frags, atom_frags):
        q = 0
        for leftover_idx in idxs:
            orig = leftover_idx if leftover_idx < tm_idx else leftover_idx + 1
            q += int(mol.GetAtomWithIdx(orig).GetFormalCharge())
            q += x_delta.get(orig, 0)
        out.append((heavy_formula(frag), q))
    return {"ox": ox, "ligands": out, "q_lig_sum": sum(q for _, q in out)}


def xyz_huckel_frags(xyz_path: str, charge: int):
    mol = get_basic_mol(xyz_path, charge)
    tmc_idx = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() in TM)
    mdis = rdMolStandardize.MetalDisconnector(params)
    mdis.SetMetalNon(Chem.MolFromSmarts(MetalNon_Hg))
    frags = mdis.Disconnect(mol)
    frag_mols = Chem.GetMolFrags(frags, asMols=True)
    ligands = []
    for f in frag_mols:
        if any(a.GetAtomicNum() in TM for a in f.GetAtoms()):
            continue
        qh = int(get_proposed_ligand_charge(f))
        atoms = [a.GetAtomicNum() for a in f.GetAtoms()]
        ac = Chem.rdmolops.GetAdjacencyMatrix(f)
        ac_ok = AC2mol(f, ac, atoms, qh, allow_charged_fragments=True, use_atom_maps=False)
        ligands.append(
            {
                "formula": heavy_formula(f),
                "n_atoms": f.GetNumAtoms(),
                "q_huckel": qh,
                "ac2mol_huckel": ac_ok is not None,
            }
        )
    return ligands


def match_and_test(xyz_ligs, ilp_ligs):
    unused = list(ilp_ligs)
    rows = []
    for lig in xyz_ligs:
        pick = None
        for i, (form, q) in enumerate(unused):
            if form == lig["formula"]:
                pick = i
                break
        q_ilp = None
        matched = False
        if pick is not None:
            _, q_ilp = unused.pop(pick)
            matched = True
        ac_ok_ilp = None
        rows.append({**lig, "q_ilp": q_ilp, "matched": matched, "ac2mol_ilp": ac_ok_ilp})
    return rows, unused


def ac2mol_with_charge(xyz_path: str, charge: int, formula, q):
    mol = get_basic_mol(xyz_path, charge)
    mdis = rdMolStandardize.MetalDisconnector(params)
    mdis.SetMetalNon(Chem.MolFromSmarts(MetalNon_Hg))
    frag_mols = Chem.GetMolFrags(mdis.Disconnect(mol), asMols=True)
    for f in frag_mols:
        if any(a.GetAtomicNum() in TM for a in f.GetAtoms()):
            continue
        if heavy_formula(f) != formula:
            continue
        atoms = [a.GetAtomicNum() for a in f.GetAtoms()]
        ac = Chem.rdmolops.GetAdjacencyMatrix(f)
        return AC2mol(f, ac, atoms, q, allow_charged_fragments=True, use_atom_maps=False) is not None
    return None


def analyze_one(row):
    cid = row["IDs"]
    xyz = XYZ_DIR / f"{cid}.xyz"
    if not xyz.exists():
        return {"IDs": cid, "status": "missing_xyz"}
    try:
        ilp = ilp_ligand_charges(row["smiles_ilp"])
        xyz_ligs = xyz_huckel_frags(str(xyz), int(row["charge"]))
    except Exception as exc:
        return {"IDs": cid, "status": f"error:{type(exc).__name__}"}
    if ilp is None:
        return {"IDs": cid, "status": "ilp_parse_fail"}
    q_complex = int(row["charge"])
    q_huckel_sum = sum(x["q_huckel"] for x in xyz_ligs)
    q_ilp_sum = ilp["q_lig_sum"]
    unused = list(ilp["ligands"])
    ligand_rows = []
    n_wrong = n_right = n_unmatched = 0
    n_huckel_cause = n_both_fail = n_down = n_huckel_ok_ac_ok = 0
    details = []
    for lig in xyz_ligs:
        pick = None
        for i, (form, q) in enumerate(unused):
            if form == lig["formula"]:
                pick = i
                break
        if pick is None:
            n_unmatched += 1
            details.append(f"{formula_str(lig['formula'])}:H={lig['q_huckel']}/ILP=? acH={int(lig['ac2mol_huckel'])}")
            continue
        _, q_ilp = unused.pop(pick)
        same = lig["q_huckel"] == q_ilp
        if same:
            n_right += 1
        else:
            n_wrong += 1
        ac_ilp = None
        if not same:
            try:
                ac_ilp = ac2mol_with_charge(str(xyz), q_complex, lig["formula"], q_ilp)
            except Exception:
                ac_ilp = False
        if same and not lig["ac2mol_huckel"]:
            n_down += 1
            tag = "huckel_ok_ac_fail"
        elif same and lig["ac2mol_huckel"]:
            n_huckel_ok_ac_ok += 1
            tag = "huckel_ok_ac_ok"
        elif (not same) and (not lig["ac2mol_huckel"]) and ac_ilp:
            n_huckel_cause += 1
            tag = "huckel_wrong_ac_ok_if_correct"
        elif (not same) and (not lig["ac2mol_huckel"]) and not ac_ilp:
            n_both_fail += 1
            tag = "huckel_wrong_ac_fail_anyway"
        elif (not same) and lig["ac2mol_huckel"]:
            tag = "huckel_wrong_ac_ok_anyway"
        else:
            tag = "other"
        details.append(
            f"{formula_str(lig['formula'])}:H={lig['q_huckel']}/ILP={q_ilp} "
            f"acH={int(lig['ac2mol_huckel'])} acILP={ac_ilp} {tag}"
        )
    struct_huckel_ok = q_huckel_sum == q_ilp_sum
    return {
        "IDs": cid,
        "status": "ok",
        "error_type": row.get("error_type", ""),
        "n_atoms": row.get("n_atoms", ""),
        "q_complex": q_complex,
        "ox_ilp": ilp["ox"],
        "q_huckel_sum": q_huckel_sum,
        "q_ilp_sum": q_ilp_sum,
        "sum_match": int(struct_huckel_ok),
        "n_lig": len(xyz_ligs),
        "n_charge_wrong": n_wrong,
        "n_charge_right": n_right,
        "n_unmatched": n_unmatched,
        "n_huckel_cause": n_huckel_cause,
        "n_both_fail": n_both_fail,
        "n_downstream": n_down,
        "n_both_ok": n_huckel_ok_ac_ok,
        "detail": " | ".join(details),
    }


def load_jobs():
    fail = {}
    with (ROOT / "ilp_ok_xyz2mol_fail.csv").open() as fh:
        for r in csv.DictReader(fh):
            if r["xyz2mol_status"] == "timeout":
                continue
            fail[r["IDs"]] = r
    with (ROOT / "tmqmg_ilp_benchmark_output" / "per_structure_results.csv").open() as fh:
        for r in csv.DictReader(fh):
            if r["IDs"] in fail:
                fail[r["IDs"]]["smiles_ilp"] = r["smiles_ilp"]
                fail[r["IDs"]]["charge"] = r["charge"]
    return [fail[k] for k in sorted(fail) if fail[k].get("smiles_ilp")]


def main():
    jobs = load_jobs()
    if len(sys.argv) > 1:
        want = set(sys.argv[1:])
        jobs = [j for j in jobs if j["IDs"] in want]
    print(f"jobs {len(jobs)}", flush=True)
    out_path = ROOT / "huckel_vs_downstream.csv"
    rows = []
    with ProcessPoolExecutor(max_workers=8) as ex:
        futs = [ex.submit(analyze_one, j) for j in jobs]
        for i, fut in enumerate(futs, 1):
            try:
                rows.append(fut.result(timeout=45))
            except FutTimeout:
                rows.append({"IDs": jobs[i - 1]["IDs"], "status": "timeout45"})
            except Exception as exc:
                rows.append({"IDs": jobs[i - 1]["IDs"], "status": f"error:{type(exc).__name__}"})
            if i % 50 == 0 or i == len(futs):
                print(f"  {i}/{len(futs)}", flush=True)
    keys = [
        "IDs", "status", "error_type", "n_atoms", "q_complex", "ox_ilp",
        "q_huckel_sum", "q_ilp_sum", "sum_match", "n_lig",
        "n_charge_wrong", "n_charge_right", "n_unmatched",
        "n_huckel_cause", "n_both_fail", "n_downstream", "n_both_ok", "detail",
    ]
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    done = [r for r in rows if r.get("status") == "ok"]
    print(f"wrote {out_path}  parsed {len(done)}/{len(rows)}")
    print("struct Hückel sum == ILP ligand sum:", sum(int(r["sum_match"]) for r in done), "/", len(done))
    print("ligand Hückel != ILP:", sum(int(r["n_charge_wrong"]) for r in done))
    print("ligand Hückel == ILP:", sum(int(r["n_charge_right"]) for r in done))
    print("Hückel wrong, AC2mol OK if correct charge:", sum(int(r["n_huckel_cause"]) for r in done))
    print("Hückel wrong, AC2mol fails anyway:", sum(int(r["n_both_fail"]) for r in done))
    print("Hückel right, AC2mol fails:", sum(int(r["n_downstream"]) for r in done))
    print("Hückel right, AC2mol OK (fail later):", sum(int(r["n_both_ok"]) for r in done))


if __name__ == "__main__":
    main()
