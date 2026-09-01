#!/usr/bin/env python3
"""Stage-trace xyz2mol failures where Hückel and AC2mol already succeeded."""

from __future__ import annotations

import csv
import logging
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, TimeoutError as FutTimeout
from itertools import combinations
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdmolops
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")
logging.disable(logging.CRITICAL)

sys.path.insert(0, "/home/zhujingyuan/xyz2mol_tm")
from xyz2mol_tm.huckel_to_smiles.xyz2mol_tmc import (  # noqa: E402
    MetalNon_Hg,
    TRANSITION_METALS_NUM,
    fix_equivalent_Os,
    fix_NO2,
    get_basic_mol,
    get_lig_mol,
    get_proposed_ligand_charge,
    params,
)

ROOT = Path("/home/zhujingyuan/TMC_Lewis_Structure")
XYZ_DIR = Path("/data/jingyuan_data/tmqmg")
TM = set(TRANSITION_METALS_NUM)


def load_e_ids():
    rows = list(csv.DictReader((ROOT / "huckel_vs_downstream.csv").open()))
    ids = []
    for r in rows:
        if r["status"] != "ok":
            continue
        if (
            int(r["n_charge_wrong"] or 0) == 0
            and int(r["n_unmatched"] or 0) == 0
            and int(r["n_downstream"] or 0) == 0
            and int(r["n_both_fail"] or 0) == 0
        ):
            ids.append(r["IDs"])
    charge = {}
    with (ROOT / "tmqmg_ilp_benchmark_output" / "per_structure_results.csv").open() as fh:
        for r in csv.DictReader(fh):
            if r["IDs"] in set(ids):
                charge[r["IDs"]] = int(r["charge"])
    return [(cid, charge[cid]) for cid in ids if cid in charge]


def trace_one(item):
    cid, overall_charge = item
    xyz = XYZ_DIR / f"{cid}.xyz"
    try:
        mol = get_basic_mol(str(xyz), overall_charge)
        tmc_idx = None
        for a in mol.GetAtoms():
            a.SetIntProp("__origIdx", a.GetIdx())
            if a.GetAtomicNum() in TM:
                tmc_idx = a.GetIdx()
        coordinating_atoms = set(
            int(i) for i in rdmolops.GetAdjacencyMatrix(mol)[tmc_idx, :].nonzero()[0]
        )
        mdis = rdMolStandardize.MetalDisconnector(params)
        mdis.SetMetalNon(Chem.MolFromSmarts(MetalNon_Hg))
        frag_mols = rdmolops.GetMolFrags(mdis.Disconnect(mol), asMols=True)
        lig_list = []
        total_lig_charge = 0
        tm_idx = None
        for i, f in enumerate(frag_mols):
            if any(a.GetAtomicNum() in TM for a in f.GetAtoms()):
                tm_idx = i
                continue
            qh = get_proposed_ligand_charge(f)
            coord = [
                a.GetIdx()
                for a in f.GetAtoms()
                if a.GetIntProp("__origIdx") in coordinating_atoms
            ]
            try:
                lig_mol, lig_charge = get_lig_mol(f, qh, coord)
            except Exception as exc:
                return cid, "get_lig_mol", type(exc).__name__, str(exc).split("\n")[0][:160]
            if lig_mol is None:
                return cid, "get_lig_mol_none", "None", f"charge={qh}"
            total_lig_charge += lig_charge
            lig_list.append(lig_mol)
        tm = Chem.RWMol(frag_mols[tm_idx])
        tm_ox = overall_charge - total_lig_charge
        for a in tm.GetAtoms():
            if a.GetAtomicNum() in TM:
                a.SetFormalCharge(tm_ox)
        for lmol in lig_list:
            tm = Chem.CombineMols(tm, lmol)
        emol = Chem.RWMol(tm)
        coordinating_atoms_idx = [
            a.GetIdx()
            for a in emol.GetAtoms()
            if a.HasProp("__origIdx") and a.GetIntProp("__origIdx") in coordinating_atoms
        ]
        tm_new = [
            a.GetIdx()
            for a in emol.GetAtoms()
            if a.HasProp("__origIdx") and a.GetIntProp("__origIdx") == tmc_idx
        ][0]
        dMat = Chem.Get3DDistanceMatrix(emol)
        cut_atoms = []
        for i, j in combinations(coordinating_atoms_idx, 2):
            bond = emol.GetBondBetweenAtoms(int(i), int(j))
            if bond and abs(dMat[i, tm_new] - dMat[j, tm_new]) >= 0.4:
                if dMat[i, tm_new] > dMat[j, tm_new] and i in coordinating_atoms_idx:
                    coordinating_atoms_idx.remove(i)
                    cut_atoms.append(i)
                if dMat[j, tm_new] > dMat[i, tm_new] and j in coordinating_atoms_idx:
                    coordinating_atoms_idx.remove(j)
                    cut_atoms.append(j)
        for j in cut_atoms:
            for i in list(coordinating_atoms_idx):
                bond = emol.GetBondBetweenAtoms(int(i), int(j))
                if bond and dMat[i, tm_new] - dMat[j, tm_new] >= -0.1:
                    coordinating_atoms_idx.remove(i)
        for i in coordinating_atoms_idx:
            if emol.GetBondBetweenAtoms(i, tm_new):
                continue
            emol.AddBond(i, tm_new, Chem.BondType.DATIVE)
        try:
            smiles = Chem.MolToSmiles(emol.GetMol())
        except Exception as exc:
            return cid, "mol_to_smiles", type(exc).__name__, str(exc).split("\n")[0][:160]
        try:
            smiles = fix_equivalent_Os(smiles)
        except Exception as exc:
            return cid, "fix_equivalent_Os", type(exc).__name__, str(exc).split("\n")[0][:160]
        try:
            smiles = fix_NO2(smiles)
        except Exception as exc:
            return cid, "fix_NO2", type(exc).__name__, str(exc).split("\n")[0][:160]
        parsed = Chem.MolFromSmiles(smiles)
        if parsed is None:
            return cid, "reparse_smiles", "None", smiles[:120]
        try:
            Chem.SanitizeMol(parsed)
        except Exception as exc:
            return cid, "final_sanitize", type(exc).__name__, str(exc).split("\n")[0][:160]
        return cid, "unexpected_success", "", ""
    except Exception as exc:
        return cid, "other", type(exc).__name__, str(exc).split("\n")[0][:160]


def main():
    jobs = load_e_ids()
    if len(sys.argv) > 1:
        want = set(sys.argv[1:])
        jobs = [j for j in jobs if j[0] in want]
    print(f"jobs {len(jobs)}", flush=True)
    out = ROOT / "late_fail_stages.csv"
    rows = []
    with ProcessPoolExecutor(max_workers=8) as ex:
        futs = [ex.submit(trace_one, j) for j in jobs]
        for i, fut in enumerate(futs, 1):
            try:
                rows.append(fut.result(timeout=40))
            except FutTimeout:
                rows.append((jobs[i - 1][0], "timeout40", "", ""))
            if i % 40 == 0 or i == len(futs):
                print(f"  {i}/{len(futs)}", flush=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["IDs", "stage", "exc", "msg"])
        w.writerows(rows)
    print("stage", Counter(r[1] for r in rows))
    print("exc  ", Counter(r[2] for r in rows if r[1] != "unexpected_success"))
    print("wrote", out)


if __name__ == "__main__":
    main()
