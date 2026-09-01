#!/usr/bin/env python3
"""Pairwise full-connectivity comparison for four TMC SMILES sources.

The graph key preserves atom elements and heavy-atom adjacency while ignoring
bond order, aromaticity, formal charge, radicals, stereochemistry, and bond
direction/type (including covalent versus dative metal-ligand bonds).
"""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent
DEFAULT_ILP = ROOT / "ilp_mismatch_rerun_current" / "merged_per_structure_results.csv"
DEFAULT_REF = ROOT / "tmqmg_smiles.csv"
DEFAULT_OUT = ROOT / "connectivity_graph_comparison"

REF_COLUMNS = {
    "CSD": "smiles_CSD_fix",
    "NBO(DFT)": "smiles_NBO_DFT_xyz",
    "Hückel(DFT)": "smiles_huckel_DFT_xyz",
}
METHODS = ("ILP", "CSD", "NBO(DFT)", "Hückel(DFT)")
PAIRS = [
    (METHODS[i], METHODS[j])
    for i in range(len(METHODS))
    for j in range(i + 1, len(METHODS))
]
FAIL_TOKENS = {"", "fail", "API_smiles_missing", "not_in_database", "nan"}


def valid_cell(value: str | None) -> bool:
    return value is not None and str(value).strip().lower() not in FAIL_TOKENS


def parse_mol(smiles: str) -> Chem.Mol | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return mol
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return mol


def connectivity_key(smiles: str) -> tuple[tuple[int, ...], tuple[tuple[int, int], ...]] | None:
    """Return a canonical element-labelled heavy-atom adjacency certificate."""
    mol = parse_mol(smiles)
    if mol is None:
        return None
    try:
        mol = Chem.RemoveHs(mol)
        rw = Chem.RWMol(mol)
        for atom in rw.GetAtoms():
            atom.SetFormalCharge(0)
            atom.SetNumRadicalElectrons(0)
            atom.SetIsAromatic(False)
            atom.SetIsotope(0)
            atom.SetAtomMapNum(0)
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
            atom.SetNumExplicitHs(0)
            atom.SetNoImplicit(True)
        for bond in rw.GetBonds():
            bond.SetBondType(Chem.BondType.SINGLE)
            bond.SetIsAromatic(False)
            bond.SetBondDir(Chem.BondDir.NONE)
            bond.SetStereo(Chem.BondStereo.STEREONONE)
        graph = rw.GetMol()
        graph.UpdatePropertyCache(strict=False)
        ranks = list(
            Chem.CanonicalRankAtoms(
                graph,
                breakTies=True,
                includeChirality=False,
                includeIsotopes=False,
            )
        )
    except Exception:
        return None

    order = sorted(range(graph.GetNumAtoms()), key=lambda idx: ranks[idx])
    canonical_index = {old: new for new, old in enumerate(order)}
    atoms = tuple(graph.GetAtomWithIdx(idx).GetAtomicNum() for idx in order)
    edges = tuple(
        sorted(
            (
                min(canonical_index[bond.GetBeginAtomIdx()], canonical_index[bond.GetEndAtomIdx()]),
                max(canonical_index[bond.GetBeginAtomIdx()], canonical_index[bond.GetEndAtomIdx()]),
            )
            for bond in graph.GetBonds()
        )
    )
    return atoms, edges


def load_ilp(path: Path) -> dict[str, str]:
    output = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row.get("status") != "ok":
                continue
            if str(row.get("rdkit_parse_ok", "")).strip().lower() not in {"1", "true"}:
                continue
            if valid_cell(row.get("smiles_ilp")):
                output[row["IDs"]] = row["smiles_ilp"]
    return output


def load_references(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            row["IDs"]: {method: row.get(column, "") for method, column in REF_COLUMNS.items()}
            for row in csv.DictReader(handle)
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ilp-csv", type=Path, default=DEFAULT_ILP)
    parser.add_argument("--ref-csv", type=Path, default=DEFAULT_REF)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    ilp = load_ilp(args.ilp_csv)
    references = load_references(args.ref_csv)
    ids = sorted(set(ilp) | set(references))

    smiles_by_id: dict[str, dict[str, str]] = {}
    unique_smiles: list[str] = []
    seen = set()
    for csd_id in ids:
        values = {"ILP": ilp.get(csd_id, ""), **references.get(csd_id, {})}
        smiles_by_id[csd_id] = values
        for smiles in values.values():
            if valid_cell(smiles) and smiles not in seen:
                seen.add(smiles)
                unique_smiles.append(smiles)

    if args.workers <= 1:
        keys = [connectivity_key(smiles) for smiles in unique_smiles]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            keys = list(pool.map(connectivity_key, unique_smiles, chunksize=32))
    cache = dict(zip(unique_smiles, keys))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    detail_path = args.out_dir / "per_structure.csv"
    summary_path = args.out_dir / "summary.csv"
    fields = [
        "IDs",
        "pair",
        "both_parsed",
        "connectivity_match",
        "atoms_a",
        "atoms_b",
        "edges_a",
        "edges_b",
    ]
    tallies = {f"{a} / {b}": [0, 0] for a, b in PAIRS}

    with detail_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for csd_id in ids:
            values = smiles_by_id[csd_id]
            features = {
                method: cache.get(values.get(method, ""))
                if valid_cell(values.get(method, ""))
                else None
                for method in METHODS
            }
            for method_a, method_b in PAIRS:
                label = f"{method_a} / {method_b}"
                key_a, key_b = features[method_a], features[method_b]
                both = key_a is not None and key_b is not None
                match = both and key_a == key_b
                if both:
                    tallies[label][1] += 1
                    tallies[label][0] += int(match)
                writer.writerow(
                    {
                        "IDs": csd_id,
                        "pair": label,
                        "both_parsed": int(both),
                        "connectivity_match": int(match) if both else "",
                        "atoms_a": len(key_a[0]) if key_a else "",
                        "atoms_b": len(key_b[0]) if key_b else "",
                        "edges_a": len(key_a[1]) if key_a else "",
                        "edges_b": len(key_b[1]) if key_b else "",
                    }
                )

    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["pair", "matches", "denominator", "match_percent"])
        for method_a, method_b in PAIRS:
            label = f"{method_a} / {method_b}"
            matches, denominator = tallies[label]
            percent = 100.0 * matches / denominator if denominator else 0.0
            writer.writerow([label, matches, denominator, f"{percent:.4f}"])
            print(f"{label}: {matches:,} / {denominator:,} ({percent:.2f}%)")

    print(f"Wrote {summary_path}")
    print(f"Wrote {detail_path}")


if __name__ == "__main__":
    main()
