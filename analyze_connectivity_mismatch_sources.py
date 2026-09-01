#!/usr/bin/env python3
"""Decompose full-graph mismatches into ligand-internal and metal-shell terms."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from rdkit import Chem, RDLogger

from compare_smiles_connectivity_graphs import (
    DEFAULT_ILP,
    DEFAULT_REF,
    METHODS,
    PAIRS,
    load_ilp,
    load_references,
    parse_mol,
    valid_cell,
)

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "connectivity_mismatch_sources"
TM_ATOMIC_NUMS = {
    21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
    39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
    57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,
    72, 73, 74, 75, 76, 77, 78, 79, 80,
}


def graph_key(graph: Chem.Mol) -> tuple[tuple[int, ...], tuple[tuple[int, int], ...]]:
    ranks = list(
        Chem.CanonicalRankAtoms(
            graph,
            breakTies=True,
            includeChirality=False,
            includeIsotopes=False,
        )
    )
    order = sorted(range(graph.GetNumAtoms()), key=lambda idx: ranks[idx])
    canonical_index = {old: new for new, old in enumerate(order)}
    atoms = tuple(graph.GetAtomWithIdx(idx).GetAtomicNum() for idx in order)
    edges = tuple(
        sorted(
            (
                min(canonical_index[b.GetBeginAtomIdx()], canonical_index[b.GetEndAtomIdx()]),
                max(canonical_index[b.GetBeginAtomIdx()], canonical_index[b.GetEndAtomIdx()]),
            )
            for b in graph.GetBonds()
        )
    )
    return atoms, edges


def normalized_heavy_graph(smiles: str) -> Chem.Mol | None:
    mol = parse_mol(smiles)
    if mol is None:
        return None
    try:
        rw = Chem.RWMol(Chem.RemoveHs(mol))
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
        return graph
    except Exception:
        return None


def connectivity_features(smiles: str) -> dict | None:
    graph = normalized_heavy_graph(smiles)
    if graph is None:
        return None
    metal_indices = [
        atom.GetIdx()
        for atom in graph.GetAtoms()
        if atom.GetAtomicNum() in TM_ATOMIC_NUMS
    ]
    if not metal_indices:
        return None
    shells = tuple(
        sorted(
            (
                graph.GetAtomWithIdx(idx).GetAtomicNum(),
                tuple(
                    sorted(
                        neighbor.GetAtomicNum()
                        for neighbor in graph.GetAtomWithIdx(idx).GetNeighbors()
                    )
                ),
            )
            for idx in metal_indices
        )
    )
    ligand = Chem.RWMol(graph)
    for idx in sorted(metal_indices, reverse=True):
        ligand.RemoveAtom(idx)
    ligand_graph = ligand.GetMol()
    ligand_graph.UpdatePropertyCache(strict=False)
    try:
        return {
            "full": graph_key(graph),
            "ligand": graph_key(ligand_graph),
            "metal_shell": shells,
        }
    except Exception:
        return None


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
    smiles_by_id = {
        csd_id: {"ILP": ilp.get(csd_id, ""), **references.get(csd_id, {})}
        for csd_id in ids
    }
    unique_smiles = list(
        dict.fromkeys(
            smiles
            for values in smiles_by_id.values()
            for smiles in values.values()
            if valid_cell(smiles)
        )
    )
    if args.workers <= 1:
        extracted = [connectivity_features(smiles) for smiles in unique_smiles]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            extracted = list(pool.map(connectivity_features, unique_smiles, chunksize=32))
    cache = dict(zip(unique_smiles, extracted))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    detail_path = args.out_dir / "per_structure.csv"
    summary_path = args.out_dir / "summary.csv"
    fields = [
        "IDs",
        "pair",
        "both_parsed",
        "full_match",
        "ligand_internal_match",
        "metal_shell_match",
        "mismatch_class",
    ]
    counts = {f"{a} / {b}": Counter() for a, b in PAIRS}

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
                feat_a, feat_b = features[method_a], features[method_b]
                both = feat_a is not None and feat_b is not None
                row = {"IDs": csd_id, "pair": label, "both_parsed": int(both)}
                if not both:
                    row.update(
                        full_match="",
                        ligand_internal_match="",
                        metal_shell_match="",
                        mismatch_class="",
                    )
                    writer.writerow(row)
                    continue

                full_match = feat_a["full"] == feat_b["full"]
                ligand_match = feat_a["ligand"] == feat_b["ligand"]
                metal_match = feat_a["metal_shell"] == feat_b["metal_shell"]
                if full_match:
                    mismatch_class = "full_match"
                elif ligand_match:
                    mismatch_class = "metal_only"
                elif metal_match:
                    mismatch_class = "ligand_only_or_same-element-metal_ambiguity"
                else:
                    mismatch_class = "ligand_and_metal"
                counts[label][mismatch_class] += 1
                counts[label]["denominator"] += 1
                row.update(
                    full_match=int(full_match),
                    ligand_internal_match=int(ligand_match),
                    metal_shell_match=int(metal_match),
                    mismatch_class=mismatch_class,
                )
                writer.writerow(row)

    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "pair",
            "denominator",
            "full_match",
            "metal_only",
            "ligand_only_or_same-element-metal_ambiguity",
            "ligand_and_metal",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for method_a, method_b in PAIRS:
            label = f"{method_a} / {method_b}"
            row = {"pair": label, **{name: counts[label][name] for name in fieldnames[1:]}}
            writer.writerow(row)
            mismatches = row["denominator"] - row["full_match"]
            print(
                f"{label}: mismatch {mismatches:,}; "
                f"metal-only {row['metal_only']:,}; "
                f"ligand-only/coarse-metal-match "
                f"{row['ligand_only_or_same-element-metal_ambiguity']:,}; "
                f"both {row['ligand_and_metal']:,}"
            )
    print(f"Wrote {summary_path}")
    print(f"Wrote {detail_path}")


if __name__ == "__main__":
    main()
