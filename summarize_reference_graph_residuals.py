#!/usr/bin/env python3
"""Summarize residual errors after full reference-graph ILP treatment."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger

from evaluate_ilp_vs_xyz2mol_smiles import TM_ATOMIC_NUMS

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent
EXPERIMENT = ROOT / "reference_graph_ilp_experiment"


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_counter(path: Path, fields: list[str], counts: Counter) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([*fields, "count"])
        for key, count in counts.most_common():
            values = key if isinstance(key, tuple) else (key,)
            writer.writerow([*values, count])


def metal_symbol(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    metals = [
        atom.GetSymbol()
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() in TM_ATOMIC_NUMS
    ]
    return metals[0] if len(metals) == 1 else ""


def main() -> None:
    references = {
        row["IDs"]: row
        for row in load_rows(ROOT / "tmqmg_smiles.csv")
    }
    metal_by_id = {
        csd_id: metal_symbol(row["smiles_CSD_fix"])
        for csd_id, row in references.items()
    }
    validated_ids = {
        row["IDs"]
        for row in load_rows(EXPERIMENT / "reference_edges.csv")
        if row["validated"] == "1" and row["in_graph_consensus"] == "1"
    }

    residual_connectivity = load_rows(
        EXPERIMENT / "full_treatment_residual_connectivity.csv"
    )
    metal_counts = Counter()
    donor_counts = Counter()
    enriched = []
    for row in residual_connectivity:
        csd_id = row["IDs"]
        metal = metal_by_id.get(csd_id, "")
        reference_neighbors = Counter(row["neighbors_reference"].split())
        output_neighbors = Counter(row["neighbors_ilp"].split())
        missing = reference_neighbors - output_neighbors
        extra = output_neighbors - reference_neighbors
        metal_counts[metal] += 1
        for donor, count in missing.items():
            donor_counts[(metal, donor)] += count
        enriched.append(
            {
                **row,
                "metal": metal,
                "missing_reference_neighbors": " ".join(missing.elements()),
                "extra_output_neighbors": " ".join(extra.elements()),
            }
        )

    if enriched:
        with (EXPERIMENT / "full_treatment_residual_connectivity_enriched.csv").open(
            "w",
            newline="",
            encoding="utf-8",
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=list(enriched[0]))
            writer.writeheader()
            writer.writerows(enriched)
    write_counter(
        EXPERIMENT / "residual_connectivity_by_metal.csv",
        ["metal"],
        metal_counts,
    )
    write_counter(
        EXPERIMENT / "residual_missing_donor_by_metal.csv",
        ["metal", "missing_donor"],
        donor_counts,
    )

    agreement = load_rows(EXPERIMENT / "agreement_per_structure.csv")
    oxidation_counts = Counter()
    ligand_counts = Counter()
    for row in agreement:
        if (
            row["treatment"] != "full"
            or row["IDs"] not in validated_ids
            or row["both_parsed"] != "1"
        ):
            continue
        metal = metal_by_id.get(row["IDs"], "")
        if row["metal_oxidation"] == "0":
            delta = int(row["ox_ilp"]) - int(row["ox_reference"])
            oxidation_counts[(metal, delta)] += 1
        if row["ligand_equiv"] == "0":
            ligand_counts[metal] += 1
    write_counter(
        EXPERIMENT / "residual_oxidation_by_metal_delta.csv",
        ["metal", "ilp_minus_reference_oxidation"],
        oxidation_counts,
    )
    write_counter(
        EXPERIMENT / "residual_ligand_by_metal.csv",
        ["metal"],
        ligand_counts,
    )


if __name__ == "__main__":
    main()
