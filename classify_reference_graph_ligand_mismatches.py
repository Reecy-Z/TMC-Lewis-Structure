#!/usr/bin/env python3
"""Classify ligand mismatches after the full reference-graph ILP treatment."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger

from evaluate_ilp_vs_xyz2mol_smiles import extract_features, tm_atoms

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent
EXPERIMENT = ROOT / "reference_graph_ilp_experiment"


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def coordination_motif(smiles: str) -> tuple[str, str]:
    """Return metal and a mutually exclusive, priority-ordered motif class."""
    mol = Chem.MolFromSmiles(smiles)
    metals = tm_atoms(mol) if mol is not None else []
    if len(metals) != 1:
        return "", "unparsed"
    metal = metals[0]
    neighbors = [
        atom for atom in metal.GetNeighbors() if atom.GetAtomicNum() != 1
    ]
    carbon_donors = [atom for atom in neighbors if atom.GetAtomicNum() == 6]
    carbon_indices = {atom.GetIdx() for atom in carbon_donors}
    has_carbonyl = False
    has_cn = False
    has_carbene = False
    has_haptic = False
    for carbon in carbon_donors:
        ligand_neighbors = [
            atom
            for atom in carbon.GetNeighbors()
            if atom.GetIdx() != metal.GetIdx() and atom.GetAtomicNum() != 1
        ]
        symbols = [atom.GetSymbol() for atom in ligand_neighbors]
        has_carbonyl |= len(ligand_neighbors) == 1 and "O" in symbols
        has_cn |= len(ligand_neighbors) == 1 and "N" in symbols
        has_carbene |= carbon.IsInRing() and "N" in symbols
        has_haptic |= any(
            atom.GetIdx() in carbon_indices for atom in ligand_neighbors
        )

    if has_haptic:
        category = "haptic_or_adjacent_multi_C"
    elif has_carbene:
        category = "heterocyclic_carbene_like"
    elif has_carbonyl:
        category = "carbonyl_like"
    elif has_cn:
        category = "C_N_donor_cyanide_or_isocyanide"
    elif carbon_donors:
        category = "other_C_bound_ligand"
    else:
        donor_symbols = {atom.GetSymbol() for atom in neighbors}
        if donor_symbols <= {"N", "O"}:
            category = "N_O_only_coordination"
        elif "P" in donor_symbols:
            category = "P_donor_no_M_C"
        elif donor_symbols & {"S", "Se"}:
            category = "S_Se_donor_no_M_C"
        elif donor_symbols & {"F", "Cl", "Br", "I"}:
            category = "halide_containing_no_M_C"
        else:
            category = "other_non_C_donor"
    return metal.GetSymbol(), category


def main() -> None:
    reference_rows = {
        row["IDs"]: row for row in load_rows(ROOT / "tmqmg_smiles.csv")
    }
    full_results = {
        row["IDs"]: row
        for row in load_rows(EXPERIMENT / "per_structure_results_full.csv")
    }
    valid_ids = {
        row["IDs"]
        for row in load_rows(EXPERIMENT / "reference_edges.csv")
        if row["validated"] == "1" and row["in_graph_consensus"] == "1"
    }
    agreement = {
        row["IDs"]: row
        for row in load_rows(EXPERIMENT / "agreement_per_structure.csv")
        if row["treatment"] == "full" and row["IDs"] in valid_ids
    }

    motif_totals = Counter()
    motif_mismatches = Counter()
    metal_totals = Counter()
    metal_mismatches = Counter()
    classified: list[dict[str, str | int]] = []
    mechanism_counts = Counter()
    oxidation_delta_counts = Counter()

    for csd_id in sorted(valid_ids):
        reference_smiles = reference_rows[csd_id]["smiles_CSD_fix"]
        metal, motif = coordination_motif(reference_smiles)
        motif_totals[motif] += 1
        metal_totals[metal] += 1
        row = agreement[csd_id]
        if row["ligand_equiv"] != "0":
            continue
        motif_mismatches[motif] += 1
        metal_mismatches[metal] += 1

        if row["both_parsed"] != "1":
            mechanism = "feature_extraction_failure"
            oxidation_delta: str | int = ""
            ilp_ligands = ""
            reference_ligands = ""
        else:
            oxidation_delta = int(row["ox_ilp"]) - int(row["ox_reference"])
            oxidation_delta_counts[oxidation_delta] += 1
            ilp_features = extract_features(
                full_results[csd_id]["smiles_ilp"],
                convert_x=True,
            )
            reference_features = extract_features(reference_smiles, convert_x=False)
            if ilp_features is None or reference_features is None:
                mechanism = "feature_extraction_failure"
                ilp_ligands = ""
                reference_ligands = ""
            else:
                ilp_ligands = " || ".join(ilp_features["ligands"])
                reference_ligands = " || ".join(reference_features["ligands"])
                if oxidation_delta != 0:
                    mechanism = "coupled_ligand_charge_and_oxidation"
                elif (
                    ilp_features["ligands_skeleton"]
                    != reference_features["ligands_skeleton"]
                ):
                    mechanism = "metal_disconnection_or_ligand_partition"
                else:
                    mechanism = "unresolved_positional_or_equivalence"
        mechanism_counts[mechanism] += 1
        classified.append(
            {
                "IDs": csd_id,
                "metal": metal,
                "motif_class": motif,
                "mechanism_class": mechanism,
                "oxidation_delta_ilp_minus_reference": oxidation_delta,
                "full_connectivity_match": row["full_connectivity"],
                "metal_connectivity_match": row["metal_connectivity"],
                "ilp_oxidation": row["ox_ilp"],
                "reference_oxidation": row["ox_reference"],
                "ilp_ligands": ilp_ligands,
                "reference_ligands": reference_ligands,
            }
        )

    write_rows(EXPERIMENT / "ligand_mismatch_classification.csv", classified)
    write_rows(
        EXPERIMENT / "ligand_mismatch_mechanism_summary.csv",
        [
            {"mechanism_class": key, "count": count}
            for key, count in mechanism_counts.most_common()
        ],
    )
    write_rows(
        EXPERIMENT / "ligand_mismatch_oxidation_delta_summary.csv",
        [
            {"oxidation_delta_ilp_minus_reference": key, "count": count}
            for key, count in oxidation_delta_counts.most_common()
        ],
    )
    motif_summary = [
        {
            "motif_class": motif,
            "cohort_count": total,
            "mismatch_count": motif_mismatches[motif],
            "mismatch_percent": 100.0 * motif_mismatches[motif] / total,
        }
        for motif, total in motif_totals.items()
    ]
    motif_summary.sort(key=lambda row: float(row["mismatch_percent"]), reverse=True)
    write_rows(EXPERIMENT / "ligand_mismatch_by_motif.csv", motif_summary)

    metal_summary = [
        {
            "metal": metal,
            "cohort_count": total,
            "mismatch_count": metal_mismatches[metal],
            "mismatch_percent": 100.0 * metal_mismatches[metal] / total,
        }
        for metal, total in metal_totals.items()
    ]
    metal_summary.sort(key=lambda row: float(row["mismatch_percent"]), reverse=True)
    write_rows(EXPERIMENT / "ligand_mismatch_by_metal.csv", metal_summary)

    canonical_mismatches = sum(
        row["ligand_canonical"] == "0" for row in agreement.values()
    )
    equivalent_mismatches = len(classified)
    summary = {
        "cohort_size": len(valid_ids),
        "canonical_ligand_mismatches": canonical_mismatches,
        "canonical_ligand_mismatch_percent": 100.0
        * canonical_mismatches
        / len(valid_ids),
        "chemically_equivalent_ligand_mismatches": equivalent_mismatches,
        "chemically_equivalent_ligand_mismatch_percent": 100.0
        * equivalent_mismatches
        / len(valid_ids),
        "canonical_only_differences_resolved_by_equivalence": canonical_mismatches
        - equivalent_mismatches,
        "mechanism_counts": dict(mechanism_counts),
    }
    (EXPERIMENT / "ligand_mismatch_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
