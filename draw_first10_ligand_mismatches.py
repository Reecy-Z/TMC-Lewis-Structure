#!/usr/bin/env python3
"""Draw the first ten ILP/reference ligand mismatches as a two-column grid."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor

ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = (
    ROOT
    / "reference_graph_ilp_experiment"
    / "ligand_mismatch_classification.csv"
)
DEFAULT_OUTPUT = (
    ROOT
    / "reference_graph_ilp_experiment"
    / "first10_ligand_mismatch_ILP_vs_reference.png"
)


def split_ligands(value: str) -> list[str]:
    return [item.strip() for item in value.split(" || ") if item.strip()]


def unmatched_ligands(ilp_value: str, reference_value: str) -> tuple[list[str], list[str]]:
    ilp = Counter(split_ligands(ilp_value))
    reference = Counter(split_ligands(reference_value))
    common = ilp & reference
    return list((ilp - common).elements()), list((reference - common).elements())


def combined_mol(smiles_list: list[str]) -> Chem.Mol:
    mol = Chem.MolFromSmiles(".".join(smiles_list))
    if mol is None:
        raise ValueError(f"Could not parse ligand SMILES: {smiles_list}")
    rdDepictor.Compute2DCoords(mol)
    return mol


def main() -> None:
    with DEFAULT_INPUT.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))[:10]

    molecules: list[Chem.Mol] = []
    legends: list[str] = []
    for row in rows:
        ilp_ligands, reference_ligands = unmatched_ligands(
            row["ilp_ligands"],
            row["reference_ligands"],
        )
        molecules.extend(
            [
                combined_mol(ilp_ligands),
                combined_mol(reference_ligands),
            ]
        )
        delta = row["oxidation_delta_ilp_minus_reference"]
        legends.extend(
            [
                f"{row['IDs']}  |  ILP unmatched ligand(s)  |  Δox={delta}",
                f"{row['IDs']}  |  Reference unmatched ligand(s)",
            ]
        )

    image = Draw.MolsToGridImage(
        molecules,
        molsPerRow=2,
        subImgSize=(720, 300),
        legends=legends,
        useSVG=False,
        returnPNG=False,
    )
    image.save(DEFAULT_OUTPUT)
    print(f"Wrote {DEFAULT_OUTPUT}")


if __name__ == "__main__":
    main()
