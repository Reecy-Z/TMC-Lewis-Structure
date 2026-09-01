#!/usr/bin/env python3
"""Draw remaining ILP/reference ligand mismatches in the first10 two-column style."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor

ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = (
    ROOT
    / "reference_graph_ilp_experiment"
    / "hypervalent_degree_rerun"
    / "still_mismatch.csv"
)
DEFAULT_OUTPUT_DIR = (
    ROOT
    / "reference_graph_ilp_experiment"
    / "hypervalent_degree_rerun"
    / "first100_pages"
)
DEFAULT_PDF = (
    ROOT
    / "reference_graph_ilp_experiment"
    / "hypervalent_degree_rerun"
    / "first100_ligand_mismatch_ILP_vs_reference.pdf"
)


def split_ligands(value: str) -> list[str]:
    return [item.strip() for item in value.split(" || ") if item.strip()]


def unmatched_ligands(ilp_value: str, reference_value: str) -> tuple[list[str], list[str]]:
    ilp = Counter(split_ligands(ilp_value))
    reference = Counter(split_ligands(reference_value))
    common = ilp & reference
    return list((ilp - common).elements()), list((reference - common).elements())


def oxidation_delta(row: dict[str, str]) -> str:
    try:
        return str(int(row["ox_ilp"]) - int(row["ox_reference"]))
    except (KeyError, TypeError, ValueError):
        return "?"


def combined_mol(smiles_list: list[str]) -> Chem.Mol:
    if not smiles_list:
        mol = Chem.MolFromSmiles("*")
    else:
        mol = Chem.MolFromSmiles(".".join(smiles_list))
    if mol is None:
        mol = Chem.MolFromSmiles("*")
    rdDepictor.Compute2DCoords(mol)
    return mol


def page_image(
    rows: list[dict[str, str]],
    *,
    sub_img_size: tuple[int, int] = (720, 300),
    legend_font_size: int = 16,
):
    molecules: list[Chem.Mol] = []
    legends: list[str] = []
    for row in rows:
        ilp_ligands, reference_ligands = unmatched_ligands(
            row.get("ilp_ligands", ""),
            row.get("reference_ligands", ""),
        )
        molecules.extend(
            [
                combined_mol(ilp_ligands),
                combined_mol(reference_ligands),
            ]
        )
        delta = oxidation_delta(row)
        legends.extend(
            [
                f"{row['IDs']}  |  ILP unmatched ligand(s)  |  Δox={delta}",
                f"{row['IDs']}  |  Reference unmatched ligand(s)",
            ]
        )
    options = Draw.MolDrawOptions()
    options.legendFontSize = legend_font_size
    options.minFontSize = max(8, legend_font_size // 2)
    options.bondLineWidth = 2.0
    return Draw.MolsToGridImage(
        molecules,
        molsPerRow=2,
        subImgSize=sub_img_size,
        legends=legends,
        useSVG=False,
        returnPNG=False,
        drawOptions=options,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--pdf", type=Path, default=DEFAULT_PDF)
    parser.add_argument("--limit", type=int, default=100)
    parser.add_argument("--page-size", type=int, default=10)
    parser.add_argument("--width", type=int, default=720, help="Per-cell width in pixels")
    parser.add_argument("--height", type=int, default=300, help="Per-cell height in pixels")
    parser.add_argument("--legend-font-size", type=int, default=16)
    parser.add_argument("--dpi", type=float, default=72.0, help="PDF raster DPI")
    args = parser.parse_args()

    with args.input.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))[: args.limit]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    images = []
    for page_index, start in enumerate(range(0, len(rows), args.page_size), start=1):
        page_rows = rows[start : start + args.page_size]
        image = page_image(
            page_rows,
            sub_img_size=(args.width, args.height),
            legend_font_size=args.legend_font_size,
        )
        images.append(image)
        path = args.output_dir / f"p{page_index:02d}.png"
        image.save(path, dpi=(args.dpi, args.dpi))
        print(f"Wrote {path}", flush=True)

    if images:
        rgb_images = [image.convert("RGB") for image in images]
        rgb_images[0].save(
            args.pdf,
            save_all=True,
            append_images=rgb_images[1:],
            resolution=args.dpi,
            dpi=(args.dpi, args.dpi),
        )
        print(f"Wrote {args.pdf} ({len(rows)} complexes, {len(images)} pages)")


if __name__ == "__main__":
    main()
