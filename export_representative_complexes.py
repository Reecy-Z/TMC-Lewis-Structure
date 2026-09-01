#!/usr/bin/env python3
"""Export representative tmQM-G complexes into representative_complex/ by category."""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "tmqmg_ilp_benchmark_output" / "per_structure_results.csv"
XYZ_DIR = Path("/data/jingyuan_data/tmqmg")
OUT_ROOT = ROOT / "representative_complex"

CATEGORIES: list[dict] = [
    {
        "folder": "01_werner_type",
        "label": "Werner-type coordination complexes",
        "entries": [
            {
                "id": "ACUZEG",
                "description": "trans-[Co(en)2Cl2]+ (CoIII dichloro bis(ethylenediamine); diamine chelate)",
            },
            {
                "id": "GOLPAZ",
                "description": "Ni(II) square-planar bis(salicylaldimine) (O,N chelate; Werner-type Ni)",
            },
            {
                "id": "JAQPIB",
                "description": "[Cd(NH3)4]2+ with phenanthroline (Werner CdII ammine)",
            },
        ],
    },
    {
        "folder": "02_metal_carbonyl",
        "label": "Metal carbonyl complexes",
        "entries": [
            {
                "id": "KINXEK",
                "description": "Cr(CO)5 plus methyl isocyanide (pentacarbonyl chromium)",
            },
            {
                "id": "DAMKUY",
                "description": "Mo(CO)5 with phenyl isocyanide",
            },
            {
                "id": "HAHRUD",
                "description": "W(CO)5 with thioformate ligand",
            },
        ],
    },
    {
        "folder": "03_metal_alkene_alkyne",
        "label": "Metal-alkene/alkyne complexes",
        "entries": [
            {
                "id": "KCEYPT",
                "description": "[PtCl3(eta2-C2H4)]- Zeise-type anion",
            },
            {
                "id": "DOXPEL",
                "description": "Pt(II) eta2-alkene with norbornene framework",
            },
            {
                "id": "XUDXUQ",
                "description": "Pt(II) complex with terminal alkyne ligands",
            },
        ],
    },
    {
        "folder": "04_metallocenes_sandwich",
        "label": "Metallocenes and (half-)sandwich complexes",
        "entries": [
            {
                "id": "FONGEV",
                "subtype": "sandwich",
                "description": "Fe(II) bis(eta5-indolyl) sandwich",
            },
            {
                "id": "BEGXOB",
                "subtype": "sandwich",
                "description": "Fe(II) bis(eta5-indolyl) sandwich with Se coligands",
            },
            {
                "id": "NUCZUI",
                "subtype": "half-sandwich",
                "description": "Cr(eta5-Cp) piano-stool with CO",
            },
            {
                "id": "IWOLIP",
                "subtype": "half-sandwich",
                "description": "Cr(eta6-C6H5F) piano-stool with three CO ligands",
            },
        ],
    },
    {
        "folder": "05_metal_carbene_alkylidene",
        "label": "Metal carbene and alkylidene complexes",
        "entries": [
            {
                "id": "KUNROA",
                "subtype": "Fischer",
                "description": "W(CO)4 Fischer alkoxycarbene (minimal 38-atom COQMOM-style M=C in tmQM-G benchmark)",
            },
            {
                "id": "DARFEH",
                "subtype": "Fischer",
                "description": "Cr(CO)5 Fischer alkoxycarbene (30-atom textbook (CO)5M=C(OMe)R prototype in tmQM-G benchmark)",
            },
            {
                "id": "COQMOM",
                "subtype": "Fischer",
                "description": "Re Fischer-type alkoxycarbene with four CO ligands",
            },
            {
                "id": "COQNAZ",
                "subtype": "Fischer",
                "description": "Re Fischer-type alkoxycarbene (para-tolyl variant)",
            },
            {
                "id": "COQNIH",
                "subtype": "Fischer",
                "description": "Re Fischer-type alkoxycarbene with fluorenyl substituent",
            },
            {
                "id": "EBOWIB",
                "subtype": "Fischer",
                "description": "Mo(VI) difluoro oxo alkoxycarbene chelate (minimal 19-atom Fischer; terminal M=O)",
            },
            {
                "id": "ZOGCON",
                "subtype": "Schrock",
                "description": "Mo(VI) Cp benzylidene alkylidene (42-atom Schrock M=C; ILP M–C double bond)",
            },
            {
                "id": "REDWAY",
                "subtype": "Schrock",
                "description": "Mo(VI) Cp dialkyl alkylidene (32-atom Schrock M=C; I + Cp + two M=C alkyl groups)",
            },
            {
                "id": "RIPXIA",
                "subtype": "NHC",
                "description": "Cu(I) chloro complex with dialkyl imidazol-2-ylidene (IPr-type NHC)",
            },
            {
                "id": "ECIHOO",
                "subtype": "NHC",
                "description": "Au(I) chloro imidazol-2-ylidene (minimal 17-atom NHC in tmQM-G benchmark)",
            },
            {
                "id": "ABICUK",
                "subtype": "NHC",
                "description": "Rh(I) bromo carbonyl with mesityl imidazol-2-ylidene (IMes-type NHC)",
            },
            {
                "id": "PEHJEU",
                "subtype": "NHC",
                "description": "Re(I) pentacarbonyl bromo with benzimidazol-2-ylidene NHC",
            },
            {
                "id": "KOQBAT",
                "subtype": "CAAC",
                "description": "Cr(0) pentacarbonyl with cyclic alkylamino carbene (minimal 23-atom CAAC in tmQM-G benchmark)",
            },
        ],
    },
    {
        "folder": "06_oxo_nitrido_imido",
        "label": "Oxo, nitrido, and imido complexes",
        "entries": [
            {
                "id": "ABIKOM",
                "subtype": "oxo",
                "description": "Mo(VI) cis-dioxo dichloride (oxo)",
            },
            {
                "id": "BULNED",
                "subtype": "nitrido",
                "description": "Mo(VI) terminal nitrido complex",
            },
            {
                "id": "BOFTOH",
                "subtype": "imido",
                "description": "Mo(VI) oxo-imido with phenylimido ligand",
            },
        ],
    },
    {
        "folder": "07_large_chelate",
        "label": "Large chelate ring complexes",
        "entries": [
            {
                "id": "AKIKUB",
                "description": "Ni(II) tetraaza macrocycle (large N4 chelate ring)",
            },
            {
                "id": "HIWJIH",
                "description": "Cd(II) tetraaza macrocycle",
            },
            {
                "id": "YOQKAR",
                "description": "Ru(II) complex with expanded N4 chelate framework",
            },
        ],
    },
    {
        "folder": "other",
        "label": "Other coordination motifs",
        "entries": [
            {
                "id": "AKOSEZ",
                "subtype": "agostic",
                "description": "Cr(CO)5 with B–H agostic donation from NHC–BH3 (30-atom; CBC B–H σ-complex L)",
            },
            {
                "id": "SAJCIP",
                "subtype": "oxo",
                "description": "W(VI) oxo tetrafluoride pyridine adduct (17-atom; terminal M=O → X2)",
            },
            {
                "id": "PAXDEA",
                "subtype": "hypervalent",
                "description": "Pt(II) trichloro DMSO (14-atom; sulfoxide S 10e expanded octet → L)",
            },
            {
                "id": "ENSPAU",
                "subtype": "hypervalent",
                "description": "Au(III) bis(sulfonate) (21-atom; sulfonate S 12e expanded octet → L)",
            },
            {
                "id": "CEFHAW",
                "subtype": "hypervalent",
                "description": "Pt(II) chelate with sulfonyl/nitro S (25-atom; S 12e expanded octet → X)",
            },
            {
                "id": "TAYBED",
                "subtype": "hexadentate",
                "description": "Co(III) hexadentate phenolic macrocycle (85-atom single ligand; 4N+2O → L4X2)",
            },
            {
                "id": "NOEPOR",
                "subtype": "tetradentate",
                "description": "Ni(II) expanded porphyrinoid macrocycle (85-atom single ligand; 4×N → X4)",
            },
            {
                "id": "ALOKOB",
                "subtype": "tetradentate",
                "description": "Zn(II) porphyrinoid macrocycle (85-atom single ligand; 4×N → L2X2)",
            },
            {
                "id": "JAWCAO",
                "subtype": "tetradentate",
                "description": "Ag porphyrinoid macrocycle (84-atom single ligand; 4×N → LX3; ALOKOB Ag analog)",
            },
            {
                "id": "BAYSUR",
                "subtype": "tetradentate",
                "description": "Ag porphyrinoid macrocycle (65-atom F-substituted; 4×N → LX3)",
            },
            {
                "id": "FIQPIG",
                "subtype": "tetradentate",
                "description": "Au carbaporphyrinoid macrocycle (81-atom single ligand; 3×N+1×C → LX3)",
            },
            {
                "id": "BAQROB",
                "subtype": "tetradentate",
                "description": "Ag carbaporphyrinoid macrocycle (56-atom single ligand; 3×N+1×C → X4)",
            },
            {
                "id": "EJIZUS",
                "subtype": "hexadentate",
                "description": "Y(III) hexadentate amine macrocycle with Si arms (85-atom single ligand; 6×N)",
            },
        ],
    },
]


def charge_tag(charge: int) -> str:
    if charge > 0:
        return f"plus{charge}"
    if charge < 0:
        return f"minus{abs(charge)}"
    return "0"


def main() -> None:
    rows = {
        row["IDs"]: row
        for row in csv.DictReader(CSV_PATH.open(encoding="utf-8"))
        if row["status"] == "ok"
    }

    manifest: list[dict] = []
    missing_xyz: list[str] = []
    missing_csv: list[str] = []

    for category in CATEGORIES:
        folder = OUT_ROOT / category["folder"]
        folder.mkdir(parents=True, exist_ok=True)

        for entry in category["entries"]:
            csd_id = entry["id"]
            row = rows.get(csd_id)
            if row is None:
                missing_csv.append(csd_id)
                continue

            charge = int(row["charge"])
            tag = charge_tag(charge)
            xyz_name = f"{csd_id}_charge_{tag}.xyz"
            src = XYZ_DIR / f"{csd_id}.xyz"
            dst = folder / xyz_name

            if not src.is_file():
                missing_xyz.append(csd_id)
                continue

            shutil.copy2(src, dst)

            manifest.append(
                {
                    "id": csd_id,
                    "category_folder": category["folder"],
                    "category": category["label"],
                    "subtype": entry.get("subtype"),
                    "charge": charge,
                    "charge_tag": tag,
                    "xyz": str(dst.relative_to(OUT_ROOT)),
                    "description": entry["description"],
                    "status": row["status"],
                    "smiles_ilp": row["smiles_ilp"],
                }
            )

    manifest_path = OUT_ROOT / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print(f"Exported {len(manifest)} structures to {OUT_ROOT}")
    if missing_csv:
        print("Missing from benchmark CSV:", ", ".join(missing_csv))
    if missing_xyz:
        print("Missing XYZ files:", ", ".join(missing_xyz))


if __name__ == "__main__":
    main()
