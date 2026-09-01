#!/usr/bin/env python3
"""Export Cr L4X4 tmQMg XYZ files (charge-free MLX: n=4, x=4)."""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from plot_mo_mlx_heatmap import (  # noqa: E402
    CSV_PATH,
    DEFAULT_XYZ_DIR,
    METAL_VALENCE,
    classify_metal_centres,
    load_engine,
    load_skip_ids,
    mlx_label,
)

OUT_DIR = ROOT / "crl4x4_xyz"
METAL = "Cr"


def write_xyz(src: Path, dest: Path, comment: str) -> None:
    lines = src.read_text(encoding="utf-8").splitlines()
    n = lines[0].strip()
    rest = lines[2:] if len(lines) >= 2 else []
    dest.write_text("\n".join([n, comment, *rest]) + "\n", encoding="utf-8")


def main() -> int:
    xyz_dir = DEFAULT_XYZ_DIR
    if not xyz_dir.is_dir():
        raise SystemExit(f"XYZ directory not found: {xyz_dir}")
    skip = load_skip_ids()
    engine = load_engine()
    OUT_DIR.mkdir(exist_ok=True)

    rows = []
    for row in csv.DictReader(CSV_PATH.open(encoding="utf-8")):
        if row["IDs"] in skip or row["status"] != "ok":
            continue
        if not re.search(rf"\[{METAL}", row.get("smiles_ilp", "")):
            continue
        rows.append((row["IDs"], int(float(row["charge"] or 0))))

    manifest = ["id\tcharge\tn\tx\ten\tvn\tlabel\tfilename"]
    kept = 0
    for i, (csd_id, charge) in enumerate(rows, start=1):
        if i % 100 == 0:
            print(f"  scanned {i}/{len(rows)}, kept {kept}", flush=True)
        src = xyz_dir / f"{csd_id}.xyz"
        if not src.is_file():
            continue
        try:
            centres = classify_metal_centres(engine, src, charge, METAL)
        except Exception as exc:
            print(f"  skip {csd_id}: {type(exc).__name__}: {exc}", flush=True)
            continue
        for l_count, x_count, _en_old, _vn in centres:
            if l_count != 4 or x_count != 4:
                continue
            en = METAL_VALENCE + 2 * l_count + x_count
            vn = x_count
            tag = f"{charge:+d}"
            dest = OUT_DIR / f"{csd_id}_charge_{tag}_EN{en}_VN{vn}.xyz"
            comment = (
                f"{csd_id} {METAL} total_charge={tag} "
                f"EN={en} VN={vn} {METAL}{mlx_label(l_count, x_count)}"
            )
            write_xyz(src, dest, comment)
            kept += 1
            manifest.append(
                f"{csd_id}\t{tag}\t{l_count}\t{x_count}\t{en}\t{vn}\t"
                f"{mlx_label(l_count, x_count)}\t{dest.name}"
            )
            print(f"  {dest.name}", flush=True)

    (OUT_DIR / "manifest.tsv").write_text("\n".join(manifest) + "\n", encoding="utf-8")
    print(f"Wrote {kept} CrL4X4 files to {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
