#!/usr/bin/env python3
"""Copy tmQMg XYZ files that land on heatmap gray cells; filename includes total charge."""

from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from plot_mo_mlx_heatmap import (  # noqa: E402
    CSV_PATH,
    DEFAULT_XYZ_DIR,
    EN_MAX,
    EN_MIN,
    METAL_ORDER,
    VN_MAX,
    VN_MIN,
    classify_metal_centres,
    is_forbidden_mlx_cell,
    load_engine,
    load_skip_ids,
    mlx_label,
)

OUT_DIR = ROOT / "gray_cell_xyz"
LIMITS = {"Cr": None, "Mo": 15, "W": 15}


def charge_tag(charge: int) -> str:
    return f"{charge:+d}"


def iter_candidates(metal: str):
    skip = load_skip_ids()
    for row in csv.DictReader(CSV_PATH.open(encoding="utf-8")):
        if row["IDs"] in skip or row["status"] != "ok":
            continue
        if not re.search(rf"\[{metal}", row.get("smiles_ilp", "")):
            continue
        charge = int(float(row["charge"] or 0))
        if charge % 2 == 0:
            continue
        yield row["IDs"], charge


def write_xyz(src: Path, dest: Path, comment: str) -> None:
    lines = src.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValueError(f"empty xyz: {src}")
    n = lines[0].strip()
    rest = lines[2:] if len(lines) >= 2 else []
    dest.write_text("\n".join([n, comment, *rest]) + "\n", encoding="utf-8")


def main() -> int:
    xyz_dir = DEFAULT_XYZ_DIR
    if not xyz_dir.is_dir():
        raise SystemExit(f"XYZ directory not found: {xyz_dir}")
    engine = load_engine()
    OUT_DIR.mkdir(exist_ok=True)
    manifest = []

    for metal in METAL_ORDER:
        metal_dir = OUT_DIR / metal
        metal_dir.mkdir(exist_ok=True)
        limit = LIMITS[metal]
        kept = 0
        seen_cells: dict[tuple[int, int], int] = defaultdict(int)

        for csd_id, charge in iter_candidates(metal):
            if limit is not None and kept >= limit:
                break
            src = xyz_dir / f"{csd_id}.xyz"
            if not src.is_file():
                continue
            try:
                centres = classify_metal_centres(engine, src, charge, metal)
            except Exception as exc:
                print(f"  skip {csd_id}: {type(exc).__name__}: {exc}")
                continue
            gray_centres = [
                (l, x, en, vn)
                for l, x, en, vn in centres
                if EN_MIN <= en <= EN_MAX
                and VN_MIN <= vn <= VN_MAX
                and is_forbidden_mlx_cell(en, vn)
            ]
            if not gray_centres:
                continue
            l, x, en, vn = gray_centres[0]
            tag = charge_tag(charge)
            dest = metal_dir / f"{csd_id}_charge_{tag}.xyz"
            comment = (
                f"{csd_id} {metal} total_charge={tag} "
                f"EN={en} VN={vn} {metal}{mlx_label(l, x)} heatmap_gray_cell"
            )
            write_xyz(src, dest, comment)
            kept += 1
            seen_cells[(en, vn)] += 1
            manifest.append(
                f"{metal}\t{csd_id}\t{tag}\t{en}\t{vn}\t{mlx_label(l, x)}\t{dest.name}"
            )
            print(f"  {dest.relative_to(ROOT)}")

        print(f"{metal}: wrote {kept} gray-cell xyz (limit={limit})")

    (OUT_DIR / "manifest.tsv").write_text(
        "metal\tid\tcharge\ten\tvn\tlabel\tfilename\n" + "\n".join(manifest) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {OUT_DIR / 'manifest.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
