#!/usr/bin/env python3
"""Pick a few tmQMg XYZ files outside the EN 12–18 / VN 0–6 heatmap window.

EN is the charge-free count: EN = 6 + 2n + VN.
"""

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
    METAL_VALENCE,
    VN_MAX,
    VN_MIN,
    classify_metal_centres,
    load_engine,
    load_skip_ids,
    mlx_label,
)

OUT_DIR = ROOT / "out_of_window_xyz"
# A few examples per reason, across metals.
QUOTA = {"EN>18": 3, "EN<12": 3, "VN>6": 3}


def charge_tag(charge: int) -> str:
    return f"{charge:+d}"


def reason_for(en: int, vn: int) -> str | None:
    tags = []
    if en > EN_MAX:
        tags.append("EN>18")
    elif en < EN_MIN:
        tags.append("EN<12")
    if vn > VN_MAX or vn < VN_MIN:
        tags.append("VN>6")
    if not tags:
        return None
    return tags[0]


def iter_rows(metal: str, *, reverse_size: bool = False):
    skip = load_skip_ids()
    rows = []
    for row in csv.DictReader(CSV_PATH.open(encoding="utf-8")):
        if row["IDs"] in skip or row["status"] != "ok":
            continue
        if not re.search(rf"\[{metal}", row.get("smiles_ilp", "")):
            continue
        rows.append(row)
    rows.sort(
        key=lambda r: int(float(r.get("n_atoms") or 999)),
        reverse=reverse_size,
    )
    for row in rows:
        yield row["IDs"], int(float(row["charge"] or 0))


def write_xyz(src: Path, dest: Path, comment: str) -> None:
    lines = src.read_text(encoding="utf-8").splitlines()
    n = lines[0].strip()
    rest = lines[2:] if len(lines) >= 2 else []
    dest.write_text("\n".join([n, comment, *rest]) + "\n", encoding="utf-8")


def main() -> int:
    xyz_dir = DEFAULT_XYZ_DIR
    if not xyz_dir.is_dir():
        raise SystemExit(f"XYZ directory not found: {xyz_dir}")
    engine = load_engine()
    OUT_DIR.mkdir(exist_ok=True)

    kept = defaultdict(int)
    seen_label = set()
    seen_id = set()
    manifest = ["metal\tid\tcharge\tn\tx\ten\tvn\treason\tfilename"]
    search = (("EN<12", False), ("EN>18", True), ("VN>6", True))

    for reason_wanted, reverse_size in search:
        if kept[reason_wanted] >= QUOTA[reason_wanted]:
            continue
        for metal in METAL_ORDER:
            if kept[reason_wanted] >= QUOTA[reason_wanted]:
                break
            for csd_id, charge in iter_rows(metal, reverse_size=reverse_size):
                if kept[reason_wanted] >= QUOTA[reason_wanted]:
                    break
                if csd_id in seen_id:
                    continue
                src = xyz_dir / f"{csd_id}.xyz"
                if not src.is_file():
                    continue
                try:
                    centres = classify_metal_centres(engine, src, charge, metal)
                except Exception as exc:
                    print(f"  skip {csd_id}: {type(exc).__name__}: {exc}")
                    continue
                for l_count, x_count, _en_old, _vn_old in centres:
                    en = METAL_VALENCE + 2 * l_count + x_count
                    vn = x_count
                    reason = reason_for(en, vn)
                    if reason != reason_wanted:
                        continue
                    label = mlx_label(l_count, x_count) or "bare"
                    key = (reason, metal, label)
                    if key in seen_label:
                        continue
                    seen_label.add(key)
                    seen_id.add(csd_id)
                    tag = charge_tag(charge)
                    dest = OUT_DIR / f"{metal}_{csd_id}_charge_{tag}_EN{en}_VN{vn}.xyz"
                    comment = (
                        f"{csd_id} {metal} total_charge={tag} "
                        f"EN={en} VN={vn} {metal}{label} reason={reason} "
                        f"formula EN=6+2n+VN"
                    )
                    write_xyz(src, dest, comment)
                    kept[reason] += 1
                    manifest.append(
                        f"{metal}\t{csd_id}\t{tag}\t{l_count}\t{x_count}\t{en}\t{vn}\t{reason}\t{dest.name}"
                    )
                    print(f"  {dest.name}  {reason}  {metal}{label}")
                    break

    (OUT_DIR / "manifest.tsv").write_text("\n".join(manifest) + "\n", encoding="utf-8")
    print("kept", dict(kept))
    print(f"Wrote {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
