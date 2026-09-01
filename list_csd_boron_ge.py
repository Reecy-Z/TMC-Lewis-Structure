#!/usr/bin/env python3
"""
List CSD codes whose heaviest-fragment formula has >= N boron atoms.

Parses ``formula_heaviest_fragment`` in ``tmqmg_smiles.csv`` (space-separated
element tokens like ``B10``, ``Br1``). Bromine ``Br*`` is not counted as boron.

Example:
  python list_csd_boron_ge.py --min-b 6 -o csd_boron_ge6.txt
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

DEFAULT_CSV = Path(__file__).resolve().parent / "tmqmg_smiles.csv"
# Element token: symbol + optional count (e.g. C40, B10, Cl2, Br1).
_TOKEN_RE = re.compile(r"^([A-Za-z]+)(\d*)$")


def boron_count_from_formula(formula: str) -> int:
    """Count boron atoms from a Hill-style fragment string."""
    if not formula or not str(formula).strip():
        return 0
    total = 0
    for raw in str(formula).split():
        m = _TOKEN_RE.match(raw.strip())
        if not m:
            continue
        sym, count_s = m.group(1), m.group(2)
        if sym == "Br":
            continue
        if sym == "B":
            total += int(count_s) if count_s else 1
    return total


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--csv",
        type=Path,
        default=DEFAULT_CSV,
        help=f"tmQM-G style CSV (default: {DEFAULT_CSV})",
    )
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "csd_codes_boron_ge6.txt",
        help="Output .txt path (one CSD code per line)",
    )
    ap.add_argument(
        "--min-b",
        type=int,
        default=6,
        help="Minimum boron count (inclusive)",
    )
    args = ap.parse_args()

    if args.min_b < 0:
        print("--min-b must be >= 0", file=sys.stderr)
        return 1
    if not args.csv.is_file():
        print(f"CSV not found: {args.csv}", file=sys.stderr)
        return 1

    codes: list[str] = []
    with args.csv.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None or "IDs" not in reader.fieldnames:
            print("CSV must have an 'IDs' column.", file=sys.stderr)
            return 1
        if "formula_heaviest_fragment" not in reader.fieldnames:
            print("CSV must have 'formula_heaviest_fragment'.", file=sys.stderr)
            return 1
        for row in reader:
            csd = (row.get("IDs") or "").strip()
            if not csd:
                continue
            formula = row.get("formula_heaviest_fragment") or ""
            if boron_count_from_formula(formula) >= args.min_b:
                codes.append(csd)

    codes.sort()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(codes) + ("\n" if codes else ""), encoding="utf-8")
    print(f"Wrote {len(codes)} CSD codes to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
