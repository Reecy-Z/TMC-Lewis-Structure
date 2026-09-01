#!/usr/bin/env python3
"""
Re-test Lewis-engine-ILP on XYZ files in tmqmg_ilp_benchmark_output/error/.

Skips CSD codes in list_error_geometry.txt and csd_codes_boron_ge6.txt.
Copies still-failing XYZ to tmqmg_ilp_benchmark_output/error_cannot_be_saved/.
"""

from __future__ import annotations

import importlib.util
import json
import re
import shutil
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ENGINE_PATH = ROOT / "Lewis-engine-ILP.py"
ERROR_DIR = ROOT / "tmqmg_ilp_benchmark_output" / "error"
OUT_FAIL_DIR = ROOT / "tmqmg_ilp_benchmark_output" / "error_cannot_be_saved"
SKIP_PATHS = [
    ROOT / "list_error_geometry.txt",
    ROOT / "csd_codes_boron_ge6.txt",
]

FNAME_RE = re.compile(r"^(.+)_(-?\d+)_(.+)\.xyz$", re.IGNORECASE)


@dataclass
class Stats:
    total_xyz: int = 0
    skipped_csd: int = 0
    tested: int = 0
    ilp_ok: int = 0
    ilp_failed: int = 0
    other_error: int = 0
    saved_fail: int = 0
    failures: list[dict] = field(default_factory=list)


def load_skip_set() -> frozenset[str]:
    out: set[str] = set()
    for path in SKIP_PATHS:
        if not path.is_file():
            continue
        for line in path.read_text(encoding="utf-8").splitlines():
            code = line.split("#", 1)[0].strip()
            if code:
                out.add(code.upper())
    return frozenset(out)


def load_engine():
    spec = importlib.util.spec_from_file_location("lewis_engine_ilp", ENGINE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {ENGINE_PATH}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def parse_xyz_entry(path: Path) -> tuple[str, int] | None:
    m = FNAME_RE.match(path.name)
    if not m:
        return None
    return m.group(1).upper(), int(m.group(2))


def run_one(engine, xyz_path: Path, charge: int) -> tuple[str, str]:
    try:
        atoms = engine.read_xyz(str(xyz_path))
        raw = engine.connectivity(atoms)
        aromatic_systems = engine.aromatic_candidate_systems(atoms, raw)
        engine.solve_bond_orders(
            atoms,
            raw,
            aromatic_systems,
            mol_charge=charge,
            metal_adjacency_edges=raw,
        )
        return "ok", ""
    except RuntimeError as exc:
        msg = str(exc)
        if "ILP failed" in msg:
            return "ilp_failed", msg
        return "error", msg
    except Exception as exc:
        return "error", f"{type(exc).__name__}: {exc}"


def main() -> int:
    if not ERROR_DIR.is_dir():
        print(f"Error dir not found: {ERROR_DIR}", file=sys.stderr)
        return 1
    if not ENGINE_PATH.is_file():
        print(f"Engine not found: {ENGINE_PATH}", file=sys.stderr)
        return 1

    skip = load_skip_set()
    engine = load_engine()
    OUT_FAIL_DIR.mkdir(parents=True, exist_ok=True)

    xyz_files = sorted(ERROR_DIR.glob("*.xyz"))
    stats = Stats(total_xyz=len(xyz_files))
    t0 = time.perf_counter()

    try:
        from tqdm import tqdm

        iterator = tqdm(xyz_files, desc="Retest error XYZ", unit="file")
    except ImportError:
        iterator = xyz_files

    for xyz_path in iterator:
        parsed = parse_xyz_entry(xyz_path)
        if parsed is None:
            stats.other_error += 1
            stats.failures.append(
                {"file": xyz_path.name, "status": "bad_filename", "error": "unparsed name"}
            )
            dest = OUT_FAIL_DIR / xyz_path.name
            shutil.copy2(xyz_path, dest)
            stats.saved_fail += 1
            continue

        csd_id, charge = parsed
        if csd_id in skip:
            stats.skipped_csd += 1
            continue

        stats.tested += 1
        status, err = run_one(engine, xyz_path, charge)
        if status == "ok":
            stats.ilp_ok += 1
            continue

        if status == "ilp_failed":
            stats.ilp_failed += 1
        else:
            stats.other_error += 1

        stats.failures.append(
            {
                "file": xyz_path.name,
                "csd_id": csd_id,
                "charge": charge,
                "status": status,
                "error": err[:500],
            }
        )
        dest = OUT_FAIL_DIR / xyz_path.name
        shutil.copy2(xyz_path, dest)
        stats.saved_fail += 1

    wall = time.perf_counter() - t0
    rate = (100.0 * stats.ilp_ok / stats.tested) if stats.tested else 0.0

    summary = {
        "error_dir": str(ERROR_DIR),
        "out_fail_dir": str(OUT_FAIL_DIR),
        "skip_lists": [str(p) for p in SKIP_PATHS],
        "skip_csd_count": len(skip),
        "total_xyz_in_error_dir": stats.total_xyz,
        "skipped_csd": stats.skipped_csd,
        "tested": stats.tested,
        "ilp_ok": stats.ilp_ok,
        "ilp_failed": stats.ilp_failed,
        "other_error": stats.other_error,
        "success_rate_percent": round(rate, 2),
        "saved_to_error_cannot_be_saved": stats.saved_fail,
        "wall_seconds": round(wall, 2),
    }

    summary_json = OUT_FAIL_DIR / "retest_summary.json"
    summary_txt = OUT_FAIL_DIR / "retest_summary.txt"
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    lines = [
        "Retest Lewis-engine-ILP on tmqmg_ilp_benchmark_output/error/",
        "=" * 60,
        f"Total XYZ in error/:           {stats.total_xyz}",
        f"Skipped (CSD skip lists):      {stats.skipped_csd}",
        f"Tested:                        {stats.tested}",
        f"ILP success:                   {stats.ilp_ok}",
        f"ILP failed:                    {stats.ilp_failed}",
        f"Other errors:                  {stats.other_error}",
        f"Success rate (of tested):      {rate:.2f}%",
        f"Saved to error_cannot_be_saved: {stats.saved_fail}",
        f"Wall time (s):                 {wall:.2f}",
        "",
        f"Summary JSON: {summary_json}",
    ]
    summary_txt.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print()
    print(summary_txt.read_text(encoding="utf-8"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
