#!/usr/bin/env python3
"""
Batch-run Lewis-engine-ILP on tmQM-G structures and validate RDKit SMILES.

Reads CSD IDs and molecular charge from tmqmg_smiles.csv, loads matching
``{ID}.xyz`` from a directory (default: /data/jingyuan_data/tmqmg), runs the
aromatic ILP + RDKit SMILES export, writes per-structure results and a summary,
and copies input XYZ for failed cases into ``tmqmg_ilp_benchmark_output/error/``
as ``{ID}_{charge}_{reason}.xyz`` (e.g. ``ABAXUX_0_rdkit_parse_failed.xyz``).

Example:
  python batch_tmqmg_ilp_benchmark.py --limit 100
  python batch_tmqmg_ilp_benchmark.py  # full ~60k set (may take hours)

Progress bar (tqdm if installed) and incremental writes:
``per_structure_results.csv``, ``summary.json``, and ``summary.txt`` are updated
during the run (each row to CSV; summary every ``--summary-every`` structures).

CSD codes listed in skip lists are not run: they do not count toward ``--limit``,
batch totals, or failure statistics, and are omitted from
``per_structure_results.csv``.

By default, if present, the script will load:
- ``csd_codes_boron_ge6.txt`` (e.g. boron-rich exclude list)
- ``list_error_geometry.txt`` (user-curated geometry/problem exclude list)
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import shutil
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

DEFAULT_CSV = Path(__file__).with_name("tmqmg_smiles.csv")
DEFAULT_XYZ_DIR = Path("/data/jingyuan_data/tmqmg")
DEFAULT_ENGINE = Path(__file__).with_name("Lewis-engine-ILP.py")
DEFAULT_OUT_DIR = Path(__file__).with_name("tmqmg_ilp_benchmark_output")
DEFAULT_SKIP_CSD_LIST = Path(__file__).with_name("csd_codes_boron_ge6.txt")
DEFAULT_SKIP_ERROR_GEOMETRY_LIST = Path(__file__).with_name("list_error_geometry.txt")

DETAIL_CSV_COLUMNS = [
    "IDs",
    "charge",
    "status",
    "ilp_seconds",
    "n_atoms",
    "rdkit_parse_ok",
    "error",
    "smiles_ilp",
]


@dataclass
class RowResult:
    csd_id: str
    charge: int
    status: str
    ilp_seconds: float = 0.0
    smiles: str = ""
    rdkit_parse_ok: bool = False
    error: str = ""
    n_atoms: int = 0


@dataclass
class OxSigmaFlag:
    csd_id: str
    charge: int
    flags: list[str]



@dataclass
class BatchCounters:
    total_rows: int = 0
    skipped_excluded: int = 0
    missing_xyz: int = 0
    ilp_success: int = 0
    ilp_failed: int = 0
    smiles_generated: int = 0
    rdkit_parse_ok: int = 0
    rdkit_parse_fail: int = 0
    other_errors: int = 0
    error_xyz_saved: int = 0


def load_lewis_engine(engine_path: Path):
    spec = importlib.util.spec_from_file_location("lewis_engine_ilp", engine_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load engine from {engine_path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def parse_charge(raw: str) -> int:
    if raw is None or str(raw).strip() == "":
        return 0
    return int(round(float(raw)))


def load_csd_skip_set(path: Path | None) -> frozenset[str]:
    """One CSD code per line; ``#`` starts an end-of-line comment. Empty if missing."""
    if path is None:
        return frozenset()
    return load_csd_skip_set_from_paths([path])


def load_csd_skip_set_from_paths(paths: list[Path]) -> frozenset[str]:
    """Union of one-code-per-line lists; silently ignores missing files."""
    out: set[str] = set()
    for path in paths:
        if not path.is_file():
            continue
        for line in path.read_text(encoding="utf-8").splitlines():
            s = line.split("#", 1)[0].strip()
            if s:
                out.add(s)
    return frozenset(out)


def run_ilp_for_xyz(
    engine: Any,
    xyz_path: Path,
    mol_charge: int,
    *,
    rdkit_Chem: Any | None,
) -> RowResult:
    csd_id = xyz_path.stem
    t0 = time.perf_counter()
    try:
        atoms = engine.read_xyz(str(xyz_path))
        raw = engine.connectivity(atoms)
        aromatic_systems = engine.aromatic_candidate_systems(atoms, raw)
        bonds, lp_out, fc_out = engine.solve_bond_orders(
            atoms,
            raw,
            aromatic_systems,
            mol_charge=mol_charge,
            metal_adjacency_edges=raw,
        )
        bonds, lp_out, fc_out, _carbene_labels = engine.apply_heterocyclic_carbene_corrections(
            atoms,
            bonds,
            lp_out,
            fc_out,
            aromatic_systems,
            raw,
            mol_charge=mol_charge,
        )
        bonds, lp_out, fc_out = engine.apply_eta_covalent_pi_corrections(
            atoms,
            bonds,
            lp_out,
            fc_out,
            mol_charge=mol_charge,
            metal_adjacency_edges=raw,
        )
        ilp_sec = time.perf_counter() - t0
        coords = [[a[2], a[3], a[4]] for a in atoms]
        atom_syms = [a[1] for a in atoms]
        fc_list = [fc_out.get(a[0], 0) for a in atoms]
        idx_to_pos = {a[0]: k for k, a in enumerate(atoms)}
        metal_adj_0 = [
            (idx_to_pos[tm], idx_to_pos[lig], ei, ej)
            for tm, lig, ei, ej in raw
            if engine.is_TM(ei) ^ engine.is_TM(ej)
        ]
        dative_ml = engine.infer_dative_ml_pairs_cbc(
            atom_syms,
            coords,
            bonds,
            lp_out,
            fc_list,
            metal_adjacency_edges=metal_adj_0,
        )
        smiles = engine.ilp_to_smiles(
            atoms,
            bonds,
            fc_out,
            dative_ml_pairs=dative_ml,
            edges=raw,
        )
        rdkit_ok = False
        if rdkit_Chem is not None and smiles:
            rdkit_ok = rdkit_Chem.MolFromSmiles(smiles) is not None
        return RowResult(
            csd_id=csd_id,
            charge=mol_charge,
            status="ok",
            ilp_seconds=ilp_sec,
            smiles=smiles,
            rdkit_parse_ok=rdkit_ok,
            n_atoms=len(atoms),
        )
    except RuntimeError as exc:
        ilp_sec = time.perf_counter() - t0
        msg = str(exc)
        if "ILP failed" in msg:
            return RowResult(
                csd_id=csd_id,
                charge=mol_charge,
                status="ilp_failed",
                ilp_seconds=ilp_sec,
                error=msg,
                n_atoms=0,
            )
        return RowResult(
            csd_id=csd_id,
            charge=mol_charge,
            status="error",
            ilp_seconds=ilp_sec,
            error=msg,
            n_atoms=0,
        )
    except Exception as exc:
        return RowResult(
            csd_id=csd_id,
            charge=mol_charge,
            status="error",
            ilp_seconds=time.perf_counter() - t0,
            error=f"{type(exc).__name__}: {exc}",
            n_atoms=0,
        )


def ox_lt_sigma_flags(engine: Any, atoms: list, bonds: list[tuple[int, int, int]], fc_out: dict[int, int]) -> list[str]:
    """
    Return a list of strings describing TM centres where ox(TM) < Σ b_tm(TM–L).

    ox(TM) is taken from fc_out entries for TM atoms (solve_bond_orders convention).
    Σ b_tm is computed from ILP bond list 'bonds' on TM–nonmetal edges with order>0.
    """
    atom_el = {i: el for i, el, *_ in atoms}
    sigma_by_tm: dict[int, int] = {}
    for a, b, order in bonds:
        ea, eb = atom_el[a], atom_el[b]
        if not engine._is_tm_nm_edge(ea, eb):
            continue
        tm = a if engine.is_TM(ea) else b
        sigma_by_tm[tm] = sigma_by_tm.get(tm, 0) + int(order)

    flags: list[str] = []
    for i, el in atom_el.items():
        if not engine.is_TM(el):
            continue
        ox = fc_out.get(i)
        if ox is None:
            continue
        sigma = sigma_by_tm.get(i, 0)
        if int(ox) < int(sigma):
            flags.append(f"{el}{i}: ox={int(ox)} < sigma={int(sigma)}")
    return flags


def run_ilp_for_xyz_probe_remote_c_relax(
    engine: Any,
    xyz_path: Path,
    mol_charge: int,
    *,
    rdkit_Chem: Any | None,
) -> tuple[RowResult, bool]:
    """
    Two-pass probe for remote-C lp penalty:

    - Pass 1: default ILP_WEIGHT_REMOTE_C_LP_VIOLATION (elastic prefer lp=0)
    - Pass 2: weight=0 (no remote-C lp penalty) if pass 1 is ILP-infeasible

    Returns (RowResult, used_second_pass).
    """
    csd_id = xyz_path.stem
    t0 = time.perf_counter()

    atoms = engine.read_xyz(str(xyz_path))
    raw = engine.connectivity(atoms)
    aromatic_systems = engine.aromatic_candidate_systems(atoms, raw)

    def _solve(*, remote_c_lp_violation_weight: float | None):
        return engine.solve_bond_orders(
            atoms,
            raw,
            aromatic_systems,
            mol_charge=mol_charge,
            metal_adjacency_edges=raw,
            remote_c_lp_violation_weight=remote_c_lp_violation_weight,
        )

    try:
        bonds, lp_out, fc_out = _solve(remote_c_lp_violation_weight=None)
        used_second = False
    except RuntimeError as exc:
        msg = str(exc)
        if "ILP failed" not in msg:
            ilp_sec = time.perf_counter() - t0
            return (
                RowResult(
                    csd_id=csd_id,
                    charge=mol_charge,
                    status="error",
                    ilp_seconds=ilp_sec,
                    error=msg,
                    n_atoms=0,
                ),
                False,
            )
        # Pass 2
        try:
            bonds, lp_out, fc_out = _solve(remote_c_lp_violation_weight=0.0)
            used_second = True
        except RuntimeError as exc2:
            ilp_sec = time.perf_counter() - t0
            msg2 = str(exc2)
            status = "ilp_failed" if "ILP failed" in msg2 else "error"
            return (
                RowResult(
                    csd_id=csd_id,
                    charge=mol_charge,
                    status=status,
                    ilp_seconds=ilp_sec,
                    error=msg2,
                    n_atoms=0,
                ),
                False,
            )

    bonds, lp_out, fc_out, _carbene_labels = engine.apply_heterocyclic_carbene_corrections(
        atoms,
        bonds,
        lp_out,
        fc_out,
        aromatic_systems,
        raw,
        mol_charge=mol_charge,
    )
    bonds, lp_out, fc_out = engine.apply_eta_covalent_pi_corrections(
        atoms,
        bonds,
        lp_out,
        fc_out,
        mol_charge=mol_charge,
        metal_adjacency_edges=raw,
    )

    ilp_sec = time.perf_counter() - t0
    coords = [[a[2], a[3], a[4]] for a in atoms]
    atom_syms = [a[1] for a in atoms]
    fc_list = [fc_out.get(a[0], 0) for a in atoms]
    idx_to_pos = {a[0]: k for k, a in enumerate(atoms)}
    metal_adj_0 = [
        (idx_to_pos[tm], idx_to_pos[lig], ei, ej)
        for tm, lig, ei, ej in raw
        if engine.is_TM(ei) ^ engine.is_TM(ej)
    ]
    dative_ml = engine.infer_dative_ml_pairs_cbc(
        atom_syms,
        coords,
        bonds,
        lp_out,
        fc_list,
        metal_adjacency_edges=metal_adj_0,
    )
    smiles = engine.ilp_to_smiles(
        atoms,
        bonds,
        fc_out,
        dative_ml_pairs=dative_ml,
        edges=raw,
    )
    rdkit_ok = False
    if rdkit_Chem is not None and smiles:
        rdkit_ok = rdkit_Chem.MolFromSmiles(smiles) is not None
    return (
        RowResult(
            csd_id=csd_id,
            charge=mol_charge,
            status="ok",
            ilp_seconds=ilp_sec,
            smiles=smiles,
            rdkit_parse_ok=rdkit_ok,
            n_atoms=len(atoms),
        ),
        used_second,
    )


def read_csv_rows(
    csv_path: Path,
    *,
    start: int,
    limit: int | None,
    exclude_ids: frozenset[str],
) -> tuple[list[tuple[str, int]], int]:
    """Return ``(rows_to_run, n_skipped_excluded)``.

    Rows whose ``IDs`` appear in ``exclude_ids`` are omitted: they do not fill
    ``--limit`` and are not run.
    """
    rows: list[tuple[str, int]] = []
    skipped_excluded = 0
    with csv_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"No header in {csv_path}")
        if "IDs" not in reader.fieldnames or "charge" not in reader.fieldnames:
            raise ValueError(
                f"CSV must contain columns 'IDs' and 'charge'; got {reader.fieldnames}"
            )
        for i, row in enumerate(reader):
            if i < start:
                continue
            if limit is not None and len(rows) >= limit:
                break
            csd = (row.get("IDs") or "").strip()
            if not csd:
                continue
            if csd in exclude_ids:
                skipped_excluded += 1
                continue
            rows.append((csd, parse_charge(row.get("charge", "0"))))
    return rows, skipped_excluded


def is_error_result(r: RowResult) -> bool:
    """True when the run should archive the input XYZ under out_dir/error/."""
    if r.status == "missing_xyz":
        return False
    if r.status != "ok":
        return True
    return not r.rdkit_parse_ok


def error_xyz_reason(r: RowResult) -> str:
    """Fail-reason tag: ilp_failed | rdkit_parse_failed | error."""
    if r.status == "ilp_failed":
        return "ilp_failed"
    if r.status == "ok" and not r.rdkit_parse_ok:
        return "rdkit_parse_failed"
    return "error"


def error_xyz_filename(csd_id: str, r: RowResult) -> str:
    """{CSD}_{charge}_{reason}.xyz"""
    return f"{csd_id}_{r.charge}_{error_xyz_reason(r)}.xyz"


def save_error_xyz(xyz_path: Path, error_dir: Path, r: RowResult) -> Path:
    error_dir.mkdir(parents=True, exist_ok=True)
    dest = error_dir / error_xyz_filename(r.csd_id, r)
    shutil.copy2(xyz_path, dest)
    return dest


def update_counters(c: BatchCounters, r: RowResult) -> None:
    c.total_rows += 1
    if r.status == "missing_xyz":
        c.missing_xyz += 1
    elif r.status == "ilp_failed":
        c.ilp_failed += 1
    elif r.status == "ok":
        c.ilp_success += 1
        if r.smiles:
            c.smiles_generated += 1
        if r.rdkit_parse_ok:
            c.rdkit_parse_ok += 1
        else:
            c.rdkit_parse_fail += 1
    else:
        c.other_errors += 1


def detail_row_values(r: RowResult) -> list[Any]:
    return [
        r.csd_id,
        r.charge,
        r.status,
        f"{r.ilp_seconds:.4f}",
        r.n_atoms,
        int(r.rdkit_parse_ok),
        r.error,
        r.smiles,
    ]


def write_detail_csv(path: Path, results: list[RowResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(DETAIL_CSV_COLUMNS)
        for r in results:
            w.writerow(detail_row_values(r))


class LiveDetailCsv:
    """Append one result row at a time (flush after each write)."""

    def __init__(self, path: Path) -> None:
        self.path = path
        self._fh = None
        self._writer: csv.writer | None = None

    def __enter__(self) -> LiveDetailCsv:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._fh = self.path.open("w", newline="", encoding="utf-8")
        self._writer = csv.writer(self._fh)
        self._writer.writerow(DETAIL_CSV_COLUMNS)
        self._fh.flush()
        return self

    def __exit__(self, *exc: object) -> None:
        if self._fh is not None:
            self._fh.flush()
            self._fh.close()
            self._fh = None
            self._writer = None

    def write_row(self, r: RowResult) -> None:
        if self._writer is None or self._fh is None:
            raise RuntimeError("LiveDetailCsv not opened")
        self._writer.writerow(detail_row_values(r))
        self._fh.flush()


def make_progress(
    rows: list[tuple[str, int]],
    *,
    disable: bool,
    desc: str = "ILP batch",
) -> tuple[Any, Any | None]:
    """Return ``(row_iterator, tqdm_bar_or_none)``."""
    if disable:
        return rows, None
    try:
        from tqdm import tqdm
    except ImportError:
        return rows, None
    bar = tqdm(rows, total=len(rows), desc=desc, unit="struct", dynamic_ncols=True)
    return bar, bar


def write_summary(
    path_txt: Path,
    path_json: Path,
    counters: BatchCounters,
    *,
    wall_seconds: float,
    csv_path: Path,
    xyz_dir: Path,
    engine_path: Path,
    error_xyz_dir: Path,
    start: int,
    limit: int | None,
    skip_list_paths: list[Path],
) -> None:
    skip_list_paths_norm = [str(p) for p in skip_list_paths]
    summary = {
        "input_csv": str(csv_path),
        "xyz_directory": str(xyz_dir),
        "engine": str(engine_path),
        "start_offset": start,
        "limit": limit,
        "skip_csd_lists": skip_list_paths_norm,
        "skipped_excluded": counters.skipped_excluded,
        "wall_time_seconds": round(wall_seconds, 3),
        "total_processed": counters.total_rows,
        "missing_xyz": counters.missing_xyz,
        "ilp_success": counters.ilp_success,
        "ilp_failed": counters.ilp_failed,
        "other_errors": counters.other_errors,
        "smiles_generated": counters.smiles_generated,
        "rdkit_parse_ok": counters.rdkit_parse_ok,
        "rdkit_parse_fail": counters.rdkit_parse_fail,
        "rdkit_parse_ok_fraction_of_ilp_success": (
            round(counters.rdkit_parse_ok / counters.ilp_success, 4)
            if counters.ilp_success
            else None
        ),
        "error_xyz_dir": str(error_xyz_dir),
        "error_xyz_saved": counters.error_xyz_saved,
    }
    path_json.parent.mkdir(parents=True, exist_ok=True)
    path_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    lines = [
        "tmQM-G ILP + RDKit SMILES batch summary",
        "=" * 60,
        f"Input CSV:        {csv_path}",
        f"XYZ directory:    {xyz_dir}",
        f"Engine:           {engine_path}",
        f"Skip lists:       {', '.join(skip_list_paths_norm) if skip_list_paths_norm else '(none)'}",
        f"Skipped (list):   {counters.skipped_excluded}",
        f"Rows processed:   {counters.total_rows}",
        f"Wall time (s):    {wall_seconds:.2f}",
        f"Wall time (min):  {wall_seconds / 60:.2f}",
        "",
        f"Missing XYZ:      {counters.missing_xyz}",
        f"ILP success:      {counters.ilp_success}",
        f"ILP failed:       {counters.ilp_failed}",
        f"Other errors:     {counters.other_errors}",
        "",
        f"SMILES generated: {counters.smiles_generated}",
        f"RDKit parse OK:   {counters.rdkit_parse_ok}",
        f"RDKit parse fail: {counters.rdkit_parse_fail}",
        "",
        f"Error XYZ dir:    {error_xyz_dir}",
        f"Error XYZ saved:  {counters.error_xyz_saved}",
    ]
    if counters.ilp_success:
        pct = 100.0 * counters.rdkit_parse_ok / counters.ilp_success
        lines.append(f"RDKit OK / ILP OK:  {pct:.2f}%")
    path_txt.write_text("\n".join(lines) + "\n", encoding="utf-8")


def progress_postfix(counters: BatchCounters) -> dict[str, int]:
    return {
        "ok": counters.ilp_success,
        "rdkit": counters.rdkit_parse_ok,
        "fail": counters.ilp_failed,
        "miss": counters.missing_xyz,
    }


def set_progress_postfix(pbar: Any, counters: BatchCounters) -> None:
    if pbar is None:
        return
    try:
        pbar.set_postfix(**progress_postfix(counters), refresh=False)
    except Exception:
        pass


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    ap.add_argument("--xyz-dir", type=Path, default=DEFAULT_XYZ_DIR)
    ap.add_argument("--engine", type=Path, default=DEFAULT_ENGINE)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    ap.add_argument("--start", type=int, default=0, help="Skip first N CSV data rows")
    ap.add_argument("--limit", type=int, default=None, help="Max rows to process")
    ap.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable tqdm progress bar",
    )
    ap.add_argument(
        "--summary-every",
        type=int,
        default=1,
        help="Rewrite summary.json / summary.txt every N structures (0 = only at end)",
    )
    ap.add_argument(
        "--progress-every",
        type=int,
        default=0,
        help="Legacy: print a line every N structures (0 = off; tqdm used by default)",
    )
    ap.add_argument(
        "--skip-csd-list",
        type=Path,
        default=None,
        metavar="PATH",
        help=(
            "Text file: one CSD code per line to skip (not run, not in totals/failures). "
            f"If omitted, use {DEFAULT_SKIP_CSD_LIST.name} and/or "
            f"{DEFAULT_SKIP_ERROR_GEOMETRY_LIST.name} when those files exist."
        ),
    )
    ap.add_argument(
        "--no-skip-csd-list",
        action="store_true",
        help="Do not load any CSD skip list.",
    )
    ap.add_argument(
        "--report-remote-c-relax",
        action="store_true",
        help=(
            "Probe which structures require disabling ILP_WEIGHT_REMOTE_C_LP_VIOLATION: "
            "run ILP with default remote-C lp penalty, then rerun with weight=0 "
            "only if needed. Writes IDs to out-dir/remote_c_lp_relax_needed.txt."
        ),
    )
    ap.add_argument(
        "--report-ox-lt-sigma",
        action="store_true",
        help=(
            "Report structures where ox(TM) < Σ b_tm(TM–L) using ILP bond orders. "
            "Writes out-dir/ox_lt_sigma.json and copies matching XYZ with charge comment "
            "to out-dir/ox_lt_sigma_xyz/."
        ),
    )
    args = ap.parse_args()

    if not args.csv.is_file():
        print(f"CSV not found: {args.csv}", file=sys.stderr)
        return 1
    if not args.engine.is_file():
        print(f"Engine not found: {args.engine}", file=sys.stderr)
        return 1
    if not args.xyz_dir.is_dir():
        print(f"XYZ directory not found: {args.xyz_dir}", file=sys.stderr)
        return 1

    try:
        from rdkit import Chem as rdkit_Chem
    except ImportError:
        print("RDKit required for SMILES validation.", file=sys.stderr)
        return 1

    skip_list_paths: list[Path] = []
    if not args.no_skip_csd_list:
        if args.skip_csd_list is not None:
            if args.skip_csd_list.is_file():
                skip_list_paths.append(args.skip_csd_list)
            else:
                print(f"Skip list not found (ignored): {args.skip_csd_list}", file=sys.stderr)
        else:
            if DEFAULT_SKIP_CSD_LIST.is_file():
                skip_list_paths.append(DEFAULT_SKIP_CSD_LIST)
            if DEFAULT_SKIP_ERROR_GEOMETRY_LIST.is_file():
                skip_list_paths.append(DEFAULT_SKIP_ERROR_GEOMETRY_LIST)

    exclude_ids = load_csd_skip_set_from_paths(skip_list_paths)
    if skip_list_paths and exclude_ids:
        print(
            f"Excluding {len(exclude_ids)} CSD code(s) from "
            + ", ".join(str(p) for p in skip_list_paths)
            + " (not counted in --limit, totals, or failures).",
            flush=True,
        )

    rows, skipped_excluded = read_csv_rows(
        args.csv,
        start=args.start,
        limit=args.limit,
        exclude_ids=exclude_ids,
    )
    if not rows:
        print(
            "No rows to process (check --start/--limit or all rows excluded).",
            file=sys.stderr,
        )
        return 1

    print(f"Loading {args.engine.name} ...")
    engine = load_lewis_engine(args.engine)
    counters = BatchCounters()
    counters.skipped_excluded = skipped_excluded
    args.out_dir.mkdir(parents=True, exist_ok=True)
    detail_path = args.out_dir / "per_structure_results.csv"
    summary_txt = args.out_dir / "summary.txt"
    summary_json = args.out_dir / "summary.json"
    error_xyz_dir = args.out_dir / "error"
    remote_c_relax_needed: list[str] = []
    ox_lt_sigma: list[OxSigmaFlag] = []
    t_wall = time.perf_counter()

    def flush_summary() -> None:
        write_summary(
            summary_txt,
            summary_json,
            counters,
            wall_seconds=time.perf_counter() - t_wall,
            csv_path=args.csv,
            xyz_dir=args.xyz_dir,
            engine_path=args.engine,
            error_xyz_dir=error_xyz_dir,
            start=args.start,
            limit=args.limit,
            skip_list_paths=skip_list_paths,
        )

    flush_summary()
    print(f"Live output: {detail_path}", flush=True)
    print(f"Live summary: {summary_json}", flush=True)

    row_iter, pbar = make_progress(rows, disable=args.no_progress)
    n_total = len(rows)

    with LiveDetailCsv(detail_path) as live_csv:
        for n, (csd_id, charge) in enumerate(row_iter, start=1):
            xyz_path = args.xyz_dir / f"{csd_id}.xyz"
            if not xyz_path.is_file():
                r = RowResult(
                    csd_id=csd_id,
                    charge=charge,
                    status="missing_xyz",
                    error=f"not found: {xyz_path}",
                )
            else:
                if args.report_remote_c_relax:
                    r, used_second = run_ilp_for_xyz_probe_remote_c_relax(
                        engine, xyz_path, charge, rdkit_Chem=rdkit_Chem
                    )
                    if used_second and r.status == "ok":
                        remote_c_relax_needed.append(csd_id)
                else:
                    r = run_ilp_for_xyz(
                        engine, xyz_path, charge, rdkit_Chem=rdkit_Chem
                    )
                if args.report_ox_lt_sigma and r.status == "ok":
                    atoms0 = engine.read_xyz(str(xyz_path))
                    raw0 = engine.connectivity(atoms0)
                    arom0 = engine.aromatic_candidate_systems(atoms0, raw0)
                    bonds0, lp0, fc0 = engine.solve_bond_orders(
                        atoms0, raw0, arom0, mol_charge=charge, metal_adjacency_edges=raw0
                    )
                    flags = ox_lt_sigma_flags(engine, atoms0, bonds0, fc0)
                    if flags:
                        ox_lt_sigma.append(
                            OxSigmaFlag(csd_id=csd_id, charge=charge, flags=flags)
                        )
                        out_xyz_dir = args.out_dir / "ox_lt_sigma_xyz"
                        out_xyz_dir.mkdir(parents=True, exist_ok=True)
                        dst = out_xyz_dir / f"{csd_id}_{charge}.xyz"
                        lines = xyz_path.read_text(encoding="utf-8").splitlines(True)
                        if len(lines) >= 2:
                            lines[1] = f"charge={charge}\n"
                        dst.write_text("".join(lines), encoding="utf-8")

            live_csv.write_row(r)
            update_counters(counters, r)
            if is_error_result(r) and xyz_path.is_file():
                save_error_xyz(xyz_path, error_xyz_dir, r)
                counters.error_xyz_saved += 1

            se = args.summary_every
            if se > 0 and (n % se == 0 or n == n_total):
                flush_summary()

            set_progress_postfix(pbar, counters)

            if args.progress_every and n % args.progress_every == 0:
                elapsed = time.perf_counter() - t_wall
                print(
                    f"[{n}/{n_total}] elapsed {elapsed:.0f}s | "
                    f"ILP ok {counters.ilp_success} | RDKit ok {counters.rdkit_parse_ok} | "
                    f"fail {counters.ilp_failed} | missing {counters.missing_xyz}"
                    f" | skipped(list) {counters.skipped_excluded}",
                    flush=True,
                )

    if pbar is not None:
        pbar.close()

    flush_summary()

    print()
    print(summary_txt.read_text(encoding="utf-8"))
    print(f"Detail CSV: {detail_path}")
    print(f"Summary JSON: {summary_json}")
    if counters.error_xyz_saved:
        print(f"Error XYZ:  {error_xyz_dir} ({counters.error_xyz_saved} file(s))")
    if args.report_remote_c_relax:
        out_path = args.out_dir / "remote_c_lp_relax_needed.txt"
        out_path.write_text("\n".join(remote_c_relax_needed) + ("\n" if remote_c_relax_needed else ""), encoding="utf-8")
        print(f"Remote-C lp relax needed: {len(remote_c_relax_needed)}")
        print(f"List saved: {out_path}")
    if args.report_ox_lt_sigma:
        out_json = args.out_dir / "ox_lt_sigma.json"
        out_json.write_text(
            json.dumps([o.__dict__ for o in ox_lt_sigma], indent=2) + "\n",
            encoding="utf-8",
        )
        print(f"Ox < sigma: {len(ox_lt_sigma)}")
        print(f"List saved: {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
