#!/usr/bin/env python3
"""Prepare XYZ-indexed Hückel connectivity for the strict consensus cohort."""

from __future__ import annotations

import argparse
import csv
import json
import signal
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from rdkit import Chem, RDLogger

from compare_smiles_connectivity_graphs import connectivity_key

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent
DEFAULT_EVALUATION = ROOT / "ilp_mismatch_rerun_current" / "evaluation" / "per_structure.csv"
DEFAULT_GRAPH_COMPARISON = ROOT / "connectivity_graph_comparison" / "per_structure.csv"
DEFAULT_REFERENCE = ROOT / "tmqmg_smiles.csv"
DEFAULT_XYZ_DIR = Path("/data/jingyuan_data/tmqmg")
DEFAULT_XYZ2MOL_ROOT = Path("/home/zhujingyuan/xyz2mol_tm")
DEFAULT_OUT_DIR = ROOT / "reference_graph_ilp_experiment"

REFERENCE_PAIRS = (
    "CSD / NBO(DFT)",
    "CSD / Hückel(DFT)",
    "NBO(DFT) / Hückel(DFT)",
)
STRICT_METRICS = ("ligand_canonical", "metal_connectivity", "metal_oxidation")


def load_consensus_cohorts(
    evaluation_path: Path,
    graph_comparison_path: Path,
) -> tuple[set[str], set[str]]:
    metric_rows: dict[str, dict[str, dict[str, str]]] = {}
    with evaluation_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["pair"] in REFERENCE_PAIRS:
                metric_rows.setdefault(row["IDs"], {})[row["pair"]] = row

    metric_consensus = {
        csd_id
        for csd_id, rows in metric_rows.items()
        if all(
            pair in rows
            and rows[pair]["both_parsed"] == "1"
            and all(rows[pair][metric] == "1" for metric in STRICT_METRICS)
            for pair in REFERENCE_PAIRS
        )
    }

    graph_rows: dict[str, dict[str, dict[str, str]]] = {}
    with graph_comparison_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["IDs"] in metric_consensus and row["pair"] in REFERENCE_PAIRS:
                graph_rows.setdefault(row["IDs"], {})[row["pair"]] = row

    graph_consensus = {
        csd_id
        for csd_id in metric_consensus
        if all(
            graph_rows.get(csd_id, {}).get(pair, {}).get("both_parsed") == "1"
            and graph_rows[csd_id][pair]["connectivity_match"] == "1"
            for pair in REFERENCE_PAIRS
        )
    }
    return metric_consensus, graph_consensus


def load_reference_rows(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["IDs"]: row for row in csv.DictReader(handle)}


def parse_charge(value: str) -> int:
    return int(round(float(value or 0)))


def encode_edges(edges: list[tuple[int, int]]) -> str:
    return ";".join(f"{i}-{j}" for i, j in edges)


def _heavy_connectivity_key(
    symbols: list[str],
    edges: list[tuple[int, int]],
) -> tuple[tuple[int, ...], tuple[tuple[int, int], ...]] | None:
    periodic_table = Chem.GetPeriodicTable()
    heavy_indices = [idx for idx, symbol in enumerate(symbols) if symbol != "H"]
    old_to_new = {old: new for new, old in enumerate(heavy_indices)}
    editable = Chem.RWMol()
    for old_idx in heavy_indices:
        atom = Chem.Atom(periodic_table.GetAtomicNumber(symbols[old_idx]))
        atom.SetNoImplicit(True)
        editable.AddAtom(atom)
    for begin, end in edges:
        if begin in old_to_new and end in old_to_new:
            editable.AddBond(old_to_new[begin], old_to_new[end], Chem.BondType.SINGLE)
    graph = editable.GetMol()
    graph.UpdatePropertyCache(strict=False)
    try:
        ranks = list(
            Chem.CanonicalRankAtoms(
                graph,
                breakTies=True,
                includeChirality=False,
                includeIsotopes=False,
            )
        )
    except Exception:
        return None
    order = sorted(range(graph.GetNumAtoms()), key=lambda idx: ranks[idx])
    canonical_index = {old: new for new, old in enumerate(order)}
    atoms = tuple(graph.GetAtomWithIdx(idx).GetAtomicNum() for idx in order)
    canonical_edges = tuple(
        sorted(
            (
                min(canonical_index[bond.GetBeginAtomIdx()], canonical_index[bond.GetEndAtomIdx()]),
                max(canonical_index[bond.GetBeginAtomIdx()], canonical_index[bond.GetEndAtomIdx()]),
            )
            for bond in graph.GetBonds()
        )
    )
    return atoms, canonical_edges


def _timeout_handler(_signum, _frame) -> None:
    raise TimeoutError("Hückel reference extraction exceeded 300 seconds")


def extract_reference_graph(
    job: tuple[str, int, str, str, str],
) -> dict[str, str | int]:
    csd_id, charge, huckel_smiles, xyz_dir, xyz2mol_root = job
    xyz_path = Path(xyz_dir) / f"{csd_id}.xyz"
    result: dict[str, str | int] = {
        "IDs": csd_id,
        "charge": charge,
        "status": "error",
        "validated": 0,
        "extraction_method": "",
        "n_atoms": 0,
        "n_edges": 0,
        "n_heavy_edges": 0,
        "edges_0": "",
        "error": "",
    }
    try:
        if not xyz_path.is_file():
            raise FileNotFoundError(xyz_path)
        sys.path.insert(0, xyz2mol_root)
        from xyz2mol_tm.huckel_to_smiles.xyz2mol_tmc import (
            get_tmc_mol,
            get_tmc_xyz_edges,
        )

        signal.signal(signal.SIGALRM, _timeout_handler)
        signal.alarm(300)
        try:
            edges = get_tmc_xyz_edges(str(xyz_path), charge)
        finally:
            signal.alarm(0)

        lines = xyz_path.read_text(encoding="utf-8").splitlines()
        symbols = [line.split()[0] for line in lines[2:] if line.strip()]
        edges = sorted(set((min(int(i), int(j)), max(int(i), int(j))) for i, j in edges))
        if any(i < 0 or j >= len(symbols) or i >= j for i, j in edges):
            raise ValueError("reference edge index outside XYZ atom range")

        mapped_key = _heavy_connectivity_key(symbols, edges)
        stored_key = connectivity_key(huckel_smiles)
        validated = mapped_key is not None and mapped_key == stored_key
        extraction_method = "fast_binary"
        if not validated:
            signal.alarm(300)
            try:
                _mol, edges = get_tmc_mol(
                    str(xyz_path),
                    charge,
                    return_xyz_edges=True,
                )
            finally:
                signal.alarm(0)
            edges = sorted(
                set((min(int(i), int(j)), max(int(i), int(j))) for i, j in edges)
            )
            mapped_key = _heavy_connectivity_key(symbols, edges)
            validated = mapped_key is not None and mapped_key == stored_key
            extraction_method = "full_huckel_fallback"
        result.update(
            {
                "status": "ok" if validated else "validation_failed",
                "validated": int(validated),
                "extraction_method": extraction_method,
                "n_atoms": len(symbols),
                "n_edges": len(edges),
                "n_heavy_edges": sum(
                    symbols[i] != "H" and symbols[j] != "H" for i, j in edges
                ),
                "edges_0": encode_edges(edges),
                "error": "" if validated else "mapped graph differs from stored Hückel SMILES",
            }
        )
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


def write_cohort(path: Path, ids: set[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["IDs"])
        writer.writerows((csd_id,) for csd_id in sorted(ids))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluation", type=Path, default=DEFAULT_EVALUATION)
    parser.add_argument("--graph-comparison", type=Path, default=DEFAULT_GRAPH_COMPARISON)
    parser.add_argument("--reference-csv", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--xyz-dir", type=Path, default=DEFAULT_XYZ_DIR)
    parser.add_argument("--xyz2mol-root", type=Path, default=DEFAULT_XYZ2MOL_ROOT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    metric_consensus, graph_consensus = load_consensus_cohorts(
        args.evaluation,
        args.graph_comparison,
    )
    references = load_reference_rows(args.reference_csv)
    selected = sorted(metric_consensus)
    if args.limit:
        selected = selected[: args.limit]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_cohort(args.out_dir / "cohort_metric_consensus_29110.csv", metric_consensus)
    write_cohort(args.out_dir / "cohort_graph_consensus_29068.csv", graph_consensus)
    write_cohort(
        args.out_dir / "cohort_metric_only_ambiguous_42.csv",
        metric_consensus - graph_consensus,
    )

    jobs = [
        (
            csd_id,
            parse_charge(references[csd_id]["charge"]),
            references[csd_id]["smiles_huckel_DFT_xyz"],
            str(args.xyz_dir),
            str(args.xyz2mol_root),
        )
        for csd_id in selected
    ]
    if args.workers <= 1:
        results = []
        for index, job in enumerate(jobs, 1):
            results.append(extract_reference_graph(job))
            if index % 250 == 0:
                print(f"prepared {index} / {len(jobs)}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            results = []
            for index, result in enumerate(
                pool.map(extract_reference_graph, jobs, chunksize=4),
                1,
            ):
                results.append(result)
                if index % 250 == 0:
                    print(f"prepared {index} / {len(jobs)}", flush=True)

    fields = [
        "IDs",
        "charge",
        "in_graph_consensus",
        "status",
        "validated",
        "extraction_method",
        "n_atoms",
        "n_edges",
        "n_heavy_edges",
        "edges_0",
        "error",
    ]
    output_path = args.out_dir / "reference_edges.csv"
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in results:
            row["in_graph_consensus"] = int(row["IDs"] in graph_consensus)
            writer.writerow(row)

    status_counts: dict[str, int] = {}
    for row in results:
        status = str(row["status"])
        status_counts[status] = status_counts.get(status, 0) + 1
    summary = {
        "metric_consensus": len(metric_consensus),
        "graph_consensus": len(graph_consensus),
        "metric_only_ambiguous": len(metric_consensus - graph_consensus),
        "processed": len(results),
        "status_counts": status_counts,
    }
    (args.out_dir / "reference_graph_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
