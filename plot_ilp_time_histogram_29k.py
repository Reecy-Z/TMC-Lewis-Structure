#!/usr/bin/env python3
"""ILP time histogram for the 28,948 reference-graph complexes (no ILP rerun)."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from plot_ilp_xyz2mol_time_histograms import (
    FIG_DPI,
    ILP_CSV,
    _add_mean_median_lines_ilp,
    _apply_nature_style,
    _plot_ilp_histogram,
)

ROOT = Path(__file__).resolve().parent
EDGE_CSV = ROOT / "reference_graph_ilp_experiment" / "reference_edges.csv"
FULL_CSV = ROOT / "reference_graph_ilp_experiment" / "per_structure_results_full.csv"
OUT_PNG = ROOT / "ilp_time_histogram_29k.png"
OUT_JSON = ROOT / "ilp_time_histogram_29k_summary.json"


def load_cohort_ids() -> set[str]:
    return {
        row["IDs"]
        for row in csv.DictReader(EDGE_CSV.open(encoding="utf-8"))
        if row.get("validated") == "1" and row.get("in_graph_consensus") == "1"
    }


def load_full_times() -> dict[str, float]:
    return {
        row["IDs"]: float(row["ilp_seconds"])
        for row in csv.DictReader(FULL_CSV.open(encoding="utf-8"))
    }


def main() -> int:
    cohort = load_cohort_ids()
    by_id: dict[str, float] = {}
    for row in csv.DictReader(ILP_CSV.open(encoding="utf-8")):
        if row["status"] == "missing_xyz":
            continue
        by_id[row["IDs"]] = float(row["ilp_seconds"])
    full_times = load_full_times()
    missing = sorted(cohort - set(by_id))
    for csd_id in missing:
        if csd_id not in full_times:
            raise SystemExit(f"No stored ILP time for {csd_id}")
        by_id[csd_id] = full_times[csd_id]

    ids = sorted(cohort)
    arr = np.asarray([by_id[csd_id] for csd_id in ids], dtype=float)

    _apply_nature_style()
    fig, ax = plt.subplots(figsize=(3.5, 2.55))
    _plot_ilp_histogram(ax, arr, panel_title=None)
    _add_mean_median_lines_ilp(ax, arr, symbols_only=True, indexed_x=False)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=FIG_DPI, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)

    summary = {
        "n_complexes": int(arr.size),
        "times_from_benchmark": int(arr.size - len(missing)),
        "times_from_reference_graph_full": missing,
        "ilp_seconds": {
            "median": float(np.median(arr)),
            "mean": float(np.mean(arr)),
            "p95": float(np.percentile(arr, 95)),
            "p99": float(np.percentile(arr, 99)),
            "max": float(np.max(arr)),
        },
        "output_png": str(OUT_PNG),
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
