#!/usr/bin/env python3
"""Bar chart of ligand / metal-connectivity / oxidation match on the 28,948 cohort."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from evaluate_ilp_vs_xyz2mol_smiles import compare_pair, extract_features
from plot_ilp_xyz2mol_time_histograms import (
    COLOR_AXIS,
    COLOR_ILP,
    FIG_DPI,
    _apply_nature_style,
    _style_nature_axes,
)

ROOT = Path(__file__).resolve().parent
EDGE_CSV = ROOT / "reference_graph_ilp_experiment" / "reference_edges.csv"
AGREE_CSV = ROOT / "reference_graph_ilp_experiment" / "agreement_per_structure.csv"
REF_SMILES = ROOT / "tmqmg_smiles.csv"
RERUN_CSVS = (
    ROOT / "reference_graph_ilp_experiment" / "hypervalent_degree_rerun" / "per_structure_results.csv",
    ROOT / "reference_graph_ilp_experiment" / "carbene_o2_rerun" / "per_structure_results.csv",
    ROOT / "reference_graph_ilp_experiment" / "ox_min_rerun" / "per_structure_results.csv",
)
OUT_PNG = ROOT / "reference_connectivity_match_bars.png"
OUT_JSON = ROOT / "reference_connectivity_match_bars_summary.json"

METRICS = (
    ("ligand_equiv", "Ligand"),
    ("metal_connectivity", "Metal connectivity"),
    ("metal_oxidation", "Metal oxidation"),
)


def load_cohort() -> set[str]:
    return {
        row["IDs"]
        for row in csv.DictReader(EDGE_CSV.open(encoding="utf-8"))
        if row.get("validated") == "1" and row.get("in_graph_consensus") == "1"
    }


def load_original_flags(cohort: set[str]) -> dict[str, dict[str, int]]:
    out: dict[str, dict[str, int]] = {}
    for row in csv.DictReader(AGREE_CSV.open(encoding="utf-8")):
        if row.get("treatment") != "full" or row["IDs"] not in cohort:
            continue
        out[row["IDs"]] = {
            key: int(row.get(key, "0") or 0)
            for key, _label in METRICS
        }
    return out


def load_latest_smiles() -> dict[str, str]:
    latest: dict[str, str] = {}
    for path in RERUN_CSVS:
        for row in csv.DictReader(path.open(encoding="utf-8")):
            if row.get("status") == "ok" and row.get("smiles_ilp"):
                latest[row["IDs"]] = row["smiles_ilp"]
    return latest


def main() -> int:
    cohort = load_cohort()
    flags = load_original_flags(cohort)
    refs = {
        row["IDs"]: row["smiles_CSD_fix"]
        for row in csv.DictReader(REF_SMILES.open(encoding="utf-8"))
    }
    latest = {cid: smi for cid, smi in load_latest_smiles().items() if cid in cohort}

    n_reparsed = 0
    n_reparse_fail = 0
    for cid, smiles in latest.items():
        ilp = extract_features(smiles, convert_x=True)
        ref = extract_features(refs[cid], convert_x=False)
        if ilp is None or ref is None:
            n_reparse_fail += 1
            flags[cid] = {key: 0 for key, _label in METRICS}
            continue
        cmp = compare_pair(ilp, ref)
        flags[cid] = {key: int(bool(cmp[key])) for key, _label in METRICS}
        n_reparsed += 1

    n = len(cohort)
    rates = {}
    for key, label in METRICS:
        n_match = sum(flags[cid][key] for cid in cohort)
        rates[key] = {
            "label": label,
            "n_match": n_match,
            "n_mismatch": n - n_match,
            "percent": 100.0 * n_match / n,
        }

    _apply_nature_style()
    fig, ax = plt.subplots(figsize=(3.6, 2.55))
    x = np.arange(len(METRICS))
    heights = [rates[key]["percent"] for key, _label in METRICS]
    ax.bar(
        x,
        heights,
        width=0.62,
        color=COLOR_ILP,
        edgecolor="white",
        linewidth=0.35,
        alpha=0.72,
    )
    ax.set_xticks(x)
    ax.set_xticklabels([label for _key, label in METRICS])
    ax.set_ylabel("Match rate (%)")
    ax.set_ylim(0.0, 100.0)
    ax.set_yticks([0, 25, 50, 75, 100])
    for xi, height in zip(x, heights):
        ax.text(
            xi,
            min(height + 1.8, 97.5),
            f"{height:.2f}%",
            ha="center",
            va="bottom",
            fontsize=7,
            color=COLOR_AXIS,
        )
    _style_nature_axes(ax)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=FIG_DPI, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)

    summary = {
        "n_complexes": n,
        "n_reparsed_from_reruns": n_reparsed,
        "n_reparse_fail": n_reparse_fail,
        "rates": rates,
        "output_png": str(OUT_PNG),
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
