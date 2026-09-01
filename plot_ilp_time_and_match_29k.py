#!/usr/bin/env python3
"""Two-panel figure: (A) 29k ILP times, (B) reference-connectivity match rates."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from plot_ilp_time_histogram_29k import load_cohort_ids, load_full_times
from plot_ilp_xyz2mol_time_histograms import (
    COLOR_AXIS,
    COLOR_ILP,
    FIG_DPI,
    ILP_CSV,
    _add_mean_median_lines_ilp,
    _apply_nature_style,
    _plot_ilp_histogram,
    _style_nature_axes,
)

ROOT = Path(__file__).resolve().parent
RATES_JSON = ROOT / "reference_connectivity_match_bars_summary.json"
OUT_PNG = ROOT / "ilp_time_and_match_29k.png"

BAR_LABELS = (
    ("ligand_equiv", "Ligand"),
    ("metal_connectivity", "Metal\nconnectivity"),
    ("metal_oxidation", "Metal\noxidation"),
)


def load_times() -> np.ndarray:
    cohort = load_cohort_ids()
    by_id: dict[str, float] = {}
    for row in csv.DictReader(ILP_CSV.open(encoding="utf-8")):
        if row["status"] == "missing_xyz":
            continue
        by_id[row["IDs"]] = float(row["ilp_seconds"])
    full_times = load_full_times()
    for csd_id in cohort:
        if csd_id not in by_id:
            by_id[csd_id] = full_times[csd_id]
    return np.asarray([by_id[csd_id] for csd_id in sorted(cohort)], dtype=float)


def _add_panel_letter(ax, letter: str) -> None:
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    ylab = ax.yaxis.label.get_window_extent(renderer)
    x0, _ = ax.transAxes.inverted().transform((ylab.x0, ylab.y0))
    ax.text(
        x0,
        1.03,
        letter,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontweight="bold",
        fontsize=9,
        color=COLOR_AXIS,
        clip_on=False,
    )


def plot_match_bars(ax, rates: dict) -> None:
    x = np.arange(len(BAR_LABELS))
    heights = [rates[key]["percent"] for key, _label in BAR_LABELS]
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
    ax.set_xticklabels([label for _key, label in BAR_LABELS])
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


def main() -> int:
    times = load_times()
    rates = json.loads(RATES_JSON.read_text(encoding="utf-8"))["rates"]

    _apply_nature_style()
    fig, (ax_a, ax_b) = plt.subplots(
        1,
        2,
        figsize=(7.2, 2.65),
        gridspec_kw={"width_ratios": [1.05, 1.0], "wspace": 0.32},
    )
    _plot_ilp_histogram(ax_a, times, panel_title=None)
    _add_mean_median_lines_ilp(ax_a, times, symbols_only=True, indexed_x=False)
    plot_match_bars(ax_b, rates)
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.22, top=0.88, wspace=0.32)
    _add_panel_letter(ax_a, "A")
    _add_panel_letter(ax_b, "B")
    fig.savefig(OUT_PNG, dpi=FIG_DPI, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    print(OUT_PNG)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
