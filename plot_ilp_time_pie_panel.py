#!/usr/bin/env python3
"""Combined panel: ILP runtime histogram (left) and failure-mode pie (right)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

from plot_ilp_success_pie import (
    COLOR_ILP_FAIL,
    COLOR_PARSE_FAIL,
    PIE_LEGEND_LABELS,
    VALUE_COLOR,
    VALUE_FONTSIZE,
    load_ilp_outcome_counts,
    plot_ilp_failure_pie,
)
from plot_ilp_xyz2mol_time_histograms import (
    N_COMPLEXES,
    _add_mean_median_lines_ilp,
    _apply_nature_style,
    _plot_ilp_histogram,
    load_ilp_keys_and_times,
    load_skip_ids,
)

ROOT = Path(__file__).resolve().parent
OUT_PNG = ROOT / "ilp_time_pie_panel.png"
OUT_JSON = ROOT / "ilp_time_pie_panel_summary.json"
FIG_DPI = 300
# Multiplier on the auto-computed pie axes side (1.0 = as large as fits below legend).
PIE_SIZE_SCALE = 1.0
# Fraction of the right column height reserved for the legend row.
PIE_LEGEND_HEIGHT_RATIO = 0.28


def _apply_panel_value_typography() -> None:
    plt.rcParams.update(
        {
            "legend.fontsize": VALUE_FONTSIZE,
            "xtick.labelsize": VALUE_FONTSIZE,
            "ytick.labelsize": VALUE_FONTSIZE,
        }
    )


def _add_panel_pie(
    fig,
    ax_legend,
    ax_pie_area,
    n_parse_fail: int,
    n_ilp_fail: int,
    n_total: int,
    *,
    pie_size_scale: float = PIE_SIZE_SCALE,
) -> None:
    """Legend in the top row; pie centered on the legend's horizontal midpoint below it."""
    if pie_size_scale <= 0:
        raise ValueError("pie_size_scale must be positive")

    ax_legend.axis("off")
    ax_pie_area.axis("off")
    legend_handles = [
        Patch(facecolor=COLOR_PARSE_FAIL, edgecolor="white", linewidth=0.8, alpha=0.88),
        Patch(facecolor=COLOR_ILP_FAIL, edgecolor="white", linewidth=0.8, alpha=0.88),
    ]
    legend = ax_legend.legend(
        legend_handles,
        PIE_LEGEND_LABELS,
        loc="upper right",
        frameon=False,
        fontsize=VALUE_FONTSIZE,
    )
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    leg_bbox = legend.get_window_extent(renderer).transformed(fig.transFigure.inverted())
    leg_center_x = 0.5 * (leg_bbox.x0 + leg_bbox.x1)

    pie_box = ax_pie_area.get_position()
    max_side = pie_size_scale * min(
        pie_box.height,
        2.0 * (leg_center_x - pie_box.x0),
        2.0 * (pie_box.x1 - leg_center_x),
    )
    pie_x0 = leg_center_x - max_side / 2
    pie_y0 = pie_box.y0 + (pie_box.height - max_side) / 2
    ax_pie = fig.add_axes([pie_x0, pie_y0, max_side, max_side])
    plot_ilp_failure_pie(
        ax_pie,
        n_parse_fail,
        n_ilp_fail,
        n_total,
        panel=True,
        show_legend=False,
        value_fontsize=VALUE_FONTSIZE,
        value_color=VALUE_COLOR,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pie-size-scale",
        type=float,
        default=PIE_SIZE_SCALE,
        metavar="S",
        help="scale factor for pie axes side relative to auto layout (default: %(default)g)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    pie_size_scale = float(args.pie_size_scale)
    if pie_size_scale <= 0:
        raise SystemExit("--pie-size-scale must be positive")
    skip_ids = load_skip_ids()
    keys, ilp_times = load_ilp_keys_and_times()
    if len(keys) != N_COMPLEXES:
        raise SystemExit(f"Expected {N_COMPLEXES} complexes, got {len(keys)}")
    if any(k[0] in skip_ids for k in keys):
        raise SystemExit("Skip-list IDs found in plotting set")

    counts = load_ilp_outcome_counts()
    n_total = counts["n_total"]
    n_parse_fail = counts["ilp_ok_smiles_not_parseable"]
    n_ilp_fail = counts["ilp_failed"]
    if n_parse_fail + n_ilp_fail == 0:
        raise SystemExit("No failure cases to plot")

    ilp_arr = np.asarray(ilp_times, dtype=float)
    _apply_nature_style()
    _apply_panel_value_typography()

    fig = plt.figure(figsize=(7.2, 2.8))
    gs = fig.add_gridspec(
        1,
        2,
        width_ratios=[1.35, 0.95],
        wspace=0.22,
        left=0.08,
        right=0.98,
        bottom=0.16,
        top=0.96,
    )
    ax_hist = fig.add_subplot(gs[0])
    gs_pie = gs[1].subgridspec(
        2,
        1,
        height_ratios=[PIE_LEGEND_HEIGHT_RATIO, 1.0 - PIE_LEGEND_HEIGHT_RATIO],
        hspace=0.03,
    )
    ax_legend = fig.add_subplot(gs_pie[0])
    ax_pie_area = fig.add_subplot(gs_pie[1])

    _plot_ilp_histogram(ax_hist, ilp_arr, panel_title=None)
    _add_mean_median_lines_ilp(ax_hist, ilp_arr, symbols_only=True, indexed_x=False)
    _add_panel_pie(
        fig,
        ax_legend,
        ax_pie_area,
        n_parse_fail,
        n_ilp_fail,
        n_total,
        pie_size_scale=pie_size_scale,
    )

    fig.savefig(OUT_PNG, dpi=FIG_DPI, pad_inches=0.04)
    plt.close(fig)

    summary = {
        "n_complexes": len(keys),
        "n_total_outcomes": n_total,
        "ilp_seconds": {
            "median": float(np.median(ilp_arr)),
            "mean": float(np.mean(ilp_arr)),
        },
        "pie_slices": {
            "smiles_not_parseable": n_parse_fail,
            "not_solved": n_ilp_fail,
        },
        "pie_size_scale": pie_size_scale,
        "pie_legend_height_ratio": PIE_LEGEND_HEIGHT_RATIO,
        "output_png": str(OUT_PNG),
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
