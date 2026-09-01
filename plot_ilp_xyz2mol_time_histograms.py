#!/usr/bin/env python3
"""Compare ILP vs xyz2mol wall-time distributions on the same 59,596 complexes."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

ROOT = Path(__file__).resolve().parent
ILP_CSV = ROOT / "tmqmg_ilp_benchmark_output" / "per_structure_results.csv"
XYZ2MOL_CSV = Path("/home/zhujingyuan/xyz2mol_tm/results/tmqmg_csd_xyz2mol_results.csv")
SKIP_PATHS = (
    ROOT / "list_error_geometry.txt",
    ROOT / "csd_codes_boron_ge6.txt",
)
OUT_PNG = ROOT / "ilp_xyz2mol_time_histograms.png"
OUT_PANEL_PNG = ROOT / "ilp_heuristic_time_histograms.png"
OUT_ILP_PNG = ROOT / "ilp_time_histogram.png"
OUT_XYZ2MOL_PNG = ROOT / "xyz2mol_time_histogram.png"
OUT_JSON = ROOT / "ilp_xyz2mol_time_histograms_summary.json"
FIG_DPI = 300
N_COMPLEXES = 59_596
FINE_BIN_WIDTH = 0.1
FINE_BIN_COUNT = 20
FINE_XMAX = FINE_BIN_WIDTH * FINE_BIN_COUNT
STAT_LINE_LW = 1.6
ILP_X_TICKS = (0.0, 0.5, 1.0, 1.5, 2.0)
XYZ2MOL_COARSE_EDGES = (5.0, 100.0, 300.0)
# Visual width of each coarse bar (2–5, 5–100, 100–300 s) in plot x-units.
# One unit equals the width of one 0.1 s fine bar; increase to widen coarse bars.
XYZ2MOL_COARSE_BAR_WIDTH = 1.0
XYZ2MOL_X_TICKS = (0.0, 0.5, 1.0, 1.5, 2.0, 5.0, 100.0, 300.0)
# Wavy x-axis baseline for non-uniform segments (2–5, 5–100, 100–300 s): "zigzag" or "wave".
XYZ2MOL_AXIS_BREAK_STYLE = "zigzag"

# Light cool-tone palette with semi-transparent bars.
BAR_ALPHA = 0.72
NATURE_FONT_SANS = ["Arial", "Helvetica", "DejaVu Sans"]
COLOR_ILP = "#7EB6D7"
COLOR_HEURISTIC = "#94C9D4"
COLOR_MEAN = "#C07878"
COLOR_MEDIAN = "#566B7A"
COLOR_AXIS = "#444444"


def _apply_nature_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": NATURE_FONT_SANS,
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 9,
            "axes.titleweight": "normal",
            "axes.linewidth": 0.6,
            "axes.edgecolor": COLOR_AXIS,
            "axes.labelcolor": COLOR_AXIS,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "xtick.color": COLOR_AXIS,
            "ytick.color": COLOR_AXIS,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "legend.fontsize": 7,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.dpi": FIG_DPI,
        }
    )


def _style_nature_axes(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(COLOR_AXIS)
    ax.spines["bottom"].set_color(COLOR_AXIS)
    ax.tick_params(colors=COLOR_AXIS, width=0.6, length=3)
    ax.grid(False)


def _add_panel_title(ax, title: str) -> None:
    ax.set_title(title, loc="center", pad=6)


def _set_xyz2mol_xticklabels(ax) -> None:
    labels = [f"{t:g}" for t in XYZ2MOL_X_TICKS]
    ax.set_xticklabels(labels)
    for tick_val, label in zip(XYZ2MOL_X_TICKS, ax.get_xticklabels()):
        if tick_val in (100.0, 300.0):
            label.set_rotation(40)
            label.set_ha("right")
            label.set_rotation_mode("anchor")


def _fine_bin_edges() -> np.ndarray:
    fine = np.arange(0.0, FINE_XMAX + FINE_BIN_WIDTH * 0.5, FINE_BIN_WIDTH)
    assert len(fine) == FINE_BIN_COUNT + 1
    assert fine[-1] == FINE_XMAX
    return fine


def _xyz2mol_histogram_bin_edges() -> np.ndarray:
    """0–2 s: 20 bins of 0.1 s; then (2,5], (5,100], (100,300]."""
    return np.concatenate([_fine_bin_edges(), np.array(XYZ2MOL_COARSE_EDGES)])


def _xyz2mol_coarse_layout(n_fine: int, coarse_bar_width: float) -> tuple[float, float]:
    coarse_start = float(n_fine)
    xmax = coarse_start + 3.0 * coarse_bar_width
    return coarse_start, xmax


def _time_to_hybrid_x(
    t: float,
    bin_edges: np.ndarray,
    *,
    n_fine: int,
    coarse_start: float,
    fine_xmax: float,
    coarse_bar_width: float,
) -> float:
    if t <= fine_xmax:
        idx = int(np.searchsorted(bin_edges, t, side="right") - 1)
        idx = int(np.clip(idx, 0, n_fine - 1))
        left = float(bin_edges[idx])
        right = float(bin_edges[idx + 1])
        frac = (t - left) / (right - left) if right > left else 0.5
        return idx + frac

    coarse_bounds = (fine_xmax, *XYZ2MOL_COARSE_EDGES)
    for i in range(3):
        left_t = float(coarse_bounds[i])
        right_t = float(coarse_bounds[i + 1])
        if t <= right_t or i == 2:
            frac = (min(t, right_t) - left_t) / (right_t - left_t) if right_t > left_t else 0.5
            return coarse_start + (i + frac) * coarse_bar_width
    return coarse_start + 3.0 * coarse_bar_width


def _segment_axis_path(
    x0: float,
    x1: float,
    *,
    style: str,
    y0: float,
    y_amp: float,
) -> tuple[np.ndarray, np.ndarray]:
    if x1 <= x0:
        return np.array([x0]), np.array([y0])
    n_teeth = max(2, min(5, int(round((x1 - x0) * 1.5)) + 2))
    if style == "wave":
        t = np.linspace(0.0, 1.0, max(24, int(30 * (x1 - x0)) + 12))
        xs = x0 + t * (x1 - x0)
        ys = y0 + y_amp * np.sin(2.0 * np.pi * n_teeth * t)
        return xs, ys
    xs = np.linspace(x0, x1, 2 * n_teeth + 1)
    ys = np.array([y0 + (i % 2) * y_amp for i in range(2 * n_teeth + 1)], dtype=float)
    return xs, ys


def _draw_xyz2mol_custom_xaxis(
    ax,
    tick_x: list[float],
    *,
    style: str = XYZ2MOL_AXIS_BREAK_STYLE,
    axis_color: str = COLOR_AXIS,
) -> None:
    """Straight baseline for 0–2 s; zigzag/wave baseline between 2–5, 5–100, 100–300 s ticks."""
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    y0 = 0.0
    y_amp = 0.03
    tick_len = 0.024
    color = axis_color
    lw = 1.0

    ax.spines["bottom"].set_visible(False)
    ax.tick_params(axis="x", length=0, pad=7)

    # Uniform segment: tick 0 … tick 2 s.
    ax.plot(
        [tick_x[0], tick_x[4]],
        [y0, y0],
        transform=trans,
        color=color,
        lw=lw,
        clip_on=False,
        zorder=10,
    )
    # Non-uniform segments: 2–5, 5–100, 100–300 s.
    for i in range(4, len(tick_x) - 1):
        xs, ys = _segment_axis_path(
            tick_x[i],
            tick_x[i + 1],
            style=style,
            y0=y0,
            y_amp=y_amp,
        )
        ax.plot(xs, ys, transform=trans, color=color, lw=lw, clip_on=False, zorder=10)

    for x in tick_x:
        ax.plot(
            [x, x],
            [y0 - tick_len, y0],
            transform=trans,
            color=color,
            lw=lw,
            clip_on=False,
            zorder=11,
        )


def _time_to_fine_index_x(t: float) -> float:
    if t <= 0.0:
        return 0.0
    if t >= FINE_XMAX:
        return float(FINE_BIN_COUNT)
    return t / FINE_BIN_WIDTH


def _draw_straight_xaxis(ax, x0: float, x1: float, *, axis_color: str = COLOR_AXIS) -> None:
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    y0 = 0.0
    tick_len = 0.024
    ax.plot([x0, x1], [y0, y0], transform=trans, color=axis_color, lw=1.0, clip_on=False, zorder=10)
    for x in ax.get_xticks():
        if x0 - 1e-9 <= x <= x1 + 1e-9:
            ax.plot(
                [x, x],
                [y0 - tick_len, y0],
                transform=trans,
                color=axis_color,
                lw=1.0,
                clip_on=False,
                zorder=11,
            )


def _align_panel_xlabels(axes, y: float = -0.21) -> None:
    for ax in axes:
        ax.set_xlabel("Time (s)")
        ax.xaxis.set_label_coords(0.5, y)


def _plot_ilp_histogram(
    ax,
    arr: np.ndarray,
    *,
    color: str = COLOR_ILP,
    show_ylabel: bool = True,
    show_xlabel: bool = True,
    panel_title: str | None = "ILP",
    indexed_x: bool = False,
) -> None:
    ilp_bins = _fine_bin_edges()
    counts, _ = np.histogram(arr, bins=ilp_bins)

    if indexed_x:
        ax.bar(
            np.arange(FINE_BIN_COUNT, dtype=float),
            counts,
            width=1.0,
            align="edge",
            color=color,
            edgecolor="white",
            linewidth=0.35,
            alpha=BAR_ALPHA,
        )
        ax.set_xlim(0.0, float(FINE_BIN_COUNT))
        tick_x = [_time_to_fine_index_x(t) for t in ILP_X_TICKS]
        ax.set_xticks(tick_x)
        ax.set_xticklabels([f"{t:g}" for t in ILP_X_TICKS])
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(axis="x", length=0, pad=7)
        _draw_straight_xaxis(ax, 0.0, float(FINE_BIN_COUNT))
    else:
        ax.hist(
            arr,
            bins=ilp_bins,
            color=color,
            edgecolor="white",
            linewidth=0.35,
            alpha=BAR_ALPHA,
        )
        ax.set_xlim(0.0, FINE_XMAX)
        ax.set_xticks(ILP_X_TICKS)
        ax.set_xticklabels([f"{t:g}" for t in ILP_X_TICKS])

    if show_xlabel:
        ax.set_xlabel("Time (s)")
    if show_ylabel:
        ax.set_ylabel("Logarithmic count")
    ax.set_yscale("log")
    if panel_title:
        _add_panel_title(ax, panel_title)
    _style_nature_axes(ax)


def _add_mean_median_lines_ilp(
    ax,
    arr: np.ndarray,
    *,
    symbols_only: bool = False,
    indexed_x: bool = False,
) -> None:
    mu = float(np.mean(arr))
    eta = float(np.median(arr))
    if indexed_x:
        x_mu = _time_to_fine_index_x(mu)
        x_eta = _time_to_fine_index_x(eta)
        ax.axvline(x_mu, color=COLOR_MEAN, ls="-.", lw=STAT_LINE_LW)
        ax.axvline(x_eta, color=COLOR_MEDIAN, ls="--", lw=STAT_LINE_LW)
    else:
        ax.axvline(mu, color=COLOR_MEAN, ls="-.", lw=STAT_LINE_LW)
        ax.axvline(eta, color=COLOR_MEDIAN, ls="--", lw=STAT_LINE_LW)
    if symbols_only:
        handles = [
            Line2D([0], [0], color=COLOR_MEAN, ls="-.", lw=STAT_LINE_LW, label=rf"$\mu$ = {mu:.3f} s"),
            Line2D([0], [0], color=COLOR_MEDIAN, ls="--", lw=STAT_LINE_LW, label=rf"$\eta$ = {eta:.3f} s"),
        ]
    else:
        handles = [
            Line2D([0], [0], color=COLOR_MEAN, ls="-.", lw=STAT_LINE_LW, label=rf"mean $\mu$ = {mu:.3f} s"),
            Line2D([0], [0], color=COLOR_MEDIAN, ls="--", lw=STAT_LINE_LW, label=rf"median $\eta$ = {eta:.3f} s"),
        ]
    ax.legend(handles=handles, frameon=False, loc="upper right")


def _plot_xyz2mol_histogram(
    ax,
    arr: np.ndarray,
    bin_edges: np.ndarray,
    *,
    coarse_bar_width: float,
    axis_break_style: str = XYZ2MOL_AXIS_BREAK_STYLE,
    color: str = COLOR_HEURISTIC,
    show_ylabel: bool = False,
    show_xlabel: bool = True,
    panel_title: str | None = "Heuristic",
) -> None:
    counts, _ = np.histogram(arr, bins=bin_edges)
    n_fine = FINE_BIN_COUNT
    coarse_start, xmax = _xyz2mol_coarse_layout(n_fine, coarse_bar_width)

    ax.bar(
        np.arange(n_fine, dtype=float),
        counts[:n_fine],
        width=1.0,
        align="edge",
        color=color,
        edgecolor="white",
        linewidth=0.35,
        alpha=BAR_ALPHA,
    )
    for i, count in enumerate(counts[n_fine : n_fine + 3]):
        left = coarse_start + i * coarse_bar_width
        ax.bar(
            left,
            count,
            width=coarse_bar_width,
            align="edge",
            color=color,
            edgecolor="white",
            linewidth=0.35,
            alpha=BAR_ALPHA,
        )

    ax.set_xlim(0.0, xmax)
    if show_xlabel:
        ax.set_xlabel("Time (s)")
    if show_ylabel:
        ax.set_ylabel("Logarithmic count")
    ax.set_yscale("log")

    tick_x = [
        _time_to_hybrid_x(
            t,
            bin_edges,
            n_fine=n_fine,
            coarse_start=coarse_start,
            fine_xmax=FINE_XMAX,
            coarse_bar_width=coarse_bar_width,
        )
        for t in XYZ2MOL_X_TICKS
    ]
    ax.set_xticks(tick_x)
    _set_xyz2mol_xticklabels(ax)
    _draw_xyz2mol_custom_xaxis(ax, tick_x, style=axis_break_style)
    if panel_title:
        _add_panel_title(ax, panel_title)
    _style_nature_axes(ax)


def _add_mean_median_lines_hybrid(
    ax,
    arr: np.ndarray,
    bin_edges: np.ndarray,
    *,
    coarse_bar_width: float,
    symbols_only: bool = False,
) -> None:
    n_fine = FINE_BIN_COUNT
    coarse_start, _ = _xyz2mol_coarse_layout(n_fine, coarse_bar_width)
    mu = float(np.mean(arr))
    eta = float(np.median(arr))
    x_mu = _time_to_hybrid_x(
        mu,
        bin_edges,
        n_fine=n_fine,
        coarse_start=coarse_start,
        fine_xmax=FINE_XMAX,
        coarse_bar_width=coarse_bar_width,
    )
    x_eta = _time_to_hybrid_x(
        eta,
        bin_edges,
        n_fine=n_fine,
        coarse_start=coarse_start,
        fine_xmax=FINE_XMAX,
        coarse_bar_width=coarse_bar_width,
    )
    ax.axvline(x_mu, color=COLOR_MEAN, ls="-.", lw=STAT_LINE_LW)
    ax.axvline(x_eta, color=COLOR_MEDIAN, ls="--", lw=STAT_LINE_LW)
    if symbols_only:
        handles = [
            Line2D([0], [0], color=COLOR_MEAN, ls="-.", lw=STAT_LINE_LW, label=rf"$\mu$ = {mu:.3f} s"),
            Line2D([0], [0], color=COLOR_MEDIAN, ls="--", lw=STAT_LINE_LW, label=rf"$\eta$ = {eta:.3f} s"),
        ]
    else:
        handles = [
            Line2D([0], [0], color=COLOR_MEAN, ls="-.", lw=STAT_LINE_LW, label=rf"mean $\mu$ = {mu:.3f} s"),
            Line2D([0], [0], color=COLOR_MEDIAN, ls="--", lw=STAT_LINE_LW, label=rf"median $\eta$ = {eta:.3f} s"),
        ]
    ax.legend(handles=handles, frameon=False, loc="upper right")


def _add_mean_median_lines(ax, arr: np.ndarray, *, symbols_only: bool = False) -> None:
    mu = float(np.mean(arr))
    eta = float(np.median(arr))
    ax.axvline(mu, color=COLOR_MEAN, ls="-.", lw=STAT_LINE_LW)
    ax.axvline(eta, color=COLOR_MEDIAN, ls="--", lw=STAT_LINE_LW)
    if symbols_only:
        handles = [
            Line2D([0], [0], color=COLOR_MEAN, ls="-.", lw=STAT_LINE_LW, label=rf"$\mu$ = {mu:.3f} s"),
            Line2D([0], [0], color=COLOR_MEDIAN, ls="--", lw=STAT_LINE_LW, label=rf"$\eta$ = {eta:.3f} s"),
        ]
    else:
        handles = [
            Line2D([0], [0], color=COLOR_MEAN, ls="-.", lw=STAT_LINE_LW, label=rf"mean $\mu$ = {mu:.3f} s"),
            Line2D([0], [0], color=COLOR_MEDIAN, ls="--", lw=STAT_LINE_LW, label=rf"median $\eta$ = {eta:.3f} s"),
        ]
    ax.legend(handles=handles, frameon=False, loc="upper right")


def _plot_combined_panel_figure(
    ilp_arr: np.ndarray,
    xyz_arr: np.ndarray,
    *,
    coarse_bar_width: float,
    axis_break_style: str,
) -> plt.Figure:
    xyz_bins = _xyz2mol_histogram_bin_edges()
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(8.2, 2.6),
        sharey=True,
        gridspec_kw={"width_ratios": [1, 1], "wspace": 0.03},
    )
    ax_ilp, ax_heur = axes

    _plot_ilp_histogram(
        ax_ilp,
        ilp_arr,
        show_ylabel=True,
        show_xlabel=False,
        panel_title="ILP",
        indexed_x=True,
    )
    _add_mean_median_lines_ilp(ax_ilp, ilp_arr, symbols_only=True, indexed_x=True)

    _plot_xyz2mol_histogram(
        ax_heur,
        xyz_arr,
        xyz_bins,
        coarse_bar_width=coarse_bar_width,
        axis_break_style=axis_break_style,
        show_ylabel=False,
        show_xlabel=False,
        panel_title="Heuristic",
    )
    _add_mean_median_lines_hybrid(
        ax_heur, xyz_arr, xyz_bins, coarse_bar_width=coarse_bar_width, symbols_only=True
    )

    fig.subplots_adjust(left=0.06, right=0.998, bottom=0.22, top=0.86, wspace=0.03)
    _align_panel_xlabels(axes)
    return fig


def norm_charge(raw: str | None) -> int:
    if raw is None or str(raw).strip() == "":
        return 0
    return int(round(float(raw)))


def load_skip_ids() -> frozenset[str]:
    out: set[str] = set()
    for path in SKIP_PATHS:
        if not path.is_file():
            continue
        for line in path.read_text(encoding="utf-8").splitlines():
            code = line.split("#", 1)[0].strip()
            if code:
                out.add(code)
    return frozenset(out)


def load_ilp_keys_and_times() -> tuple[list[tuple[str, int]], list[float]]:
    rows: list[tuple[str, int, float]] = []
    for row in csv.DictReader(ILP_CSV.open(encoding="utf-8")):
        if row["status"] == "missing_xyz":
            continue
        key = (row["IDs"], norm_charge(row["charge"]))
        rows.append((key[0], key[1], float(row["ilp_seconds"])))
    keys = [(a, b) for a, b, _ in rows]
    times = [t for _, _, t in rows]
    return keys, times


def load_xyz2mol_times(keys: list[tuple[str, int]]) -> list[float]:
    by_key: dict[tuple[str, int], float] = {}
    for row in csv.DictReader(XYZ2MOL_CSV.open(encoding="utf-8")):
        key = (row["ID"], norm_charge(row["charge"]))
        raw = str(row.get("elapsed_s", "")).strip()
        by_key[key] = float(raw) if raw else 0.0
    missing = [k for k in keys if k not in by_key]
    if missing:
        raise RuntimeError(f"xyz2mol missing {len(missing)} keys, e.g. {missing[:5]}")
    return [by_key[k] for k in keys]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--coarse-bar-width",
        type=float,
        default=XYZ2MOL_COARSE_BAR_WIDTH,
        metavar="W",
        help=(
            "visual width of each coarse bar (2–5, 5–100, 100–300 s) in plot x-units; "
            "one unit = one 0.1 s fine bar (default: %(default)g)"
        ),
    )
    parser.add_argument(
        "--axis-break-style",
        choices=("zigzag", "wave"),
        default=XYZ2MOL_AXIS_BREAK_STYLE,
        help="wavy x-axis baseline style for 2–5, 5–100, 100–300 s segments (default: %(default)s)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    coarse_bar_width = float(args.coarse_bar_width)
    axis_break_style = str(args.axis_break_style)
    if coarse_bar_width <= 0:
        raise SystemExit("--coarse-bar-width must be positive")

    if not ILP_CSV.is_file():
        raise SystemExit(f"ILP CSV not found: {ILP_CSV}")
    if not XYZ2MOL_CSV.is_file():
        raise SystemExit(f"xyz2mol CSV not found: {XYZ2MOL_CSV}")

    skip_ids = load_skip_ids()
    keys, ilp_times = load_ilp_keys_and_times()
    xyz_times = load_xyz2mol_times(keys)

    if len(keys) != N_COMPLEXES:
        raise SystemExit(f"Expected {N_COMPLEXES} complexes, got {len(keys)}")
    if any(k[0] in skip_ids for k in keys):
        raise SystemExit("Skip-list IDs found in plotting set")

    ilp_arr = np.asarray(ilp_times, dtype=float)
    xyz_arr = np.asarray(xyz_times, dtype=float)
    _apply_nature_style()

    n_gt_5_ilp = int((ilp_arr > 5.0).sum())
    n_gt_5_xyz = int((xyz_arr > 5.0).sum())
    xmax = float(np.ceil(max(ilp_arr.max(), xyz_arr.max())))

    fig_panel = _plot_combined_panel_figure(
        ilp_arr,
        xyz_arr,
        coarse_bar_width=coarse_bar_width,
        axis_break_style=axis_break_style,
    )
    fig_panel.savefig(OUT_PANEL_PNG, dpi=FIG_DPI, bbox_inches="tight", pad_inches=0.04)
    fig_panel.savefig(OUT_PNG, dpi=FIG_DPI, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig_panel)

    ilp_bins = _fine_bin_edges()
    fig_ilp, ax_ilp = plt.subplots(figsize=(3.5, 2.55))
    _plot_ilp_histogram(ax_ilp, ilp_arr, panel_title=None)
    _add_mean_median_lines_ilp(ax_ilp, ilp_arr, symbols_only=True, indexed_x=False)
    fig_ilp.tight_layout()
    fig_ilp.savefig(OUT_ILP_PNG, dpi=FIG_DPI, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig_ilp)

    xyz_bins = _xyz2mol_histogram_bin_edges()
    fig_xyz, ax_xyz = plt.subplots(figsize=(7.08, 2.55))
    _plot_xyz2mol_histogram(
        ax_xyz,
        xyz_arr,
        xyz_bins,
        coarse_bar_width=coarse_bar_width,
        axis_break_style=axis_break_style,
        show_ylabel=True,
        panel_title=None,
    )
    _add_mean_median_lines_hybrid(
        ax_xyz, xyz_arr, xyz_bins, coarse_bar_width=coarse_bar_width, symbols_only=True
    )
    fig_xyz.subplots_adjust(bottom=0.22)
    fig_xyz.tight_layout()
    fig_xyz.savefig(OUT_XYZ2MOL_PNG, dpi=FIG_DPI, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig_xyz)

    summary = {
        "n_complexes": len(keys),
        "skip_list_size": len(skip_ids),
        "ilp_seconds": {
            "median": float(np.median(ilp_arr)),
            "mean": float(np.mean(ilp_arr)),
            "p95": float(np.percentile(ilp_arr, 95)),
            "p99": float(np.percentile(ilp_arr, 99)),
            "max": float(np.max(ilp_arr)),
        },
        "xyz2mol_elapsed_s": {
            "median": float(np.median(xyz_arr)),
            "mean": float(np.mean(xyz_arr)),
            "p95": float(np.percentile(xyz_arr, 95)),
            "p99": float(np.percentile(xyz_arr, 99)),
            "max": float(np.max(xyz_arr)),
            "n_ge_300s": int((xyz_arr >= 300).sum()),
            "n_gt_5s": n_gt_5_xyz,
        },
        "ilp_n_gt_5s": n_gt_5_ilp,
        "histogram_xmax_s": xmax,
        "output_png": str(OUT_PNG),
        "output_panel_png": str(OUT_PANEL_PNG),
        "output_ilp_png": str(OUT_ILP_PNG),
        "output_xyz2mol_png": str(OUT_XYZ2MOL_PNG),
        "xyz2mol_histogram_bins": (
            f"{FINE_BIN_COUNT} fine bins of {FINE_BIN_WIDTH:g} s in "
            f"0–{FINE_XMAX:g} s + 3 wider coarse bins (2–5, 5–100, 100–300 s)"
        ),
        "ilp_histogram_bins": f"{FINE_BIN_COUNT} bins of {FINE_BIN_WIDTH:g} s in 0–{FINE_XMAX:g} s",
        "xyz2mol_coarse_bar_width": coarse_bar_width,
        "xyz2mol_axis_break_style": axis_break_style,
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
