#!/usr/bin/env python3
"""Group-6 MLX count tables (L rows × X columns) from heatmap summary JSON."""

from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parent
DEFAULT_SUMMARY_JSON = {
    "Cr": ROOT / "cr_mlx_heatmap_summary.json",
    "Mo": ROOT / "mo_mlx_heatmap_summary.json",
    "W": ROOT / "w_mlx_heatmap_summary.json",
}
METAL_ORDER = ("Cr", "Mo", "W")
OUT_PNG = ROOT / "group6_mlx_table.png"
FIG_DPI = 300

L_MAX = 6
X_MAX = 6
L_VALUES = list(range(0, L_MAX + 1))  # 0 → 6, top to bottom
X_VALUES = list(range(0, X_MAX + 1))

NATURE_FONT_SANS = ["Arial", "Helvetica", "DejaVu Sans"]
COLOR_AXIS = "#444444"
COLOR_GRID = "#666666"
COLOR_HEADER_FILL = "#f4f4f4"
NATURE_HEATMAP_CMAP = LinearSegmentedColormap.from_list(
    "mlx_red",
    ["#ffffff", "#fee5d9", "#fcae91", "#fb6a4a", "#de2d26", "#a50f15"],
    N=256,
)
SUB_TO_ASCII = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")
LX_RE = re.compile(r"(?:L(\d+))?(?:X(\d+))?$")


def _apply_nature_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": NATURE_FONT_SANS,
            "font.size": 8,
            "axes.linewidth": 0.6,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.dpi": FIG_DPI,
        }
    )


def parse_lx_label(label: str) -> tuple[int, int]:
    text = (label or "").translate(SUB_TO_ASCII).strip()
    match = LX_RE.fullmatch(text)
    if not match or not text:
        raise ValueError(f"Cannot parse MLX label: {label!r}")
    l = int(match.group(1)) if match.group(1) else 0
    x = int(match.group(2)) if match.group(2) else 0
    return l, x


def counts_by_lx(summary: dict) -> Counter:
    labels = summary.get("cell_labels", {})
    out: Counter = Counter()
    for key, count in summary["cell_counts"].items():
        label = labels.get(key, "")
        l, x = parse_lx_label(label)
        out[(l, x)] += int(count)
    return out


def load_metal_counts() -> dict[str, Counter]:
    metal_counts = {}
    for metal in METAL_ORDER:
        path = DEFAULT_SUMMARY_JSON[metal]
        summary = json.loads(path.read_text(encoding="utf-8"))
        metal_counts[metal] = counts_by_lx(summary)
    return metal_counts


def _draw_split_corner(ax, x0: float, y0: float, w: float, h: float, metal: str) -> None:
    """Header cell: two rays from bottom-right → n (rows), No. of MLₙXₘ, m (columns)."""
    ax.add_patch(
        Rectangle(
            (x0, y0),
            w,
            h,
            facecolor="#ffffff",
            edgecolor=COLOR_GRID,
            lw=0.7,
            zorder=3,
        )
    )
    # y increases downward: y0 is the visual top. Both lines start at bottom-right.
    bottom_right = (x0 + w, y0 + h)
    left_hit = (x0, y0 + 0.58 * h)
    top_hit = (x0 + 0.58 * w, y0)
    ax.plot(
        [bottom_right[0], left_hit[0]],
        [bottom_right[1], left_hit[1]],
        color=COLOR_GRID,
        lw=0.7,
        zorder=4,
    )
    ax.plot(
        [bottom_right[0], top_hit[0]],
        [bottom_right[1], top_hit[1]],
        color=COLOR_GRID,
        lw=0.7,
        zorder=4,
    )

    ax.text(
        x0 + 0.22 * w,
        y0 + 0.84 * h,
        r"$n$",
        ha="center",
        va="center",
        fontsize=9,
        color="#111111",
        zorder=5,
    )
    ax.text(
        x0 + 0.86 * w,
        y0 + 0.20 * h,
        r"$m$",
        ha="center",
        va="center",
        fontsize=9,
        color="#111111",
        zorder=5,
    )
    ax.text(
        x0 + 0.36 * w,
        y0 + 0.34 * h,
        f"No. of\n$\\mathrm{{{metal}L}}_n\\mathrm{{X}}_m$",
        ha="center",
        va="center",
        fontsize=6.2,
        color="#111111",
        linespacing=1.05,
        zorder=5,
    )


def _plot_one_table(ax, cell_counts: Counter, metal: str) -> None:
    n_row = len(L_VALUES)
    n_col = len(X_VALUES)
    hw, hh = 1.85, 1.40

    grid = np.zeros((n_row, n_col), dtype=float)
    for i, l in enumerate(L_VALUES):
        for j, x in enumerate(X_VALUES):
            grid[i, j] = float(cell_counts.get((l, x), 0))

    ax.set_xlim(0, hw + n_col)
    ax.set_ylim(0, hh + n_row)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.axis("off")

    for i, l in enumerate(L_VALUES):
        for j, x in enumerate(X_VALUES):
            count = int(grid[i, j])
            ax.add_patch(
                Rectangle(
                    (hw + j, hh + i),
                    1,
                    1,
                    facecolor="#ffffff",
                    edgecolor=COLOR_GRID,
                    lw=0.55,
                    zorder=2,
                )
            )
            if count:
                ax.text(
                    hw + j + 0.5,
                    hh + i + 0.5,
                    str(count),
                    ha="center",
                    va="center",
                    fontsize=10,
                    color="#111111",
                    zorder=5,
                )

    for j, x in enumerate(X_VALUES):
        ax.add_patch(
            Rectangle(
                (hw + j, 0),
                1,
                hh,
                facecolor=COLOR_HEADER_FILL,
                edgecolor=COLOR_GRID,
                lw=0.7,
                zorder=3,
            )
        )
        ax.text(
            hw + j + 0.5,
            hh * 0.5,
            str(x),
            ha="center",
            va="center",
            fontsize=10,
            color=COLOR_AXIS,
            zorder=5,
        )

    for i, l in enumerate(L_VALUES):
        ax.add_patch(
            Rectangle(
                (0, hh + i),
                hw,
                1,
                facecolor=COLOR_HEADER_FILL,
                edgecolor=COLOR_GRID,
                lw=0.7,
                zorder=3,
            )
        )
        ax.text(
            hw * 0.5,
            hh + i + 0.5,
            str(l),
            ha="center",
            va="center",
            fontsize=10,
            color=COLOR_AXIS,
            zorder=5,
        )

    _draw_split_corner(ax, 0.0, 0.0, hw, hh, metal)


def plot_tables(metal_counts: dict[str, Counter], out_png: Path) -> None:
    _apply_nature_style()
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.6))
    fig.subplots_adjust(wspace=0.08, left=0.02, right=0.99, top=0.97, bottom=0.03)
    for ax, metal in zip(axes, METAL_ORDER):
        _plot_one_table(ax, metal_counts[metal], metal)
    fig.savefig(out_png, dpi=FIG_DPI, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


def main() -> int:
    metal_counts = load_metal_counts()
    plot_tables(metal_counts, OUT_PNG)
    print(f"Wrote {OUT_PNG}")
    for metal in METAL_ORDER:
        total = sum(
            c
            for (l, x), c in metal_counts[metal].items()
            if 0 <= l <= L_MAX and 0 <= x <= X_MAX
        )
        outside = sum(
            c
            for (l, x), c in metal_counts[metal].items()
            if not (0 <= l <= L_MAX and 0 <= x <= X_MAX)
        )
        print(f"  {metal}: in-grid {total}, L/X outside 0–6: {outside}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
