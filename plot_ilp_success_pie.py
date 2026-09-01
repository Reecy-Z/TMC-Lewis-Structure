#!/usr/bin/env python3
"""Pie chart of ILP failure modes on tmQM-G (excl. skip lists & missing_xyz)."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
ILP_CSV = ROOT / "tmqmg_ilp_benchmark_output" / "per_structure_results.csv"
SKIP_PATHS = (
    ROOT / "list_error_geometry.txt",
    ROOT / "csd_codes_boron_ge6.txt",
)
OUT_PNG = ROOT / "ilp_success_pie.png"
OUT_JSON = ROOT / "ilp_success_pie_summary.json"
FIG_DPI = 300

COLOR_PARSE_FAIL = "#F2D0A9"
COLOR_ILP_FAIL = "#A8BEC9"
VALUE_FONTSIZE = 8
VALUE_COLOR = "#444444"
PIE_LEGEND_LABELS = (
    "Solved ILP with unparsable SMILES",
    "Unsolved ILP instances",
)
FONT_SANS = ["Arial", "Helvetica", "DejaVu Sans"]


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


def load_ilp_outcome_counts() -> dict[str, int]:
    skip_ids = load_skip_ids()
    n_total = 0
    n_ilp_ok_parse_ok = 0
    n_ilp_ok_parse_fail = 0
    n_ilp_failed = 0

    for row in csv.DictReader(ILP_CSV.open(encoding="utf-8")):
        if row["IDs"] in skip_ids:
            continue
        if row["status"] == "missing_xyz":
            continue
        n_total += 1
        status = row["status"]
        parse_ok = str(row.get("rdkit_parse_ok", "")).strip().lower() in ("1", "true")
        if status == "ok":
            if parse_ok:
                n_ilp_ok_parse_ok += 1
            else:
                n_ilp_ok_parse_fail += 1
        elif status == "ilp_failed":
            n_ilp_failed += 1
        else:
            raise RuntimeError(f"Unexpected status: {status!r} for {row['IDs']}")

    return {
        "n_total": n_total,
        "ilp_ok_smiles_parseable": n_ilp_ok_parse_ok,
        "ilp_ok_smiles_not_parseable": n_ilp_ok_parse_fail,
        "ilp_failed": n_ilp_failed,
    }


def _pct(n: int, total: int) -> float:
    return 100.0 * n / total if total else 0.0


def plot_ilp_failure_pie(
    ax,
    n_parse_fail: int,
    n_ilp_fail: int,
    n_total: int,
    *,
    legend_loc: str = "upper right",
    legend_bbox: tuple[float, float] | None = None,
    panel: bool = False,
    show_legend: bool = True,
    value_fontsize: float = VALUE_FONTSIZE,
    value_color: str = VALUE_COLOR,
) -> list:
    slices = [n_parse_fail, n_ilp_fail]
    colors = [COLOR_PARSE_FAIL, COLOR_ILP_FAIL]
    n_failures = sum(slices)
    if n_failures == 0:
        raise ValueError("No failure cases to plot")

    def _slice_autopct(pct: float) -> str:
        count = int(round(pct * n_failures / 100.0))
        frac_total = 100.0 * count / n_total
        return f"{count} ({frac_total:.2f}%)"

    pie_kwargs: dict = {
        "colors": colors,
        "startangle": 90,
        "counterclock": False,
        "autopct": _slice_autopct,
        "pctdistance": 0.55,
        "wedgeprops": {"edgecolor": "white", "linewidth": 2.0, "alpha": 0.88},
        "textprops": {"fontsize": value_fontsize, "color": value_color},
    }
    if panel:
        pie_kwargs["radius"] = 0.95

    wedges, _, autotexts = ax.pie(slices, **pie_kwargs)
    for t in autotexts:
        t.set_fontsize(value_fontsize)
        t.set_color(value_color)

    if panel:
        ax.set_xlim(-1.05, 1.05)
        ax.set_ylim(-1.05, 1.05)
        ax.set_aspect("equal")
        ax.axis("off")
    else:
        ax.set_aspect("equal")

    if show_legend:
        legend_kwargs: dict = {"frameon": False, "fontsize": value_fontsize}
        if legend_bbox is not None:
            ax.legend(
                wedges,
                PIE_LEGEND_LABELS,
                loc=legend_loc,
                bbox_to_anchor=legend_bbox,
                **legend_kwargs,
            )
        else:
            ax.legend(wedges, PIE_LEGEND_LABELS, loc=legend_loc, **legend_kwargs)
    return wedges


def main() -> int:
    if not ILP_CSV.is_file():
        raise SystemExit(f"ILP CSV not found: {ILP_CSV}")

    counts = load_ilp_outcome_counts()
    n_total = counts["n_total"]
    n_parse_ok = counts["ilp_ok_smiles_parseable"]
    n_parse_fail = counts["ilp_ok_smiles_not_parseable"]
    n_ilp_fail = counts["ilp_failed"]
    if n_parse_fail + n_ilp_fail == 0:
        raise SystemExit("No failure cases to plot")

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": FONT_SANS,
            "font.size": VALUE_FONTSIZE,
            "legend.fontsize": VALUE_FONTSIZE,
            "savefig.dpi": FIG_DPI,
            "figure.facecolor": "white",
        }
    )

    fig, ax = plt.subplots(figsize=(4.2, 3.2))
    plot_ilp_failure_pie(ax, n_parse_fail, n_ilp_fail, n_total)

    fig.savefig(OUT_PNG, dpi=FIG_DPI, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)

    n_ok = n_parse_ok + n_parse_fail
    n_failures = n_parse_fail + n_ilp_fail
    summary = {
        **counts,
        "pie_slices": {"smiles_not_parseable": n_parse_fail, "not_solved": n_ilp_fail},
        "ilp_solved": n_ok,
        "ilp_solved_fraction": _pct(n_ok, n_total) / 100.0,
        "ilp_failed_fraction": _pct(n_ilp_fail, n_total) / 100.0,
        "smiles_parseable_fraction_of_solved": _pct(n_parse_ok, n_ok) / 100.0 if n_ok else 0.0,
        "not_parseable_fraction_of_failures": _pct(n_parse_fail, n_failures) / 100.0,
        "not_solved_fraction_of_failures": _pct(n_ilp_fail, n_failures) / 100.0,
        "output_png": str(OUT_PNG),
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
