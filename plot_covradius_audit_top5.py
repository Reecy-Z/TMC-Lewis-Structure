#!/usr/bin/env python3
"""Top-5 M–L pair bar charts across alternative P99 distance margins.

Top-5 pairs are ranked from tmQM_covradius_audit.pkl (CCDC R_cov + 0.45 Å).
Each pair compares R_cov+0.45 Å with P99, P99+0.05 Å, and P99+0.15 Å.
"""

from __future__ import annotations

import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
AUDIT_COV_PKL = ROOT / "tmQM_covradius_audit.pkl"
AUDIT_P99_PKLS = (
    ("P99", ROOT / "tmQM_covradius_audit_cutoff_P99.pkl"),
    ("P99 + 0.05 Å", ROOT / "tmQM_covradius_audit_cutoff_P99_0p05.pkl"),
    ("P99 + 0.15 Å", ROOT / "tmQM_covradius_audit_cutoff.pkl"),
)
MIN_REF_BONDS = 100
MIN_REF_NONBONDS = 100
FIG_DPI = 300

LABEL_COV = r"$R_{\mathrm{cov}}(M)+R_{\mathrm{cov}}(L)+0.45$ Å"
FN_COLORS = ("#8E2C2C", "#C44E52", "#E98B79", "#F6C1A7")
FP_COLORS = ("#2F5597", "#4C72B0", "#79A6D2", "#B9D4EA")


def _rows_from_bep(bep: dict) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for pair, rec in bep.items():
        tp, fp, fn, tn = rec["tp"], rec["fp"], rec["fn"], rec["tn"]
        n_bond = tp + fn
        n_nonbond = fp + tn
        out[pair] = {
            "fn": fn,
            "fp": fp,
            "n_bond": n_bond,
            "n_nonbond": n_nonbond,
            "miss_rate": fn / n_bond if n_bond else 0.0,
            "fp_rate": fp / n_nonbond if n_nonbond else 0.0,
        }
    return out


def _load_audits() -> tuple[dict[str, dict], list[tuple[str, dict[str, dict]]]]:
    with open(AUDIT_COV_PKL, "rb") as fh:
        cov_bep = pickle.load(fh)["by_element_pair"]
    alternatives = []
    for label, path in AUDIT_P99_PKLS:
        with open(path, "rb") as fh:
            alternatives.append((label, _rows_from_bep(pickle.load(fh)["by_element_pair"])))
    return _rows_from_bep(cov_bep), alternatives


def top5_fn_pairs(cov_rows: dict[str, dict]) -> list[str]:
    eligible = [
        p
        for p, r in cov_rows.items()
        if r["n_bond"] >= MIN_REF_BONDS and r["miss_rate"] > 0
    ]
    return sorted(eligible, key=lambda p: cov_rows[p]["miss_rate"], reverse=True)[:5]


def top5_fp_pairs(cov_rows: dict[str, dict]) -> list[str]:
    eligible = [
        p
        for p, r in cov_rows.items()
        if r["n_bond"] >= MIN_REF_BONDS
        and r["n_nonbond"] >= MIN_REF_NONBONDS
        and r["fp_rate"] > 0
    ]
    return sorted(eligible, key=lambda p: cov_rows[p]["fp_rate"], reverse=True)[:5]


def plot_compare_top5(
    pairs: list[str],
    cov_rows: dict[str, dict],
    alternatives: list[tuple[str, dict[str, dict]]],
    *,
    rate_key: str,
    num_key: str,
    den_key: str,
    title: str,
    ylabel: str,
    out_path: Path,
    colors: tuple[str, ...],
) -> None:
    n = len(pairs)
    series = [(LABEL_COV, cov_rows), *alternatives]
    all_rates = [
        [rows[p][rate_key] * 100 for p in pairs]
        for _label, rows in series
    ]
    ymax = max(max(rates) for rates in all_rates) * 1.42

    fig, ax = plt.subplots(figsize=(11.5, 6.0))
    x = np.arange(n)
    width = 0.19
    offsets = (np.arange(len(series)) - (len(series) - 1) / 2) * width
    bar_groups = []
    for offset, (label, _rows), rates, color in zip(offsets, series, all_rates, colors):
        bars = ax.bar(
            x + offset,
            rates,
            width,
            label=label,
            color=color,
            edgecolor="#333333",
            linewidth=0.7,
        )
        bar_groups.append(bars)

    ax.set_ylim(0, ymax)
    ax.set_xticks(x)
    ax.set_xticklabels(pairs, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=13, pad=12)
    ax.yaxis.grid(True, linestyle="--", alpha=0.35)
    ax.set_axisbelow(True)
    ax.legend(loc="upper right", fontsize=9, framealpha=0.92, ncol=2)

    for bars, (_label, rows) in zip(bar_groups, series):
        for bar, pair in zip(bars, pairs):
            row = rows[pair]
            pct = row[rate_key] * 100
            num, den = row[num_key], row[den_key]
            h = bar.get_height()
            label = "0%" if num == 0 else f"{pct:.2f}%\n({num:,}/{den:,})"
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                h + ymax * 0.015 if h > 0 else ymax * 0.01,
                label,
                ha="center",
                va="bottom",
                fontsize=7,
                fontweight="medium",
            )

    fig.tight_layout()
    fig.savefig(out_path, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    cov_rows, alternatives = _load_audits()
    top_fn = top5_fn_pairs(cov_rows)
    top_fp = top5_fp_pairs(cov_rows)

    plot_compare_top5(
        top_fn,
        cov_rows,
        alternatives,
        rate_key="miss_rate",
        num_key="fn",
        den_key="n_bond",
        title="M–L pairs with the highest missed bonds",
        ylabel="Missed bond rate (%)",
        out_path=ROOT / "audit_top5_FN.png",
        colors=FN_COLORS,
    )
    plot_compare_top5(
        top_fp,
        cov_rows,
        alternatives,
        rate_key="fp_rate",
        num_key="fp",
        den_key="n_nonbond",
        title="M–L pairs with the highest spurious bonds",
        ylabel="Spurious bond rate (%)",
        out_path=ROOT / "audit_top5_FP.png",
        colors=FP_COLORS,
    )

    print("Wrote", ROOT / "audit_top5_FN.png")
    for p in top_fn:
        values = [(LABEL_COV, cov_rows[p]), *[(label, rows[p]) for label, rows in alternatives]]
        print("  " + p + ": " + " | ".join(
            f"{label} {row['miss_rate']*100:.2f}% ({row['fn']}/{row['n_bond']})"
            for label, row in values
        ))
    print("Wrote", ROOT / "audit_top5_FP.png")
    for p in top_fp:
        values = [(LABEL_COV, cov_rows[p]), *[(label, rows[p]) for label, rows in alternatives]]
        print("  " + p + ": " + " | ".join(
            f"{label} {row['fp_rate']*100:.2f}% ({row['fp']}/{row['n_nonbond']})"
            for label, row in values
        ))


if __name__ == "__main__":
    main()
