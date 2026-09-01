#!/usr/bin/env python3
"""Plot whole-dataset missed and spurious M–L bonds across P99 margins."""

from __future__ import annotations

import csv
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
FIG_DPI = 300
REFERENCE_BONDS = 491_502
REFERENCE_NONBONDS = 5_532_159

AUDITS = (
    ("P99", ROOT / "tmQM_covradius_audit_cutoff_P99.pkl"),
    ("P99 + 0.05 Å", ROOT / "tmQM_covradius_audit_cutoff_P99_0p05.pkl"),
    ("P99 + 0.15 Å", ROOT / "tmQM_covradius_audit_cutoff.pkl"),
)


def load_pair_records(path: Path) -> dict[str, dict]:
    with path.open("rb") as fh:
        return pickle.load(fh)["by_element_pair"]


def apply_fixed_ni_n_override(records: dict[str, dict], fixed_record: dict) -> dict[str, dict]:
    """Keep the engine's 2.620 Å Mercury GUI Ni–N exception fixed across margins."""
    out = {pair: dict(rec) for pair, rec in records.items()}
    out["Ni-N"] = dict(fixed_record)
    return out


def totals(records: dict[str, dict]) -> dict[str, int]:
    return {
        key: sum(rec[key] for rec in records.values())
        for key in ("tp", "fp", "fn", "tn")
    }


def collect_rows() -> list[dict]:
    loaded = [(label, load_pair_records(path)) for label, path in AUDITS]
    fixed_ni_n = loaded[-1][1]["Ni-N"]
    rows = []
    for label, records in loaded:
        result = totals(apply_fixed_ni_n_override(records, fixed_ni_n))
        tp, fp, fn, tn = (result[key] for key in ("tp", "fp", "fn", "tn"))
        rows.append(
            {
                "cutoff": label,
                "reference_bonds": tp + fn,
                "missed_bonds": fn,
                "missed_bond_rate_pct": 100 * fn / (tp + fn),
                "reference_nonbonds": fp + tn,
                "spurious_bonds": fp,
                "spurious_reference_bond_ratio_pct": 100 * fp / (tp + fn),
                "spurious_nonbond_rate_pct": 100 * fp / (fp + tn),
                "predicted_bonds": tp + fp,
                "spurious_predicted_fraction_pct": 100 * fp / (tp + fp),
                "total_errors": fn + fp,
            }
        )
    return rows


def write_csv(rows: list[dict]) -> Path:
    out = ROOT / "audit_margin_overall_FN_FP.csv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return out


def plot(rows: list[dict]) -> Path:
    labels = [row["cutoff"] for row in rows]
    x = np.arange(len(rows))
    colors = ("#4C72B0", "#79A6D2", "#B9D4EA")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.9))

    panels = (
        ("missed_bonds", "missed_bond_rate_pct", "Missed reference bonds", "Missed bonds (count)"),
        (
            "spurious_bonds",
            "spurious_reference_bond_ratio_pct",
            "Spurious predicted bonds",
            "Spurious bonds (count)",
        ),
    )
    for ax, (count_key, rate_key, title, ylabel) in zip(axes, panels):
        counts = [row[count_key] for row in rows]
        bars = ax.bar(x, counts, width=0.62, color=colors, edgecolor="#333333", linewidth=0.8)
        ax.set_xticks(x, labels)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.yaxis.grid(True, linestyle="--", alpha=0.35)
        ax.set_axisbelow(True)
        ax.set_ylim(0, max(counts) * 1.23)
        for bar, row in zip(bars, rows):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(counts) * 0.025,
                f"{row[count_key]:,}\n({row[rate_key]:.3f}%)",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    fig.suptitle(
        "Whole-dataset M–L connectivity sensitivity to the P99 distance margin",
        fontsize=13,
        y=1.01,
    )
    fig.text(
        0.5,
        -0.015,
        (
            f"Reference: {REFERENCE_BONDS:,} CSD M–L bonds and "
            f"{REFERENCE_NONBONDS:,} nonbonded metal–nonmetal pairs. "
            "Percentages in both panels are relative to the reference-bond count. "
            "Ni–N fixed at the 2.620 Å Mercury GUI limit."
        ),
        ha="center",
        fontsize=8.5,
    )
    fig.tight_layout()
    out = ROOT / "audit_margin_overall_FN_FP.png"
    fig.savefig(out, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    rows = collect_rows()
    csv_path = write_csv(rows)
    fig_path = plot(rows)
    print("Wrote", csv_path)
    print("Wrote", fig_path)
    for row in rows:
        print(
            f"{row['cutoff']}: missed {row['missed_bonds']:,}/{row['reference_bonds']:,} "
            f"({row['missed_bond_rate_pct']:.3f}%), spurious "
            f"{row['spurious_bonds']:,}/{row['reference_bonds']:,} reference bonds "
            f"({row['spurious_reference_bond_ratio_pct']:.3f}%)"
        )


if __name__ == "__main__":
    main()
