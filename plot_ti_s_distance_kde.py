#!/usr/bin/env python3
"""KDE of reference M–L bond distances from tmQM_metal_bond_data_CIF.pkl."""

from __future__ import annotations

import json
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
PKL = ROOT / "tmQM_metal_bond_data_CIF.pkl"
CCDC_JSON = ROOT / "ccdc_covalent_radii.json"
LIMITS_JSON = ROOT / "tm_nonmetal_bond_limits.json"
COV_MARGIN = 0.45
P99_MARGIN = 0.05
FIG_DPI = 300

# Match plot_covradius_audit_top5.py
LABEL_COV = r"$R_{\mathrm{cov}}(M)+R_{\mathrm{cov}}(L)+0.45$ Å"
LABEL_P99 = "99th %ile of M–L distances + 0.05 Å"
COLOR_COV = "#C44E52"
COLOR_P99 = "#F4A582"

PAIRS: tuple[tuple[str, str], ...] = (("Hg", "O"), ("La", "N"))


def gaussian_kde_1d(samples: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Silverman bandwidth Gaussian KDE without scipy."""
    x = np.asarray(samples, dtype=float)
    std = x.std(ddof=1) if len(x) > 1 else 1.0
    bw = 1.06 * std * len(x) ** (-1 / 5) if std > 0 else 0.1
    diffs = (grid[:, None] - x[None, :]) / bw
    return np.exp(-0.5 * diffs**2).sum(axis=1) / (len(x) * bw * np.sqrt(2 * np.pi))


def load_pair_distances(metal: str, ligand: str) -> np.ndarray:
    with open(PKL, "rb") as fh:
        data = pickle.load(fh)
    dp = data["analysis"]["distances_by_pair"]
    for key in ((metal, ligand, "TM"), (ligand, metal, "TM")):
        if key in dp:
            return np.asarray(dp[key], dtype=float)
    raise KeyError(f"{metal}-{ligand} distances not found in analysis.distances_by_pair")


def load_cov_cutoff(metal: str, ligand: str) -> tuple[float, float, float]:
    with open(CCDC_JSON, encoding="utf-8") as fh:
        ccdc = json.load(fh)["by_symbol"]
    r_m, r_l = float(ccdc[metal]), float(ccdc[ligand])
    return r_m, r_l, r_m + r_l + COV_MARGIN


def load_p99_cutoff(metal: str, ligand: str) -> float:
    pair = f"{metal}-{ligand}"
    with open(LIMITS_JSON, encoding="utf-8") as fh:
        stats = json.load(fh)["stats"][pair]
    return float(stats["p99_A"]) + P99_MARGIN


def plot_pair_kde(metal: str, ligand: str, out_path: Path) -> None:
    dists = load_pair_distances(metal, ligand)
    r_m, r_l, cov_cut = load_cov_cutoff(metal, ligand)
    p99_cut = load_p99_cutoff(metal, ligand)

    x_min = min(dists.min() - 0.08, cov_cut - 0.35)
    x_max = max(cov_cut + 0.12, p99_cut + 0.12, dists.max() + 0.08)
    xs = np.linspace(x_min, x_max, 400)
    ys = gaussian_kde_1d(dists, xs)

    fig, ax = plt.subplots(figsize=(8.5, 5))
    ax.fill_between(xs, ys, alpha=0.35, color="#4C72B0")
    ax.plot(
        xs,
        ys,
        color="#4C72B0",
        linewidth=2,
        label=f"CSD {metal}–{ligand} bonds (n={len(dists):,})",
    )

    ax.axvline(cov_cut, color=COLOR_COV, linestyle="--", linewidth=1.8, label=LABEL_COV)
    ax.axvline(p99_cut, color=COLOR_P99, linestyle="-.", linewidth=1.8, label=LABEL_P99)

    ax.set_xlabel(f"{metal}–{ligand} distance (Å)", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.set_title(f"{metal}–{ligand} Bond Length Distribution", fontsize=13, pad=10)
    ax.set_xlim(x_min, x_max)
    ax.yaxis.grid(True, linestyle="--", alpha=0.3)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", fontsize=9.5, framealpha=0.95)

    fig.tight_layout()
    fig.savefig(out_path, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")
    print(f"  n={len(dists)}, range [{dists.min():.3f}, {dists.max():.3f}] Å")
    print(f"  R_cov({metal})={r_m:.2f}, R_cov({ligand})={r_l:.2f}, cov+{COV_MARGIN}={cov_cut:.2f} Å")
    print(f"  P99+{P99_MARGIN}={p99_cut:.3f} Å")


def main() -> None:
    for metal, ligand in PAIRS:
        out = ROOT / f"{metal}-{ligand}_distance_KDE.png"
        plot_pair_kde(metal, ligand, out)


if __name__ == "__main__":
    main()
