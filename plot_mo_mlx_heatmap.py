#!/usr/bin/env python3
"""MLX heatmaps for Group-6 metals (Cr, Mo, W) in tmQM-G (ILP + CBC classification)."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

ROOT = Path(__file__).resolve().parent
ENGINE_PATH = ROOT / "Lewis-engine-ILP.py"
CSV_PATH = ROOT / "tmqmg_ilp_benchmark_output" / "per_structure_results.csv"
DEFAULT_XYZ_DIR = Path("/data/jingyuan_data/tmqmg")
SKIP_PATHS = (
    ROOT / "list_error_geometry.txt",
    ROOT / "csd_codes_boron_ge6.txt",
)
METAL_ORDER = ("Cr", "Mo", "W")
METAL_VALENCE = 6  # neutral Cr/Mo/W d+s electron count (Group 6)
METAL_TITLES = {
    "Cr": "chromium",
    "Mo": "molybdenum",
    "W": "tungsten",
}
DEFAULT_SUMMARY_JSON = {
    "Cr": ROOT / "cr_mlx_heatmap_summary.json",
    "Mo": ROOT / "mo_mlx_heatmap_summary.json",
    "W": ROOT / "w_mlx_heatmap_summary.json",
}
OUT_COMBINED_PNG = ROOT / "group6_mlx_heatmap.png"
OUT_REVISED_PNG = ROOT / "group6_mlx_heatmap_revised.png"
_LX_LABEL_RE = re.compile(r"(?:L(\d+))?(?:X(\d+))?$")
_SUB_TO_ASCII = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")
FIG_DPI = 300
# Horizontal gap between Cr / Mo / W panels (matplotlib wspace; larger = wider gap).
PANEL_WSPACE = 0.10
# Font size for count numbers inside heatmap cells.
COUNT_FONTSIZE = 6
# Font size for MLX formula labels (e.g. MoL₃X₂) inside heatmap cells.
MLX_FONTSIZE = 6.5
# Font sizes for axis tick numbers and axis titles (E.N. / V.N.).
TICK_FONTSIZE = 7
AXIS_LABEL_FONTSIZE = 8

EN_MIN, EN_MAX = 12, 18
VN_MIN, VN_MAX = 0, 6
CHARGE_RANGE = (-2, -1, 0, 1, 2)

NATURE_FONT_SANS = ["Arial", "Helvetica", "DejaVu Sans"]
COLOR_AXIS = "#444444"
COLOR_FORBIDDEN = "#8a8a8a"
COLOR_FORBIDDEN_EDGE = "#707070"
NATURE_HEATMAP_CMAP = LinearSegmentedColormap.from_list(
    "mlx_red",
    ["#ffffff", "#fee5d9", "#fcae91", "#fb6a4a", "#de2d26", "#a50f15"],
    N=256,
)
FORBIDDEN_BASE_CMAP = LinearSegmentedColormap.from_list(
    "forbidden_base", [COLOR_FORBIDDEN, COLOR_FORBIDDEN]
)


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
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.dpi": FIG_DPI,
        }
    )


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


def load_engine():
    if "rdkit" not in sys.modules:
        import types

        rdkit_stub = types.ModuleType("rdkit")
        chem_stub = types.ModuleType("rdkit.Chem")
        rdkit_stub.Chem = chem_stub
        sys.modules["rdkit"] = rdkit_stub
        sys.modules["rdkit.Chem"] = chem_stub

    spec = importlib.util.spec_from_file_location("lewis_engine_ilp", ENGINE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {ENGINE_PATH}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _subscript(n: int) -> str:
    return "".join("₀₁₂₃₄₅₆₇₈₉"[int(c)] for c in str(n))


def mlx_mathtext_for_cell(en: int, vn: int, metal: str) -> str:
    """Matplotlib mathtext for canonical ML_lX_x at grid (EN, VN), q=0."""
    l = (en - vn - METAL_VALENCE) // 2
    x = vn
    s = rf"\mathrm{{{metal}}}"
    if l:
        s += rf"\mathrm{{L}}_{{{l}}}" if l > 1 else r"\mathrm{L}"
    if x:
        s += rf"\mathrm{{X}}_{{{x}}}" if x > 1 else r"\mathrm{X}"
    return f"${s}$"


def mlx_label(l: int, x: int) -> str:
    parts = []
    if l:
        parts.append(f"L{_subscript(l)}")
    if x:
        parts.append(f"X{_subscript(x)}")
    return "".join(parts) if parts else "—"


def parse_lx_label(label: str) -> tuple[int, int]:
    text = (label or "").translate(_SUB_TO_ASCII).strip()
    match = _LX_LABEL_RE.fullmatch(text)
    if not match or not text:
        raise ValueError(f"Cannot parse MLX label: {label!r}")
    l = int(match.group(1)) if match.group(1) else 0
    x = int(match.group(2)) if match.group(2) else 0
    return l, x


def rebin_counts_en_no_charge(
    cell_counts: Counter,
    cell_labels: dict[tuple[int, int], str],
) -> Counter:
    """Re-bin (EN, VN) using EN = m + 2n + VN (drop molecular charge q)."""
    out: Counter = Counter()
    for (en_old, vn), count in cell_counts.items():
        label = cell_labels.get((en_old, vn), "")
        n, x = parse_lx_label(label)
        en_new = METAL_VALENCE + 2 * n + x
        out[(en_new, x)] += int(count)
    return out


def is_forbidden_mlx_cell(en: int, vn: int) -> bool:
    """Neutral MLX grid (m=6): forbidden unless EN = 6 + 2l + VN (checkerboard)."""
    if not (EN_MIN <= en <= EN_MAX and VN_MIN <= vn <= VN_MAX):
        return True
    delta = en - vn - METAL_VALENCE
    return delta < 0 or delta % 2 != 0


def is_valid_mlx_cell(en: int, vn: int) -> bool:
    if not (EN_MIN <= en <= EN_MAX and VN_MIN <= vn <= VN_MAX):
        return False
    for q in CHARGE_RANGE:
        l_num = en - vn - METAL_VALENCE + q
        if l_num < 0 or l_num % 2 != 0:
            continue
        l = l_num // 2
        if METAL_VALENCE + 2 * l + vn - q == en:
            return True
    return False


def mlx_formula_for_cell(en: int, vn: int) -> str:
    for q in (0, 1, -1, 2, -2):
        l_num = en - vn - METAL_VALENCE + q
        if l_num < 0 or l_num % 2 != 0:
            continue
        l = l_num // 2
        if METAL_VALENCE + 2 * l + vn - q == en:
            return mlx_label(l, vn)
    return ""


def classify_metal_centres(
    engine, xyz_path: Path, charge: int, metal: str
) -> list[tuple[int, int, int, int]]:
    """Return list of (l, x, en, vn) for each centre of the given metal."""
    atoms = engine.read_xyz(str(xyz_path))
    if not any(a[1] == metal for a in atoms):
        return []

    raw = engine.connectivity(atoms)
    aromatic_systems = engine.aromatic_candidate_systems(atoms, raw)
    bonds, lp_out, fc_out = engine.solve_bond_orders(
        atoms,
        raw,
        aromatic_systems,
        mol_charge=charge,
        metal_adjacency_edges=raw,
    )
    bonds, lp_out, fc_out, _carbene_labels = engine.apply_heterocyclic_carbene_corrections(
        atoms,
        bonds,
        lp_out,
        fc_out,
        aromatic_systems,
        raw,
        mol_charge=charge,
    )
    bonds, lp_out, fc_out = engine.apply_eta_covalent_pi_corrections(
        atoms,
        bonds,
        lp_out,
        fc_out,
        mol_charge=charge,
        metal_adjacency_edges=raw,
    )

    coords = [[a[2], a[3], a[4]] for a in atoms]
    atom_syms = [a[1] for a in atoms]
    fc_list = [fc_out.get(a[0], 0) for a in atoms]
    idx_to_pos = {a[0]: k for k, a in enumerate(atoms)}
    metal_adj_0 = [
        (idx_to_pos[tm], idx_to_pos[lig], ei, ej)
        for tm, lig, ei, ej in raw
        if engine.is_TM(ei) ^ engine.is_TM(ej)
    ]

    bo0: dict[tuple[int, int], int] = {}
    for i, j, order in bonds:
        ii, jj = i - 1, j - 1
        bo0[(min(ii, jj), max(ii, jj))] = int(order)

    lp0 = {k - 1: v for k, v in lp_out.items()}

    cbc_bundle = engine.classify_cbc_ligands(
        atom_syms,
        coords,
        bo0,
        lp0,
        fc_list,
        metal_adjacency_edges=metal_adj_0,
    )
    results, _ = cbc_bundle

    out: list[tuple[int, int, int, int]] = []
    for metal_idx, records in results.items():
        if atom_syms[metal_idx] != metal:
            continue
        l_count, x_count = engine.mlx_lx_counts_from_cbc_records(metal_idx, records, bo0)
        en, vn = engine.mlx_en_vn(METAL_VALENCE, l_count, x_count, charge)
        out.append((l_count, x_count, en, vn))
    return out


def collect_mlx(
    engine,
    metal: str,
    *,
    xyz_dir: Path,
    limit: int | None = None,
) -> tuple[Counter, dict[tuple[int, int], str], list[dict], dict]:
    skip = load_skip_ids()
    rows: list[tuple[str, int]] = []
    for row in csv.DictReader(CSV_PATH.open(encoding="utf-8")):
        if row["IDs"] in skip or row["status"] != "ok":
            continue
        if not re.search(rf"\[{metal}", row.get("smiles_ilp", "")):
            continue
        rows.append((row["IDs"], int(float(row["charge"] or 0))))

    if limit is not None:
        rows = rows[:limit]

    cell_counts: Counter = Counter()
    cell_label_votes: dict[tuple[int, int], Counter] = defaultdict(Counter)
    records: list[dict] = []
    errors: Counter = Counter()

    try:
        from tqdm import tqdm

        iterator = tqdm(rows, desc=f"{metal} MLX", unit="struct")
    except ImportError:
        iterator = rows

    for csd_id, charge in iterator:
        xyz_path = xyz_dir / f"{csd_id}.xyz"
        if not xyz_path.is_file():
            errors["missing_xyz"] += 1
            continue
        try:
            centres = classify_metal_centres(engine, xyz_path, charge, metal)
        except Exception as exc:
            errors[type(exc).__name__] += 1
            continue
        for l_count, x_count, en, vn in centres:
            key = (en, vn)
            label = mlx_label(l_count, x_count)
            cell_counts[key] += 1
            cell_label_votes[key][label] += 1
            records.append(
                {
                    "id": csd_id,
                    "charge": charge,
                    "l": l_count,
                    "x": x_count,
                    "en": en,
                    "vn": vn,
                    "label": label,
                }
            )

    cell_labels = {k: votes.most_common(1)[0][0] for k, votes in cell_label_votes.items()}
    return cell_counts, cell_labels, records, dict(errors)


def _build_grids(cell_counts: Counter) -> tuple[np.ndarray, np.ndarray, list[int], list[int]]:
    en_vals = list(range(EN_MIN, EN_MAX + 1))
    vn_vals = list(range(VN_MAX, VN_MIN - 1, -1))

    count_grid = np.zeros((len(vn_vals), len(en_vals)), dtype=float)
    forbidden = np.zeros((len(vn_vals), len(en_vals)), dtype=bool)

    for i, vn in enumerate(vn_vals):
        for j, en in enumerate(en_vals):
            if is_forbidden_mlx_cell(en, vn):
                forbidden[i, j] = True
                continue
            count_grid[i, j] = float(cell_counts.get((en, vn), 0))

    return count_grid, forbidden, en_vals, vn_vals


def _plot_heatmap_on_ax(
    ax,
    cell_counts: Counter,
    metal: str,
    *,
    vmax: float,
    show_ylabel: bool = True,
    count_fontsize: float = COUNT_FONTSIZE,
    mlx_fontsize: float = MLX_FONTSIZE,
    tick_fontsize: float = TICK_FONTSIZE,
    axis_label_fontsize: float = AXIS_LABEL_FONTSIZE,
) -> plt.cm.ScalarMappable:
    count_grid, forbidden, en_vals, vn_vals = _build_grids(cell_counts)

    ax.imshow(
        np.ones((len(vn_vals), len(en_vals))),
        cmap=FORBIDDEN_BASE_CMAP,
        vmin=0,
        vmax=1,
        aspect="equal",
        origin="lower",
        interpolation="nearest",
        zorder=1,
    )

    masked_counts = np.ma.masked_where(forbidden, count_grid)
    im = ax.imshow(
        masked_counts,
        cmap=NATURE_HEATMAP_CMAP,
        vmin=0.0,
        vmax=vmax,
        aspect="equal",
        origin="lower",
        interpolation="nearest",
        zorder=2,
    )

    for i in range(len(vn_vals)):
        for j in range(len(en_vals)):
            if forbidden[i, j]:
                ax.add_patch(
                    plt.Rectangle(
                        (j - 0.5, i - 0.5),
                        1,
                        1,
                        facecolor=COLOR_FORBIDDEN,
                        edgecolor=COLOR_FORBIDDEN_EDGE,
                        lw=0.35,
                        zorder=3,
                    )
                )
                continue
            en, vn = en_vals[j], vn_vals[i]
            count = int(count_grid[i, j])
            mlx_txt = mlx_mathtext_for_cell(en, vn, metal)
            count_color = "#111111" if count > 0 else COLOR_AXIS
            ax.text(
                j,
                i + (0.16 if count > 0 else 0.0),
                mlx_txt,
                ha="center",
                va="center",
                fontsize=mlx_fontsize,
                color="#111111",
                zorder=10,
            )
            if count > 0:
                ax.text(
                    j,
                    i - 0.24,
                    f"{count}",
                    ha="center",
                    va="center",
                    fontsize=count_fontsize,
                    color=count_color,
                    zorder=10,
                )

    ax.set_xticks(range(len(en_vals)))
    ax.set_xticklabels([str(v) for v in en_vals])
    ax.set_yticks(range(len(vn_vals)))
    ax.set_yticklabels([str(v) for v in vn_vals])
    ax.set_xlabel("E.N.", fontsize=axis_label_fontsize)
    if show_ylabel:
        ax.set_ylabel("V.N.", fontsize=axis_label_fontsize)
    else:
        ax.set_ylabel("")
        ax.tick_params(labelleft=True)
    ax.tick_params(
        colors=COLOR_AXIS,
        width=0.6,
        length=3,
        labelsize=tick_fontsize,
    )
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(COLOR_AXIS)
        ax.spines[spine].set_linewidth(0.6)

    return im


def _vmax_for_counts(cell_counts: Counter) -> float:
    count_grid, forbidden, _, _ = _build_grids(cell_counts)
    allowed = count_grid[~forbidden]
    return max(float(allowed.max()) if allowed.size else 1.0, 1.0)


def plot_combined_heatmap(
    metal_counts: dict[str, Counter],
    out_png: Path,
    *,
    wspace: float = PANEL_WSPACE,
    count_fontsize: float = COUNT_FONTSIZE,
    mlx_fontsize: float = MLX_FONTSIZE,
    tick_fontsize: float = TICK_FONTSIZE,
    axis_label_fontsize: float = AXIS_LABEL_FONTSIZE,
) -> None:
    _apply_nature_style()

    n = len(METAL_ORDER)
    fig, axes = plt.subplots(1, n, figsize=(3.75 * n, 4.5))
    if n == 1:
        axes = [axes]
    fig.subplots_adjust(wspace=wspace, left=0.07, right=0.99, top=0.96, bottom=0.14)

    for ax, metal in zip(axes, METAL_ORDER):
        vmax = _vmax_for_counts(metal_counts[metal])
        im = _plot_heatmap_on_ax(
            ax,
            metal_counts[metal],
            metal,
            vmax=vmax,
            show_ylabel=(metal == METAL_ORDER[0]),
            count_fontsize=count_fontsize,
            mlx_fontsize=mlx_fontsize,
            tick_fontsize=tick_fontsize,
            axis_label_fontsize=axis_label_fontsize,
        )
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="4.5%", pad=0.05)
        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label(
            f"Number of {metal} complex",
            fontsize=6.5,
            color=COLOR_AXIS,
            labelpad=3,
        )
        cbar.ax.tick_params(labelsize=6, colors=COLOR_AXIS, width=0.6, length=2.0)
        cbar.outline.set_linewidth(0.6)
        cbar.outline.set_edgecolor(COLOR_AXIS)

    fig.savefig(out_png, dpi=FIG_DPI, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


def plot_single_heatmap(
    cell_counts: Counter,
    metal: str,
    out_png: Path,
    *,
    count_fontsize: float = COUNT_FONTSIZE,
    mlx_fontsize: float = MLX_FONTSIZE,
    tick_fontsize: float = TICK_FONTSIZE,
    axis_label_fontsize: float = AXIS_LABEL_FONTSIZE,
) -> None:
    _apply_nature_style()
    vmax = _vmax_for_counts(cell_counts)

    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    im = _plot_heatmap_on_ax(
        ax,
        cell_counts,
        metal,
        vmax=vmax,
        show_ylabel=True,
        count_fontsize=count_fontsize,
        mlx_fontsize=mlx_fontsize,
        tick_fontsize=tick_fontsize,
        axis_label_fontsize=axis_label_fontsize,
    )
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label(
        f"Number of {metal} complex",
        fontsize=8,
        color=COLOR_AXIS,
        labelpad=6,
    )
    cbar.ax.tick_params(labelsize=7, colors=COLOR_AXIS, width=0.6, length=2.5)
    cbar.outline.set_linewidth(0.6)
    cbar.outline.set_edgecolor(COLOR_AXIS)
    fig.savefig(out_png, dpi=FIG_DPI, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


def load_summary_json(path: Path) -> tuple[Counter, dict[tuple[int, int], str]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    cell_counts = Counter({tuple(map(int, k.split(","))): v for k, v in data["cell_counts"].items()})
    cell_labels = {
        tuple(map(int, k.split(","))): v for k, v in data.get("cell_labels", {}).items()
    }
    return cell_counts, cell_labels


def save_summary_json(
    path: Path,
    metal: str,
    cell_counts: Counter,
    cell_labels: dict[tuple[int, int], str],
    records: list[dict],
    errors: dict,
) -> None:
    summary = {
        "metal": metal,
        "m_valence": METAL_VALENCE,
        "en_formula": (
            "EN = m + 2l + x - q; VN = x; "
            f"x = sum of ILP {metal}-L bond orders on CBC-X records "
            "(single/double/triple -> 1/2/3)"
        ),
        "n_structures_processed": len({r["id"] for r in records}),
        f"n_{metal.lower()}_centres": sum(cell_counts.values()),
        "cell_counts": {f"{en},{vn}": c for (en, vn), c in sorted(cell_counts.items())},
        "cell_labels": {f"{en},{vn}": lab for (en, vn), lab in sorted(cell_labels.items())},
        "top_cells": [
            {
                "en": en,
                "vn": vn,
                "count": c,
                "label": cell_labels.get((en, vn), mlx_formula_for_cell(en, vn)),
            }
            for (en, vn), c in cell_counts.most_common(12)
        ],
        "errors": errors,
    }
    path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--xyz-dir", type=Path, default=DEFAULT_XYZ_DIR)
    p.add_argument("--limit", type=int, default=None, metavar="N")
    p.add_argument(
        "--metal",
        choices=[*METAL_ORDER, "all"],
        default="all",
        help="single metal plot, or all three combined (default: all)",
    )
    p.add_argument("--out-png", type=Path, default=None)
    p.add_argument(
        "--en-no-charge",
        action="store_true",
        help="plot EN = m + 2n + VN (no −q); rebin JSON counts onto the q=0 checkerboard",
    )
    p.add_argument(
        "--from-json",
        type=Path,
        nargs="*",
        default=None,
        help="load summaries (Cr Mo W order, or one file for --metal)",
    )
    p.add_argument(
        "--recompute",
        action="store_true",
        help="recompute even if summary JSON exists",
    )
    p.add_argument(
        "--wspace",
        type=float,
        default=PANEL_WSPACE,
        metavar="W",
        help=f"horizontal gap between panels (default: {PANEL_WSPACE})",
    )
    p.add_argument(
        "--count-fontsize",
        type=float,
        default=COUNT_FONTSIZE,
        metavar="FS",
        help=f"font size for count numbers in cells (default: {COUNT_FONTSIZE})",
    )
    p.add_argument(
        "--mlx-fontsize",
        type=float,
        default=MLX_FONTSIZE,
        metavar="FS",
        help=f"font size for MLX formula labels in cells (default: {MLX_FONTSIZE})",
    )
    p.add_argument(
        "--tick-fontsize",
        type=float,
        default=TICK_FONTSIZE,
        metavar="FS",
        help=f"font size for X/Y axis tick numbers (default: {TICK_FONTSIZE})",
    )
    p.add_argument(
        "--axis-label-fontsize",
        type=float,
        default=AXIS_LABEL_FONTSIZE,
        metavar="FS",
        help=f"font size for E.N. / V.N. axis titles (default: {AXIS_LABEL_FONTSIZE})",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if args.metal == "all":
        metals = METAL_ORDER
        if args.out_png:
            out_png = args.out_png
        elif args.en_no_charge:
            out_png = OUT_REVISED_PNG
        else:
            out_png = OUT_COMBINED_PNG
    else:
        metals = (args.metal,)
        if args.out_png:
            out_png = args.out_png
        elif args.en_no_charge:
            out_png = ROOT / f"{args.metal.lower()}_mlx_heatmap_revised.png"
        else:
            out_png = ROOT / f"{args.metal.lower()}_mlx_heatmap.png"

    from_json = list(args.from_json) if args.from_json is not None else None
    metal_counts: dict[str, Counter] = {}

    engine = None
    for idx, metal in enumerate(metals):
        summary_path = DEFAULT_SUMMARY_JSON[metal]
        loaded_from_json = False

        if from_json is not None:
            json_path = from_json[idx] if len(from_json) > idx else from_json[0]
            cell_counts, cell_labels = load_summary_json(json_path)
            if args.en_no_charge:
                cell_counts = rebin_counts_en_no_charge(cell_counts, cell_labels)
            metal_counts[metal] = cell_counts
            loaded_from_json = True
        elif not args.recompute and summary_path.is_file():
            cell_counts, cell_labels = load_summary_json(summary_path)
            if args.en_no_charge:
                cell_counts = rebin_counts_en_no_charge(cell_counts, cell_labels)
            metal_counts[metal] = cell_counts
            loaded_from_json = True
            print(f"Loaded {metal} from {summary_path}")

        if loaded_from_json:
            continue

        if engine is None:
            if not CSV_PATH.is_file():
                raise SystemExit(f"CSV not found: {CSV_PATH}")
            if not args.xyz_dir.is_dir():
                raise SystemExit(f"XYZ directory not found: {args.xyz_dir}")
            engine = load_engine()

        cell_counts, cell_labels, records, errors = collect_mlx(
            engine, metal, xyz_dir=args.xyz_dir, limit=args.limit
        )
        save_summary_json(summary_path, metal, cell_counts, cell_labels, records, errors)
        if args.en_no_charge:
            cell_counts = rebin_counts_en_no_charge(cell_counts, cell_labels)
        metal_counts[metal] = cell_counts

    if args.metal == "all":
        plot_combined_heatmap(
            metal_counts,
            out_png,
            wspace=args.wspace,
            count_fontsize=args.count_fontsize,
            mlx_fontsize=args.mlx_fontsize,
            tick_fontsize=args.tick_fontsize,
            axis_label_fontsize=args.axis_label_fontsize,
        )
    else:
        plot_single_heatmap(
            metal_counts[metals[0]],
            metals[0],
            out_png,
            count_fontsize=args.count_fontsize,
            mlx_fontsize=args.mlx_fontsize,
            tick_fontsize=args.tick_fontsize,
            axis_label_fontsize=args.axis_label_fontsize,
        )

    print(f"Wrote {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
