#!/usr/bin/env python3
"""Compare tmQM M-L pairs in /data/jingyuan_data/tmqmg vs tm_nonmetal_bond_limits.json stats."""

import json
import math
import os
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.abspath(__file__))
LIMITS_JSON = os.path.join(ROOT, "tm_nonmetal_bond_limits.json")
CCDC_JSON = os.path.join(ROOT, "ccdc_covalent_radii.json")
TMQMG_DIR = "/data/jingyuan_data/tmqmg"

TM_SET = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu",
}
S_BLOCK_SYMS = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "Tl", "Pb", "Bi", "Po", "At", "Rn"}
COV_BOND_MARGIN = 0.45
COV_BOND_MARGIN_S_BLOCK = 0.40
TM_NONMETAL_P99_MARGIN_A = 0.05
PROBE_MARGIN_A = 0.55


def load_cov_radii():
    with open(CCDC_JSON, encoding="utf-8") as fh:
        data = json.load(fh)
    by = data.get("by_symbol") or data
    return {k: float(v) for k, v in by.items() if isinstance(v, (int, float))}


def load_tm_limits():
    with open(LIMITS_JSON, encoding="utf-8") as fh:
        data = json.load(fh)
    stats = data.get("stats") or {}
    limits = {}
    for pair, rec in stats.items():
        if isinstance(rec, dict) and rec.get("p99_A") is not None:
            limits[str(pair)] = float(rec["p99_A"]) + TM_NONMETAL_P99_MARGIN_A
    return limits, set(stats.keys())


COV_R = load_cov_radii()
TM_LIMITS, STATS_KEYS = load_tm_limits()


def read_xyz(path):
    with open(path, encoding="utf-8") as f:
        lines = f.readlines()
    n = int(lines[0].strip().split()[0])
    atoms = []
    for k in range(n):
        p = lines[2 + k].split()
        raw = p[0]
        sym = raw[0].upper() + raw[1:].lower() if len(raw) > 1 else raw.upper()
        x, y, z = float(p[1]), float(p[2]), float(p[3])
        atoms.append((k + 1, sym, x, y, z))
    return atoms


def dist(a, b):
    return math.sqrt((a[2] - b[2]) ** 2 + (a[3] - b[3]) ** 2 + (a[4] - b[4]) ** 2)


def bond_cutoff_cov(ei, ej):
    margin = COV_BOND_MARGIN_S_BLOCK if ei in S_BLOCK_SYMS or ej in S_BLOCK_SYMS else COV_BOND_MARGIN
    ri = COV_R.get(ei)
    rj = COV_R.get(ej)
    if ri is None or rj is None:
        return None
    return ri + rj + margin


def tm_nonmetal_cutoff(metal, ligand):
    return TM_LIMITS.get(f"{metal}-{ligand}")


def probe_cutoff(metal, ligand):
    lim = tm_nonmetal_cutoff(metal, ligand)
    if lim is not None:
        return lim
    cov = bond_cutoff_cov(metal, ligand)
    if cov is None:
        return 4.0
    return cov + PROBE_MARGIN_A - COV_BOND_MARGIN


def main():
    pair_counts = Counter()
    pair_distances = defaultdict(list)
    structures_with_tm = 0
    n_files = 0
    n_skip_no_cov = 0

    for name in sorted(os.listdir(TMQMG_DIR)):
        if not name.endswith(".xyz"):
            continue
        path = os.path.join(TMQMG_DIR, name)
        if not os.path.isfile(path):
            continue
        n_files += 1
        try:
            atoms = read_xyz(path)
        except Exception:
            continue
        tms = [(i, el) for i, el, *_ in atoms if el in TM_SET]
        if not tms:
            continue
        structures_with_tm += 1
        non_tm = [(i, el) for i, el, *_ in atoms if el not in TM_SET]
        atom_by_idx = {i: a for i, a in enumerate(atoms, start=1)}
        for tm_idx, tm_sym in tms:
            tm_a = atom_by_idx[tm_idx]
            for lig_idx, lig_sym in non_tm:
                cov = bond_cutoff_cov(tm_sym, lig_sym)
                if cov is None:
                    n_skip_no_cov += 1
                    continue
                d = dist(tm_a, atom_by_idx[lig_idx])
                if d >= probe_cutoff(tm_sym, lig_sym):
                    continue
                key = f"{tm_sym}-{lig_sym}"
                pair_counts[key] += 1
                pair_distances[key].append(d)
        if n_files % 10000 == 0:
            print(f"  scanned {n_files} ...", flush=True)

    observed = set(pair_counts.keys())
    missing = sorted(observed - STATS_KEYS)

    print(f"JSON stats pairs: {len(STATS_KEYS)}")
    print(f"Scanned xyz: {n_files}, with TM: {structures_with_tm}")
    print(f"Observed M-L pairs (probe): {len(observed)}")
    print(f"Missing from JSON: {len(missing)}")

    print("\n=== Missing pairs (sorted by count) ===")
    for key in sorted(missing, key=lambda k: (-pair_counts[k], k)):
        ds = sorted(pair_distances[key])
        p99 = ds[int(0.99 * (len(ds) - 1))]
        cov_fb = bond_cutoff_cov(*key.split("-", 1))
        print(
            f"  {key:12s}  n={pair_counts[key]:6d}  "
            f"d=[{ds[0]:.3f},{ds[-1]:.3f}]  p99~{p99:.3f}  cov_fb={cov_fb:.3f}"
        )

    by_metal = defaultdict(list)
    for key in missing:
        m, l = key.split("-", 1)
        by_metal[m].append(l)
    print("\n=== By metal (missing ligand symbols) ===")
    for m in sorted(by_metal, key=lambda x: -len(by_metal[x])):
        print(f"  {m}: {sorted(by_metal[m])}")

    out_path = os.path.join(ROOT, "tmqmg_ml_pair_coverage_gap.json")
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(
            {
                "n_xyz": n_files,
                "n_structures_with_tm": structures_with_tm,
                "missing_pairs": [
                    {
                        "pair": k,
                        "count": pair_counts[k],
                        "d_min": min(pair_distances[k]),
                        "d_max": max(pair_distances[k]),
                        "cov_fallback_A": bond_cutoff_cov(*k.split("-", 1)),
                    }
                    for k in missing
                ],
            },
            fh,
            indent=2,
        )
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
