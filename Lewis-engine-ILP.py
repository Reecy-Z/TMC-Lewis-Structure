#!/usr/bin/env python3
"""
Lewis-engine-ILP.py

"""

from __future__ import annotations

import itertools
import json
import math
import os
import sys
from collections import defaultdict

LOCAL_VENDOR = os.path.join(os.path.dirname(__file__), "_vendor")
if os.path.isdir(LOCAL_VENDOR) and LOCAL_VENDOR not in sys.path:
    sys.path.insert(0, LOCAL_VENDOR)

import pulp
from rdkit import Chem

# ---------------------------------------------------------------------------
# ILP configuration (solve_bond_orders)
#
#   ILP_HARD_*    — true hard constraints; False drops them (may become infeasible
#                   if nothing else compensates).
#   ILP_WEIGHT_*  — objective penalties only (0 = off).
# ---------------------------------------------------------------------------

# --- Hard constraints (must hold exactly) ---
# Hard octet: oct_plus == oct_minus == 0.
ILP_HARD_OCTET = True
# Σfc(ligands) + Σox(TM) − Σb_tm = mol_charge.
ILP_HARD_MOL_CHARGE_BALANCE = True
# Σox(TM) ≥ Σb_tm on non-η M–L only (η orders excluded). Halogen (F/Cl/Br/I) ligands are omitted.
ILP_HARD_OX_GE_SIGMA = True
# C with no TM neighbor in connectivity: prefer lp = 0 (soft penalty when enabled).
ILP_HARD_C_LP_ONLY_TM_NEIGHBORS = True
# η-fragment carbons to a TM (≥ETA_MIN_COORDINATING_GROUP_SIZE): lp = 0.
ILP_HARD_ETA_CARBON_LP_ZERO = True

# --- Soft objective weights (0 = disabled) ---
ILP_WEIGHT_FORMAL_CHARGE = 100.0
ILP_WEIGHT_AROMATIC_DEVIATION = 100.0
ILP_WEIGHT_ENEG_NEGATIVE_FC = 10.0
ILP_WEIGHT_ML_DISTANCE_CLASS = 50.0
ILP_WEIGHT_ETA_GROUP_MAX_DOUBLE_BONDS = 25.0
ILP_WEIGHT_REMOTE_C_LP_VIOLATION = 100.0

# Aromatic 4n+2 target: fixed n in pi_target = 4*n + 2 (default n=1 → 6 pi e per system).
# Set to None in solve_bond_orders(..., aromatic_huckel_n=None) to restore variable k.
AROMATIC_HUCKEL_N = 1

# Soft tie-break: similar M–L contact distances on the same ligand → same σ/dative class (z_cov).
ML_DISTANCE_CLASS_EPSILON = 0.15  # Å; w_ij = max(0, ε − |d_i − d_j|)

# η fragments: ≥2 contiguous TM-bound atoms on one ligand (see _eta_carbon_atom_ids).
ETA_MIN_COORDINATING_GROUP_SIZE = 2

# CBC wall-time limit for prob.solve inside solve_bond_orders (0 = no limit).
SOLVE_BOND_ORDERS_CBC_TIME_LIMIT_SEC = 10.0


# Covalent radii: ccdc_covalent_radii.json (CCDC ChemistryLib; Z=1–118 + D). TM–L cutoffs use
# p99_A + margin; other pairs use R_cov sum + COV_BOND_MARGIN.
COV_BOND_MARGIN = 0.45
COV_BOND_MARGIN_S_BLOCK = 0.40
S_BLOCK_SYMS = frozenset({
    "Li", "Be", "Na", "Mg", "K", "Ca", "Rb", "Sr", "Cs", "Ba",
})

_CCDC_COV_RADII_JSON = os.path.join(
    os.path.dirname(__file__), "ccdc_covalent_radii.json"
)

_TM_NONMETAL_LIMITS_JSON = os.path.join(
    os.path.dirname(__file__), "tmQM_tm_nonmetal_bond_limits.json"
)
# TM–L connectivity cutoff: empirical p99_A + this margin (Å)
TM_NONMETAL_P99_MARGIN_A = 0.05

TM_SET = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu",
}

TM_COMMON_OXIDATION_STATES = {
    "Sc": [3],
    "Ti": [3, 4],
    "V": [2, 3, 4, 5],
    "Cr": [0, 2, 3, 4, 6],
    "Mn": [1, 2, 3, 4, 6, 7],
    "Fe": [2, 3],
    "Co": [1, 2, 3],
    "Ni": [2],
    "Cu": [1, 2],
    "Zn": [2],
    "Y": [3],
    "Zr": [4],
    "Nb": [3, 4, 5],
    "Mo": [0, 2, 3, 4, 5, 6],
    "Tc": [2, 3, 4, 5, 6, 7],
    "Ru": [2, 3, 4, 5, 6, 7, 8],
    "Rh": [1, 3],
    "Pd": [2, 4],
    "Ag": [1],
    "Cd": [1, 2],
    "La": [3],
    "Hf": [4],
    "Ta": [3, 4, 5],
    "W": [0, 2, 3, 4, 5, 6],
    "Re": [1, 2, 3, 4, 5, 6, 7],
    "Os": [2, 3, 4, 5, 6, 7, 8],
    "Ir": [1, 3],
    "Pt": [2, 4],
    "Au": [1, 3],
    "Hg": [1, 2],
}

# Standard covalent bond capacities
STD_CAP = {
    "H": 1, "He": 0, "Li": 1, "Be": 2, "B": 3, "C": 4, "N": 3, "O": 2, "F": 1, "Ne": 0,
    "Na": 1, "Mg": 2, "Al": 3, "Si": 4, "P": 3, "S": 2, "Cl": 1, "Ar": 0,
    "K": 1, "Ca": 2, "Ga": 3, "Ge": 4, "As": 3, "Se": 2, "Br": 1, "Kr": 0,
    "Rb": 1, "Sr": 2, "In": 3, "Sn": 4, "Sb": 3, "Te": 2, "I": 1, "Xe": 0,
    "Cs": 1, "Ba": 2, "Tl": 3, "Pb": 4, "Bi": 3,
}

# Neutral-atom valence electron totals 
VALENCE = {
    "H": 1, "He": 2,
    "Li": 1, "Be": 2, "B": 3, "C": 4, "N": 5, "O": 6, "F": 7, "Ne": 8,
    "Na": 1, "Mg": 2, "Al": 3, "Si": 4, "P": 5, "S": 6, "Cl": 7, "Ar": 8,
    "K": 1, "Ca": 2, "Ga": 3, "Ge": 4, "As": 5, "Se": 6, "Br": 7, "Kr": 8,
    "Rb": 1, "Sr": 2, "In": 3, "Sn": 4, "Sb": 5, "Te": 6, "I": 7, "Xe": 8,
    "Cs": 1, "Ba": 2, "Tl": 3, "Pb": 4, "Bi": 5,
    # TM: d+s electrons of neutral ground-state atom
    "Sc": 3, "Ti": 4, "V": 5, "Cr": 6, "Mn": 7, "Fe": 8, "Co": 9, "Ni": 10, "Cu": 11, "Zn": 12,
    "Y": 3, "Zr": 4, "Nb": 5, "Mo": 6, "Tc": 7, "Ru": 8, "Rh": 9, "Pd": 10, "Ag": 11, "Cd": 12,
    "Hf": 4, "Ta": 5, "W": 6, "Re": 7, "Os": 8, "Ir": 9, "Pt": 10, "Au": 11, "Hg": 12,
    "La": 3, "Ce": 4, "Pr": 5, "Nd": 6, "Pm": 7, "Sm": 8, "Eu": 9, "Gd": 10, "Tb": 11,
    "Dy": 12, "Ho": 13, "Er": 14, "Tm": 15, "Yb": 16, "Lu": 17,
}

# Same keys as STD_CAP;
VALENCE_TARGET = STD_CAP

ALKALI_HEAVY_IONIC = frozenset({"Na", "K", "Rb", "Cs", "Fr"})
ALKALINE_EARTH_IONIC = frozenset({"Mg", "Ca", "Sr", "Ba", "Ra"})

# Pauling electronegativity
ENEG = {
    "F": 3.98,
    "O": 3.44,
    "Cl": 3.16,
    "N": 3.04,
    "Br": 2.96,
    "I": 2.66,
    "S": 2.58,
    "Se": 2.55,
    "C": 2.55,
    "P": 2.19,
    "H": 2.20,
    "As": 2.18,
    "Te": 2.10,
    "Si": 1.90,
    "B": 2.04,
    "Ge": 2.01,
    "Sn": 1.96,
    "Sb": 2.05,
    "Pb": 2.33,
    "Al": 1.61,
    "Ga": 1.81,
    "In": 1.78,
    "Tl": 2.04,
}

# Neutral valence electron counts for formal charge.
VALENCE_ELECTRONS = dict(VALENCE)

# Halides at M: always forced to covalent M–X single bond (MLX X-type).
TM_MONATOMIC_COV_LIGANDS = frozenset({"F", "Cl", "Br", "I"})

def _load_ccdc_covalent_radii(path: str = _CCDC_COV_RADII_JSON) -> dict[str, float]:
    """CCDC ChemistryLib Element.covalent_radius() values (Å) from ccdc_covalent_radii.json."""
    try:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
    except OSError:
        return {}
    by_sym = data.get("by_symbol")
    if not isinstance(by_sym, dict):
        return {}
    out: dict[str, float] = {}
    for sym, val in by_sym.items():
        try:
            out[str(sym)] = float(val)
        except (TypeError, ValueError):
            continue
    return out

COV_R_CCDC = _load_ccdc_covalent_radii()
CORE_E = {k: 0 for k in COV_R_CCDC}

def _load_tm_nonmetal_bond_limits(path: str = _TM_NONMETAL_LIMITS_JSON) -> dict[str, float]:
    """TM–ligand cutoff (Å): tmQM p99_A + margin, or Mercury GUI limit_A (no margin)."""
    try:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
    except OSError:
        return {}
    stats = data.get("stats")
    if not isinstance(stats, dict):
        return {}
    out: dict[str, float] = {}
    margin = TM_NONMETAL_P99_MARGIN_A
    for pair, rec in stats.items():
        if not isinstance(rec, dict):
            continue
        if rec.get("source") == "mercury_gui":
            lim = rec.get("limit_A")
            if lim is None:
                continue
            try:
                out[str(pair)] = float(lim)
            except (TypeError, ValueError):
                continue
            continue
        p99 = rec.get("p99_A")
        if p99 is None:
            continue
        try:
            out[str(pair)] = float(p99) + margin
        except (TypeError, ValueError):
            continue
    return out


TM_NONMETAL_BOND_LIMITS = _load_tm_nonmetal_bond_limits()


def _cov_radius(sym: str) -> float | None:
    return COV_R_CCDC.get(sym)


def _has_cov_radius(sym: str) -> bool:
    return _cov_radius(sym) is not None


def is_ionlike_s_block_metal(sym: str) -> bool:
    return sym in ALKALI_HEAVY_IONIC or sym in ALKALINE_EARTH_IONIC


def validate_atom_symbols(atom_syms: list[str]) -> None:
    """Raise ValueError if any symbol is missing from VALENCE or COV_R_CCDC."""
    missing_val = sorted({s for s in atom_syms if s not in VALENCE})
    missing_cov = sorted({s for s in atom_syms if not _has_cov_radius(s)})
    if not missing_val and not missing_cov:
        return
    lines: list[str] = []
    for sym in missing_val:
        idxs = ", ".join(str(i + 1) for i, s in enumerate(atom_syms) if s == sym)
        lines.append(f"{sym} (atom index: {idxs}): not in VALENCE")
    for sym in missing_cov:
        if sym in missing_val:
            continue
        idxs = ", ".join(str(i + 1) for i, s in enumerate(atom_syms) if s == sym)
        lines.append(f"{sym} (atom index: {idxs}): no covalent radius in COV_R_CCDC")
    raise ValueError(
        "Unsupported element symbol(s) in XYZ:\n" + "\n".join(f"  - {ln}" for ln in lines)
    )


def _looks_like_xyz_atom_line(line: str) -> bool:
    """True if stripped line has element + three floats (common coord row)."""
    p = line.strip().split()
    if len(p) < 4:
        return False
    try:
        float(p[1])
        float(p[2])
        float(p[3])
    except ValueError:
        return False
    head = p[0]
    if not head:
        return False
    # reject pure numbers as element
    try:
        float(head)
        return False
    except ValueError:
        pass
    return True


def read_xyz(path: str):
    """
    XYZ reader — always uses the same layout (1-based line numbers):

      line 1: integer n (number of atoms)
      line 2: comment / title (ignored for coordinates; may be empty)
      lines 3 .. n+2: element + x y z

    Line 2 must exist. Coordinates always start at line 3, even for minimal files:
    put a placeholder comment on line 2 if you have nothing to say.
    """
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    if not lines:
        raise ValueError(f"empty XYZ file: {path}")

    try:
        n = int(lines[0].strip().split()[0])
    except (ValueError, IndexError) as exc:
        raise ValueError(
            f"{path}: line 1 must be a positive integer atom count, got {lines[0]!r}"
        ) from exc
    if n <= 0:
        raise ValueError(f"{path}: atom count must be positive, got {n}")
    if len(lines) < 2:
        raise ValueError(
            f"{path}: need line 2 (comment/title). Coordinates start at line 3."
        )
    if _looks_like_xyz_atom_line(lines[1]):
        raise ValueError(
            f"{path}: line 2 looks like a coordinate row, but this program always "
            f"treats line 2 as the comment line. Insert a title/comment on line 2 "
            f"and move coordinates to lines 3–{n + 2}."
        )
    need = 2 + n
    if len(lines) < need:
        raise ValueError(
            f"{path}: for n={n}, need at least {need} lines "
            f"(1 count + 1 comment + {n} coords starting at line 3); "
            f"got {len(lines)}"
        )

    atoms = []
    for k in range(n):
        file_line_no = 3 + k
        ln = lines[2 + k]
        p = ln.split()
        if len(p) < 4:
            raise ValueError(
                f"{path}: line {file_line_no}: expected element + x y z, got {ln!r}"
            )
        raw = p[0]
        sym = raw[0].upper() + raw[1:].lower() if len(raw) > 1 else raw.upper()
        try:
            x, y, z = float(p[1]), float(p[2]), float(p[3])
        except ValueError as exc:
            raise ValueError(
                f"{path}: line {file_line_no}: non-numeric coordinates in {ln!r}"
            ) from exc
        atoms.append((k + 1, sym, x, y, z))
    validate_atom_symbols([a[1] for a in atoms])
    return atoms


def dist(a, b) -> float:
    def _xyz(p):
        if len(p) >= 5:
            return p[2], p[3], p[4]
        return p[0], p[1], p[2]

    ax, ay, az = _xyz(a)
    bx, by, bz = _xyz(b)
    return math.sqrt((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2)


def is_tm(symbol: str) -> bool:
    return symbol in TM_SET


def is_TM(symbol: str) -> bool:
    """Compatibility alias for viewer bridge code."""
    return is_tm(symbol)


def _prune_h_to_closest_nonmetal_neighbor(atoms, edges):
    """
    If H is within cutoff of more than one non-metal, keep only the closest contact.
    TM–H edges are unchanged (not counted toward the one-nonmetal limit).
    """
    by_h: dict[int, list[tuple[tuple[int, int, str, str], float]]] = defaultdict(list)
    for i, j, ei, ej in edges:
        if ei != "H" and ej != "H":
            continue
        if ei == "H":
            h_idx, other_idx, other_sym = i, j, ej
        else:
            h_idx, other_idx, other_sym = j, i, ei
        if is_tm(other_sym):
            continue
        d = dist(atoms[h_idx - 1], atoms[other_idx - 1])
        by_h[h_idx].append(((i, j, ei, ej), d))

    drop = set()
    for h_idx, candidates in by_h.items():
        if len(candidates) <= 1:
            continue
        candidates.sort(key=lambda item: item[1])
        for edge, _d in candidates[1:]:
            drop.add(edge)

    if not drop:
        return edges
    return [e for e in edges if e not in drop]


def _bond_cutoff_cov(ei, ej):
    margin = (
        COV_BOND_MARGIN_S_BLOCK
        if ei in S_BLOCK_SYMS or ej in S_BLOCK_SYMS
        else COV_BOND_MARGIN
    )
    return _cov_radius(ei) + _cov_radius(ej) + margin


def _tm_nonmetal_bond_cutoff(metal: str, ligand: str) -> float | None:
    """P99+0.05 Å limit for TM–nonmetal, or None to fall back to covalent radii."""
    if is_tm(metal) and not is_tm(ligand):
        return TM_NONMETAL_BOND_LIMITS.get(f"{metal}-{ligand}")
    return None


def _bond_cutoff(ei, ej):
    """Distance cutoff (Å): TM–nonmetal uses tmQM P99+0.05; all other pairs use COV+margin."""
    tm_lim = _tm_nonmetal_bond_cutoff(ei, ej)
    if tm_lim is None:
        tm_lim = _tm_nonmetal_bond_cutoff(ej, ei)
    if tm_lim is not None:
        return tm_lim
    return _bond_cutoff_cov(ei, ej)


def _within_cov_bond_cutoff(d, ei, ej):
    if not _has_cov_radius(ei) or not _has_cov_radius(ej):
        return False
    return d < _bond_cutoff(ei, ej)


def connectivity(atoms, factor=None):
    """Raw connectivity: TM–nonmetal uses tmQM P99+0.05 Å; other pairs COV+0.45/0.40 Å."""
    _ = factor  # deprecated; kept for call-site compatibility
    edges = []
    n = len(atoms)
    for i in range(n):
        ai = atoms[i]
        for j in range(i + 1, n):
            aj = atoms[j]
            ei, ej = ai[1], aj[1]
            if not _has_cov_radius(ei) or not _has_cov_radius(ej):
                continue
            cutoff = _bond_cutoff(ei, ej)
            if dist(ai, aj) < cutoff:
                edges.append((ai[0], aj[0], ei, ej))
    return _prune_h_to_closest_nonmetal_neighbor(atoms, edges)


def check_octet_violations(
    atoms, bo, adj, lp=None, fc=None, *, metal_adjacency_edges=None, coords=None
):
    """
    Valence check using ILP electron accounting only:
    bond_e + 2×lp vs VALENCE[sym] − fc (no extra dative-LP term).
    """
    violations = []
    check_syms = {
        "C", "N", "O", "S", "P", "B", "Si", "Se", "Te", "F", "Cl", "Br", "I",
    }

    for i, sym in enumerate(atoms):
        if sym not in check_syms:
            continue
        bond_e = sum(bo.get((min(i, k), max(i, k)), 0) for k in adj[i])
        lp_cnt = lp.get(i, 0) if lp else 0
        fc_val = fc[i] if fc else 0
        total_e = bond_e + 2 * lp_cnt
        exp_val = VALENCE.get(sym, 4) - fc_val

        if total_e == exp_val:
            continue

        neutral_val = VALENCE.get(sym, 4)
        if fc_val != 0 and total_e == neutral_val:
            continue

        if sym in ("B", "Si", "P", "Ge", "Sn", "As", "Sb", "Se", "Te", "Pb", "Bi"):
            if bond_e >= 3:
                continue

        violations.append(
            (i, f"{sym}{i+1}", bond_e, lp_cnt, 0, lp_cnt, exp_val, fc_val)
        )
    return violations


def print_summary_and_choose_ilp(atom_syms, bonds, lp_out, mol_charge):
    """
    Lewis-engine.py print_summary_and_choose layout for ILP results.
    All valence is reported as Lewis (non-Lewis rows zero); NBO-style table only.
    """
    te = float(sum(VALENCE.get(s, 0) for s in atom_syms) - mol_charge)
    if te <= 0:
        te = 1.0
    le, vnl, ryl = te, 0.0, 0.0
    stats_charge = float(mol_charge)
    print("          -------------------------------")
    print(f"                 Total Lewis{le:10.5f}  ({100 * le / te:8.4f}%)")
    print(f"           Valence non-Lewis{vnl:10.5f}  ({100 * vnl / te:8.4f}%)")
    print(f"           Rydberg non-Lewis{ryl:10.5f}  ({100 * ryl / te:8.4f}%)")
    print("          -------------------------------")
    print(f"               Total unit  1{te:10.5f}  ({100.0:8.4f}%)")
    print(f"              Charge unit  1{stats_charge:10.5f}")
    print()
    print_choose_block(bonds, lp_out)


def print_octet_report(
    atom_syms, bo0, lp_full, fc0, *, metal_adjacency_edges=None, coords=None
):
    """Lewis-engine.py-style octet / valence block after $CHOOSE."""
    adj_check = defaultdict(list)
    for i, j in bo0:
        adj_check[i].append(j)
        adj_check[j].append(i)
    violations = check_octet_violations(
        atom_syms,
        bo0,
        adj_check,
        lp=lp_full,
        fc=fc0,
        metal_adjacency_edges=metal_adjacency_edges,
        coords=coords,
    )
    if violations:
        print()
        print("=" * 72)
        print("  ⚠  Octet / Valence Violations in Best Lewis Structure")
        print("=" * 72)
        print(
            "  Atoms where bond_e + 2×LP ≠ expected valence electrons (ILP lp only)."
        )
        print("  These indicate a suboptimal Kekulé pattern or unusual geometry.")
        print()
        print(
            f"  {'Atom':<10} {'Bond-e':>6}  {'LP':>4}  {'Total-e':>7}  {'Exp-e':>6}  {'FC':>4}  Notes"
        )
        print(f"  {'-'*68}")
        for idx, label, be, lp_free, lp_dat, lp_cnt, exp_val, fc_val in violations:
            total_e = be + 2 * lp_cnt
            adj_syms = ", ".join(f"{atom_syms[k]}{k + 1}" for k in adj_check[idx])
            lp_str = f"{lp_cnt}" if lp_dat == 0 else f"{lp_cnt}({lp_free}+{lp_dat}d)"
            print(
                f"  {label:<10} {be:>6}  {lp_str:>4}  {total_e:>7}  {exp_val:>6}  {fc_val:>+4}  bonded to: {adj_syms}"
            )
        print()
        print(f"  → {len(violations)} atom(s) violate expected valence in this Lewis structure.")
    else:
        print()
        print("  ✓  All atoms satisfy their expected valence (octet rule OK).")

    nonzero_fc = [
        (i, f"{atom_syms[i]}{i + 1}", int(fc0[i]))
        for i in range(len(atom_syms))
        if i < len(fc0) and fc0[i] != 0 and not is_TM(atom_syms[i])
    ]
    if nonzero_fc:
        print()
        print("  Atoms with non-zero formal charge:")
        for _i, label, fc_val in nonzero_fc:
            print(f"    {label}:  fc = {fc_val:+d}")
    else:
        print()
        print("  No non-zero formal charges on non-metal atoms.")


def print_choose_block(bonds, lp_out):
    print(" $CHOOSE")
    lone_items = sorted((idx, v) for idx, v in lp_out.items() if v > 0)
    if lone_items:
        print("   LONE " + " ".join(f"{idx} {v}" for idx, v in lone_items) + " END")
    line = "   BOND"
    for i, j, bo in bonds:
        if bo == 1:
            tok = f"S {i} {j}"
        elif bo == 2:
            tok = f"D {i} {j}"
        elif bo == 3:
            tok = f"T {i} {j}"
        else:
            tok = f"S {i} {j}"
        if len(line) + 1 + len(tok) > 72:
            print(line)
            line = "       " + tok
        else:
            line += " " + tok
    print(line + " END")
    print(" $END")


_ORGANIC = {"C", "H", "N", "O", "F", "Cl"}


def is_inorganic(sym):
    return sym not in _ORGANIC


_NORMAL_BONDS = {
    "H": 1,
    "B": 3, "C": 4, "N": 3, "O": 2, "F": 1,
    "Si": 4, "P": 3, "S": 2, "Cl": 1,
    "Ge": 4, "As": 3, "Se": 2, "Br": 1,
    "Sn": 4, "Sb": 3, "Te": 2, "I": 1,
    "Pb": 4, "Bi": 3,
    "Al": 3, "Ga": 3, "In": 3, "Tl": 3,
}


def normal_bonds(sym):
    return _NORMAL_BONDS.get(sym, 4)


_OXSTATE_ENEG = dict(ENEG)
# Atoms eligible for CBC η/π cluster BFS (phospholyl, aza-Cp, etc.).
PI_CLUSTER_SYMS = frozenset({"C", "N", "P"})

_OXSTATE_ENEG.update({
    "Sc": 1.36, "Ti": 1.54, "V": 1.63, "Cr": 1.66, "Mn": 1.55, "Fe": 1.83, "Co": 1.88,
    "Ni": 1.91, "Cu": 1.90, "Zn": 1.65, "Y": 1.22, "Zr": 1.33, "Nb": 1.60, "Mo": 2.16,
    "Tc": 1.90, "Ru": 2.20, "Rh": 2.28, "Pd": 2.20, "Ag": 1.93, "Cd": 1.69,
    "Hf": 1.30, "Ta": 1.50, "W": 2.36, "Re": 1.90, "Os": 2.20, "Ir": 2.20,
    "Pt": 2.28, "Au": 2.54, "Hg": 2.00,
})

# Minimum Lewis bond order inside an η fragment to count as one CBC L (η²) pair.
ILP_ETA_PI_BOND_MIN_ORDER = 2


def _ligand_component_by_array_index(atom_syms, bo, metal_adjacency_edges=None):
    """Array index → ligand fragment id (non-TM subgraph without counting TM–L as connectivity)."""
    n = len(atom_syms)
    pseudo_atoms = [(i + 1, atom_syms[i]) for i in range(n)]
    seen = set()
    edges_id = []
    for (i, j), order in bo.items():
        if order <= 0 or is_TM(atom_syms[i]) or is_TM(atom_syms[j]):
            continue
        key = (min(i, j), max(i, j))
        if key in seen:
            continue
        seen.add(key)
        edges_id.append((i + 1, j + 1, atom_syms[i], atom_syms[j]))
    if metal_adjacency_edges:
        for i, j, ei, ej in metal_adjacency_edges:
            if is_TM(ei) ^ is_TM(ej):
                continue
            if is_TM(atom_syms[i]) or is_TM(atom_syms[j]):
                continue
            key = (min(i, j), max(i, j))
            if key in seen:
                continue
            seen.add(key)
            edges_id.append((i + 1, j + 1, atom_syms[i], atom_syms[j]))
    by_atom_id = _nonmetal_ligand_components(pseudo_atoms, edges_id)
    return {i: by_atom_id.get(i + 1) for i in range(n)}


def _ligand_adjacency_from_edges(atoms_packed, edges):
    """Non-TM ligand connectivity from raw edges (no M–L edges)."""
    adj = defaultdict(set)
    for i, j, ei, ej in edges:
        if _is_tm_nm_edge(ei, ej):
            continue
        if is_TM(ei) or is_TM(ej):
            continue
        # η/π skeleton adjacency: allow C/N/P connectivity (e.g. phospholyl, N–P–N),
        # but keep other inorganic-heavy edges out of the ligand graph.
        if (is_inorganic(ei) or is_inorganic(ej)) and not (ei in PI_CLUSTER_SYMS and ej in PI_CLUSTER_SYMS):
            continue
        adj[i].add(j)
        adj[j].add(i)
    return adj


def _metal_ligand_contacts_from_edges(atoms_packed, edges):
    """(metal_atom_id, ligand_atom_id) for each inorganic–ligand contact in *edges*."""
    pairs = []
    for i, j, ei, ej in edges:
        if is_TM(ei) and not is_TM(ej):
            pairs.append((i, j))
        elif is_TM(ej) and not is_TM(ei):
            pairs.append((j, i))
        elif is_inorganic(ei) and not is_TM(ej) and not is_inorganic(ej):
            pairs.append((i, j))
        elif is_inorganic(ej) and not is_TM(ei) and not is_inorganic(ei):
            pairs.append((j, i))
    return pairs


def _eta_coordinating_groups_atom_ids(
    metal_id, coord_ids, lig_comp, lig_adj, *, min_group_size=1
):
    """
    Same as _eta_coordinating_groups but atom ids (XYZ indices); drops groups
    smaller than *min_group_size* when that is > 1.
    """
    by_ligand = defaultdict(set)
    for i in coord_ids:
        lid = lig_comp.get(i)
        if lid is None:
            continue
        by_ligand[lid].add(i)
    groups = []
    for lid, coord_set in by_ligand.items():
        ligand_atoms = {i for i, l in lig_comp.items() if l == lid}
        visited_coord = set()
        for seed in sorted(coord_set):
            if seed in visited_coord:
                continue
            queue = [seed]
            seen = {seed}
            comp_coord = set()
            while queue:
                u = queue.pop()
                if u in coord_set:
                    comp_coord.add(u)
                # η candidates require CONTIGUOUS coordinating atoms on the ligand
                # skeleton: traversal is restricted to the coordinating-atom subgraph.
                for v in lig_adj[u]:
                    if v in coord_set and v not in seen:
                        seen.add(v)
                        queue.append(v)
            visited_coord |= comp_coord
            if len(comp_coord) >= min_group_size:
                groups.append(sorted(comp_coord))
    if min_group_size <= 1:
        for i in sorted(coord_ids):
            if lig_comp.get(i) is None:
                groups.append([i])
    return groups


def _eta_carbon_atom_ids(
    atoms_packed,
    edges,
    *,
    min_group_size=ETA_MIN_COORDINATING_GROUP_SIZE,
):
    """Carbon atom ids on a geometric η fragment to a transition-metal centre (is_TM)."""
    atom_el = {i: el for i, el, *_ in atoms_packed}
    lig_comp = _nonmetal_ligand_components(atoms_packed, edges)
    lig_adj = _ligand_adjacency_from_edges(atoms_packed, edges)
    metals = {i for i, el, *_ in atoms_packed if is_TM(el)}
    eta_c = set()
    by_metal = defaultdict(set)
    for tm, lig in _metal_ligand_contacts_from_edges(atoms_packed, edges):
        if tm in metals:
            by_metal[tm].add(lig)
    for tm, coord in by_metal.items():
        for group in _eta_coordinating_groups_atom_ids(
            tm, coord, lig_comp, lig_adj, min_group_size=min_group_size
        ):
            for aid in group:
                if atom_el.get(aid) == "C":
                    eta_c.add(aid)
    return eta_c


def _eta_ligand_atoms_by_metal(
    atoms_packed,
    edges,
    *,
    min_group_size=ETA_MIN_COORDINATING_GROUP_SIZE,
):
    """
    Per transition-metal centre: ligand atom ids in geometric η fragments
    (≥ *min_group_size* contiguous TM-bound atoms on one ligand skeleton).
    """
    lig_comp = _nonmetal_ligand_components(atoms_packed, edges)
    lig_adj = _ligand_adjacency_from_edges(atoms_packed, edges)
    metals = {i for i, el, *_ in atoms_packed if is_TM(el)}
    out = {}
    by_metal = defaultdict(set)
    for tm, lig in _metal_ligand_contacts_from_edges(atoms_packed, edges):
        if tm in metals:
            by_metal[tm].add(lig)
    for tm, coord in by_metal.items():
        eta_atoms = set()
        for group in _eta_coordinating_groups_atom_ids(
            tm, coord, lig_comp, lig_adj, min_group_size=min_group_size
        ):
            eta_atoms.update(group)
        if eta_atoms:
            out[tm] = eta_atoms
    return out


def _eta_coordinating_groups(metal_idx, coord_indices, lig_comp, adj_lewis):
    """
    Within each ligand fragment, partition M-coordinating atoms into connected
    groups via the ligand skeleton (paths may use non-coordinating atoms).
    """
    by_ligand = defaultdict(set)
    for i in coord_indices:
        lid = lig_comp.get(i)
        if lid is None:
            continue
        by_ligand[lid].add(i)
    groups = []
    for lid, coord_set in by_ligand.items():
        ligand_atoms = {i for i, l in lig_comp.items() if l == lid}
        visited_coord = set()
        for seed in sorted(coord_set):
            if seed in visited_coord:
                continue
            queue = [seed]
            seen = {seed}
            comp_coord = set()
            while queue:
                u = queue.pop()
                if u in coord_set:
                    comp_coord.add(u)
                # η candidates require CONTIGUOUS coordinating atoms on the ligand
                # skeleton: traversal is restricted to the coordinating-atom subgraph.
                for v in adj_lewis[u]:
                    if v in coord_set and v not in seen:
                        seen.add(v)
                        queue.append(v)
            visited_coord |= comp_coord
            if comp_coord:
                groups.append(sorted(comp_coord))
    for i in sorted(coord_indices):
        if lig_comp.get(i) is None:
            groups.append([i])
    return groups


def _ilp_cbc_records_for_eta_group(
    metal_idx,
    group,
    bo_ij,
    adj_lewis,
    lp_lewis=None,
    *,
    pi_min_order=None,
):
    """
    η fragment: ILP bo(M–L)>0 → X; dative ends with bo(L–L')≥pi_min_order → one L pair;
    remaining dative coordinators → single-atom L only if lp > 0 (non-η dative).
    """
    if pi_min_order is None:
        pi_min_order = ILP_ETA_PI_BOND_MIN_ORDER
    records = []
    dative = set()
    for i in group:
        if bo_ij(i, metal_idx) > 0:
            records.append(((i,), "X"))
        else:
            dative.add(i)
    in_pair = set()
    for a in sorted(dative):
        for b in adj_lewis[a]:
            if b not in dative or b <= a:
                continue
            if bo_ij(a, b) >= pi_min_order:
                records.append(((a, b), "L"))
                in_pair.add(a)
                in_pair.add(b)
    for i in sorted(dative):
        if i not in in_pair:
            lp_cnt = lp_lewis.get(i, 0) if lp_lewis else 0
            if lp_cnt > 0:
                records.append(((i,), "L"))
    return records


def _eta_group_internal_lp_edge_keys(
    atoms_packed,
    edges,
    *,
    min_group_size=ETA_MIN_COORDINATING_GROUP_SIZE,
):
    """
    Ligand-skeleton edge keys (min, max) with both ends in the same geometric η
    group (≥ *min_group_size* TM-bound atoms) for any transition metal centre.
    """
    lig_comp = _nonmetal_ligand_components(atoms_packed, edges)
    lig_adj = _ligand_adjacency_from_edges(atoms_packed, edges)
    metals = {i for i, el, *_ in atoms_packed if is_TM(el)}
    keys = set()
    for tm in metals:
        coord = {
            lig
            for t, lig in _metal_ligand_contacts_from_edges(atoms_packed, edges)
            if t == tm
        }
        for group in _eta_coordinating_groups_atom_ids(
            tm, coord, lig_comp, lig_adj, min_group_size=min_group_size
        ):
            gset = set(group)
            for a in group:
                for b in lig_adj.get(a, ()):
                    if b not in gset or a >= b:
                        continue
                    keys.add((a, b))
    return sorted(keys)


_SIGMA_AGOSTIC_H_PARENTS = frozenset({"B", "C", "Si", "Al", "Ga"})
_DIHYDROGEN_HH_MAX = 1.15


def _cbc_record_for_h_neighbor(
    metal_idx,
    h_idx,
    atoms,
    coords,
    adj_lewis,
    adj_full,
    bo_ij,
    *,
    seen_hh_pairs=None,
):
    """
    CBC for H in the metal coordination sphere (aligned with Lewis-engine.py):

    M–H in Lewis (bo>0) → X; agostic X–H (X = B,C,…) with b_tm=0 → (X,H) L;
    η²-H₂ → (H,H) L; protic H → X; geometric M···H only → None (skip).
    """
    h_nbrs = list(adj_lewis[h_idx])
    if bo_ij(metal_idx, h_idx) > 0 or any(is_TM(atoms[k]) for k in h_nbrs):
        return ((h_idx,), "X")

    if len(h_nbrs) == 0:
        for other_h in sorted(adj_full[metal_idx]):
            if atoms[other_h] != "H" or other_h == h_idx:
                continue
            other_nbrs = list(adj_lewis[other_h])
            if any(is_TM(atoms[k]) for k in other_nbrs):
                continue
            non_tm = [k for k in other_nbrs if not is_TM(atoms[k])]
            if len(non_tm) > 1:
                continue
            if len(non_tm) == 1 and non_tm[0] != h_idx:
                continue
            if dist(coords[h_idx], coords[other_h]) >= _DIHYDROGEN_HH_MAX:
                continue
            pair = tuple(sorted((h_idx, other_h)))
            if seen_hh_pairs is not None:
                if pair in seen_hh_pairs:
                    return None
                seen_hh_pairs.add(pair)
            return (pair, "L")
    elif (
        len(h_nbrs) == 1
        and atoms[h_nbrs[0]] == "H"
        and not any(is_TM(atoms[k]) for k in h_nbrs)
    ):
        other_h = h_nbrs[0]
        if other_h in adj_full[metal_idx]:
            pair = tuple(sorted((h_idx, other_h)))
            if seen_hh_pairs is not None:
                if pair in seen_hh_pairs:
                    return None
                seen_hh_pairs.add(pair)
            return (pair, "L")

    if len(h_nbrs) == 1 and atoms[h_nbrs[0]] in _SIGMA_AGOSTIC_H_PARENTS:
        return ((h_nbrs[0], h_idx), "L")

    if len(h_nbrs) == 1 and atoms[h_nbrs[0]] in ("N", "O", "S", "F", "Cl", "Br", "I"):
        return ((h_idx,), "X")

    if bo_ij(metal_idx, h_idx) == 0:
        return None
    return ((h_idx,), "X")


def classify_cbc_ligands(atoms, coords, bo, lp, fc, charge=0, *, metal_adjacency_edges=None):
    n = len(atoms)

    def _full_connectivity_dative(atoms, coords, bo_lewis, lp_lewis, dative_ml_pairs):
        # ILP lp[] is non-donated LPs only; each dative M–L still has one pair on the bond.
        dative_lp_by_lig = defaultdict(int)
        for _tm, lig in dative_ml_pairs:
            dative_lp_by_lig[lig] += 1

        def _effective_lp(lig_idx):
            free = lp_lewis.get(lig_idx, 0) if lp_lewis else 0
            return free + dative_lp_by_lig.get(lig_idx, 0)

        lewis_adj = defaultdict(set)
        for (a, b) in bo_lewis:
            lewis_adj[a].add(b)
            lewis_adj[b].add(a)

        tm_sigma_partners = defaultdict(set)
        for a, b in bo_lewis:
            sa, sb = atoms[a], atoms[b]
            if is_TM(sa) and not is_TM(sb):
                tm_sigma_partners[a].add(b)
            elif is_TM(sb) and not is_TM(sa):
                tm_sigma_partners[b].add(a)

        def is_terminal_co_o(lig_idx, tm_idx):
            sym = atoms[lig_idx]
            if sym not in ("O", "N"):
                return False
            nbrs = list(lewis_adj[lig_idx])
            if len(nbrs) != 1:
                return False
            c_nbr = nbrs[0]
            if atoms[c_nbr] != "C":
                return False
            return tm_idx in lewis_adj[c_nbr]

        def _donor_lp_on_c(c_idx, tm_idx):
            """LP on C: free lone pairs and/or one pair on a dative M–C bond (ILP lp may be 0)."""
            lp_free = lp_lewis.get(c_idx, 0)
            if lp_free > 0:
                return lp_free
            if (tm_idx, c_idx) in ml_bo_zero:
                return 1
            return 0

        def is_lp_carbene_c_to_tm(c_idx, tm_idx):
            """C coordinating to TM (σ, dative, or lp); substituents on this C are backbone."""
            if atoms[c_idx] != "C":
                return False
            if not _within_cov_bond_cutoff(
                dist(coords[c_idx], coords[tm_idx]), atoms[tm_idx], "C"
            ):
                return False
            if c_idx in tm_sigma_partners[tm_idx]:
                return True
            key = (min(c_idx, tm_idx), max(c_idx, tm_idx))
            if bo_lewis.get(key, 0) > 0:
                return True
            return _donor_lp_on_c(c_idx, tm_idx) > 0

        def is_backbone_atom(lig_idx, tm_idx):
            for nbr in lewis_adj[lig_idx]:
                if nbr in tm_sigma_partners[tm_idx]:
                    return True
                if is_lp_carbene_c_to_tm(nbr, tm_idx):
                    return True
            return False

        def is_saturated_no_lp(lig_idx):
            sym = atoms[lig_idx]
            if sym == "C" or is_TM(sym) or sym == "H":
                return False
            if _effective_lp(lig_idx) > 0:
                return False
            cap = STD_CAP.get(sym, 4)
            bond_e = sum(
                bo_lewis.get((min(lig_idx, k), max(lig_idx, k)), 0)
                for k in lewis_adj[lig_idx]
            )
            return bond_e >= cap

        bonds = []
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                si, sj = atoms[i], atoms[j]
                tm_i, tm_j = is_TM(si), is_TM(sj)
                d_ij = dist(coords[i], coords[j])

                if (tm_i and not tm_j) or (tm_j and not tm_i):
                    tm = i if tm_i else j
                    lig = j if tm_i else i
                    sym_tm, sym_lig = atoms[tm], atoms[lig]
                    key_ml = (min(tm, lig), max(tm, lig))

                    if is_terminal_co_o(lig, tm):
                        continue
                    # Keep explicit Lewis M–L σ bonds even if the ligand atom is saturated
                    # and has no lone pairs (e.g. M–SiR3). The saturated-no-LP filter is
                    # only meant to suppress *purely geometric* M···L contacts.
                    if bo_lewis.get(key_ml, 0) <= 0 and is_saturated_no_lp(lig):
                        continue
                    if not _within_cov_bond_cutoff(d_ij, sym_tm, sym_lig):
                        continue
                    bonds.append((i, j))
                elif _within_cov_bond_cutoff(d_ij, si, sj):
                    bonds.append((i, j))

        # Geometric M–L with b_tm=0: keep C (η/haptic graph); other elements need lp > 0.
        for tm, lig in ml_bo_zero:
            if atoms[lig] != "C" and _effective_lp(lig) <= 0:
                continue
            bonds.append((tm, lig))

        tm_indices = [i for i in range(len(atoms)) if is_TM(atoms[i])]
        bond_set = set(map(tuple, [sorted(b) for b in bonds]))
        adj_tmp = defaultdict(set)
        for i, j in bonds:
            adj_tmp[i].add(j)
            adj_tmp[j].add(i)

        all_C_adj = defaultdict(set)
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                if atoms[i] == "C" and atoms[j] == "C":
                    if _within_cov_bond_cutoff(dist(coords[i], coords[j]), "C", "C"):
                        all_C_adj[i].add(j)
                        all_C_adj[j].add(i)

        for tm in tm_indices:
            tm_C_all = {j for j in adj_tmp[tm] if atoms[j] == "C"}
            tm_C_dative = {c for c in tm_C_all if bo_lewis.get((min(tm, c), max(tm, c)), 0) == 0}
            if not tm_C_dative:
                continue
            visited_c = set(tm_C_all)
            frontier = set(tm_C_dative)
            for _ in range(3):
                new_frontier = set()
                for c in frontier:
                    for nb in all_C_adj[c]:
                        if nb not in visited_c:
                            d_nb = dist(coords[tm], coords[nb])
                            if _within_cov_bond_cutoff(
                                d_nb, atoms[tm], atoms[nb]
                            ):
                                key = (min(tm, nb), max(tm, nb))
                                if key not in bond_set:
                                    bonds.append((tm, nb))
                                    bond_set.add(key)
                                    adj_tmp[tm].add(nb)
                                    adj_tmp[nb].add(tm)
                                visited_c.add(nb)
                                new_frontier.add(nb)
                frontier = new_frontier

        return bonds

    # M–L with no Lewis bond (b_tm=0): used for η expansion / coordination graph only.
    ml_bo_zero = set()
    if metal_adjacency_edges:
        for i, j, ei, ej in metal_adjacency_edges:
            if not is_TM(ei) ^ is_TM(ej):
                continue
            tm_i, lig_i = (i, j) if is_TM(ei) else (j, i)
            key = (min(tm_i, lig_i), max(tm_i, lig_i))
            if bo.get(key, 0) == 0:
                ml_bo_zero.add((tm_i, lig_i))
    else:
        for i in range(n):
            if not is_TM(atoms[i]):
                continue
            for j in range(n):
                if i == j or is_TM(atoms[j]):
                    continue
                key = (min(i, j), max(i, j))
                if bo.get(key, 0) != 0:
                    continue
                if _within_cov_bond_cutoff(dist(coords[i], coords[j]), atoms[i], atoms[j]):
                    ml_bo_zero.add((i, j))

    lp = dict(lp) if lp else {}
    all_bonds_list = _full_connectivity_dative(atoms, coords, bo, lp, ml_bo_zero)
    adj_full = defaultdict(set)
    for i, j in all_bonds_list:
        adj_full[i].add(j)
        adj_full[j].add(i)

    adj_lewis = defaultdict(set)
    for (i, j) in bo:
        adj_lewis[i].add(j)
        adj_lewis[j].add(i)

    def bo_ij(i, j):
        return bo.get((min(i, j), max(i, j)), 0)

    def _neighbor_primary_from_records(atom_idx, records):
        for tup, t in records:
            if atom_idx in tup and t == "X":
                return "X"
        for tup, t in records:
            if atom_idx in tup:
                return t
        return None

    def _fix_symmetric_chelate_O_records(metal_idx, records):
        """Promote one chelating O from L to X when two O donors share a π-linked C backbone."""
        nbr_cls = {}
        for tup, t in records:
            for a in tup:
                nbr_cls[a] = [t]
        o_L = [
            idx for idx, types in nbr_cls.items()
            if atoms[idx] == "O" and types == ["L"] and bo_ij(idx, metal_idx) == 0
        ]
        if len(o_L) < 2:
            return records
        for i in range(len(o_L)):
            for j in range(i + 1, len(o_L)):
                oi, oj = o_L[i], o_L[j]
                ci = next(
                    (k for k in adj_lewis[oi] if atoms[k] == "C" and bo_ij(oi, k) >= 2),
                    None,
                )
                cj = next(
                    (k for k in adj_lewis[oj] if atoms[k] == "C" and bo_ij(oj, k) >= 2),
                    None,
                )
                if ci is None or cj is None:
                    continue
                connected = ci == cj or cj in adj_lewis[ci] or ci in adj_lewis[cj]
                if not connected:
                    for mid in adj_lewis[ci]:
                        if cj in adj_lewis[mid]:
                            connected = True
                            break
                if connected:
                    flip = max(oi, oj)
                    out = []
                    for tup, t in records:
                        if tup == (flip,) and t == "L":
                            out.append((tup, "X"))
                        else:
                            out.append((tup, t))
                    return out
        return records

    lig_comp = _ligand_component_by_array_index(atoms, bo, metal_adjacency_edges)

    results = {}
    neighbor_cbc: dict[int, dict[int, str]] = {}
    for metal_idx in range(n):
        sym_m = atoms[metal_idx]
        if not is_inorganic(sym_m):
            continue
        if sym_m in ("Br", "I"):
            n_sub_m = sum(bo_ij(metal_idx, k) for k in adj_lewis[metal_idx])
            tm_nbrs = [j for j in adj_full[metal_idx] if is_TM(atoms[j])]
            if n_sub_m == 0 and tm_nbrs:
                continue
        neighbours = sorted(adj_full[metal_idx])
        if not neighbours:
            continue

        interaction_records = []
        coord_for_eta = set()
        seen_hh_pairs = set()
        for nbr_idx in neighbours:
            sym_n = atoms[nbr_idx]
            if is_TM(sym_n):
                en_m = _OXSTATE_ENEG.get(sym_m, 2.0)
                en_n = _OXSTATE_ENEG.get(sym_n, 2.0)
                interaction_records.append(((nbr_idx,), "L" if en_n < en_m else "Z"))
                continue
            if sym_n == "H":
                h_rec = _cbc_record_for_h_neighbor(
                    metal_idx,
                    nbr_idx,
                    atoms,
                    coords,
                    adj_lewis,
                    adj_full,
                    bo_ij,
                    seen_hh_pairs=seen_hh_pairs,
                )
                if h_rec is not None:
                    interaction_records.append(h_rec)
                continue
            coord_for_eta.add(nbr_idx)

        for group in _eta_coordinating_groups(
            metal_idx, coord_for_eta, lig_comp, adj_lewis
        ):
            interaction_records.extend(
                _ilp_cbc_records_for_eta_group(
                    metal_idx, group, bo_ij, adj_lewis, lp
                )
            )

        interaction_records = _fix_symmetric_chelate_O_records(
            metal_idx, interaction_records
        )
        if not interaction_records:
            continue

        neighbor_cbc[metal_idx] = {
            nbr_idx: _neighbor_primary_from_records(nbr_idx, interaction_records)
            for nbr_idx in neighbours
            if _neighbor_primary_from_records(nbr_idx, interaction_records) is not None
        }
        results[metal_idx] = interaction_records

    return results, neighbor_cbc


def cbc_x_ligands_from_interaction_records(cbc_interaction_records):
    """
    Per-TM set of ligand indices that appear in a CBC record classified **X**
    (same records as ``print_cbc_report``).
    """
    out = {}
    for metal_idx, records in (cbc_interaction_records or {}).items():
        xs = set()
        for atoms_t, typ in records:
            if typ == "X":
                xs.update(atoms_t)
        if xs:
            out[metal_idx] = xs
    return out


def print_cbc_report(
    atoms,
    coords,
    bo,
    lp,
    fc,
    charge=0,
    *,
    metal_adjacency_edges=None,
    cbc_bundle=None,
):
    """Print the CBC classification table for every inorganic centre."""
    if cbc_bundle is None:
        cbc_bundle = classify_cbc_ligands(
            atoms, coords, bo, lp, fc, charge, metal_adjacency_edges=metal_adjacency_edges
        )
    results, _neighbor_cbc = cbc_bundle

    if not results:
        print("\n  (No inorganic atoms found — no ligand classification.)")
        return

    print()
    print("=" * 72)
    print("  Ligand Classification  (L = dative, X = covalent, Z = Lewis-acid)")
    print("=" * 72)

    adj_lewis = defaultdict(set)
    for (i, j) in bo:
        adj_lewis[i].add(j)
        adj_lewis[j].add(i)

    def bo_ij(i, j):
        return bo.get((min(i, j), max(i, j)), 0)

    def avg_dist(metal_idx, atom_tuple):
        if len(atom_tuple) == 1:
            a = atom_tuple[0]
            return math.sqrt(sum((coords[metal_idx][k] - coords[a][k]) ** 2 for k in range(3)))
        a, b = atom_tuple
        if atoms[b] == "H" and atoms[a] in ("B", "C", "Si", "Al", "Ga"):
            return math.sqrt(sum((coords[metal_idx][k] - coords[b][k]) ** 2 for k in range(3)))
        dists = [
            math.sqrt(sum((coords[metal_idx][k] - coords[x][k]) ** 2 for k in range(3)))
            for x in atom_tuple
        ]
        return sum(dists) / len(dists)

    def row_label(atom_tuple):
        parts = [f"{atoms[a]}{a+1}" for a in atom_tuple]
        return "-".join(parts)

    def row_note(atom_tuple, cbc_char, metal_idx):
        if cbc_char == "Z":
            return "Lewis acid — metal donates electrons to ligand"
        rep = atom_tuple[0]
        sym = atoms[rep]
        n_val = _NORMAL_BONDS.get(sym, 4)
        n_sub = sum(bo_ij(rep, k) for k in adj_lewis[rep] if not is_inorganic(atoms[k]))
        lp_cnt = lp.get(rep, 0) if lp else 0

        if len(atom_tuple) == 2:
            a, b = atom_tuple
            if atoms[a] == "H" and atoms[b] == "H":
                return "η²-H₂ dihydrogen — σ(H-H) dative donor"
            if atoms[b] == "H" and atoms[a] in ("B", "C", "Si", "Al", "Ga"):
                return f"agostic {atoms[a]}-H σ-bond → σ-complex dative donation to M"
            bond_type = "C=C" if (atoms[a] == "C" and atoms[b] == "C") else (
                "N=N" if (atoms[a] == "N" and atoms[b] == "N") else "π"
            )
            return f"η² {bond_type} π-bond → dative π donation to M"

        if sym == "H":
            h_lewis_nbrs = list(adj_lewis[rep])
            if cbc_char == "X":
                if h_lewis_nbrs:
                    return "hydride H — covalent M-H bond"
                return "hydride H — covalent M-H bond"
            if cbc_char == "L":
                if len(h_lewis_nbrs) == 1:
                    parent = h_lewis_nbrs[0]
                    parent_sym = atoms[parent]
                    return f"agostic {parent_sym}-H → σ-complex dative donation to M"
                for other_h in adj_lewis[rep]:
                    if atoms[other_h] == "H":
                        return "η²-H₂ dihydrogen — σ(H-H) dative donor"
                return "H σ-complex — dative donor"

        if cbc_char == "L":
            is_co = (sym == "C" and any(atoms[k] == "O" and bo_ij(rep, k) >= 2 for k in adj_lewis[rep]))
            if is_co:
                return "CO — σ-donor via C lone pair"
            if sym in ("C", "N") and n_sub <= 2:
                return f"{sym} carbene/NHC — lone-pair σ-donor"
            if sym == "O" and n_sub < n_val and lp_cnt > 0:
                fc_rep = fc[rep] if fc and rep < len(fc) else 0
                if fc_rep < 0:
                    return f"O⁻ (oxyanion) — lone-pair dative donor ({lp_cnt} LP)"
                return f"O=N/S (resonance) — lone-pair dative donor ({lp_cnt} LP)"
            if lp_cnt > 0:
                return f"{sym}: {n_sub}/{n_val} bonds → lone-pair donor ({lp_cnt} LP)"
            return f"{sym}: {n_sub}/{n_val} bonds → dative donor"

        if cbc_char == "X":
            n_sub_total = sum(bo_ij(rep, k) for k in adj_lewis[rep] if not is_TM(atoms[k]))
            deficit_real = n_val - n_sub_total
            if deficit_real <= 0:
                return f"{sym}: {n_sub_total}/{n_val} bonds + covalent σ bond to M"
            return f"{sym}: {n_sub_total}/{n_val} bonds → {deficit_real} short → covalent σ to M"
        return ""

    for metal_idx in sorted(results):
        sym_m = atoms[metal_idx]
        interaction_records = results[metal_idx]

        total_L = sum(1 for _, t in interaction_records if t == "L")
        total_X = sum(1 for _, t in interaction_records if t == "X")
        total_Z = sum(1 for _, t in interaction_records if t == "Z")

        parts = []
        if total_L:
            parts.append(f"L{_subscript(total_L)}")
        if total_X:
            parts.append(f"X{_subscript(total_X)}")
        if total_Z:
            parts.append(f"Z{_subscript(total_Z)}")
        designation = "".join(parts) if parts else "—"

        print(f"\n  {sym_m}{metal_idx+1}   [{designation}]")
        print(f"  {'Neighbour':<16} {'Avg.Dist':>9}  {'Type':<5}  Notes")
        print(f"  {'-'*65}")

        seen_multi = {}
        collapsed = []
        for atom_tuple, cbc_char in interaction_records:
            key = (atom_tuple, cbc_char)
            if key in seen_multi:
                seen_multi[key] += 1
            else:
                seen_multi[key] = 1
                collapsed.append((atom_tuple, cbc_char))

        for atom_tuple, cbc_char in collapsed:
            key = (atom_tuple, cbc_char)
            count = seen_multi[key]
            lbl = row_label(atom_tuple)
            if count > 1:
                lbl = f"{lbl}(×{count})"
            d_str = f"{avg_dist(metal_idx, atom_tuple):.2f} Å"
            note = row_note(atom_tuple, cbc_char, metal_idx)
            if count > 1:
                sym_rep = atoms[atom_tuple[0]]
                bond_names = {2: "double (2×)", 3: "triple (3×)"}
                bond_desc = bond_names.get(count, f"{count}×")
                note = f"{sym_rep}: {bond_desc} covalent M={sym_rep} bond (nitrido/oxo/imido)"
            print(f"  {lbl:<16} {d_str:>9}  {cbc_char:<5}  {note}")

        print(f"\n  Total:  {total_L}L + {total_X}X + {total_Z}Z  →  {designation}")


def _subscript(n):
    subs = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    return str(n).translate(subs) if n > 1 else ""


# ----- ilp_bond_order_aromatic_workflow_V2.py -----
# Aromatic-aware ILP workflow V2: TM–nonmetal bonds are optimized inside the ILP
#
# Key ideas:
# 1) Build raw connectivity from XYZ.
# 2) Remove metal-nonmetal edges for ring/fused-ring detection.
# 3) Keep only planar ring/fused-ring systems as aromatic candidates.
# 4) Soft (optional): aromatic 6π via ILP_WEIGHT_AROMATIC_DEVIATION (>0 enables block).
# 4b) Soft: similar M–L distances on a ligand → same z_cov (ILP_WEIGHT_ML_DISTANCE_CLASS).
# 4c) Hard (optional): η-fragment carbons → lp = 0 (ILP_HARD_ETA_CARBON_LP_ZERO).
# 4d) Hard (always): terminal CO (2-atom C+O fragment, M–L via C) → C≡O triple, M–C dative (b_tm=0).
# 5) For O/N/S/P in aromatic systems:
#    - if no double bond around that atom, one lone pair (2e) can contribute;
#    - if double-bonded, pi contribution comes from double bonds (2 per double).
# 6) For ring carbons not bonded to a transition metal: if unsaturated (any incident
#    multiple bond or formal charge <= -1), the same optional lone-pair pi term applies.
# 7) Hard (optional): remote C → lp = 0 (ILP_HARD_C_LP_ONLY_TM_NEIGHBORS).
# 8) Hard (optional): mol_charge = Σfc(ligands) + Σox(TM) − Σb_tm (ILP_HARD_MOL_CHARGE_BALANCE).
# 9) Hard (optional): Σox(TM) ≥ Σb_tm on non-η M–L only (ILP_HARD_OX_GE_SIGMA); F/Cl/Br/I ligands omitted from RHS.
# 10) Soft: minimize Σox(TM) (ILP_WEIGHT_TM_OX_MINIMIZE).
# See ILP_HARD_* / ILP_WEIGHT_* at top of file.
#
# ILP uses full XYZ connectivity (raw): every TM–nonmetal contact is either
# covalent order 1/2/3 (b_tm) or no Lewis M–L bond (b_tm=0), except monatomic
# Cl/Br/I/F/H — fixed covalent single. Octet uses lp + bond_sum only (no dative
# electron term in ILP). CBC/SMILES assign dative (L) when b_tm=0 and lp/fc/η rules apply.
# Ring detection strips M–L edges before cycle search.
# $CHOOSE omits b_tm=0 M–L; covalent M–L prints as S/D/T.
# """

base = sys.modules[__name__]

pulp = base.pulp


def _force_monatomic_ml_single_cov(lig_idx, atom_el, edges):
    """
    Monatomic ML single covalent rule: F/Cl/Br/I always; H only if isolated
    (no non-TM neighbor).

    F/Cl/Br/I: always. H: only if it has no non-TM neighbour in *edges* (isolated
    hydride); agostic B-H/C-H (H bonded to another atom) may use b_tm=0.
    """
    sym = atom_el[lig_idx]
    if sym in TM_MONATOMIC_COV_LIGANDS:
        return True
    if sym != "H":
        return False
    for i, j, _ei, _ej in edges:
        if lig_idx not in (i, j):
            continue
        other = j if i == lig_idx else i
        if not base.is_TM(atom_el[other]):
            return False
    return True


def _non_tm_neighbors_in_edges(lig_idx, atom_el, edges):
    out = []
    for i, j, _ei, _ej in edges:
        if lig_idx not in (i, j):
            continue
        other = j if i == lig_idx else i
        if not base.is_TM(atom_el[other]):
            out.append(other)
    return out


def _filter_tm_nm_keys_agostic_shadow_ligands(tm_nm_keys, atom_el, edges):
    """
    Drop M–X from ILP tm_nm_keys when bridging H has X–H and both M–H and M–X exist
    in connectivity (e.g. TM–H–B with B on the coordination list). Keeps M–H only
    so X is not forced to satisfy separate M–X b_tm / lp dative rules; H stays
    bridging/agostic (not isolated hydride).
    """
    key_set = set(tm_nm_keys)
    drop: set[tuple[int, int]] = set()
    for tm, lig in tm_nm_keys:
        if atom_el[lig] != "H":
            continue
        if _force_monatomic_ml_single_cov(lig, atom_el, edges):
            continue
        for x in _non_tm_neighbors_in_edges(lig, atom_el, edges):
            if x == lig or (tm, x) not in key_set:
                continue
            if atom_el[x] not in _SIGMA_AGOSTIC_H_PARENTS:
                continue
            drop.add((tm, x))
    if not drop:
        return tm_nm_keys
    return [(tm, lig) for tm, lig in tm_nm_keys if (tm, lig) not in drop]


def _filter_tm_nm_keys_oxo_bridge_shadow_center(tm_nm_keys, atom_el, edges):
    """
    General O–X–O rule:

    If a metal M coordinates to both O atoms in an O–X–O motif (O and X are directly
    bonded in the ligand skeleton), then drop the geometric M···X contact from ILP
    tm_nm_keys. This prevents the central atom X (often a borderline cutoff contact)
    from becoming a separate M–L variable that overconstrains charge/octet bookkeeping.

    Trigger condition (per metal M and center X):
      - X is a non-TM, non-O atom in tm_nm_keys (i.e., there is a geometric M–X edge)
      - X has at least two O neighbours in *edges* (non-TM neighbours)
      - at least two of those O atoms are also in tm_nm_keys for the same metal M
    """
    key_set = set(tm_nm_keys)
    drop: set[tuple[int, int]] = set()

    lig_comp = _nonmetal_ligand_components([(i, atom_el[i], 0.0, 0.0, 0.0) for i in atom_el], edges)

    for tm, x in tm_nm_keys:
        sx = atom_el.get(x)
        if sx is None or base.is_TM(sx) or sx == "O":
            continue
        o_nbrs = [
            nb
            for nb in _non_tm_neighbors_in_edges(x, atom_el, edges)
            if atom_el.get(nb) == "O"
        ]
        if len(o_nbrs) < 2:
            continue
        # Require O–X–O to be on the same ligand fragment (non-TM component).
        lid_x = lig_comp.get(x)
        if lid_x is None:
            continue
        o_coord = [o for o in o_nbrs if (tm, o) in key_set and lig_comp.get(o) == lid_x]
        if len(o_coord) >= 2:
            drop.add((tm, x))

    if not drop:
        return tm_nm_keys
    return [(tm, lig) for tm, lig in tm_nm_keys if (tm, lig) not in drop]


def _remove_metal_nonmetal_edges(raw_edges):
    """For ring detection only: remove edges where exactly one endpoint is TM."""
    keep = []
    for i, j, ei, ej in raw_edges:
        if base.is_TM(ei) ^ base.is_TM(ej):
            continue
        keep.append((i, j, ei, ej))
    return keep


def _is_tm_nm_edge(ei, ej):
    """True iff exactly one endpoint is a transition metal."""
    return base.is_TM(ei) ^ base.is_TM(ej)


def _tm_nm_orient(i, j, ei, ej):
    """
    For a TM–nonmetal edge, return (tm_idx, lig_idx).
    Precondition: _is_tm_nm_edge(ei, ej).
    """
    if base.is_TM(ei) and not base.is_TM(ej):
        return i, j
    return j, i


def _atoms_bonded_to_tm(edge_list):
    """All atom indices with at least one TM neighbour in edge_list."""
    out = set()
    for i, j, ei, ej in edge_list:
        if base.is_TM(ei):
            out.add(j)
        if base.is_TM(ej):
            out.add(i)
    return out


def _tm_oxidation_ilp_vars(prob, atoms, *, use_discrete_states=True):
    """
    Per-TM oxidation-state variable constrained to TM_COMMON_OXIDATION_STATES[sym].
    Returns {tm_index: int constant or LpVariable}.
    """
    out = {}
    common = base.TM_COMMON_OXIDATION_STATES
    for i, el, *_ in atoms:
        if not base.is_TM(el):
            continue
        if not use_discrete_states:
            raise RuntimeError(
                "TM oxidation states must be discrete and restricted to "
                "TM_COMMON_OXIDATION_STATES in this build."
            )

        allowed = common.get(el)
        if not allowed:
            raise ValueError(
                f"TM element {el!r} not found in TM_COMMON_OXIDATION_STATES; "
                f"cannot constrain oxidation state."
            )
        states = sorted({int(s) for s in allowed})
        if len(states) == 1:
            out[i] = states[0]
            continue
        picks = []
        for s in states:
            y = pulp.LpVariable(f"ox_{i}_{s}", cat="Binary")
            picks.append((y, s))
        prob += pulp.lpSum(y for y, _ in picks) == 1
        ox = pulp.LpVariable(f"ox_{i}", cat="Integer")
        prob += ox == pulp.lpSum(s * y for y, s in picks)
        out[i] = ox
    return out


def _build_adj_from_edges(edges):
    adj = defaultdict(set)
    for i, j, *_ in edges:
        adj[i].add(j)
        adj[j].add(i)
    return adj


def _find_simple_rings(adj, atom_symbol, max_size=12):
    """Simple cycle finder (set-unique), adapted from Lewis-engine style."""
    nodes = sorted(adj.keys())
    rings = []
    seen = set()
    for start in nodes:
        if atom_symbol[start] == "H":
            continue
        stack = [(start, [start], {start})]
        while stack:
            node, path, visited = stack.pop()
            for nb in adj[node]:
                if atom_symbol[nb] == "H":
                    continue
                if nb == start and len(path) >= 4:
                    key = frozenset(path)
                    if key not in seen:
                        seen.add(key)
                        rings.append(sorted(path))
                elif nb not in visited and len(path) < max_size:
                    stack.append((nb, path + [nb], visited | {nb}))
    return [frozenset(r) for r in rings]


def _fused_groups(rings):
    """Group rings that share an edge (>=2 common atoms)."""
    rings = list(rings)
    if not rings:
        return []
    parent = list(range(len(rings)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i in range(len(rings)):
        for j in range(i + 1, len(rings)):
            if len(rings[i] & rings[j]) >= 2:
                union(i, j)

    g = defaultdict(list)
    for i in range(len(rings)):
        g[find(i)].append(i)
    groups = []
    for idxs in g.values():
        s = set()
        for k in idxs:
            s |= set(rings[k])
        groups.append(tuple(sorted(s)))
    return groups


def _minimal_rings(rings):
    """
    Keep only smallest rings in an inclusion sense:
    drop any ring that strictly contains another ring.
    """
    rings = [frozenset(r) for r in rings]
    out = []
    for r in rings:
        has_strict_subset = any((s < r) for s in rings if s is not r)
        if not has_strict_subset:
            out.append(r)
    # de-dup preserve stable order by (size, atoms)
    uniq = {}
    for r in sorted(out, key=lambda x: (len(x), tuple(sorted(x)))):
        uniq[r] = None
    return list(uniq.keys())


def _plane_rmsd(points):
    """
    Heuristic planarity check without external deps:
    try all point triplets as candidate planes and take minimal RMS distance.
    """
    if len(points) <= 3:
        return 0.0
    best = float("inf")
    for a, b, c in itertools.combinations(points, 3):
        ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
        vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
        nx, ny, nz = (
            uy * vz - uz * vy,
            uz * vx - ux * vz,
            ux * vy - uy * vx,
        )
        norm = math.sqrt(nx * nx + ny * ny + nz * nz)
        if norm < 1e-8:
            continue
        d0 = -(nx * a[0] + ny * a[1] + nz * a[2])
        sq = 0.0
        for p in points:
            dist = (nx * p[0] + ny * p[1] + nz * p[2] + d0) / norm
            sq += dist * dist
        rms = math.sqrt(sq / len(points))
        if rms < best:
            best = rms
    return best if best < float("inf") else 999.0


def aromatic_candidate_systems(atoms_packed, raw_edges, max_ring_size=12, planarity_rmsd=0.12):
    atom_symbol = {idx: el for idx, el, *_ in atoms_packed}
    atom_xyz = {idx: (x, y, z) for idx, _el, x, y, z in atoms_packed}
    cut_edges = _remove_metal_nonmetal_edges(raw_edges)
    adj = _build_adj_from_edges(cut_edges)
    rings_all = _find_simple_rings(adj, atom_symbol, max_size=max_ring_size)
    rings_min = _minimal_rings(rings_all)

    # Step 1: smallest planar rings
    planar_min_rings = []
    for r in rings_min:
        syms = [atom_symbol[i] for i in r]
        if any(base.is_TM(s) or s == "H" for s in syms):
            continue
        if sum(1 for s in syms if s in ("C", "N", "O", "S", "P")) < 4:
            continue
        pts = [atom_xyz[i] for i in sorted(r)]
        if _plane_rmsd(pts) <= planarity_rmsd:
            planar_min_rings.append(frozenset(r))

    # Step 2: assemble larger fused systems from planar minimal rings
    groups = _fused_groups(planar_min_rings)
    out = []
    for g in groups:
        out.append(tuple(sorted(g)))
    return out


def _nonmetal_ligand_components(atoms, edges):
    """
    Map each non-TM atom id → ligand component id (connected via non–M–L edges only).
    Intended for single-TM complexes; multiple TM still get separate ligand graphs.
    """
    non_tm = {a[0] for a in atoms if not base.is_TM(a[1])}
    adj = defaultdict(set)
    for i, j, ei, ej in edges:
        if _is_tm_nm_edge(ei, ej):
            continue
        if i in non_tm and j in non_tm:
            adj[i].add(j)
            adj[j].add(i)
    lig_id = {}
    next_id = 0
    for start in non_tm:
        if start in lig_id:
            continue
        stack = [start]
        lig_id[start] = next_id
        while stack:
            cur = stack.pop()
            for nb in adj[cur]:
                if nb not in lig_id:
                    lig_id[nb] = next_id
                    stack.append(nb)
        next_id += 1
    return lig_id


def ligand_ml_contacts_by_component(atoms, edges):
    """
    After connectivity: ligand_id → list of M–L contacts on that ligand fragment.

    Each contact: (tm_idx, lig_idx, distance_Å, metal_sym, lig_sym).
    """
    atom_el = {a[0]: a[1] for a in atoms}
    xyz = {a[0]: (a[2], a[3], a[4]) for a in atoms}
    lig_comp = _nonmetal_ligand_components(atoms, edges)
    out = defaultdict(list)
    for i, j, ei, ej in edges:
        if not _is_tm_nm_edge(ei, ej):
            continue
        tm, lig = _tm_nm_orient(i, j, ei, ej)
        if lig not in lig_comp:
            continue
        lid = lig_comp[lig]
        d = dist(xyz[tm], xyz[lig])
        out[lid].append((tm, lig, d, atom_el[tm], atom_el[lig]))
    return dict(out)


def _terminal_co_tm_c_o_triples(atoms, edges, tm_nm_keys, atom_el):
    """
    One (tm, c, o) per TM–C contact that is a terminal carbonyl: nonmetal component is
    exactly {C,O}, coordination through C, and C–O exists in *edges*.
    """
    lig_comp = _nonmetal_ligand_components(atoms, edges)
    out = []
    for tm, lig in tm_nm_keys:
        if atom_el.get(lig) != "C":
            continue
        lid = lig_comp.get(lig)
        if lid is None:
            continue
        comp = {aid for aid, l in lig_comp.items() if l == lid}
        if len(comp) != 2:
            continue
        syms = {atom_el.get(aid) for aid in comp}
        if syms != {"C", "O"}:
            continue
        c_idx = lig
        o_idx = next(aid for aid in comp if atom_el.get(aid) == "O")
        has_co = False
        for i, j, ei, ej in edges:
            if {i, j} == {c_idx, o_idx} and "C" in (ei, ej) and "O" in (ei, ej):
                has_co = True
                break
        if not has_co:
            continue
        out.append((tm, c_idx, o_idx))
    return out


def _ligand_ml_zcov_pair_weights(
    atoms,
    edges,
    *,
    epsilon=ML_DISTANCE_CLASS_EPSILON,
    same_element_pair_only=True,
):
    """
    Pairwise weights for soft z_cov alignment: similar contact length → same σ vs dative.

    Returns [(weight, (tm, lig_a), (tm, lig_b)), ...].
    """
    by_ligand = ligand_ml_contacts_by_component(atoms, edges)
    pairs = []
    for contacts in by_ligand.values():
        buckets = defaultdict(list)
        for tm, lig, d, sym_m, sym_l in contacts:
            key = (sym_m, sym_l) if same_element_pair_only else 0
            buckets[key].append(((tm, lig), d))
        for items in buckets.values():
            n = len(items)
            for i in range(n):
                (key_a, d_a) = items[i]
                for j in range(i + 1, n):
                    (key_b, d_b) = items[j]
                    w = max(0.0, float(epsilon) - abs(d_a - d_b))
                    if w > 0.0:
                        pairs.append((w, key_a, key_b))
    return pairs


def solve_bond_orders(
    atoms,
    edges,
    aromatic_systems,
    mol_charge=0,
    *,
    metal_adjacency_edges=None,
    tm_ox_minimize_weight=None,
    aromatic_huckel_n=AROMATIC_HUCKEL_N,
    ml_distance_class_weight=None,
    c_lp_tm_neighbor_hard=None,
    retry_relax_c_lp_tm=True,
    __from_c_lp_relaxed_retry=False,
    solve_time_limit_sec=SOLVE_BOND_ORDERS_CBC_TIME_LIMIT_SEC,
):
    apply_c_lp_tm = (
        ILP_HARD_C_LP_ONLY_TM_NEIGHBORS
        if c_lp_tm_neighbor_hard is None
        else bool(c_lp_tm_neighbor_hard)
    )

    prob = pulp.LpProblem("BondOrderAssignmentAromatic", pulp.LpMinimize)

    atom_el = {i: el for i, el, *_ in atoms}

    ilp_edges = [t for t in edges if not _is_tm_nm_edge(t[2], t[3])]
    tm_nm_keys = []
    seen_tm_nm = set()
    for i, j, ei, ej in edges:
        if not _is_tm_nm_edge(ei, ej):
            continue
        tm, lig = _tm_nm_orient(i, j, ei, ej)
        key = (tm, lig)
        if key not in seen_tm_nm:
            seen_tm_nm.add(key)
            tm_nm_keys.append(key)
    tm_nm_keys = _filter_tm_nm_keys_agostic_shadow_ligands(tm_nm_keys, atom_el, edges)
    tm_nm_keys = _filter_tm_nm_keys_oxo_bridge_shadow_center(tm_nm_keys, atom_el, edges)

    u2 = {}
    u3 = {}
    for i, j, _ei, _ej in ilp_edges:
        key = (i, j) if i < j else (j, i)
        u2[key] = pulp.LpVariable(f"u2_{key[0]}_{key[1]}", cat="Binary")
        u3[key] = pulp.LpVariable(f"u3_{key[0]}_{key[1]}", cat="Binary")
        prob += u3[key] <= u2[key]

    # TM–nonmetal: covalent order b_tm in {0,1,2,3}; b_tm=0 means no Lewis M–L bond.
    u2_tm = {}
    u3_tm = {}
    b_tm = {}
    z_cov_tm = {}
    for tm, lig in tm_nm_keys:
        z = pulp.LpVariable(f"zcov_{tm}_{lig}", cat="Binary")
        u2t = pulp.LpVariable(f"u2tm_{tm}_{lig}", cat="Binary")
        u3t = pulp.LpVariable(f"u3tm_{tm}_{lig}", cat="Binary")
        prob += u3t <= u2t
        prob += u2t <= z
        prob += u3t <= z
        z_cov_tm[(tm, lig)] = z
        u2_tm[(tm, lig)] = u2t
        u3_tm[(tm, lig)] = u3t
        s = 1 + u2t + u3t
        b = pulp.LpVariable(f"btm_{tm}_{lig}", lowBound=0, upBound=3, cat="Integer")
        prob += b <= s
        prob += b <= 3 * z
        prob += b >= s - 3 * (1 - z)
        prob += b >= 0
        b_tm[(tm, lig)] = b

    inc_edges = defaultdict(list)
    for i, j, *_ in edges:
        inc_edges[i].append((i, j))
        inc_edges[j].append((i, j))

    def _incident_bond_order_sum(idx):
        terms = []
        for a, b in inc_edges.get(idx, []):
            ea, eb = atom_el[a], atom_el[b]
            if _is_tm_nm_edge(ea, eb):
                tm, lig = _tm_nm_orient(a, b, ea, eb)
                if lig != idx:
                    continue
                key = (tm, lig)
                if key not in b_tm:
                    continue  # agostic-shadowed M–X (ILP keeps M–H only)
                terms.append(b_tm[key])
            else:
                kk = (min(a, b), max(a, b))
                terms.append(1 + u2[kk] + u3[kk])
        return pulp.lpSum(terms) if terms else 0

    # Atoms bonded to a TM in the full structure (e.g. Cl–Au), for aromatic C rules
    # and for excluding M-bound ligand atoms from mol_charge summation.
    tm_neighbor_atoms = _atoms_bonded_to_tm(
        metal_adjacency_edges if metal_adjacency_edges is not None else edges
    )

    lp, oct_plus, oct_minus, q, abs_q, q_neg = {}, {}, {}, {}, {}, {}
    b_oct8_choice = {}
    c_oct8_choice = {}
    si_oct8_choice = {}
    s_oct_choice = {}
    p_oct_choice = {}
    for i, el, *_ in atoms:
        if base.is_tm(el) or el not in base.VALENCE_ELECTRONS:
            continue
        q[i] = pulp.LpVariable(f"q_{i}", lowBound=-4, upBound=4, cat="Integer")
        abs_q[i] = pulp.LpVariable(f"absq_{i}", lowBound=0, cat="Integer")
        q_neg[i] = pulp.LpVariable(f"qneg_{i}", lowBound=0, cat="Integer")
        if base.is_ionlike_s_block_metal(el):
            continue
        lp[i] = pulp.LpVariable(f"lp_{i}", lowBound=0, cat="Integer")
        oct_plus[i] = pulp.LpVariable(f"octp_{i}", lowBound=0, cat="Integer")
        oct_minus[i] = pulp.LpVariable(f"octm_{i}", lowBound=0, cat="Integer")
        if el == "B":
            # For boron: choose 6e (y=0) vs 8e (y=1) local-electron target.
            b_oct8_choice[i] = pulp.LpVariable(f"b8_{i}", cat="Binary")
        if el == "C" and i in tm_neighbor_atoms:
            # For TM-bound carbon (carbene / related): allow 6e vs 8e local-electron targets.
            c_oct8_choice[i] = pulp.LpVariable(f"c8_{i}", cat="Binary")
        if (
            el == "Si"
            and i in tm_neighbor_atoms
        ):
            # For TM-bound silicon: choose 6e (y=0) vs 8e (y=1) local-electron target.
            # (Models silylene-like :SiR2 donation while preserving hard-octet.)
            si_oct8_choice[i] = pulp.LpVariable(f"si8_{i}", cat="Binary")
        if el == "S":
            # For sulfur: allow 8e/10e/12e expanded-octet targets (equal preference).
            y10 = pulp.LpVariable(f"s10_{i}", cat="Binary")
            y12 = pulp.LpVariable(f"s12_{i}", cat="Binary")
            prob += y10 + y12 <= 1
            s_oct_choice[i] = (y10, y12)
        if el == "P":
            # For phosphorus: allow 8e/10e/12e expanded-octet targets (equal preference).
            y10 = pulp.LpVariable(f"p10_{i}", cat="Binary")
            y12 = pulp.LpVariable(f"p12_{i}", cat="Binary")
            prob += y10 + y12 <= 1
            p_oct_choice[i] = (y10, y12)

    for i, el, *_ in atoms:
        if i not in q:
            continue
        bond_sum = _incident_bond_order_sum(i)
        ve = base.VALENCE_ELECTRONS[el]
        if base.is_ionlike_s_block_metal(el):
            prob += q[i] == ve - bond_sum
            prob += abs_q[i] >= q[i]
            prob += abs_q[i] >= -q[i]
            prob += q_neg[i] >= -q[i]
            continue
        if el == "B" and i in b_oct8_choice:
            oct_target = 6 + 2 * b_oct8_choice[i]
        elif el == "C" and i in c_oct8_choice:
            oct_target = 6 + 2 * c_oct8_choice[i]
        elif el == "Si" and i in si_oct8_choice:
            oct_target = 6 + 2 * si_oct8_choice[i]
        elif el == "S" and i in s_oct_choice:
            y10, y12 = s_oct_choice[i]
            oct_target = 8 + 2 * y10 + 4 * y12
        elif el == "P" and i in p_oct_choice:
            y10, y12 = p_oct_choice[i]
            oct_target = 8 + 2 * y10 + 4 * y12
        else:
            oct_target = 2 if el in ("H", "Li") else 8
        local_e = 2 * lp[i] + 2 * bond_sum
        assigned_e = 2 * lp[i] + bond_sum
        prob += local_e - oct_target == oct_plus[i] - oct_minus[i]
        if ILP_HARD_OCTET:
            prob += oct_plus[i] == 0
            prob += oct_minus[i] == 0
        prob += q[i] == ve - assigned_e
        prob += abs_q[i] >= q[i]
        prob += abs_q[i] >= -q[i]
        prob += q_neg[i] >= -q[i]

    for i, j, ei, ej in ilp_edges:
        if "H" in (ei, ej) or "Cl" in (ei, ej):
            key = (min(i, j), max(i, j))
            prob += u2[key] == 0
            prob += u3[key] == 0

    heavy_deg = defaultdict(int)
    for i, j, ei, ej in edges:
        if ei != "H":
            heavy_deg[i] += 1
        if ej != "H":
            heavy_deg[j] += 1

    for i, j, ei, ej in ilp_edges:
        key = (min(i, j), max(i, j))
        if (ei == "H" or ej == "H"):
            prob += u3[key] == 0
            continue
        if (base.is_TM(ei) or base.is_TM(ej)):
            prob += u3[key] == 0
            continue
        if (heavy_deg[i] > 2 or heavy_deg[j] > 2):
            prob += u3[key] == 0

    # Monatomic ML single covalent rule:
    # - F/Cl/Br/I at TM are forced to a covalent single bond (b_tm=1).
    # - H is forced only when isolated (no non-TM neighbor).
    for tm, lig in tm_nm_keys:
        if not _force_monatomic_ml_single_cov(lig, atom_el, edges):
            continue
        prob += z_cov_tm[(tm, lig)] == 1
        prob += u2_tm[(tm, lig)] == 0
        prob += u3_tm[(tm, lig)] == 0

    # Terminal CO (2-atom C+O fragment, M–L through C): C≡O triple; M–C dative (b_tm=0).
    terminal_cos = _terminal_co_tm_c_o_triples(atoms, edges, tm_nm_keys, atom_el)
    co_triple_done = set()
    for tm, c_idx, o_idx in terminal_cos:
        kk = (min(c_idx, o_idx), max(c_idx, o_idx))
        if kk not in u2:
            raise RuntimeError(
                f"Terminal CO: C–O edge {kk} not in ILP edge set (check connectivity)"
            )
        if kk not in co_triple_done:
            prob += u2[kk] == 1
            prob += u3[kk] == 1
            co_triple_done.add(kk)
        prob += b_tm[(tm, c_idx)] == 0

    # --- Aromatic 4n+2 penalty ---
    aromatic_dev_terms = []
    if ILP_WEIGHT_AROMATIC_DEVIATION > 0:
        for sys_idx, sys_atoms in enumerate(aromatic_systems):
            # pi_e = 2 * (#internal double-bond EDGES) + 2 per ring atom with LP not on a ring double bond.
            # Count double bonds by EDGE (once per bond), not by per-atom incidence.
            sys_atom_set = set(sys_atoms)
            sys_edges = sorted(
                {
                    (min(a, b), max(a, b))
                    for i in sys_atoms
                    for a, b in inc_edges.get(i, [])
                    if (a in sys_atom_set and b in sys_atom_set)
                }
            )
            internal_double_edges = pulp.lpSum(
                u2[(a, b)] for a, b in sys_edges
            ) if sys_edges else 0
            internal_triple_edges = pulp.lpSum(
                u3[(a, b)] for a, b in sys_edges
            ) if sys_edges else 0

            double_present = {}
            lp_pair = {}
            for i in sys_atoms:
                dsum = pulp.lpSum(
                    u2[(min(a, b), max(a, b))] for a, b in inc_edges.get(i, [])
                    if (a in sys_atom_set and b in sys_atom_set)
                ) if inc_edges.get(i) else 0
                dp = pulp.LpVariable(f"sys{sys_idx}_dblp_{i}", cat="Binary")
                double_present[i] = dp
                deg = max(
                    1,
                    len([
                        1
                        for a, b in inc_edges.get(i, [])
                        if a in sys_atom_set and b in sys_atom_set
                    ]),
                )
                prob += dsum >= dp
                prob += dsum <= deg * dp

                if i in lp:
                    lpb = pulp.LpVariable(f"sys{sys_idx}_lppair_{i}", cat="Binary")
                    lp_pair[i] = lpb
                    prob += lp[i] >= lpb
                    prob += lpb <= 1 - dp
                else:
                    lp_pair[i] = 0

            pi_e = (
                2 * (internal_double_edges - internal_triple_edges)
                + 4 * internal_triple_edges
                + pulp.lpSum(2 * lp_pair[i] for i in sys_atoms)
            )
            dev_p = pulp.LpVariable(f"sys{sys_idx}_devp", lowBound=0, cat="Integer")
            dev_m = pulp.LpVariable(f"sys{sys_idx}_devm", lowBound=0, cat="Integer")
            if aromatic_huckel_n is not None:
                pi_target = 4 * int(aromatic_huckel_n) + 2
                prob += pi_e - pi_target == dev_p - dev_m
            else:
                max_pi = max(2, 4 * len(sys_atoms))
                kmax = max(0, (max_pi - 2) // 4)
                k = pulp.LpVariable(
                    f"sys{sys_idx}_k", lowBound=0, upBound=kmax, cat="Integer"
                )
                prob += pi_e - (4 * k + 2) == dev_p - dev_m
            aromatic_dev_terms.append(dev_p + dev_m)

    tm_ox_vars = _tm_oxidation_ilp_vars(prob, atoms, use_discrete_states=True)
    sigma_cov_ml_all = (
        pulp.lpSum(b_tm[k] for k in tm_nm_keys) if tm_nm_keys else 0
    )
    eta_lig_by_metal = _eta_ligand_atoms_by_metal(atoms, edges)
    sigma_cov_ml_non_eta_terms = [
        b_tm[(tm, lig)]
        for tm, lig in tm_nm_keys
        if lig not in eta_lig_by_metal.get(tm, ())
    ]
    sigma_cov_ml_non_eta = (
        pulp.lpSum(sigma_cov_ml_non_eta_terms) if sigma_cov_ml_non_eta_terms else 0
    )

    # Hard: dative-only M–L (b_tm == 0) requires lp >= 1 on the ligand atom, except:
    #   • η/haptic atoms (lp == 0 allowed);
    #   • bridging/agostic H (non-TM neighbour in edges).
    # M–X dropped from tm_nm_keys when X is agostic parent of bridging H (M–H–X).
    for tm, lig in tm_nm_keys:
        if lig in eta_lig_by_metal.get(tm, ()):
            continue
        if lig not in lp:
            continue
        if atom_el[lig] == "H" and not _force_monatomic_ml_single_cov(lig, atom_el, edges):
            continue
        # b_tm ∈ {0,1,2,3}. Enforce lp>=1 when b_tm==0, else lp>=0.
        prob += lp[lig] >= 1 - b_tm[(tm, lig)]

    ligand_q_sum = pulp.lpSum(q[i] for i in q) if q else 0
    tm_ox_sum = pulp.lpSum(
        v if not isinstance(v, int) else v for v in tm_ox_vars.values()
    )
    if ILP_HARD_MOL_CHARGE_BALANCE:
        prob += ligand_q_sum + tm_ox_sum - sigma_cov_ml_all == mol_charge

    if ILP_HARD_OX_GE_SIGMA and tm_ox_vars:
        prob += tm_ox_sum >= sigma_cov_ml_non_eta

    octet_penalty = pulp.lpSum(oct_plus.values()) + pulp.lpSum(oct_minus.values())
    formal_charge_penalty = pulp.lpSum(abs_q.values()) if abs_q else 0
    tm_nm_high_order = (
        pulp.lpSum(u2_tm[k] + u3_tm[k] for k in tm_nm_keys) if tm_nm_keys else 0
    )
    eneg_neg_charge_penalty = pulp.lpSum(
        max(0.0, 4.0 - base.ENEG.get(el, 2.0)) * q_neg[i]
        for i, el, *_ in atoms if i in q_neg
    )
    aromatic_penalty = pulp.lpSum(aromatic_dev_terms) if aromatic_dev_terms else 0

    eta_internal_keys = _eta_group_internal_lp_edge_keys(atoms, edges)
    eta_group_u2_sum = (
        pulp.lpSum(u2[k] for k in eta_internal_keys if k in u2)
        if eta_internal_keys
        else 0
    )

    ml_zcov_terms = []
    for pidx, (w_ij, key_a, key_b) in enumerate(
        _ligand_ml_zcov_pair_weights(atoms, edges)
    ):
        if key_a not in z_cov_tm or key_b not in z_cov_tm:
            continue
        za, zb = z_cov_tm[key_a], z_cov_tm[key_b]
        diff = pulp.LpVariable(f"mlzc_{pidx}", cat="Binary")
        prob += diff >= za - zb
        prob += diff >= zb - za
        ml_zcov_terms.append(w_ij * diff)
    ml_zcov_penalty = pulp.lpSum(ml_zcov_terms) if ml_zcov_terms else 0

    ml_zcov_w = (
        ILP_WEIGHT_ML_DISTANCE_CLASS
        if ml_distance_class_weight is None
        else ml_distance_class_weight
    )
    # Note: we intentionally do NOT add a "minimize TM oxidation state" soft objective
    # (tm_ox_penalty). TM oxidation variables remain present for hard constraints.

    remote_c_lp_violation = {}
    if apply_c_lp_tm:
        for i, el, *_ in atoms:
            if el == "C" and i in lp and i not in tm_neighbor_atoms:
                # Elastic remote-C lp=0: allow lp>0 but penalize it heavily.
                v = pulp.LpVariable(f"rcv_{i}", lowBound=0, cat="Integer")
                remote_c_lp_violation[i] = v
                prob += lp[i] <= v

    if ILP_HARD_ETA_CARBON_LP_ZERO:
        for i in _eta_carbon_atom_ids(atoms, edges):
            if i in lp:
                prob += lp[i] == 0

    objective_terms = []
    if ILP_WEIGHT_FORMAL_CHARGE > 0:
        objective_terms.append(ILP_WEIGHT_FORMAL_CHARGE * formal_charge_penalty)
    if ILP_WEIGHT_ENEG_NEGATIVE_FC > 0:
        objective_terms.append(ILP_WEIGHT_ENEG_NEGATIVE_FC * eneg_neg_charge_penalty)
    if ILP_WEIGHT_AROMATIC_DEVIATION > 0:
        objective_terms.append(ILP_WEIGHT_AROMATIC_DEVIATION * aromatic_penalty)
    if ILP_WEIGHT_REMOTE_C_LP_VIOLATION > 0 and remote_c_lp_violation:
        objective_terms.append(
            ILP_WEIGHT_REMOTE_C_LP_VIOLATION
            * pulp.lpSum(remote_c_lp_violation.values())
        )
    if ml_zcov_w > 0:
        objective_terms.append(ml_zcov_w * ml_zcov_penalty)
    if ILP_WEIGHT_ETA_GROUP_MAX_DOUBLE_BONDS > 0 and eta_internal_keys:
        objective_terms.append(
            -ILP_WEIGHT_ETA_GROUP_MAX_DOUBLE_BONDS * eta_group_u2_sum
        )

    if not objective_terms:
        raise RuntimeError("ILP objective is empty: enable at least one ILP_WEIGHT_* term")
    prob += pulp.lpSum(objective_terms)

    cbc_opts: dict = {"msg": False}
    if solve_time_limit_sec is not None and float(solve_time_limit_sec) > 0:
        cbc_opts["timeLimit"] = float(solve_time_limit_sec)
    status = prob.solve(pulp.PULP_CBC_CMD(**cbc_opts))
    st = pulp.LpStatus[status]
    if st not in ("Optimal", "Integer Feasible"):
        # Infeasible / Unbounded are model issues, not wall-clock timeout.
        _timeout_like = {"Not Solved", "Undefined"}
        if (
            solve_time_limit_sec
            and float(solve_time_limit_sec) > 0
            and st in _timeout_like
        ):
            raise RuntimeError(
                f"ILP timeout: solve_bond_orders hit CBC timeLimit "
                f"{float(solve_time_limit_sec):g}s (status {st})"
            )
        raise RuntimeError(
            f"ILP failed: {st} "
            "(check mol_charge, η-carbon lp=0 vs ring anion/aromatic 6π, connectivity, "
            "remote-C lp preference, and LP/octet balance on donors)"
        )

    bonds = []
    for i, j, _ei, _ej in ilp_edges:
        key = (min(i, j), max(i, j))
        order = 1 + int(round(pulp.value(u2[key]))) + int(round(pulp.value(u3[key])))
        bonds.append((key[0], key[1], order))
    for tm, lig in tm_nm_keys:
        b_val = int(round(pulp.value(b_tm[(tm, lig)])))
        if b_val <= 0:
            continue
        a, b = (tm, lig) if tm < lig else (lig, tm)
        bonds.append((a, b, b_val))
    bonds.sort()
    lp_out = {i: int(round(pulp.value(v))) for i, v in lp.items()}
    fc_out = {i: int(round(pulp.value(v))) for i, v in q.items()}

    # Stash remote-C lp violations for printing by the CLI.
    base.LAST_REMOTE_C_LP_VIOLATIONS = []
    if remote_c_lp_violation:
        viol = sorted(
            (
                i,
                atom_el.get(i, "?"),
                lp_out.get(i, 0),
                int(round(pulp.value(remote_c_lp_violation[i]))),
            )
            for i in remote_c_lp_violation
            if lp_out.get(i, 0) > 0
        )
        base.LAST_REMOTE_C_LP_VIOLATIONS = viol
    for tm_i, ox_v in tm_ox_vars.items():
        if isinstance(ox_v, int):
            fc_out[tm_i] = ox_v
        else:
            fc_out[tm_i] = int(round(pulp.value(ox_v)))
    return bonds, lp_out, fc_out


def _array_idx(atom_id: int, n_atoms: int) -> int:
    """Map 1-based XYZ atom numbers to 0-based array indices when needed."""
    if atom_id >= n_atoms and atom_id >= 1:
        return atom_id - 1
    return atom_id


def _ml_bo_zero_from_adj(bo0, metal_adjacency_edges, n_atoms):
    """M–L contacts with ILP Lewis bond order 0; keys are 0-based (tm, lig)."""
    out = set()
    if not metal_adjacency_edges:
        return out
    for i, j, ei, ej in metal_adjacency_edges:
        if not is_TM(ei) ^ is_TM(ej):
            continue
        tm_id, lig_id = (i, j) if is_TM(ei) else (j, i)
        tm_i = _array_idx(tm_id, n_atoms)
        lig_i = _array_idx(lig_id, n_atoms)
        if bo0.get((min(tm_i, lig_i), max(tm_i, lig_i)), 0) == 0:
            out.add((tm_i, lig_i))
    return out


def _expand_eta_dative_ligands(
    metal_idx,
    l_seeds,
    ml_bo_zero,
    atoms,
    adj_lewis,
    *,
    coord_indices=None,
    lig_comp=None,
):
    """
    All M–L (b_tm=0) contacts in the same η coordinating group as any CBC L seed.
    """
    metal_contacts = {lig for (tm, lig) in ml_bo_zero if tm == metal_idx}
    seeds = [a for a in l_seeds if a in metal_contacts]
    if not seeds:
        return set()
    if coord_indices is not None and lig_comp is not None:
        out = set()
        for group in _eta_coordinating_groups(
            metal_idx, coord_indices, lig_comp, adj_lewis
        ):
            if any(s in group for s in seeds):
                out.update(i for i in group if (metal_idx, i) in ml_bo_zero)
        return out
    visited = set(seeds)
    queue = list(seeds)
    while queue:
        cur = queue.pop()
        for nb in adj_lewis[cur]:
            if nb not in metal_contacts:
                continue
            if nb not in visited:
                visited.add(nb)
                queue.append(nb)
    return visited


def infer_dative_ml_pairs_cbc(
    atoms,
    coords,
    bonds,
    lp,
    fc,
    *,
    metal_adjacency_edges=None,
    edges=None,
):
    """
    CBC-stage dative M–L for SMILES: ILP b_tm=0 and CBC type L.

    Single-atom L donors (e.g. N with lp) → one dative arrow each.
    η-type L (including paired atom tuples) → one dative per ring/contact atom
    with b_tm=0 (e.g. Cp η⁵ → five C->M arrows).
    Returns (tm, lig) with 1-based atom indices (ligand donates to TM).
    """
    adjacency = metal_adjacency_edges if metal_adjacency_edges is not None else edges
    if adjacency is None:
        return []

    n_atoms = len(atoms)
    bo0 = {}
    for i, j, o in bonds:
        bo0[(min(i - 1, j - 1), max(i - 1, j - 1))] = int(o)

    adj_lewis = defaultdict(set)
    for (i, j), order in bo0.items():
        if order > 0:
            adj_lewis[i].add(j)
            adj_lewis[j].add(i)

    lp0 = {k - 1: v for k, v in lp.items()} if lp else {}
    if isinstance(fc, list):
        fc0 = fc
    else:
        fc0 = [fc.get(k, 0) for k in range(n_atoms)]

    ml_bo_zero = _ml_bo_zero_from_adj(bo0, adjacency, n_atoms)
    lig_comp = _ligand_component_by_array_index(atoms, bo0, adjacency)

    cbc, _neighbor_cbc = classify_cbc_ligands(
        atoms,
        coords,
        bo0,
        lp0,
        fc0,
        metal_adjacency_edges=adjacency,
    )

    def bo_ij(m, lig):
        return bo0.get((min(m, lig), max(m, lig)), 0)

    coord_by_metal = defaultdict(set)
    for tm, lig in ml_bo_zero:
        coord_by_metal[tm].add(lig)
    for tm, lig, ei, ej in adjacency:
        if not is_TM(ei) ^ is_TM(ej):
            continue
        tmi = _array_idx(tm, n_atoms)
        ligi = _array_idx(lig, n_atoms)
        coord_by_metal[tmi].add(ligi)

    out = []
    seen = set()
    for metal_idx, records in cbc.items():
        l_seeds = set()
        for atom_tuple, cbc_char in records:
            if cbc_char != "L":
                continue
            for a in atom_tuple:
                l_seeds.add(a)

        dative_ligs = set(l_seeds)
        if l_seeds:
            dative_ligs |= _expand_eta_dative_ligands(
                metal_idx,
                l_seeds,
                ml_bo_zero,
                atoms,
                adj_lewis,
                coord_indices=coord_by_metal.get(metal_idx, ()),
                lig_comp=lig_comp,
            )

        for lig in sorted(dative_ligs):
            if (metal_idx, lig) not in ml_bo_zero:
                continue
            if bo_ij(metal_idx, lig) != 0:
                continue
            pair = (metal_idx + 1, lig + 1)
            if pair not in seen:
                seen.add(pair)
                out.append(pair)
    return out


def infer_dative_ml_pairs(atoms, edges, bonds, *, coords=None, lp=None, fc=None):
    """
    Dative M–L for SMILES: use infer_dative_ml_pairs_cbc when coords/lp/fc supplied;
    otherwise only geometric candidates with ILP order 0 (no CBC filter).
    """
    if coords is not None and lp is not None and fc is not None:
        return infer_dative_ml_pairs_cbc(
            atoms,
            coords,
            bonds,
            lp,
            fc,
            metal_adjacency_edges=edges,
            edges=edges,
        )
    orders = {}
    for i, j, o in bonds:
        orders[(min(i, j), max(i, j))] = int(o)
    out = []
    for i, j, ei, ej in edges:
        if not _is_tm_nm_edge(ei, ej):
            continue
        tm, lig = _tm_nm_orient(i, j, ei, ej)
        if orders.get((min(tm, lig), max(tm, lig)), 0) == 0:
            out.append((tm, lig))
    return out


def _require_rdkit():
    if Chem is None:
        raise ImportError(
            "RDKit is required for SMILES export. Install with: pip install rdkit"
        )


def _ilp_bond_order_to_rdkit(order: int):
    _require_rdkit()
    if order <= 1:
        return Chem.BondType.SINGLE
    if order == 2:
        return Chem.BondType.DOUBLE
    if order >= 3:
        return Chem.BondType.TRIPLE
    return Chem.BondType.SINGLE


def ilp_to_rdkit_mol(atoms, bonds, fc_out, *, dative_ml_pairs=None, edges=None):
    """
  Build an RDKit Mol from ILP results.

  *atoms*: (idx, symbol, x, y, z) with 1-based *idx*.
  *bonds*: (i, j, order) for covalent edges (including covalent M–L); 1-based.
  *fc_out*: formal charge per idx; TM entries are oxidation states for [M+n] labels.
  *dative_ml_pairs*: (tm, lig) 1-based from ``infer_dative_ml_pairs_cbc``; η rings
  get one ``lig->M`` per contact atom (e.g. five arrows for Cp η⁵).
    """
    _require_rdkit()
    if dative_ml_pairs is None:
        raise ValueError(
            "ilp_to_rdkit_mol: provide dative_ml_pairs from infer_dative_ml_pairs_cbc"
        )

    pt = Chem.GetPeriodicTable()
    idx_to_rd = {}
    rw = Chem.RWMol()
    for idx, sym, *_ in atoms:
        z = pt.GetAtomicNumber(sym)
        if z <= 0:
            raise ValueError(f"Unknown element {sym!r} for atom {idx}")
        atom = Chem.Atom(z)
        if is_tm(sym):
            if idx in fc_out:
                atom.SetFormalCharge(int(fc_out[idx]))
        elif idx in fc_out:
            atom.SetFormalCharge(int(fc_out[idx]))
        idx_to_rd[idx] = rw.AddAtom(atom)

    added = set()

    def _add_edge(i, j, btype):
        key = (min(i, j), max(i, j))
        if key in added:
            return
        ri, rj = idx_to_rd[i], idx_to_rd[j]
        rw.AddBond(ri, rj, btype)
        added.add(key)

    for i, j, order in bonds:
        if int(order) <= 0:
            continue
        _add_edge(i, j, _ilp_bond_order_to_rdkit(int(order)))

    for tm, lig in dative_ml_pairs:
        _add_edge(lig, tm, Chem.BondType.DATIVE)

    mol = rw.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Chem.rdchem.MolSanitizeException:
        mol.UpdatePropertyCache(strict=False)
    return mol


def ilp_to_smiles(
    atoms,
    bonds,
    fc_out,
    *,
    dative_ml_pairs=None,
    edges=None,
    canonical=False,
):
    """RDKit SMILES: covalent orders as ``=`` / ``#``; dative M–L as ``->`` / ``<-``."""
    mol = ilp_to_rdkit_mol(
        atoms,
        bonds,
        fc_out,
        dative_ml_pairs=dative_ml_pairs,
        edges=edges,
    )
    return Chem.MolToSmiles(mol, canonical=canonical)


def print_ilp_smiles_report(
    atoms,
    bonds,
    fc_out,
    *,
    dative_ml_pairs=None,
    edges=None,
    canonical=False,
):
    """Print SMILES derived from ILP bond orders and dative M–L contacts."""
    print()
    print("=" * 72)
    print("  SMILES Generation")
    print("=" * 72)
    try:
        smiles = ilp_to_smiles(
            atoms,
            bonds,
            fc_out,
            dative_ml_pairs=dative_ml_pairs,
            edges=edges,
            canonical=canonical,
        )
        n_dative = len(
            dative_ml_pairs
            if dative_ml_pairs is not None
            else (infer_dative_ml_pairs(atoms, edges, bonds) if edges else [])
        )
        print(f"  {smiles}")
    except ImportError as exc:
        print(f"  (skipped: {exc})")
    except Exception as exc:
        print(f"  ⚠ SMILES generation failed: {exc}")
    print()


def tm_oxidation_from_charge_balance(
    atom_syms,
    fc,
    mol_charge: int = 0,
    bo=None,
    *,
    cbc_interaction_records=None,
    cbc_neighbor_types=None,
) -> dict[int, tuple[int, int, int]]:
    """
    For each TM:
      ox(M) = mol_charge − Σfc(ligand atoms) + σ.

    σ is Σ(ILP M–L bond orders) when no CBC filter is given (ILP Lewis only).
    After CBC, pass *cbc_interaction_records* (from ``classify_cbc_ligands``):
    only ligand atoms in a record classified **X** add their ILP bond order;
    **L** and **Z** contribute 0.  Prefer this over *cbc_neighbor_types*.

    Returns {tm_idx: (ox_state, sum_fc_ligands, sum_sigma)}.
    """
    n = len(atom_syms)
    if len(fc) != n:
        raise ValueError(f"fc length {len(fc)} != atom count {n}")
    cbc_x_ligands = None
    if cbc_interaction_records is not None:
        cbc_x_ligands = cbc_x_ligands_from_interaction_records(cbc_interaction_records)
    elif cbc_neighbor_types is not None:
        cbc_x_ligands = {
            m: {lig for lig, t in nbrs.items() if t == "X"}
            for m, nbrs in cbc_neighbor_types.items()
        }
    out: dict[int, tuple[int, int, int]] = {}
    for m, sym in enumerate(atom_syms):
        if not base.is_TM(sym):
            continue
        ligands_sum = sum(
            fc[j] for j in range(n) if j != m and not base.is_TM(atom_syms[j])
        )
        sigma_sum = 0
        if bo is not None:
            x_ligs = cbc_x_ligands.get(m, ()) if cbc_x_ligands is not None else None
            for lig in range(n):
                if lig == m or base.is_TM(atom_syms[lig]):
                    continue
                order = bo.get((min(m, lig), max(m, lig)), 0)
                if order <= 0:
                    continue
                if x_ligs is not None and lig not in x_ligs:
                    continue
                sigma_sum += order
        ox = int(mol_charge) - int(ligands_sum) + int(sigma_sum)
        out[m] = (ox, int(ligands_sum), int(sigma_sum))
    return out


def print_tm_oxidation_sigma_report(
    atom_syms,
    fc,
    mol_charge: int = 0,
    bo=None,
    *,
    cbc_interaction_records=None,
    cbc_neighbor_types=None,
) -> None:
    """Report TM oxidation from charge / σ balance (σ uses CBC X records when provided)."""
    common = base.TM_COMMON_OXIDATION_STATES
    results = tm_oxidation_from_charge_balance(
        atom_syms,
        fc,
        mol_charge,
        bo=bo,
        cbc_interaction_records=cbc_interaction_records,
        cbc_neighbor_types=cbc_neighbor_types,
    )
    if not results:
        return

    print()
    print("=" * 72)
    print("  TM Oxidation State (from formal-charge balance)")
    print("=" * 72)
    if cbc_interaction_records is not None or cbc_neighbor_types is not None:
        print(
            "  Rule (after CBC): mol_charge = Σfc(ligands) + ox(M)"
            " − σ, where σ = Σ(ILP bond order) on ligand atoms in CBC records"
            " classified X (same as the CBC table above)."
        )
    else:
        print(
            "  Rule (ILP): mol_charge = Σfc(ligands) + ox(M)"
            " − Σ(covalent M–L bond orders in Lewis *bo*)."
        )
    print("  Compare ox(M) to TM_COMMON_OXIDATION_STATES.")
    print()

    for m in sorted(results):
        sym = atom_syms[m]
        ox_state, ligands_sum, sigma_sum = results[m]
        label = f"{sym}{m + 1}"

        ref = common.get(sym)
        if ref is None:
            status = "—  (no TM_COMMON_OXIDATION_STATES entry)"
        elif ox_state in ref:
            status = f"√  (common for {sym}: {ref})"
        else:
            status = f"⚠ WARNING: oxidation state {ox_state} not in common oxidation states {ref}"

        print(
            f"  {label}:  oxidation state = {ox_state} = "
            f"[{mol_charge} − ({ligands_sum}) + (σ={sigma_sum})]"
        )
        print(f"           {status}")
        print()


def main():
    if len(sys.argv) not in (2, 3):
        raise SystemExit(
            "Usage: python Lewis-engine-ILP.py <molecule.xyz> [molecular_charge]"
        )
    xyz = sys.argv[1]
    charge = int(sys.argv[2]) if len(sys.argv) == 3 else 0
    atoms = base.read_xyz(xyz)
    raw = base.connectivity(atoms)
    aromatic_systems = aromatic_candidate_systems(atoms, raw)

    bonds, lp_out, fc_out = solve_bond_orders(
        atoms, raw, aromatic_systems, mol_charge=charge, metal_adjacency_edges=raw
    )
    atom_syms = [a[1] for a in atoms]
    print(f"Read {len(atoms)} atoms from {xyz} (charge={charge})")
    print(f"Aromatic candidate systems (assembled from minimal planar rings): {len(aromatic_systems)}")
    for i, s in enumerate(aromatic_systems, start=1):
        print(f"  System {i}: atoms {list(s)}")
    print()
    base.print_summary_and_choose_ilp(atom_syms, bonds, lp_out, charge)

    coords = [[a[2], a[3], a[4]] for a in atoms]
    bo0 = {(i - 1, j - 1): o for i, j, o in bonds}
    lp_full = {i - 1: lp_out.get(i, 0) for i in range(1, len(atom_syms) + 1)}
    fc0 = [fc_out.get(a[0], 0) for a in atoms]
    idx_to_pos = {a[0]: k for k, a in enumerate(atoms)}
    metal_adj_0 = [
        (idx_to_pos[tm], idx_to_pos[lig], ei, ej)
        for tm, lig, ei, ej in raw
        if base.is_TM(ei) ^ base.is_TM(ej)
    ]
    base.print_octet_report(
        atom_syms, bo0, lp_full, fc0,
        metal_adjacency_edges=metal_adj_0, coords=coords,
    )
    viol = getattr(base, "LAST_REMOTE_C_LP_VIOLATIONS", None) or []
    if viol:
        items = ", ".join(f"{sym}{i}(lp={lpv},v={vv})" for i, sym, lpv, vv in viol[:25])
        more = "" if len(viol) <= 25 else f" ... +{len(viol) - 25} more"
        print(f"  ⚠ remote-C lp>0 used (elastic): {items}{more}")

    lp_by_arr = {k: lp_out.get(a[0], 0) for k, a in enumerate(atoms)}
    cbc_bundle = classify_cbc_ligands(
        atom_syms, coords, bo0, lp_by_arr, fc0, charge, metal_adjacency_edges=metal_adj_0
    )
    base.print_cbc_report(
        atom_syms, coords, bo0, lp_by_arr, fc0, charge,
        metal_adjacency_edges=metal_adj_0, cbc_bundle=cbc_bundle,
    )
    cbc_interaction_records, _cbc_neighbor_types = cbc_bundle
    print_tm_oxidation_sigma_report(
        atom_syms, fc0, charge, bo=bo0, cbc_interaction_records=cbc_interaction_records
    )
    dative_ml_pairs = infer_dative_ml_pairs_cbc(
        atom_syms,
        coords,
        bonds,
        lp_out,
        fc0,
        metal_adjacency_edges=metal_adj_0,
    )
    print_ilp_smiles_report(
        atoms, bonds, fc_out, dative_ml_pairs=dative_ml_pairs, edges=raw
    )


if __name__ == "__main__":
    main()

