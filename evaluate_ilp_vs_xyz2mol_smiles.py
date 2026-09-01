#!/usr/bin/env python3
"""All pairwise SMILES comparisons among ILP and xyz2mol_tm methods.

Methods (4 → 6 unordered pairs = 3+2+1):
  ILP, CSD, NBO(DFT), Hückel(DFT)

Metrics:
  1. ligand canonical SMILES after strip_stereo / isomericSmiles=False
     (ILP covalent X rewritten as anions); plus resonance-aware ligand match
     and NHC/CAAC + isocyanide equivalence
  2. metal neighbor element multiset (bond order ignored)
  3. metal oxidation state (formal charge on the TM atom)
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent
DEFAULT_ILP = ROOT / "tmqmg_ilp_benchmark_output" / "per_structure_results.csv"
DEFAULT_REF = ROOT / "tmqmg_smiles.csv"
DEFAULT_OUT = ROOT / "evaluate_ilp_vs_xyz2mol"

FAIL_TOKENS = {
    "",
    "fail",
    "API_smiles_missing",
    "not_in_database",
}

TM_ATOMIC_NUMS = {
    21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
    39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
    57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,
    72, 73, 74, 75, 76, 77, 78, 79, 80,
}

REF_COLUMNS = {
    "CSD": "smiles_CSD_fix",
    "NBO(DFT)": "smiles_NBO_DFT_xyz",
    "Hückel(DFT)": "smiles_huckel_DFT_xyz",
}

METHODS = ("ILP", "CSD", "NBO(DFT)", "Hückel(DFT)")
PAIRS = [
    (METHODS[i], METHODS[j])
    for i in range(len(METHODS))
    for j in range(i + 1, len(METHODS))
]

_METAL_NON = (
    "[#3,#11,#12,#13,#19,#21,#22,#23,#24,#25,#26,#27,#28,#29,#30,"
    "#39,#40,#41,#42,#43,#44,#45,#46,#47,#48,"
    "#57,#58,#59,#60,#61,#62,#63,#64,#65,#66,#67,#68,#69,#70,#71,"
    "#72,#73,#74,#75,#76,#77,#78,#79,#80]~"
    "[#1,B,#6,#14,#15,#33,#51,#16,#34,#52,Cl,Br,I,#85]"
)


def _disconnector():
    params = rdMolStandardize.MetalDisconnectorOptions()
    params.splitAromaticC = True
    params.splitGrignards = True
    params.adjustCharges = False
    mdis = rdMolStandardize.MetalDisconnector(params)
    mdis.SetMetalNon(Chem.MolFromSmarts(_METAL_NON))
    return mdis


MDIS = None


def _ensure_mdis():
    global MDIS
    if MDIS is None:
        MDIS = _disconnector()
    return MDIS


def is_valid_smiles_cell(smi: str | None) -> bool:
    if smi is None:
        return False
    text = str(smi).strip()
    if text in FAIL_TOKENS or text.lower() == "nan":
        return False
    return True


def parse_mol(smi: str):
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        return rewrite_z_boron_dative_for_match(mol)
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    if mol is None:
        return None
    mol = rewrite_z_boron_dative_for_match(mol)
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return mol


def rewrite_z_boron_dative_for_match(mol):
    """Reverse Z-type M→B datives to B→M for RDKit valence / ligand matching.

    ILP SMILES keep Cu->B (Z). RDKit then counts B as valence 4 and refuses
    to parse. Reference writes B->Cu. This flip is match-only; stored ILP
    SMILES are unchanged.
    """
    flips = []
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.DATIVE:
            continue
        begin = bond.GetBeginAtom()
        end = bond.GetEndAtom()
        if begin.GetAtomicNum() in TM_ATOMIC_NUMS and end.GetAtomicNum() == 5:
            flips.append((begin.GetIdx(), end.GetIdx()))
    if not flips:
        return mol
    rw = Chem.RWMol(mol)
    for tm_idx, b_idx in flips:
        rw.RemoveBond(tm_idx, b_idx)
        rw.AddBond(b_idx, tm_idx, Chem.BondType.DATIVE)
    out = rw.GetMol()
    try:
        Chem.SanitizeMol(out)
    except Exception:
        out.UpdatePropertyCache(strict=False)
    return out


def tm_atoms(mol):
    return [a for a in mol.GetAtoms() if a.GetAtomicNum() in TM_ATOMIC_NUMS]


def strip_stereo_mol(mol):
    """Drop atom/bond stereo so NBO @SP3 and /C= do not affect comparison."""
    m = Chem.Mol(mol)
    Chem.RemoveStereochemistry(m)
    for atom in m.GetAtoms():
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    for bond in m.GetBonds():
        bond.SetStereo(Chem.BondStereo.STEREONONE)
        bond.SetBondDir(Chem.BondDir.NONE)
    return m


def mol_to_smi(mol) -> str:
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)


def ligand_resonance_set(frag) -> frozenset[str]:
    """Enumerate charge-separated resonance forms for small ligands.

    Large conjugated ligands (porphyrins, etc.) can hang inside
    ResonanceMolSupplier; keep the stereo-stripped SMILES only.
    """
    forms = {mol_to_smi(frag)}
    if frag.GetNumHeavyAtoms() > 24:
        return frozenset(forms)
    try:
        supplier = Chem.ResonanceMolSupplier(
            frag, flags=Chem.ALLOW_CHARGE_SEPARATION, maxStructs=64
        )
        for res in supplier:
            if res is None:
                continue
            forms.add(mol_to_smi(strip_stereo_mol(res)))
    except Exception:
        pass
    return frozenset(forms)


def resonance_multiset_match(sets_a, sets_b) -> bool:
    if len(sets_a) != len(sets_b):
        return False
    unused = list(sets_b)
    for sa in sets_a:
        hit = None
        for i, sb in enumerate(unused):
            if sa & sb:
                hit = i
                break
        if hit is None:
            return False
        unused.pop(hit)
    return True


def _heavy_neighbors(atom):
    return [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]


def _is_pyrazolin_5_ylidene_atom(atom, mol) -> bool:
    """C in a 5-ring with adjacent N–N, bonded to one ring N and one ring C."""
    if not atom.IsInRingSize(5):
        return False
    heavy = _heavy_neighbors(atom)
    if len(heavy) != 2:
        return False
    if {n.GetAtomicNum() for n in heavy} != {6, 7}:
        return False
    for ring in mol.GetRingInfo().AtomRings():
        if atom.GetIdx() not in ring or len(ring) != 5:
            continue
        nitrogens = [
            mol.GetAtomWithIdx(idx)
            for idx in ring
            if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7
        ]
        if len(nitrogens) != 2:
            continue
        if mol.GetBondBetweenAtoms(nitrogens[0].GetIdx(), nitrogens[1].GetIdx()) is None:
            continue
        return True
    return False


def _is_imidazolidin_ylidene_atom(atom, mol) -> bool:
    """C in a saturated 5-ring with two non-adjacent N (imidazolidine).

    Covers imidazolidin-2-ylidene (C between two N is handled separately) and
    C4/C5 connectivity: carbene C bonded to one ring N and one ring C.
    """
    if not atom.IsInRingSize(5):
        return False
    heavy = _heavy_neighbors(atom)
    if len(heavy) != 2:
        return False
    if {n.GetAtomicNum() for n in heavy} != {6, 7}:
        return False
    for ring in mol.GetRingInfo().AtomRings():
        if atom.GetIdx() not in ring or len(ring) != 5:
            continue
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if sum(a.GetAtomicNum() == 6 for a in ring_atoms) != 3:
            continue
        nitrogens = [a for a in ring_atoms if a.GetAtomicNum() == 7]
        if len(nitrogens) != 2:
            continue
        if mol.GetBondBetweenAtoms(nitrogens[0].GetIdx(), nitrogens[1].GetIdx()) is not None:
            continue
        if any(a.GetIsAromatic() for a in ring_atoms):
            continue
        ring_set = set(ring)
        if any(n.GetIdx() not in ring_set for n in heavy):
            continue
        saturated = True
        nring = len(ring)
        for i in range(nring):
            bond = mol.GetBondBetweenAtoms(ring[i], ring[(i + 1) % nring])
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                saturated = False
                break
        if not saturated:
            continue
        return True
    return False


def normalize_nhc_caac_isocyanide(mol):
    """Map ILP vs xyz2mol NHC/CAAC and isocyanide drawings to one form.

    NHC/CAAC/OXA/ADC/PYZ/IMD/SNHC: CH2 or [C-2] between N/N (ring or acyclic ADC),
    N/C(quaternary), N/O (oxazolidin-2-ylidene), N/C in a pyrazolin-5-ylidene
    5-ring (adjacent N–N), N/C in a saturated imidazolidine 5-ring
    (two non-adjacent N), or N/N in a saturated 6/7-membered NHC -> [C]
    Isocyanide: C=N-R or [C-2]=N-R -> [C-]#[N+]R
    """
    rw = Chem.RWMol(mol)
    for atom in rw.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        heavy = _heavy_neighbors(atom)
        if len(heavy) != 1 or heavy[0].GetAtomicNum() != 7:
            continue
        n_atom = heavy[0]
        if len(_heavy_neighbors(n_atom)) != 2:
            continue
        bond = rw.GetBondBetweenAtoms(atom.GetIdx(), n_atom.GetIdx())
        if bond is None:
            continue
        if (
            bond.GetBondType() == Chem.BondType.TRIPLE
            and atom.GetFormalCharge() == -1
            and n_atom.GetFormalCharge() == 1
        ):
            continue
        is_imine = bond.GetBondType() == Chem.BondType.DOUBLE
        is_c2m = atom.GetFormalCharge() == -2
        if not (is_imine or is_c2m):
            continue
        bond.SetBondType(Chem.BondType.TRIPLE)
        atom.SetFormalCharge(-1)
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(True)
        atom.SetNumRadicalElectrons(0)
        n_atom.SetFormalCharge(1)
        n_atom.SetNumExplicitHs(0)
    for atom in rw.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.GetIsAromatic():
            continue
        heavy = _heavy_neighbors(atom)
        if len(heavy) != 2:
            continue
        n_n = sum(n.GetAtomicNum() == 7 for n in heavy)
        n_c = sum(n.GetAtomicNum() == 6 for n in heavy)
        n_o = sum(n.GetAtomicNum() == 8 for n in heavy)
        in_nhc_ring = (
            atom.IsInRingSize(5) or atom.IsInRingSize(6) or atom.IsInRingSize(7)
        )
        if n_n == 2 and (in_nhc_ring or not atom.IsInRing()):
            pass
        elif n_n == 1 and n_c == 1 and in_nhc_ring:
            c_nbr = next(n for n in heavy if n.GetAtomicNum() == 6)
            if (
                len(_heavy_neighbors(c_nbr)) < 3
                and not _is_pyrazolin_5_ylidene_atom(atom, rw)
                and not _is_imidazolidin_ylidene_atom(atom, rw)
            ):
                continue
        elif n_n == 1 and n_o == 1 and in_nhc_ring:
            pass
        else:
            continue
        n_h = atom.GetTotalNumHs()
        chg = atom.GetFormalCharge()
        already = n_h == 0 and chg == 0 and atom.GetNumRadicalElectrons() >= 1
        if not (already or n_h == 2 or chg == -2):
            continue
        atom.SetFormalCharge(0)
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(True)
        atom.SetNumRadicalElectrons(2)
    out = rw.GetMol()
    try:
        Chem.SanitizeMol(out)
    except Exception:
        pass
    return strip_stereo_mol(out)


def ligand_equiv_set(frag) -> frozenset[str]:
    """Resonance forms plus NHC/CAAC/OXA/IMD/ADC/PYZ and isocyanide-normalized forms."""
    forms = set(ligand_resonance_set(frag))
    try:
        norm = normalize_nhc_caac_isocyanide(frag)
        if mol_to_smi(norm) != mol_to_smi(frag):
            forms.update(ligand_resonance_set(norm))
            forms.add(mol_to_smi(norm))
    except Exception:
        pass
    return frozenset(forms)


_BOND_ORDER = {
    Chem.BondType.SINGLE: 1,
    Chem.BondType.DOUBLE: 2,
    Chem.BondType.TRIPLE: 3,
}


def convert_ilp_x_to_charged_dative(mol):
    """Rewrite covalent M–L (ILP X) as dative M←L with anion charge on L.

    ILP SMILES keep X as covalent/neutral (Pd(Cl), phenolate as Pd–O). xyz2mol
    writes all M–L as dative, so disconnect yields [Cl-], [O-], etc. Metal
    oxidation is left unchanged.
    """
    rw = Chem.RWMol(mol)
    edits = []
    for atom in rw.GetAtoms():
        if atom.GetAtomicNum() not in TM_ATOMIC_NUMS:
            continue
        tm_idx = atom.GetIdx()
        for nbr in atom.GetNeighbors():
            bond = rw.GetBondBetweenAtoms(tm_idx, nbr.GetIdx())
            if bond is None or bond.GetBondType() == Chem.BondType.DATIVE:
                continue
            order = _BOND_ORDER.get(
                bond.GetBondType(), int(round(bond.GetBondTypeAsDouble() or 1))
            )
            if order < 1:
                order = 1
            edits.append((tm_idx, nbr.GetIdx(), order))
    for tm_idx, lig_idx, order in edits:
        rw.RemoveBond(tm_idx, lig_idx)
        rw.AddBond(lig_idx, tm_idx, Chem.BondType.DATIVE)
        lig = rw.GetAtomWithIdx(lig_idx)
        lig.SetFormalCharge(int(lig.GetFormalCharge()) - order)
        # Keep explicit hydrides (e.g. PhSiH2–Ni → [SiH2-] not [SiH-]).
        # Charge already accounts for the former covalent M–L bond.
    out = rw.GetMol()
    try:
        Chem.SanitizeMol(out)
    except Exception:
        return None
    try:
        out = Chem.RemoveHs(out)
    except Exception:
        pass
    return out


def _ligand_skeleton_smiles(frag) -> str:
    """Canonical heavy-atom SMILES with formal charges stripped.

    ILP writes X ligands as covalent/neutral (Cl, ArOH) while xyz2mol writes
    anions after disconnect ([Cl-], ArO-). The skeleton still identifies the
    same ligand connectivity.
    """
    rw = Chem.RWMol(frag)
    for idx in range(rw.GetNumAtoms() - 1, -1, -1):
        if rw.GetAtomWithIdx(idx).GetAtomicNum() == 1:
            rw.RemoveAtom(idx)
    for atom in rw.GetAtoms():
        atom.SetFormalCharge(0)
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(True)
    return Chem.MolToSmiles(rw.GetMol(), canonical=True)


def extract_features(smi: str, convert_x: bool = False) -> dict | None:
    """Return ligand canonical SMILES, metal neighbor symbols, TM oxidation."""
    if not is_valid_smiles_cell(smi):
        return None
    mol = parse_mol(str(smi).strip())
    if mol is None:
        return None
    try:
        mol = Chem.RemoveHs(mol)
    except Exception:
        return None
    mol = strip_stereo_mol(mol)

    metals = tm_atoms(mol)
    if len(metals) != 1:
        return None
    tm = metals[0]
    ox = int(tm.GetFormalCharge())
    neighbors = tuple(sorted(n.GetSymbol() for n in tm.GetNeighbors()))

    lig_mol = mol
    if convert_x:
        converted = convert_ilp_x_to_charged_dative(mol)
        if converted is not None:
            lig_mol = converted

    try:
        disconnected = _ensure_mdis().Disconnect(Chem.Mol(lig_mol))
        frags = Chem.GetMolFrags(disconnected, asMols=True, sanitizeFrags=True)
    except Exception:
        return None

    ligands: list[str] = []
    ligands_resonance: list[frozenset[str]] = []
    ligands_equiv: list[frozenset[str]] = []
    ligands_skeleton: list[str] = []
    for frag in frags:
        if any(a.GetAtomicNum() in TM_ATOMIC_NUMS for a in frag.GetAtoms()):
            continue
        try:
            frag = strip_stereo_mol(frag)
            ligands.append(mol_to_smi(frag))
            ligands_resonance.append(ligand_resonance_set(frag))
            ligands_equiv.append(ligand_equiv_set(frag))
            ligands_skeleton.append(_ligand_skeleton_smiles(frag))
        except Exception:
            return None
    return {
        "ligands": tuple(sorted(ligands)),
        "ligands_resonance": tuple(ligands_resonance),
        "ligands_equiv": tuple(ligands_equiv),
        "ligands_skeleton": tuple(sorted(ligands_skeleton)),
        "neighbors": neighbors,
        "ox": ox,
    }


def load_ilp(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("status") != "ok":
                continue
            if str(row.get("rdkit_parse_ok", "")).strip() not in {"1", "True", "true"}:
                continue
            smi = row.get("smiles_ilp", "")
            if is_valid_smiles_cell(smi):
                out[row["IDs"]] = smi
    return out


def load_ref(path: Path) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh):
            out[row["IDs"]] = {key: row.get(col, "") for key, col in REF_COLUMNS.items()}
    return out


def pair_label(a: str, b: str) -> str:
    return f"{a} / {b}"


def compare_pair(feat_a: dict, feat_b: dict) -> dict[str, bool]:
    return {
        "ligand_canonical": feat_a["ligands"] == feat_b["ligands"],
        "ligand_resonance": resonance_multiset_match(
            feat_a["ligands_resonance"], feat_b["ligands_resonance"]
        ),
        "ligand_equiv": resonance_multiset_match(
            feat_a["ligands_equiv"], feat_b["ligands_equiv"]
        ),
        "ligand_canonical_ignore_charge": feat_a["ligands_skeleton"]
        == feat_b["ligands_skeleton"],
        "metal_connectivity": feat_a["neighbors"] == feat_b["neighbors"],
        "metal_oxidation": feat_a["ox"] == feat_b["ox"],
    }


def _extract_one(item: tuple[str, str]) -> tuple[tuple[str, str], dict | None]:
    kind, smi = item
    return (kind, smi), extract_features(smi, convert_x=(kind == "ilp"))


def fmt_frac(num: int, den: int) -> str:
    if den == 0:
        return "0 / 0"
    return f"{num} / {den} ({100.0 * num / den:.1f}%)"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ilp-csv", type=Path, default=DEFAULT_ILP)
    parser.add_argument("--ref-csv", type=Path, default=DEFAULT_REF)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    ilp = load_ilp(args.ilp_csv)
    refs = load_ref(args.ref_csv)
    ids = sorted(set(refs) | set(ilp))
    if args.limit:
        ids = ids[: args.limit]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    detail_path = args.out_dir / "per_structure.csv"
    summary_path = args.out_dir / "summary.csv"

    unique: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for csd_id in ids:
        if csd_id in ilp:
            job = ("ilp", ilp[csd_id])
            if job not in seen:
                seen.add(job)
                unique.append(job)
        for ref_name in REF_COLUMNS:
            ref_smi = refs.get(csd_id, {}).get(ref_name, "")
            if not is_valid_smiles_cell(ref_smi):
                continue
            job = ("ref", ref_smi)
            if job not in seen:
                seen.add(job)
                unique.append(job)

    print(f"extracting features for {len(unique)} unique SMILES ...", flush=True)
    print(f"pairs ({len(PAIRS)}): " + ", ".join(pair_label(a, b) for a, b in PAIRS), flush=True)
    feat_cache: dict[tuple[str, str], dict | None] = {}
    if args.workers <= 1:
        for i, job in enumerate(unique, 1):
            _, feat = _extract_one(job)
            feat_cache[job] = feat
            if i % 5000 == 0:
                print(f"  {i} / {len(unique)}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            for i, (job, feat) in enumerate(pool.map(_extract_one, unique, chunksize=16), 1):
                feat_cache[job] = feat
                if i % 5000 == 0:
                    print(f"  {i} / {len(unique)}", flush=True)

    def feat_of(csd_id: str, method: str):
        if method == "ILP":
            smi = ilp.get(csd_id, "")
            if not smi:
                return None
            return feat_cache.get(("ilp", smi))
        smi = refs.get(csd_id, {}).get(method, "")
        if not is_valid_smiles_cell(smi):
            return None
        return feat_cache.get(("ref", smi))

    tallies = {pair_label(a, b): Counter() for a, b in PAIRS}
    metrics = (
        "ligand_canonical",
        "ligand_resonance",
        "ligand_equiv",
        "ligand_canonical_ignore_charge",
        "metal_connectivity",
        "metal_oxidation",
    )
    detail_fields = [
        "IDs",
        "pair",
        "both_parsed",
        "ligand_canonical",
        "ligand_resonance",
        "ligand_equiv",
        "ligand_canonical_ignore_charge",
        "metal_connectivity",
        "metal_oxidation",
        "ox_a",
        "ox_b",
        "neighbors_a",
        "neighbors_b",
    ]

    n_done = 0
    with detail_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=detail_fields)
        writer.writeheader()
        for csd_id in ids:
            feats = {m: feat_of(csd_id, m) for m in METHODS}
            for a, b in PAIRS:
                label = pair_label(a, b)
                fa, fb = feats[a], feats[b]
                both = fa is not None and fb is not None
                row = {
                    "IDs": csd_id,
                    "pair": label,
                    "both_parsed": int(both),
                    "ligand_canonical": "",
                    "ligand_resonance": "",
                    "ligand_equiv": "",
                    "ligand_canonical_ignore_charge": "",
                    "metal_connectivity": "",
                    "metal_oxidation": "",
                    "ox_a": fa["ox"] if fa else "",
                    "ox_b": fb["ox"] if fb else "",
                    "neighbors_a": " ".join(fa["neighbors"]) if fa else "",
                    "neighbors_b": " ".join(fb["neighbors"]) if fb else "",
                }
                if both:
                    tallies[label]["n_both"] += 1
                    cmp = compare_pair(fa, fb)
                    for key, val in cmp.items():
                        row[key] = int(val)
                        if val:
                            tallies[label][key] += 1
                writer.writerow(row)
            n_done += 1
            if n_done % 5000 == 0:
                print(f"processed {n_done} / {len(ids)}", flush=True)

    pair_headers = [pair_label(a, b) for a, b in PAIRS]
    with summary_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["metric", *pair_headers])
        for metric in metrics:
            writer.writerow(
                [metric]
                + [
                    f"{tallies[h][metric]} / {tallies[h]['n_both']}"
                    for h in pair_headers
                ]
            )
        writer.writerow(
            ["n_both_parsed"] + [str(tallies[h]["n_both"]) for h in pair_headers]
        )

    print(f"structures considered: {len(ids)}")
    print()
    print("metric," + ",".join(pair_headers))
    for metric in metrics:
        cells = [fmt_frac(tallies[h][metric], tallies[h]["n_both"]) for h in pair_headers]
        print(metric + "," + ",".join(cells))
    print()
    print("Wrote", summary_path)
    print("Wrote", detail_path)


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        sys.exit(0)
