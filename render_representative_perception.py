#!/usr/bin/env python3
"""Perceive representative TMCs with Open Babel and xyz2mol_tmc; save 2D images."""

from __future__ import annotations

import json
import re
import subprocess
import sys
import traceback
from collections import defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point3D

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent
XYZ_ROOT = ROOT / "representative_complex"
OUT_ROOT = XYZ_ROOT / "perception_images"
XYZ2MOL_ROOT = Path("/home/zhujingyuan/xyz2mol_tm")

sys.path.insert(0, str(XYZ2MOL_ROOT))
from xyz2mol_tm.huckel_to_smiles.xyz2mol_tmc import get_tmc_mol  # noqa: E402

CHARGE_RE = re.compile(r"_charge_(plus|minus)?(\d+)\.xyz$", re.I)
BOND_LABEL = {
    Chem.BondType.SINGLE: "1",
    Chem.BondType.DOUBLE: "2",
    Chem.BondType.TRIPLE: "3",
    Chem.BondType.AROMATIC: "a",
    Chem.BondType.DATIVE: "dative",
    Chem.BondType.UNSPECIFIED: "?",
}


def parse_charge(path: Path) -> int:
    m = CHARGE_RE.search(path.name)
    if not m:
        return 0
    sign, n = m.group(1), int(m.group(2))
    if sign == "minus":
        return -n
    return n


def csd_id(path: Path) -> str:
    return path.name.split("_charge_")[0]


def mol_from_openbabel(xyz: Path) -> Chem.Mol:
    proc = subprocess.run(
        ["obabel", str(xyz), "-omol", "--gen2d"],
        check=False,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0 or not proc.stdout.strip():
        raise RuntimeError(proc.stderr.strip() or "Open Babel produced no MOL block")
    mol = Chem.MolFromMolBlock(proc.stdout, sanitize=False, removeHs=False)
    if mol is None:
        raise RuntimeError("RDKit failed to parse Open Babel MOL block")
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass
    try:
        mol = Chem.RemoveHs(mol, sanitize=False)
    except Exception:
        rw = Chem.RWMol(mol)
        for idx in sorted(
            (a.GetIdx() for a in rw.GetAtoms() if a.GetAtomicNum() == 1),
            reverse=True,
        ):
            rw.RemoveAtom(idx)
        mol = rw.GetMol()
    # Open Babel XYZ perception writes M–L as MOL type=1. RDKit SanitizeMol
    # rewrites some of those to DATIVE; restore single bonds for drawing.
    dative = {
        Chem.BondType.DATIVE,
        Chem.BondType.DATIVER,
        Chem.BondType.DATIVEL,
    }
    for bond in mol.GetBonds():
        if bond.GetBondType() in dative:
            bond.SetBondType(Chem.BondType.SINGLE)
    mol.UpdatePropertyCache(strict=False)
    return mol


def mol_from_xyz2mol(xyz: Path, charge: int) -> Chem.Mol:
    mol = get_tmc_mol(xyz, charge, with_stereo=False)
    if mol is None:
        raise RuntimeError("xyz2mol_tmc returned None")
    return mol


_TM_NUMS = (
    set(range(21, 31))
    | set(range(39, 49))
    | set(range(72, 81))
    | set(range(104, 113))
    | {57, 89}
)


def _is_tm_atom(atom: Chem.Atom) -> bool:
    return atom.GetAtomicNum() in _TM_NUMS


def _haptic_units(mol: Chem.Mol) -> list[tuple[int, list[int]]]:
    """Connected metal-bound heavy atoms: Cp, arene, allyl, η²-olefin, …"""
    ligand_adj: dict[int, list[int]] = defaultdict(list)
    metal_nbrs: dict[int, set[int]] = defaultdict(set)
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        ai, aj = mol.GetAtomWithIdx(i), mol.GetAtomWithIdx(j)
        ti, tj = _is_tm_atom(ai), _is_tm_atom(aj)
        if ti ^ tj:
            metal, lig = (i, j) if ti else (j, i)
            if mol.GetAtomWithIdx(lig).GetAtomicNum() > 1:
                metal_nbrs[metal].add(lig)
        elif not ti and not tj:
            ligand_adj[i].append(j)
            ligand_adj[j].append(i)
    units: list[tuple[int, list[int]]] = []
    for metal, ligs in metal_nbrs.items():
        seen: set[int] = set()
        for start in ligs:
            if start in seen:
                continue
            stack = [start]
            seen.add(start)
            comp: list[int] = []
            while stack:
                u = stack.pop()
                comp.append(u)
                for v in ligand_adj[u]:
                    if v in ligs and v not in seen:
                        seen.add(v)
                        stack.append(v)
            if len(comp) >= 2:
                units.append((metal, comp))
    return units


def _attach_haptic_centroids(mol: Chem.Mol) -> tuple[Chem.Mol, list[tuple[int, int, list[int]]]]:
    """Replace η contacts with a dummy at each haptic set (depiction only)."""
    units = _haptic_units(mol)
    if not units:
        return mol, []
    rw = Chem.RWMol(mol)
    for metal, comp in units:
        for atom_idx in comp:
            if rw.GetBondBetweenAtoms(metal, atom_idx) is not None:
                rw.RemoveBond(metal, atom_idx)
    records: list[tuple[int, int, list[int]]] = []
    for metal, comp in units:
        dummy = Chem.Atom(0)
        dummy.SetProp("dummyLabel", "")
        dummy.SetNoImplicit(True)
        d_idx = rw.AddAtom(dummy)
        for atom_idx in comp:
            rw.AddBond(d_idx, atom_idx, Chem.BondType.SINGLE)
        rw.AddBond(metal, d_idx, Chem.BondType.SINGLE)
        records.append((d_idx, metal, list(comp)))
    out = rw.GetMol()
    out.UpdatePropertyCache(strict=False)
    return out, records


def _finalize_haptic_centroids(
    mol: Chem.Mol, records: list[tuple[int, int, list[int]]]
) -> Chem.Mol:
    """Park each dummy at the ring centroid; keep one dashed metal–centroid bond."""
    if not records or mol.GetNumConformers() == 0:
        return mol
    conf = mol.GetConformer()
    conf.Set3D(False)
    for d_idx, _metal, comp in records:
        pts = [conf.GetAtomPosition(i) for i in comp]
        cx = sum(p.x for p in pts) / len(pts)
        cy = sum(p.y for p in pts) / len(pts)
        conf.SetAtomPosition(d_idx, Point3D(cx, cy, 0.0))
    rw = Chem.RWMol(mol)
    for d_idx, metal, comp in records:
        for atom_idx in comp:
            if rw.GetBondBetweenAtoms(d_idx, atom_idx) is not None:
                rw.RemoveBond(d_idx, atom_idx)
        bond = rw.GetBondBetweenAtoms(metal, d_idx)
        if bond is not None:
            bond.SetBondType(Chem.BondType.HYDROGEN)
    out = rw.GetMol()
    src = mol.GetConformer()
    dst = Chem.Conformer(out.GetNumAtoms())
    dst.Set3D(False)
    for i in range(out.GetNumAtoms()):
        dst.SetAtomPosition(i, src.GetAtomPosition(i))
    out.RemoveAllConformers()
    out.AddConformer(dst, assignId=True)
    out.UpdatePropertyCache(strict=False)
    return out


def draw_mol(mol: Chem.Mol, title: str, png_path: Path) -> None:
    mol = Chem.Mol(mol)
    for atom in mol.GetAtoms():
        atom.ClearProp("atomNote")
    mol, haptic_recs = _attach_haptic_centroids(mol)
    try:
        rdDepictor.SetPreferCoordGen(True)
    except Exception:
        pass
    try:
        rdDepictor.Compute2DCoords(mol)
    except Exception:
        if mol.GetNumConformers() == 0:
            AllChem.Compute2DCoords(mol)
    mol = _finalize_haptic_centroids(mol, haptic_recs)

    n = max(mol.GetNumAtoms(), 8)
    w = min(1600, 420 + 18 * n)
    h = min(1200, 320 + 14 * n)
    drawer = rdMolDraw2D.MolDraw2DCairo(w, h)
    opts = drawer.drawOptions()
    opts.addAtomIndices = False
    opts.annotationFontScale = 0.7
    opts.legendFontSize = 18
    opts.dummyIsotopeLabels = False
    opts.dummiesAreAttachments = False
    drawer.DrawMolecule(mol, legend=title)
    drawer.FinishDrawing()
    png_path.parent.mkdir(parents=True, exist_ok=True)
    png_path.write_bytes(drawer.GetDrawingText())


def draw_error(message: str, png_path: Path) -> None:
    from PIL import Image, ImageDraw, ImageFont

    img = Image.new("RGB", (900, 360), "white")
    draw = ImageDraw.Draw(img)
    try:
        font = ImageFont.truetype("DejaVuSans.ttf", 18)
    except Exception:
        font = ImageFont.load_default()
    draw.multiline_text((24, 24), message[:1200], fill="black", font=font, spacing=4)
    png_path.parent.mkdir(parents=True, exist_ok=True)
    img.save(png_path)


def mosaic(left: Path, right: Path, out: Path, labels: tuple[str, str]) -> None:
    from PIL import Image, ImageDraw, ImageFont

    im_l = Image.open(left).convert("RGB")
    im_r = Image.open(right).convert("RGB")
    pad, header = 16, 40
    h = max(im_l.height, im_r.height)
    canvas = Image.new("RGB", (im_l.width + im_r.width + 3 * pad, h + header + 2 * pad), "white")
    canvas.paste(im_l, (pad, header + pad))
    canvas.paste(im_r, (2 * pad + im_l.width, header + pad))
    draw = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.truetype("DejaVuSans.ttf", 22)
    except Exception:
        font = ImageFont.load_default()
    draw.text((pad, 10), labels[0], fill="black", font=font)
    draw.text((2 * pad + im_l.width, 10), labels[1], fill="black", font=font)
    out.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out)


def mol_summary(mol: Chem.Mol) -> dict:
    bonds = []
    for b in mol.GetBonds():
        a = b.GetBeginAtom()
        c = b.GetEndAtom()
        bonds.append(
            {
                "a": f"{a.GetSymbol()}{a.GetIdx()}",
                "b": f"{c.GetSymbol()}{c.GetIdx()}",
                "order": BOND_LABEL.get(b.GetBondType(), str(b.GetBondType())),
            }
        )
    charges = [
        {"atom": f"{a.GetSymbol()}{a.GetIdx()}", "fc": a.GetFormalCharge()}
        for a in mol.GetAtoms()
        if a.GetFormalCharge() != 0
    ]
    tm = [
        {"atom": f"{a.GetSymbol()}{a.GetIdx()}", "fc": a.GetFormalCharge()}
        for a in mol.GetAtoms()
        if a.GetAtomicNum() in {
            21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
            39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
            57, 72, 73, 74, 75, 76, 77, 78, 79, 80,
        }
    ]
    try:
        smiles = Chem.MolToSmiles(mol)
    except Exception as exc:
        smiles = f"(MolToSmiles failed: {exc})"
    return {
        "n_atoms": mol.GetNumAtoms(),
        "n_bonds": mol.GetNumBonds(),
        "smiles": smiles,
        "tm": tm,
        "formal_charges": charges,
        "bonds": bonds,
    }


def main() -> None:
    xyz_files = sorted(XYZ_ROOT.glob("*/*.xyz"))
    if not xyz_files:
        raise SystemExit(f"No XYZ files under {XYZ_ROOT}")

    records = []
    for xyz in xyz_files:
        cid = csd_id(xyz)
        charge = parse_charge(xyz)
        rec = {
            "id": cid,
            "xyz": str(xyz.relative_to(ROOT)),
            "charge": charge,
            "category": xyz.parent.name,
        }
        print(f"=== {cid}  charge={charge} ===", flush=True)

        ob_png = OUT_ROOT / "openbabel" / f"{cid}.png"
        tm_png = OUT_ROOT / "xyz2mol_tm" / f"{cid}.png"
        cmp_png = OUT_ROOT / "comparison" / f"{cid}.png"

        try:
            ob_mol = mol_from_openbabel(xyz)
            draw_mol(ob_mol, f"{cid}  Open Babel  q={charge}", ob_png)
            rec["openbabel"] = mol_summary(ob_mol)
            rec["openbabel"]["png"] = str(ob_png.relative_to(ROOT))
            print(f"  Open Babel OK  bonds={ob_mol.GetNumBonds()}", flush=True)
        except Exception as exc:
            rec["openbabel"] = {"error": str(exc), "traceback": traceback.format_exc()}
            draw_error(f"{cid}  Open Babel failed\n{exc}", ob_png)
            print(f"  Open Babel FAIL: {exc}", flush=True)

        try:
            tm_mol = mol_from_xyz2mol(xyz, charge)
            draw_mol(tm_mol, f"{cid}  xyz2mol_tmc  q={charge}", tm_png)
            rec["xyz2mol_tm"] = mol_summary(tm_mol)
            rec["xyz2mol_tm"]["png"] = str(tm_png.relative_to(ROOT))
            print(f"  xyz2mol_tmc OK  {rec['xyz2mol_tm']['smiles'][:80]}", flush=True)
        except Exception as exc:
            rec["xyz2mol_tm"] = {"error": str(exc), "traceback": traceback.format_exc()}
            draw_error(f"{cid}  xyz2mol_tmc failed\n{exc}", tm_png)
            print(f"  xyz2mol_tmc FAIL: {exc}", flush=True)

        try:
            mosaic(ob_png, tm_png, cmp_png, (f"{cid}  Open Babel", f"{cid}  xyz2mol_tmc"))
            rec["comparison_png"] = str(cmp_png.relative_to(ROOT))
        except Exception as exc:
            rec["comparison_error"] = str(exc)

        records.append(rec)

    summary_path = OUT_ROOT / "summary.json"
    summary_path.write_text(json.dumps(records, indent=2, ensure_ascii=False) + "\n")
    print(f"\nWrote {len(records)} records to {summary_path}")


if __name__ == "__main__":
    main()
