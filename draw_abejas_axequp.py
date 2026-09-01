#!/usr/bin/env python3
"""Clearer RDKit depictions of ABEJAS and AXEQUP (ILP SMILES)."""

import csv
from pathlib import Path

import numpy as np
from PIL import Image
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

rdDepictor.SetPreferCoordGen(True)

ROOT = Path(__file__).resolve().parent
ILP_CSV = ROOT / "tmqmg_ilp_benchmark_output" / "per_structure_results.csv"
TM = {
    21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
    39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
    57, 72, 73, 74, 75, 76, 77, 78, 79, 80,
}
PH_SMARTS = Chem.MolFromSmarts("[#15,#21,#22,#23,#24,#25,#26,#27,#28,#29,#30,#39,#40,#41,#42,#43,#44,#45,#46,#47,#48,#57,#72,#73,#74,#75,#76,#77,#78,#79,#80]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1")


def load_ilp_smiles():
    out = {}
    with ILP_CSV.open() as fh:
        for row in csv.DictReader(fh):
            if row["IDs"] in ("ABEJAS", "AXEQUP"):
                out[row["IDs"]] = row["smiles_ilp"]
    return out


def parse_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    try:
        mol = Chem.RemoveHs(mol)
    except Exception:
        pass
    return mol


def strip_haptic_bonds(mol):
    """Drop all M–C bonds to a 5-membered ring so later 2D can use a centroid link."""
    Chem.FastFindRings(mol)
    tm = next(atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in TM)
    carbon_nbrs = [n for n in tm.GetNeighbors() if n.GetAtomicNum() == 6]
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) != 5:
            continue
        haptic = [n for n in carbon_nbrs if n.GetIdx() in ring]
        if len(haptic) < 4:
            continue
        rw = Chem.RWMol(mol)
        for nbr in haptic:
            rw.RemoveBond(tm.GetIdx(), nbr.GetIdx())
        out = rw.GetMol()
        Chem.FastFindRings(out)
        return out
    return mol


def neutralize_dative_and_co(mol):
    rw = Chem.RWMol(mol)
    for bond in rw.GetBonds():
        if bond.GetBondType() == Chem.BondType.DATIVE:
            bond.SetBondType(Chem.BondType.SINGLE)
    for bond in rw.GetBonds():
        a = bond.GetBeginAtom()
        b = bond.GetEndAtom()
        if {a.GetAtomicNum(), b.GetAtomicNum()} != {6, 8}:
            continue
        if {a.GetFormalCharge(), b.GetFormalCharge()} != {-1, 1}:
            continue
        a.SetFormalCharge(0)
        b.SetFormalCharge(0)
        bond.SetBondType(Chem.BondType.TRIPLE)
    return rw.GetMol()


def _is_unsubstituted_ph(mol, match):
    """Keep backbone aryls; only collapse rings whose outer carbons are degree 2."""
    return all(mol.GetAtomWithIdx(idx).GetDegree() == 2 for idx in match[2:])


def collapse_simple_phenyls(mol):
    matches = [
        m for m in mol.GetSubstructMatches(PH_SMARTS) if _is_unsubstituted_ph(mol, m)
    ]
    ipso = []
    to_delete = set()
    for match in matches:
        ipso.append(match[1])
        to_delete.update(match[2:])
    rw = Chem.RWMol(mol)
    for idx in ipso:
        atom = rw.GetAtomWithIdx(idx)
        atom.SetIsAromatic(False)
        atom.SetAtomicNum(6)
        atom.SetFormalCharge(0)
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(True)
        atom.SetProp("atomLabel", "Ph")
        for bond in atom.GetBonds():
            bond.SetIsAromatic(False)
            bond.SetBondType(Chem.BondType.SINGLE)
    for idx in sorted(to_delete, reverse=True):
        rw.RemoveAtom(idx)
    out = rw.GetMol()
    Chem.SanitizeMol(out, catchErrors=True)
    Chem.FastFindRings(out)
    return out


def _ring_tether_idx(mol, ring):
    ring_set = set(ring)
    for idx in ring:
        for nbr in mol.GetAtomWithIdx(idx).GetNeighbors():
            if nbr.GetIdx() in ring_set:
                continue
            if nbr.GetAtomicNum() == 6:
                return idx
    return ring[0]


def add_cp_centroid_bond(mol):
    tm = next(atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in TM)
    Chem.FastFindRings(mol)
    ring = next((r for r in mol.GetRingInfo().AtomRings() if len(r) == 5), None)
    if ring is None:
        return mol
    conf = mol.GetConformer()
    cx = float(np.mean([conf.GetAtomPosition(i).x for i in ring]))
    cy = float(np.mean([conf.GetAtomPosition(i).y for i in ring]))
    radii = [
        float(np.hypot(conf.GetAtomPosition(i).x - cx, conf.GetAtomPosition(i).y - cy))
        for i in ring
    ]
    radius = float(np.mean(radii)) or 1.0
    ring_set = set(ring)
    rest_x = []
    rest_y = []
    for atom in mol.GetAtoms():
        if atom.GetIdx() in ring_set or atom.GetIdx() == tm.GetIdx():
            continue
        pos = conf.GetAtomPosition(atom.GetIdx())
        rest_x.append(pos.x)
        rest_y.append(pos.y)
    ox = float(np.mean(rest_x)) if rest_x else cx
    oy = float(np.mean(rest_y)) if rest_y else cy
    vx, vy = cx - ox, cy - oy
    nrm = float(np.hypot(vx, vy)) or 1.0
    old = conf.GetAtomPosition(tm.GetIdx())
    mx = cx + vx / nrm * (radius + 1.25)
    my = cy + vy / nrm * (radius + 1.25)
    dx, dy = mx - old.x, my - old.y
    conf.SetAtomPosition(tm.GetIdx(), (mx, my, 0.0))
    # Keep carbonyls attached to the metal after the move.
    for nbr in tm.GetNeighbors():
        if nbr.GetAtomicNum() != 6:
            continue
        o_nbrs = [n for n in nbr.GetNeighbors() if n.GetAtomicNum() == 8]
        if not o_nbrs:
            continue
        cpos = conf.GetAtomPosition(nbr.GetIdx())
        conf.SetAtomPosition(nbr.GetIdx(), (cpos.x + dx, cpos.y + dy, 0.0))
        for oxygen in o_nbrs:
            opos = conf.GetAtomPosition(oxygen.GetIdx())
            conf.SetAtomPosition(oxygen.GetIdx(), (opos.x + dx, opos.y + dy, 0.0))
    rw = Chem.RWMol(mol)
    dummy_idx = rw.AddAtom(Chem.Atom(0))
    dummy = rw.GetAtomWithIdx(dummy_idx)
    dummy.SetNoImplicit(True)
    dummy.SetProp("atomLabel", " ")
    rw.AddBond(tm.GetIdx(), dummy_idx, Chem.BondType.SINGLE)
    out = rw.GetMol()
    new_conf = Chem.Conformer(out.GetNumAtoms())
    old = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        new_conf.SetAtomPosition(i, old.GetAtomPosition(i))
    new_conf.SetAtomPosition(dummy_idx, (cx, cy, 0.0))
    out.RemoveAllConformers()
    out.AddConformer(new_conf, assignId=True)
    return out


def clean_labels(mol):
    rw = Chem.RWMol(mol)
    for atom in rw.GetAtoms():
        if atom.GetAtomicNum() in TM:
            atom.SetFormalCharge(0)
            if atom.HasProp("atomLabel"):
                atom.ClearProp("atomLabel")
        if atom.GetAtomicNum() == 0 and atom.HasProp("atomLabel"):
            if atom.GetProp("atomLabel").strip() == "Ph":
                atom.SetProp("atomLabel", "Ph")
        if atom.GetAtomicNum() != 6:
            continue
        nbr_nums = {n.GetAtomicNum() for n in atom.GetNeighbors()}
        if 8 in nbr_nums and nbr_nums & TM:
            atom.SetProp("atomLabel", "C")
    return rw.GetMol()


def spread_ph_labels(mol):
    conf = mol.GetConformer()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:
            continue
        phs = [
            n
            for n in atom.GetNeighbors()
            if n.HasProp("atomLabel") and n.GetProp("atomLabel") == "Ph"
        ]
        others = [
            n
            for n in atom.GetNeighbors()
            if not (n.HasProp("atomLabel") and n.GetProp("atomLabel") == "Ph")
        ]
        if len(phs) < 2 or not others:
            continue
        p = conf.GetAtomPosition(atom.GetIdx())
        o = conf.GetAtomPosition(others[0].GetIdx())
        vx, vy = p.x - o.x, p.y - o.y
        nrm = float(np.hypot(vx, vy)) or 1.0
        px, py = -vy / nrm, vx / nrm
        for i, ph in enumerate(phs):
            sign = -1.0 if i == 0 else 1.0
            conf.SetAtomPosition(
                ph.GetIdx(),
                (
                    p.x + px * sign * 1.35 + vx / nrm * 0.85,
                    p.y + py * sign * 1.35 + vy / nrm * 0.85,
                    0.0,
                ),
            )
    return mol


def prepare(smiles, haptic=False):
    mol = parse_mol(smiles)
    if haptic:
        mol = strip_haptic_bonds(mol)
    mol = neutralize_dative_and_co(mol)
    mol = collapse_simple_phenyls(mol)
    mol = clean_labels(mol)
    rdDepictor.Compute2DCoords(mol)
    try:
        rdDepictor.StraightenDepiction(mol)
    except Exception:
        pass
    if haptic:
        mol = add_cp_centroid_bond(mol)
    return spread_ph_labels(mol)


def crop_white(image, pad=40):
    arr = np.array(image.convert("RGB"))
    mask = arr.min(axis=2) < 250
    ys, xs = np.where(mask)
    box = (
        max(0, int(xs.min()) - pad),
        max(0, int(ys.min()) - pad),
        min(image.width, int(xs.max()) + pad + 1),
        min(image.height, int(ys.max()) + pad + 1),
    )
    return image.crop(box)


def draw_one(mol, path, title):
    drawer = rdMolDraw2D.MolDraw2DCairo(1200, 1000)
    opt = drawer.drawOptions()
    opt.bondLineWidth = 2.6
    opt.padding = 0.08
    opt.minFontSize = 20
    opt.maxFontSize = 28
    opt.legendFontSize = 26
    opt.additionalAtomLabelPadding = 0.16
    opt.explicitMethyl = True
    opt.dummyIsotopeLabels = False
    drawer.DrawMolecule(mol, legend=title)
    drawer.FinishDrawing()
    raw = Path("/tmp") / f"{path.stem}_raw.png"
    raw.write_bytes(drawer.GetDrawingText())
    image = crop_white(Image.open(raw))
    image.save(path)
    print(path.name, "atoms", mol.GetNumAtoms(), "size", image.size)
    return image


def stitch(left, right, dest, gap=36):
    height = max(left.height, right.height)

    def resize_h(im):
        width = int(im.width * height / im.height)
        return im.resize((width, height), Image.Resampling.LANCZOS)

    left, right = resize_h(left), resize_h(right)
    canvas = Image.new("RGB", (left.width + right.width + gap, height), "white")
    canvas.paste(left, (0, 0))
    canvas.paste(right, (left.width + gap, 0))
    canvas.save(dest)
    print(dest.name, canvas.size)


def main():
    smiles = load_ilp_smiles()
    abejas = prepare(smiles["ABEJAS"], haptic=False)
    axequp = prepare(smiles["AXEQUP"], haptic=True)
    im1 = draw_one(abejas, ROOT / "ABEJAS_ilp_rdkit.png", "ABEJAS")
    im2 = draw_one(axequp, ROOT / "AXEQUP_ilp_rdkit.png", "AXEQUP")
    stitch(im1, im2, ROOT / "ABEJAS_AXEQUP_ilp_rdkit.png")


if __name__ == "__main__":
    main()
