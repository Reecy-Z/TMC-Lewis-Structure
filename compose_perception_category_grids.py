#!/usr/bin/env python3
"""Per-category grids: rows = complexes, columns = Open Babel | xyz2mol_tmc | ILP."""

from __future__ import annotations

import importlib.util
import json
import sys
import traceback
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from render_representative_perception import (  # noqa: E402
    OUT_ROOT,
    XYZ_ROOT,
    csd_id,
    draw_error,
    draw_mol,
    parse_charge,
)

CELL_W = 560
CELL_H = 420
LABEL_W = 130
HEADER_H = 52
PAD = 10
COL_TITLES = ("Open Babel", "xyz2mol_tmc", "ILP")


def load_engine():
    spec = importlib.util.spec_from_file_location(
        "lewis_engine_ilp", str(ROOT / "Lewis-engine-ILP.py")
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def ilp_mol(engine, xyz: Path, charge: int) -> Chem.Mol:
    atoms = engine.read_xyz(str(xyz))
    raw = engine.connectivity(atoms)
    aromatic = engine.aromatic_candidate_systems(atoms, raw)
    bonds, lp_out, fc_out = engine.solve_bond_orders(
        atoms, raw, aromatic, mol_charge=charge, metal_adjacency_edges=raw
    )
    bonds, lp_out, fc_out, _ = engine.apply_heterocyclic_carbene_corrections(
        atoms, bonds, lp_out, fc_out, aromatic, raw, mol_charge=charge
    )
    bonds, lp_out, fc_out = engine.apply_eta_covalent_pi_corrections(
        atoms, bonds, lp_out, fc_out, mol_charge=charge, metal_adjacency_edges=raw
    )
    atom_syms = [a[1] for a in atoms]
    coords = [[a[2], a[3], a[4]] for a in atoms]
    fc0 = [fc_out.get(a[0], 0) for a in atoms]
    idx_to_pos = {a[0]: k for k, a in enumerate(atoms)}
    metal_adj_0 = [
        (idx_to_pos[tm], idx_to_pos[lig], ei, ej)
        for tm, lig, ei, ej in raw
        if engine.is_TM(ei) ^ engine.is_TM(ej)
    ]
    dative = engine.infer_dative_ml_pairs_cbc(
        atom_syms,
        coords,
        bonds,
        lp_out,
        fc0,
        metal_adjacency_edges=metal_adj_0,
    )
    mol = engine.ilp_to_rdkit_mol(
        atoms, bonds, fc_out, dative_ml_pairs=dative, edges=raw
    )
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
    return mol


def fit_cell(im: Image.Image, w: int, h: int) -> Image.Image:
    im = im.convert("RGB")
    scale = min(w / im.width, h / im.height, 1.0)
    nw, nh = max(1, int(im.width * scale)), max(1, int(im.height * scale))
    im = im.resize((nw, nh), Image.Resampling.LANCZOS)
    canvas = Image.new("RGB", (w, h), "white")
    canvas.paste(im, ((w - nw) // 2, (h - nh) // 2))
    return canvas


def font(size: int):
    try:
        return ImageFont.truetype("DejaVuSans.ttf", size)
    except Exception:
        return ImageFont.load_default()


def compose_category(title: str, rows: list[tuple[str, Path, Path, Path]], out: Path) -> None:
    n = len(rows)
    width = LABEL_W + 3 * CELL_W + 4 * PAD
    height = HEADER_H + n * (CELL_H + PAD) + PAD
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)
    f_h, f_l = font(22), font(18)
    draw.text((PAD, 12), title, fill="black", font=f_h)
    for i, col in enumerate(COL_TITLES):
        x = LABEL_W + PAD + i * (CELL_W + PAD) + 8
        draw.text((x, 16), col, fill="black", font=f_h)

    for r, (cid, p0, p1, p2) in enumerate(rows):
        y = HEADER_H + r * (CELL_H + PAD)
        draw.text((PAD, y + CELL_H // 2 - 10), cid, fill="black", font=f_l)
        for c, path in enumerate((p0, p1, p2)):
            x = LABEL_W + PAD + c * (CELL_W + PAD)
            if path.is_file():
                cell = fit_cell(Image.open(path), CELL_W, CELL_H)
            else:
                cell = Image.new("RGB", (CELL_W, CELL_H), "white")
                ImageDraw.Draw(cell).text((20, 20), "missing", fill="red", font=f_l)
            canvas.paste(cell, (x, y))
    out.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out, "PNG")
    print(f"wrote {out} ({n} rows)")


def main() -> None:
    skip_ilp = "--skip-ilp" in sys.argv
    engine = None if skip_ilp else load_engine()
    manifest = json.loads((XYZ_ROOT / "manifest.json").read_text())
    by_cat: dict[str, list] = defaultdict(list)
    for rec in manifest:
        by_cat[rec["category_folder"]].append(rec)
    for recs in by_cat.values():
        recs.sort(key=lambda r: r.get("panel", "Z"))

    ilp_dir = OUT_ROOT / "ilp"
    ilp_dir.mkdir(parents=True, exist_ok=True)

    if not skip_ilp:
        for rec in manifest:
            cid = rec["id"]
            xyz = XYZ_ROOT / rec["xyz"]
            if not xyz.is_file():
                matches = list(XYZ_ROOT.glob(f"*/{cid}_*.xyz")) + list(
                    XYZ_ROOT.glob(f"*/{cid}.xyz")
                )
                if not matches:
                    print(f"ILP skip {cid}: xyz missing")
                    continue
                xyz = matches[0]
            charge = rec.get("charge")
            if charge is None:
                charge = parse_charge(xyz)
            png = ilp_dir / f"{cid}.png"
            print(f"ILP {cid} q={charge}", flush=True)
            try:
                mol = ilp_mol(engine, xyz, int(charge))
                draw_mol(mol, f"{cid}  ILP  q={charge}", png)
                print(f"  ILP OK atoms={mol.GetNumAtoms()}", flush=True)
            except Exception as exc:
                draw_error(f"{cid}  ILP failed\n{exc}", png)
                print(f"  ILP FAIL: {exc}\n{traceback.format_exc()}", flush=True)

    cat_titles = {
        "01_werner_type": "Figure 2  Werner-type coordination complexes",
        "02_eta_pi_donor": "Figure 3  eta-coordinated pi-donor ligands",
        "03_metal_carbene": "Figure 4  Metal-carbene complexes",
        "04_additional_motifs": "Figure 5  Additional structural and bonding motifs",
    }
    grid_dir = OUT_ROOT / "category_grids"
    for folder, title in cat_titles.items():
        rows = []
        for rec in by_cat[folder]:
            cid = rec["id"]
            rows.append(
                (
                    cid,
                    OUT_ROOT / "openbabel" / f"{cid}.png",
                    OUT_ROOT / "xyz2mol_tm" / f"{cid}.png",
                    OUT_ROOT / "ilp" / f"{cid}.png",
                )
            )
        compose_category(title, rows, grid_dir / f"{folder}.png")


if __name__ == "__main__":
    main()
