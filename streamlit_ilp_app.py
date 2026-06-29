#!/usr/bin/env python3
"""
Streamlit app for ILP-based Lewis structure analysis (Lewis-engine-ILP backend).

Run:
    streamlit run streamlit_ilp_app.py --server.address 0.0.0.0 --server.port 8501
"""

from __future__ import annotations

import contextlib
import io
import importlib.util
import json
import math
import sys
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components


APP_TITLE = "Orbis QC — TMC Lewis Analyzer"
ROOT_DIR = Path(__file__).resolve().parent
ENGINE_FILE = ROOT_DIR / "Lewis-engine-ILP.py"
VIEWER_FILE = ROOT_DIR / "XYZ_Viewer.html"
# Streamlit iframe height for the embedded XYZ viewer (px). Increase for a taller 3D panel.
VIEWER_IFRAME_HEIGHT = 1100
VIEWER_MIN_HEIGHT_PX = 820

# Built-in XYZ examples (paths relative to repo root).
DEMO_CASES = (
    {
        "id": "abasec",
        "label": "ABASEC (+1)",
        "path": ROOT_DIR / "test" / "ABASEC_charge_plus_1.xyz",
        "charge": 1,
    },
    {
        "id": "abelok",
        "label": "ABELOK (0)",
        "path": ROOT_DIR / "test" / "ABELOK_charge_0.xyz",
        "charge": 0,
    },
    {
        "id": "kadxao",
        "label": "KADXAO (0)",
        "path": ROOT_DIR / "test" / "KADXAO_charge_0.xyz",
        "charge": 0,
    },
)

# d-/f-block metals treated as transition metals for upload validation (minus engine.TM_SET).
_TM_LIKE_ELEMENTS = frozenset({
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
})


class UnsupportedTransitionMetalError(ValueError):
    """XYZ contains a transition-metal element outside Lewis-engine-ILP TM_SET."""


class IlpSolveError(RuntimeError):
    """ILP solver did not return an optimal / integer-feasible solution."""


def inject_styles():
    st.markdown(
        """
<style>
@import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;500;600&family=Sora:wght@300;400;500;600;700&display=swap');
html, body, [class*="css"] { font-family: 'Sora', sans-serif; }
.stApp { background: #f0f4f8; }
[data-testid="stHeader"] { display:none; }
#MainMenu, footer { visibility:hidden; }
section[data-testid="stSidebar"],
[data-testid="stSidebarCollapsedControl"]{
  display:none !important;
}
[data-testid="stAppViewContainer"]{
  margin-left:0 !important;
}

.orbis-topbar{
  background: linear-gradient(135deg,#111827,#1f2937 55%,#374151);
  border-radius: 12px;
  padding: 10px 14px;
  color:#fff;
  margin-bottom: 10px;
  display:flex;
  align-items:center;
  justify-content:space-between;
}
.orbis-left{ display:flex; align-items:center; gap:10px; }
.orbis-topbar .orbis-title{ color:#fff; }
.orbis-topbar .orbis-sub{ color:#fff; opacity:.68; }
.orbis-icon{
  width:34px; height:34px; border-radius:10px;
  background:linear-gradient(135deg,rgba(45,106,79,.5),rgba(27,67,50,.6));
  border:1px solid rgba(116,198,157,.35);
  box-shadow:0 2px 10px rgba(0,0,0,.25),0 0 0 1px rgba(116,198,157,.25);
  display:flex; align-items:center; justify-content:center;
  font-size:16px;
}
.orbis-title{ font-size: 18px; font-weight: 700; line-height: 1.05; letter-spacing:-.2px; }
.orbis-sub{ font-size: 11px; opacity: .68; margin-top: 1px; font-style: italic; }
.status-chip{
  font-family:'JetBrains Mono', monospace;
  font-size:10px;
  padding:4px 10px;
  border-radius:999px;
  border:1px solid rgba(34,197,94,.45);
  color:#bbf7d0;
  background:rgba(22,163,74,.2);
}
.mono-label{ font-family: 'JetBrains Mono', monospace; letter-spacing: .4px; color:#64748b; font-size:12px; text-transform: uppercase; }
.tiny-label{ font-family:'JetBrains Mono', monospace; color:#718096; font-size:10px; text-transform:uppercase; letter-spacing:.5px; margin: 2px 0 6px 0; }
.field-label{
  font-family:'JetBrains Mono', monospace;
  color:#64748b;
  font-size:11px;
  text-transform:uppercase;
  letter-spacing:.4px;
  margin:0 0 4px 2px;
  line-height:1.2;
}

.sidebar-card{
  background:#ffffff;
  border:1px solid #d0dae8;
  border-radius:12px;
  padding:4px 10px 8px 10px;
}
.control-panel-marker,
.control-label-row{ display:none !important; }
[data-testid="stVerticalBlock"]:has(.control-panel-marker) [data-testid="stMarkdown"]:has(.control-panel-marker){
  display:none !important;
  height:0 !important;
  margin:0 !important;
  padding:0 !important;
  overflow:hidden !important;
}
[data-testid="stVerticalBlock"]:has(.control-panel-marker){
  margin-bottom:8px !important;
  gap:0 !important;
  padding-top:0 !important;
}
[data-testid="stVerticalBlock"]:has(.control-panel-marker) [data-testid="stMarkdownContainer"]{
  margin:0 !important;
  padding:0 !important;
}
[data-testid="stVerticalBlock"]:has(.control-panel-marker) [data-testid="stHorizontalBlock"]:has(p.field-label):not(:has([data-testid="stFileUploader"])){
  align-items:flex-end !important;
  margin-bottom:0 !important;
}
[data-testid="stVerticalBlock"]:has(.control-panel-marker) [data-testid="stHorizontalBlock"]:has(p.field-label):not(:has([data-testid="stFileUploader"])) > [data-testid="column"]{
  display:flex !important;
  align-items:flex-end !important;
  justify-content:flex-start !important;
}
[data-testid="stVerticalBlock"]:has(.control-panel-marker) [data-testid="stHorizontalBlock"]:has(p.field-label):not(:has([data-testid="stFileUploader"])) p.field-label{
  margin:0 0 4px 2px !important;
}
[data-testid="stVerticalBlock"]:has(.control-panel-marker) [data-testid="stHorizontalBlock"]:has([data-testid="stFileUploader"]){
  margin-top:0 !important;
  margin-bottom:10px !important;
  align-items:stretch !important;
}
[data-testid="stVerticalBlock"]:has(.control-panel-marker) [data-testid="stMarkdown"] p.tiny-label{
  margin:2px 0 6px 0 !important;
}
[data-testid="stVerticalBlock"]:has(.control-panel-marker) [data-testid="stHorizontalBlock"]:has([data-testid="stFileUploader"]) > [data-testid="column"]{
  display:flex !important;
  flex-direction:column !important;
  justify-content:flex-end !important;
}
[data-testid="stVerticalBlock"]:has(.control-panel-marker) [data-testid="stHorizontalBlock"]:has([data-testid="stFileUploader"]) [data-testid="stFileUploader"]{
  margin:0 !important;
  flex:0 0 40px !important;
}
[data-testid="stVerticalBlock"]:has(.control-panel-marker) [data-testid="stHorizontalBlock"]:has([data-testid="stFileUploader"]) [data-testid="stFileUploader"] > div{
  gap:0 !important;
  height:100% !important;
}
[data-testid="stVerticalBlock"]:has(.control-panel-marker) [data-testid="stHorizontalBlock"]:has([data-testid="stFileUploader"]) [data-testid="stFileUploader"] section{
  box-sizing:border-box !important;
  height:40px !important;
  min-height:40px !important;
  max-height:40px !important;
  margin:0 !important;
  padding:0 8px !important;
  display:flex !important;
  align-items:center !important;
  justify-content:flex-start !important;
  overflow:hidden !important;
}
[data-testid="stVerticalBlock"]:has(.control-panel-marker) [data-testid="stHorizontalBlock"]:has([data-testid="stFileUploader"]) [data-testid="stFileUploaderDropzone"]{
  width:100% !important;
  height:100% !important;
  min-height:0 !important;
  margin:0 !important;
  padding:0 !important;
  display:flex !important;
  align-items:center !important;
  gap:0 !important;
}
[data-testid="stVerticalBlock"]:has(.control-panel-marker) [data-testid="stHorizontalBlock"]:has([data-testid="stFileUploader"]) [data-testid="stFileUploaderDropzoneInstructions"],
[data-testid="stVerticalBlock"]:has(.control-panel-marker) [data-testid="stHorizontalBlock"]:has([data-testid="stFileUploader"]) [data-testid="stFileUploader"] small{
  display:none !important;
}
[data-testid="stVerticalBlock"]:has(.control-panel-marker) [data-testid="stHorizontalBlock"]:has([data-testid="stFileUploader"]) [data-testid="stFileUploader"] button{
  min-height:30px !important;
  height:30px !important;
  padding:0 10px !important;
  font-size:12px !important;
  line-height:1 !important;
  margin:0 !important;
}
[data-testid="stVerticalBlock"]:has(.control-panel-marker) [data-testid="stHorizontalBlock"]:has([data-testid="stFileUploader"]) [data-testid="column"]:nth-child(2) [data-testid="stNumberInput"] div[data-baseweb="input"]{
  min-height:40px !important;
  height:40px !important;
}
[data-testid="stVerticalBlock"]:has(.control-panel-marker) [data-testid="stHorizontalBlock"]:has([data-testid="stFileUploader"]) [data-testid="column"]:nth-child(3) .stButton > button{
  min-height:40px !important;
  height:40px !important;
}
.result-card{
  background:#ffffff;
  border:1px solid #d0dae8;
  border-radius:12px;
  padding:2px 8px 10px 8px;
}
.section-head{
  border:1px solid #e2eaf4;
  background:#f7f9fc;
  border-radius:9px;
  padding:7px 10px;
  font-family:'JetBrains Mono', monospace;
  color:#475569;
  font-size:11px;
  text-transform:uppercase;
  letter-spacing:.5px;
  margin:6px 0 7px 0;
}

div[data-testid="stMetric"]{
  border: 1px solid #d0dae8; border-radius: 10px; background: #fff; padding: 8px 10px;
}
div[data-testid="stMetricLabel"] p{
  font-family:'JetBrains Mono', monospace;
  text-transform:uppercase;
  font-size:11px;
  color:#64748b;
}
div[data-testid="stFileUploader"] section{
  border: 1.5px dashed #d0dae8; border-radius: 10px;
  padding: 0.15rem 0.35rem;
  background:#f8fafc;
}
div[data-testid="stFileUploader"] small{ color:#718096 !important; }
div[data-testid="stFileUploaderDropzone"]{ padding: 0.1rem 0.25rem; }
div[data-baseweb="input"] input { font-family:'JetBrains Mono', monospace; }
div[data-baseweb="input"]{ min-height:34px; }

div.stButton > button {
  background: linear-gradient(135deg,#1f2937,#374151);
  color:#fff;
  border:1px solid #1f2937;
  border-radius:9px;
  font-weight:600;
  min-height: 38px;
  font-size: 13px;
}
div.stButton > button:hover {
  border-color:#111827;
  background: linear-gradient(135deg,#374151,#1f2937);
}

div[data-testid="stCodeBlock"]{
  border:1px solid #d0dae8;
  border-radius:10px;
}
[data-testid="stCustomComponentV1"]{
  margin-bottom:0 !important;
}
</style>
        """,
        unsafe_allow_html=True,
    )


@st.cache_resource(show_spinner=False)
def load_engine():
    spec = importlib.util.spec_from_file_location("lewis_engine_ilp", str(ENGINE_FILE))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load engine from {ENGINE_FILE}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_xyz_text(xyz_text: str):
    lines = [ln.strip() for ln in xyz_text.replace("\r\n", "\n").split("\n")]
    while lines and not lines[-1]:
        lines.pop()
    if len(lines) < 2:
        raise ValueError("XYZ content is too short.")
    n = int(lines[0])
    atom_lines = lines[2 : 2 + n]
    if len(atom_lines) != n:
        raise ValueError(f"Expected {n} atom lines, got {len(atom_lines)}.")
    atoms, coords = [], []
    for ln in atom_lines:
        parts = ln.split()
        if len(parts) < 4:
            raise ValueError(f"Malformed atom line: {ln}")
        raw = parts[0]
        sym = raw[0].upper() + raw[1:].lower() if len(raw) > 1 else raw.upper()
        atoms.append(sym)
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return atoms, coords


def unsupported_transition_metals(engine, atoms: list[str]) -> list[str]:
    """Symbols present in XYZ that look like TMs but are not in engine.TM_SET."""
    tm_set = getattr(engine, "TM_SET", frozenset())
    return sorted({sym for sym in atoms if sym in _TM_LIKE_ELEMENTS and sym not in tm_set})


def validate_transition_metals(engine, atoms: list[str]) -> None:
    bad = unsupported_transition_metals(engine, atoms)
    if not bad:
        return
    supported = ", ".join(sorted(getattr(engine, "TM_SET", ())))
    lines = []
    for sym in bad:
        idxs = [str(i + 1) for i, a in enumerate(atoms) if a == sym]
        lines.append(f"**{sym}** (atom index: {', '.join(idxs)})")
    detail = "\n".join(f"- {line}" for line in lines)
    raise UnsupportedTransitionMetalError(
        f"Unsupported transition metal(s) in XYZ: {', '.join(bad)}.\n\n"
        f"{detail}\n\n"
        f"Supported transition metals: {supported}"
    )


def show_unsupported_tm_error(exc: UnsupportedTransitionMetalError) -> None:
    st.error("Unsupported transition metal in uploaded XYZ")
    st.markdown(str(exc))
    st.info(
        "Replace the metal center with a supported transition metal, or remove it from the "
        "structure before running ILP analysis."
    )


def show_ilp_solve_error(
    exc: IlpSolveError | RuntimeError,
    *,
    atoms=None,
    coords=None,
    raw_edges=None,
) -> None:
    st.error(
        "ILP could not find a valid Lewis structure — only initial connectivity is shown below."
    )
    st.warning(
        "Try adjusting the **total molecular charge** (Mol Charge) and/or the **XYZ geometry** "
        "(bond lengths, connectivity, metal–ligand contacts), then run the analysis again."
    )
    if atoms is not None and coords is not None and raw_edges is not None:
        show_connectivity_preview(atoms, coords, raw_edges, show_captions=False)


def run_aromatic_workflow_in_memory(engine, atoms, coords, mol_charge):
    """
    Streamlit-compatible wrapper for Lewis-engine-ILP.py aromatic ILP workflow.
    Returns bo/lp/fc in the same shapes as find_best_lewis().
    """
    atoms_packed = [
        (i + 1, atoms[i], coords[i][0], coords[i][1], coords[i][2])
        for i in range(len(atoms))
    ]
    raw_edges = engine.connectivity(atoms_packed)
    aromatic_systems = engine.aromatic_candidate_systems(atoms_packed, raw_edges)

    try:
        bonds, lp_out, fc_out = engine.solve_bond_orders(
            atoms_packed,
            raw_edges,
            aromatic_systems,
            mol_charge=int(mol_charge),
            metal_adjacency_edges=raw_edges,
        )
        bonds, lp_out, fc_out, _carbene_labels = engine.apply_heterocyclic_carbene_corrections(
            atoms_packed,
            bonds,
            lp_out,
            fc_out,
            aromatic_systems,
            raw_edges,
            mol_charge=int(mol_charge),
        )
        bonds, lp_out, fc_out = engine.apply_eta_covalent_pi_corrections(
            atoms_packed,
            bonds,
            lp_out,
            fc_out,
            mol_charge=int(mol_charge),
            metal_adjacency_edges=raw_edges,
        )
    except RuntimeError as exc:
        if "ILP failed" in str(exc):
            raise IlpSolveError(str(exc)) from exc
        raise
    bo = {(i - 1, j - 1): o for i, j, o in bonds}
    lp = {i - 1: v for i, v in lp_out.items() if v > 0}
    fc = [0] * len(atoms)
    for i in range(len(atoms)):
        fc[i] = fc_out.get(i + 1, 0)
    stats = {
        "mode": "aromatic-aware",
        "aromatic_system_count": len(aromatic_systems),
        "aromatic_systems": [list(s) for s in aromatic_systems],
    }
    coords = [[a[2], a[3], a[4]] for a in atoms_packed]
    idx_to_pos = {a[0]: k for k, a in enumerate(atoms_packed)}
    metal_adj_0 = [
        (idx_to_pos[tm], idx_to_pos[lig], ei, ej)
        for tm, lig, ei, ej in raw_edges
        if engine.is_TM(ei) ^ engine.is_TM(ej)
    ]
    dative_ml_pairs = engine.infer_dative_ml_pairs_cbc(
        atoms,
        coords,
        bonds,
        lp_out,
        fc,
        metal_adjacency_edges=metal_adj_0,
    )
    return bo, lp, fc, stats, aromatic_systems, bonds, fc_out, dative_ml_pairs, raw_edges


def format_choose_block(bo: dict, lp: dict) -> str:
    out = io.StringIO()
    print(" $CHOOSE", file=out)
    lone_items = [(i + 1, v) for i, v in sorted(lp.items()) if v > 0]
    if lone_items:
        print("   LONE " + " ".join(f"{idx} {v}" for idx, v in lone_items) + " END", file=out)
    sym = {1: "S", 2: "D", 3: "T"}
    tokens = [f"{sym.get(o, 'S')} {i+1} {j+1}" for (i, j), o in sorted(bo.items())]
    line = "   BOND"
    for tok in tokens:
        if len(line) + 1 + len(tok) > 72:
            print(line, file=out)
            line = "       " + tok
        else:
            line += " " + tok
    print(line + " END", file=out)
    print(" $END", file=out)
    return out.getvalue()


def metal_adjacency_array_indices(backend, atoms, coords):
    """TM–ligand edges as 0-based array indices (for CBC dative M–C LP)."""
    atoms_packed = [
        (i + 1, atoms[i], coords[i][0], coords[i][1], coords[i][2])
        for i in range(len(atoms))
    ]
    idx_to_pos = {i + 1: i for i in range(len(atoms))}
    out = []
    for i, j, ei, ej in backend.connectivity(atoms_packed):
        if not (backend.is_TM(ei) ^ backend.is_TM(ej)):
            continue
        tm, lig = (i, j) if backend.is_TM(ei) else (j, i)
        out.append((idx_to_pos[tm], idx_to_pos[lig], ei, ej))
    return out


def format_octet_report(
    engine,
    atoms,
    coords,
    bo,
    lp,
    fc,
    *,
    metal_adjacency_edges=None,
) -> str:
    """Octet check + non-metal formal charges (same text as Lewis-engine-ILP CLI)."""
    lp_full = {i: lp.get(i, 0) for i in range(len(atoms))}
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        engine.print_octet_report(
            atoms,
            bo,
            lp_full,
            fc,
            metal_adjacency_edges=metal_adjacency_edges,
            coords=coords,
        )
    return out.getvalue().strip()


def format_cbc_report(
    backend, atoms, coords, bo, lp, fc, charge, *, metal_adjacency_edges=None
):
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        cbc_bundle = backend.classify_cbc_ligands(
            atoms, coords, bo, lp, fc, charge, metal_adjacency_edges=metal_adjacency_edges
        )
        backend.print_cbc_report(
            atoms, coords, bo, lp, fc, charge,
            metal_adjacency_edges=metal_adjacency_edges, cbc_bundle=cbc_bundle,
        )
    return out.getvalue().strip(), cbc_bundle


def format_tm_oxidation_report(
    engine,
    atoms,
    fc,
    mol_charge: int,
    bo=None,
    *,
    cbc_interaction_records=None,
    cbc_neighbor_types=None,
) -> str:
    """Same text as Lewis-engine-ILP.print_tm_oxidation_sigma_report."""
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        engine.print_tm_oxidation_sigma_report(
            atoms,
            fc,
            mol_charge,
            bo=bo,
            cbc_interaction_records=cbc_interaction_records,
            cbc_neighbor_types=cbc_neighbor_types,
        )
    return out.getvalue().strip()


def format_tm_oxidation_summary(
    engine,
    atoms,
    fc,
    mol_charge: int,
    bo=None,
    *,
    cbc_interaction_records=None,
    cbc_neighbor_types=None,
) -> str:
    """One-line TM oxidation summary for the Streamlit UI."""
    common = engine.TM_COMMON_OXIDATION_STATES
    results = engine.tm_oxidation_from_charge_balance(
        atoms,
        fc,
        mol_charge,
        bo=bo,
        cbc_interaction_records=cbc_interaction_records,
        cbc_neighbor_types=cbc_neighbor_types,
    )
    lines: list[str] = []
    for m in sorted(results):
        sym = atoms[m]
        ox_state, _, _ = results[m]
        label = f"{sym}{m + 1}"
        ref = common.get(sym)
        if ref is None:
            status = "—  (no TM_COMMON_OXIDATION_STATES entry)"
        elif ox_state in ref:
            status = f"√  (common for {sym}: {ref})"
        else:
            status = (
                f"⚠ WARNING: oxidation state {ox_state} "
                f"not in common oxidation states {ref}"
            )
        lines.append(f"{label}:  oxidation state = {ox_state} {status}")
    return "\n".join(lines)


def build_viewer_payload(
    backend, atoms, coords, bo, lp, fc, charge, *, metal_adjacency_edges=None
):
    results, _neighbor_cbc = backend.classify_cbc_ligands(
        atoms, coords, bo, lp, fc, charge, metal_adjacency_edges=metal_adjacency_edges
    )

    l_pairs = set()
    for metal_idx in sorted(results):
        if not backend.is_TM(atoms[metal_idx]):
            continue
        for atom_tuple, cbc_char in results[metal_idx]:
            if cbc_char == "L" and len(atom_tuple) == 1:
                lig = atom_tuple[0]
                l_pairs.add((min(metal_idx, lig), max(metal_idx, lig)))

    dative_bonds = []
    for metal_idx in sorted(results):
        if not backend.is_TM(atoms[metal_idx]):
            continue
        for atom_tuple, cbc_char in results[metal_idx]:
            if cbc_char != "L":
                continue
            donor_atoms = list(atom_tuple)
            pi_bond_pair = None
            if len(atom_tuple) == 2:
                a, b = atom_tuple
                if bo.get((min(a, b), max(a, b)), 0) >= 2:
                    pi_bond_pair = [a, b]
            dative_bonds.append(
                {
                    "donor_atoms": donor_atoms,
                    "acceptor": metal_idx,
                    "pi_bond_pair": pi_bond_pair,
                }
            )

    bonds_list = []
    bo_keys = set()
    for (i, j), order in sorted(bo.items()):
        key = (min(i, j), max(i, j))
        bo_keys.add(key)
        if key in l_pairs:
            continue
        bonds_list.append({"i": i, "j": j, "order": order})

    x_counts = {}
    for metal_idx in sorted(results):
        if not backend.is_TM(atoms[metal_idx]):
            continue
        for atom_tuple, cbc_char in results[metal_idx]:
            if cbc_char not in ("X", "Z") or len(atom_tuple) != 1:
                continue
            lig = atom_tuple[0]
            if lig == metal_idx:
                continue
            key = (min(metal_idx, lig), max(metal_idx, lig))
            if key not in bo_keys:
                x_counts[key] = x_counts.get(key, 0) + 1

    for (i, j), count in x_counts.items():
        bonds_list.append({"i": i, "j": j, "order": min(count, 3)})

    return {
        "atoms": atoms,
        "coords": coords,
        "bonds": bonds_list,
        "dative_bonds": dative_bonds,
        "formal_charges": list(fc),
    }


def build_connectivity_viewer_payload(atoms, coords, raw_edges) -> dict:
    """Step-1 connectivity only: every raw edge as order=1; M–L as normal bond lines."""
    bonds_list = []
    seen: set[tuple[int, int]] = set()
    for i, j, _ei, _ej in raw_edges:
        a, b = int(i) - 1, int(j) - 1
        key = (min(a, b), max(a, b))
        if key in seen:
            continue
        seen.add(key)
        bonds_list.append({"i": a, "j": b, "order": 1})
    return {
        "atoms": atoms,
        "coords": coords,
        "bonds": bonds_list,
        "dative_bonds": [],
        "formal_charges": [0] * len(atoms),
        "view_mode": "connectivity",
    }


def show_connectivity_preview(
    atoms,
    coords,
    raw_edges,
    *,
    show_captions: bool = True,
    collapsible: bool = False,
) -> None:
    connectivity_payload = build_connectivity_viewer_payload(atoms, coords, raw_edges)

    def render_body() -> None:
        if show_captions:
            st.caption(
                "Distance-based `connectivity()` only: each contact is a single bond line, "
                "including M–L. No ILP bond orders, CBC, or dative arrows."
            )
            n_conn = len(connectivity_payload["bonds"])
            st.caption(f"{n_conn} connectivity edge(s) from Lewis-engine step 1.")
        show_3d_preview(connectivity_payload)

    if collapsible:
        with st.expander("Connectivity Preview", expanded=False):
            render_body()
    else:
        st.markdown(
            '<div class="section-head">Connectivity Preview</div>',
            unsafe_allow_html=True,
        )
        render_body()


def show_3d_preview(payload: dict):
    if not VIEWER_FILE.exists():
        st.error(f"Required viewer file missing: `{VIEWER_FILE.name}`")
        return

    html_text = VIEWER_FILE.read_text(encoding="utf-8")
    payload_literal = json.dumps(payload)
    inject_js = f"""
<style>
  html, body {{
    width: 100% !important;
    height: 100% !important;
    margin: 0 !important;
    overflow: hidden !important;
  }}
  #app {{
    width: 100% !important;
    height: 100% !important;
  }}
  #header {{
    display: none !important;
  }}
  #app {{
    grid-template-columns: 1fr !important;
    grid-template-rows: auto 1fr var(--status-h) !important;
    grid-template-areas: "panel" "viewer" "status" !important;
  }}
  /* Keep original viewer controls, hide only upload/python bridge UI */
  #py-init-bar, #file-strip {{
    display: none !important;
  }}
  #panel {{
    border-right: none !important;
    border-bottom: 1px solid var(--border) !important;
    max-height: none !important;
    overflow: visible !important;
  }}
  #panel-body {{
    overflow: visible !important;
  }}
  #viewer {{
    position: relative !important;
    min-height: {VIEWER_MIN_HEIGHT_PX}px !important;
  }}
  #plot {{
    width: 100% !important;
    height: 100% !important;
  }}
</style>
<script>
(function() {{
  const payload = {payload_literal};

  function prepEmbeddedMode() {{
    const fileInfo = document.getElementById('file-info');
    if (fileInfo) fileInfo.textContent = '';
    const visInfo = document.getElementById('vis-info');
    const isConn = payload.view_mode === 'connectivity';
    if (visInfo) {{
      visInfo.textContent = isConn
        ? 'Step-1 connectivity (all single bonds; M–L shown as lines).'
        : 'ILP Lewis structure loaded from Streamlit.';
    }}
    const dativeSection = document.getElementById('dative-section');
    if (dativeSection) dativeSection.style.display = isConn ? 'none' : '';
    const placeholderHint = document.querySelector('#viewer-placeholder .placeholder-hint');
    if (placeholderHint) placeholderHint.textContent = 'Structure is provided by Streamlit ILP backend.';
  }}

  function setDataAndRender() {{
    if (typeof state === 'undefined' || typeof visualize !== 'function') return false;
    state.atoms = payload.atoms || [];
    state.coords = payload.coords || [];
    state.bonds = payload.bonds || [];
    state.dativeBonds = payload.dative_bonds || [];
    state.loaded = true;
    state._newFile = true;

    const placeholder = document.getElementById('viewer-placeholder');
    if (placeholder) placeholder.style.display = 'none';
    const modeChip = document.getElementById('mode-chip');
    if (modeChip) {{
      modeChip.textContent = payload.view_mode === 'connectivity'
        ? `Connectivity · ${{state.atoms.length}} atoms`
        : `Lewis · ${{state.atoms.length}} atoms`;
      modeChip.className = 'ok';
    }}
    const fileInfo = document.getElementById('file-info');
    if (fileInfo) {{
      fileInfo.textContent = `${{state.bonds.length}} bonds`
        + (payload.view_mode === 'connectivity' ? ' (step 1)' : ` · ${{state.dativeBonds.length}} dative`);
    }}
    if (typeof updateDativeInfo === 'function') updateDativeInfo();
    visualize();
    setTimeout(() => {{
      window.dispatchEvent(new Event('resize'));
      if (window.Plotly && typeof Plotly.Plots?.resize === 'function') {{
        const p = document.getElementById('plot');
        if (p) Plotly.Plots.resize(p);
      }}
    }}, 60);
    return true;
  }}

  function boot() {{
    prepEmbeddedMode();
    if (!setDataAndRender()) {{
      const timer = setInterval(() => {{
        prepEmbeddedMode();
        if (setDataAndRender()) clearInterval(timer);
      }}, 100);
    }}
  }}
  boot();
}})();
</script>
"""
    components.html(
        html_text.replace("</body>", inject_js + "\n</body>"),
        height=VIEWER_IFRAME_HEIGHT,
        scrolling=False,
    )


def _xyz_export_stem(xyz_file, *, demo_xyz_name: str | None = None) -> str:
    """Base filename for download (upload or quick-example)."""
    if xyz_file is not None:
        return Path(xyz_file.name).stem
    if demo_xyz_name:
        return Path(demo_xyz_name).stem
    return "ilp_result"


def _load_pending_demo_case(case: dict) -> bool:
    """Apply a queued demo before any widget with key ``mol_charge`` is drawn."""
    path = Path(case["path"])
    if not path.is_file():
        st.session_state.demo_load_error = f"Example file not found: {path}"
        return False
    st.session_state.demo_xyz_text = path.read_text(encoding="utf-8")
    st.session_state.demo_xyz_name = path.name
    st.session_state.mol_charge = int(case["charge"])
    st.session_state.pop("demo_load_error", None)
    return True


def run_analyzer_app() -> None:
    pending_demo = st.session_state.pop("pending_demo_case", None)
    if pending_demo is not None:
        _load_pending_demo_case(pending_demo)

    st.markdown(
        """
        <div class="orbis-topbar">
          <div class="orbis-left">
            <div class="orbis-icon">✦</div>
            <div>
              <div class="orbis-title">Orbis QC</div>
              <div class="orbis-sub">XYZ Structure Viewer + Lewis Analysis</div>
            </div>
          </div>
          <div class="status-chip">XYZ · Ready</div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    engine = load_engine()

    with st.container():
        st.markdown('<span class="control-panel-marker"></span>', unsafe_allow_html=True)
        lh1, lh2, lh3 = st.columns([2.4, 0.9, 1.2], gap="small")
        with lh1:
            st.markdown('<p class="field-label">XYZ</p>', unsafe_allow_html=True)
        with lh2:
            st.markdown('<p class="field-label">Charge</p>', unsafe_allow_html=True)
        with lh3:
            st.markdown('<p class="field-label">Action</p>', unsafe_allow_html=True)

        i1, i2, i3 = st.columns([2.4, 0.9, 1.2], gap="small")
        with i1:
            xyz_file = st.file_uploader("XYZ", type=["xyz"], label_visibility="collapsed")
            if xyz_file is not None:
                st.session_state.pop("demo_xyz_text", None)
                st.session_state.pop("demo_xyz_name", None)
        with i2:
            mol_charge = st.number_input(
                "Charge",
                step=1,
                format="%d",
                label_visibility="collapsed",
                key="mol_charge",
                help="Total molecular charge (e.g. −1 for anions, +1 for cations).",
            )
        with i3:
            has_xyz = xyz_file is not None or bool(st.session_state.get("demo_xyz_text"))
            run_btn = st.button(
                "Run Analysis",
                type="primary",
                disabled=not has_xyz,
                use_container_width=True,
            )

        st.markdown('<p class="tiny-label">Quick examples</p>', unsafe_allow_html=True)
        ex1, ex2, ex3 = st.columns(3, gap="small")
        for col, case in zip((ex1, ex2, ex3), DEMO_CASES):
            with col:
                if st.button(
                    case["label"],
                    key=f"demo_{case['id']}",
                    use_container_width=True,
                ):
                    st.session_state.pending_demo_case = case
                    st.rerun()
        if st.session_state.get("demo_load_error"):
            st.error(st.session_state.demo_load_error)
        if st.session_state.get("demo_xyz_name") and xyz_file is None:
            st.caption(f"Loaded example: `{st.session_state.demo_xyz_name}`")

    if not has_xyz:
        st.info("Upload an `.xyz` file or pick a quick example above to preview and analyze.")
        return

    try:
        if xyz_file is not None:
            xyz_text = xyz_file.getvalue().decode("utf-8", errors="replace")
        else:
            xyz_text = st.session_state.demo_xyz_text
        atoms, coords = parse_xyz_text(xyz_text)
        validate_transition_metals(engine, atoms)
        engine.validate_atom_symbols(atoms)

        atoms_packed = [
            (i + 1, atoms[i], coords[i][0], coords[i][1], coords[i][2])
            for i in range(len(atoms))
        ]
        raw_edges = engine.connectivity(atoms_packed)
        aromatic_systems_preview = engine.aromatic_candidate_systems(atoms_packed, raw_edges)
        st.markdown('<div class="section-head">Aromatic Candidate Systems</div>', unsafe_allow_html=True)
        if not aromatic_systems_preview:
            st.caption("No planar aromatic candidate systems detected.")
        else:
            for i, s in enumerate(aromatic_systems_preview, start=1):
                st.caption(f"System {i}: atoms {list(s)}")
        if not run_btn:
            return

        with st.spinner("Running ILP solver..."):
            (
                bo,
                lp,
                fc,
                stats,
                aromatic_systems,
                bonds_1b,
                fc_out_1b,
                dative_ml_pairs,
                raw_edges,
            ) = run_aromatic_workflow_in_memory(
                engine, atoms, coords, int(mol_charge)
            )
            choose_block = format_choose_block(bo, lp)
            metal_adj_0 = metal_adjacency_array_indices(engine, atoms, coords)
            octet_report = format_octet_report(
                engine,
                atoms,
                coords,
                bo,
                lp,
                fc,
                metal_adjacency_edges=metal_adj_0,
            )
            cbc_report, cbc_bundle = format_cbc_report(
                engine, atoms, coords, bo, lp, fc, int(mol_charge),
                metal_adjacency_edges=metal_adj_0,
            )
            cbc_interaction_records, _cbc_neighbor_types = cbc_bundle
            ox_report = format_tm_oxidation_report(
                engine, atoms, fc, int(mol_charge), bo=bo,
                cbc_interaction_records=cbc_interaction_records,
            )
            ox_summary = format_tm_oxidation_summary(
                engine, atoms, fc, int(mol_charge), bo=bo,
                cbc_interaction_records=cbc_interaction_records,
            )
            viewer_payload = build_viewer_payload(
                engine, atoms, coords, bo, lp, fc, int(mol_charge),
                metal_adjacency_edges=metal_adj_0,
            )

        show_3d_preview(viewer_payload)

        show_connectivity_preview(atoms, coords, raw_edges, collapsible=True)

        st.markdown('<div class="result-card">', unsafe_allow_html=True)
        c1, c2, c3 = st.columns(3)
        c1.metric("Atoms", len(atoms))
        c2.metric("Charge Input", int(mol_charge))
        c3.metric("Lewis Bonds", len(bo))

        st.markdown('<div class="section-head">Aromatic Summary</div>', unsafe_allow_html=True)
        st.caption(f"Detected aromatic systems: {len(aromatic_systems)}")
        for i, s in enumerate(aromatic_systems, start=1):
            st.caption(f"System {i}: atoms {list(s)}")

        st.markdown('<div class="section-head">$CHOOSE</div>', unsafe_allow_html=True)
        st.code(choose_block, language="text")

        st.markdown('<div class="section-head">Octet / Valence</div>', unsafe_allow_html=True)
        st.code(octet_report, language="text")

        st.markdown('<div class="section-head">SMILES</div>', unsafe_allow_html=True)
        try:
            atoms_packed = [
                (i + 1, atoms[i], coords[i][0], coords[i][1], coords[i][2])
                for i in range(len(atoms))
            ]
            smiles_txt = engine.ilp_to_smiles(
                atoms_packed,
                bonds_1b,
                fc_out_1b,
                dative_ml_pairs=dative_ml_pairs,
                edges=raw_edges,
            )
            st.code(smiles_txt, language="text")
        except ImportError as exc:
            st.caption(f"SMILES skipped: {exc}")
        except Exception as exc:
            st.warning(f"SMILES generation failed: {exc}")

        st.markdown('<div class="section-head">Coordination</div>', unsafe_allow_html=True)
        st.code(cbc_report, language="text")

        if ox_summary:
            st.markdown(
                '<div class="section-head">TM Oxidation State</div>',
                unsafe_allow_html=True,
            )
            st.code(ox_summary, language="text")

        export_text = (
            f"Read {len(atoms)} atoms  (charge={int(mol_charge)})\n\n"
            f"{choose_block}\n\n"
            f"{octet_report}\n"
        )
        try:
            export_text += f"\nSMILES:\n{smiles_txt}\n"
        except NameError:
            pass
        export_text += f"\n{cbc_report}\n"
        if ox_report:
            export_text += f"\n{ox_report}\n"
        export_text += f"\nStats:\n{stats}\n"
        st.download_button(
            "Download Result (.txt)",
            data=export_text,
            file_name=f"{_xyz_export_stem(xyz_file, demo_xyz_name=st.session_state.get('demo_xyz_name'))}_ilp_result.txt",
            mime="text/plain",
        )
        st.markdown("</div>", unsafe_allow_html=True)

    except UnsupportedTransitionMetalError as exc:
        show_unsupported_tm_error(exc)
    except IlpSolveError as exc:
        show_ilp_solve_error(exc, atoms=atoms, coords=coords, raw_edges=raw_edges)
    except ValueError as exc:
        msg = str(exc)
        if msg.startswith("Unsupported element symbol"):
            st.error("Unsupported element symbol in XYZ")
            st.markdown(msg)
        else:
            st.error(f"Invalid XYZ input: {exc}")
    except Exception as exc:
        if "ILP failed" in str(exc):
            show_ilp_solve_error(exc, atoms=atoms, coords=coords, raw_edges=raw_edges)
        else:
            st.error(f"Analysis failed: {exc}")


def main() -> None:
    st.set_page_config(page_title=APP_TITLE, layout="wide", initial_sidebar_state="collapsed")
    inject_styles()
    run_analyzer_app()


if __name__ == "__main__":
    main()

