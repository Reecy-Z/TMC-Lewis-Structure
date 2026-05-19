#!/usr/bin/env python3
"""
Streamlit app for ILP-based Lewis structure analysis (Lewis-engine-ILP backend).

Run:
    streamlit run streamlit_ilp_app.py --server.address 0.0.0.0 --server.port 8501

Auth credentials: copy `.streamlit/secrets.toml.example` to `.streamlit/secrets.toml`.
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
VIEWER_FILE = ROOT_DIR / "XYZ_Viewer_fixed_V2.html"
# Streamlit iframe height for the embedded XYZ viewer (px). Increase for a taller 3D panel.
VIEWER_IFRAME_HEIGHT = 1100
VIEWER_MIN_HEIGHT_PX = 820

# Built-in XYZ examples (paths relative to repo root).
DEMO_CASES = (
    {
        "id": "abasec",
        "label": "ABASEC (+1)",
        "path": ROOT_DIR / "test_5_22" / "ABASEC_charge_plus_1.xyz",
        "charge": 1,
    },
    {
        "id": "abelok",
        "label": "ABELOK (0)",
        "path": ROOT_DIR / "test_5_22" / "ABELOK_charge_0.xyz",
        "charge": 0,
    },
    {
        "id": "kadxao",
        "label": "KADXAO (0)",
        "path": ROOT_DIR / "test_5_22" / "KADXAO_charge_0.xyz",
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
[data-testid="stAppViewContainer"]{
  margin-left:0 !important;
}
[data-testid="stToolbar"],
[data-testid="stDecoration"],
[data-testid="stStatusWidget"]{
  display:none !important;
}
section.main{
  padding-top:0 !important;
  padding-left:1rem !important;
  padding-right:1rem !important;
}
section.main .block-container,
section.main [data-testid="stMainBlockContainer"]{
  padding-top:0 !important;
  padding-bottom:1rem !important;
  max-width:100% !important;
}
section.main > div:first-child{
  padding-top:0 !important;
  margin-top:0 !important;
}
section.main [data-testid="stVerticalBlock"] > div:first-child{
  padding-top:0 !important;
  margin-top:0 !important;
}
section.main [data-testid="stElementContainer"]{
  padding-top:0 !important;
  margin-top:0 !important;
}
section.main [data-testid="stVerticalBlockBorderWrapper"]{
  margin-top:0 !important;
}
[data-testid="stHeader"] { display:none; }
#MainMenu, footer { visibility:hidden; }
section[data-testid="stSidebar"], [data-testid="stSidebarCollapsedControl"] {
  display: none !important;
}

/* Top bar container */
[data-testid="stVerticalBlockBorderWrapper"]:has(.orbis-topbar-inner){
  background: linear-gradient(135deg,#111827,#1f2937 55%,#374151) !important;
  border:none !important;
  border-radius:12px !important;
  padding:8px 14px !important;
  margin:0 0 8px 0 !important;
}
[data-testid="stVerticalBlockBorderWrapper"]:has(.orbis-topbar-inner) [data-testid="stHorizontalBlock"]{
  align-items:center !important;
}
[data-testid="stVerticalBlockBorderWrapper"]:has(.orbis-topbar-inner) .orbis-title{ color:#fff; }
[data-testid="stVerticalBlockBorderWrapper"]:has(.orbis-topbar-inner) .orbis-sub{ color:#fff; opacity:.68; }
[data-testid="stVerticalBlockBorderWrapper"]:has(.orbis-topbar-inner) [data-testid="column"]:last-child{
  display:flex !important;
  justify-content:flex-end !important;
  align-items:center !important;
}
[data-testid="stVerticalBlockBorderWrapper"]:has(.orbis-topbar-inner) [data-testid="column"]:last-child .stButton{
  margin:0 !important;
}
[data-testid="stVerticalBlockBorderWrapper"]:has(.orbis-topbar-inner) [data-testid="column"]:last-child .stButton > button{
  min-height:34px !important;
  padding:6px 16px !important;
  font-size:12px !important;
  font-weight:600 !important;
  background:#f1f5f9 !important;
  color:#1e293b !important;
  border:1px solid #cbd5e1 !important;
  border-radius:8px !important;
  box-shadow:0 1px 2px rgba(0,0,0,.06) !important;
}
[data-testid="stVerticalBlockBorderWrapper"]:has(.orbis-topbar-inner) [data-testid="column"]:last-child .stButton > button:hover{
  background:#fff !important;
  border-color:#94a3b8 !important;
  color:#0f172a !important;
}

/* Control card container */
[data-testid="stVerticalBlockBorderWrapper"]:has(.control-card-inner){
  background:#ffffff !important;
  border:1px solid #d0dae8 !important;
  border-radius:12px !important;
  padding:8px 10px 6px 10px !important;
  margin:0 0 8px 0 !important;
}
[data-testid="stVerticalBlockBorderWrapper"]:has(.control-card-inner) [data-testid="stHorizontalBlock"]{
  align-items:stretch !important;
}
[data-testid="stVerticalBlockBorderWrapper"]:has(.control-card-inner) [data-testid="stHorizontalBlock"] > [data-testid="column"]{
  display:flex !important;
  flex-direction:column !important;
}
[data-testid="stVerticalBlockBorderWrapper"]:has(.control-card-inner) [data-testid="stFileUploader"] section{
  flex:1 1 auto;
  min-height:88px;
  display:flex;
  flex-direction:column;
  justify-content:center;
}
[data-testid="stVerticalBlockBorderWrapper"]:has(.control-card-inner) [data-testid="column"]:nth-child(2) [data-testid="stNumberInput"]{
  margin-top:auto;
  margin-bottom:auto;
}
[data-testid="stVerticalBlockBorderWrapper"]:has(.control-card-inner) [data-testid="column"]:nth-child(3) .stButton{
  flex:1 1 auto;
  display:flex !important;
  flex-direction:column !important;
  margin-top:0 !important;
}
[data-testid="stVerticalBlockBorderWrapper"]:has(.control-card-inner) [data-testid="column"]:nth-child(3) .stButton > button{
  flex:1 1 auto;
  min-height:88px !important;
  height:100% !important;
}

.orbis-topbar-inner,
.control-card-inner{ display:none !important; height:0 !important; margin:0 !important; padding:0 !important; }

.orbis-left{ display:flex; align-items:center; gap:10px; }
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
.status-line{
  margin-top:8px;
  border:1px solid #d0dae8;
  background:#f7f9fc;
  border-radius:10px;
  padding:6px 10px;
  font-family:'JetBrains Mono', monospace;
  font-size:10px;
  color:#64748b;
}
.login-wrap{
  max-width:400px;
  margin:12vh auto 0 auto;
  padding:0 12px;
}
.login-card{
  background:#ffffff;
  border:1px solid #d0dae8;
  border-radius:14px;
  padding:28px 26px 22px 26px;
  box-shadow:0 8px 28px rgba(15,23,42,.08);
}
.login-brand{
  text-align:center;
  margin-bottom:22px;
}
.login-brand .orbis-title{ font-size:22px; }
.login-brand .orbis-sub{ font-size:12px; margin-top:4px; }
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


def ensure_unknown_elements(backend, atoms):
    unknown = sorted(set(a for a in atoms if a not in backend.VALENCE))
    for u in unknown:
        backend.VALENCE[u] = 4
        backend.CORE_E[u] = 0
        backend.COV_R[u] = 0.77
        backend.STD_CAP[u] = 4
        backend.ENEG[u] = 2.55


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


def show_ilp_solve_error(exc: IlpSolveError | RuntimeError) -> None:
    st.error("ILP could not find a valid Lewis structure")
    st.warning(
        "Try adjusting the **total molecular charge** (Mol Charge) and/or the **XYZ geometry** "
        "(bond lengths, connectivity, metal–ligand contacts), then run the analysis again."
    )
    st.caption(f"Solver detail: {exc}")


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
        bonds, lp_out, fc_out = engine.solve_bond_orders_aromatic(
            atoms_packed,
            raw_edges,
            aromatic_systems,
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


_VIEWER_SUPPLEMENTAL_DONORS = frozenset({"O", "N", "S"})
_SUPPLEMENTAL_ANGLE_MIN_DEG = 170.0
_SUPPLEMENTAL_ANGLE_MAX_DEG = 190.0


def _bo_ij(bo, i, j):
    return bo.get((min(i, j), max(i, j)), 0)


def _angle_at_center_deg(coords, a_idx, center_idx, c_idx):
    """Angle at center_idx (degrees) between vectors center→a and center→c."""
    ax, ay, az = coords[a_idx]
    bx, by, bz = coords[center_idx]
    cx, cy, cz = coords[c_idx]
    v1 = (ax - bx, ay - by, az - bz)
    v2 = (cx - bx, cy - by, cz - bz)
    n1 = math.sqrt(v1[0] ** 2 + v1[1] ** 2 + v1[2] ** 2)
    n2 = math.sqrt(v2[0] ** 2 + v2[1] ** 2 + v2[2] ** 2)
    if n1 < 1e-10 or n2 < 1e-10:
        return None
    cosang = max(-1.0, min(1.0, (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (n1 * n2)))
    return math.degrees(math.acos(cosang))


def should_draw_supplemental_dative_arrow(
    atoms, coords, bo, lig, tm, ml_order, backend, *, lp=None
):
    """
    O/N/S covalent M–L with lp >= 1: draw extra arrow if
      (1) M–L bond order is 2, or
      (2) M–L is single AND some non-metal neighbor of lig is double-bonded to lig
          AND ∠(TM–lig–neighbor) is in [170°, 190°].
    """
    if (lp or {}).get(lig, 0) < 1:
        return False
    if ml_order == 2:
        return True
    if ml_order != 1:
        return False
    n = len(atoms)
    for nbr in range(n):
        if nbr in (lig, tm):
            continue
        if backend.is_TM(atoms[nbr]):
            continue
        if _bo_ij(bo, lig, nbr) < 2:
            continue
        ang = _angle_at_center_deg(coords, tm, lig, nbr)
        if ang is not None and _SUPPLEMENTAL_ANGLE_MIN_DEG <= ang <= _SUPPLEMENTAL_ANGLE_MAX_DEG:
            return True
    return False


def supplemental_dative_on_covalent_ml(
    atoms, coords, bo, lp, backend, *, pure_l_pairs=None
):
    """Build lp_on_covalent_ml dative entries for the 3D viewer (see should_draw_*)."""
    pure_l_pairs = pure_l_pairs or set()
    out = []
    seen = set()
    for (i, j), order in bo.items():
        if order <= 0:
            continue
        si, sj = atoms[i], atoms[j]
        if backend.is_TM(si) and not backend.is_TM(sj):
            tm, lig = i, j
        elif backend.is_TM(sj) and not backend.is_TM(si):
            tm, lig = j, i
        else:
            continue
        if atoms[lig] not in _VIEWER_SUPPLEMENTAL_DONORS:
            continue
        key = (min(tm, lig), max(tm, lig))
        if key in pure_l_pairs or key in seen:
            continue
        if not should_draw_supplemental_dative_arrow(
            atoms, coords, bo, lig, tm, int(order), backend, lp=lp
        ):
            continue
        out.append(
            {
                "donor_atoms": [lig],
                "acceptor": tm,
                "pi_bond_pair": None,
                "kind": "lp_on_covalent_ml",
                "covalent_order": int(order),
            }
        )
        seen.add(key)
    return out


def _ml_bundle_bond_keys(dative_bonds):
    keys = set()
    for d in dative_bonds:
        if d.get("kind") != "lp_on_covalent_ml":
            continue
        donors = d.get("donor_atoms") or []
        if len(donors) != 1:
            continue
        lig, tm = donors[0], d["acceptor"]
        keys.add((min(lig, tm), max(lig, tm)))
    return keys


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


def tm_oxidation_warning_messages(
    engine,
    atoms,
    fc,
    mol_charge: int,
    bo=None,
    *,
    cbc_interaction_records=None,
    cbc_neighbor_types=None,
) -> list[str]:
    common = engine.TM_COMMON_OXIDATION_STATES
    results = engine.tm_oxidation_from_charge_balance(
        atoms,
        fc,
        mol_charge,
        bo=bo,
        cbc_interaction_records=cbc_interaction_records,
        cbc_neighbor_types=cbc_neighbor_types,
    )
    msgs: list[str] = []
    for m in sorted(results):
        sym = atoms[m]
        ox_state, _others_sum, _sigma_sum = results[m]
        ref = common.get(sym)
        if ref is not None and ox_state not in ref:
            msgs.append(
                f"{sym}{m + 1}: oxidation state {ox_state} not in common oxidation states {ref}"
            )
    return msgs


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

    dative_bonds.extend(
        supplemental_dative_on_covalent_ml(
            atoms, coords, bo, lp, backend, pure_l_pairs=l_pairs
        )
    )
    bundle_keys = _ml_bundle_bond_keys(dative_bonds)

    bonds_list = []
    bo_keys = set()
    for (i, j), order in sorted(bo.items()):
        key = (min(i, j), max(i, j))
        bo_keys.add(key)
        if key in l_pairs or key in bundle_keys:
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
        key = (min(i, j), max(i, j))
        if key in bundle_keys:
            continue
        bonds_list.append({"i": i, "j": j, "order": min(count, 3)})

    return {
        "atoms": atoms,
        "coords": coords,
        "bonds": bonds_list,
        "dative_bonds": dative_bonds,
        "formal_charges": list(fc),
    }


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
    if (visInfo) visInfo.textContent = 'ILP backend data loaded from Streamlit.';
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
      modeChip.textContent = `XYZ · ${{state.atoms.length}} atoms`;
      modeChip.className = 'ok';
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
    st.caption("Viewer uses original HTML renderer with ILP-derived bond/dative payload.")


def init_auth_session() -> None:
    if "authenticated" not in st.session_state:
        st.session_state.authenticated = False


def get_auth_credentials() -> tuple[str, str]:
    """Read login credentials from .streamlit/secrets.toml [auth] section."""
    try:
        auth = st.secrets["auth"]
        username = auth["username"]
        password = auth["password"]
    except (KeyError, TypeError, AttributeError) as exc:
        raise RuntimeError(
            "Missing auth secrets. Create `.streamlit/secrets.toml` with:\n\n"
            "[auth]\n"
            'username = "your_username"\n'
            'password = "your_password"'
        ) from exc
    return str(username), str(password)


def show_login_page(auth_username: str, auth_password: str) -> None:
    st.markdown(
        """
        <div class="login-wrap">
          <div class="login-card">
            <div class="login-brand">
              <div class="orbis-title">Orbis QC</div>
              <div class="orbis-sub">TMC Lewis Analyzer — sign in to continue</div>
            </div>
          </div>
        </div>
        """,
        unsafe_allow_html=True,
    )
    _left, center, _right = st.columns([1, 1.2, 1])
    with center:
        with st.form("orbis_login_form", clear_on_submit=False):
            username = st.text_input("Username", placeholder="Enter username")
            password = st.text_input("Password", type="password", placeholder="Enter password")
            submitted = st.form_submit_button("Sign in", use_container_width=True, type="primary")
        if submitted:
            if username == auth_username and password == auth_password:
                st.session_state.authenticated = True
                st.rerun()
            else:
                st.error("Invalid username or password.")


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

    with st.container(border=True):
        st.markdown('<span class="orbis-topbar-inner"></span>', unsafe_allow_html=True)
        top_left, top_right = st.columns([5.4, 0.8], gap="small")
        with top_left:
            st.markdown(
                """
                <div class="orbis-left">
                  <div class="orbis-icon">✦</div>
                  <div>
                    <div class="orbis-title">Orbis QC</div>
                    <div class="orbis-sub">XYZ Structure Viewer + Lewis Analysis</div>
                  </div>
                </div>
                """,
                unsafe_allow_html=True,
            )
        with top_right:
            if st.button("Log out", key="logout_btn"):
                st.session_state.authenticated = False
                st.rerun()

    engine = load_engine()

    with st.container(border=True):
        st.markdown('<span class="control-card-inner"></span>', unsafe_allow_html=True)
        i1, i2, i3 = st.columns([2.4, 0.9, 1.2], gap="small")
        with i1:
            xyz_file = st.file_uploader("XYZ", type=["xyz"], label_visibility="collapsed")
            if xyz_file is not None:
                st.session_state.pop("demo_xyz_text", None)
                st.session_state.pop("demo_xyz_name", None)
        with i2:
            st.markdown('<p class="field-label">Charge</p>', unsafe_allow_html=True)
            mol_charge = st.number_input(
                "Charge",
                step=1,
                format="%d",
                label_visibility="collapsed",
                key="mol_charge",
                help="Total molecular charge (e.g. −1 for anions, +1 for cations).",
            )
        with i3:
            st.markdown('<p class="field-label">&nbsp;</p>', unsafe_allow_html=True)
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
        ensure_unknown_elements(engine, atoms)

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
            cbc_report, cbc_bundle = format_cbc_report(
                engine, atoms, coords, bo, lp, fc, int(mol_charge),
                metal_adjacency_edges=metal_adj_0,
            )
            cbc_interaction_records, _cbc_neighbor_types = cbc_bundle
            ox_report = format_tm_oxidation_report(
                engine, atoms, fc, int(mol_charge), bo=bo,
                cbc_interaction_records=cbc_interaction_records,
            )
            ox_warnings = tm_oxidation_warning_messages(
                engine, atoms, fc, int(mol_charge), bo=bo,
                cbc_interaction_records=cbc_interaction_records,
            )
            viewer_payload = build_viewer_payload(
                engine, atoms, coords, bo, lp, fc, int(mol_charge),
                metal_adjacency_edges=metal_adj_0,
            )

        st.markdown('<div class="result-card">', unsafe_allow_html=True)
        show_3d_preview(viewer_payload)
        st.markdown("</div>", unsafe_allow_html=True)

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

        st.markdown('<div class="section-head">CBC Report</div>', unsafe_allow_html=True)
        st.code(cbc_report, language="text")

        if ox_report:
            st.markdown(
                '<div class="section-head">TM Oxidation State (charge balance)</div>',
                unsafe_allow_html=True,
            )
            for msg in ox_warnings:
                st.warning(msg)
            st.code(ox_report, language="text")

        export_text = (
            f"Read {len(atoms)} atoms  (charge={int(mol_charge)})\n\n"
            f"{choose_block}\n"
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

        st.markdown(
            '<div class="status-line">Drag → rotate · Scroll → zoom · Right-drag → pan · ILP backend active</div>',
            unsafe_allow_html=True,
        )

    except UnsupportedTransitionMetalError as exc:
        show_unsupported_tm_error(exc)
    except IlpSolveError as exc:
        show_ilp_solve_error(exc)
    except ValueError as exc:
        st.error(f"Invalid XYZ input: {exc}")
    except Exception as exc:
        if "ILP failed" in str(exc):
            show_ilp_solve_error(exc)
        else:
            st.error(f"Analysis failed: {exc}")


def main() -> None:
    st.set_page_config(
        page_title=APP_TITLE,
        layout="wide",
        initial_sidebar_state="collapsed",
    )
    inject_styles()
    init_auth_session()
    try:
        auth_username, auth_password = get_auth_credentials()
    except RuntimeError as exc:
        st.error(str(exc))
        st.stop()
    if not st.session_state.authenticated:
        show_login_page(auth_username, auth_password)
        return
    run_analyzer_app()


if __name__ == "__main__":
    main()

