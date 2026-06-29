# TMC-Lewis-Structure

Code repository for the manuscript **Automated Generation of Lewis-like Structures for Transition-Metal Complexes**.

ILP-based Lewis structures for transition-metal complexes: bond orders, formal charges, CBC (L/X/Z), oxidation states, 3D viewer, and SMILES from XYZ.

## Repository layout

| Path | Role |
|------|------|
| `streamlit_ilp_app.py` | Web UI (Streamlit) |
| `Lewis-engine-ILP.py` | ILP Lewis / CBC / SMILES backend |
| `XYZ_Viewer.html` | Embedded 3D structure viewer |
| `test/` | Built-in example XYZ files |

## Run locally

```bash
pip install -r requirements.txt
streamlit run streamlit_ilp_app.py
```

Optional flags (e.g. bind to all interfaces):

```bash
streamlit run streamlit_ilp_app.py --server.address 0.0.0.0 --server.port 8501
```

## Use `Lewis-engine-ILP.py` directly

Command-line usage:

```bash
python Lewis-engine-ILP.py <molecule.xyz> [molecular_charge]
```

Examples:

```bash
# Neutral example
python Lewis-engine-ILP.py test/ABELOK_charge_0.xyz

# Charged example (+1)
python Lewis-engine-ILP.py test/ABASEC_charge_plus_1.xyz 1
```

The script prints:
- Lewis assignment summary (bond orders, lone pairs, formal charges)
- Octet check and CBC ligand classification
- TM oxidation-state / sigma-bond report
- ILP-derived SMILES (RDKit parse check)
