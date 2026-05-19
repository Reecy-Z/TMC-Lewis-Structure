# TMC-Lewis-Structure

ILP-based Lewis structures for transition-metal complexes: bond orders, formal charges, CBC (L/X/Z), oxidation states, 3D viewer, and SMILES from XYZ.

## Repository layout

| Path | Role |
|------|------|
| `streamlit_ilp_app.py` | Web UI (Streamlit) |
| `Lewis-engine-ILP.py` | ILP Lewis / CBC / SMILES backend |
| `XYZ_Viewer_fixed_V2.html` | Embedded 3D structure viewer |
| `test_5_22/` | Built-in example XYZ files |

## Run locally

```bash
pip install -r requirements.txt
cp .streamlit/secrets.toml.example .streamlit/secrets.toml
# Edit secrets.toml with [auth] username / password
streamlit run streamlit_ilp_app.py
```

## Deploy on [Streamlit Community Cloud](https://share.streamlit.io)

1. Push this repo to GitHub: `https://github.com/Reecy-Z/TMC-Lewis-Structure`
2. Sign in at [share.streamlit.io](https://share.streamlit.io) → **Create app**
3. **Repository:** `Reecy-Z/TMC-Lewis-Structure`  
   **Branch:** `main`  
   **Main file path:** `streamlit_ilp_app.py`
4. Open **Advanced settings** → **Secrets** and paste (set your own credentials):

```toml
[auth]
username = "your_username"
password = "your_password"
```

5. **Deploy**. The app URL will look like `https://tmc-lewis-structure-xxxx.streamlit.app`.

### Notes

- **RDKit** is optional: SMILES export is skipped if `rdkit-pypi` fails to install; Lewis analysis and the 3D viewer still work.
- Example structures load from `test_5_22/` via the **Quick examples** buttons.
- Do not commit `.streamlit/secrets.toml` (it is in `.gitignore`).
