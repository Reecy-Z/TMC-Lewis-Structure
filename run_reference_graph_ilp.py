#!/usr/bin/env python3
"""Run ILP with baseline/reference ligand and metal connectivity combinations."""

from __future__ import annotations

import argparse
import csv
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

from batch_tmqmg_ilp_benchmark import (
    DETAIL_CSV_COLUMNS,
    RowResult,
    load_lewis_engine,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_ENGINE = ROOT / "Lewis-engine-ILP.py"
DEFAULT_XYZ_DIR = Path("/data/jingyuan_data/tmqmg")
DEFAULT_REFERENCE_EDGES = ROOT / "reference_graph_ilp_experiment" / "reference_edges.csv"
DEFAULT_BASELINE = ROOT / "ilp_mismatch_rerun_current" / "merged_per_structure_results.csv"
DEFAULT_OUT_DIR = ROOT / "reference_graph_ilp_experiment"

TREATMENTS = ("ml_only", "ligand_only", "full")
OUTPUT_FIELDS = DETAIL_CSV_COLUMNS + [
    "treatment",
    "edge_changed",
    "result_source",
]

_ENGINE: Any | None = None
_RDKIT_CHEM: Any | None = None
_XYZ_DIR: Path | None = None


def decode_edges(value: str) -> list[tuple[int, int]]:
    if not value:
        return []
    return [
        tuple(int(part) for part in token.split("-", 1))
        for token in value.split(";")
        if token
    ]


def encode_edges(edges: set[tuple[int, int]] | list[tuple[int, int]]) -> str:
    return ";".join(f"{i}-{j}" for i, j in sorted(edges))


def row_result_dict(result: RowResult) -> dict[str, str | int | float]:
    return {
        "IDs": result.csd_id,
        "charge": result.charge,
        "status": result.status,
        "ilp_seconds": f"{result.ilp_seconds:.6f}",
        "n_atoms": result.n_atoms,
        "rdkit_parse_ok": int(result.rdkit_parse_ok),
        "error": result.error,
        "smiles_ilp": result.smiles,
    }


def init_worker(engine_path: str, xyz_dir: str) -> None:
    global _ENGINE, _RDKIT_CHEM, _XYZ_DIR
    from rdkit import Chem

    _ENGINE = load_lewis_engine(Path(engine_path))
    _RDKIT_CHEM = Chem
    _XYZ_DIR = Path(xyz_dir)


def pairs_to_engine_edges(
    atoms: list,
    pairs: list[tuple[int, int]],
) -> list[tuple[int, int, str, str]]:
    edges = []
    for begin, end in pairs:
        atom_i = atoms[begin]
        atom_j = atoms[end]
        edges.append((atom_i[0], atom_j[0], atom_i[1], atom_j[1]))
    return edges


def run_ilp_with_edges(
    csd_id: str,
    charge: int,
    edge_pairs: list[tuple[int, int]],
) -> dict[str, str | int | float]:
    assert _ENGINE is not None and _XYZ_DIR is not None
    xyz_path = _XYZ_DIR / f"{csd_id}.xyz"
    started = time.perf_counter()
    try:
        atoms = _ENGINE.read_xyz(str(xyz_path))
        edges = pairs_to_engine_edges(atoms, edge_pairs)
        aromatic_systems = _ENGINE.aromatic_candidate_systems(atoms, edges)
        bonds, lp_out, fc_out = _ENGINE.solve_bond_orders(
            atoms,
            edges,
            aromatic_systems,
            mol_charge=charge,
            metal_adjacency_edges=edges,
        )
        bonds, lp_out, fc_out, _labels = _ENGINE.apply_heterocyclic_carbene_corrections(
            atoms,
            bonds,
            lp_out,
            fc_out,
            aromatic_systems,
            edges,
            mol_charge=charge,
        )
        bonds, lp_out, fc_out = _ENGINE.apply_eta_covalent_pi_corrections(
            atoms,
            bonds,
            lp_out,
            fc_out,
            mol_charge=charge,
            metal_adjacency_edges=edges,
        )
        coords = [[atom[2], atom[3], atom[4]] for atom in atoms]
        atom_symbols = [atom[1] for atom in atoms]
        formal_charges = [fc_out.get(atom[0], 0) for atom in atoms]
        index_to_position = {atom[0]: idx for idx, atom in enumerate(atoms)}
        metal_edges_0 = [
            (
                index_to_position[begin],
                index_to_position[end],
                symbol_i,
                symbol_j,
            )
            for begin, end, symbol_i, symbol_j in edges
            if _ENGINE.is_TM(symbol_i) ^ _ENGINE.is_TM(symbol_j)
        ]
        dative_pairs = _ENGINE.infer_dative_ml_pairs_cbc(
            atom_symbols,
            coords,
            bonds,
            lp_out,
            formal_charges,
            metal_adjacency_edges=metal_edges_0,
        )
        smiles = _ENGINE.ilp_to_smiles(
            atoms,
            bonds,
            fc_out,
            dative_ml_pairs=dative_pairs,
            edges=edges,
        )
        parsed = _RDKIT_CHEM.MolFromSmiles(smiles) if _RDKIT_CHEM is not None else None
        parse_ok = parsed is not None
        return row_result_dict(
            RowResult(
                csd_id=csd_id,
                charge=charge,
                status="ok",
                ilp_seconds=time.perf_counter() - started,
                smiles=smiles,
                rdkit_parse_ok=parse_ok,
                n_atoms=len(atoms),
            )
        )
    except RuntimeError as exc:
        status = "ilp_failed" if "ILP failed" in str(exc) else "error"
        return row_result_dict(
            RowResult(
                csd_id=csd_id,
                charge=charge,
                status=status,
                ilp_seconds=time.perf_counter() - started,
                error=str(exc),
            )
        )
    except Exception as exc:
        return row_result_dict(
            RowResult(
                csd_id=csd_id,
                charge=charge,
                status="error",
                ilp_seconds=time.perf_counter() - started,
                error=f"{type(exc).__name__}: {exc}",
            )
        )


def split_edges(
    atoms: list,
    edges: set[tuple[int, int]],
) -> tuple[set[tuple[int, int]], set[tuple[int, int]]]:
    assert _ENGINE is not None
    ligand_edges: set[tuple[int, int]] = set()
    metal_edges: set[tuple[int, int]] = set()
    for begin, end in edges:
        is_metal_edge = _ENGINE.is_TM(atoms[begin][1]) ^ _ENGINE.is_TM(atoms[end][1])
        (metal_edges if is_metal_edge else ligand_edges).add((begin, end))
    return ligand_edges, metal_edges


def prepare_variants(
    job: tuple[str, int, str],
) -> dict[str, str | int]:
    assert _ENGINE is not None and _XYZ_DIR is not None
    csd_id, charge, reference_encoded = job
    result: dict[str, str | int] = {
        "IDs": csd_id,
        "charge": charge,
        "status": "error",
        "raw_edges_0": "",
        "ml_only_edges_0": "",
        "ligand_only_edges_0": "",
        "full_edges_0": reference_encoded,
        "ml_only_changed": 0,
        "ligand_only_changed": 0,
        "full_changed": 0,
        "raw_edge_count": 0,
        "reference_edge_count": 0,
        "raw_ml_count": 0,
        "reference_ml_count": 0,
        "error": "",
    }
    try:
        xyz_path = _XYZ_DIR / f"{csd_id}.xyz"
        atoms = _ENGINE.read_xyz(str(xyz_path))
        raw_engine_edges = _ENGINE.connectivity(atoms)
        atom_id_to_position = {atom[0]: idx for idx, atom in enumerate(atoms)}
        raw = {
            tuple(
                sorted(
                    (
                        atom_id_to_position[begin],
                        atom_id_to_position[end],
                    )
                )
            )
            for begin, end, _symbol_i, _symbol_j in raw_engine_edges
        }
        reference = set(decode_edges(reference_encoded))
        if any(i < 0 or j >= len(atoms) or i >= j for i, j in reference):
            raise ValueError("reference edge index outside XYZ atom range")

        raw_ligand, raw_metal = split_edges(atoms, raw)
        reference_ligand, reference_metal = split_edges(atoms, reference)
        variants = {
            "ml_only": raw_ligand | reference_metal,
            "ligand_only": reference_ligand | raw_metal,
            "full": reference,
        }
        result.update(
            {
                "status": "ok",
                "raw_edges_0": encode_edges(raw),
                "ml_only_edges_0": encode_edges(variants["ml_only"]),
                "ligand_only_edges_0": encode_edges(variants["ligand_only"]),
                "full_edges_0": encode_edges(variants["full"]),
                "ml_only_changed": int(variants["ml_only"] != raw),
                "ligand_only_changed": int(variants["ligand_only"] != raw),
                "full_changed": int(variants["full"] != raw),
                "raw_edge_count": len(raw),
                "reference_edge_count": len(reference),
                "raw_ml_count": len(raw_metal),
                "reference_ml_count": len(reference_metal),
                "error": "",
            }
        )
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


def solve_changed_variants(
    row: dict[str, str | int],
) -> tuple[str, dict[str, dict[str, str | int | float]]]:
    csd_id = str(row["IDs"])
    charge = int(row["charge"])
    raw_encoded = str(row["raw_edges_0"])
    unique_results: dict[str, dict[str, str | int | float]] = {}
    mode_results: dict[str, dict[str, str | int | float]] = {}
    force_solve = bool(int(row.get("force_solve", 0)))
    for treatment in TREATMENTS:
        encoded = str(row[f"{treatment}_edges_0"])
        if encoded == raw_encoded and not force_solve:
            continue
        if encoded not in unique_results:
            unique_results[encoded] = run_ilp_with_edges(
                csd_id,
                charge,
                decode_edges(encoded),
            )
        mode_results[treatment] = unique_results[encoded]
    return csd_id, mode_results


def load_csv_by_id(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["IDs"]: row for row in csv.DictReader(handle)}


def normalize_baseline_row(
    csd_id: str,
    charge: int,
    baseline: dict[str, str] | None,
) -> dict[str, str | int | float]:
    if baseline is None:
        return {
            "IDs": csd_id,
            "charge": charge,
            "status": "missing_baseline",
            "ilp_seconds": "0.000000",
            "n_atoms": 0,
            "rdkit_parse_ok": 0,
            "error": "No current ILP baseline row",
            "smiles_ilp": "",
        }
    return {field: baseline.get(field, "") for field in DETAIL_CSV_COLUMNS}


def write_treatment_results(
    out_dir: Path,
    reference_rows: list[dict[str, str]],
    variant_rows: dict[str, dict[str, str | int]],
    baseline_rows: dict[str, dict[str, str]],
    solved: dict[str, dict[str, dict[str, str | int | float]]],
) -> dict[str, dict[str, int]]:
    summaries: dict[str, dict[str, int]] = {}
    for treatment in TREATMENTS:
        counts: dict[str, int] = {}
        path = out_dir / f"per_structure_results_{treatment}.csv"
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS)
            writer.writeheader()
            for reference in reference_rows:
                csd_id = reference["IDs"]
                charge = int(reference["charge"])
                variant = variant_rows.get(csd_id)
                changed = bool(
                    variant
                    and variant.get("status") == "ok"
                    and int(variant[f"{treatment}_changed"]) == 1
                )
                force_solved = bool(
                    variant
                    and int(variant.get("force_solve", 0)) == 1
                )
                if reference.get("validated") != "1" or not variant or variant.get("status") != "ok":
                    base = {
                        "IDs": csd_id,
                        "charge": charge,
                        "status": "reference_graph_unavailable",
                        "ilp_seconds": "0.000000",
                        "n_atoms": reference.get("n_atoms", 0),
                        "rdkit_parse_ok": 0,
                        "error": reference.get("error", "") or (variant or {}).get("error", ""),
                        "smiles_ilp": "",
                    }
                    source = "unavailable"
                elif changed or force_solved:
                    base = solved.get(csd_id, {}).get(treatment)
                    if base is None:
                        base = {
                            "IDs": csd_id,
                            "charge": charge,
                            "status": "missing_treatment_result",
                            "ilp_seconds": "0.000000",
                            "n_atoms": 0,
                            "rdkit_parse_ok": 0,
                            "error": "Changed graph was not solved",
                            "smiles_ilp": "",
                        }
                    source = "treatment_rerun"
                else:
                    base = normalize_baseline_row(
                        csd_id,
                        charge,
                        baseline_rows.get(csd_id),
                    )
                    source = "baseline_reused"
                output = dict(base)
                output.update(
                    {
                        "treatment": treatment,
                        "edge_changed": int(changed),
                        "result_source": source,
                    }
                )
                writer.writerow(output)
                status = str(output["status"])
                counts[status] = counts.get(status, 0) + 1
        summaries[treatment] = counts
        print(f"Wrote {path}", flush=True)
    return summaries


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", type=Path, default=DEFAULT_ENGINE)
    parser.add_argument("--xyz-dir", type=Path, default=DEFAULT_XYZ_DIR)
    parser.add_argument("--reference-edges", type=Path, default=DEFAULT_REFERENCE_EDGES)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    reference_rows = list(load_csv_by_id(args.reference_edges).values())
    reference_rows.sort(key=lambda row: row["IDs"])
    if args.limit:
        reference_rows = reference_rows[: args.limit]
    baseline_rows = load_csv_by_id(args.baseline)
    jobs = [
        (row["IDs"], int(row["charge"]), row["edges_0"])
        for row in reference_rows
        if row.get("validated") == "1"
    ]
    args.out_dir.mkdir(parents=True, exist_ok=True)

    variants: list[dict[str, str | int]] = []
    solved: dict[str, dict[str, dict[str, str | int | float]]] = {}
    with ProcessPoolExecutor(
        max_workers=args.workers,
        initializer=init_worker,
        initargs=(str(args.engine), str(args.xyz_dir)),
    ) as pool:
        for index, row in enumerate(pool.map(prepare_variants, jobs, chunksize=8), 1):
            row["force_solve"] = int(str(row["IDs"]) not in baseline_rows)
            variants.append(row)
            if index % 500 == 0:
                print(f"compared graph inputs {index} / {len(jobs)}", flush=True)

        changed_rows = [
            row
            for row in variants
            if row["status"] == "ok"
            and (
                int(row.get("force_solve", 0)) == 1
                or any(
                    int(row[f"{treatment}_changed"]) == 1
                    for treatment in TREATMENTS
                )
            )
        ]
        print(f"structures requiring at least one rerun: {len(changed_rows)}", flush=True)
        for index, (csd_id, results) in enumerate(
            pool.map(solve_changed_variants, changed_rows, chunksize=1),
            1,
        ):
            solved[csd_id] = results
            if index % 100 == 0:
                print(f"solved changed structures {index} / {len(changed_rows)}", flush=True)

    audit_fields = [
        "IDs",
        "charge",
        "status",
        "raw_edge_count",
        "reference_edge_count",
        "raw_ml_count",
        "reference_ml_count",
        "ml_only_changed",
        "ligand_only_changed",
        "full_changed",
        "force_solve",
        "raw_edges_0",
        "ml_only_edges_0",
        "ligand_only_edges_0",
        "full_edges_0",
        "error",
    ]
    audit_path = args.out_dir / "edge_input_audit.csv"
    with audit_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=audit_fields)
        writer.writeheader()
        writer.writerows(variants)

    variant_by_id = {str(row["IDs"]): row for row in variants}
    treatment_summaries = write_treatment_results(
        args.out_dir,
        reference_rows,
        variant_by_id,
        baseline_rows,
        solved,
    )
    edge_change_counts = {
        treatment: sum(int(row[f"{treatment}_changed"]) for row in variants)
        for treatment in TREATMENTS
    }
    summary = {
        "reference_rows": len(reference_rows),
        "validated_reference_rows": len(jobs),
        "edge_change_counts": edge_change_counts,
        "structures_rerun": len(changed_rows),
        "treatment_status_counts": treatment_summaries,
    }
    (args.out_dir / "run_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    print(f"Wrote {audit_path}")


if __name__ == "__main__":
    main()
