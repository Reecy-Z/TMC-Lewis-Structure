#!/usr/bin/env python3
"""Re-solve prior ligand mismatches with updated ILP + reference connectivity."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from compare_smiles_connectivity_graphs import connectivity_key
from evaluate_ilp_vs_xyz2mol_smiles import compare_pair, extract_features
from run_reference_graph_ilp import (
    decode_edges,
    init_worker,
    load_csv_by_id,
    run_ilp_with_edges,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_MISMATCHES = (
    ROOT / "reference_graph_ilp_experiment" / "ligand_mismatch_classification.csv"
)
DEFAULT_REFERENCE_EDGES = (
    ROOT / "reference_graph_ilp_experiment" / "reference_edges.csv"
)
DEFAULT_REFERENCE_SMILES = ROOT / "tmqmg_smiles.csv"
DEFAULT_OUT_DIR = (
    ROOT / "reference_graph_ilp_experiment" / "hypervalent_degree_rerun"
)
DEFAULT_ENGINE = ROOT / "Lewis-engine-ILP.py"
DEFAULT_XYZ_DIR = Path("/data/jingyuan_data/tmqmg")


def solve_one(job: tuple[str, int, str]) -> dict[str, str | int | float]:
    csd_id, charge, encoded = job
    return run_ilp_with_edges(csd_id, charge, decode_edges(encoded))


def write_rows(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    if fields is None:
        fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mismatches", type=Path, default=DEFAULT_MISMATCHES)
    parser.add_argument("--reference-edges", type=Path, default=DEFAULT_REFERENCE_EDGES)
    parser.add_argument("--reference-smiles", type=Path, default=DEFAULT_REFERENCE_SMILES)
    parser.add_argument("--engine", type=Path, default=DEFAULT_ENGINE)
    parser.add_argument("--xyz-dir", type=Path, default=DEFAULT_XYZ_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    mismatches = list(load_csv_by_id(args.mismatches).values())
    mismatches.sort(key=lambda row: row["IDs"])
    if args.limit:
        mismatches = mismatches[: args.limit]
    reference_edges = load_csv_by_id(args.reference_edges)
    reference_smiles = load_csv_by_id(args.reference_smiles)

    jobs: list[tuple[str, int, str]] = []
    skipped: list[dict[str, str]] = []
    for row in mismatches:
        csd_id = row["IDs"]
        edge_row = reference_edges.get(csd_id)
        if (
            edge_row is None
            or edge_row.get("validated") != "1"
            or not edge_row.get("edges_0")
        ):
            skipped.append(
                {
                    "IDs": csd_id,
                    "reason": "reference_graph_unavailable",
                }
            )
            continue
        jobs.append((csd_id, int(edge_row["charge"]), edge_row["edges_0"]))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    print(f"solving {len(jobs)} mismatches; skipped {len(skipped)}", flush=True)

    solved: dict[str, dict[str, str | int | float]] = {}
    with ProcessPoolExecutor(
        max_workers=args.workers,
        initializer=init_worker,
        initargs=(str(args.engine), str(args.xyz_dir)),
    ) as pool:
        for index, row in enumerate(pool.map(solve_one, jobs, chunksize=1), 1):
            solved[str(row["IDs"])] = row
            if index % 50 == 0 or index == len(jobs):
                print(f"solved {index} / {len(jobs)}", flush=True)

    result_fields = [
        "IDs",
        "charge",
        "status",
        "ilp_seconds",
        "n_atoms",
        "rdkit_parse_ok",
        "error",
        "smiles_ilp",
    ]
    write_rows(
        args.out_dir / "per_structure_results.csv",
        [solved[csd_id] for csd_id, _charge, _edges in jobs],
        result_fields,
    )

    details: list[dict[str, str | int]] = []
    motif_fixed: Counter[str] = Counter()
    motif_remain: Counter[str] = Counter()
    mechanism_fixed: Counter[str] = Counter()
    mechanism_remain: Counter[str] = Counter()
    counts = Counter()

    for old in mismatches:
        csd_id = old["IDs"]
        ref_smi = reference_smiles[csd_id]["smiles_CSD_fix"]
        result = solved.get(csd_id)
        motif = old.get("motif_class") or old.get("old_motif_class", "")
        mechanism = old.get("mechanism_class") or old.get("old_mechanism_class", "")
        if result is None:
            status = "skipped"
            smiles = ""
        else:
            status = str(result["status"])
            smiles = str(result.get("smiles_ilp") or "")
        ilp_feat = extract_features(smiles, convert_x=True) if smiles else None
        ref_feat = extract_features(ref_smi, convert_x=False)
        parsed = ilp_feat is not None and ref_feat is not None
        comparisons = compare_pair(ilp_feat, ref_feat) if parsed else {}
        ligand_equiv = int(comparisons.get("ligand_equiv", False))
        metal_oxidation = int(comparisons.get("metal_oxidation", False))
        metal_connectivity = int(comparisons.get("metal_connectivity", False))
        full_connectivity = int(
            parsed
            and connectivity_key(smiles) is not None
            and connectivity_key(smiles) == connectivity_key(ref_smi)
        )
        if ligand_equiv:
            counts["ligand_fixed"] += 1
            motif_fixed[motif] += 1
            mechanism_fixed[mechanism] += 1
        else:
            counts["ligand_remain"] += 1
            motif_remain[motif] += 1
            mechanism_remain[mechanism] += 1
        if metal_oxidation:
            counts["oxidation_fixed"] += 1
        else:
            counts["oxidation_remain"] += 1
        details.append(
            {
                "IDs": csd_id,
                "status": status,
                "both_parsed": int(parsed),
                "ligand_equiv": ligand_equiv,
                "ligand_canonical": int(comparisons.get("ligand_canonical", False)),
                "metal_connectivity": metal_connectivity,
                "metal_oxidation": metal_oxidation,
                "full_connectivity": full_connectivity,
                "ox_ilp": ilp_feat["ox"] if ilp_feat else "",
                "ox_reference": ref_feat["ox"] if ref_feat else old.get(
                    "reference_oxidation", ""
                ),
                "old_ilp_oxidation": old.get("ilp_oxidation", ""),
                "old_mechanism_class": mechanism,
                "old_motif_class": motif,
                "ilp_ligands": (
                    " || ".join(ilp_feat["ligands"]) if ilp_feat else ""
                ),
                "reference_ligands": (
                    " || ".join(ref_feat["ligands"]) if ref_feat else ""
                ),
                "old_ilp_ligands": old.get("ilp_ligands", ""),
                "smiles_ilp": smiles,
                "error": (result or {}).get("error", ""),
            }
        )

    write_rows(args.out_dir / "agreement.csv", details)
    write_rows(
        args.out_dir / "fixed.csv",
        [row for row in details if int(row["ligand_equiv"]) == 1],
    )
    write_rows(
        args.out_dir / "still_mismatch.csv",
        [row for row in details if int(row["ligand_equiv"]) == 0],
    )

    total = len(details)
    motif_rows = []
    for motif in sorted(set(motif_fixed) | set(motif_remain)):
        n_old = motif_fixed[motif] + motif_remain[motif]
        motif_rows.append(
            {
                "motif_class": motif,
                "old_mismatches": n_old,
                "fixed": motif_fixed[motif],
                "remain": motif_remain[motif],
                "fixed_percent": 100.0 * motif_fixed[motif] / n_old if n_old else 0.0,
            }
        )
    motif_rows.sort(key=lambda row: (-int(row["fixed"]), str(row["motif_class"])))
    write_rows(args.out_dir / "fixed_by_motif.csv", motif_rows)

    mechanism_rows = []
    for mechanism in sorted(set(mechanism_fixed) | set(mechanism_remain)):
        n_old = mechanism_fixed[mechanism] + mechanism_remain[mechanism]
        mechanism_rows.append(
            {
                "mechanism_class": mechanism,
                "old_mismatches": n_old,
                "fixed": mechanism_fixed[mechanism],
                "remain": mechanism_remain[mechanism],
                "fixed_percent": (
                    100.0 * mechanism_fixed[mechanism] / n_old if n_old else 0.0
                ),
            }
        )
    write_rows(args.out_dir / "fixed_by_mechanism.csv", mechanism_rows)

    status_counts = Counter(str(solved[csd_id]["status"]) for csd_id, _, _ in jobs)
    summary = {
        "old_ligand_mismatches": total,
        "jobs_solved": len(jobs),
        "skipped_unavailable_graph": len(skipped),
        "ilp_status_counts": dict(status_counts),
        "ligand_equiv_fixed": counts["ligand_fixed"],
        "ligand_equiv_remain": counts["ligand_remain"],
        "ligand_equiv_fixed_percent": (
            100.0 * counts["ligand_fixed"] / total if total else 0.0
        ),
        "oxidation_now_match": counts["oxidation_fixed"],
        "oxidation_still_mismatch": counts["oxidation_remain"],
        "fixed_by_motif": {
            row["motif_class"]: {
                "fixed": row["fixed"],
                "remain": row["remain"],
            }
            for row in motif_rows
        },
        "fixed_by_mechanism": {
            row["mechanism_class"]: {
                "fixed": row["fixed"],
                "remain": row["remain"],
            }
            for row in mechanism_rows
        },
    }
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    print(f"Wrote {args.out_dir}")


if __name__ == "__main__":
    main()
