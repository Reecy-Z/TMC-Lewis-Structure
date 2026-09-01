#!/usr/bin/env python3
"""Evaluate fixed-cohort agreement for reference-graph ILP treatments."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from compare_smiles_connectivity_graphs import connectivity_key
from evaluate_ilp_vs_xyz2mol_smiles import (
    _extract_one,
    compare_pair,
    is_valid_smiles_cell,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_EXPERIMENT = ROOT / "reference_graph_ilp_experiment"
DEFAULT_BASELINE = ROOT / "ilp_mismatch_rerun_current" / "merged_per_structure_results.csv"
DEFAULT_REFERENCE = ROOT / "tmqmg_smiles.csv"

METHOD_FILES = {
    "baseline": None,
    "ml_only": "per_structure_results_ml_only.csv",
    "ligand_only": "per_structure_results_ligand_only.csv",
    "full": "per_structure_results_full.csv",
}
METRICS = (
    "ligand_canonical",
    "ligand_equiv",
    "metal_connectivity",
    "metal_oxidation",
    "full_connectivity",
)


def load_ids(path: Path) -> set[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["IDs"] for row in csv.DictReader(handle)}


def load_reference(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["IDs"]: row for row in csv.DictReader(handle)}


def load_smiles_results(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    smiles: dict[str, str] = {}
    statuses: dict[str, str] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            csd_id = row["IDs"]
            statuses[csd_id] = row.get("status", "")
            if (
                row.get("status") == "ok"
                and row.get("rdkit_parse_ok", "").strip().lower() in {"1", "true"}
                and is_valid_smiles_cell(row.get("smiles_ilp"))
            ):
                smiles[csd_id] = row["smiles_ilp"]
    return smiles, statuses


def extract_caches(
    methods: dict[str, dict[str, str]],
    references: dict[str, dict[str, str]],
    ids: set[str],
    workers: int,
) -> tuple[dict[tuple[str, str], dict | None], dict[str, object | None]]:
    feature_jobs: list[tuple[str, str]] = []
    seen_features: set[tuple[str, str]] = set()
    unique_smiles: list[str] = []
    seen_smiles: set[str] = set()
    for csd_id in sorted(ids):
        ref_smiles = references[csd_id]["smiles_CSD_fix"]
        job = ("ref", ref_smiles)
        if job not in seen_features:
            seen_features.add(job)
            feature_jobs.append(job)
        if ref_smiles not in seen_smiles:
            seen_smiles.add(ref_smiles)
            unique_smiles.append(ref_smiles)
        for smiles_by_id in methods.values():
            smiles = smiles_by_id.get(csd_id, "")
            if not smiles:
                continue
            job = ("ilp", smiles)
            if job not in seen_features:
                seen_features.add(job)
                feature_jobs.append(job)
            if smiles not in seen_smiles:
                seen_smiles.add(smiles)
                unique_smiles.append(smiles)

    if workers <= 1:
        feature_items = [_extract_one(job) for job in feature_jobs]
        graph_keys = [connectivity_key(smiles) for smiles in unique_smiles]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            feature_items = list(pool.map(_extract_one, feature_jobs, chunksize=16))
            graph_keys = list(pool.map(connectivity_key, unique_smiles, chunksize=32))
    return dict(feature_items), dict(zip(unique_smiles, graph_keys))


def evaluate_details(
    metric_ids: set[str],
    graph_ids: set[str],
    methods: dict[str, dict[str, str]],
    statuses: dict[str, dict[str, str]],
    references: dict[str, dict[str, str]],
    feature_cache: dict[tuple[str, str], dict | None],
    graph_cache: dict[str, object | None],
) -> list[dict[str, str | int]]:
    details: list[dict[str, str | int]] = []
    for csd_id in sorted(metric_ids):
        reference_smiles = references[csd_id]["smiles_CSD_fix"]
        reference_features = feature_cache.get(("ref", reference_smiles))
        reference_graph = graph_cache.get(reference_smiles)
        for treatment, smiles_by_id in methods.items():
            smiles = smiles_by_id.get(csd_id, "")
            features = feature_cache.get(("ilp", smiles)) if smiles else None
            parsed = features is not None and reference_features is not None
            comparisons = compare_pair(features, reference_features) if parsed else {}
            full_connectivity = bool(
                parsed
                and graph_cache.get(smiles) is not None
                and graph_cache.get(smiles) == reference_graph
            )
            details.append(
                {
                    "IDs": csd_id,
                    "in_graph_consensus": int(csd_id in graph_ids),
                    "treatment": treatment,
                    "status": statuses[treatment].get(csd_id, "missing"),
                    "both_parsed": int(parsed),
                    "ligand_canonical": int(comparisons.get("ligand_canonical", False)),
                    "ligand_equiv": int(comparisons.get("ligand_equiv", False)),
                    "metal_connectivity": int(
                        comparisons.get("metal_connectivity", False)
                    ),
                    "metal_oxidation": int(comparisons.get("metal_oxidation", False)),
                    "full_connectivity": int(full_connectivity),
                    "ox_ilp": features["ox"] if features else "",
                    "ox_reference": reference_features["ox"] if reference_features else "",
                    "neighbors_ilp": " ".join(features["neighbors"]) if features else "",
                    "neighbors_reference": (
                        " ".join(reference_features["neighbors"])
                        if reference_features
                        else ""
                    ),
                }
            )
    return details


def write_detail(path: Path, details: list[dict[str, str | int]]) -> None:
    fields = [
        "IDs",
        "in_graph_consensus",
        "treatment",
        "status",
        "both_parsed",
        *METRICS,
        "ox_ilp",
        "ox_reference",
        "neighbors_ilp",
        "neighbors_reference",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(details)


def summarize(
    details: list[dict[str, str | int]],
    metric_ids: set[str],
    graph_ids: set[str],
    validated_graph_ids: set[str],
) -> list[dict[str, str | int | float]]:
    rows: list[dict[str, str | int | float]] = []
    cohorts = {
        "metric_consensus_29110": metric_ids,
        "graph_consensus_29068": graph_ids,
        "validated_graph_consensus": validated_graph_ids,
    }
    by_key = {(str(row["IDs"]), str(row["treatment"])): row for row in details}
    for cohort_name, cohort_ids in cohorts.items():
        total = len(cohort_ids)
        for treatment in METHOD_FILES:
            comparable = sum(
                int(by_key[(csd_id, treatment)]["both_parsed"])
                for csd_id in cohort_ids
            )
            for metric in METRICS:
                matches = sum(
                    int(by_key[(csd_id, treatment)][metric])
                    for csd_id in cohort_ids
                )
                rows.append(
                    {
                        "cohort": cohort_name,
                        "treatment": treatment,
                        "metric": metric,
                        "matches": matches,
                        "total": total,
                        "comparable": comparable,
                        "fixed_denominator_percent": 100.0 * matches / total,
                        "comparable_percent": (
                            100.0 * matches / comparable if comparable else 0.0
                        ),
                    }
                )
    return rows


def transition_rows(
    details: list[dict[str, str | int]],
    cohort_ids: set[str],
) -> list[dict[str, str | int]]:
    indexed = {(str(row["IDs"]), str(row["treatment"])): row for row in details}
    output: list[dict[str, str | int]] = []
    for treatment in ("ml_only", "ligand_only", "full"):
        for metric in METRICS:
            counts = Counter()
            for csd_id in cohort_ids:
                baseline_match = bool(int(indexed[(csd_id, "baseline")][metric]))
                treatment_match = bool(int(indexed[(csd_id, treatment)][metric]))
                if baseline_match and treatment_match:
                    counts["stayed_match"] += 1
                elif baseline_match:
                    counts["regressed"] += 1
                elif treatment_match:
                    counts["fixed"] += 1
                else:
                    counts["stayed_mismatch"] += 1
            output.append(
                {
                    "treatment": treatment,
                    "metric": metric,
                    "fixed": counts["fixed"],
                    "regressed": counts["regressed"],
                    "net_gain": counts["fixed"] - counts["regressed"],
                    "stayed_match": counts["stayed_match"],
                    "stayed_mismatch": counts["stayed_mismatch"],
                }
            )
    return output


def factorial_rows(
    summary_rows: list[dict[str, str | int | float]],
) -> list[dict[str, str | float]]:
    rates = {
        (str(row["cohort"]), str(row["treatment"]), str(row["metric"])): float(
            row["fixed_denominator_percent"]
        )
        for row in summary_rows
    }
    output: list[dict[str, str | float]] = []
    cohorts = sorted({str(row["cohort"]) for row in summary_rows})
    for cohort in cohorts:
        for metric in METRICS:
            baseline = rates[(cohort, "baseline", metric)]
            ml_only = rates[(cohort, "ml_only", metric)]
            ligand_only = rates[(cohort, "ligand_only", metric)]
            full = rates[(cohort, "full", metric)]
            output.append(
                {
                    "cohort": cohort,
                    "metric": metric,
                    "baseline_percent": baseline,
                    "ml_edge_effect_pp": ml_only - baseline,
                    "ligand_edge_effect_pp": ligand_only - baseline,
                    "full_effect_pp": full - baseline,
                    "interaction_pp": full - ml_only - ligand_only + baseline,
                    "residual_mismatch_percent": 100.0 - full,
                }
            )
    return output


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def residual_connectivity_rows(
    details: list[dict[str, str | int]],
    validated_graph_ids: set[str],
    edge_audit: dict[str, dict[str, str]],
) -> tuple[list[dict[str, str | int]], Counter]:
    output: list[dict[str, str | int]] = []
    relations = Counter()
    for row in details:
        if (
            row["treatment"] != "full"
            or row["IDs"] not in validated_graph_ids
            or int(row["full_connectivity"]) == 1
        ):
            continue
        ilp_neighbors = str(row["neighbors_ilp"]).split()
        reference_neighbors = str(row["neighbors_reference"]).split()
        if len(ilp_neighbors) < len(reference_neighbors):
            relation = "fewer_output_metal_contacts"
        elif len(ilp_neighbors) > len(reference_neighbors):
            relation = "more_output_metal_contacts"
        elif ilp_neighbors != reference_neighbors:
            relation = "different_elements_same_count"
        else:
            relation = "same_shell_graph_difference"
        relations[relation] += 1
        audit = edge_audit.get(str(row["IDs"]), {})
        output.append(
            {
                **row,
                "residual_class": relation,
                "coordination_delta": len(ilp_neighbors) - len(reference_neighbors),
                "raw_ml_count": audit.get("raw_ml_count", ""),
                "reference_ml_count": audit.get("reference_ml_count", ""),
            }
        )
    return output, relations


def plot_summary(
    path: Path,
    summary_rows: list[dict[str, str | int | float]],
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    cohort = "validated_graph_consensus"
    metrics = list(METRICS)
    labels = [
        "Ligand canonical",
        "Ligand equivalent",
        "Metal connectivity",
        "Metal oxidation",
        "Full connectivity",
    ]
    rates = {
        (str(row["treatment"]), str(row["metric"])): float(
            row["fixed_denominator_percent"]
        )
        for row in summary_rows
        if row["cohort"] == cohort
    }
    cohort_size = next(
        int(row["total"]) for row in summary_rows if row["cohort"] == cohort
    )
    x = np.arange(len(metrics))
    width = 0.19
    fig, axis = plt.subplots(figsize=(11.5, 5.8))
    for index, treatment in enumerate(METHOD_FILES):
        values = [rates[(treatment, metric)] for metric in metrics]
        bars = axis.bar(
            x + (index - 1.5) * width,
            values,
            width,
            label=treatment.replace("_", " "),
        )
        axis.bar_label(bars, fmt="%.1f", fontsize=7, padding=2)
    axis.set_ylabel("Agreement on fixed cohort (%)")
    axis.set_xlabel("Agreement metric")
    axis.set_title("ILP agreement after controlled reference-graph replacement")
    axis.set_xticks(x, labels, rotation=12, ha="right")
    axis.set_ylim(75, 101)
    axis.legend(ncol=4, frameon=False, loc="lower center")
    axis.grid(axis="y", alpha=0.25)
    fig.text(
        0.01,
        0.01,
        f"Cohort: {cohort_size:,} structures with exact CSD/NBO/Hückel full-graph "
        "consensus and validated XYZ-indexed reference graphs; failures count as mismatches.",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(path, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment-dir", type=Path, default=DEFAULT_EXPERIMENT)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    metric_ids = load_ids(args.experiment_dir / "cohort_metric_consensus_29110.csv")
    graph_ids = load_ids(args.experiment_dir / "cohort_graph_consensus_29068.csv")
    with (args.experiment_dir / "reference_edges.csv").open(
        newline="",
        encoding="utf-8",
    ) as handle:
        validated_ids = {
            row["IDs"] for row in csv.DictReader(handle) if row.get("validated") == "1"
        }
    validated_graph_ids = graph_ids & validated_ids
    references = load_reference(args.reference)

    methods: dict[str, dict[str, str]] = {}
    statuses: dict[str, dict[str, str]] = {}
    for treatment, filename in METHOD_FILES.items():
        path = args.baseline if filename is None else args.experiment_dir / filename
        methods[treatment], statuses[treatment] = load_smiles_results(path)

    feature_cache, graph_cache = extract_caches(
        methods,
        references,
        metric_ids,
        args.workers,
    )
    details = evaluate_details(
        metric_ids,
        graph_ids,
        methods,
        statuses,
        references,
        feature_cache,
        graph_cache,
    )
    detail_path = args.experiment_dir / "agreement_per_structure.csv"
    write_detail(detail_path, details)

    summary_rows = summarize(details, metric_ids, graph_ids, validated_graph_ids)
    transition_summary = transition_rows(details, validated_graph_ids)
    factorial_summary = factorial_rows(summary_rows)
    write_csv(args.experiment_dir / "agreement_summary.csv", summary_rows)
    write_csv(args.experiment_dir / "transition_summary.csv", transition_summary)
    write_csv(args.experiment_dir / "factorial_attribution.csv", factorial_summary)

    residual_patterns = Counter()
    for row in details:
        if row["treatment"] != "full" or row["IDs"] not in validated_graph_ids:
            continue
        pattern = "|".join(
            f"{metric}={int(row[metric])}"
            for metric in (
                "ligand_equiv",
                "metal_connectivity",
                "metal_oxidation",
                "full_connectivity",
            )
        )
        residual_patterns[pattern] += 1
    residual_rows = [
        {"pattern": pattern, "count": count}
        for pattern, count in residual_patterns.most_common()
    ]
    write_csv(args.experiment_dir / "full_treatment_residual_patterns.csv", residual_rows)
    edge_audit = load_reference(args.experiment_dir / "edge_input_audit.csv")
    residual_connectivity, residual_relations = residual_connectivity_rows(
        details,
        validated_graph_ids,
        edge_audit,
    )
    write_csv(
        args.experiment_dir / "full_treatment_residual_connectivity.csv",
        residual_connectivity,
    )

    plot_path = args.experiment_dir / "reference_graph_agreement.png"
    plot_summary(plot_path, summary_rows)
    machine_summary = {
        "metric_cohort_size": len(metric_ids),
        "graph_cohort_size": len(graph_ids),
        "validated_graph_cohort_size": len(validated_graph_ids),
        "agreement_rows": len(details),
        "primary_summary": [
            row for row in summary_rows if row["cohort"] == "validated_graph_consensus"
        ],
        "primary_transitions": transition_summary,
        "primary_factorial_attribution": [
            row
            for row in factorial_summary
            if row["cohort"] == "validated_graph_consensus"
        ],
        "full_treatment_residual_patterns": residual_rows,
        "full_treatment_residual_connectivity_classes": dict(residual_relations),
    }
    (args.experiment_dir / "analysis_summary.json").write_text(
        json.dumps(machine_summary, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(machine_summary, indent=2))
    print(f"Wrote {detail_path}")
    print(f"Wrote {plot_path}")


if __name__ == "__main__":
    main()
