# Reference-graph ILP diagnostic

## Cohorts

- 29,110 structures: CSD, NBO(DFT), and Hückel(DFT) agree on strict ligand
  canonical structure, metal-neighbor element multiset, and metal oxidation state.
- 29,068 structures: the three references also have identical full heavy-atom
  connectivity.
- 28,948 structures: the reference graph was reproduced and mapped to DFT XYZ atom
  indices with the current Hückel workflow. This is the primary causal cohort.
- 120 exact-consensus structures were excluded from the causal cohort: 113 historical
  Hückel graphs were not reproduced by the current code and 7 XYZ files were absent.

## Controlled intervention

The experiment holds geometry, total charge, ILP parameters, aromaticity detection,
post-corrections, CBC classification, and SMILES generation fixed. It independently
replaces ligand-internal and metal–ligand (M–L) edges with the validated reference
graph.

Across the primary cohort, 918 structures have different M–L input edges, while zero
have different ligand-internal input edges. Consequently, the M–L-only and full-graph
treatments are identical. Ten unchanged-graph controls reproduced the baseline exact
SMILES and all semantic metrics, supporting reuse of baseline results when an input
edge set is unchanged.

## Main results

| Metric | Baseline ILP | Full reference graph | Change |
|---|---:|---:|---:|
| Ligand canonical | 85.021% | 85.205% | +0.183 pp |
| Chemically equivalent ligand | 95.703% | 95.882% | +0.180 pp |
| Metal connectivity | 96.608% | 99.717% | +3.109 pp |
| Metal oxidation state | 96.390% | 96.476% | +0.086 pp |
| Full heavy-atom connectivity | 96.601% | 99.710% | +3.109 pp |

Reference M–L edges fix 900 metal-connectivity mismatches with no regressions. This
shows that initial M–L bond perception explains most connectivity disagreement.

Only 25 net oxidation-state mismatches and 52 net chemically equivalent ligand
mismatches are corrected. After graph replacement, 1,020 oxidation mismatches remain;
1,012 also have a ligand mismatch. The dominant residual problem is therefore ILP
bond-order/formal-charge and covalent-versus-dative assignment, not connectivity.

## Residual connectivity

After full graph replacement, 84 full-connectivity mismatches remain:

- 82 outputs contain fewer metal contacts;
- 2 have the same metal-neighbor element multiset but a different donor identity;
- 83 of the 84 already had the reference graph as the original ILP input;
- the single changed-input residual (`PUJNIR`) becomes ILP-infeasible.

The reduced shells lose 85 C, 8 O, 4 Br, and 2 N contacts. Main metal hotspots are
Cu (20 structures), Ni (14), Pt (13), and Pd (6). These cases localize the remaining
connectivity problem to downstream ILP bond variables, CBC/dative inference, or
SMILES export.

## Key files

- `agreement_summary.csv`: fixed-denominator rates for all cohorts and treatments.
- `transition_summary.csv`: mismatch-to-match and regression counts.
- `factorial_attribution.csv`: M–L, ligand-edge, full, and interaction effects.
- `agreement_per_structure.csv`: per-structure metrics.
- `edge_input_audit.csv`: baseline and treatment edge sets.
- `full_treatment_residual_connectivity_enriched.csv`: residual connectivity cases.
- `residual_oxidation_by_metal_delta.csv`: residual oxidation-state hotspots.
- `ligand_mismatch_classification.csv`: per-structure ligand mismatch mechanism and
  coordination-motif class.
- `ligand_mismatch_by_motif.csv` and `ligand_mismatch_by_metal.csv`: normalized
  mismatch rates for ligand categories.
- `reference_graph_agreement.png`: summary plot.
