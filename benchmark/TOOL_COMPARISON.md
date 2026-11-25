# JESS vs Side-chain matcher (v12) under identical inputs
Same pairs, same annotated residues on both sides, same decoy residue sets (identical generators and seeds), same role-concordance evaluation.
JESS: templates from query catalytic residues (Cα + functional atom; chemical-class residue matching), searched only against the target's annotated or decoy residues; subtemplates of 2–4 residues; RMSD ≤ 2.0 Å, distance cut-off 1.5 Å; score = −log10 E-value of the best (largest, then lowest-E) hit.

## Benchmark (training, 691 pairs): 691 pairs
### Specificity (score-based AUC, real vs decoy)
| decoy family   | v12 AUC (95% CI)    | JESS AUC (95% CI)   |
|:---------------|:--------------------|:--------------------|
| random         | 0.818 (0.800–0.836) | 0.588 (0.573–0.603) |
| surface        | 0.826 (0.809–0.844) | 0.589 (0.575–0.603) |
| composition    | 0.810 (0.791–0.827) | 0.580 (0.567–0.595) |
| compact        | 0.685 (0.664–0.706) | 0.586 (0.571–0.602) |
| spatial        | 0.690 (0.668–0.711) | 0.579 (0.565–0.593) |

### Calls at matched false-positive rate (threshold = 95th / 99th percentile of compact+spatial decoy scores in this set)
| tool   |   decoy pct | decoy FP rate   | real pairs called   | role concordance of calls   | random mapping, same pairs   |       p |
|:-------|------------:|:----------------|:--------------------|:----------------------------|:-----------------------------|--------:|
| v12    |          95 | 5.0%            | 160 (23.2%)         | 0.543 (0.484–0.599)         | 0.439 (0.399–0.476)          | 8.4e-06 |
| JESS   |          95 | 2.1%            | 126 (18.2%)         | 0.548 (0.487–0.605)         | 0.407 (0.365–0.449)          | 3.2e-06 |
| v12    |          99 | 1.0%            | 58 (8.4%)           | 0.606 (0.523–0.698)         | 0.446 (0.387–0.507)          | 0.0004  |
| JESS   |          99 | 1.0%            | 75 (10.9%)          | 0.561 (0.479–0.640)         | 0.380 (0.335–0.427)          | 1.7e-05 |

- JESS returns any hit for 18.2% of real pairs vs 2.0% of compact/spatial decoys; v12 aligns 81.9% of real pairs.
- On pairs where JESS hits (N=125): JESS roles 0.553 (0.489–0.622), v12 roles 0.544 (0.485–0.608), random 0.408 (0.362–0.449)
- JESS hit sizes: 2 res: 102, 3 res: 17, 4 res: 7; non-sequence-ordered JESS hits: 35.7%
- Overlap: pairs called by both 90; JESS only 36; v12 (95th pct) only 70

## Independent validation (BLAST-filtered): 420 pairs
### Specificity (score-based AUC, real vs decoy)
| decoy family   | v12 AUC (95% CI)    | JESS AUC (95% CI)   |
|:---------------|:--------------------|:--------------------|
| random         | 0.815 (0.792–0.839) | 0.558 (0.544–0.574) |
| surface        | 0.827 (0.804–0.850) | 0.557 (0.541–0.571) |
| composition    | 0.790 (0.767–0.815) | 0.546 (0.532–0.561) |
| compact        | 0.696 (0.670–0.721) | 0.551 (0.535–0.568) |
| spatial        | 0.701 (0.675–0.728) | 0.551 (0.536–0.565) |

### Calls at matched false-positive rate (threshold = 95th / 99th percentile of compact+spatial decoy scores in this set)
| tool   |   decoy pct | decoy FP rate   | real pairs called   | role concordance of calls   | random mapping, same pairs   |     p |
|:-------|------------:|:----------------|:--------------------|:----------------------------|:-----------------------------|------:|
| v12    |          95 | 5.0%            | 94 (22.4%)          | 0.416 (0.333–0.496)         | 0.367 (0.313–0.421)          | 0.07  |
| JESS   |          95 | 1.3%            | 49 (11.7%)          | 0.369 (0.258–0.469)         | 0.312 (0.248–0.378)          | 0.22  |
| v12    |          99 | 1.0%            | 37 (8.8%)           | 0.557 (0.418–0.680)         | 0.440 (0.349–0.524)          | 0.014 |
| JESS   |          99 | 1.0%            | 44 (10.5%)          | 0.366 (0.252–0.474)         | 0.318 (0.254–0.378)          | 0.32  |

- JESS returns any hit for 11.7% of real pairs vs 1.2% of compact/spatial decoys; v12 aligns 77.9% of real pairs.
- On pairs where JESS hits (N=48): JESS roles 0.366 (0.257–0.477), v12 roles 0.339 (0.241–0.442), random 0.308 (0.240–0.375)
- JESS hit sizes: 2 res: 46, 3 res: 2, 4 res: 1; non-sequence-ordered JESS hits: 28.6%
- Overlap: pairs called by both 25; JESS only 24; v12 (95th pct) only 69
