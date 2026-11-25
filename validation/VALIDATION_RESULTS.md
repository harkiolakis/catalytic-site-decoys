# Independent validation of Side-chain matcher (v12) — results against pre-registered endpoints
- hash check v12.py (preregistration/frozen_files): unchanged
- hash check run_v12.py (preregistration/frozen_files): unchanged
- hash check select_validation.py (preregistration/frozen_files): unchanged
- hash check validation_pairs.csv (preregistration/frozen_files): unchanged
- hash check validation_entries.csv (preregistration/frozen_files): unchanged
- hash check validation_sites.json (preregistration/frozen_files): unchanged

## BLAST E-value leakage filter (declared deviation; recommended primary)
Pairs analysable: 420 (v12 alignment returned for 77.9%)

### Primary endpoints
| endpoint                                                     |   pairs | value (95% CI)                                               | criterion         | met   |
|:-------------------------------------------------------------|--------:|:-------------------------------------------------------------|:------------------|:------|
| P1 AUC vs compact decoys                                     |     398 | 0.696 (0.668–0.723)                                          | lower CI > 0.60   | YES   |
| P1 AUC vs spatial decoys                                     |     416 | 0.701 (0.674–0.725)                                          | lower CI > 0.60   | YES   |
| P2 role concordance of significant matches vs random mapping |      82 | 0.428 (0.341–0.519) vs 0.368 (0.315–0.425); Δ=0.060; p=0.047 | p<0.05 and Δ≥0.05 | YES   |
| P3 false-positive rate, compact decoys                       |     420 | 4.7%                                                         | ≤10%              | YES   |
| P3 false-positive rate, spatial decoys                       |     420 | 3.7%                                                         | ≤10%              | YES   |
Significant matches (score ≥ 0.248651): 82 of 420 pairs (19.5%)

### Secondary endpoints
- S1 role concordance, pairs aligned by both (N=324): v12 0.363 (0.322–0.403) | v11 0.352 (0.316–0.388) | random 0.340 (0.311–0.368); v12 vs v11 p=0.47
- S1 v12 alignments with ≥3 residues (N=36): v12 0.334 (0.225–0.455) | v11 0.307 (0.223–0.389) | random 0.313 (0.239–0.393); p(v12>random)=0.74
- S2 other decoys:
| decoy       |   pairs | AUC (95% CI)        |
|:------------|--------:|:--------------------|
| random      |     420 | 0.815 (0.793–0.838) |
| surface     |     420 | 0.827 (0.805–0.848) |
| composition |     395 | 0.790 (0.768–0.814) |
- S3 by stratum:
| stratum                  |   pairs |   significant % | role (significant)   | random (same pairs)   |
|:-------------------------|--------:|----------------:|:---------------------|:----------------------|
| ec_stratum=diff_ec_class |     161 |            14.9 | 0.326 (0.188–0.500)  | 0.300 (0.211–0.397)   |
| ec_stratum=same_ec3      |      83 |            31.3 | 0.649 (0.477–0.796)  | 0.461 (0.365–0.548)   |
| ec_stratum=same_ec_class |     176 |            18.2 | 0.326 (0.203–0.448)  | 0.343 (0.251–0.441)   |
| tm_stratum=TM0.3-0.5     |     128 |            19.5 | 0.303 (0.170–0.450)  | 0.241 (0.159–0.335)   |
| tm_stratum=TM<0.3        |     282 |            19.1 | 0.469 (0.358–0.577)  | 0.424 (0.358–0.498)   |
| tm_stratum=TM>=0.5       |      10 |            30   | 0.733 (0.200–1.000)  | 0.423 (0.064–0.667)   |
- S3 significant matches between structurally dissimilar proteins (TM<0.3): 54
- S4 significant alignments not sequence-ordered: 50.0% (N=82)

Top significant matches between structurally dissimilar proteins (TM<0.3):
| query   | target   | ec_stratum    |    tm |   v12_k |   v12_score |   v12_rmsd_fc |   v12_role | v12_ordered   | v12_corr                |
|:--------|:---------|:--------------|------:|--------:|------------:|--------------:|-----------:|:--------------|:------------------------|
| 1k30    | 1cw0     | diff_ec_class | 0.234 |       2 |       0.876 |         0.913 |      1     | True          | 139-69;144-97           |
| 1mud    | 1s2k     | same_ec_class | 0.222 |       2 |       0.802 |         1.014 |      0     | False         | 37-136;138-53           |
| 1jm6    | 1cw0     | diff_ec_class | 0.206 |       2 |       0.724 |         1.537 |      0     | True          | 1239-69;1243-97         |
| 1mud    | 1dpx     | same_ec_class | 0.233 |       2 |       0.671 |         0.58  |      0     | True          | 37-35;138-52            |
| 1ru4    | 2abk     | same_ec_class | 0.212 |       2 |       0.581 |         1.246 |      0     | False         | 268-138;273-120         |
| 1cjy    | 1qd6     | same_ec3      | 0.174 |       3 |       0.519 |         0.937 |      0.667 | False         | 197-146;228-144;549-142 |
| 1jm6    | 1mj9     | same_ec_class | 0.252 |       2 |       0.516 |         0.702 |      0.5   | True          | 1239-304;1243-338       |
| 1k32    | 1j00     | same_ec_class | 0.119 |       2 |       0.498 |         0.169 |      1     | False         | 746-157;918-44          |
| 1ru4    | 1qaz     | same_ec3      | 0.287 |       2 |       0.465 |         1.561 |      1     | True          | 268-191;273-239         |
| 1dki    | 1qd6     | same_ec_class | 0.222 |       2 |       0.418 |         1.782 |      0.5   | False         | 47-144;195-142          |

## Pre-registered %-identity leakage filter (as registered)
Pairs analysable: 160 (v12 alignment returned for 77.5%)

### Primary endpoints
| endpoint                                                     |   pairs | value (95% CI)                                               | criterion         | met   |
|:-------------------------------------------------------------|--------:|:-------------------------------------------------------------|:------------------|:------|
| P1 AUC vs compact decoys                                     |     144 | 0.690 (0.641–0.735)                                          | lower CI > 0.60   | YES   |
| P1 AUC vs spatial decoys                                     |     156 | 0.713 (0.670–0.755)                                          | lower CI > 0.60   | YES   |
| P2 role concordance of significant matches vs random mapping |      30 | 0.486 (0.344–0.636) vs 0.385 (0.302–0.470); Δ=0.101; p=0.079 | p<0.05 and Δ≥0.05 | NO    |
| P3 false-positive rate, compact decoys                       |     160 | 5.5%                                                         | ≤10%              | YES   |
| P3 false-positive rate, spatial decoys                       |     160 | 2.2%                                                         | ≤10%              | YES   |
Significant matches (score ≥ 0.248651): 30 of 160 pairs (18.8%)

### Secondary endpoints
- S1 role concordance, pairs aligned by both (N=124): v12 0.404 (0.339–0.475) | v11 0.389 (0.333–0.445) | random 0.374 (0.329–0.416); v12 vs v11 p=0.82
- S1 v12 alignments with ≥3 residues (N=13): v12 0.276 (0.128–0.436) | v11 0.376 (0.237–0.505) | random 0.302 (0.196–0.402); p(v12>random)=0.68
- S2 other decoys:
| decoy       |   pairs | AUC (95% CI)        |
|:------------|--------:|:--------------------|
| random      |     160 | 0.814 (0.777–0.853) |
| surface     |     160 | 0.826 (0.785–0.865) |
| composition |     157 | 0.813 (0.777–0.851) |
- S3 by stratum:
| stratum                  |   pairs |   significant % | role (significant)   | random (same pairs)   |
|:-------------------------|--------:|----------------:|:---------------------|:----------------------|
| ec_stratum=diff_ec_class |      72 |            15.3 | 0.318 (0.091–0.545)  | 0.291 (0.192–0.389)   |
| ec_stratum=same_ec3      |      27 |            33.3 | 0.815 (0.667–0.944)  | 0.509 (0.380–0.648)   |
| ec_stratum=same_ec_class |      61 |            16.4 | 0.375 (0.125–0.600)  | 0.377 (0.216–0.567)   |
| tm_stratum=TM0.3-0.5     |      47 |            12.8 | 0.375 (0.042–0.708)  | 0.145 (0.037–0.257)   |
| tm_stratum=TM<0.3        |     111 |            20.7 | 0.493 (0.341–0.652)  | 0.439 (0.351–0.530)   |
| tm_stratum=TM>=0.5       |       2 |            50   | n/a                  | n/a                   |
- S3 significant matches between structurally dissimilar proteins (TM<0.3): 23
- S4 significant alignments not sequence-ordered: 43.3% (N=30)
