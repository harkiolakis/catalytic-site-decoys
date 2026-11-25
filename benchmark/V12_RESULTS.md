# Side-chain matcher (v12) (permutation-tolerant, chemistry-aware) vs v11 Mode C — same 691 pairs, current M-CSA
Pairs: 691. v12 returns an alignment (≥2 matched residues within 3 Å on both Cα and functional group) for 81.9%; v11 Mode C for 99.4% (v11 has no rejection threshold).
v12 matched residues: median 2; k≥3 in 9.8% of pairs. Non-sequence-ordered alignments: 44.7%.

## 1. Functional correctness (share ≥1 M-CSA mechanistic role per aligned residue pair)
Pairs aligned by both v12 and v11 with role data: 566
| method                        |   N | mean role concordance (95% CI)   |   vs v11 (paired Wilcoxon p) |
|:------------------------------|----:|:---------------------------------|-----------------------------:|
| v12                           | 566 | 0.424 (0.390–0.456)              |                      0.1     |
| v12 + connectivity descriptor | 566 | 0.421 (0.389–0.452)              |                      0.13    |
| v12 without chemistry         | 566 | 0.382 (0.352–0.412)              |                      0.047   |
| v12 Cα only                   | 566 | 0.431 (0.400–0.459)              |                      0.02    |
| v11 Mode C                    | 566 | 0.409 (0.381–0.435)              |                              |
| Lowest-RMSD mapping (Cα)      | 566 | 0.391 (0.365–0.416)              |                      0.076   |
| Random mapping                | 566 | 0.380 (0.356–0.403)              |                      0.00033 |
- same_h: v12 0.534 | v11 0.501 | random 0.416 (N=132)
- same_t_diff_h: v12 0.388 | v11 0.400 | random 0.390 (N=112)
- diff_t: v12 0.391 | v11 0.375 | random 0.361 (N=322)

All v12 alignments (N=566): role concordance 0.424 (0.389–0.452); for k≥3 0.607 (0.529–0.684)

## 2. Geometry of v12 alignments
- v12 Cα RMSD: median 1.09 Å, <2 Å 95.9% (N=566)
- v12 functional-group RMSD: median 1.22 Å, <2 Å 92.4% (N=566)
- v11 Cα RMSD (same pairs): median 3.89 Å, <2 Å 24.7% (N=566)

## 3. Specificity: real site vs decoy target sites (score = normalised alignment score; decoys with no alignment score 0)
| method   | decoy       | stratum          |   pairs |   real mean score |   decoy mean score | AUC (95% CI)        |   Wilcoxon p |   decoy aligned % |   real aligned % |
|:---------|:------------|:-----------------|--------:|------------------:|-------------------:|:--------------------|-------------:|------------------:|-----------------:|
| v12      | random      | all              |     691 |             0.165 |              0.036 | 0.818 (0.800–0.836) |        6e-90 |              36.8 |             81.9 |
| v12      | random      | min site size ≥3 |     450 |             0.202 |              0.05  | 0.890 (0.872–0.907) |        2e-69 |              51.3 |             96   |
| v12      | random      | min site size 2  |     241 |             0.096 |              0.009 | 0.741 (0.709–0.774) |        1e-21 |               9.8 |             55.6 |
| v12      | surface     | all              |     691 |             0.165 |              0.032 | 0.826 (0.809–0.844) |        4e-92 |              30.9 |             81.9 |
| v12      | surface     | min site size ≥3 |     450 |             0.202 |              0.046 | 0.891 (0.875–0.907) |        3e-70 |              43.7 |             96   |
| v12      | surface     | min site size 2  |     241 |             0.096 |              0.006 | 0.753 (0.722–0.787) |        2e-23 |               6.9 |             55.6 |
| v12      | composition | all              |     650 |             0.163 |              0.037 | 0.810 (0.791–0.826) |        3e-85 |              31.5 |             81.1 |
| v12      | composition | min site size ≥3 |     413 |             0.202 |              0.052 | 0.875 (0.857–0.892) |        6e-65 |              44.6 |             95.6 |
| v12      | composition | min site size 2  |     237 |             0.097 |              0.011 | 0.741 (0.708–0.772) |        3e-21 |               8.8 |             55.7 |
| v12      | compact     | all              |     642 |             0.17  |              0.09  | 0.685 (0.665–0.706) |        6e-49 |              71.3 |             83   |
| v12      | compact     | min site size ≥3 |     418 |             0.207 |              0.117 | 0.738 (0.711–0.764) |        2e-41 |              89   |             96.4 |
| v12      | compact     | min site size 2  |     224 |             0.102 |              0.041 | 0.638 (0.597–0.681) |        2e-10 |              38.2 |             58   |
| v12      | spatial     | all              |     682 |             0.166 |              0.084 | 0.690 (0.671–0.711) |        1e-51 |              67.7 |             81.7 |
| v12      | spatial     | min site size ≥3 |     442 |             0.204 |              0.11  | 0.747 (0.722–0.771) |        8e-44 |              85.7 |             95.9 |
| v12      | spatial     | min site size 2  |     240 |             0.097 |              0.038 | 0.639 (0.602–0.675) |        1e-10 |              34.6 |             55.4 |
| v12_desc | random      | all              |     691 |             0.212 |              0.048 | 0.822 (0.804–0.839) |        5e-91 |              36.8 |             81.9 |
| v12_desc | random      | min site size ≥3 |     450 |             0.258 |              0.066 | 0.896 (0.878–0.913) |        5e-70 |              51.3 |             96   |
| v12_desc | random      | min site size 2  |     241 |             0.125 |              0.013 | 0.741 (0.708–0.776) |        1e-22 |               9.8 |             55.6 |
| v12_desc | surface     | all              |     691 |             0.212 |              0.04  | 0.836 (0.819–0.853) |        9e-94 |              30.9 |             81.9 |
| v12_desc | surface     | min site size ≥3 |     450 |             0.258 |              0.057 | 0.905 (0.890–0.921) |        2e-71 |              43.7 |             96   |
| v12_desc | surface     | min site size 2  |     241 |             0.125 |              0.008 | 0.754 (0.721–0.788) |        4e-24 |               6.9 |             55.6 |
| v12_desc | composition | all              |     650 |             0.209 |              0.046 | 0.818 (0.800–0.834) |        1e-86 |              31.5 |             81.1 |
| v12_desc | composition | min site size ≥3 |     413 |             0.258 |              0.065 | 0.887 (0.871–0.902) |        1e-65 |              44.6 |             95.6 |
| v12_desc | composition | min site size 2  |     237 |             0.126 |              0.014 | 0.743 (0.708–0.774) |        1e-22 |               8.8 |             55.7 |
| v12_desc | compact     | all              |     642 |             0.218 |              0.116 | 0.693 (0.670–0.714) |        1e-51 |              71.3 |             83   |
| v12_desc | compact     | min site size ≥3 |     418 |             0.264 |              0.149 | 0.751 (0.725–0.778) |        9e-45 |              89   |             96.4 |
| v12_desc | compact     | min site size 2  |     224 |             0.132 |              0.055 | 0.640 (0.600–0.680) |        8e-11 |              38.2 |             58   |
| v12_desc | spatial     | all              |     682 |             0.213 |              0.112 | 0.692 (0.672–0.714) |        1e-53 |              67.7 |             81.7 |
| v12_desc | spatial     | min site size ≥3 |     442 |             0.26  |              0.144 | 0.750 (0.724–0.774) |        4e-46 |              85.7 |             95.9 |
| v12_desc | spatial     | min site size 2  |     240 |             0.125 |              0.052 | 0.639 (0.602–0.675) |        5e-11 |              34.6 |             55.4 |

For reference, v11 Mode C (RMSD-based): compact AUC 0.505, spatial 0.525, random 0.846.

## 4. Ablations (real pairs)
| variant                           |   aligned % |   median k |   k≥3 % | role concordance    |   median Cα RMSD |
|:----------------------------------|------------:|-----------:|--------:|:--------------------|-----------------:|
| v12                               |        81.9 |          2 |     9.8 | 0.424 (0.388–0.453) |             1.09 |
| + connectivity descriptor (w=0.5) |        81.9 |          2 |    12.4 | 0.421 (0.392–0.452) |             1.08 |
| no chemistry                      |        81.9 |          2 |    25.9 | 0.382 (0.350–0.413) |             1.1  |
| Cα only (no functional group)     |        91.9 |          2 |    23.2 | 0.425 (0.397–0.457) |             0.4  |

## 5. Evolutionary relationship (v12 score)
Mean score: same_h 0.239, same_t_diff_h 0.151, diff_t 0.141; Kruskal–Wallis p=3.9e-07
Aligned %: same_h 85%, same_t_diff_h 81%, diff_t 81%

## 6. Top different-topology v12 alignments with ≥3 matched residues (case-study candidates)
| query   | target   |   v12_k |   v12_score |   v12_rmsd_ca |   v12_rmsd_fc |   v12_role | v12_ordered   | v12_corr                     |
|:--------|:---------|--------:|------------:|--------------:|--------------:|-----------:|:--------------|:-----------------------------|
| 1qh5    | 4x9e     |       4 |       0.445 |         0.959 |         1.36  |      1     | False         | 58-118;59-117;134-268;173-69 |
| 1bs1    | 1p3d     |       3 |       0.354 |         0.601 |         0.357 |      1     | True          | 15-129;16-130;115-173        |
| 1qib    | 1uw8     |       3 |       0.339 |         0.563 |         1.009 |      0.667 | True          | 201-95;205-97;211-162        |
| 1itq    | 2ybp     |       3 |       0.32  |         1.235 |         1.199 |      0.667 | False         | 125-190;152-188;198-276      |
| 1rhs    | 1pq5     |       3 |       0.279 |         0.567 |         1.065 |      1     | True          | 250-193;251-194;252-195      |
| 2phk    | 1fr8     |       3 |       0.262 |         1.342 |         1.594 |      0     | False         | 154-318;167-317;186-252      |
| 1ca2    | 1nia     |       3 |       0.254 |         1.305 |         1.5   |      1     | False         | 94-150;96-95;119-145         |
| 1jnr    | 1oth     |       3 |       0.25  |         1.346 |         1.749 |      1     | False         | 141-263;265-330;398-141      |
| 1vr7    | 1azw     |       3 |       0.229 |         0.898 |         1.268 |      0     | False         | 62-111;63-110;83-43          |
| 1djx    | 1dqs     |       3 |       0.228 |         1.564 |         1.404 |      0.667 | False         | 312-268;356-275;390-194      |
| 1mpx    | 4hpn     |       3 |       0.22  |         1.332 |         1.156 |      0.333 | True          | 82-166;307-269;340-296       |
| 1zoi    | 1l6p     |       3 |       0.217 |         1.012 |         1.533 |      0.667 | False         | 97-109;227-68;256-70         |
