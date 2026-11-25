# Benchmark results (current M-CSA annotations)
Benchmark: 924 original pairs; 233 lack ≥2 annotated catalytic residues on both chains in current M-CSA; 4 produce no alignment; **687 aligned pairs** analysed (failures excluded from all RMSD statistics).
Site size: query median 4 (IQR 3–6), target median 4; aligned correspondences median 3.

## Table 1. Mode C (v11) site RMSD overall, by CATH relationship and by number of aligned residues
| set              |   N | median RMSD (95% CI)   | IQR       |   <2 Å % |   <3 Å % |   <5 Å % |
|:-----------------|----:|:-----------------------|:----------|---------:|---------:|---------:|
| All pairs        | 687 | 3.69 (3.37–3.98)       | 1.89–5.62 |     26.3 |     39.9 |     67.5 |
| Same HSF         | 156 | 2.42 (2.14–3.01)       | 0.97–4.31 |     39.7 |     57.1 |     84   |
| Same T, diff HSF | 138 | 4.19 (3.40–4.75)       | 1.98–5.77 |     25.4 |     36.2 |     64.5 |
| Different T      | 393 | 4.04 (3.68–4.28)       | 2.25–5.90 |     21.4 |     34.4 |     62.1 |
| 2 aligned        | 247 | 1.83 (1.43–2.12)       | 0.63–3.08 |     54.3 |     72.1 |     93.1 |
| 3 aligned        | 162 | 3.58 (3.12–4.11)       | 2.16–5.11 |     20.4 |     40.1 |     72.8 |
| 4 aligned        | 135 | 4.78 (4.50–5.32)       | 3.79–6.43 |      3.7 |      8.9 |     53.3 |
| >4 aligned       | 143 | 6.02 (5.61–6.62)       | 4.41–8.33 |      6.3 |     13.3 |     30.8 |
| ≥3 aligned       | 440 | 4.81 (4.51–5.13)       | 3.18–6.52 |     10.7 |     21.8 |     53.2 |

Kruskal–Wallis across CATH groups: H=34.56, p=3.13e-08, ε²=0.048
- Different T vs Same HSF: Δmedian 1.62 Å (95% CI 0.93–2.00), Cliff's δ 0.318, MW p=6.35e-09; 26% of Different T below Same HSF median
- Same T, diff HSF vs Same HSF: Δmedian 1.77 Å (95% CI 0.83–2.39), Cliff's δ 0.283, MW p=2.9e-05; 30% of Same T, diff HSF below Same HSF median
- Different T vs Same T, diff HSF: Δmedian -0.15 Å (95% CI -0.85–0.71), Cliff's δ 0.022, MW p=0.706; 52% of Different T below Same T, diff HSF median

## Table 2. Same-input baselines (same catalytic residues, same number of correspondences k as Mode C)
| set                                                |   N | median RMSD (95% CI)   | IQR       |   <2 Å % |   <3 Å % |   <5 Å % |
|:---------------------------------------------------|----:|:-----------------------|:----------|---------:|---------:|---------:|
| Octant-descriptor aligner (v11)                    | 687 | 3.69 (3.37–3.94)       | 1.89–5.62 |     26.3 |     39.9 |     67.5 |
| Random one-to-one mapping (median of 200)          | 687 | 4.31 (4.05–4.65)       | 2.48–6.38 |     18   |     31.6 |     58.8 |
| Random order-preserving mapping (median of 200)    | 687 | 4.02 (3.70–4.26)       | 2.29–5.75 |     19.5 |     34.6 |     64.6 |
| Best order-preserving mapping (exhaustive)         | 687 | 1.78 (1.58–2.11)       | 0.54–3.71 |     52.1 |     67.2 |     87.2 |
| Best unrestricted mapping (exhaustive/100k sample) | 687 | 1.30 (1.11–1.46)       | 0.47–2.51 |     66.4 |     82.1 |     94.5 |

Mode C vs random mapping: paired Wilcoxon p=1.6e-18; Mode C lower in 63.5% of pairs.
Mode C finds the best order-preserving mapping in 29.5% of pairs; median percentile of Mode C among all mappings: 32.6% beat it.
The best unrestricted mapping is NOT sequence-order-preserving and ≥0.5 Å better than Mode C in 49.8% of pairs (53.4% of different-topology pairs) — the order-preserving Smith–Waterman cannot recover permuted sites.

## Table 3. Ablation of Mode C components (same pairs; each configuration's own k)
| configuration                                                             |   aligned |   median k |   median RMSD |   <2 Å % |   same k as full: N |   identical correspondence RMSD % | median RMSD at same k (config vs full)   |
|:--------------------------------------------------------------------------|----------:|-----------:|--------------:|---------:|--------------------:|----------------------------------:|:-----------------------------------------|
| Full v11 (72-slot conn. + 24-slot exact + BLOSUM + distance + refinement) |       687 |          3 |          3.69 |     26.3 |                 687 |                             100   | 3.69 vs 3.69                             |
| 72-slot connectivity channel only (no exact-offset channel)               |       687 |          3 |          3.71 |     27.2 |                 666 |                              77.6 | 3.57 vs 3.63                             |
| v10 exact-offset descriptor (2-sphere)                                    |       613 |          3 |          3.33 |     27.2 |                 477 |                              59.1 | 3.33 vs 3.21                             |
| v11 without BLOSUM                                                        |       685 |          3 |          3.61 |     25.3 |                 662 |                              79   | 3.50 vs 3.63                             |
| v11 without distance-compatibility term                                   |       687 |          3 |          3.74 |     24.6 |                 677 |                              91.7 | 3.71 vs 3.68                             |
| v11 without iterative refinement                                          |       687 |          3 |          3.69 |     26.2 |                 686 |                              99.7 | 3.69 vs 3.69                             |
| No descriptor (BLOSUM + distance)                                         |       509 |          3 |          2.92 |     24.6 |                 326 |                              47.2 | 2.98 vs 3.04                             |
| Distance term only                                                        |       633 |          3 |          2.94 |     35.5 |                 546 |                              34.4 | 2.86 vs 3.41                             |
| BLOSUM only                                                               |       312 |          2 |          2.54 |     16.9 |                 129 |                              57.4 | 2.72 vs 2.33                             |

## Table 4. Functional correctness: fraction of aligned residue pairs sharing ≥1 M-CSA mechanistic role
| method                           |   N pairs | mean role concordance (95% CI)   |   vs random (paired Wilcoxon p) |
|:---------------------------------|----------:|:---------------------------------|--------------------------------:|
| Mode C (v11)                     |       687 | 0.414 (0.387–0.440)              |                         2.9e-05 |
| v10 descriptor                   |       613 | 0.429 (0.402–0.455)              |                         1.3e-07 |
| Distance term only               |       633 | 0.396 (0.370–0.423)              |                         0.1     |
| BLOSUM only                      |       312 | 0.460 (0.420–0.500)              |                         1.2e-07 |
| Random mapping (mean of 50)      |       687 | 0.382 (0.361–0.403)              |                                 |
| Lowest-RMSD unrestricted mapping |       687 | 0.395 (0.368–0.421)              |                         0.088   |
- Same HSF: Mode C 0.488 vs random 0.410 (N=156)
- Same T, diff HSF: Mode C 0.415 vs random 0.402 (N=138)
- Different T: Mode C 0.384 vs random 0.363 (N=393)

## Side-chain (functional-group centroid) RMSD under the Mode C correspondence
| set                      |   N | median RMSD (95% CI)   | IQR       |   <2 Å % |   <3 Å % |   <5 Å % |
|:-------------------------|----:|:-----------------------|:----------|---------:|---------:|---------:|
| Cα RMSD                  | 687 | 3.69 (3.40–3.94)       | 1.89–5.62 |     26.3 |     39.9 |     67.5 |
| Side-chain centroid RMSD | 687 | 3.09 (2.94–3.29)       | 1.64–4.77 |     31.9 |     47.5 |     77   |
Spearman(Cα, side-chain RMSD) = 0.87

## Table 5. Decoy controls (target catalytic residues replaced; 5 trials per pair per mode; pair-level bootstrap CIs)
| decoy       | stratum   |   pairs |   real median |   decoy median | AUC (95% CI)        |   Wilcoxon p (per pair) |   real <2 Å % |   decoy <2 Å % | Rg real / decoy (Å)   |
|:------------|:----------|--------:|--------------:|---------------:|:--------------------|------------------------:|--------------:|---------------:|:----------------------|
| random      | all       |     687 |          3.69 |          10.53 | 0.846 (0.833–0.859) |                  2e-109 |          26.3 |            6.3 | 6.8 / 16.7            |
| random      | k=2       |     247 |          1.83 |           7.78 | 0.840 (0.816–0.862) |                  3e-39  |          54.3 |           16.3 | 6.0 / 14.2            |
| random      | k≥3       |     440 |          4.81 |          12.08 | 0.887 (0.870–0.903) |                  8e-72  |          10.7 |            1   | 7.1 / 17.8            |
| surface     | all       |     687 |          3.69 |          12.89 | 0.891 (0.881–0.902) |                  7e-113 |          26.3 |            4   | 6.8 / 19.1            |
| surface     | k=2       |     247 |          1.83 |          10.1  | 0.888 (0.869–0.905) |                  1e-41  |          54.3 |            9.8 | 6.0 / 17.0            |
| surface     | k≥3       |     440 |          4.81 |          14.5  | 0.927 (0.914–0.939) |                  3e-73  |          10.7 |            0.9 | 7.1 / 20.1            |
| sequence    | all       |     687 |          3.69 |           9.11 | 0.774 (0.757–0.790) |                  5e-99  |          26.3 |           10.6 | 6.8 / 14.8            |
| sequence    | k=2       |     247 |          1.83 |           5.73 | 0.745 (0.716–0.775) |                  7e-33  |          54.3 |           26.4 | 6.0 / 13.1            |
| sequence    | k≥3       |     440 |          4.81 |          10.32 | 0.827 (0.809–0.848) |                  6e-67  |          10.7 |            2   | 7.1 / 15.7            |
| composition | all       |     646 |          3.62 |          11.07 | 0.859 (0.847–0.872) |                  8e-104 |          27.1 |            4.9 | 6.7 / 17.0            |
| composition | k=2       |     242 |          1.79 |           8.34 | 0.862 (0.838–0.885) |                  4e-39  |          54.5 |           12.1 | 5.9 / 15.1            |
| composition | k≥3       |     404 |          4.8  |          12.66 | 0.900 (0.884–0.916) |                  2e-66  |          10.6 |            0.7 | 7.1 / 18.0            |
| compact     | all       |     635 |          3.46 |           3.51 | 0.505 (0.490–0.522) |                  0.7    |          28.2 |           26.1 | 6.6 / 5.9             |
| compact     | k=2       |     229 |          1.62 |           1.7  | 0.504 (0.466–0.541) |                  1      |          57.6 |           59.3 | 5.8 / 5.1             |
| compact     | k≥3       |     406 |          4.6  |           4.61 | 0.497 (0.475–0.521) |                  0.6    |          11.6 |            8.3 | 6.9 / 6.3             |
| spatial     | all       |     686 |          3.69 |           3.91 | 0.525 (0.511–0.538) |                  0.003  |          26.4 |           22.7 | 6.7 / 6.7             |
| spatial     | k=2       |     247 |          1.83 |           1.98 | 0.523 (0.494–0.555) |                  0.07   |          54.3 |           51.9 | 6.0 / 6.0             |
| spatial     | k≥3       |     439 |          4.81 |           4.99 | 0.527 (0.507–0.546) |                  0.02   |          10.7 |            6.5 | 7.0 / 7.0             |
| spatial_seq | all       |     509 |          3.32 |           3.45 | 0.515 (0.498–0.533) |                  0.06   |          30.3 |           26.9 | 6.3 / 6.6             |
| spatial_seq | k=2       |     217 |          1.87 |           2.14 | 0.521 (0.488–0.550) |                  0.06   |          53   |           50.2 | 5.8 / 5.8             |
| spatial_seq | k≥3       |     292 |          4.51 |           4.43 | 0.511 (0.484–0.535) |                  0.3    |          13.4 |            9.4 | 6.7 / 6.9             |

Pairs with ≥1 valid decoy: compact 92%, composition 94%, random 100%, sequence 100%, spatial 100%, spatial_seq 74%, surface 100%
Spatial decoys: distance tolerance used — 1 Å: 83%, 1.5 Å: 11%, 2 Å: 5%, 3 Å: 1%

## EC analysis (current M-CSA EC numbers)
Pairs with EC on both sides: 687/687
|   ec_shared |   N |   median_RMSD |   pct_lt2 |   role_conc |
|------------:|----:|--------------:|----------:|------------:|
|           0 | 494 |         3.88  |      23.5 |       0.395 |
|           1 | 119 |         3.065 |      31.1 |       0.44  |
|           2 |  26 |         3.413 |      19.2 |       0.543 |
|           3 |  46 |         2.302 |      47.8 |       0.465 |
|           4 |   2 |         1.328 |      50   |       0.7   |
ec_shared        0    1   2   3  4  All
cath_relation                          
same_h          81   27   9  37  2  156
same_t_diff_h  110   22   4   2  0  138
diff_t         303   70  13   7  0  393
All            494  119  26  46  2  687
Different-topology pairs sharing EC first digit: 90 of 393 (manuscript §2.4 states none)
Spearman(shared EC digits, RMSD) = -0.148, p=0.0001
Sub-2 Å pairs: 1% identical EC, 64% different EC class; of the latter, 80% are 2-residue alignments
Identical-EC pairs NOT reaching <5 Å: 0% (N=2)

## Global TM-score vs site RMSD
Spearman(TM-score, site RMSD) = -0.194
Sub-2 Å different-topology pairs: 84; with k=2: 69; TM-score ≥0.5: 1

Best different-topology alignments with ≥4 residues (case-study candidates):
| query   | target   |   k |   v11_full_rmsd |   v11_sc_rmsd |   v11_full_role |    tm |   ec_shared |
|:--------|:---------|----:|----------------:|--------------:|----------------:|------:|------------:|
| 1s9c    | 1b3r     |   4 |           1.264 |         1.568 |            1    | 0.345 |           0 |
| 1deu    | 1ey3     |   4 |           1.565 |         0.766 |            0.25 | 0.227 |           0 |
| 1g64    | 1mkb     |   4 |           2.612 |         2.078 |            0.75 | 0.307 |           0 |
| 1eug    | 1bt3     |   4 |           2.633 |         1.873 |            0    | 0.326 |           0 |
| 5fit    | 1gcb     |   4 |           2.871 |         3.568 |            0.75 | 0.339 |           1 |
| 1x9d    | 1bzy     |   4 |           2.876 |         1.673 |            0.25 | 0.205 |           0 |
| 1or8    | 1zp3     |   4 |           2.895 |         2.198 |            0    | 0.266 |           0 |
| 5cnv    | 1mrq     |   4 |           2.954 |         1.813 |            0.25 | 0.231 |           1 |
| 1exp    | 1q0n     |   4 |           3.032 |         2.161 |            0.25 | 0.216 |           0 |
| 1lok    | 1gcb     |   4 |           3.084 |         3.059 |            0    | 0.278 |           2 |
