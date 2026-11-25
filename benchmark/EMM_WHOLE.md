# Enzyme Motif Miner in native whole-structure mode

EnzyMM 0.5.0 with its default Jess parameters and templates of ≥3 side-chain residues, searched against the complete target structure (all chains, first model). 'Filtered' = matches the EnzyMM logistic-regression filter predicts correct. A match is on-site when any matched residue is an annotated M-CSA catalytic residue of the target (residue number, any chain).

## Benchmark

### (A) Pair task: query-entry templates vs whole target structure (691 pairs; 572 with ≥1 usable template)

| outcome                              | query-entry templates   | unrelated-entry templates (2 per pair, different EC class)   |
|:-------------------------------------|:------------------------|:-------------------------------------------------------------|
| Any match (unfiltered)               | 2.4% [1.5–4.1] (14/572) | 0.3% [0.1–0.8] (3/1144)                                      |
| Filtered match                       | 1.7% [1.0–3.2] (10/572) | 0.0% [0.0–0.3] (0/1144)                                      |
| Filtered match on the catalytic site | 1.7% [1.0–3.2] (10/572) | 0.0% [0.0–0.3] (0/1144)                                      |

Filtered matches: median 4 residues, RMSD 0.41 Å; precision (share of matched residues annotated) 0.82; recall of the annotated site 0.71.

Filtered on-site match rate by CATH relation: diff_t 0.0% [0.0–1.2] (0/323), same_h 7.1% [3.9–12.6] (10/141), same_t_diff_h 0.0% [0.0–3.4] (0/108)

### (B) Annotation task: whole library minus the target's own entry vs whole target structure (435 structures)

| outcome                                                    | value                       |
|:-----------------------------------------------------------|:----------------------------|
| ≥1 filtered match                                          | 41.4% [36.8–46.1] (180/435) |
| Top filtered match on the catalytic site                   | 33.3% [29.1–37.9] (145/435) |
| Any filtered match on the catalytic site                   | 33.3% [29.1–37.9] (145/435) |
| Top filtered match off-site (residue-level false positive) | 8.0% [5.8–11.0] (35/435)    |

Top filtered match: precision 0.60, recall 0.54; matching template's entry shares EC to level ≥3 with the target in 38.3% [31.5–45.6] (69/180); level 1 in 65.0% [57.8–71.6] (117/180).

## Validation (BLAST-filtered set)

### (A) Pair task: query-entry templates vs whole target structure (451 pairs; 337 with ≥1 usable template)

| outcome                              | query-entry templates   | unrelated-entry templates (2 per pair, different EC class)   |
|:-------------------------------------|:------------------------|:-------------------------------------------------------------|
| Any match (unfiltered)               | 1.2% [0.5–3.0] (4/337)  | 0.3% [0.1–1.1] (2/674)                                       |
| Filtered match                       | 0.6% [0.2–2.1] (2/337)  | 0.0% [0.0–0.6] (0/674)                                       |
| Filtered match on the catalytic site | 0.6% [0.2–2.1] (2/337)  | 0.0% [0.0–0.6] (0/674)                                       |

Filtered matches: median 6 residues, RMSD 0.15 Å; precision (share of matched residues annotated) 0.65; recall of the annotated site 1.00.

Filtered on-site match rate by EC stratum: diff_ec_class 0.0% [0.0–2.6] (0/142), same_ec3 3.3% [0.9–11.4] (2/60), same_ec_class 0.0% [0.0–2.8] (0/135)

### (B) Annotation task: whole library minus the target's own entry vs whole target structure (124 structures)

| outcome                                                    | value                      |
|:-----------------------------------------------------------|:---------------------------|
| ≥1 filtered match                                          | 22.6% [16.1–30.7] (28/124) |
| Top filtered match on the catalytic site                   | 16.1% [10.7–23.6] (20/124) |
| Any filtered match on the catalytic site                   | 16.1% [10.7–23.6] (20/124) |
| Top filtered match off-site (residue-level false positive) | 6.5% [3.3–12.2] (8/124)    |

Top filtered match: precision 0.58, recall 0.46; matching template's entry shares EC to level ≥3 with the target in 28.6% [15.3–47.1] (8/28); level 1 in 64.3% [45.8–79.3] (18/28).
