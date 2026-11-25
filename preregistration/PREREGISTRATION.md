# CS-TOURIST v12 — independent validation: pre-registered analysis plan

Written **before** any validation structure was downloaded or any validation alignment was run. The frozen parameters, threshold and SHA-256 hashes of code and pair lists are in `FROZEN_CONFIG.json`. The method, parameters and threshold will not be changed after the validation results are seen. Any deviation will be reported as such.

## Validation set
- **Source:** M-CSA curated data (`mcsa/curated_data.csv`, downloaded 24 Nov 2026).
- **Independence:** an entry is excluded if any of its PDB entries **or** UniProt accessions occurs in the 924-pair benchmark (training set). This leaves 255 entries, of which **202** have a chain with ≥2 annotated catalytic residues.
- **Chain:** one per entry, the chain with the most annotated residues.
- **Pairs (seed 2026):** **1,118**:
  - all 218 pairs sharing ≥3 EC digits;
  - 400 pairs sharing the EC class only;
  - 500 pairs from different EC classes;
  - each entry used in at most 15 pairs.
- **Post-download leakage filter (fixed now):** drop any validation chain with >30% sequence identity (global alignment, BLOSUM62) to any benchmark chain, together with all of its pairs. Report how many were dropped.
- **Structural strata:** global TM-score (TM-align, normalised by the query) in place of CATH, with cut-offs < 0.3, 0.3–0.5 and ≥ 0.5.

## Methods compared
1. **v12** (frozen, `rerun/v12.py`, PARAMS as in `FROZEN_CONFIG.json`). This is the primary method.
2. **v11 Mode C** (original method, same code and settings as the rerun) for comparison.
3. **Same-input baselines:** random one-to-one mapping and lowest-RMSD unrestricted mapping.
- **Decoys:** random, surface, composition, compact (Rg-matched) and spatial (Cα-distance-matched), 5 trials per pair, as in the benchmark.

## Primary endpoints
- **P1. Specificity.** AUC of the v12 score, real vs **compact** decoys and real vs **spatial** decoys. Confidence intervals by pair-level bootstrap. *Success:* lower 95% CI bound > 0.60 for both.
- **P2. Functional correctness of significant matches.** Among pairs with v12 score ≥ **0.2487** (the frozen threshold), the mean fraction of aligned residue pairs sharing ≥1 M-CSA mechanistic role, compared with a random mapping of the same residues (paired Wilcoxon). *Success:* significantly greater (p < 0.05) **and** a difference ≥ 0.05.
- **P3. False-positive control.** The fraction of compact and spatial decoys scoring ≥ threshold. *Expected:* ≈ 5% (calibration check; *success* if ≤ 10%).

## Secondary endpoints
- **S1.** Role concordance of v12 vs v11 on pairs both align; and for v12 alignments with ≥3 matched residues.
- **S2.** Specificity vs random, surface and composition decoys.
- **S3.** Rate of significant matches and role concordance by EC stratum and TM-score stratum, including the count of significant matches between structurally dissimilar proteins (TM < 0.3).
- **S4.** Proportion of significant v12 alignments that are not sequence-ordered.

All endpoints are reported regardless of outcome.
