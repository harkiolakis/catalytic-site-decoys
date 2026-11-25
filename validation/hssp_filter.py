"""Sensitivity analysis (declared deviation): homology filter using the HSSP curve (Rost 1999):
a validation chain is dropped if, for any benchmark chain, local-alignment identity exceeds
pI(L) = 480 * L^(-0.32 * (1 + exp(-L/1000))) percent (L = aligned length; L>=11), with >=50% coverage of the shorter chain."""
import os, sys, json, time, pickle
import numpy as np, pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
cache = pickle.load(open(f"{HERE}/trace_cache.pkl", "rb"))
ent = pd.read_csv(f"{HERE}/validation_entries.csv").dropna(subset=["pdb"])
bp = pd.read_csv(f"{ROOT}/benchmark/pairs_924.csv")
bench = {(r.query.lower(), r.q_chain) for r in bp.itertuples()} | {(r.target.lower(), r.t_chain) for r in bp.itertuples()}
bseq = {k: cache[(f"{ROOT}/structures/{k[0]}.pdb", k[1])][1] for k in sorted(bench) if (f"{ROOT}/structures/{k[0]}.pdb", k[1]) in cache}
val = {(r.pdb, r.chain): cache[(f"{ROOT}/structures_validation/{r.pdb}.pdb", r.chain)][1] for r in ent.itertuples() if (f"{ROOT}/structures_validation/{r.pdb}.pdb", r.chain) in cache}
al = PairwiseAligner(); al.substitution_matrix = substitution_matrices.load("BLOSUM62"); al.open_gap_score = -10; al.extend_gap_score = -0.5; al.mode = "local"
def hssp(L): return 480 * L ** (-0.32 * (1 + np.exp(-L / 1000))) / 100 if L >= 11 else 1.0
def kmers(s, k=3): return {s[i:i + k] for i in range(len(s) - k + 1)}
bk = {k: kmers(v) for k, v in bseq.items()}
C = f"{HERE}/hssp_cache.json"; lc = json.load(open(C)) if os.path.exists(C) else {}; t0 = time.time()
for key, s in val.items():
    sk = f"{key[0]}|{key[1]}"
    if sk in lc: continue
    if time.time() - t0 > 150: break
    ks = kmers(s); best = (-9, 0, 0, None)
    for bkey, bs in bseq.items():
        if len(ks & bk[bkey]) / max(1, min(len(ks), len(bk[bkey]))) < 0.10: continue
        aln = al.align(s, bs)[0]; bl = list(zip(*aln.aligned))
        if not bl: continue
        same = sum(1 for (a1, e1), (b1, _) in bl for i in range(e1 - a1) if s[a1 + i] == bs[b1 + i])
        L = max(bl[-1][0][1] - bl[0][0][0], bl[-1][1][1] - bl[0][1][0]); cov = min(bl[-1][0][1] - bl[0][0][0], bl[-1][1][1] - bl[0][1][0]) / min(len(s), len(bs))
        if cov < 0.5: continue
        margin = same / L - hssp(L)
        if margin > best[0]: best = (float(margin), float(same / L), int(L), list(bkey))
    lc[sk] = best
json.dump(lc, open(C, "w"))
if len(lc) < len(val): print("done", len(lc), "/", len(val)); sys.exit(0)
drop = {tuple(k.split("|")) for k, v in lc.items() if v[0] > 0}
pairs = pd.read_csv(f"{HERE}/validation_pairs_final_preregistered.csv")  # has TM; rebuild from full list
full = pd.read_csv(f"{HERE}/validation_pairs.csv")
ok = full[[((r.query, r.q_chain) in val and (r.target, r.t_chain) in val and (r.query, r.q_chain) not in drop and (r.target, r.t_chain) not in drop) for r in full.itertuples()]].copy()
ok.to_csv(f"{HERE}/validation_pairs_hssp_noTM.csv", index=False)
with open(f"{HERE}/leakage_report_hssp.txt", "w") as f:
    f.write(f"HSSP-curve filter: chains dropped {len(drop)} of {len(val)}; pairs {len(full)} -> {len(ok)}\n")
    for k in sorted(drop): f.write(f"  {k} -> {lc[f'{k[0]}|{k[1]}'][3]} identity {lc[f'{k[0]}|{k[1]}'][1]:.2f} over {lc[f'{k[0]}|{k[1]}'][2]} aa\n")
print(open(f"{HERE}/leakage_report_hssp.txt").read()[:1500])
