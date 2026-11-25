"""Leakage filter (declared deviation from the pre-registered %-identity rule, which proved to sit at the
best-of-765 noise floor): BLAST-style local alignment (BLOSUM62, gap open 11 / extend 1) of each validation
chain against all benchmark chains; Karlin-Altschul E-value with database = all benchmark chains
(lambda=0.267, K=0.041, gapped BLOSUM62 11/1). Drop chains whose best hit has E < 1e-3."""
import os, sys, json, time, pickle, math
import numpy as np, pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
cache = pickle.load(open(f"{HERE}/trace_cache.pkl", "rb"))
ent = pd.read_csv(f"{HERE}/validation_entries.csv").dropna(subset=["pdb"])
bp = pd.read_csv(f"{ROOT}/benchmark/pairs_924.csv")
bench = sorted({(r.query.lower(), r.q_chain) for r in bp.itertuples()} | {(r.target.lower(), r.t_chain) for r in bp.itertuples()})
bseq = {k: cache[(f"{ROOT}/structures/{k[0]}.pdb", k[1])][1] for k in bench if (f"{ROOT}/structures/{k[0]}.pdb", k[1]) in cache}
val = {(r.pdb, r.chain): cache[(f"{ROOT}/structures_validation/{r.pdb}.pdb", r.chain)][1] for r in ent.itertuples() if (f"{ROOT}/structures_validation/{r.pdb}.pdb", r.chain) in cache}
DB = sum(len(s) for s in bseq.values()); LAM, K = 0.267, 0.041
al = PairwiseAligner(); al.substitution_matrix = substitution_matrices.load("BLOSUM62"); al.open_gap_score = -11; al.extend_gap_score = -1; al.mode = "local"
def kmers(s, k=3): return {s[i:i + k] for i in range(len(s) - k + 1)}
bk = {k: kmers(v) for k, v in bseq.items()}
C = f"{HERE}/blast_cache.json"
try: lc = json.load(open(C))
except Exception: lc = {}
t0 = time.time()
for key, s in val.items():
    sk = f"{key[0]}|{key[1]}"
    if sk in lc: continue
    if time.time() - t0 > 150: break
    ks = kmers(s); best = [0.0, None]
    for bkey, bs in bseq.items():
        if len(ks & bk[bkey]) / max(1, min(len(ks), len(bk[bkey]))) < 0.05: continue
        sc = al.score(s, bs)
        if sc > best[0]: best = [float(sc), list(bkey)]
    E = K * len(val[key]) * DB * math.exp(-LAM * best[0]) if best[0] > 0 else 1e9
    lc[sk] = [best[0], E, best[1]]
json.dump(lc, open(C, "w"))
if len(lc) < len(val): print("done", len(lc), "/", len(val)); sys.exit(0)
drop = {tuple(k.split("|")) for k, v in lc.items() if v[1] < 1e-3}
full = pd.read_csv(f"{HERE}/validation_pairs.csv")
ok = full[[((r.query, r.q_chain) in val and (r.target, r.t_chain) in val and (r.query, r.q_chain) not in drop and (r.target, r.t_chain) not in drop) for r in full.itertuples()]].copy()
ok.to_csv(f"{HERE}/validation_pairs_blast_noTM.csv", index=False)
Es = np.array([v[1] for v in lc.values()])
with open(f"{HERE}/leakage_report_blast.txt", "w") as f:
    f.write(f"database: {len(bseq)} benchmark chains, {DB} residues; validation chains: {len(val)}\n")
    f.write(f"best-hit E-value distribution: <1e-10 {int((Es<1e-10).sum())}, 1e-10..1e-3 {int(((Es>=1e-10)&(Es<1e-3)).sum())}, 1e-3..1 {int(((Es>=1e-3)&(Es<1)).sum())}, >=1 {int((Es>=1).sum())}\n")
    f.write(f"chains dropped (E<1e-3): {len(drop)}; pairs {len(full)} -> {len(ok)}\n")
    for k in sorted(drop, key=lambda x: lc[f'{x[0]}|{x[1]}'][1]): f.write(f"  {k} -> {lc[f'{k[0]}|{k[1]}'][2]} E={lc[f'{k[0]}|{k[1]}'][1]:.1e}\n")
print(open(f"{HERE}/leakage_report_blast.txt").read()[:2500])
