"""Post-download steps fixed in PREREGISTRATION.md: (1) drop validation chains with >30% identity to any benchmark
chain; (2) TM-score strata for validation pairs. Writes validation_pairs_final.csv and leakage_report.txt."""
import os, sys, csv, json
import numpy as np, pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices
from tmtools import tm_align
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
sys.path.insert(0, os.path.join(ROOT, "benchmark"))
from csd.io import parse_chain_trace
import pickle
CACHE = os.path.join(HERE, "trace_cache.pkl")
_cache = pickle.load(open(CACHE, "rb")) if os.path.exists(CACHE) else {}
def trace(path, ch):
    k = (path, ch)
    if k not in _cache:
        t = parse_chain_trace(path, ch); _cache[k] = (np.asarray(t.coords, float), "".join(t.sequence))
    return _cache[k]
STAGE = sys.argv[1] if len(sys.argv) > 1 else "all"
ent = pd.read_csv(f"{HERE}/validation_entries.csv"); pairs = pd.read_csv(f"{HERE}/validation_pairs.csv")
bp = pd.read_csv(f"{ROOT}/benchmark/pairs_924.csv")
bench = {(r.query.lower(), r.q_chain) for r in bp.itertuples()} | {(r.target.lower(), r.t_chain) for r in bp.itertuples()}
val, missing = {}, []
for r in ent.itertuples():
    p = f"{ROOT}/structures_validation/{r.pdb}.pdb"
    try: val[(r.pdb, r.chain)] = trace(p, r.chain)
    except Exception as e: missing.append((r.pdb, r.chain, str(e)[:40]))
bseq = {}
for pdb, ch in sorted(bench):
    try: bseq[(pdb, ch)] = trace(f"{ROOT}/structures/{pdb}.pdb", ch)[1]
    except Exception: pass
if STAGE == "parse":
    pickle.dump(_cache, open(CACHE, "wb")); print("cached", len(_cache)); sys.exit(0)
al = PairwiseAligner(); al.substitution_matrix = substitution_matrices.load("BLOSUM62"); al.open_gap_score = -10; al.extend_gap_score = -0.5; al.mode = "local"
def ident(a, b):
    """Local alignment; identity over aligned columns (incl. internal gaps); 0 if coverage of the shorter chain < 50%."""
    aln = al.align(a, b)[0]; blocks = list(zip(*aln.aligned))
    if not blocks: return 0.0
    same = sum(1 for (s1, e1), (s2, e2) in blocks for i in range(e1 - s1) if a[s1 + i] == b[s2 + i])
    ncol = max(blocks[-1][0][1] - blocks[0][0][0], blocks[-1][1][1] - blocks[0][1][0])
    cov = min(blocks[-1][0][1] - blocks[0][0][0], blocks[-1][1][1] - blocks[0][1][0]) / min(len(a), len(b))
    return same / ncol if cov >= 0.5 else 0.0
def kmers(s, k=3): return {s[i:i + k] for i in range(len(s) - k + 1)}
bk = {k: kmers(v) for k, v in bseq.items()}
import time
LC = os.path.join(HERE, "leak_cache.json")
lc = json.load(open(LC)) if os.path.exists(LC) else {}
t0 = time.time()
for key, (_, s) in val.items():
    sk = f"{key[0]}|{key[1]}"
    if sk in lc: continue
    if time.time() - t0 > 140: break
    ks = kmers(s); best = 0.0; hit = None
    for bkey, bs in bseq.items():
        if len(ks & bk[bkey]) / max(1, min(len(ks), len(bk[bkey]))) < 0.10: continue
        x = ident(s, bs)
        if x > best: best, hit = x, bkey
    lc[sk] = [best, list(hit) if hit else None]
json.dump(lc, open(LC, "w"))
if len(lc) < len(val): print("identity done for", len(lc), "/", len(val)); sys.exit(0)
leak = {tuple(k.split("|")): (v[0], v[1]) for k, v in lc.items()}
drop = {k for k, (x, _) in leak.items() if x > 0.30}
keep = pairs[[((r.query, r.q_chain) in val and (r.target, r.t_chain) in val and (r.query, r.q_chain) not in drop and (r.target, r.t_chain) not in drop) for r in pairs.itertuples()]].copy()
tms = []
for r in keep.itertuples():
    (ca, sa), (cb, sb) = val[(r.query, r.q_chain)], val[(r.target, r.t_chain)]
    tms.append(tm_align(ca, cb, sa, sb).tm_norm_chain1)
keep["tm"] = tms; keep["tm_stratum"] = pd.cut(keep.tm, [-1, 0.3, 0.5, 2], labels=["TM<0.3", "TM0.3-0.5", "TM>=0.5"])
keep.to_csv(f"{HERE}/validation_pairs_final.csv", index=False)
with open(f"{HERE}/leakage_report.txt", "w") as f:
    f.write(f"validation chains parsed: {len(val)} / {len(ent)}; unparsed: {missing}\n")
    f.write(f"chains dropped for >30% identity to a benchmark chain: {len(drop)}\n")
    for k in sorted(drop): f.write(f"  {k} -> {leak[k][1]} identity {leak[k][0]:.2f}\n")
    f.write(f"pairs: {len(pairs)} -> {len(keep)}\n{keep.groupby(['ec_stratum','tm_stratum'], observed=True).size().to_string()}\n")
print(open(f"{HERE}/leakage_report.txt").read())
