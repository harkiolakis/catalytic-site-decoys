"""Global TM-score (TM-align, normalised by the query chain) for every benchmark pair -> benchmark_tm.csv
python compute_benchmark_tm.py <pairs.csv> <structures_dir> <out.csv>"""
import sys, os, csv
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from tmtools import tm_align
from csd.io import parse_chain_trace
PAIRS, SDIR, OUT = sys.argv[1:4]; cache = {}
def trace(p, c):
    k = (p, c)
    if k not in cache:
        t = parse_chain_trace(os.path.join(SDIR, p + ".pdb"), c); cache[k] = (t.coords, "".join(t.sequence))
    return cache[k]
rows = []
for r in csv.DictReader(open(PAIRS)):
    try:
        (a, sa), (b, sb) = trace(r["query"].lower(), r["q_chain"]), trace(r["target"].lower(), r["t_chain"])
        rows.append(dict(pair_idx=r["pair_idx"], tm=tm_align(a, b, sa, sb).tm_norm_chain1))
    except Exception as e: rows.append(dict(pair_idx=r["pair_idx"], tm=""))
with open(OUT, "w", newline="") as f: w = csv.DictWriter(f, ["pair_idx", "tm"]); w.writeheader(); w.writerows(rows)
print("wrote", len(rows))
