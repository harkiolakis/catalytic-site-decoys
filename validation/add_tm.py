"""Add global TM-scores (TM-align, normalised by the query chain) to a validation pair list.
python add_tm.py <in.csv> <out.csv>   (needs trace_cache.pkl from leakage_and_strata.py)"""
import pandas as pd, pickle, os, sys
from tmtools import tm_align
HERE = os.path.dirname(os.path.abspath(__file__)); c = pickle.load(open(f"{HERE}/trace_cache.pkl", "rb"))
inp, out = sys.argv[1], sys.argv[2]
p = pd.read_csv(inp, dtype=str)
def g(pdb, ch): return c[(f"{HERE}/../structures_validation/{pdb}.pdb", ch)]
tm = []
for r in p.itertuples():
    (a, sa), (b, sb) = g(r.query, r.q_chain), g(r.target, r.t_chain); tm.append(tm_align(a, b, sa, sb).tm_norm_chain1)
p["tm"] = tm; p["tm_stratum"] = pd.cut(p.tm, [-1, 0.3, 0.5, 2], labels=["TM<0.3", "TM0.3-0.5", "TM>=0.5"])
p.to_csv(out, index=False); print(len(p)); print(p.groupby(["ec_stratum", "tm_stratum"], observed=True).size().to_string())
