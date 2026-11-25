"""Run the octant-descriptor aligner: v10 (exact-offset) and v11 (connectivity-aware dual-channel) on a pair list.
python run_modec.py <sites.json> <pairs.csv> <structures_dir> <out.csv> [workers]"""
import sys, os, csv, json, warnings
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from multiprocessing import Pool
from csd.io import parse_chain_trace
SITES, PAIRS, SDIR, OUT = sys.argv[1:5]; W = int(sys.argv[5]) if len(sys.argv) > 5 else 2
_g = {}
def init():
    from dp_fast import warmup; warmup()
    _g["sites"] = json.load(open(SITES))
def idx_set(trace, lst):
    m = {trace.res_ids[i][1]: i for i in range(len(trace.res_ids))}
    return sorted({m[r] for r in lst if r in m})
def work(p):
    from octant_aligner import align_octant_mode_c
    out = []
    s = _g["sites"]
    ql = s.get(p["query"].lower(), {}).get(p["q_chain"], []); tl = s.get(p["target"].lower(), {}).get(p["t_chain"], [])
    base = dict(p, q_listed=len(ql), t_listed=len(tl))
    try:
        A = parse_chain_trace(os.path.join(SDIR, p["query"].lower() + ".pdb"), p["q_chain"])
        Bt = parse_chain_trace(os.path.join(SDIR, p["target"].lower() + ".pdb"), p["t_chain"])
    except Exception as e:
        return [dict(base, method="parse", error=str(e)[:80])]
    qa, ta = idx_set(A, ql), idx_set(Bt, tl)
    base.update(q_cat_total=len(qa), t_cat_total=len(ta))
    if len(qa) < 2 or len(ta) < 2:
        return [dict(base, method="skip", error="too_few_cat")]
    for name, kw in [("v10", dict(use_connectivity=False)), ("v11", dict(use_connectivity=True, w_conn=1.0))]:
        try:
            pairs, rmsd, n, tm, info = align_octant_mode_c(A, Bt, set(qa), set(ta), **kw)
            out.append(dict(base, method=name, n_aligned=n, rmsd=rmsd, tm_like=tm,
                            corr=";".join(f"{A.res_ids[a][1]}-{Bt.res_ids[b][1]}" for a, b in pairs),
                            error=info.get("error", "") if isinstance(info, dict) else ""))
        except Exception as e:
            out.append(dict(base, method=name, error=str(e)[:80]))
    return out
if __name__ == "__main__":
    pairs = list(csv.DictReader(open(PAIRS)))
    rows = []
    with Pool(W, initializer=init) as pool:
        for i, r in enumerate(pool.imap_unordered(work, pairs, chunksize=4)):
            rows += r
            if i % 100 == 0: print(i, flush=True)
    keys = sorted({k for r in rows for k in r})
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(rows)
    print("wrote", OUT, len(rows))
