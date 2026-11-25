"""Descriptor-free baselines on the same decoys as run_decoys.py (identical seeds).
For the real target site and each decoy: best unrestricted mapping RMSD (oracle_any), best order-preserving
mapping RMSD (oracle_mono), median RMSD of random mappings (random_any), all at k=min(m,n); and the
distance-only aligner (octant-descriptor aligner machinery with descriptor and BLOSUM weights set to 0).
python run_decoy_baselines.py <sites.json> <pairs.csv> <structures_dir> <out.csv> [trials] [workers] [start] [end]"""
import sys, os, csv, json, warnings
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
SITES, PAIRS, SDIR, OUT = sys.argv[1:5]
NT = int(sys.argv[5]) if len(sys.argv) > 5 else 5; W = int(sys.argv[6]) if len(sys.argv) > 6 else 2
S0 = int(sys.argv[7]) if len(sys.argv) > 7 else 0; S1 = int(sys.argv[8]) if len(sys.argv) > 8 else 10**9
_argv = sys.argv; sys.argv = [_argv[0], SITES, PAIRS, SDIR, "/dev/null"]
import run_decoys as RD           # reuse make_decoy and seeding exactly
sys.argv = [_argv[0], SITES, "x", PAIRS, SDIR, "/dev/null"]
import run_experiments as RE      # batch_rmsd, mappings
sys.argv = ["x", "a", "b", "c", "d", "e"]
import run_v12 as RV              # the decoy generator used for v12 / JESS / EnzyMM
sys.argv = _argv
import numpy as np
from multiprocessing import Pool
GEN = os.environ.get("GEN", "v11")   # v11: run_decoys generator; v12: run_v12 generator
MODES = os.environ.get("DECOY_MODES", "random,surface,composition,compact,spatial").split(",")
RE.CAP = 20000
_g = {}
def init():
    from dp_fast import warmup; warmup(); _g["sites"] = json.load(open(SITES))

def measures(A, B, XA, qa, idx, rng):
    from octant_aligner_ablate import align_octant_mode_c
    XB = B.coords[idx]; m, n = len(qa), len(idx); k = min(m, n); out = {}
    for mono in (True, False):
        Q, T, _ = RE.mappings(m, n, k, mono, rng); r = RE.batch_rmsd(XA[Q], XB[T])
        out["oracle_mono" if mono else "oracle_any"] = float(r.min())
        if not mono:
            sel = rng.choice(len(r), min(200, len(r)), replace=False); out["random_any"] = float(np.median(r[sel]))
    try:
        pairs, rmsd, nal, tm, info = align_octant_mode_c(A, B, set(qa), set(idx), use_connectivity=True, w_conn=1.0, w_desc=0.0, beta=0.0)
        out["dist_only"] = rmsd if nal else np.nan; out["dist_only_n"] = nal
    except Exception: out["dist_only"] = np.nan; out["dist_only_n"] = 0
    if GEN == "v12":
        try:
            pairs, rmsd, nal, tm, info = align_octant_mode_c(A, B, set(qa), set(idx), use_connectivity=True, w_conn=1.0)
            out["v11"] = rmsd if nal else np.nan; out["v11_n"] = nal
        except Exception: out["v11"] = np.nan; out["v11_n"] = 0
    out["k"] = k; return out

def work(p):
    from csd.io import parse_chain_trace
    s = _g["sites"]; qp, tp = p["query"].lower(), p["target"].lower()
    ql, tl = s.get(qp, {}).get(p["q_chain"], []), s.get(tp, {}).get(p["t_chain"], [])
    try:
        A = parse_chain_trace(os.path.join(SDIR, qp + ".pdb"), p["q_chain"]); B = parse_chain_trace(os.path.join(SDIR, tp + ".pdb"), p["t_chain"])
    except Exception: return []
    mA = {A.res_ids[i][1]: i for i in range(len(A.res_ids))}; mB = {B.res_ids[i][1]: i for i in range(len(B.res_ids))}
    qa = sorted({mA[r] for r in ql if r in mA}); ta = sorted({mB[r] for r in tl if r in mB})
    if len(qa) < 2 or len(ta) < 2: return []
    X = B.coords; Dall = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1)); nbr = (Dall < 10).sum(1) - 1
    XA = A.coords[qa]; prng = np.random.default_rng(int(p["pair_idx"]) + 7)
    base = dict(pair_idx=p["pair_idx"], cath_relation=p["cath_relation"], m=len(qa), n=len(ta))
    out = [dict(base, mode="real", trial=-1, **measures(A, B, XA, qa, ta, prng))]
    for mode in MODES:
        for trial in range(NT):
            if GEN == "v12":
                rng = np.random.default_rng(int(p["pair_idx"]) * 1000 + trial * 10 + RV.DMODES.index(mode))
                d = RV.make_decoy(mode, B.coords, "".join(B.sequence), ta, rng)
            else:
                rng = np.random.default_rng(int(p["pair_idx"]) * 1000 + trial * 10 + RD.ALL_MODES.index(mode))
                d, tol = RD.make_decoy(mode, B, ta, Dall, nbr, rng)
            row = dict(base, mode=mode, trial=trial)
            if d is None: row["error"] = "no_decoy"
            else: row.update(measures(A, B, XA, qa, d, prng))
            out.append(row)
    return out

if __name__ == "__main__":
    pairs = list(csv.DictReader(open(PAIRS)))[S0:S1]; rows = []
    with Pool(W, initializer=init) as pool:
        for i, r in enumerate(pool.imap_unordered(work, pairs, chunksize=1)):
            rows += r
            if i % 25 == 0: print(i, flush=True)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(rows)
    print("wrote", OUT, len(rows))
