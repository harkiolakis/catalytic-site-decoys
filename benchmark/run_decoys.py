"""Decoy controls for the octant-descriptor aligner (v11, full configuration).
Modes: random, surface (as original), spatial (all pairwise Cα distances within tol of the catalytic set),
sequence (same sequence spacing ±2), spatial_seq (both; spacing ±5), composition (same residue types).
python run_decoys.py <sites.json> <pairs.csv> <structures_dir> <out.csv> [trials] [workers]"""
import sys, os, csv, json, random, warnings
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from multiprocessing import Pool
SITES, PAIRS, SDIR, OUT = sys.argv[1:5]
NT = int(sys.argv[5]) if len(sys.argv) > 5 else 5; W = int(sys.argv[6]) if len(sys.argv) > 6 else 2
MODES = os.environ.get("DECOY_MODES", "random,surface,spatial,sequence,spatial_seq,composition").split(",")
ALL_MODES = ["random", "surface", "spatial", "sequence", "spatial_seq", "composition", "compact"]
_g = {}
def init():
    from dp_fast import warmup; warmup(); _g["sites"] = json.load(open(SITES))

def backtrack(Dall, T, avail, tol, rng, seq=None, max_nodes=20000):
    """Find decoy list matching pairwise distances of T (and optionally sequence spacing ±seq)."""
    DT = Dall[np.ix_(T, T)]; n = len(T); nodes = [0]
    def rec(chosen):
        if len(chosen) == n: return list(chosen)
        j = len(chosen); mask = avail.copy(); mask[chosen] = False
        for k, c in enumerate(chosen):
            mask &= np.abs(Dall[c] - DT[k, j]) <= tol
        if seq is not None:
            idx = np.arange(len(avail)); want = chosen[0] + (T[j] - T[0])
            mask &= np.abs(idx - want) <= seq
        cand = np.flatnonzero(mask); rng.shuffle(cand)
        for c in cand[:20]:
            nodes[0] += 1
            if nodes[0] > max_nodes: return None
            r = rec(chosen + [int(c)])
            if r: return r
        return None
    starts = np.flatnonzero(avail); rng.shuffle(starts)
    for s in starts[:200]:
        r = rec([int(s)])
        if r: return r
        if nodes[0] > max_nodes: break
    return None

def make_decoy(mode, B, T, Dall, nbr, rng):
    nres = len(B.coords); avail = np.ones(nres, bool); avail[T] = False; avail[:2] = False  # descriptor needs i>=2
    n = len(T)
    if mode == "random":
        c = np.flatnonzero(avail); return sorted(rng.choice(c, n, replace=False).tolist()) if len(c) >= n else None, 0
    if mode == "surface":
        c = np.flatnonzero(avail); c = c[np.argsort(nbr[c])][:max(n, len(c) // 2)]
        return sorted(rng.choice(c, n, replace=False).tolist()) if len(c) >= n else None, 0
    if mode == "sequence":
        span = T[-1] - T[0]
        for _ in range(500):
            s0 = int(rng.integers(2, max(3, nres - span - 2)))
            d = [s0] + [s0 + (t - T[0]) + int(rng.integers(-2, 3)) for t in T[1:]]
            if all(2 <= x < nres and avail[x] for x in d) and len(set(d)) == n: return sorted(d), 0
        return None, 0
    if mode == "composition":
        seq = B.sequence; out = []
        for t in T:
            c = [i for i in np.flatnonzero(avail) if seq[i] == seq[t] and i not in out]
            if not c: return None, 0
            out.append(int(rng.choice(c)))
        return sorted(out), 0
    if mode == "compact":
        Y = B.coords[T]; rg0 = float(np.sqrt(((Y - Y.mean(0)) ** 2).sum(1).mean()))
        c = np.flatnonzero(avail)
        for _ in range(300):
            ctr = int(rng.choice(c)); near = c[np.argsort(Dall[ctr, c])][:3 * n]
            d = rng.choice(near, n, replace=False); Z = B.coords[d]
            rg = float(np.sqrt(((Z - Z.mean(0)) ** 2).sum(1).mean()))
            if abs(rg - rg0) <= 0.15 * rg0 + 0.5: return sorted(d.tolist()), 0
        return None, 0
    if mode in ("spatial", "spatial_seq"):
        for tol in (1.0, 1.5, 2.0, 3.0):
            r = backtrack(Dall, list(T), avail, tol, rng, seq=5 if mode == "spatial_seq" else None)
            if r: return sorted(r), tol
        return None, np.nan

def work(p):
    from csd.io import parse_chain_trace
    from octant_aligner import align_octant_mode_c
    s = _g["sites"]; qp, tp = p["query"].lower(), p["target"].lower()
    ql, tl = s.get(qp, {}).get(p["q_chain"], []), s.get(tp, {}).get(p["t_chain"], [])
    try:
        A = parse_chain_trace(os.path.join(SDIR, qp + ".pdb"), p["q_chain"]); B = parse_chain_trace(os.path.join(SDIR, tp + ".pdb"), p["t_chain"])
    except Exception: return []
    mA = {A.res_ids[i][1]: i for i in range(len(A.res_ids))}; mB = {B.res_ids[i][1]: i for i in range(len(B.res_ids))}
    qa = sorted({mA[r] for r in ql if r in mA}); ta = sorted({mB[r] for r in tl if r in mB})
    if len(qa) < 2 or len(ta) < 2: return []
    X = B.coords; Dall = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1))
    nbr = (Dall < 10).sum(1) - 1
    def rg(ix): Y = X[ix]; return float(np.sqrt(((Y - Y.mean(0)) ** 2).sum(1).mean()))
    pr, rr, nr, tm, info = align_octant_mode_c(A, B, set(qa), set(ta), use_connectivity=True, w_conn=1.0)
    base = dict(pair_idx=p["pair_idx"], query=p["query"], target=p["target"], cath_relation=p["cath_relation"],
                m=len(qa), n=len(ta), real_n=nr, real_rmsd=rr if nr else np.nan, real_rg=rg(ta))
    out = []
    for mode in MODES:
        for trial in range(NT):
            rng = np.random.default_rng(int(p["pair_idx"]) * 1000 + trial * 10 + ALL_MODES.index(mode))
            d, tol = make_decoy(mode, B, ta, Dall, nbr, rng)
            row = dict(base, mode=mode, trial=trial)
            if d is None: row["error"] = "no_decoy"; out.append(row); continue
            try:
                pairs, r, nal, tm, info = align_octant_mode_c(A, B, set(qa), set(d), use_connectivity=True, w_conn=1.0)
                row.update(decoy_n=nal, decoy_rmsd=r if nal else np.nan, decoy_rg=rg(d), tol=tol,
                           decoy_nbr=float(nbr[d].mean()), real_nbr=float(nbr[ta].mean()))
            except Exception as e: row["error"] = str(e)[:50]
            out.append(row)
    return out

if __name__ == "__main__":
    pairs = list(csv.DictReader(open(PAIRS))); rows = []
    with Pool(W, initializer=init) as pool:
        for i, r in enumerate(pool.imap_unordered(work, pairs, chunksize=2)):
            rows += r
            if i % 50 == 0: print(i, flush=True)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(rows)
    print("wrote", OUT, len(rows))
