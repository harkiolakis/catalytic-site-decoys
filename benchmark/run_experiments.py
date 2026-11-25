"""Same-input baselines, ablations, side-chain RMSD and role concordance for the octant-descriptor aligner (v11).
python run_experiments.py <sites.json> <roles.json> <pairs.csv> <structures_dir> <out.csv> [workers]"""
import sys, os, csv, json, math, itertools, warnings
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from multiprocessing import Pool
SITES, ROLES, PAIRS, SDIR, OUT = sys.argv[1:6]; W = int(sys.argv[6]) if len(sys.argv) > 6 else 2
CAP = 100000; NRAND = 200
CONFIGS = {
    "v11_full":        dict(use_connectivity=True, w_conn=1.0),
    "v11_conn_only":   dict(use_connectivity=True, w_conn=0.0),
    "v10_exact":       dict(use_connectivity=False),
    "v11_no_blosum":   dict(use_connectivity=True, w_conn=1.0, beta=0.0),
    "v11_no_dist":     dict(use_connectivity=True, w_conn=1.0, dist_constraint_weight=0.0),
    "v11_no_refine":   dict(use_connectivity=True, w_conn=1.0, refined=False),
    "blosum_dist_only":dict(use_connectivity=True, w_conn=1.0, w_desc=0.0),
    "dist_only":       dict(use_connectivity=True, w_conn=1.0, w_desc=0.0, beta=0.0),
    "blosum_only":     dict(use_connectivity=True, w_conn=1.0, w_desc=0.0, dist_constraint_weight=0.0),
}
_g = {}
def init():
    from dp_fast import warmup; warmup()
    _g["sites"] = json.load(open(SITES)); _g["roles"] = json.load(open(ROLES))
    _g["rng"] = np.random.default_rng(os.getpid())

def batch_rmsd(A, B):
    """A,B: (N,k,3) paired sets -> RMSD after optimal superposition (B onto A)."""
    a = A - A.mean(1, keepdims=True); b = B - B.mean(1, keepdims=True)
    H = np.einsum("nki,nkj->nij", b, a)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(np.einsum("nij,njk->nik", np.transpose(Vt, (0, 2, 1)), np.transpose(U, (0, 2, 1)))))
    S[:, 2] *= d
    e = (a * a).sum((1, 2)) + (b * b).sum((1, 2)) - 2 * S.sum(1)
    return np.sqrt(np.maximum(e, 0) / A.shape[1])

def kabsch_rmsd(X, Y): return float(batch_rmsd(X[None], Y[None])[0])

def mappings(m, n, k, mono, rng):
    """Index arrays (N,k) for query and target of all (or CAP sampled) injective k-mappings."""
    if mono: total = math.comb(m, k) * math.comb(n, k)
    else: total = math.comb(m, k) * math.perm(n, k)
    if total <= CAP:
        qs = list(itertools.combinations(range(m), k))
        ts = list(itertools.combinations(range(n), k)) if mono else list(itertools.permutations(range(n), k))
        Q = np.repeat(np.array(qs), len(ts), 0); T = np.tile(np.array(ts), (len(qs), 1))
        return Q, T, True
    Q = np.array([np.sort(rng.choice(m, k, replace=False)) for _ in range(CAP)])
    T = np.array([np.sort(rng.choice(n, k, replace=False)) if mono else rng.choice(n, k, replace=False) for _ in range(CAP)])
    return Q, T, False

def sidechain_centroids(pdb, chain, resids):
    out = {}
    for l in open(pdb, errors="ignore"):
        if l.startswith("ENDMDL"): break
        if not l.startswith("ATOM") or l[21] != chain: continue
        try: r = int(l[22:26])
        except ValueError: continue
        if r not in resids or l[26] != " ": continue
        name = l[12:16].strip(); alt = l[16]
        if alt not in " A": continue
        xyz = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
        d = out.setdefault(r, {"sc": [], "ca": None})
        if name == "CA": d["ca"] = xyz
        elif name not in ("N", "C", "O", "OXT") and not name.startswith("H"): d["sc"].append(xyz)
    return {r: np.mean(v["sc"], 0) if v["sc"] else np.array(v["ca"]) for r, v in out.items() if v["ca"] is not None}

def role_share(roles, qp, qc, qr, tp, tc, tr):
    a = set(roles.get(f"{qp}|{qc}|{qr}", [])); b = set(roles.get(f"{tp}|{tc}|{tr}", []))
    if not a or not b: return np.nan
    return float(len(a & b) > 0)

def work(p):
    from csd.io import parse_chain_trace
    from octant_aligner_ablate import align_octant_mode_c
    s, roles, rng = _g["sites"], _g["roles"], _g["rng"]
    qp, tp, qc, tc = p["query"].lower(), p["target"].lower(), p["q_chain"], p["t_chain"]
    ql, tl = s.get(qp, {}).get(qc, []), s.get(tp, {}).get(tc, [])
    row = dict(p)
    try:
        A = parse_chain_trace(os.path.join(SDIR, qp + ".pdb"), qc); Bt = parse_chain_trace(os.path.join(SDIR, tp + ".pdb"), tc)
    except Exception as e:
        row["error"] = "parse:" + str(e)[:60]; return row
    mA = {A.res_ids[i][1]: i for i in range(len(A.res_ids))}; mB = {Bt.res_ids[i][1]: i for i in range(len(Bt.res_ids))}
    qa = sorted({mA[r] for r in ql if r in mA}); ta = sorted({mB[r] for r in tl if r in mB})
    row.update(m=len(qa), n=len(ta))
    if len(qa) < 2 or len(ta) < 2:
        row["error"] = "too_few_cat"; return row
    XA, XB = A.coords[qa], Bt.coords[ta]
    row["rg_q"] = float(np.sqrt(((XA - XA.mean(0)) ** 2).sum(1).mean())); row["rg_t"] = float(np.sqrt(((XB - XB.mean(0)) ** 2).sum(1).mean()))
    resA = [A.res_ids[i][1] for i in qa]; resB = [Bt.res_ids[i][1] for i in ta]
    posA = {v: k for k, v in enumerate(qa)}; posB = {v: k for k, v in enumerate(ta)}
    full_corr = None
    for name, kw in CONFIGS.items():
        try:
            pairs, rmsd, nal, tm, info = align_octant_mode_c(A, Bt, set(qa), set(ta), **kw)
        except Exception as e:
            row[f"{name}_err"] = str(e)[:40]; continue
        row[f"{name}_n"] = nal; row[f"{name}_rmsd"] = rmsd if nal else np.nan
        if nal and name in ("v11_full", "v10_exact", "blosum_only", "dist_only"):
            sh = [role_share(roles, qp, qc, A.res_ids[a][1], tp, tc, Bt.res_ids[b][1]) for a, b in pairs]
            row[f"{name}_role"] = np.nanmean(sh) if not all(np.isnan(sh)) else np.nan
        if name == "v11_full" and nal: full_corr = [(posA[a], posB[b]) for a, b in pairs]
    if full_corr is None:
        row["error"] = "v11_no_alignment"; return row
    k = len(full_corr); row["k"] = k
    row["v11_corr"] = ";".join(f"{resA[i]}-{resB[j]}" for i, j in full_corr)
    row["v11_monotone"] = all(full_corr[x][1] < full_corr[x + 1][1] for x in range(k - 1))
    # side-chain centroid RMSD under v11 correspondence
    scA = sidechain_centroids(os.path.join(SDIR, qp + ".pdb"), qc, set(resA)); scB = sidechain_centroids(os.path.join(SDIR, tp + ".pdb"), tc, set(resB))
    try:
        row["v11_sc_rmsd"] = kabsch_rmsd(np.array([scA[resA[i]] for i, j in full_corr]), np.array([scB[resB[j]] for i, j in full_corr]))
    except KeyError: row["v11_sc_rmsd"] = np.nan
    m, n = len(qa), len(ta)
    for mono in (True, False):
        tag = "mono" if mono else "any"
        Q, T, exact = mappings(m, n, k, mono, rng)
        r = batch_rmsd(XA[Q], XB[T])
        row[f"oracle_{tag}_rmsd"] = float(r.min()); row[f"oracle_{tag}_exact"] = exact
        best = np.argmin(r)
        sel = rng.choice(len(r), min(NRAND, len(r)), replace=False)
        row[f"random_{tag}_median"] = float(np.median(r[sel])); row[f"random_{tag}_mean"] = float(r[sel].mean())
        row[f"v11_pct_vs_{tag}"] = float((r < row["v11_full_rmsd"] - 1e-9).mean())  # fraction of mappings beating v11
        if not mono:
            sh = [role_share(roles, qp, qc, resA[i], tp, tc, resB[j]) for i, j in zip(Q[best], T[best])]
            row["oracle_any_role"] = np.nanmean(sh) if not all(np.isnan(sh)) else np.nan
            row["oracle_any_monotone"] = bool(all(T[best][x] < T[best][x + 1] for x in range(k - 1)))
            rs = []
            for x in sel[:50]:
                sh = [role_share(roles, qp, qc, resA[i], tp, tc, resB[j]) for i, j in zip(Q[x], T[x])]
                if not all(np.isnan(sh)): rs.append(np.nanmean(sh))
            row["random_any_role"] = float(np.mean(rs)) if rs else np.nan
    return row

if __name__ == "__main__":
    pairs = list(csv.DictReader(open(PAIRS)))
    rows = []
    with Pool(W, initializer=init) as pool:
        for i, r in enumerate(pool.imap_unordered(work, pairs, chunksize=2)):
            rows.append(r)
            if i % 50 == 0: print(i, flush=True)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(rows)
    print("wrote", OUT, len(rows))
