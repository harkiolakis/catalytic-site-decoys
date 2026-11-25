"""Benchmark v12 on real pairs + decoys.  python run_v12.py <sites> <roles> <pairs.csv> <structdir> <out_prefix> [trials] [workers]"""
import sys, os, csv, json, warnings
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from multiprocessing import Pool
SITES, ROLES, PAIRS, SDIR, OUTP = sys.argv[1:6]; NT = int(sys.argv[6]) if len(sys.argv) > 6 else 5; W = int(sys.argv[7]) if len(sys.argv) > 7 else 2
DMODES = ["random", "surface", "compact", "spatial", "composition"]
VARIANTS = {"v12": dict(), "v12_desc": dict(w_desc=0.5), "v12_nochem": dict(nochem=True), "v12_ca_only": dict(ca_only=True)}
_g = {}
def init():
    _g["sites"] = json.load(open(SITES)); _g["roles"] = json.load(open(ROLES))

def desc_sets(trace, idxs):
    from tourist.descriptor import compute_signature_connectivity_aware, relative_offsets_multi, CONNECTIVITY_MATCH_SLOTS
    out = []
    for i in idxs:
        sig = compute_signature_connectivity_aware(trace.coords, i, 12.0, 4.0, 8.0) if i >= 2 else None
        rel = relative_offsets_multi(sig, i) if sig is not None else None
        out.append({s for s in CONNECTIVITY_MATCH_SLOTS if rel is not None and s < len(rel) and len(rel[s]) > 0})
    return out

def jacc(a, b): return np.array([[len(x & y) / len(x | y) if (x | y) else 0.0 for y in b] for x in a])

def run_variant(q, t, dq, dt, name):
    import v12
    p = dict(v12.PARAMS); kw = VARIANTS[name]
    if kw.get("ca_only"): q = [(a, c, c) for a, c, f in q]; t = [(a, c, c) for a, c, f in t]
    if kw.get("nochem"): q = [("ALA", c, f) for a, c, f in q]; t = [("ALA", c, f) for a, c, f in t]
    if "w_desc" in kw: p["w_desc"] = kw["w_desc"]
    return v12.align(q, t, desc_sim=jacc(dq, dt) if p["w_desc"] > 0 else None, p=p)

def make_decoy(mode, coords, seq, T, rng):
    nres = len(coords); avail = np.ones(nres, bool); avail[T] = False; avail[:2] = False; n = len(T)
    Dall = np.sqrt(((coords[:, None] - coords[None]) ** 2).sum(-1)); nbr = (Dall < 10).sum(1) - 1
    c = np.flatnonzero(avail)
    if len(c) < n: return None
    if mode == "random": return sorted(rng.choice(c, n, replace=False).tolist())
    if mode == "surface":
        c2 = c[np.argsort(nbr[c])][:max(n, len(c) // 2)]; return sorted(rng.choice(c2, n, replace=False).tolist())
    if mode == "composition":
        out = []
        for t in T:
            cc = [i for i in c if seq[i] == seq[t] and i not in out]
            if not cc: return None
            out.append(int(rng.choice(cc)))
        return sorted(out)
    Y = coords[T]; rg0 = float(np.sqrt(((Y - Y.mean(0)) ** 2).sum(1).mean()))
    if mode == "compact":
        for _ in range(300):
            ctr = int(rng.choice(c)); near = c[np.argsort(Dall[ctr, c])][:3 * n]; d = rng.choice(near, n, replace=False)
            Z = coords[d]
            if abs(float(np.sqrt(((Z - Z.mean(0)) ** 2).sum(1).mean())) - rg0) <= 0.15 * rg0 + 0.5: return sorted(d.tolist())
        return None
    if mode == "spatial":  # all pairwise Ca distances within tol of the catalytic set (backtracking)
        DT = Dall[np.ix_(T, T)]
        for tol in (1.0, 1.5, 2.0):
            nodes = [0]
            def rec(ch):
                if len(ch) == n: return ch
                mask = avail.copy(); mask[ch] = False
                for k, x in enumerate(ch): mask &= np.abs(Dall[x] - DT[k, len(ch)]) <= tol
                cand = np.flatnonzero(mask); rng.shuffle(cand)
                for x in cand[:20]:
                    nodes[0] += 1
                    if nodes[0] > 20000: return None
                    r = rec(ch + [int(x)])
                    if r: return r
                return None
            st = c.copy(); rng.shuffle(st)
            for s0 in st[:200]:
                r = rec([int(s0)])
                if r: return sorted(r)
                if nodes[0] > 20000: break
        return None

def work(p):
    import v12
    from tourist.io import parse_chain_trace
    s, roles = _g["sites"], _g["roles"]
    qp, tp, qc, tc = p["query"].lower(), p["target"].lower(), p["q_chain"], p["t_chain"]
    try:
        A = parse_chain_trace(os.path.join(SDIR, qp + ".pdb"), qc); B = parse_chain_trace(os.path.join(SDIR, tp + ".pdb"), tc)
        RA = v12.parse_residues(os.path.join(SDIR, qp + ".pdb"), qc); RB = v12.parse_residues(os.path.join(SDIR, tp + ".pdb"), tc)
    except Exception as e: return [], []
    mA = {A.res_ids[i][1]: i for i in range(len(A.res_ids))}; mB = {B.res_ids[i][1]: i for i in range(len(B.res_ids))}
    qa = sorted({mA[r] for r in s.get(qp, {}).get(qc, []) if r in mA and r in RA}); ta = sorted({mB[r] for r in s.get(tp, {}).get(tc, []) if r in mB and r in RB})
    if len(qa) < 2 or len(ta) < 2: return [], []
    qres = [A.res_ids[i][1] for i in qa]; q = [RA[r] for r in qres]; dq = desc_sets(A, qa)
    def target(idx):
        rr = [B.res_ids[i][1] for i in idx]; return rr, [RB[r] for r in rr], desc_sets(B, idx)
    tres, t, dt = target(ta)
    real = dict(p, m=len(q), n=len(t))
    for v in VARIANTS:
        r = run_variant(q, t, dq, dt, v)
        if r is None: real[f"{v}_k"] = 0; real[f"{v}_score"] = 0.0; continue
        mt = r["matches"]
        real.update({f"{v}_k": len(mt), f"{v}_score": r["score_norm"], f"{v}_rmsd_ca": r["rmsd_ca"], f"{v}_rmsd_fc": r["rmsd_fc"],
                     f"{v}_ordered": all(mt[x][1] < mt[x + 1][1] for x in range(len(mt) - 1))})
        sh = [len(set(roles.get(f"{qp}|{qc}|{qres[i]}", [])) & set(roles.get(f"{tp}|{tc}|{tres[j]}", []))) > 0
              for i, j in mt if roles.get(f"{qp}|{qc}|{qres[i]}") and roles.get(f"{tp}|{tc}|{tres[j]}")]
        real[f"{v}_role"] = float(np.mean(sh)) if sh else np.nan
        if v == "v12": real["v12_corr"] = ";".join(f"{qres[i]}-{tres[j]}" for i, j in mt)
    decs = []
    seq = "".join(B.sequence)
    for mode in DMODES:
        for trial in range(NT):
            rng = np.random.default_rng(int(p["pair_idx"]) * 1000 + trial * 10 + DMODES.index(mode))
            d = make_decoy(mode, B.coords, seq, ta, rng)
            row = dict(pair_idx=p["pair_idx"], cath_relation=p["cath_relation"], mode=mode, trial=trial, m=len(q), n=len(t))
            if d is None or any(B.res_ids[i][1] not in RB for i in d): row["error"] = "no_decoy"; decs.append(row); continue
            _, td, ddt = target(d)
            for v in ("v12", "v12_desc"):
                r = run_variant(q, td, dq, ddt, v)
                row[f"{v}_k"] = len(r["matches"]) if r else 0; row[f"{v}_score"] = r["score_norm"] if r else 0.0
                row[f"{v}_rmsd_fc"] = r["rmsd_fc"] if r else np.nan
            decs.append(row)
    return [real], decs

if __name__ == "__main__":
    pairs = list(csv.DictReader(open(PAIRS))); R, Dd = [], []
    with Pool(W, initializer=init) as pool:
        for i, (r, d) in enumerate(pool.imap_unordered(work, pairs, chunksize=2)):
            R += r; Dd += d
    for rows, suf in ((R, "_real.csv"), (Dd, "_decoys.csv")):
        keys = list(dict.fromkeys(k for r in rows for k in r))
        with open(OUTP + suf, "w", newline="") as f:
            w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(rows)
    print("wrote", len(R), len(Dd))
