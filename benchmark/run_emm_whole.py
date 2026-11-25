"""Enzyme Motif Miner in its native whole-structure mode (default Jess parameters, templates with >=3 side-chain
residues, logistic-regression filter recorded as predicted_correct).
(A) pair task: templates of the query's M-CSA entry (none built from the target PDB) searched against the whole target
    structure; negative control: templates of 2 random M-CSA entries whose EC class (1st digit) differs from both
    query and target, searched against the same structure.
(B) annotation task: the whole template library minus the target's own M-CSA entry (and any template built from the
    target PDB) searched against the whole target structure.
A match is 'on site' if any matched target residue is an annotated catalytic residue (residue number, any chain);
precision = fraction of matched residues that are annotated.
python run_emm_whole.py <sites.json> <pairs.csv> <structdir> <out_prefix> [workers]"""
import sys, os, csv, json, collections, warnings
warnings.filterwarnings("ignore")
import numpy as np
SITES, PAIRS, SDIR, OUTP = sys.argv[1:5]; W = int(sys.argv[5]) if len(sys.argv) > 5 else 2
HERE = os.path.dirname(os.path.abspath(__file__))
_g = {}
def init():
    from enzymm.template import load_templates
    _g["s"] = json.load(open(SITES)); T = list(load_templates(warn=False, verbose=False))
    by = collections.defaultdict(list)
    for t in T: by[t.mcsa_id].append(t)
    _g["T"], _g["by"] = T, by
    pm, ec = {}, {}
    for r in csv.reader(open(os.path.join(HERE, "..", "mcsa", "curated_data.csv"))):
        if r[0] == "M-CSA ID": continue
        pm.setdefault(r[2].lower(), int(r[0][1:])); ec.setdefault(int(r[0][1:]), r[3])
    _g["pm"], _g["ec"] = pm, ec

def matcher(temps):
    import rich.console
    from enzymm.jess_run import Matcher
    return Matcher(temps, cpus=1, filter_matches=False, verbose=False, warn=False,
                   console=rich.console.Console(file=open(os.devnull, "w"), width=100))

def summarise(ms, site, prefix):
    """ms: list of matches. Returns dict of outcome measures for best match (EMM ranking) among all / predicted-correct."""
    out = {}
    def key(m): return (bool(getattr(m, "predicted_correct", False)), len(m.matched_residues), -m.hit.rmsd)
    for tag, sub in (("any", ms), ("pc", [m for m in ms if getattr(m, "predicted_correct", False)])):
        out[f"{prefix}_{tag}_n"] = len(sub)
        if not sub: out[f"{prefix}_{tag}_hit"] = False; continue
        b = max(sub, key=key); res = [int(x[2]) for x in b.matched_residues]
        on = [r in site for r in res]
        out.update({f"{prefix}_{tag}_hit": True, f"{prefix}_{tag}_k": len(res), f"{prefix}_{tag}_rmsd": b.hit.rmsd,
                    f"{prefix}_{tag}_onsite": any(on), f"{prefix}_{tag}_prec": float(np.mean(on)),
                    f"{prefix}_{tag}_recall": len(set(res) & site) / len(site) if site else np.nan,
                    f"{prefix}_{tag}_mcsa": getattr(b.hit.template(), "mcsa_id", None)})
        # any predicted-correct/any match on site at all
        out[f"{prefix}_{tag}_anyonsite"] = any(any(int(x[2]) in site for x in m.matched_residues) for m in sub)
    return out

def ec_shared(a, b):
    k = 0
    for x, y in zip((a or "").split("."), (b or "").split(".")):
        if x == y and x not in ("-", "", "n"): k += 1
        else: break
    return k

def work(p):
    import pyjess
    s, pm, by = _g["s"], _g["pm"], _g["by"]
    qp, tp = p["query"].lower(), p["target"].lower()
    site = {int(r) for ch, rs in s.get(tp, {}).items() for r in rs}
    if not site: return None
    try: mol = pyjess.Molecule.load(os.path.join(SDIR, tp + ".pdb"))
    except Exception as e: return dict(pair_idx=p["pair_idx"], error=str(e)[:40])
    mq = int(p["q_mcsa"][1:]) if p.get("q_mcsa") else pm.get(qp); mt = int(p["t_mcsa"][1:]) if p.get("t_mcsa") else pm.get(tp)
    row = dict(pair_idx=p["pair_idx"], query=qp, target=tp, q_mcsa=mq, t_mcsa=mt, site_n=len(site),
               ec_q=_g["ec"].get(mq), ec_t=_g["ec"].get(mt))
    ok = lambda t: (t.pdb_id or "").lower() != tp and t.effective_size >= 3
    tq = [t for t in by.get(mq, []) if ok(t)] if mq != mt else []
    row["A_ntemp"] = len(tq)
    if tq: row.update(summarise(matcher(tq).run_single(mol), site, "A"))
    # negative control: 2 entries with different EC class from both
    rng = np.random.default_rng(int(p["pair_idx"]) + 99)
    c1q, c1t = (row["ec_q"] or "x").split(".")[0], (row["ec_t"] or "x").split(".")[0]
    cands = sorted(e for e in by if e not in (mq, mt) and (_g["ec"].get(e) or "y").split(".")[0] not in (c1q, c1t)
                   and any(ok(t) for t in by[e]))
    for j, e in enumerate(rng.choice(cands, 2, replace=False)):
        row.update(summarise(matcher([t for t in by[e] if ok(t)]).run_single(mol), site, f"C{j}"))
    return row

def work_lib(tp_mt):
    """(B) one target structure against the whole library minus its own entry."""
    import pyjess
    tp, mt = tp_mt; s = _g["s"]; site = {int(r) for ch, rs in s.get(tp, {}).items() for r in rs}
    if not site: return None
    try: mol = pyjess.Molecule.load(os.path.join(SDIR, tp + ".pdb"))
    except Exception as e: return dict(target=tp, error=str(e)[:40])
    temps = [t for t in _g["T"] if t.mcsa_id != mt and (t.pdb_id or "").lower() != tp and t.effective_size >= 3]
    ms = matcher(temps).run_single(mol)
    row = dict(target=tp, t_mcsa=mt, ec_t=_g["ec"].get(mt), site_n=len(site)); row.update(summarise(ms, site, "B"))
    for tag in ("any", "pc"):
        e = row.get(f"B_{tag}_mcsa")
        row[f"B_{tag}_ec_shared"] = ec_shared(_g["ec"].get(e), row["ec_t"]) if e else np.nan
    return row

if __name__ == "__main__":
    from multiprocessing import Pool
    pairs = list(csv.DictReader(open(PAIRS), ))
    with Pool(W, initializer=init) as pool:
        R = [r for r in pool.imap_unordered(work, pairs, chunksize=4) if r]
        init(); tg = sorted({(p["target"].lower(), int(p["t_mcsa"][1:]) if p.get("t_mcsa") else _g["pm"].get(p["target"].lower())) for p in pairs})
        L = [r for r in pool.imap_unordered(work_lib, tg, chunksize=2) if r]
    for rows, suf in ((R, "_pairs.csv"), (L, "_library.csv")):
        keys = list(dict.fromkeys(k for r in rows for k in r))
        with open(OUTP + suf, "w", newline="") as f: w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(rows)
    print("wrote", len(R), len(L))
