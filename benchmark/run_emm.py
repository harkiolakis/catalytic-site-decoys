"""Enzyme Motif Miner (Hackett et al.) under the same-input protocol.
For pair (q, t): EMM templates of q's M-CSA entry (excluding any built from t's PDB) are searched, with EMM defaults
(size-specific Jess parameters, best_match, logistic-regression filter on RMSD + residue orientation, templates with
>=3 side-chain residues), against a molecule containing only t's annotated residues (or decoy residues, same generators/seeds).
Score = 10*[predicted correct] + matched residues - RMSD/10 for the best match (0 if none).
python run_emm.py <sites> <roles> <pairs.csv> <structdir> <out_prefix> [trials]"""
import sys, os, csv, json, warnings, collections
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
SITES, ROLES, PAIRS, SDIR, OUTP = sys.argv[1:6]; NT = int(sys.argv[6]) if len(sys.argv) > 6 else 5
_argv = sys.argv; sys.argv = ["x", "a", "b", "c", "d", "e"]; import run_v12 as RV; sys.argv = _argv
import sidechain_matcher as v12, run_jess as RJ
import enzymm
PERMISSIVE = os.environ.get("EMM_PERMISSIVE") == "1"
from enzymm.template import load_templates
from enzymm.jess_run import Matcher
D = os.path.join(os.path.dirname(enzymm.__file__), "data")
_g = {}
def init():
    _g["s"] = json.load(open(SITES)); _g["r"] = json.load(open(ROLES))
    H = json.load(open(os.path.join(D, "catalytic_residue_homologs_information.json")))
    _g["hroles"] = {(int(m), pc): {int(k): set(v.get("roles_summary", [])) for k, v in d["residues"].items()} for m, chains in H.items() for pc, d in chains.items()}
    by = collections.defaultdict(list)
    for t in load_templates(warn=False, verbose=False): by[t.mcsa_id].append(t)
    _g["T"] = by
    pm = {}
    for r in csv.reader(open(os.path.join(HERE, "..", "mcsa", "curated_data.csv"))):
        if r[0] != "M-CSA ID": pm[r[2].lower()] = int(r[0][1:])
    _g["pm"] = pm
def best(templates, mol):
    if not templates: return None
    import rich.console
    M = Matcher(templates, cpus=1, filter_matches=False, match_small_templates=PERMISSIVE, verbose=False, warn=False, console=rich.console.Console(file=open(os.devnull, 'w'), width=100, force_terminal=False))
    try: ms = M.run_single(mol)
    except ZeroDivisionError: return None
    b = None
    for m in ms:
        try: ok = bool(m.predicted_correct)
        except Exception: ok = False
        n = len(m.atom_triplets); sc = 10 * ok + n - m.hit.rmsd / 10
        if b is None or sc > b[0]: b = (sc, n, m.hit.rmsd, ok, m)
    return b
def work(p):
    from csd.io import parse_chain_trace
    s, roles = _g["s"], _g["r"]
    qp, tp, qc, tc = p["query"].lower(), p["target"].lower(), p["q_chain"], p["t_chain"]
    try:
        A = parse_chain_trace(os.path.join(SDIR, qp + ".pdb"), qc); B = parse_chain_trace(os.path.join(SDIR, tp + ".pdb"), tc)
        RA = v12.parse_residues(os.path.join(SDIR, qp + ".pdb"), qc); RB = v12.parse_residues(os.path.join(SDIR, tp + ".pdb"), tc)
        T = RJ.read_atoms(os.path.join(SDIR, tp + ".pdb"), tc)
    except Exception: return [], []
    mA = {A.res_ids[i][1]: i for i in range(len(A.res_ids))}; mB = {B.res_ids[i][1]: i for i in range(len(B.res_ids))}
    qa = sorted({mA[r] for r in s.get(qp, {}).get(qc, []) if r in mA and r in RA}); ta = sorted({mB[r] for r in s.get(tp, {}).get(tc, []) if r in mB and r in RB})
    if len(qa) < 2 or len(ta) < 2: return [], []
    mq = _g["pm"].get(qp); temps = [t for t in _g["T"].get(mq, []) if (t.pdb_id or "").lower() != tp and (PERMISSIVE or t.effective_size >= 3)]  # EMM default: >=3 side-chain residues; permissive: all sizes
    row = dict(pair_idx=p["pair_idx"], m=len(qa), n=len(ta), q_mcsa=mq, t_mcsa=_g["pm"].get(tp), n_templates=len(temps),
               max_template_size=max((t.effective_size for t in temps), default=0))
    def run(tidx): return best(temps, RJ.molecule(T, [B.res_ids[i][1] for i in tidx]))
    b = run(ta)
    if b:
        m = b[4]; tm = m.hit.template(); tatoms = list(tm)
        tres = [(tatoms[3 * i].chain_id, tatoms[3 * i].residue_number) for i in range(len(m.atom_triplets))]
        qres = [int(x[2]) for x in m.matched_residues]
        sh = []
        for (ch, rn), tr in zip(tres, qres):
            a = _g["hroles"].get((mq, f"{(tm.pdb_id or '').lower()}{ch}"), {}).get(rn); bb = roles.get(f"{tp}|{tc}|{tr}")
            if a and bb: sh.append(len(a & set(bb)) > 0)
        row.update(emm_score=b[0], emm_k=b[1], emm_rmsd=b[2], emm_ok=b[3], emm_role=float(np.mean(sh)) if sh else np.nan)
    else: row.update(emm_score=0.0, emm_k=0, emm_ok=False)
    decs = []; seq = "".join(B.sequence)
    for mode in RV.DMODES:
        for trial in range(NT):
            r2 = np.random.default_rng(int(p["pair_idx"]) * 1000 + trial * 10 + RV.DMODES.index(mode))
            d = RV.make_decoy(mode, B.coords, seq, ta, r2)
            dr = dict(pair_idx=p["pair_idx"], mode=mode, trial=trial)
            if d is None or any(B.res_ids[i][1] not in RB for i in d): dr["error"] = "no_decoy"
            else:
                bd = run(d); dr.update(emm_score=bd[0] if bd else 0.0, emm_k=bd[1] if bd else 0, emm_ok=bd[3] if bd else False)
            decs.append(dr)
    return [row], decs
if __name__ == "__main__":
    from multiprocessing import Pool
    pairs = list(csv.DictReader(open(PAIRS))); R, Dd = [], []
    with Pool(2, initializer=init) as pool:
        for r, d in pool.imap_unordered(work, pairs, chunksize=2): R += r; Dd += d
    for rows, suf in ((R, "_real.csv"), (Dd, "_decoys.csv")):
        keys = list(dict.fromkeys(k for r in rows for k in r))
        with open(OUTP + suf, "w", newline="") as f: w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(rows)
    print("wrote", len(R), len(Dd))
