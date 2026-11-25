"""JESS (pyjess) under the same-input protocol as v12.
Template = query catalytic residues, two atoms per residue (CA + functional atom(s)), residue names allowed
within the same chemical class as v12; molecule = ONLY the target's annotated (or decoy) residues.
Subtemplates of size 2..4 (all subsets, capped at 60 per size); best hit = largest matched size, then lowest E-value.
Score = -log10(E) of best hit (0 if no hit). Decoys regenerated with the exact seeds/generators of run_v12.py.
python run_jess.py <sites> <roles> <pairs.csv> <structdir> <out_prefix> [trials]"""
import sys, os, csv, json, itertools, warnings
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pyjess
SITES, ROLES, PAIRS, SDIR, OUTP = sys.argv[1:6]; NT = int(sys.argv[6]) if len(sys.argv) > 6 else 5
_argv = sys.argv; sys.argv = ["x", "a", "b", "c", "d", "e"]; import run_v12 as RV; sys.argv = _argv
import sidechain_matcher as v12
RMSD_T, DIST_C, MAX_DYN, CAP = 2.0, 1.5, 1.5, 60
CLASS_MEMBERS = {}
for res, cl in v12.CLASS.items(): CLASS_MEMBERS.setdefault(cl, []).append(res)
def tatoms(chain, resnum, resname, at):
    cl = v12.CLASS.get(resname); names = CLASS_MEMBERS[cl] if cl else [resname]
    fa = sorted({a for n in names for a in v12.FUNC_ATOMS.get(n, [])}) or None
    out = [pyjess.TemplateAtom(chain_id="A", residue_number=resnum, x=at["CA"][0], y=at["CA"][1], z=at["CA"][2], residue_names=names, atom_names=["CA"])]
    own = [a for a in v12.FUNC_ATOMS.get(resname, []) if a in at and a != "CA"]
    if own and fa:
        p = at[own[0]]
        out.append(pyjess.TemplateAtom(chain_id="A", residue_number=resnum, x=p[0], y=p[1], z=p[2], residue_names=names, atom_names=[a for a in fa if a != "CA"]))
    return out
def read_atoms(pdb, chain):
    res = {}
    for l in open(pdb, errors="ignore"):
        if l.startswith("ENDMDL"): break
        if not l.startswith(("ATOM", "HETATM")) or l[21] != chain or l[26] != " " or l[16] not in " A": continue
        try: r = int(l[22:26])
        except ValueError: continue
        d = res.setdefault(r, {"name": l[17:20], "at": {}, "lines": []})
        d["at"][l[12:16].strip()] = (float(l[30:38]), float(l[38:46]), float(l[46:54])); d["lines"].append(l)
    return res
def molecule(res, resnums):
    lines = ["ATOM  " + l[6:16] + " " + l[17:21] + "A" + l[22:] if l.startswith("HETATM") else l[:21] + "A" + l[22:] for r in resnums if r in res for l in res[r]["lines"]]
    return pyjess.Molecule.loads("\n".join(lines) + "\nEND\n", format="pdb", id="t")
def best_hit(Q, qres, mol, rng):
    best = None
    for k in (4, 3, 2):
        if k > len(qres): continue
        subs = list(itertools.combinations(qres, k)); rng.shuffle(subs); subs = subs[:CAP]
        temps = []
        for s in subs:
            ats = [a for r in s for a in tatoms("A", r, Q[r]["name"], Q[r]["at"])]
            temps.append(pyjess.Template(ats, id="|".join(map(str, s))))
        hits = list(pyjess.Jess(temps).query(mol, RMSD_T, DIST_C, MAX_DYN, ignore_chain="residues"))
        for h in hits:
            sc = -h.log_evalue / np.log(10)
            if best is None or sc > best[0]:
                qs = list(map(int, h.template().id.split("|")))
                tr = []
                for a in h.atoms():
                    if a.residue_number not in tr: tr.append(a.residue_number)
                best = (sc, k, h.rmsd, list(zip(qs, tr)))
        if best is not None: break  # largest matching size wins
    return best
_g = {}
def work(p):
    from csd.io import parse_chain_trace
    s, roles = _g["s"], _g["r"]
    qp, tp, qc, tc = p["query"].lower(), p["target"].lower(), p["q_chain"], p["t_chain"]
    try:
        A = parse_chain_trace(os.path.join(SDIR, qp + ".pdb"), qc); B = parse_chain_trace(os.path.join(SDIR, tp + ".pdb"), tc)
        RA = v12.parse_residues(os.path.join(SDIR, qp + ".pdb"), qc); RB = v12.parse_residues(os.path.join(SDIR, tp + ".pdb"), tc)
        Q = read_atoms(os.path.join(SDIR, qp + ".pdb"), qc); T = read_atoms(os.path.join(SDIR, tp + ".pdb"), tc)
    except Exception: return [], []
    mA = {A.res_ids[i][1]: i for i in range(len(A.res_ids))}; mB = {B.res_ids[i][1]: i for i in range(len(B.res_ids))}
    qa = sorted({mA[r] for r in s.get(qp, {}).get(qc, []) if r in mA and r in RA}); ta = sorted({mB[r] for r in s.get(tp, {}).get(tc, []) if r in mB and r in RB})
    if len(qa) < 2 or len(ta) < 2: return [], []
    qres = [A.res_ids[i][1] for i in qa if A.res_ids[i][1] in Q and "CA" in Q[A.res_ids[i][1]]["at"]]
    rng = np.random.default_rng(int(p["pair_idx"]))
    def run(tidx):
        tres = [B.res_ids[i][1] for i in tidx]
        return best_hit(Q, qres, molecule(T, tres), rng)
    b = run(ta); row = dict(pair_idx=p["pair_idx"], cath_relation=p.get("cath_relation", ""), m=len(qa), n=len(ta))
    if b:
        sh = [len(set(roles.get(f"{qp}|{qc}|{i}", [])) & set(roles.get(f"{tp}|{tc}|{j}", []))) > 0 for i, j in b[3] if roles.get(f"{qp}|{qc}|{i}") and roles.get(f"{tp}|{tc}|{j}")]
        row.update(jess_score=b[0], jess_k=b[1], jess_rmsd=b[2], jess_role=float(np.mean(sh)) if sh else np.nan,
                   jess_ordered=all(b[3][x][1] < b[3][x + 1][1] for x in range(len(b[3]) - 1)))
    else: row.update(jess_score=0.0, jess_k=0)
    decs = []; seq = "".join(B.sequence)
    for mode in RV.DMODES:
        for trial in range(NT):
            r2 = np.random.default_rng(int(p["pair_idx"]) * 1000 + trial * 10 + RV.DMODES.index(mode))
            d = RV.make_decoy(mode, B.coords, seq, ta, r2)
            dr = dict(pair_idx=p["pair_idx"], mode=mode, trial=trial, m=len(qa), n=len(ta))
            if d is None or any(B.res_ids[i][1] not in RB for i in d): dr["error"] = "no_decoy"
            else:
                bd = run(d); dr.update(jess_score=bd[0] if bd else 0.0, jess_k=bd[1] if bd else 0)
            decs.append(dr)
    return [row], decs
def init(): _g["s"] = json.load(open(SITES)); _g["r"] = json.load(open(ROLES))
if __name__ == "__main__":
    from multiprocessing import Pool
    pairs = list(csv.DictReader(open(PAIRS))); R, D = [], []
    with Pool(2, initializer=init) as pool:
        for r, d in pool.imap_unordered(work, pairs, chunksize=2): R += r; D += d
    for rows, suf in ((R, "_real.csv"), (D, "_decoys.csv")):
        keys = list(dict.fromkeys(k for r in rows for k in r))
        with open(OUTP + suf, "w", newline="") as f: w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(rows)
    print("wrote", len(R), len(D))
