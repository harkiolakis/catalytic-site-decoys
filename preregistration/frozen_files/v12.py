"""CS-TOURIST v12: permutation-tolerant, chemistry-aware catalytic-site matching.

Each catalytic residue is represented by two points: its C-alpha and the centroid of its
side-chain functional group. Alignment = best rigid superposition found from two-residue
seeds, followed by optimal (Hungarian) one-to-one assignment and iterative refinement.
No sequence-order constraint. Parameters fixed a priori (see PARAMS).
"""
import numpy as np
from scipy.optimize import linear_sum_assignment

PARAMS = dict(d_cut=3.0, d0=2.0, seed_tol=2.0, max_seeds=3000, n_refine=3, w_desc=0.0, min_chem=0.2)

FUNC_ATOMS = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"], "HIS": ["ND1", "NE2"], "LYS": ["NZ"], "ARG": ["NE", "NH1", "NH2"],
              "SER": ["OG"], "THR": ["OG1"], "CYS": ["SG"], "TYR": ["OH"], "ASN": ["OD1", "ND2"], "GLN": ["OE1", "NE2"],
              "TRP": ["NE1"], "MET": ["SD"], "GLY": ["CA"]}
CLASS = {"ASP": "acid", "GLU": "acid", "LYS": "base", "ARG": "base", "HIS": "his", "SER": "oh", "THR": "oh", "TYR": "tyr",
         "CYS": "sh", "ASN": "amide", "GLN": "amide", "PHE": "arom", "TRP": "arom", "GLY": "gly"}
RELATED = {frozenset(["acid", "amide"]): 0.5, frozenset(["his", "base"]): 0.5, frozenset(["oh", "sh"]): 0.6, frozenset(["oh", "tyr"]): 0.6,
           frozenset(["tyr", "arom"]): 0.5, frozenset(["his", "acid"]): 0.3, frozenset(["his", "amide"]): 0.3}

def chem(a, b):
    ca, cb = CLASS.get(a, "hyd"), CLASS.get(b, "hyd")
    if ca == cb: return 1.0
    return RELATED.get(frozenset([ca, cb]), PARAMS["min_chem"])

def parse_residues(pdb, chain):
    """{resid: (resname, CA xyz, functional-centroid xyz)} for one chain, first model, altloc blank/A, no insertion codes."""
    atoms = {}
    for l in open(pdb, errors="ignore"):
        if l.startswith("ENDMDL"): break
        if not (l.startswith("ATOM") or l.startswith("HETATM")) or l[21] != chain or l[26] != " " or l[16] not in " A": continue
        try: r = int(l[22:26])
        except ValueError: continue
        d = atoms.setdefault(r, {"name": l[17:20], "at": {}})
        d["at"][l[12:16].strip()] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    out = {}
    for r, d in atoms.items():
        at = d["at"]
        if "CA" not in at: continue
        fa = [at[a] for a in FUNC_ATOMS.get(d["name"], []) if a in at]
        if not fa:
            sc = [v for k, v in at.items() if k not in ("N", "CA", "C", "O", "OXT") and not k.startswith("H")]
            fa = sc if sc else [at["CA"]]
        out[r] = (d["name"], at["CA"], np.mean(fa, 0))
    return out

def kabsch(P, Q):
    """R,t minimising |R Q + t - P|."""
    pc, qc = P.mean(0), Q.mean(0); H = (Q - qc).T @ (P - pc)
    U, S, Vt = np.linalg.svd(H); d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, pc - R @ qc

def align(q, t, desc_sim=None, p=PARAMS):
    """q,t: lists of (resname, ca, fc). Returns dict with matches [(i,j)], score, rmsd_ca, rmsd_fc."""
    m, n = len(q), len(t)
    if m < 2 or n < 2: return None
    QA = np.array([x[1] for x in q]); QF = np.array([x[2] for x in q]); TA = np.array([x[1] for x in t]); TF = np.array([x[2] for x in t])
    C = np.array([[chem(a[0], b[0]) for b in t] for a in q])
    if desc_sim is not None and p["w_desc"] > 0: C = C + p["w_desc"] * desc_sim
    dQA = np.linalg.norm(QA[:, None] - QA[None], axis=-1); dTA = np.linalg.norm(TA[:, None] - TA[None], axis=-1)
    dQF = np.linalg.norm(QF[:, None] - QF[None], axis=-1); dTF = np.linalg.norm(TF[:, None] - TF[None], axis=-1)
    seeds = []
    for i in range(m):
        for i2 in range(i + 1, m):
            for j in range(n):
                for j2 in range(n):
                    if j2 == j: continue
                    if abs(dQA[i, i2] - dTA[j, j2]) > p["seed_tol"] or abs(dQF[i, i2] - dTF[j, j2]) > p["seed_tol"]: continue
                    seeds.append((C[i, j] + C[i2, j2], i, i2, j, j2))
    if not seeds: return None
    seeds.sort(reverse=True); seeds = seeds[:p["max_seeds"]]
    best = None
    for _, i, i2, j, j2 in seeds:
        R, tt = kabsch(np.vstack([QA[[i, i2]], QF[[i, i2]]]), np.vstack([TA[[j, j2]], TF[[j, j2]]]))
        for _ in range(p["n_refine"]):
            A2 = TA @ R.T + tt; F2 = TF @ R.T + tt
            dA = np.linalg.norm(QA[:, None] - A2[None], axis=-1); dF = np.linalg.norm(QF[:, None] - F2[None], axis=-1)
            W = C * 0.5 * (1 / (1 + (dA / p["d0"]) ** 2) + 1 / (1 + (dF / p["d0"]) ** 2))
            W[(dA > p["d_cut"]) | (dF > p["d_cut"])] = 0
            r, c = linear_sum_assignment(-W); keep = W[r, c] > 0; r, c = r[keep], c[keep]
            if len(r) < 2: break
            R, tt = kabsch(np.vstack([QA[r], QF[r]]), np.vstack([TA[c], TF[c]]))
        if len(r) < 2: continue
        score = float(W[r, c].sum())
        if best is None or score > best["score"] + 1e-9:
            A2 = TA[c] @ R.T + tt; F2 = TF[c] @ R.T + tt
            best = dict(matches=list(zip(r.tolist(), c.tolist())), score=score,
                        rmsd_ca=float(np.sqrt(((QA[r] - A2) ** 2).sum(1).mean())), rmsd_fc=float(np.sqrt(((QF[r] - F2) ** 2).sum(1).mean())))
    if best: best["score_norm"] = best["score"] / max(m, n)
    return best
