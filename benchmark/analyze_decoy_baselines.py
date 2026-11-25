"""Is the compactness confound method-independent?  AUC (real vs decoy) for every method on identical decoys.
Lower RMSD = more site-like for RMSD methods; higher score for v12/JESS/EnzyMM.  AUC = P(real ranks above decoy).
Pair-level bootstrap 95% CIs (1000 resamples).  Writes DECOY_BASELINES.md and T_decoy_baselines.csv."""
import numpy as np, pandas as pd, os
from scipy import stats
rng = np.random.default_rng(1)
H = os.path.dirname(os.path.abspath(__file__)); V = os.path.join(H, "..", "validation")
out = []; L = []
def say(s=""): L.append(s)
def auc(real, dec, higher):
    r, d = np.asarray(real, float), np.asarray(dec, float)
    if not higher: r, d = -r, -d
    return stats.mannwhitneyu(r, d).statistic / (len(r) * len(d))
def boot(R, D, col_r, col_d, higher):
    """R: per-pair real values (index pair_idx); D: decoy rows with pair_idx."""
    R = R.dropna(); D = D[D.pair_idx.isin(R.index)].dropna(subset=[col_d]); R = R[R.index.isin(D.pair_idx)]
    if len(R) < 20: return None
    dec = {p: g[col_d].values for p, g in D.groupby("pair_idx")}; ps = R.index.values
    a = auc(R.values, D[col_d].values, higher); v = []
    for _ in range(1000):
        s = rng.choice(ps, len(ps)); v.append(auc(R.loc[s].values, np.concatenate([dec[q] for q in s]), higher))
    lo, hi = np.percentile(v, [2.5, 97.5]); return a, lo, hi, len(R)
def fmt(x): return "—" if x is None else f"{x[0]:.3f} ({x[1]:.3f}–{x[2]:.3f})"

def block(label, B, extra, modes, keep=None):
    """B: baseline csv (real rows mode=='real'); extra: list of (name, real_series, decoy_df, col, higher)."""
    B = B.copy(); B["pair_idx"] = B.pair_idx.astype(int)
    if keep is not None: B = B[B.pair_idx.isin(keep)]
    real = B[B["mode"] == "real"].set_index("pair_idx"); dec = B[B["mode"] != "real"]
    meths = [("Random mapping (median)", "random_any", False), ("Best order-preserving mapping", "oracle_mono", False),
             ("Best unrestricted mapping", "oracle_any", False), ("Distance-only aligner", "dist_only", False)]
    if "v11" in B: meths.insert(0, ("v11 (octant descriptor)", "v11", False))
    rows = []
    for name, col, hi in meths:
        rows.append(dict(method=name, **{m: fmt(boot(real[col], dec[dec["mode"] == m], col, col, hi)) for m in modes}))
    for name, R, Dd, col, hi in extra:
        R = R.copy(); R.index = R.index.astype(int); Dd = Dd.copy(); Dd["pair_idx"] = Dd.pair_idx.astype(int)
        if keep is not None: R = R[R.index.isin(keep)]; Dd = Dd[Dd.pair_idx.isin(keep)]
        rows.append(dict(method=name, **{m: fmt(boot(R, Dd[Dd["mode"] == m], col, col, hi)) for m in modes}))
    t = pd.DataFrame(rows); t.insert(0, "set", label); out.append(t)
    say(f"### {label}"); say(t.drop(columns="set").to_markdown(index=False)); say()
    # Rg of decoys vs real (target site)
    return t

MODES = ["random", "surface", "composition", "compact", "spatial"]
say("# Is the compactness confound method-independent?\n")
say("AUC = probability that the real catalytic site ranks above a decoy (0.5 = chance). Decoys replace the target's catalytic residues; "
    "5 per pair per family. RMSD-based methods: lower RMSD ranks higher. 95% CI: pair-level bootstrap. "
    "Mapping baselines use k = min(m, n) residues; the random-mapping value is the median RMSD of 200 random injective mappings.\n")

# --- benchmark, v11 decoy generator (the decoys of the v11 decoy table) ---
B1 = pd.read_csv(os.path.join(H, "baselines", "bench.csv"))
DM = pd.read_csv(os.path.join(H, "decoys_mcsa_current.csv")); DM = DM[DM.real_n > 0]
v11r = DM.groupby("pair_idx").real_rmsd.first(); v11d = DM.rename(columns={"decoy_rmsd": "v"})[["pair_idx", "mode", "v"]]
block("Benchmark (691 pairs) — decoy set 1 (generator used for v11)", B1, [("v11 (octant descriptor)", v11r, v11d, "v", False)], MODES)

# --- benchmark, v12 decoy generator (decoys of v12 / JESS / EnzyMM tables) ---
if os.path.exists(os.path.join(H, "baselines", "bench_g12.csv")):
    B2 = pd.read_csv(os.path.join(H, "baselines", "bench_g12.csv"))
    VR = pd.read_csv(os.path.join(H, "v12_real.csv")).set_index("pair_idx"); VD = pd.read_csv(os.path.join(H, "v12_decoys.csv"))
    JR = pd.read_csv(os.path.join(H, "jess_real.csv")).set_index("pair_idx"); JD = pd.read_csv(os.path.join(H, "jess_decoys.csv"))
    ER = pd.read_csv(os.path.join(H, "emm_real.csv")).set_index("pair_idx"); ED = pd.read_csv(os.path.join(H, "emm_decoys.csv"))
    block("Benchmark — decoy set 2 (generator used for v12, JESS and EnzyMM)", B2,
          [("v12 (two-point, order-free)", VR.v12_score, VD, "v12_score", True), ("JESS (residue-restricted)", JR.jess_score, JD, "jess_score", True),
           ("EnzyMM (residue-restricted)", ER.emm_score, ED, "emm_score", True)], MODES)

# --- validation (BLAST-filtered set) ---
U = pd.read_csv(os.path.join(V, "validation_pairs_union.csv")); keep = set(U[U.set_blast == 1].pair_idx.astype(int))
for f, lab, gen in (("val.csv", "decoy set 1", 1), ("val_g12.csv", "decoy set 2", 2)):
    p = os.path.join(H, "baselines", f)
    if not os.path.exists(p): continue
    extra = []
    if gen == 2:
        VR = pd.read_csv(os.path.join(V, "val_v12_real.csv")).set_index("pair_idx"); VD = pd.read_csv(os.path.join(V, "val_v12_decoys.csv"))
        JR = pd.read_csv(os.path.join(V, "val_jess_real.csv")).set_index("pair_idx"); JD = pd.read_csv(os.path.join(V, "val_jess_decoys.csv"))
        extra = [("v12 (two-point, order-free)", VR.v12_score, VD, "v12_score", True), ("JESS (residue-restricted)", JR.jess_score, JD, "jess_score", True)]
    block(f"Validation (BLAST-filtered set) — {lab}", pd.read_csv(p), extra, MODES, keep)

say("Note: spatial decoys reproduce every pairwise Cα distance of the real site (to 1–2 Å), so any method that uses only Cα geometry "
    "is at chance against them by construction; compact (Rg-matched) decoys are the informative control.")
T = pd.concat(out); T.to_csv(os.path.join(H, "T_decoy_baselines.csv"), index=False)
open(os.path.join(H, "DECOY_BASELINES.md"), "w").write("\n".join(L)); print("\n".join(L))
