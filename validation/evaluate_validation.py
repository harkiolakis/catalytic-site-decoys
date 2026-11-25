"""Evaluate the pre-registered endpoints (P1-P3) and secondary analyses on the validation set; checks the frozen
fingerprints first. Writes VALIDATION_RESULTS.md.  python evaluate_validation.py"""
import json, hashlib, os
import numpy as np, pandas as pd
from scipy import stats
H = os.path.dirname(os.path.abspath(__file__)); R0 = os.path.join(H, "..")
cfg = json.load(open(os.path.join(R0, "preregistration", "FROZEN_CONFIG.json"))); T = cfg["significance_threshold"]
L = []
def say(s=""): print(s); L.append(str(s))
def md(df): return df.to_markdown(index=False)
say("# Independent validation of Side-chain matcher (v12) — results against pre-registered endpoints")
for f, h in cfg["sha256"].items():
    ok = hashlib.sha256(open(os.path.join(R0, "preregistration", "frozen_files", os.path.basename(f)), "rb").read()).hexdigest() == h
    say(f"- hash check {os.path.basename(f)} (preregistration/frozen_files): {'unchanged' if ok else 'CHANGED'}")
R = pd.read_csv(f"{H}/val_v12_real.csv", dtype={"query": str, "target": str}); D = pd.read_csv(f"{H}/val_v12_decoys.csv")
E = pd.read_csv(f"{H}/val_exp.csv", dtype={"query": str, "target": str}); U = pd.read_csv(f"{H}/validation_pairs_union.csv", dtype=str)
U["pair_idx"] = U.pair_idx.astype(int); U["tm"] = U.tm.astype(float)
M = R.merge(E[["pair_idx", "v11_full_rmsd", "v11_full_role", "random_any_role", "oracle_any_role", "v11_full_n"]], on="pair_idx", how="left")
M = M.merge(U[["pair_idx", "set_prereg", "set_blast", "tm", "tm_stratum", "ec_stratum"]], on="pair_idx", how="left", suffixes=("", "_u"))
rng = np.random.default_rng(11); B = 1000
def auc(r, d): return stats.mannwhitneyu(r, d).statistic / (len(r) * len(d))
def auc_ci(real, g):
    g = g[g.v12_score.notna()]; dec = {p: x.v12_score.values for p, x in g.groupby("pair_idx")}; ps = np.array([p for p in real.index if p in dec]); real = real.loc[ps]; bs = []
    for _ in range(B):
        s = rng.choice(ps, len(ps)); bs.append(auc(real.loc[s].values, np.concatenate([dec[q] for q in s])))
    return auc(real.values, np.concatenate([dec[q] for q in ps])), *np.percentile(bs, [2.5, 97.5]), len(ps)
def mci(x):
    x = np.asarray(x, float); x = x[~np.isnan(x)]
    if len(x) < 2: return "n/a"
    i = rng.integers(0, len(x), (B, len(x))); lo, hi = np.percentile(x[i].mean(1), [2.5, 97.5]); return f"{x.mean():.3f} ({lo:.3f}–{hi:.3f})"
for setname, col, note in [("BLAST E-value leakage filter (declared deviation; recommended primary)", "set_blast", ""),
                           ("Pre-registered %-identity leakage filter (as registered)", "set_prereg", "")]:
    S = M[M[col].astype(str) == "1"].copy(); Dd = D[D.pair_idx.isin(S.pair_idx)]
    real = S.set_index("pair_idx").v12_score
    say(f"\n## {setname}\nPairs analysable: {len(S)} (v12 alignment returned for {100*(S.v12_k>=2).mean():.1f}%)")
    say("\n### Primary endpoints")
    rows = []
    for mode in ["compact", "spatial"]:
        a, lo, hi, n = auc_ci(real, Dd[Dd["mode"] == mode]); rows.append({"endpoint": f"P1 AUC vs {mode} decoys", "pairs": n, "value (95% CI)": f"{a:.3f} ({lo:.3f}–{hi:.3f})", "criterion": "lower CI > 0.60", "met": "YES" if lo > 0.60 else "NO"})
    sig = S[S.v12_score >= T].dropna(subset=["v12_role", "random_any_role"])
    if len(sig) >= 5:
        w = stats.wilcoxon(sig.v12_role, sig.random_any_role); d = sig.v12_role.mean() - sig.random_any_role.mean()
        rows.append({"endpoint": "P2 role concordance of significant matches vs random mapping", "pairs": len(sig),
                     "value (95% CI)": f"{mci(sig.v12_role)} vs {mci(sig.random_any_role)}; Δ={d:.3f}; p={w.pvalue:.2g}", "criterion": "p<0.05 and Δ≥0.05",
                     "met": "YES" if (w.pvalue < 0.05 and d >= 0.05) else "NO"})
    for mode in ["compact", "spatial"]:
        x = Dd[(Dd["mode"] == mode)].v12_score.dropna(); fp = (x >= T).mean()
        rows.append({"endpoint": f"P3 false-positive rate, {mode} decoys", "pairs": Dd[Dd['mode']==mode].pair_idx.nunique(), "value (95% CI)": f"{100*fp:.1f}%", "criterion": "≤10%", "met": "YES" if fp <= 0.10 else "NO"})
    say(md(pd.DataFrame(rows)))
    say(f"Significant matches (score ≥ {T}): {(S.v12_score>=T).sum()} of {len(S)} pairs ({100*(S.v12_score>=T).mean():.1f}%)")
    say("\n### Secondary endpoints")
    b = S[(S.v12_k >= 2) & S.v11_full_role.notna() & S.v12_role.notna()]
    if len(b) > 5:
        say(f"- S1 role concordance, pairs aligned by both (N={len(b)}): v12 {mci(b.v12_role)} | v11 {mci(b.v11_full_role)} | random {mci(b.random_any_role)}; v12 vs v11 p={stats.wilcoxon(b.v12_role, b.v11_full_role).pvalue:.2g}")
    k3 = S[S.v12_k >= 3].dropna(subset=["v12_role", "random_any_role"])
    if len(k3) > 5: say(f"- S1 v12 alignments with ≥3 residues (N={len(k3)}): v12 {mci(k3.v12_role)} | v11 {mci(k3.v11_full_role)} | random {mci(k3.random_any_role)}; p(v12>random)={stats.wilcoxon(k3.v12_role, k3.random_any_role).pvalue:.2g}")
    rows = []
    for mode in ["random", "surface", "composition"]:
        a, lo, hi, n = auc_ci(real, Dd[Dd["mode"] == mode]); rows.append({"decoy": mode, "pairs": n, "AUC (95% CI)": f"{a:.3f} ({lo:.3f}–{hi:.3f})"})
    say("- S2 other decoys:"); say(md(pd.DataFrame(rows)))
    st = []
    for col2 in ["ec_stratum", "tm_stratum"]:
        for k, g in S.groupby(col2):
            sg = g[g.v12_score >= T]
            st.append({"stratum": f"{col2}={k}", "pairs": len(g), "significant %": round(100*len(sg)/len(g), 1), "role (significant)": mci(sg.v12_role), "random (same pairs)": mci(sg.random_any_role)})
    say("- S3 by stratum:"); say(md(pd.DataFrame(st)))
    say(f"- S3 significant matches between structurally dissimilar proteins (TM<0.3): {int(((S.v12_score>=T)&(S.tm<0.3)).sum())}")
    s4 = S[S.v12_score >= T]; say(f"- S4 significant alignments not sequence-ordered: {100*(s4.v12_ordered.astype(str)=='False').mean():.1f}% (N={len(s4)})")
    if col == "set_blast":
        top = S[(S.v12_score >= T) & (S.tm < 0.3)].sort_values("v12_score", ascending=False).head(10)
        say("\nTop significant matches between structurally dissimilar proteins (TM<0.3):")
        say(md(top[["query", "target", "ec_stratum", "tm", "v12_k", "v12_score", "v12_rmsd_fc", "v12_role", "v12_ordered", "v12_corr"]].round(3)))
open(f"{H}/VALIDATION_RESULTS.md", "w").write("\n".join(L) + "\n")
