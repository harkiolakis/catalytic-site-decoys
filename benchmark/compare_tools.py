"""JESS vs v12 (vs v11) under identical inputs, decoys and role-concordance evaluation."""
import numpy as np, pandas as pd
from scipy import stats
rng = np.random.default_rng(3); B = 1000
L = []
def say(s=""): print(s); L.append(str(s))
def md(df): return df.to_markdown(index=False)
def auc(r, d): return stats.mannwhitneyu(r, d).statistic / (len(r) * len(d))
def auc_ci(real, dec):
    dec = dec.dropna(subset=["score"]); g = {p: x.score.values for p, x in dec.groupby("pair_idx")}
    ps = np.array([p for p in real.index if p in g]); r = real.loc[ps]; bs = []
    for _ in range(B):
        s = rng.choice(ps, len(ps)); bs.append(auc(r.loc[s].values, np.concatenate([g[q] for q in s])))
    return f"{auc(r.values, np.concatenate([g[q] for q in ps])):.3f} ({np.percentile(bs,2.5):.3f}–{np.percentile(bs,97.5):.3f})"
def mci(x):
    x = np.asarray(x, float); x = x[~np.isnan(x)]
    if len(x) < 3: return "n/a"
    i = rng.integers(0, len(x), (B, len(x))); return f"{x.mean():.3f} ({np.percentile(x[i].mean(1),2.5):.3f}–{np.percentile(x[i].mean(1),97.5):.3f})"
SETS = {"Benchmark (training, 691 pairs)": dict(v12="v12_real.csv", v12d="v12_decoys.csv", j="jess_real.csv", jd="jess_decoys.csv", e="analysis_pairs.csv", keep=None),
        "Independent validation (BLAST-filtered)": dict(v12="../validation/val_v12_real.csv", v12d="../validation/val_v12_decoys.csv", j="../validation/val_jess_real.csv",
                                                       jd="../validation/val_jess_decoys.csv", e="../validation/val_exp.csv", keep="../validation/validation_pairs_blast.csv")}
say("# JESS vs Side-chain matcher (v12) under identical inputs")
say("Same pairs, same annotated residues on both sides, same decoy residue sets (identical generators and seeds), same role-concordance evaluation.")
say("JESS: templates from query catalytic residues (Cα + functional atom; chemical-class residue matching), searched only against the target's annotated or decoy residues; subtemplates of 2–4 residues; RMSD ≤ 2.0 Å, distance cut-off 1.5 Å; score = −log10 E-value of the best (largest, then lowest-E) hit.")
for name, f in SETS.items():
    V = pd.read_csv(f["v12"]); J = pd.read_csv(f["j"]); E = pd.read_csv(f["e"])
    VD = pd.read_csv(f["v12d"]); JD = pd.read_csv(f["jd"])
    if f["keep"]:
        k = set(pd.read_csv(f["keep"], dtype=str).pair_idx.astype(int)); V, J, E = V[V.pair_idx.isin(k)], J[J.pair_idx.isin(k)], E[E.pair_idx.isin(k)]
        VD, JD = VD[VD.pair_idx.isin(k)], JD[JD.pair_idx.isin(k)]
    common = sorted(set(V.pair_idx) & set(J.pair_idx)); V = V.set_index("pair_idx").loc[common]; J = J.set_index("pair_idx").loc[common]
    E = E.set_index("pair_idx").reindex(common)
    say(f"\n## {name}: {len(common)} pairs")
    rows = []
    for mode in ["random", "surface", "composition", "compact", "spatial"]:
        vd = VD[VD["mode"] == mode].rename(columns={"v12_score": "score"}); jd = JD[JD["mode"] == mode].rename(columns={"jess_score": "score"})
        rows.append({"decoy family": mode, "v12 AUC (95% CI)": auc_ci(V.v12_score, vd), "JESS AUC (95% CI)": auc_ci(J.jess_score, jd)})
    say("### Specificity (score-based AUC, real vs decoy)"); say(md(pd.DataFrame(rows)))
    say("\n### Calls at matched false-positive rate (threshold = 95th / 99th percentile of compact+spatial decoy scores in this set)")
    rows = []
    for q in (95, 99):
        for tool, real, dec, role, rr in [("v12", V.v12_score, pd.concat([VD[VD["mode"].isin(["compact", "spatial"])].v12_score]), V.v12_role, E.random_any_role),
                                          ("JESS", J.jess_score, pd.concat([JD[JD["mode"].isin(["compact", "spatial"])].jess_score]), J.jess_role, E.random_any_role)]:
            dec = dec.dropna(); t = np.percentile(dec, q)
            if tool == "JESS" and t <= 0: t = 1e-9  # any hit
            called = real >= t if tool == "v12" else real > t
            fp = (dec >= t).mean() if tool == "v12" else (dec > t).mean()
            x = pd.DataFrame({"role": role[called], "rand": rr[called]}).dropna()
            p = stats.wilcoxon(x.role, x.rand).pvalue if len(x) > 5 and (x.role != x.rand).any() else np.nan
            rows.append({"tool": tool, "decoy pct": q, "decoy FP rate": f"{100*fp:.1f}%", "real pairs called": f"{called.sum()} ({100*called.mean():.1f}%)",
                         "role concordance of calls": mci(x.role), "random mapping, same pairs": mci(x.rand), "p": f"{p:.2g}"})
    say(md(pd.DataFrame(rows)))
    hitJ = J.jess_k >= 2
    say(f"\n- JESS returns any hit for {100*hitJ.mean():.1f}% of real pairs vs {100*(JD[JD['mode'].isin(['compact','spatial'])].jess_k>=2).mean():.1f}% of compact/spatial decoys; v12 aligns {100*(V.v12_k>=2).mean():.1f}% of real pairs.")
    x = pd.DataFrame({"j": J.jess_role[hitJ], "r": E.random_any_role[hitJ], "v": V.v12_role[hitJ]}).dropna()
    say(f"- On pairs where JESS hits (N={len(x)}): JESS roles {mci(x.j)}, v12 roles {mci(x.v)}, random {mci(x.r)}")
    say(f"- JESS hit sizes: " + ", ".join(f"{int(k)} res: {int(n)}" for k, n in J[hitJ].jess_k.value_counts().sort_index().items()) + f"; non-sequence-ordered JESS hits: {100*(J[hitJ].jess_ordered.astype(str)=='False').mean():.1f}%")
    both = (J.jess_k >= 2) & (V.v12_score >= np.percentile(pd.concat([VD[VD['mode'].isin(['compact','spatial'])].v12_score]).dropna(), 95))
    say(f"- Overlap: pairs called by both {int(both.sum())}; JESS only {int(((J.jess_k>=2)&~both).sum())}; v12 (95th pct) only {int((~(J.jess_k>=2)&both.index.isin(V.index[V.v12_score>=np.percentile(pd.concat([VD[VD['mode'].isin(['compact','spatial'])].v12_score]).dropna(),95)])).sum())}")
open("TOOL_COMPARISON.md", "w").write("\n".join(L) + "\n")
