import numpy as np, pandas as pd
from scipy import stats
rng = np.random.default_rng(7); B = 1000
R = pd.read_csv("v12_real.csv"); D = pd.read_csv("v12_decoys.csv"); E = pd.read_csv("analysis_pairs.csv")
M = R.merge(E[["pair_idx", "k", "v11_full_rmsd", "v11_sc_rmsd", "v11_full_role", "random_any_role", "oracle_any_role", "tm"]], on="pair_idx", how="left")
L = []
def say(s=""): print(s); L.append(str(s))
def md(df): return df.to_markdown(index=False)
def mci(x):
    x = np.asarray(x, float); x = x[~np.isnan(x)]; i = rng.integers(0, len(x), (B, len(x))); lo, hi = np.percentile(x[i].mean(1), [2.5, 97.5]); return f"{x.mean():.3f} ({lo:.3f}–{hi:.3f})"
say("# Side-chain matcher (v12) (permutation-tolerant, chemistry-aware) vs v11 Mode C — same 691 pairs, current M-CSA")
say(f"Pairs: {len(R)}. v12 returns an alignment (≥2 matched residues within 3 Å on both Cα and functional group) for {100*(R.v12_k>=2).mean():.1f}%; v11 Mode C for {100*M.v11_full_rmsd.notna().mean():.1f}% (v11 has no rejection threshold).")
say(f"v12 matched residues: median {R[R.v12_k>=2].v12_k.median():.0f}; k≥3 in {100*(R.v12_k>=3).mean():.1f}% of pairs. Non-sequence-ordered alignments: {100*(R[R.v12_k>=2].v12_ordered==False).mean():.1f}%.")
say()
say("## 1. Functional correctness (share ≥1 M-CSA mechanistic role per aligned residue pair)")
c = M[(M.v12_k >= 2) & M.v11_full_role.notna() & M.v12_role.notna() & M.random_any_role.notna()]
rows = []
for col, lab in [("v12_role", "v12"), ("v12_desc_role", "v12 + connectivity descriptor"), ("v12_nochem_role", "v12 without chemistry"), ("v12_ca_only_role", "v12 Cα only"),
                 ("v11_full_role", "v11 Mode C"), ("oracle_any_role", "Lowest-RMSD mapping (Cα)"), ("random_any_role", "Random mapping")]:
    x = c[col]; rows.append({"method": lab, "N": int(x.notna().sum()), "mean role concordance (95% CI)": mci(x),
                             "vs v11 (paired Wilcoxon p)": "" if col == "v11_full_role" else f"{stats.wilcoxon(x, c.v11_full_role, nan_policy='omit').pvalue:.2g}"})
say(f"Pairs aligned by both v12 and v11 with role data: {len(c)}"); say(md(pd.DataFrame(rows)))
for cr in ["same_h", "same_t_diff_h", "diff_t"]:
    s = c[c.cath_relation == cr]; say(f"- {cr}: v12 {s.v12_role.mean():.3f} | v11 {s.v11_full_role.mean():.3f} | random {s.random_any_role.mean():.3f} (N={len(s)})")
pp = M[M.v12_k >= 2]
say(f"\nAll v12 alignments (N={pp.v12_role.notna().sum()}): role concordance {mci(pp.v12_role)}; for k≥3 {mci(pp[pp.v12_k>=3].v12_role)}")
say()
say("## 2. Geometry of v12 alignments")
for lab, x in [("v12 Cα RMSD", pp.v12_rmsd_ca), ("v12 functional-group RMSD", pp.v12_rmsd_fc), ("v11 Cα RMSD (same pairs)", pp.v11_full_rmsd)]:
    x = x.dropna(); say(f"- {lab}: median {x.median():.2f} Å, <2 Å {100*(x<2).mean():.1f}% (N={len(x)})")
say()
say("## 3. Specificity: real site vs decoy target sites (score = normalised alignment score; decoys with no alignment score 0)")
def auc(r, d): return stats.mannwhitneyu(r, d).statistic / (len(r) * len(d))
rows = []
for v in ("v12", "v12_desc"):
    for mode in ["random", "surface", "composition", "compact", "spatial"]:
        for lab, f in [("all", lambda x: x), ("min site size ≥3", lambda x: x[np.minimum(x.m, x.n) >= 3]), ("min site size 2", lambda x: x[np.minimum(x.m, x.n) == 2])]:
            g = f(D[(D["mode"] == mode) & D[f"{v}_score"].notna()]); g = g[g.pair_idx.isin(R.pair_idx)]
            real = R.set_index("pair_idx")[f"{v}_score"].loc[g.pair_idx.unique()]
            if len(real) < 20: continue
            dec = {p: x[f"{v}_score"].values for p, x in g.groupby("pair_idx")}; ps = real.index.values; bs = []
            for _ in range(B):
                s = rng.choice(ps, len(ps)); bs.append(auc(real.loc[s].values, np.concatenate([dec[q] for q in s])))
            lo, hi = np.percentile(bs, [2.5, 97.5]); dm = g.groupby("pair_idx")[f"{v}_score"].mean()
            rows.append({"method": v, "decoy": mode, "stratum": lab, "pairs": len(real), "real mean score": round(real.mean(), 3), "decoy mean score": round(dm.mean(), 3),
                         "AUC (95% CI)": f"{auc(real.values, g[f'{v}_score'].values):.3f} ({lo:.3f}–{hi:.3f})", "Wilcoxon p": f"{stats.wilcoxon(real, dm.loc[real.index]).pvalue:.1g}",
                         "decoy aligned %": round(100*(g[f'{v}_k']>=2).mean(), 1), "real aligned %": round(100*(R.set_index('pair_idx').loc[real.index][f'{v}_k']>=2).mean(), 1)})
t = pd.DataFrame(rows); t.to_csv("T_v12_decoys.csv", index=False); say(md(t))
say("\nFor reference, v11 Mode C (RMSD-based): compact AUC 0.505, spatial 0.525, random 0.846.")
say()
say("## 4. Ablations (real pairs)")
rows = []
for v, lab in [("v12", "v12"), ("v12_desc", "+ connectivity descriptor (w=0.5)"), ("v12_nochem", "no chemistry"), ("v12_ca_only", "Cα only (no functional group)")]:
    a = R[R[f"{v}_k"] >= 2]
    rows.append({"variant": lab, "aligned %": round(100*len(a)/len(R), 1), "median k": a[f"{v}_k"].median(), "k≥3 %": round(100*(R[f'{v}_k']>=3).mean(), 1),
                 "role concordance": mci(a[f"{v}_role"]), "median Cα RMSD": round(a[f"{v}_rmsd_ca"].median(), 2)})
say(md(pd.DataFrame(rows)))
say()
say("## 5. Evolutionary relationship (v12 score)")
g = {c: R[R.cath_relation == c].v12_score.values for c in ["same_h", "same_t_diff_h", "diff_t"]}
H, p = stats.kruskal(*g.values()); say(f"Mean score: " + ", ".join(f"{k} {v.mean():.3f}" for k, v in g.items()) + f"; Kruskal–Wallis p={p:.2g}")
say(f"Aligned %: " + ", ".join(f"{k} {100*(R[R.cath_relation==k].v12_k>=2).mean():.0f}%" for k in g))
say()
say("## 6. Top different-topology v12 alignments with ≥3 matched residues (case-study candidates)")
cc = R[(R.cath_relation == "diff_t") & (R.v12_k >= 3)].sort_values("v12_score", ascending=False).head(12)
say(md(cc[["query", "target", "v12_k", "v12_score", "v12_rmsd_ca", "v12_rmsd_fc", "v12_role", "v12_ordered", "v12_corr"]].round(3)))
cc.to_csv("v12_casestudy_candidates.csv", index=False)
open("V12_RESULTS.md", "w").write("\n".join(L) + "\n")
