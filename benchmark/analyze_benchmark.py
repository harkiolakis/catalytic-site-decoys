"""Benchmark tables: octant-descriptor aligner (v11) RMSD, same-input baselines, ablations, roles, decoys, EC and TM-score analyses.
Run from benchmark/ after run_experiments.py and run_decoys.py; writes BENCHMARK_RESULTS.md, T1–T5 CSVs and analysis_pairs.csv."""
import json, numpy as np, pandas as pd
from scipy import stats
rng = np.random.default_rng(42); B = 2000
E = pd.read_csv("exp_mcsa_current.csv"); D = pd.read_csv("decoys_mcsa_current.csv")
tm = pd.read_csv("benchmark_tm.csv").set_index("pair_idx").tm  # TM-align, normalised by query (compute_benchmark_tm.py)
ec = json.load(open("../mcsa/ec_by_pdb.json"))
CO = ["same_h", "same_t_diff_h", "diff_t"]; CL = {"same_h": "Same HSF", "same_t_diff_h": "Same T, diff HSF", "diff_t": "Different T"}
L = []
def say(s=""): print(s); L.append(str(s))
def md(df): return df.to_markdown(index=False)
def bci(x, f=np.median):
    x = np.asarray(x); i = rng.integers(0, len(x), (B, len(x))); v = np.apply_along_axis(f, 1, x[i]); return np.percentile(v, [2.5, 97.5])
def cliff(a, b): a, b = np.asarray(a), np.asarray(b); return ((a[:, None] > b).sum() - (a[:, None] < b).sum()) / (len(a) * len(b))
def summ(name, r):
    r = np.asarray(r, float); r = r[~np.isnan(r)]; lo, hi = bci(r)
    return {"set": name, "N": len(r), "median RMSD (95% CI)": f"{np.median(r):.2f} ({lo:.2f}–{hi:.2f})", "IQR": f"{np.percentile(r,25):.2f}–{np.percentile(r,75):.2f}",
            "<2 Å %": round(100*(r<2).mean(),1), "<3 Å %": round(100*(r<3).mean(),1), "<5 Å %": round(100*(r<5).mean(),1)}

say("# Benchmark results (current M-CSA annotations)")
say(f"Benchmark: 924 original pairs; {(E.error=='too_few_cat').sum()} lack ≥2 annotated catalytic residues on both chains in current M-CSA; "
    f"{(E.error=='v11_no_alignment').sum()} produce no alignment; **{E.error.isna().sum()} aligned pairs** analysed (failures excluded from all RMSD statistics).")
A = E[E.error.isna()].copy(); A["tm"] = A.pair_idx.map(tm)
A["nc"] = pd.cut(A.k, [1, 2, 3, 4, 99], labels=["2", "3", "4", ">4"])
say(f"Site size: query median {A.m.median():.0f} (IQR {A.m.quantile(.25):.0f}–{A.m.quantile(.75):.0f}), target median {A.n.median():.0f}; aligned correspondences median {A.k.median():.0f}.")
say()
say("## Table 1. Mode C (v11) site RMSD overall, by CATH relationship and by number of aligned residues")
rows = [summ("All pairs", A.v11_full_rmsd)] + [summ(CL[c], A[A.cath_relation == c].v11_full_rmsd) for c in CO]
rows += [summ(f"{k} aligned", A[A.nc == k].v11_full_rmsd) for k in ["2", "3", "4", ">4"]] + [summ("≥3 aligned", A[A.k >= 3].v11_full_rmsd)]
t1 = pd.DataFrame(rows); t1.to_csv("T1_rmsd.csv", index=False); say(md(t1)); say()
g = {c: A[A.cath_relation == c].v11_full_rmsd.values for c in CO}
H, p = stats.kruskal(*g.values()); say(f"Kruskal–Wallis across CATH groups: H={H:.2f}, p={p:.3g}, ε²={(H-2)/(len(A)-3):.3f}")
for a, b in [("diff_t", "same_h"), ("same_t_diff_h", "same_h"), ("diff_t", "same_t_diff_h")]:
    x, y = g[a], g[b]; i1 = rng.integers(0, len(x), (B, len(x))); i2 = rng.integers(0, len(y), (B, len(y)))
    dd = np.median(x[i1], 1) - np.median(y[i2], 1); lo, hi = np.percentile(dd, [2.5, 97.5])
    say(f"- {CL[a]} vs {CL[b]}: Δmedian {np.median(x)-np.median(y):.2f} Å (95% CI {lo:.2f}–{hi:.2f}), Cliff's δ {cliff(x,y):.3f}, MW p={stats.mannwhitneyu(x,y).pvalue:.3g}; {100*(x<np.median(y)).mean():.0f}% of {CL[a]} below {CL[b]} median")
say()

say("## Table 2. Same-input baselines (same catalytic residues, same number of correspondences k as Mode C)")
rows = []
for col, lab in [("v11_full_rmsd", "Octant-descriptor aligner (v11)"), ("random_any_median", "Random one-to-one mapping (median of 200)"),
                 ("random_mono_median", "Random order-preserving mapping (median of 200)"), ("oracle_mono_rmsd", "Best order-preserving mapping (exhaustive)"),
                 ("oracle_any_rmsd", "Best unrestricted mapping (exhaustive/100k sample)")]:
    rows.append(summ(lab, A[col]))
t2 = pd.DataFrame(rows); t2.to_csv("T2_baselines.csv", index=False); say(md(t2))
w = stats.wilcoxon(A.v11_full_rmsd, A.random_any_median)
say(f"\nMode C vs random mapping: paired Wilcoxon p={w.pvalue:.2g}; Mode C lower in {100*(A.v11_full_rmsd<A.random_any_median).mean():.1f}% of pairs.")
eq = np.isclose(A.v11_full_rmsd, A.oracle_mono_rmsd, atol=1e-3)
say(f"Mode C finds the best order-preserving mapping in {100*eq.mean():.1f}% of pairs; median percentile of Mode C among all mappings: {100*A.v11_pct_vs_any.median():.1f}% beat it.")
nm = (~A.oracle_any_monotone.astype(bool)) & (A.oracle_any_rmsd < A.v11_full_rmsd - 0.5)
say(f"The best unrestricted mapping is NOT sequence-order-preserving and ≥0.5 Å better than Mode C in {100*nm.mean():.1f}% of pairs "
    f"({100*nm[A.cath_relation=='diff_t'].mean():.1f}% of different-topology pairs) — the order-preserving Smith–Waterman cannot recover permuted sites.")
say()

say("## Table 3. Ablation of Mode C components (same pairs; each configuration's own k)")
cf = [("v11_full", "Full v11 (72-slot conn. + 24-slot exact + BLOSUM + distance + refinement)"), ("v11_conn_only", "72-slot connectivity channel only (no exact-offset channel)"),
      ("v10_exact", "v10 exact-offset descriptor (2-sphere)"), ("v11_no_blosum", "v11 without BLOSUM"), ("v11_no_dist", "v11 without distance-compatibility term"),
      ("v11_no_refine", "v11 without iterative refinement"), ("blosum_dist_only", "No descriptor (BLOSUM + distance)"), ("dist_only", "Distance term only"), ("blosum_only", "BLOSUM only")]
rows = []
for c, lab in cf:
    r = A[f"{c}_rmsd"]; n = A[f"{c}_n"]
    sub = A[r.notna() & A.v11_full_rmsd.notna()]
    same_k = sub[sub[f"{c}_n"] == sub.v11_full_n]
    rows.append({"configuration": lab, "aligned": int(r.notna().sum()), "median k": n[r.notna()].median(), "median RMSD": round(r.median(), 2),
                 "<2 Å %": round(100*(r<2).sum()/len(A), 1), "same k as full: N": len(same_k),
                 "identical correspondence RMSD %": round(100*np.isclose(same_k[f'{c}_rmsd'], same_k.v11_full_rmsd, atol=1e-3).mean(), 1) if len(same_k) else np.nan,
                 "median RMSD at same k (config vs full)": f"{same_k[f'{c}_rmsd'].median():.2f} vs {same_k.v11_full_rmsd.median():.2f}" if len(same_k) else ""})
t3 = pd.DataFrame(rows); t3.to_csv("T3_ablation.csv", index=False); say(md(t3)); say()

say("## Table 4. Functional correctness: fraction of aligned residue pairs sharing ≥1 M-CSA mechanistic role")
rows = []
for col, lab in [("v11_full_role", "Mode C (v11)"), ("v10_exact_role", "v10 descriptor"), ("dist_only_role", "Distance term only"), ("blosum_only_role", "BLOSUM only"),
                 ("random_any_role", "Random mapping (mean of 50)"), ("oracle_any_role", "Lowest-RMSD unrestricted mapping")]:
    sub = A[list(dict.fromkeys([col, "random_any_role"]))].dropna()
    if not len(sub): continue
    lo, hi = bci(sub[col].values, np.mean)
    rows.append({"method": lab, "N pairs": len(sub), "mean role concordance (95% CI)": f"{sub[col].mean():.3f} ({lo:.3f}–{hi:.3f})",
                 "vs random (paired Wilcoxon p)": "" if col == "random_any_role" else f"{stats.wilcoxon(sub[col], sub.random_any_role).pvalue:.2g}"})
t4 = pd.DataFrame(rows); t4.to_csv("T4_roles.csv", index=False); say(md(t4))
for c in CO:
    s = A[A.cath_relation == c][["v11_full_role", "random_any_role"]].dropna()
    say(f"- {CL[c]}: Mode C {s.v11_full_role.mean():.3f} vs random {s.random_any_role.mean():.3f} (N={len(s)})")
say()
say("## Side-chain (functional-group centroid) RMSD under the Mode C correspondence")
say(md(pd.DataFrame([summ("Cα RMSD", A.v11_full_rmsd), summ("Side-chain centroid RMSD", A.v11_sc_rmsd)])))
say(f"Spearman(Cα, side-chain RMSD) = {stats.spearmanr(A.v11_full_rmsd, A.v11_sc_rmsd, nan_policy='omit').statistic:.2f}")
say()

say("## Table 5. Decoy controls (target catalytic residues replaced; 5 trials per pair per mode; pair-level bootstrap CIs)")
D = D[D.real_n > 0]
def auc(r, d): return stats.mannwhitneyu(d, r).statistic / (len(r) * len(d))
def pboot(g):
    real = g.groupby("pair_idx").real_rmsd.first(); dec = {p: x.decoy_rmsd.values for p, x in g.groupby("pair_idx")}; ps = real.index.values; v = []
    for _ in range(1000):
        s = rng.choice(ps, len(ps)); v.append(auc(real.loc[s].values, np.concatenate([dec[q] for q in s])))
    return np.percentile(v, [2.5, 97.5])
rows = []
for mode in ["random", "surface", "sequence", "composition", "compact", "spatial", "spatial_seq"]:
    for lab, f in [("all", lambda x: x), ("k=2", lambda x: x[x.real_n == 2]), ("k≥3", lambda x: x[x.real_n >= 3])]:
        g = f(D[(D["mode"] == mode) & D.decoy_rmsd.notna()])
        if g.pair_idx.nunique() < 20: continue
        real = g.groupby("pair_idx").real_rmsd.first(); dm = g.groupby("pair_idx").decoy_rmsd.mean(); lo, hi = pboot(g)
        rows.append({"decoy": mode, "stratum": lab, "pairs": len(real), "real median": round(real.median(), 2), "decoy median": round(dm.median(), 2),
                     "AUC (95% CI)": f"{auc(real.values, g.decoy_rmsd.values):.3f} ({lo:.3f}–{hi:.3f})", "Wilcoxon p (per pair)": f"{stats.wilcoxon(dm.loc[real.index], real).pvalue:.1g}",
                     "real <2 Å %": round(100*(real<2).mean(), 1), "decoy <2 Å %": round(100*(g.decoy_rmsd<2).mean(), 1),
                     "Rg real / decoy (Å)": f"{g.groupby('pair_idx').real_rg.first().median():.1f} / {g.groupby('pair_idx').decoy_rg.mean().median():.1f}"})
t5 = pd.DataFrame(rows); t5.to_csv("T5_decoys.csv", index=False); say(md(t5))
cov = {m: x.groupby("pair_idx").decoy_rmsd.apply(lambda y: y.notna().any()).mean() for m, x in D.groupby("mode")}
say("\nPairs with ≥1 valid decoy: " + ", ".join(f"{m} {100*v:.0f}%" for m, v in cov.items()))
sp = D[(D["mode"] == "spatial") & D.decoy_rmsd.notna()]
say(f"Spatial decoys: distance tolerance used — " + ", ".join(f"{t:g} Å: {100*(sp.tol==t).mean():.0f}%" for t in sorted(sp.tol.dropna().unique())))
say()

say("## EC analysis (current M-CSA EC numbers)")
def shared(e1, e2):
    b = 0
    for x in e1:
        for y in e2:
            px, py = x.split("."), y.split("."); k = 0
            for i in range(min(len(px), len(py))):
                if px[i] == py[i] and px[i] not in ("-", "", "n"): k += 1
                else: break
            b = max(b, k)
    return b
A["ec_shared"] = [shared(ec[q.lower()], ec[t.lower()]) if q.lower() in ec and t.lower() in ec else np.nan for q, t in zip(A["query"], A["target"])]
say(f"Pairs with EC on both sides: {A.ec_shared.notna().sum()}/{len(A)}")
say(md(A.dropna(subset=["ec_shared"]).groupby("ec_shared").agg(N=("v11_full_rmsd", "size"), median_RMSD=("v11_full_rmsd", "median"),
    pct_lt2=("v11_full_rmsd", lambda x: round(100*(x<2).mean(), 1)), role_conc=("v11_full_role", "mean")).round(3).reset_index()))
say(pd.crosstab(A.cath_relation, A.ec_shared, margins=True).reindex(CO + ["All"]).to_string())
dt = A[A.cath_relation == "diff_t"]
say(f"Different-topology pairs sharing EC first digit: {(dt.ec_shared>=1).sum()} of {dt.ec_shared.notna().sum()} (manuscript §2.4 states none)")
e = A.dropna(subset=["ec_shared"]); r = stats.spearmanr(e.ec_shared, e.v11_full_rmsd)
say(f"Spearman(shared EC digits, RMSD) = {r.statistic:.3f}, p={r.pvalue:.2g}")
s2 = e[e.v11_full_rmsd < 2]
say(f"Sub-2 Å pairs: {100*(s2.ec_shared==4).mean():.0f}% identical EC, {100*(s2.ec_shared==0).mean():.0f}% different EC class; of the latter, {100*(s2[s2.ec_shared==0].k==2).mean():.0f}% are 2-residue alignments")
same_ec = e[e.ec_shared == 4]; say(f"Identical-EC pairs NOT reaching <5 Å: {100*(same_ec.v11_full_rmsd>=5).mean():.0f}% (N={len(same_ec)})")
say()
say("## Global TM-score vs site RMSD")
say(f"Spearman(TM-score, site RMSD) = {stats.spearmanr(A.tm, A.v11_full_rmsd, nan_policy='omit').statistic:.3f}")
s = A[(A.cath_relation == "diff_t") & (A.v11_full_rmsd < 2)]
say(f"Sub-2 Å different-topology pairs: {len(s)}; with k=2: {(s.k==2).sum()}; TM-score ≥0.5: {(s.tm>=0.5).sum()}")
cand = A[(A.cath_relation == "diff_t") & (A.k >= 4)].sort_values("v11_full_rmsd").head(10)[["query", "target", "k", "v11_full_rmsd", "v11_sc_rmsd", "v11_full_role", "tm", "ec_shared"]]
say("\nBest different-topology alignments with ≥4 residues (case-study candidates):"); say(md(cand.round(3)))
cand.to_csv("casestudy_candidates.csv", index=False)
A.to_csv("analysis_pairs.csv", index=False)
open("BENCHMARK_RESULTS.md", "w").write("\n".join(L) + "\n")
