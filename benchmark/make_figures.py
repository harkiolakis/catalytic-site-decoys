"""Figures for 'Compactness, not catalytic geometry...'. Run from benchmark/. Writes fig1_auc_by_method.png, fig2_rmsd_by_size.png, fig3_roles_matched_fpr.png to argv[1]."""
import re, os, sys, json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
H = os.path.dirname(os.path.abspath(__file__)); OUT = sys.argv[1] if len(sys.argv) > 1 else "."
C = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "grid": "#e4e3df", "ink": "#2b2a27", "muted": "#6b6a64", "light": "#b9b8b2"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": C["muted"], "axes.labelcolor": C["ink"], "xtick.color": C["muted"], "ytick.color": C["ink"],
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": C["grid"], "grid.linewidth": 0.6, "axes.axisbelow": True})
def ci(s): m = re.match(r"([\d.]+) \(([\d.]+)–([\d.]+)\)", str(s).strip()); return tuple(map(float, m.groups()))
# ---- Figure 1: dumbbell, random vs compact decoys, by method
T = pd.read_csv(os.path.join(H, "T_decoy_baselines.csv")); AB = json.load(open(os.path.join(H, "abl_results.json")))
def get(setkey, meth):
    r = T[T.set.str.contains(setkey, regex=False) & (T.method == meth)].iloc[0]; return ci(r["random"]), ci(r["compact"])
def abl(s, v): d = AB[f"{s}|{v}"]; return tuple(d["random"][k] for k in ("auc", "lo", "hi")), tuple(d["compact"][k] for k in ("auc", "lo", "hi"))
rows = [("Random one-to-one mapping", "Random mapping (median)", "cα"), ("Best order-preserving mapping (oracle)", "Best order-preserving mapping", "cα"),
        ("Best unrestricted mapping (oracle)", "Best unrestricted mapping", "cα"), ("Distance-only aligner", "Distance-only aligner", "cα"),
        ("Octant-descriptor aligner", "v11 (octant descriptor)", "cα"), ("Side-chain matcher, Cα only, no chemistry", "abl:v12_ca_nochem", "cα"),
        ("Side-chain matcher, Cα only", "abl:v12_ca_only", "sc"), ("Side-chain matcher, no chemistry", "abl:v12_nochem", "sc"), ("Side-chain matcher", "v12 (two-point, order-free)", "sc"), ("Jess (residue-restricted)", "JESS (residue-restricted)", "sc")]
fig, axs = plt.subplots(1, 2, figsize=(7.2, 5.2), sharey=True)
y = np.arange(len(rows))[::-1].astype(float); y[len(rows) - 6:] -= 0  # spacing
for ax, (title, sk, ak) in zip(axs, [("Benchmark (691 pairs)", "Benchmark — decoy set 2", "benchmark"), ("Independent validation (420 pairs)", "Validation (BLAST-filtered set) — decoy set 2", "validation")]):
    for yy, (lab, key, grp) in zip(y, rows):
        try: r, c = abl(ak, key[4:]) if key.startswith("abl:") else get(sk, key)
        except (IndexError, KeyError): continue
        ax.plot([r[0], c[0]], [yy, yy], color=C["light"], lw=1.6, zorder=1)
        ax.errorbar(r[0], yy, xerr=[[r[0] - r[1]], [r[2] - r[0]]], fmt="o", ms=6.5, mfc="white", mec=C["blue"], mew=1.8, ecolor=C["blue"], elinewidth=1.4, zorder=3)
        ax.errorbar(c[0], yy, xerr=[[c[0] - c[1]], [c[2] - c[0]]], fmt="o", ms=6.5, color=C["orange"], ecolor=C["orange"], elinewidth=1.4, zorder=3)
    ax.axvline(0.5, color=C["muted"], lw=1, ls="--", zorder=2)
    ax.axhspan(y[5] - 0.5, y[0] + 0.5, color="#f3f2ef", zorder=-1, lw=0)
    ax.set_xlim(0.4, 0.95); ax.set_title(title, fontsize=9, color=C["ink"]); ax.grid(axis="y", visible=False); ax.set_xlabel("AUC, real site vs decoy (95% CI)")
axs[0].set_yticks(y); axs[0].set_yticklabels([r[0] for r in rows]); axs[0].set_ylim(y[-1] - 0.6, y[0] + 0.6)
axs[0].text(0.405, y[0] + 0.3, "Cα coordinates only", fontsize=7.5, color=C["muted"], va="center"); axs[0].text(0.405, y[6] + 0.3, "Chemistry and/or side-chain positions", fontsize=7.5, color=C["muted"], va="center")
axs[1].text(0.505, y[-1] - 0.45, "chance", fontsize=7.5, color=C["muted"])
fig.legend(handles=[Line2D([], [], marker="o", ls="", mfc="white", mec=C["blue"], mew=1.8, ms=7, label="vs random decoys"),
                    Line2D([], [], marker="o", ls="", color=C["orange"], ms=7, label="vs compactness-matched decoys")], loc="upper center", ncol=2, frameon=False, fontsize=8.5)
fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(os.path.join(OUT, "fig1_auc_by_method.png"), dpi=300); plt.close(fig)
# ---- Figure 2: RMSD vs aligned residues
A = pd.read_csv(os.path.join(H, "analysis_pairs.csv")); A["nc"] = pd.cut(A.k, [1, 2, 3, 4, 99], labels=["2", "3", "4", ">4"])
cols = [("random_any_median", "Random one-to-one mapping", C["muted"], "o"), ("v11_full_rmsd", "Octant-descriptor aligner", C["blue"], "s"),
        ("oracle_mono_rmsd", "Best order-preserving mapping", C["aqua"], "^"), ("oracle_any_rmsd", "Best unrestricted mapping", C["orange"], "D")]
fig, ax = plt.subplots(figsize=(6.0, 3.4)); x = np.arange(4)
for col, name, c, mk in cols:
    g = A.groupby("nc", observed=True)[col]; med = g.median().values; q1 = g.quantile(.25).values; q3 = g.quantile(.75).values
    ax.plot(x, med, marker=mk, color=c, lw=2, ms=6, label=name); ax.fill_between(x, q1, q3, color=c, alpha=0.10, lw=0)
ax.set_xticks(x); ax.set_xticklabels([f"{k}\n(n={n})" for k, n in A.nc.value_counts().reindex(["2", "3", "4", ">4"]).items()])
ax.set_xlabel("Aligned catalytic residues (k)"); ax.set_ylabel("Cα RMSD after superposition (Å)\nmedian, band = IQR"); ax.set_ylim(0, None)
ax.legend(frameon=False, fontsize=7.5, loc="upper left"); ax.grid(axis="x", visible=False)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig2_rmsd_by_size.png"), dpi=300); plt.close(fig)
# ---- Figure 3: role concordance at matched 1% FPR
def tables(md):
    out, cur = [], []
    for l in open(md).read().splitlines():
        if l.startswith("|"): cur.append(l)
        elif cur: out.append(cur); cur = []
    if cur: out.append(cur)
    return [pd.DataFrame([[x.strip() for x in r.strip("|").split("|")] for r in t[2:]], columns=[x.strip() for x in t[0].strip("|").split("|")]) for t in out]
tc = tables(os.path.join(H, "TOOL_COMPARISON.md")); rows = []
for setname, t in (("Benchmark", tc[1]), ("Independent validation", tc[3])):
    for _, r in t.iterrows():
        rows.append(dict(set=setname, tool=r["tool"], pct=int(r["decoy pct"]), called=r["real pairs called"], role=ci(r["role concordance of calls"]), rand=ci(r["random mapping, same pairs"])))
D = pd.DataFrame(rows); NAME = {"v12": "Side-chain matcher", "JESS": "Jess"}
fig, axs = plt.subplots(1, 2, figsize=(6.8, 3.2), sharey=True)
for ax, s in zip(axs, ["Benchmark", "Independent validation"]):
    sub = D[(D.set == s) & (D.pct == 99)]
    for i, (_, r) in enumerate(sub.iterrows()):
        c = C["orange"] if r.tool == "v12" else C["aqua"]
        ax.errorbar(i - 0.12, r.role[0], yerr=[[r.role[0] - r.role[1]], [r.role[2] - r.role[0]]], fmt="o", color=c, ms=7, lw=2)
        ax.errorbar(i + 0.12, r.rand[0], yerr=[[r.rand[0] - r.rand[1]], [r.rand[2] - r.rand[0]]], fmt="o", color=C["muted"], mfc="white", mew=2, ms=7, lw=2)
        ax.text(i, 0.06, f"calls: {r.called}", ha="center", fontsize=7.5, color=C["ink"])
    ax.set_xticks([0, 1]); ax.set_xticklabels([NAME.get(t, t) for t in sub.tool]); ax.set_xlim(-0.6, 1.6); ax.set_ylim(0, 0.8); ax.set_title(s, fontsize=9, color=C["ink"]); ax.grid(axis="x", visible=False)
axs[0].set_ylabel("Aligned pairs sharing ≥1 M-CSA role\n(mean, 95% CI)")
fig.legend(handles=[Line2D([], [], marker="o", color=C["orange"], ls="", label="Side-chain matcher calls"), Line2D([], [], marker="o", color=C["aqua"], ls="", label="Jess calls"),
                    Line2D([], [], marker="o", color=C["muted"], mfc="white", mew=2, ls="", label="Random mapping, same pairs")], loc="upper center", ncol=3, frameon=False, fontsize=7.5)
fig.tight_layout(rect=(0, 0, 1, 0.9)); fig.savefig(os.path.join(OUT, "fig3_roles_matched_fpr.png"), dpi=300); plt.close(fig)
print("figures written")
