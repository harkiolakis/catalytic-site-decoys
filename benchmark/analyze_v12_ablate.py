"""AUC of side-chain matcher (v12) ablations vs decoys (same decoys/seeds as v12) -> ABLATION_DECOYS.md, abl_results.json"""
import pandas as pd, numpy as np, json, os
from scipy import stats
rng = np.random.default_rng(1); H = os.path.dirname(os.path.abspath(__file__))
def auc(r, d): return stats.mannwhitneyu(r, d).statistic / (len(r) * len(d))
def boot(R, D, col):
    R = R.dropna(); D = D[D.pair_idx.isin(R.index)].dropna(subset=[col]); R = R[R.index.isin(D.pair_idx)]
    dec = {p: g[col].values for p, g in D.groupby("pair_idx")}; ps = R.index.values; a = auc(R.values, D[col].values); v = []
    for _ in range(1000):
        s = rng.choice(ps, len(ps)); v.append(auc(R.loc[s].values, np.concatenate([dec[q] for q in s])))
    lo, hi = np.percentile(v, [2.5, 97.5]); return dict(auc=round(a, 3), lo=round(lo, 3), hi=round(hi, 3), pairs=int(len(R)))
U = pd.read_csv(os.path.join(H, "..", "validation", "validation_pairs_union.csv")); keep = set(U[U.set_blast == 1].pair_idx)
out = {}; L = ["# Side-chain matcher ablations against decoys\n", "| set | variant | random | compact | spatial |", "|---|---|---|---|---|"]
for lab, pre, kp in (("benchmark", "v12abl", None), ("validation", "val_v12abl", keep)):
    R = pd.read_csv(os.path.join(H, pre + "_real.csv")); D = pd.read_csv(os.path.join(H, pre + "_decoys.csv"))
    if os.path.exists(os.path.join(H, pre + "2_real.csv")):
        R2 = pd.read_csv(os.path.join(H, pre + "2_real.csv")); D2 = pd.read_csv(os.path.join(H, pre + "2_decoys.csv"))
        R = R.merge(R2[["pair_idx", "v12_ca_nochem_score"]], on="pair_idx", how="left")
        D = D.merge(D2[["pair_idx", "mode", "trial", "v12_ca_nochem_score"]], on=["pair_idx", "mode", "trial"], how="left")
    if kp is not None: R = R[R.pair_idx.isin(kp)]; D = D[D.pair_idx.isin(kp)]
    D = D[D.error.isna()] if "error" in D else D
    R = R.set_index("pair_idx")
    for v in [x for x in ("v12_ca_nochem", "v12_ca_only", "v12_nochem") if f"{x}_score" in R]:
        row = {}
        for m in ("random", "compact", "spatial"):
            row[m] = boot(R[f"{v}_score"], D[D["mode"] == m], f"{v}_score")
        out[f"{lab}|{v}"] = row
        L.append(f"| {lab} | {v} | " + " | ".join(f"{row[m]['auc']:.3f} ({row[m]['lo']:.3f}–{row[m]['hi']:.3f})" for m in ("random", "compact", "spatial")) + " |")
json.dump(out, open(os.path.join(H, "abl_results.json"), "w"), indent=1)
open(os.path.join(H, "ABLATION_DECOYS.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
