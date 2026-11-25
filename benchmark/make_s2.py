"""Build s2.json (Supplementary Table S2 blocks) and abl.json (manuscript numbers) from T_decoy_baselines.csv and abl_results.json"""
import pandas as pd, json, sys
T = pd.read_csv("T_decoy_baselines.csv"); AB = json.load(open("abl_results.json"))
NAME = {"Random mapping (median)": "Random one-to-one mapping", "Best order-preserving mapping": "Best order-preserving mapping (oracle)",
        "Best unrestricted mapping": "Best unrestricted mapping (oracle)", "Distance-only aligner": "Distance-only aligner", "v11 (octant descriptor)": "Octant-descriptor aligner",
        "v12 (two-point, order-free)": "Side-chain matcher", "JESS (residue-restricted)": "Jess (residue-restricted)", "EnzyMM (residue-restricted)": "EnzyMM (residue-restricted)"}
ABN = {"v12_ca_nochem": "Side-chain matcher, Cα only, no chemistry", "v12_ca_only": "Side-chain matcher, Cα only", "v12_nochem": "Side-chain matcher, no chemistry"}
SETN = {"Benchmark (691 pairs) — decoy set 1 (generator used for v11)": "Benchmark (691 pairs), decoy set 1",
        "Benchmark — decoy set 2 (generator used for v12, JESS and EnzyMM)": "Benchmark (691 pairs), decoy set 2",
        "Validation (BLAST-filtered set) — decoy set 1": "Validation (420 pairs), decoy set 1", "Validation (BLAST-filtered set) — decoy set 2": "Validation (420 pairs), decoy set 2"}
fam = ["random", "surface", "composition", "compact", "spatial"]; blocks = []
for s, g in T.groupby("set", sort=False):
    rows = [[NAME[m]] + [r[f] for f in fam] for m, r in zip(g.method, g.to_dict("records"))]
    if "decoy set 2" in s:
        k = "benchmark" if s.startswith("Benchmark") else "validation"
        for v, n in ABN.items():
            if f"{k}|{v}" in AB:
                d = AB[f"{k}|{v}"]; rows.append([n] + ["—", "—"] + ["—"] + [f"{d[f]['auc']:.3f} ({d[f]['lo']:.3f}–{d[f]['hi']:.3f})" if f in d else "—" for f in ("compact", "spatial")])
                rows[-1][1] = f"{d['random']['auc']:.3f} ({d['random']['lo']:.3f}–{d['random']['hi']:.3f})"
    blocks.append(dict(set=SETN[s], rows=rows))
json.dump(blocks, open(sys.argv[1] + "/s2.json", "w"), ensure_ascii=False, indent=1)
print(json.dumps({k: {f: v[f]["auc"] for f in v} for k, v in AB.items()}, indent=1))
