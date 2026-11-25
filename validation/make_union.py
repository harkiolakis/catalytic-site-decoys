"""Combine the two leakage-filtered validation pair lists into validation_pairs_union.csv, flagging membership:
set_blast  - pairs kept by the BLAST E-value filter (validation_pairs_blast.csv; primary analysis)
set_prereg - pairs kept by the pre-registered %-identity filter (validation_pairs_final_preregistered.csv)
Downstream scripts are run once on the union; analyses then select a set by its flag.  python make_union.py"""
import os, pandas as pd
H = os.path.dirname(os.path.abspath(__file__))
b = pd.read_csv(f"{H}/validation_pairs_blast.csv", dtype={"query": str, "target": str})
p = pd.read_csv(f"{H}/validation_pairs_final_preregistered.csv", dtype={"query": str, "target": str})
u = pd.concat([b, p[~p.pair_idx.isin(b.pair_idx)]], ignore_index=True)
u["set_blast"] = u.pair_idx.isin(b.pair_idx).astype(int); u["set_prereg"] = u.pair_idx.isin(p.pair_idx).astype(int)
u = u.sort_values("pair_idx")[list(b.columns) + ["set_blast", "set_prereg"]]
u.to_csv(f"{H}/validation_pairs_union.csv", index=False); print(len(u), "pairs;", int(u.set_blast.sum()), "BLAST set;", int(u.set_prereg.sum()), "pre-registered set")
