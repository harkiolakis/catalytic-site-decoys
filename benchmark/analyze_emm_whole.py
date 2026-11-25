"""Summarise EnzyMM whole-structure runs -> EMM_WHOLE.md"""
import pandas as pd, numpy as np, os, sys
from scipy import stats
H = os.path.dirname(os.path.abspath(__file__)); L = []
def say(s=""): L.append(s)
def pct(x): x = pd.Series(x).astype(float); return f"{100*x.mean():.1f}% ({int(x.sum())}/{len(x)})"
def wilson(k, n):
    if n == 0: return (np.nan, np.nan)
    z = 1.96; p = k / n; d = 1 + z*z/n; c = p + z*z/(2*n); r = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n)); return (max(0.0, (c-r)/d), min(1.0, (c+r)/d))
def rate(s):
    s = pd.Series(s).fillna(False).astype(bool); lo, hi = wilson(s.sum(), len(s)); return f"{100*s.mean():.1f}% [{100*lo:.1f}–{100*hi:.1f}] ({int(s.sum())}/{len(s)})"
say("# Enzyme Motif Miner in native whole-structure mode\n")
say("EnzyMM 0.5.0 with its default Jess parameters and templates of ≥3 side-chain residues, searched against the complete target "
    "structure (all chains, first model). 'Filtered' = matches the EnzyMM logistic-regression filter predicts correct. "
    "A match is on-site when any matched residue is an annotated M-CSA catalytic residue of the target (residue number, any chain).\n")
for lab, pref, keepf in (("Benchmark", "emmw_bench", "analysis_pairs.csv"), ("Validation (BLAST-filtered set)", "emmw_val", None)):
    if not os.path.exists(os.path.join(H, pref + "_pairs.csv")): continue
    P = pd.read_csv(os.path.join(H, pref + "_pairs.csv")); Lb = pd.read_csv(os.path.join(H, pref + "_library.csv"))
    if keepf:
        P = P.merge(pd.read_csv(os.path.join(H, "pairs_924.csv"))[["pair_idx", "cath_relation"]], on="pair_idx", how="left")
        keep = set(pd.read_csv(os.path.join(H, "decoys_mcsa_current.csv")).pair_idx.unique()); P = P[P.pair_idx.isin(keep)]
        tg = set(P.target); Lb = Lb[Lb.target.isin(tg)]
    else:
        U = pd.read_csv(os.path.join(H, "..", "validation", "validation_pairs_union.csv"), dtype={"query": str, "target": str})
        keep = set(U[U.set_blast == 1].pair_idx); P = P[P.pair_idx.isin(keep)].merge(U[["pair_idx", "ec_stratum"]], on="pair_idx", how="left"); Lb = Lb[Lb.target.isin(set(P.target))]
    say(f"## {lab}\n")
    say(f"### (A) Pair task: query-entry templates vs whole target structure ({len(P)} pairs; {int((P.A_ntemp>0).sum())} with ≥1 usable template)\n")
    A = P[P.A_ntemp > 0].reindex(columns=list(P.columns) + [c for c in ("A_pc_onsite", "A_any_onsite", "A_pc_k", "A_pc_rmsd", "A_pc_prec", "A_pc_recall") if c not in P], fill_value=np.nan)
    C = pd.concat([A.reindex(columns=[f"C{j}_any_hit", f"C{j}_pc_hit", f"C{j}_pc_onsite", f"C{j}_any_onsite"], fill_value=False).set_axis(["any_hit", "pc_hit", "pc_onsite", "any_onsite"], axis=1) for j in (0, 1)])
    rows = [dict(outcome="Any match (unfiltered)", real=rate(A.A_any_hit), control=rate(C.any_hit)),
            dict(outcome="Filtered match", real=rate(A.A_pc_hit), control=rate(C.pc_hit)),
            dict(outcome="Filtered match on the catalytic site", real=rate(A.A_pc_onsite), control=rate(C.pc_onsite))]
    say(pd.DataFrame(rows).rename(columns={"real": "query-entry templates", "control": "unrelated-entry templates (2 per pair, different EC class)"}).to_markdown(index=False))
    h = A[A.A_pc_hit == True]
    if len(h): say(f"\nFiltered matches: median {h.A_pc_k.median():.0f} residues, RMSD {h.A_pc_rmsd.median():.2f} Å; precision (share of matched residues annotated) {h.A_pc_prec.mean():.2f}; recall of the annotated site {h.A_pc_recall.mean():.2f}.")
    if "cath_relation" in P and P.cath_relation.notna().any():
        say("\nFiltered on-site match rate by CATH relation: " + ", ".join(f"{k} {rate(g.A_pc_onsite)}" for k, g in A.groupby("cath_relation")))
    if "ec_stratum" in P:
        say("\nFiltered on-site match rate by EC stratum: " + ", ".join(f"{k} {rate(g.A_pc_onsite)}" for k, g in A.groupby("ec_stratum")))
    say(f"\n### (B) Annotation task: whole library minus the target's own entry vs whole target structure ({len(Lb)} structures)\n")
    rows = [dict(outcome="≥1 filtered match", value=rate(Lb.B_pc_hit)),
            dict(outcome="Top filtered match on the catalytic site", value=rate(Lb.B_pc_onsite)),
            dict(outcome="Any filtered match on the catalytic site", value=rate(Lb.B_pc_anyonsite)),
            dict(outcome="Top filtered match off-site (residue-level false positive)", value=rate(Lb.B_pc_hit.fillna(False).astype(bool) & ~Lb.B_pc_onsite.fillna(False).astype(bool)))]
    say(pd.DataFrame(rows).to_markdown(index=False))
    hh = Lb[Lb.B_pc_hit == True]
    if len(hh):
        say(f"\nTop filtered match: precision {hh.B_pc_prec.mean():.2f}, recall {hh.B_pc_recall.mean():.2f}; matching template's entry shares EC to level ≥3 with the target in {rate(hh.B_pc_ec_shared >= 3)}; level 1 in {rate(hh.B_pc_ec_shared >= 1)}.")
    say()
open(os.path.join(H, "EMM_WHOLE.md"), "w").write("\n".join(L)); print("\n".join(L))
