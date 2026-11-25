"""Select the independent validation set: M-CSA entries sharing no PDB entry and no UniProt accession with the
924-pair benchmark; one chain per entry (most annotated catalytic residues, >=2); stratified pairs by EC relation.
Deterministic (seed 2026)."""
import csv, collections, itertools, json, random, os
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
rows = [r for r in csv.reader(open(f"{ROOT}/mcsa/curated_data.csv")) if r[0] != "M-CSA ID"]
bench = set(l.strip() for l in open(f"{ROOT}/data_prep/benchmark_pdb_ids.txt"))
ent = collections.defaultdict(lambda: {"pdb": set(), "uni": set(), "ec": set(), "res": collections.defaultdict(set)})
for r in rows:
    e = ent[r[0]]; e["pdb"].add(r[2].lower()); e["uni"].update(u.strip() for u in r[1].replace(";", ",").split(",") if u.strip())
    if r[3]: e["ec"].add(r[3].strip())
    if r[4] == "residue":
        try: e["res"][(r[2].lower(), r[6])].add(int(r[7]))
        except ValueError: pass
bench_ids = {k for k, e in ent.items() if e["pdb"] & bench}
bench_uni = set().union(*[ent[k]["uni"] for k in bench_ids])
sel = []
for k, e in sorted(ent.items()):
    if k in bench_ids or e["uni"] & bench_uni: continue
    (pdb, ch), res = max(e["res"].items(), key=lambda kv: (len(kv[1]), kv[0][1] == "A")) if e["res"] else ((None, None), set())
    if len(res) < 2: continue
    sel.append(dict(mcsa_id=k, pdb=pdb, chain=ch, n_cat=len(res), residues=" ".join(map(str, sorted(res))), ec=";".join(sorted(e["ec"])), uniprot=";".join(sorted(e["uni"]))))
def ec_rel(a, b):
    best = 0
    for x in a.split(";"):
        for y in b.split(";"):
            px, py = x.split("."), y.split("."); k = 0
            for i in range(min(len(px), len(py))):
                if px[i] == py[i] and px[i] not in ("-", ""): k += 1
                else: break
            best = max(best, k)
    return best
pairs = collections.defaultdict(list)
for a, b in itertools.combinations(range(len(sel)), 2):
    s = ec_rel(sel[a]["ec"], sel[b]["ec"]); stratum = "same_ec3" if s >= 3 else ("same_ec_class" if s >= 1 else "diff_ec_class")
    pairs[stratum].append((a, b))
rng = random.Random(2026); quota = {"same_ec3": 10**9, "same_ec_class": 400, "diff_ec_class": 500}; cap = 15
use = collections.Counter(); out = []
for stratum in ["same_ec3", "same_ec_class", "diff_ec_class"]:
    cand = pairs[stratum][:]; rng.shuffle(cand); n = 0
    for a, b in cand:
        if n >= quota[stratum]: break
        if use[a] >= cap or use[b] >= cap: continue
        q, t = (a, b) if rng.random() < 0.5 else (b, a)
        out.append(dict(pair_idx=len(out), query=sel[q]["pdb"], q_chain=sel[q]["chain"], target=sel[t]["pdb"], t_chain=sel[t]["chain"],
                        cath_relation="", ec_stratum=stratum, ec_shared=ec_rel(sel[q]["ec"], sel[t]["ec"]), q_mcsa=sel[q]["mcsa_id"], t_mcsa=sel[t]["mcsa_id"]))
        use[a] += 1; use[b] += 1; n += 1
V = f"{ROOT}/validation"
with open(f"{V}/validation_entries.csv", "w", newline="") as f:
    w = csv.DictWriter(f, sel[0].keys()); w.writeheader(); w.writerows(sel)
with open(f"{V}/validation_pairs.csv", "w", newline="") as f:
    w = csv.DictWriter(f, out[0].keys()); w.writeheader(); w.writerows(out)
json.dump({s["pdb"]: {s["chain"]: sorted(map(int, s["residues"].split()))} for s in sel}, open(f"{V}/validation_sites.json", "w"))
open(f"{V}/validation_pdb_ids.txt", "w").write("\n".join(sorted({s['pdb'] for s in sel})) + "\n")
print("entries", len(sel), "| pairs", len(out), collections.Counter(p["ec_stratum"] for p in out), "| available", {k: len(v) for k, v in pairs.items()})
print("entry usage: max", max(use.values()), "entries used", len(use))
