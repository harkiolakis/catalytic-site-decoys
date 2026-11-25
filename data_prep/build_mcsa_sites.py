"""Build catalytic-residue lists (csa_sites.json), EC numbers (ec_by_pdb.json) and mechanistic roles (roles.json)
from the M-CSA curated flat file.  python build_mcsa_sites.py [repo_root]"""
import csv, json, collections, os, sys
ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sites = collections.defaultdict(lambda: collections.defaultdict(set))
roles = collections.defaultdict(set); ec = {}
for row in csv.reader(open(os.path.join(ROOT, "mcsa/curated_data.csv"))):
    if row[0] == "M-CSA ID": continue
    pdb = row[2].lower()
    if row[3]: ec.setdefault(pdb, set()).update(e.strip() for e in row[3].split(";") if e.strip())
    if row[4] != "residue": continue
    try: resid = int(row[7])
    except ValueError: continue
    ch = row[6]; sites[pdb][ch].add(resid)
    if len(row) > 9 and row[9]: roles[(pdb, ch, resid)].add(row[9])
out = {p: {c: sorted(v) for c, v in cs.items()} for p, cs in sites.items()}
json.dump(out, open(os.path.join(ROOT, "mcsa/csa_sites.json"), "w"))
json.dump({p: sorted(v) for p, v in ec.items()}, open(os.path.join(ROOT, "mcsa/ec_by_pdb.json"), "w"))
json.dump({f"{p}|{c}|{r}": sorted(v) for (p, c, r), v in roles.items()}, open(os.path.join(ROOT, "mcsa/roles.json"), "w"))
print("PDB entries with sites:", len(out), "| with EC:", len(ec), "| residues with roles:", len(roles))
