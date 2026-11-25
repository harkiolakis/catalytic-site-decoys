"""Download PDB files for the validation set (run on a computer with internet):  python fetch_validation_structures.py"""
import os, urllib.request
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "..", "structures_validation"); os.makedirs(OUT, exist_ok=True)
ids = [l.strip() for l in open(os.path.join(HERE, "validation_pdb_ids.txt")) if l.strip()]; failed = []
for pid in ids:
    if os.path.exists(os.path.join(OUT, pid + ".pdb")): continue
    try: urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pid.upper()}.pdb", os.path.join(OUT, pid + ".pdb"))
    except Exception: failed.append(pid)
print(f"done; {len(ids)-len(failed)} present, failed (likely mmCIF-only): {failed}")
