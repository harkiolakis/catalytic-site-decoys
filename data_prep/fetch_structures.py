"""Download the PDB files for the 924-pair catalytic-site benchmark.
Run on a computer with normal internet access (standard library only):
    python fetch_structures.py
Files go to ../structures/<id>.pdb (falls back to .cif when no PDB-format file exists).
"""
import os, sys, time, urllib.request
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "structures"); os.makedirs(OUT, exist_ok=True)
ids = [l.strip() for l in open(os.path.join(HERE, "benchmark_pdb_ids.txt")) if l.strip()]
failed = []
for i, pid in enumerate(ids, 1):
    if any(os.path.exists(os.path.join(OUT, pid + e)) for e in (".pdb", ".cif")):
        continue
    for ext in (".pdb", ".cif"):
        try:
            urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pid.upper()}{ext}", os.path.join(OUT, pid + ext))
            break
        except Exception:
            continue
    else:
        failed.append(pid)
    if i % 50 == 0: print(f"{i}/{len(ids)}"); time.sleep(0.5)
print(f"done; {len(ids)-len(failed)} ok, {len(failed)} failed: {failed}")
