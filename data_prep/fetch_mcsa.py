"""Download M-CSA curated catalytic-residue data (residues + EC numbers) for the benchmark.
Run on a computer with normal internet access:  python fetch_mcsa.py
Saves ../mcsa/curated_data.csv  (the M-CSA flat file with PDB id, chain, residue number, roles, EC)."""
import os, urllib.request
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "mcsa"); os.makedirs(OUT, exist_ok=True)
url = "https://www.ebi.ac.uk/thornton-srv/m-csa/media/flat_files/curated_data.csv"
dest = os.path.join(OUT, "curated_data.csv")
urllib.request.urlretrieve(url, dest)
print("saved", dest, os.path.getsize(dest), "bytes")
