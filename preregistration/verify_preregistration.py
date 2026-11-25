"""Check the frozen files against the SHA-256 fingerprints recorded in FROZEN_CONFIG.json before validation.
The paths in FROZEN_CONFIG.json are those of the working tree at freezing time; byte-identical copies are in frozen_files/.
python verify_preregistration.py"""
import hashlib, json, os
H = os.path.dirname(os.path.abspath(__file__)); cfg = json.load(open(os.path.join(H, "FROZEN_CONFIG.json")))
print("frozen at", cfg["frozen_at"]); allok = True
for path, h in cfg["sha256"].items():
    f = os.path.join(H, "frozen_files", os.path.basename(path)); ok = hashlib.sha256(open(f, "rb").read()).hexdigest() == h
    allok &= ok; print(f"{'OK      ' if ok else 'MISMATCH'} {os.path.basename(path)}")
print("all fingerprints match" if allok else "one or more fingerprints do NOT match")
