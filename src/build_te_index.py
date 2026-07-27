"""Build a per-chromosome, strand-aware interval index of L1/Alu/SVA subfamilies.

Reads the canonical GTF and writes a pickle consumed by enrich_stranded.py.
Intervals are converted from GTF 1-based inclusive to
BED 0-based half-open so they can be compared directly against the peak BEDs.

Per chromosome we store starts/ends/subfamily-id/strand sorted by start, plus
maxlen (the longest element on that chromosome). maxlen bounds the candidate
window during overlap queries: an element starting before (peak_start - maxlen)
must end before peak_start and therefore cannot overlap.

Usage: python build_te_index.py [canonical.gtf] [out.pkl]
"""
import numpy as np, pickle, sys, os
from collections import defaultdict

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GTF = sys.argv[1] if len(sys.argv) > 1 else f"{PROJ}/hg38_L1_Alu_SVA_canonical.gtf"
OUT = sys.argv[2] if len(sys.argv) > 2 else f"{PROJ}/cache/te_index.pkl"

FAM2ELEM = {("SINE", "Alu"): "Alu", ("LINE", "L1"): "L1", ("Retroposon", "SVA"): "SVA"}

by_chrom = defaultdict(lambda: ([], [], [], []))
subfam_id, subfam_names, subfam_elem = {}, [], []

with open(GTF) as fh:
    for line in fh:
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        a = f[8].split("|")
        cls = fam = sub = ""
        for i in range(len(a) - 1):
            if a[i] == "class_id":
                cls = a[i + 1]
            elif a[i] == "family_id":
                fam = a[i + 1]
            elif a[i] == "subfamily_id":
                sub = a[i + 1]
        elem = FAM2ELEM.get((cls, fam))
        if elem is None:
            continue
        if sub not in subfam_id:
            subfam_id[sub] = len(subfam_names)
            subfam_names.append(sub)
            subfam_elem.append(elem)
        st, en, sid, sd = by_chrom[f[0]]
        st.append(int(f[3]) - 1)          # GTF 1-based inclusive -> BED 0-based
        en.append(int(f[4]))              # GTF end inclusive     -> BED half-open
        sid.append(subfam_id[sub])
        sd.append(1 if f[6] == "+" else -1)

index = {}
for c, (st, en, sid, sd) in by_chrom.items():
    st = np.asarray(st, np.int64); en = np.asarray(en, np.int64)
    sid = np.asarray(sid, np.int32); sd = np.asarray(sd, np.int8)
    o = np.argsort(st, kind="stable")
    index[c] = dict(starts=st[o], ends=en[o], sid=sid[o], strand=sd[o],
                    maxlen=int((en - st).max()))

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "wb") as fh:
    pickle.dump(dict(index=index, subfam_names=subfam_names, subfam_elem=subfam_elem), fh)

tot = sum(len(v["starts"]) for v in index.values())
nplus = sum(int((v["strand"] == 1).sum()) for v in index.values())
gtf_lines = sum(1 for _ in open(GTF) if not _.startswith("#"))
print(f"wrote {OUT}")
print(f"  chroms={len(index)} elements={tot} subfamilies={len(subfam_names)} "
      f"plus={nplus} minus={tot - nplus}")
assert tot == gtf_lines, f"indexed {tot} != {gtf_lines} GTF records"
print("  element count matches GTF: PASS")
