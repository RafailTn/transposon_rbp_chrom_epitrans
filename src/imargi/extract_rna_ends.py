#!/usr/bin/env python
"""Stream the iMARGI pairs once into compact per-chromosome RNA-end arrays.

Separated from trans_test.py because this pass reads 7.35 GB and takes ~5 min,
while the test itself must be re-run many times (nulls, strata, thresholds).

WHAT IS STORED, AND WHY IT IS SUFFICIENT.
Per RNA end, three things:

    pos      int32   genomic position of the RNA end (1-based, sorted)
    strand   uint8   1 if the RNA end aligned '+', else 0  -- RAW, see below
    dclass   uint8   0 cis-proximal (<=1 Mb), 1 cis-distal, 2 trans
    mapq     uint8   MAPQ of the RNA end (mapq1)
    mapq_dna uint8   MAPQ of the DNA end (mapq2)

WHY MAPQ IS KEPT. The 4DN pipeline ran `pairtools parse --min-mapq 1`, which
drops only genuinely tied alignments (MAPQ 0); rows at MAPQ 4-9 survive. That is
a deliberate and correct choice for repeat work -- it is why SVA coverage exists
at all -- but it leaves a specific artifact live, and only MAPQ can test it.

If an RNA end is MISMAPPED onto a TE copy from elsewhere in the genome, its DNA
end is still mapped correctly, near wherever the RNA actually came from. Measured
against the wrong anchor that DNA end looks far away, so mismapping manufactures
exactly a distal/trans signature. It is `mapq1`, the RNA end, that discriminates:
the anchor is what the mismapping corrupts. Replicate splits CANNOT test this --
mismapping is systematic, not stochastic, and reproduces perfectly across
replicates aligned to the same reference.

Both MAPQs are stored because the DNA-end quality is the relevant one for any
later question about contact TARGETS, and re-running this pass costs 5 minutes
per file.

`dclass` is precomputed here, and that is the trick that makes the shift null
affordable. The distance class is a property of the CONTACT -- the relationship
between the RNA end and its own DNA end -- and does not depend on which TE copy
we later attribute the RNA end to. So a toroidal shift of the TE anchors changes
*which* RNA ends are captured but never their distance class, and the whole null
reduces to re-slicing these sorted arrays. Recomputing overlaps per shift against
248 M pairs would be hopeless; re-slicing is milliseconds.

STRAND IS STORED RAW. check_pairs.py measured that strand1 is INVERTED relative
to the transcript strand on this library (4/4 housekeeping genes, 300-800x). That
mapping is applied downstream, not baked in here, so the arrays stay valid if the
convention is ever re-measured on a different file. trans_test.py re-derives it
and asserts.

Coordinates are 1-based throughout (.pairs is 1-based, and so is te_copy_key.tsv
and GENCODE). Nothing is converted anywhere in this stage.

Usage: extract_rna_ends.py [pairs.gz] [outdir]
"""
import os
import subprocess
import sys
from collections import defaultdict

import numpy as np

PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# Accepts SEVERAL pairs files: the last argument is the output dir, everything
# before it is input. A biorep is deposited as two pairs files and they must be
# pooled into one array set, since the whole point of a replicate split is that
# the two halves are independent of EACH OTHER, not internally subdivided.
if len(sys.argv) > 2:
    PAIRS_LIST = sys.argv[1:-1]
    OUT = sys.argv[-1]
else:
    PAIRS_LIST = [sys.argv[1]] if len(sys.argv) > 1 else [
        os.path.join(PROJ, "data", "hic", "4DNFIGDJIRV3.pairs.gz")]
    OUT = os.path.join(PROJ, "data", "imargi")

PROX = 1_000_000
FLUSH = 8_000_000     # records buffered per chromosome before hitting disk

# Only the real chromosomes. Scaffolds and chrEBV carry no TE analysis and would
# just add empty shards.
KEEP = set("chr%s" % c for c in list(range(1, 23)) + ["X", "Y"])

DT = np.dtype([("pos", np.int32), ("strand", np.uint8), ("dclass", np.uint8),
               ("mapq", np.uint8), ("mapq_dna", np.uint8)])


def main():
    for f in PAIRS_LIST:
        if not os.path.exists(f):
            sys.exit("missing %s" % f)
    raw = os.path.join(OUT, "raw")
    os.makedirs(raw, exist_ok=True)
    for f in os.listdir(raw):
        os.remove(os.path.join(raw, f))

    buf = defaultdict(list)
    nbuf = defaultdict(int)
    handles = {}

    def flush(c):
        if not buf[c]:
            return
        a = np.concatenate(buf[c])
        if c not in handles:
            handles[c] = open(os.path.join(raw, "%s.bin" % c), "wb")
        handles[c].write(a.tobytes())
        buf[c] = []
        nbuf[c] = 0

    # cols 2-6 chrom1 pos1 chrom2 pos2 strand1, then 9,10 mapq1 mapq2.
    # cut emits them in FILE order regardless of how the list is written, so the
    # field indices below are 0..6 in that order.
    cmd = "zcat %s | grep -v '^#' | cut -f2,3,4,5,6,9,10" % " ".join(
        "'%s'" % f for f in PAIRS_LIST)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True,
                         bufsize=1 << 22)
    print("streaming %s ..." % ", ".join(os.path.basename(f) for f in PAIRS_LIST),
          flush=True)
    n = kept = 0
    while True:
        lines = p.stdout.readlines(1 << 26)
        if not lines:
            break
        rec = defaultdict(lambda: ([], [], [], [], []))
        for ln in lines:
            f = ln.split("\t")
            if len(f) < 7:
                continue
            r = rec[(f[0], f[2])]
            r[0].append(int(f[1]))      # rna pos
            r[1].append(int(f[3]))      # dna pos
            r[2].append(f[4][0])        # rna strand
            r[3].append(int(f[5]))      # mapq1, RNA end
            r[4].append(int(f[6]))      # mapq2, DNA end
        n += len(lines)
        for (c1, c2), (rp, dp, st, mq, mqd) in rec.items():
            if c1 not in KEEP:
                continue
            a = np.array(rp, np.int64)
            b = np.array(dp, np.int64)
            out = np.empty(a.size, DT)
            out["pos"] = a
            out["strand"] = (np.array(st, dtype="U1") == "+").astype(np.uint8)
            # bwa caps MAPQ at 60, so uint8 cannot overflow.
            out["mapq"] = np.array(mq, np.uint8)
            out["mapq_dna"] = np.array(mqd, np.uint8)
            if c1 == c2:
                out["dclass"] = np.where(np.abs(b - a) <= PROX, 0, 1).astype(np.uint8)
            else:
                out["dclass"] = 2
            buf[c1].append(out)
            nbuf[c1] += a.size
            kept += a.size
            if nbuf[c1] >= FLUSH:
                flush(c1)
        print("\r  %d pairs read, %d kept" % (n, kept), end="", flush=True)
    p.stdout.close()
    p.wait()
    for c in list(buf):
        flush(c)
    for h in handles.values():
        h.close()
    print("\n\nsorting shards by position ...", flush=True)

    meta = {}
    for f in sorted(os.listdir(raw)):
        c = f[:-4]
        a = np.fromfile(os.path.join(raw, f), dtype=DT)
        a.sort(order="pos", kind="stable")
        np.savez(os.path.join(OUT, "%s.npz" % c),
                 pos=a["pos"], strand=a["strand"], dclass=a["dclass"],
                 mapq=a["mapq"], mapq_dna=a["mapq_dna"])
        meta[c] = int(a.size)
        # Sorted order is what every downstream searchsorted depends on.
        assert np.all(np.diff(a["pos"]) >= 0), "%s not sorted" % c
        os.remove(os.path.join(raw, f))
        print("  %-6s %10d RNA ends" % (c, a.size), flush=True)
    os.rmdir(raw)

    with open(os.path.join(OUT, "meta.tsv"), "w") as fh:
        fh.write("chrom\tn_rna_ends\n")
        for c, k in meta.items():
            fh.write("%s\t%d\n" % (c, k))
    print("\ntotal pairs read : %d" % n)
    print("RNA ends stored  : %d  (%.1f%%, rest on scaffolds/chrEBV)"
          % (kept, 100.0 * kept / max(n, 1)))
    print("wrote %s/<chrom>.npz" % OUT)


if __name__ == "__main__":
    main()
