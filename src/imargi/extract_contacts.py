#!/usr/bin/env python
"""Stream the pairs keeping the DNA END, for contacts whose RNA end is in a TE.

WHY THIS EXISTS SEPARATELY FROM extract_rna_ends.py. That script stores only the
distance CLASS of each contact, which is what makes a 1000-shift null affordable
-- the class is a property of the contact and survives moving the TE anchors, so
the null reduces to re-slicing sorted arrays. It is also, deliberately, a
throw-away of the DNA end's position, and that is exactly the field any question
about WHERE an RNA goes needs. Aggregate-fraction tests (trans_test.py,
within_intron_test.py) cannot see a focused target set by construction; this is
the extraction that lets one be looked for.

WHY IT IS SMALL. Only contacts whose RNA end falls inside an indexed TE copy are
kept -- a few percent of 248 M pairs -- so the output is a working set rather
than a re-encoding of the input.

WHAT IT DOES NOT GIVE. Pairs are (RNA, DNA), never (DNA, DNA). There is one DNA
anchor per record, so no contact matrix, no loop calls, no compartment calling
comes out of this. The DNA ends define an RNA's NEIGHBOURHOOD. DNA-DNA questions
need src/hic.

BOTH MAPQs ARE KEPT AND BOTH MATTER, for different reasons.
  mapq (RNA end) corrupts the ANCHOR: a mismapped RNA end attributes the contact
    to the wrong TE copy. Measured in this dataset to manufacture a q = 4.7e-05
    distal hit that inverts at MAPQ >= 30 -- see the README.
  mapq_dna (DNA end) corrupts the TARGET, which only starts to matter now.
    The homotypic test ("does SVA RNA land on SVA DNA") is maximally exposed:
    both ends repeat-derived is the classic Hi-C repeat artifact and produces
    precisely that signal. Filter both ends and report the comparison.

Stored per contact:

    rna_pos    int32   RNA end, 1-based
    rna_strand uint8   1 if the RNA end aligned '+', else 0 -- RAW, inverted
                       relative to transcript strand; see extract_rna_ends.py
    mapq       uint8   RNA end MAPQ
    copy       int32   row index into te_copy_key.tsv of the containing copy
    dna_chrom  uint8   index into the chroms array saved alongside
    dna_pos    int32   DNA end, 1-based
    mapq_dna   uint8   DNA end MAPQ

Copies nest (an Alu inside an L1), so a position can fall in several. The
INNERMOST-BY-START copy wins -- the last one whose start precedes the read --
which is one searchsorted and matches how the rest of the project resolves the
candidate window. Copies are not deduplicated upstream, so this is a convention,
not a fact about the annotation.

Usage: extract_contacts.py [pairs.gz ...] [outdir]
"""
import os
import subprocess
import sys
from collections import defaultdict

import numpy as np

PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
KEY = os.path.join(PROJ, "data", "chrrna", "te_copy_key.tsv")

if len(sys.argv) > 2:
    PAIRS_LIST, OUT = sys.argv[1:-1], sys.argv[-1]
else:
    PAIRS_LIST = [sys.argv[1]] if len(sys.argv) > 1 else [
        os.path.join(PROJ, "data", "hic", "4DNFIGDJIRV3.pairs.gz")]
    OUT = os.path.join(PROJ, "data", "imargi_contacts")

KEEP = set("chr%s" % c for c in list(range(1, 23)) + ["X", "Y"])
FLUSH = 4_000_000

DT = np.dtype([("rna_pos", np.int32), ("rna_strand", np.uint8),
               ("mapq", np.uint8), ("copy", np.int32),
               ("dna_chrom", np.uint8), ("dna_pos", np.int32),
               ("mapq_dna", np.uint8)])


def load_copies():
    """Per-chromosome sorted copy starts/ends plus their row index in the key.

    The row index, not the subfamily name, is what gets stored: it joins back to
    te_copy_key.tsv AND to te_copy_class.tsv, so subfamily, family, strand and
    host-gene class all stay available downstream without duplicating any of
    them into every one of tens of millions of records.
    """
    per = defaultdict(lambda: ([], [], []))
    with open(KEY) as fh:
        next(fh)
        for i, ln in enumerate(fh):
            f = ln.rstrip("\n").split("\t")
            per[f[3]][0].append(int(f[4]))
            per[f[3]][1].append(int(f[5]))
            per[f[3]][2].append(i)
    out = {}
    for c, (s, e, idx) in per.items():
        s = np.array(s, np.int64)
        o = np.argsort(s, kind="stable")
        out[c] = (s[o], np.array(e, np.int64)[o], np.array(idx, np.int64)[o],
                  int((np.array(e, np.int64) - s).max()) + 1)
    return out


def main():
    for f in PAIRS_LIST:
        if not os.path.exists(f):
            sys.exit("missing %s" % f)
    os.makedirs(OUT, exist_ok=True)
    raw = os.path.join(OUT, "raw")
    os.makedirs(raw, exist_ok=True)
    for f in os.listdir(raw):
        os.remove(os.path.join(raw, f))

    copies = load_copies()
    chroms = sorted(KEEP)
    cidx = {c: i for i, c in enumerate(chroms)}
    print("TE copies loaded for %d sequences" % len(copies), flush=True)

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
        rec = defaultdict(lambda: ([], [], [], [], [], []))
        for ln in lines:
            f = ln.split("\t")
            if len(f) < 7:
                continue
            r = rec[f[0]]
            r[0].append(int(f[1]))        # rna pos
            r[1].append(f[2])             # dna chrom
            r[2].append(int(f[3]))        # dna pos
            r[3].append(f[4][0])          # rna strand
            r[4].append(int(f[5]))        # mapq1
            r[5].append(int(f[6]))        # mapq2
        n += len(lines)
        for c1, (rp, dc, dp, st, mq, mqd) in rec.items():
            if c1 not in KEEP or c1 not in copies:
                continue
            cs, ce, cix, maxlen = copies[c1]
            a = np.array(rp, np.int64)
            # Innermost-by-start containing copy, if any. maxlen bounds how far
            # back the candidate can begin -- the same window argument the
            # searchsorted overlap idiom uses elsewhere in this project.
            j = np.searchsorted(cs, a, side="right") - 1
            ok = j >= 0
            jj = np.clip(j, 0, cs.size - 1)
            ok &= (ce[jj] >= a) & (cs[jj] >= a - maxlen)
            if not ok.any():
                continue
            k = np.nonzero(ok)[0]
            dcn = np.array([cidx.get(x, 255) for x in dc], np.int64)[k]
            m = dcn != 255
            if not m.any():
                continue
            k = k[m]
            out = np.empty(k.size, DT)
            out["rna_pos"] = a[k]
            out["rna_strand"] = (np.array(st, dtype="U1")[k] == "+")
            out["mapq"] = np.array(mq, np.uint8)[k]
            out["copy"] = cix[jj[k]]
            out["dna_chrom"] = dcn[m]
            out["dna_pos"] = np.array(dp, np.int64)[k]
            out["mapq_dna"] = np.array(mqd, np.uint8)[k]
            buf[c1].append(out)
            nbuf[c1] += k.size
            kept += k.size
            if nbuf[c1] >= FLUSH:
                flush(c1)
        print("\r  %d pairs read, %d in TE copies" % (n, kept), end="", flush=True)
    p.stdout.close()
    p.wait()
    for c in list(buf):
        flush(c)
    for h in handles.values():
        h.close()
    print("\n\nsorting shards ...", flush=True)

    tot = 0
    for f in sorted(os.listdir(raw)):
        c = f[:-4]
        arr = np.fromfile(os.path.join(raw, f), dtype=DT)
        arr.sort(order="rna_pos", kind="stable")
        np.savez(os.path.join(OUT, "%s.npz" % c),
                 **{k: arr[k] for k in DT.names})
        assert np.all(np.diff(arr["rna_pos"]) >= 0), "%s not sorted" % c
        tot += arr.size
        os.remove(os.path.join(raw, f))
        print("  %-6s %10d contacts" % (c, arr.size), flush=True)
    os.rmdir(raw)
    np.save(os.path.join(OUT, "chroms.npy"), np.array(chroms))

    print("\npairs read            : %d" % n)
    print("RNA end inside a TE   : %d  (%.2f%%)" % (tot, 100.0 * tot / max(n, 1)))
    print("wrote %s/<chrom>.npz  + chroms.npy" % OUT)


if __name__ == "__main__":
    main()
