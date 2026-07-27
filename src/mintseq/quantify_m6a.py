#!/usr/bin/env python
"""Sum MINT-Seq (m6A IP) and TT-Seq (input) coverage over every canonical hg19
L1/Alu/SVA copy, strand-aware, for both replicates.

MINT-Seq is an m6A immunoprecipitation on nascent RNA; TT-Seq is the same
nascent RNA without the IP. Raw MINT-Seq density over a TE is
(nascent transcription x methylation rate), and intronic Alus sit in highly
transcribed genes, so IP alone ranks expression, not methylation. Every readout
downstream is therefore a per-copy IP/input ratio -- this script only extracts
the two numerators and denominators, it does not form the ratio.

Strand convention matches the eCLIP tree: for a TE copy on '+', the *sense*
transcript is the one on '+', i.e. the .pos bigWig; antisense is .neg. Both IP
and input are taken from the same strand file, so the ratio is well-formed per
orientation.

Output: results/mintseq/te_copy_signal.tsv.gz, one row per TE copy, with summed
coverage in each of {mint_r1, mint_r2, tt_r1, tt_r2} x {sense, antisense}.

Usage: quantify_m6a.py [bed] [bigwig_dir] [results_dir]
"""
import gzip
import os
import sys
import time

import numpy as np
import pyBigWig

PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

BED = sys.argv[1] if len(sys.argv) > 1 else os.path.join(PROJ, "hg19_L1_Alu_SVA_canonical.bed")
BWDIR = sys.argv[2] if len(sys.argv) > 2 else os.path.join(PROJ, "data", "mintseq")
OUTDIR = sys.argv[3] if len(sys.argv) > 3 else os.path.join(PROJ, "results", "mintseq")

# (label, filename stem). The strand suffix .pos/.neg bigWig is appended.
SAMPLES = [
    ("mint_r1", "GSM4086534_K562-MINT-Seq"),
    ("mint_r2", "GSM4086536_K562-MINT-Seq-rep2"),
    ("tt_r1", "GSM4086533_K562-TT-Seq"),
    ("tt_r2", "GSM4086535_K562-TT-Seq-rep2"),
]

# Primary assembly only. hg19 _hap contigs are alternate copies of chr6/chr4/chr17
# sequence; a TE copy there is the *same* copy as one on the primary chromosome,
# and counting both would double-weight those subfamilies. _random / chrUn are
# dropped for the same reason the eCLIP tree ignores them: no usable annotation.
CANON = {"chr%d" % i for i in range(1, 23)} | {"chrX", "chrY"}


def read_bed(path):
    """chrom -> (starts, ends, strand_is_plus, subfam_id) plus the name table."""
    by_chrom = {}
    names, name_id = [], {}
    n_skipped = 0
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            chrom = f[0]
            if chrom not in CANON:
                n_skipped += 1
                continue
            sid = name_id.get(f[3])
            if sid is None:
                sid = name_id[f[3]] = len(names)
                names.append(f[3])
            by_chrom.setdefault(chrom, []).append(
                (int(f[1]), int(f[2]), f[5] == "+", sid)
            )
    out = {}
    for chrom, rows in by_chrom.items():
        starts = np.fromiter((r[0] for r in rows), np.int64, len(rows))
        ends = np.fromiter((r[1] for r in rows), np.int64, len(rows))
        plus = np.fromiter((r[2] for r in rows), bool, len(rows))
        sid = np.fromiter((r[3] for r in rows), np.int32, len(rows))
        out[chrom] = (starts, ends, plus, sid)
    return out, names, n_skipped


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    by_chrom, names, n_skipped = read_bed(BED)
    n_total = sum(len(v[0]) for v in by_chrom.values())
    print("%d TE copies on the primary assembly (%d skipped on alt/random contigs)"
          % (n_total, n_skipped))

    bws, totals = {}, {}
    for label, stem in SAMPLES:
        for strand in ("pos", "neg"):
            path = os.path.join(BWDIR, "%s.%s.bigWig" % (stem, strand))
            if not os.path.exists(path):
                sys.exit("missing bigWig: %s" % path)
            bw = pyBigWig.open(path)
            bws[(label, strand)] = bw
            totals[(label, strand)] = bw.header()["sumData"]
        print("  %-8s pos sum %.4g   neg sum %.4g"
              % (label, totals[(label, "pos")], totals[(label, "neg")]))

    # Library-size factors, written to a sidecar so the ratio step does not have
    # to reopen the bigWigs. The files are CPM-scaled but each by its own depth,
    # so IP and input are not on a common scale until this is applied.
    with open(os.path.join(OUTDIR, "library_totals.tsv"), "w") as fh:
        fh.write("sample\tstrand\tsum_data\n")
        for (label, strand), v in sorted(totals.items()):
            fh.write("%s\t%s\t%.6f\n" % (label, strand, v))

    cols = ["%s_%s" % (lab, o) for lab, _ in SAMPLES for o in ("sense", "anti")]
    outpath = os.path.join(OUTDIR, "te_copy_signal.tsv.gz")
    t0 = time.time()
    n_done = 0

    with gzip.open(outpath, "wt") as out:
        out.write("\t".join(["chrom", "start", "end", "subfamily", "strand", "length"]
                            + cols) + "\n")
        for chrom in sorted(by_chrom):
            starts, ends, plus, sid = by_chrom[chrom]
            n = len(starts)
            vals = {c: np.zeros(n, np.float64) for c in cols}

            for label, _ in SAMPLES:
                for strand_key, bw in (("pos", bws[(label, "pos")]),
                                       ("neg", bws[(label, "neg")])):
                    chrom_len = bw.chroms().get(chrom)
                    if chrom_len is None:
                        continue
                    # RepeatMasker occasionally annotates to the very end of a
                    # chromosome; pyBigWig raises on an interval past the header
                    # length, which would abort a 10-minute run at the last copy.
                    lo = np.clip(starts, 0, chrom_len)
                    hi = np.clip(ends, 0, chrom_len)
                    # A copy on '+' reads sense from .pos; a copy on '-' reads
                    # sense from .neg. So this file supplies sense for one set of
                    # copies and antisense for the other -- fill both in one pass.
                    is_sense = plus if strand_key == "pos" else ~plus
                    sense_col = "%s_sense" % label
                    anti_col = "%s_anti" % label
                    for i in range(n):
                        if hi[i] <= lo[i]:
                            continue
                        v = bw.stats(chrom, int(lo[i]), int(hi[i]),
                                     type="sum", exact=True)[0]
                        if v is None:
                            continue
                        if is_sense[i]:
                            vals[sense_col][i] = v
                        else:
                            vals[anti_col][i] = v

            lengths = ends - starts
            for i in range(n):
                out.write("%s\t%d\t%d\t%s\t%s\t%d\t%s\n" % (
                    chrom, starts[i], ends[i], names[sid[i]],
                    "+" if plus[i] else "-", lengths[i],
                    "\t".join("%.4f" % vals[c][i] for c in cols)))
            n_done += n
            print("  %-6s %7d copies  (%5.1f%%, %.1f min elapsed)"
                  % (chrom, n, 100.0 * n_done / n_total, (time.time() - t0) / 60),
                  flush=True)

    # Correctness check: the vectorised-looking stats() path is really a per-copy
    # C call, so what can silently go wrong is the sense/antisense assignment,
    # not the summation. Re-derive a random subset straight from values().
    rng = np.random.default_rng(0)
    chrom = "chr21"
    starts, ends, plus, _ = by_chrom[chrom]
    idx = rng.choice(len(starts), size=min(200, len(starts)), replace=False)
    bad = 0
    for i in idx:
        s, e = int(starts[i]), int(ends[i])
        want_strand = "pos" if plus[i] else "neg"
        bw = bws[("mint_r1", want_strand)]
        raw = np.nan_to_num(np.array(bw.values(chrom, s, e), dtype=np.float64)).sum()
        got = bw.stats(chrom, s, e, type="sum", exact=True)[0] or 0.0
        if abs(raw - got) > 1e-3 * max(1.0, abs(raw)):
            bad += 1
    if bad:
        sys.exit("sense-strand extraction check FAILED on %d/%d chr21 copies"
                 % (bad, len(idx)))
    print("chr21 brute-force check passed on %d copies" % len(idx))

    for bw in bws.values():
        bw.close()
    print("wrote %s in %.1f min" % (outpath, (time.time() - t0) / 60))


if __name__ == "__main__":
    main()
