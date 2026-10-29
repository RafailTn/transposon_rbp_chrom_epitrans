#!/usr/bin/env python
"""Label every TE copy with its genomic context relative to host genes, and write
an intron SAF for the host-intron baseline.

WHY. Chromatin-vs-cytoplasm retention over TE reads is confounded by where the
TE sits. A sense Alu in a 3'UTR rides out to the cytoplasm inside a mature mRNA;
an L1 fragment in an intron is seen in chromatin RNA because the pre-mRNA around
it is. Without splitting by context, "Alu is exported, L1 is retained" largely
restates exon-versus-intron. This script supplies the split; 10_context.py
aggregates counts with it and 11_context_deseq2.R tests within each context.

Context is decided at the TE MIDPOINT, against protein_coding and lncRNA genes
only (pseudogenes and small-RNA genes would mislabel far more than they fix).
Levels, highest priority first:

    utr3      inside a 3'UTR (GENCODE 'UTR' lying 3' of its transcript's CDS)
    exon      inside any other exon
    intron    inside a gene body but not an exon
    intergenic

Each non-intergenic level is suffixed _same or _opp: whether the host feature
is on the same strand as the TE annotation. If both strands qualify at the top
level, _same wins -- a sense-strand host is the one that can carry the TE's
sense-oriented reads. Combined with the read orientation from the two SAFs this
gives e.g. 'sense reads of intron_same TEs' = host pre-mRNA running through the
TE in its own direction.

Writes:
  data/te_context.tsv      gene_id  context
  data/gencode_intron.saf  per-strand merged intronic intervals (gene body minus
                           exons of all gene types on that strand; primary
                           chromosomes only; rDNA masked), for 09_gene_count.sh

Usage: 08_context.py [gencode.gtf.gz] [datadir]
"""
import bisect
import gzip
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)          # project root; scripts live in src/
GTF = sys.argv[1] if len(sys.argv) > 1 else (
    "/mnt/raid1/rafail/aluaca/ref/gencode.v47.primary_assembly.annotation.gtf.gz")
DATA = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "data")
KEY = os.path.join(DATA, "te_copy_key.tsv")

GENE_TYPES = {"protein_coding", "lncRNA"}

# Masked from the intron SAF: the rDNA arrays. In GENCODE v47 lncRNA gene bodies
# span the 45S repeats on chr21 and on the unplaced rDNA contigs, so their
# "introns" are rRNA -- in the cytoplasm libraries (~35% rRNA) that was ~75% of
# all intronic fragments and flattened the intron baseline to log2FC ~0.
# Unplaced/alt contigs are dropped wholesale for the same reason.
RDNA_MASK = [("chr21", 8000000, 8500000)]


def attr(attrs, name):
    """Exact-key lookup in a GTF attribute column (no substring matching)."""
    for part in attrs.split(";"):
        part = part.strip()
        if part.startswith(name + " "):
            return part[len(name) + 1:].strip('"')
    return None


def merge(ivs):
    ivs.sort()
    out = []
    for s, e in ivs:
        if out and s <= out[-1][1] + 1:
            if e > out[-1][1]:
                out[-1][1] = e
        else:
            out.append([s, e])
    return out


def subtract(a, b):
    """a minus b; both merged, sorted, 1-based inclusive."""
    out = []
    j = 0
    for s, e in a:
        cur = s
        while j < len(b) and b[j][1] < cur:
            j += 1
        k = j
        while k < len(b) and b[k][0] <= e:
            if b[k][0] > cur:
                out.append([cur, b[k][0] - 1])
            cur = max(cur, b[k][1] + 1)
            k += 1
        if cur <= e:
            out.append([cur, e])
    return out


class Index:
    """Merged intervals per (chrom, strand), point lookup by bisect."""

    def __init__(self, ivs):
        self.starts, self.ends = {}, {}
        for k, v in ivs.items():
            m = merge(v)
            self.starts[k] = [x[0] for x in m]
            self.ends[k] = [x[1] for x in m]

    def hit(self, chrom, strand, pos):
        st = self.starts.get((chrom, strand))
        if not st:
            return False
        i = bisect.bisect_right(st, pos) - 1
        return i >= 0 and self.ends[(chrom, strand)][i] >= pos


def main():
    genes = defaultdict(list)
    exons = defaultdict(list)
    all_exons = defaultdict(list)
    cds_span = {}            # transcript -> (min CDS start, max CDS end)
    utrs = []                # (transcript, chrom, strand, start, end)
    opener = gzip.open if GTF.endswith(".gz") else open
    with opener(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            feat = f[2]
            if feat not in ("gene", "exon", "CDS", "UTR"):
                continue
            chrom, s, e, strand = f[0], int(f[3]), int(f[4]), f[6]
            if feat == "exon":
                # Exons of EVERY gene type are subtracted from introns, so small
                # RNAs nested in introns (RN7SL, snoRNAs...) are not counted as
                # intronic signal.
                all_exons[(chrom, strand)].append([s, e])
            if attr(f[8], "gene_type") not in GENE_TYPES:
                continue
            if feat == "gene":
                genes[(chrom, strand)].append([s, e])
            elif feat == "exon":
                exons[(chrom, strand)].append([s, e])
            elif feat == "CDS":
                t = attr(f[8], "transcript_id")
                lo, hi = cds_span.get(t, (s, e))
                cds_span[t] = (min(lo, s), max(hi, e))
            else:
                utrs.append((attr(f[8], "transcript_id"), chrom, strand, s, e))

    utr3 = defaultdict(list)
    for t, chrom, strand, s, e in utrs:
        span = cds_span.get(t)
        if span is None:
            continue
        if (strand == "+" and s > span[1]) or (strand == "-" and e < span[0]):
            utr3[(chrom, strand)].append([s, e])

    idx_utr3, idx_exon, idx_gene = Index(utr3), Index(exons), Index(genes)

    flip = {"+": "-", "-": "+"}
    tally = defaultdict(int)
    out = open(os.path.join(DATA, "te_context.tsv"), "w")
    out.write("gene_id\tcontext\n")
    with open(KEY) as fh:
        next(fh)
        for line in fh:
            gid, _, _, chrom, s, e, strand = line.rstrip("\n").split("\t")
            mid = (int(s) + int(e)) // 2
            ctx = "intergenic"
            for name, idx in (("utr3", idx_utr3), ("exon", idx_exon),
                              ("intron", idx_gene)):
                same = idx.hit(chrom, strand, mid)
                opp = idx.hit(chrom, flip[strand], mid)
                if same or opp:
                    ctx = name + ("_same" if same else "_opp")
                    break
            tally[ctx] += 1
            out.write("%s\t%s\n" % (gid, ctx))
    out.close()

    # Introns: gene body minus every exon on the same strand, so an intron of
    # one isoform that is an exon of another counts as exon.
    saf = open(os.path.join(DATA, "gencode_intron.saf"), "w")
    saf.write("GeneID\tChr\tStart\tEnd\tStrand\n")
    n = 0
    for k in sorted(genes):
        chrom, strand = k
        if not chrom.startswith("chr") or chrom == "chrM":
            continue
        name = "intron_plus" if strand == "+" else "intron_minus"
        mask = merge(all_exons.get(k, []) + [[ms, me] for mc, ms, me in RDNA_MASK if mc == chrom])
        for s, e in subtract(merge(genes[k]), mask):
            saf.write("%s\t%s\t%d\t%d\t%s\n" % (name, chrom, s, e, strand))
            n += 1
    saf.close()

    total = sum(tally.values())
    print("TE copies by context:")
    for k in sorted(tally, key=lambda x: -tally[x]):
        print("  %-14s %9d  (%.1f%%)" % (k, tally[k], 100.0 * tally[k] / total))
    print("wrote %s and %s (%d intronic intervals)"
          % (os.path.join(DATA, "te_context.tsv"),
             os.path.join(DATA, "gencode_intron.saf"), n))


if __name__ == "__main__":
    main()
