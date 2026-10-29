#!/usr/bin/env python
"""Boundary test: are cytoplasmic young-Alu sense reads autonomous Pol III
transcripts, or pieces of longer transcripts?

WHY. Under gene normalisation (11), Alu sense reads show an age gradient --
AluY cytoplasm-enriched, AluJ chromatin-enriched -- that is sense-only and holds
in contexts with no host transcription in the Alu's direction (intergenic,
intron_opp). One explanation: young Alus keep intact internal Pol III promoters
(A/B boxes) and free Alu RNAs are cytoplasmic. Pol III transcripts start at the
Alu's own 5' end and stop at the first T-run (>= 4 T) downstream. Read-through
by a longer transcript instead covers the flanks as densely as the body.

WINDOWS, per copy, all on the Alu's sense strand:
    UP    200 bp upstream of the Alu 5' end
    BODY  the Alu
    DN1   Alu 3' end -> end of the first TTTT on the sense strand (<= 300 bp)
    DN2   200 bp past that TTTT
and junction windows (2 bp, counted with --minOverlap 2, so a read must
contain both bases to count):
    J5    last UP base + first Alu base
    J3    last Alu base + first DN1 base
    JT    last T of the terminator + next base

Pol III predicts UP << BODY, J5 ~ 0, DN1 > DN2 and JT << J3.
Read-through predicts UP ~ BODY ~ DN1 ~ DN2 and J5 ~ J3 ~ JT.

COPY FILTER. Alu >= 280 bp (near full length, so the annotated 5' end is the
real Alu start where Pol III initiates); context intergenic, intron_opp, or
intron_same (the last as a read-through POSITIVE control: host pre-mRNA runs
through those in the Alu's direction); no other TE copy overlapping any
flank window (so flank reads are not repeat reads placed at random); DN1 >= 20
bp and a TTTT found within 300 bp.

Counting is per subfamily lineage (AluJ / AluS / AluY), summed over copies: per-
copy body counts are random partitions (see 03_align.sh), lineage totals are not.
Densities are fragments per kb of window, divided by sf_gene.

Writes data/boundary/windows.saf, junctions.saf, counts in work/boundary_counts/,
and results/boundary_test.tsv.

Usage: 16_boundary.py [bam_dir] [genome.fa]
"""
import bisect
import os
import subprocess
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)          # project root; scripts live in src/
BAMDIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "work", "bam")
FASTA = sys.argv[2] if len(sys.argv) > 2 else (
    "/mnt/raid1/rafail/aluaca/ref/GRCh38.primary_assembly.genome.fa")
DATA = os.path.join(ROOT, "data", "boundary")
CNT = os.path.join(ROOT, "work", "boundary_counts")
RES = os.path.join(ROOT, "results")
THREADS = os.environ.get("THREADS", "16")

FLANK, MAXTERM, MINLEN, MINDN1 = 200, 300, 280, 20
CONTEXTS = ("intergenic", "intron_opp", "intron_same")
WINDOWS = ("UP", "BODY", "DN1", "DN2")
JUNCS = ("J5", "J3", "JT")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def lineage(sub):
    for L in ("AluY", "AluS", "AluJ"):
        if sub.startswith(L):
            return L
    return None


def read_fasta(path, keep):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name in keep:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif name in keep:
                buf.append(line.rstrip())
    if name in keep:
        seqs[name] = "".join(buf).upper()
    return seqs


def samples():
    with open(os.path.join(ROOT, "config", "runs.tsv")) as fh:
        next(fh)
        return [l.split("\t")[0] for l in fh if l.strip()]


def build():
    os.makedirs(DATA, exist_ok=True)
    ctx = {}
    with open(os.path.join(ROOT, "data", "te_context.tsv")) as fh:
        next(fh)
        for line in fh:
            g, c = line.rstrip("\n").split("\t")
            ctx[g] = c
    key, te_iv = [], defaultdict(list)
    with open(os.path.join(ROOT, "data", "te_copy_key.tsv")) as fh:
        next(fh)
        for line in fh:
            gid, sub, fam, chrom, s, e, strand = line.rstrip("\n").split("\t")
            s, e = int(s), int(e)
            te_iv[chrom].append((s, e))
            if fam == "Alu" and e - s + 1 >= MINLEN and ctx[gid] in CONTEXTS and lineage(sub):
                key.append((gid, lineage(sub), ctx[gid], chrom, s, e, strand))
    starts, ends = {}, {}
    for c, v in te_iv.items():
        v.sort()
        starts[c] = [x[0] for x in v]
        # running max of ends lets one bisect answer "does anything overlap [a,b]"
        m, run = 0, []
        for x in v:
            m = max(m, x[1])
            run.append(m)
        ends[c] = run

    def te_overlaps(chrom, a, b):
        i = bisect.bisect_right(starts[chrom], b) - 1
        return i >= 0 and ends[chrom][i] >= a

    print("loading genome ...", flush=True)
    seqs = read_fasta(FASTA, {k[3] for k in key})

    win = open(os.path.join(DATA, "windows.saf"), "w")
    jun = open(os.path.join(DATA, "junctions.saf"), "w")
    for fh in (win, jun):
        fh.write("GeneID\tChr\tStart\tEnd\tStrand\n")
    ncopy = defaultdict(int)
    wlen = defaultdict(int)
    for gid, L, c, chrom, s, e, strand in key:
        seq = seqs.get(chrom)
        if seq is None:
            continue
        # Everything in sense coordinates: walk away from the Alu 3' end.
        if strand == "+":
            down = seq[e:e + MAXTERM]                         # 0-based e = 1-based e+1
            k = down.find("TTTT")
        else:
            down = seq[max(0, s - 1 - MAXTERM):s - 1].translate(COMP)[::-1]
            k = down.find("TTTT")
        if k < 0 or k + 4 < MINDN1:
            continue
        dn1 = k + 4                                           # through the last T
        if strand == "+":
            w = {"UP": (s - FLANK, s - 1), "BODY": (s, e),
                 "DN1": (e + 1, e + dn1), "DN2": (e + dn1 + 1, e + dn1 + FLANK)}
            j = {"J5": (s - 1, s), "J3": (e, e + 1), "JT": (e + dn1, e + dn1 + 1)}
        else:
            w = {"UP": (e + 1, e + FLANK), "BODY": (s, e),
                 "DN1": (s - dn1, s - 1), "DN2": (s - dn1 - FLANK, s - dn1 - 1)}
            j = {"J5": (e, e + 1), "J3": (s - 1, s), "JT": (s - dn1 - 1, s - dn1)}
        lo = min(a for a, _ in w.values())
        hi = max(b for _, b in w.values())
        if lo < 1 or hi > len(seq):
            continue
        # The Alu never overlaps its own flanks, so any hit is another TE copy.
        if any(te_overlaps(chrom, *w[n]) for n in ("UP", "DN1", "DN2")):
            continue
        grp = "%s|%s" % (L, c)
        ncopy[grp] += 1
        for n, (a, b) in w.items():
            win.write("%s|%s\t%s\t%d\t%d\t%s\n" % (grp, n, chrom, a, b, strand))
            wlen[grp + "|" + n] += b - a + 1
        for n, (a, b) in j.items():
            jun.write("%s|%s\t%s\t%d\t%d\t%s\n" % (grp, n, chrom, a, b, strand))
    win.close()
    jun.close()
    with open(os.path.join(DATA, "copies.tsv"), "w") as fh:
        fh.write("group\tn_copies\twindow\ttotal_bp\n")
        for grp in sorted(ncopy):
            for n in WINDOWS:
                fh.write("%s\t%d\t%s\t%d\n" % (grp, ncopy[grp], n, wlen[grp + "|" + n]))
    print("copies kept per lineage|context:")
    for grp in sorted(ncopy):
        print("  %-22s %7d" % (grp, ncopy[grp]))


def count():
    os.makedirs(CNT, exist_ok=True)
    for s in samples():
        bam = os.path.join(BAMDIR, "%s.Aligned.sortedByCoord.out.bam" % s)
        for kind, extra in (("windows", []), ("junctions", ["--minOverlap", "2"])):
            out = os.path.join(CNT, "%s.%s.txt" % (s, kind))
            if os.path.exists(out) and os.path.getsize(out):
                continue
            strand = os.environ.get("STRAND_" + s, os.environ.get("STRAND_DEFAULT", "2"))
            cmd = ["featureCounts", "-M", "-O", "-F", "SAF", "-T", THREADS, "-p",
                   "--countReadPairs", "-s", strand] + extra + [
                   "-a", os.path.join(DATA, kind + ".saf"), "-o", out, bam]
            print("== %s %s" % (s, kind), flush=True)
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def read_fc(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Geneid"):
                continue
            f = line.rstrip("\n").split("\t")
            out[f[0]] = float(f[-1])
    return out


def summarise():
    sf = {}
    with open(os.path.join(RES, "size_factors.tsv")) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sf[f[0]] = float(f[2])
    ncopy, wlen = {}, {}
    with open(os.path.join(DATA, "copies.tsv")) as fh:
        next(fh)
        for line in fh:
            g, n, w, bp = line.rstrip("\n").split("\t")
            ncopy[g] = int(n)
            wlen[g + "|" + w] = int(bp)
    groups = {"NEB_chrRNA": [s for s in samples() if "chrRNA" in s],
              "Takara_chromatin": [s for s in samples() if "_chr_Takara" in s],
              "Takara_cytoplasm": [s for s in samples() if "_cyto_" in s]}
    out = open(os.path.join(RES, "boundary_test.tsv"), "w")
    cols = ["lineage", "context", "fraction", "n_copies"] + \
           ["dens_" + w for w in WINDOWS] + ["raw_" + w for w in WINDOWS] + \
           ["raw_" + j for j in JUNCS] + ["UP_over_BODY", "DN2_over_DN1", "J5_over_J3", "JT_over_J3"]
    out.write("\t".join(cols) + "\n")
    print("\n%-5s %-11s %-17s %6s | %7s %7s %7s %7s | %6s %6s %6s | %7s %7s %6s %6s" % (
        "lin", "context", "fraction", "copies", "dUP", "dBODY", "dDN1", "dDN2",
        "J5", "J3", "JT", "UP/BODY", "DN2/DN1", "J5/J3", "JT/J3"))
    for grp in sorted(ncopy, key=lambda g: (g.split("|")[1], g)):
        L, c = grp.split("|")
        for frac, ss in groups.items():
            raw = defaultdict(float)
            dens = defaultdict(float)
            for s in ss:
                w = read_fc(os.path.join(CNT, "%s.windows.txt" % s))
                jn = read_fc(os.path.join(CNT, "%s.junctions.txt" % s))
                for n in WINDOWS:
                    v = w.get(grp + "|" + n, 0.0)
                    raw[n] += v
                    dens[n] += v / sf[s] / (wlen[grp + "|" + n] / 1000.0) / len(ss)
                for n in JUNCS:
                    raw[n] += jn.get(grp + "|" + n, 0.0)

            def r(a, b):
                return a / b if b > 0 else float("nan")
            ratios = [r(dens["UP"], dens["BODY"]), r(dens["DN2"], dens["DN1"]),
                      r(raw["J5"], raw["J3"]), r(raw["JT"], raw["J3"])]
            out.write("\t".join([L, c, frac, str(ncopy[grp])]
                                + ["%.4f" % dens[n] for n in WINDOWS]
                                + ["%d" % raw[n] for n in WINDOWS]
                                + ["%d" % raw[n] for n in JUNCS]
                                + ["%.3f" % x for x in ratios]) + "\n")
            print("%-5s %-11s %-17s %6d | %7.2f %7.2f %7.2f %7.2f | %6d %6d %6d | %7.2f %7.2f %6.2f %6.2f" % (
                L, c, frac, ncopy[grp], dens["UP"], dens["BODY"], dens["DN1"], dens["DN2"],
                raw["J5"], raw["J3"], raw["JT"], *ratios))
    out.close()
    print("\nwrote %s" % os.path.join(RES, "boundary_test.tsv"))


if __name__ == "__main__":
    if not os.path.exists(os.path.join(DATA, "copies.tsv")):
        build()
    count()
    summarise()
