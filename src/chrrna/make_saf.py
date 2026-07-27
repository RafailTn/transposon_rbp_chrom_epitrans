#!/usr/bin/env python
"""Build sense and antisense SAF annotations from hg38_L1_Alu_SVA_canonical.gtf.

featureCounts SAF is 'GeneID Chr Start End Strand', tab-separated, with a header
line, and its coordinates are 1-based inclusive -- the same convention the GTF
already uses, so the columns copy straight across with no offset. (Note this is
the one place in the project where no BED conversion happens: build_te_index.py
converts GTF -> BED, this converts GTF -> SAF, and SAF is GTF-like.)

Two files are written:

  te_sense.saf       strand exactly as annotated
  te_antisense.saf   strand flipped

Both are then counted with the SAME featureCounts -s flag. Flipping the
annotation rather than flipping -s is deliberate: it makes the orientation
explicit in the file being counted instead of resting on a reading of what
'-s 1' means for a reversely-stranded paired library, which is easy to get
backwards and silently inverts the whole sense/antisense split. Same class of
trap as the .pos/.neg bigWig convention in the mintseq tree.

GeneID is a positional id 'TE<n>', with the subfamily/coordinate lookup written
to te_copy_key.tsv. Encoding the subfamily into GeneID would mean parsing it back
out of featureCounts output, and the project already has a scar from
delimiter-packed identifiers (see the pipe-delimited GTF attribute trap).

Usage: make_saf.py [gtf] [outdir]
"""
import os
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
GTF = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    PROJ, "hg38_L1_Alu_SVA_canonical.gtf")
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(PROJ, "data", "chrrna")

FLIP = {"+": "-", "-": "+"}


def subfamily_of(attrs):
    """Exact-field lookup in the pipe-packed attribute string.

    'subfamily_id' contains 'family_id' as a substring, so an unanchored search
    for family_id matches the subfamily slot instead. Split and key off exact
    names -- never a regex.
    """
    a = attrs.split("|")
    for i in range(len(a) - 1):
        if a[i] == "subfamily_id":
            return a[i + 1]
    return None


def main():
    os.makedirs(OUT, exist_ok=True)
    sense = open(os.path.join(OUT, "te_sense.saf"), "w")
    anti = open(os.path.join(OUT, "te_antisense.saf"), "w")
    key = open(os.path.join(OUT, "te_copy_key.tsv"), "w")
    for fh in (sense, anti):
        fh.write("GeneID\tChr\tStart\tEnd\tStrand\n")
    key.write("gene_id\tsubfamily\tfamily\tchrom\tstart\tend\tstrand\n")

    n = 0
    n_bad = 0
    for line in open(GTF):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        chrom, start, end, strand, attrs = f[0], f[3], f[4], f[6], f[8]
        sfam = subfamily_of(attrs)
        if sfam is None or strand not in FLIP:
            n_bad += 1
            continue
        n += 1
        gid = "TE%d" % n
        sense.write("%s\t%s\t%s\t%s\t%s\n" % (gid, chrom, start, end, strand))
        anti.write("%s\t%s\t%s\t%s\t%s\n" % (gid, chrom, start, end, FLIP[strand]))
        fam = ("Alu" if sfam.startswith("Alu")
               else "L1" if sfam.startswith("L1") else "SVA")
        key.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                  % (gid, sfam, fam, chrom, start, end, strand))

    for fh in (sense, anti, key):
        fh.close()

    if n_bad:
        print("WARNING: skipped %d records with no subfamily_id or no strand" % n_bad)
    if n == 0:
        sys.exit("no records parsed from %s -- wrong file?" % GTF)
    print("wrote %d TE copies to:" % n)
    for name in ("te_sense.saf", "te_antisense.saf", "te_copy_key.tsv"):
        print("  %s" % os.path.join(OUT, name))


if __name__ == "__main__":
    main()
