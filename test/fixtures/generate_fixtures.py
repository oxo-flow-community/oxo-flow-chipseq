#!/usr/bin/env python3
"""Generate the synthetic ChIP-seq fixtures for oxo-flow-chipseq.

The previous hand-made kit was 200 reads on a 1.8kb genome: bwa MAPQ
collapses, the dedup filter (-F 0x400, upstream keep_dups=false)
reduced everything to 7 unique reads, and SPP/macs2 halted on the
empty input. This generator emits, deterministically:

  references/genome.fa      chr1 (20kb) + chr2 (10kb)
  references/bwa_index/     BWA index of that genome (built by the
                            companion step — see README note below)
  references/chrom.sizes    UCSC chrom.sizes
  references/gene.bed       one annotated 'gene' region per peak
  raw/<SAMPLE>_R{1,2}.fastq.gz  8 samples (S1_REP1/2, S2_REP1/2,
                            C1_REP1/2, C2_REP1/2) x 5000 paired 100bp
                            reads, ~200bp inserts, 0.5% errors, ~25%
                            PCR-duplicated pairs at 2-7x total
                            multiplicity (consecutive levels — preseq
                            needs a gap-free histogram)

S1/S2 (ChIP) reads concentrate in 6 peak regions (~80% of S1/S2
pairs), C1/C2 (input) reads are uniform — macs2's
treatment-vs-control comparison has real signal. The S2* samples are
the second antibody's replicates for the multi-antibody consensus
tests (see test/fixtures/samples.tsv).

After regeneration, rebuild the BWA index (bwa index is not ported as
a rule — upstream PREPARE_GENOME is user-provided):

    bwa index test/fixtures/references/genome.fa
    mv genome.fa.* test/fixtures/references/bwa_index/
"""
import gzip
import io
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw")
REF = os.path.join(HERE, "references")
READ_LEN = 100
PAIRS_PER_SAMPLE = 5000
SEED = 13

CHR1_LEN = 20000
CHR2_LEN = 10000
# (chr, start, end) peak regions for the ChIP samples
PEAKS = [
    ("chr1", 2000, 2600),
    ("chr1", 5000, 5600),
    ("chr1", 9000, 9600),
    ("chr1", 14000, 14600),
    ("chr2", 1000, 1600),
    ("chr2", 7000, 7600),
]
# S1* (H3K4me3) and S2* (H3K27ac) are the two antibody replicate sets;
# C1*/C2* are their input controls. See test/fixtures/samples.tsv.
SAMPLES = [
    "S1_REP1",
    "S1_REP2",
    "C1_REP1",
    "C1_REP2",
    "S2_REP1",
    "S2_REP2",
    "C2_REP1",
    "C2_REP2",
]


def make_genome(rng):
    chr1 = "".join(rng.choice("ACGT") for _ in range(CHR1_LEN))
    chr2 = "".join(rng.choice("ACGT") for _ in range(CHR2_LEN))
    with open(os.path.join(REF, "genome.fa"), "w") as f:
        f.write(">chr1\n")
        for i in range(0, len(chr1), 80):
            f.write(chr1[i : i + 80] + "\n")
        f.write(">chr2\n")
        for i in range(0, len(chr2), 80):
            f.write(chr2[i : i + 80] + "\n")
    with open(os.path.join(REF, "chrom.sizes"), "w") as f:
        f.write(f"chr1\t{CHR1_LEN}\nchr2\t{CHR2_LEN}\n")
    with open(os.path.join(REF, "gene.bed"), "w") as f:
        for i, (chrom, s, e) in enumerate(PEAKS):
            f.write(f"{chrom}\t{s}\t{e}\tpeak_gene_{i}\t0\t+\t{s}\t{e}\t0\t1\t{e - s}\t0\n")
    return {"chr1": chr1, "chr2": chr2}


def mutate(seq, rng, rate=0.005):
    bases = list(seq)
    for i in range(len(bases)):
        if rng.random() < rate:
            bases[i] = rng.choice([b for b in "ACGT" if b != bases[i]])
    return "".join(bases)


def draw_start(genomes, is_chip, rng):
    if is_chip and rng.random() < 0.8:
        chrom, s, e = rng.choice(PEAKS)
        return chrom, rng.randrange(s, e - 250)
    chrom = rng.choice(["chr1", "chr2"])
    clen = len(genomes[chrom])
    return chrom, rng.randrange(0, clen - 250)


def write_sample(name, genomes, is_chip, rng):
    os.makedirs(RAW, exist_ok=True)
    r1_lines, r2_lines = [], []
    for i in range(PAIRS_PER_SAMPLE):
        chrom, start = draw_start(genomes, is_chip, rng)
        insert = rng.randint(180, 220)
        frag = genomes[chrom][start : start + insert]
        r1 = mutate(frag[:READ_LEN], rng)
        r2 = mutate(frag[-READ_LEN:][::-1].translate(str.maketrans("ACGT", "TGCA")), rng)
        rid = f"@{name}_{i}"
        r1_lines.append((rid + "/1", r1))
        r2_lines.append((rid + "/2", r2))
        # ~25% PCR duplicates at 2-7x total (consecutive levels)
        if i % 4 == 0:
            mult = rng.choice([1, 2, 3, 4, 5, 6])
            for k in range(mult):
                r1_lines.append((rid + f"_d{k}/1", r1))
                r2_lines.append((rid + f"_d{k}/2", r2))
    # mtime=0 keeps the gzip headers deterministic — regenerating the same
    # sample set reproduces byte-identical files (no git churn).
    with gzip.GzipFile(
        filename=os.path.join(RAW, f"{name}_R1.fastq.gz"), mode="wb", mtime=0
    ) as gz1, gzip.GzipFile(
        filename=os.path.join(RAW, f"{name}_R2.fastq.gz"), mode="wb", mtime=0
    ) as gz2:
        f1, f2 = io.TextIOWrapper(gz1, encoding="utf-8"), io.TextIOWrapper(
            gz2, encoding="utf-8"
        )
        for (h1, s1), (h2, s2) in zip(r1_lines, r2_lines):
            f1.write(f"{h1}\n{s1}\n+\n{'I' * READ_LEN}\n")
            f2.write(f"{h2}\n{s2}\n+\n{'I' * READ_LEN}\n")


def main():
    rng = random.Random(SEED)
    genomes = make_genome(rng)
    for name in SAMPLES:
        is_chip = name.startswith(("S1", "S2"))
        write_sample(name, genomes, is_chip, random.Random(SEED + len(name)))
    print("chipseq fixtures regenerated: 20kb/10kb genome, 8 samples x 5000 pairs, 6 peaks")
    print("next: rebuild the BWA index (see docstring)")


if __name__ == "__main__":
    main()
