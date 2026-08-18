#!/usr/bin/env python3
"""Add PCR-duplicate structure to the chipseq fixtures.

preseq lc_extrap requires duplicate-count levels (live error: 'max count
before zero is less than min required count (4) duplicates removed').
This post-processor duplicates a subset of the existing genome-derived
pairs at varying multiplicities (2x/4x/8x/16x) so the duplicate-count
curve has >=4 levels.

Usage: python3 test/fixtures/add_duplicates.py   (overwrites raw/*.fastq.gz)
"""
import gzip
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw")
SEED = 11


def read_pairs(name):
    with gzip.open(os.path.join(RAW, f"{name}_R1.fastq.gz"), "rt") as f1, gzip.open(
        os.path.join(RAW, f"{name}_R2.fastq.gz"), "rt"
    ) as f2:
        pairs = []
        while True:
            h1 = f1.readline().strip()
            if not h1:
                break
            r1 = (h1, f1.readline().strip(), f1.readline().strip(), f1.readline().strip())
            h2 = f2.readline().strip()
            r2 = (h2, f2.readline().strip(), f2.readline().strip(), f2.readline().strip())
            pairs.append((r1, r2))
    return pairs


def write_pairs(name, pairs):
    with gzip.open(os.path.join(RAW, f"{name}_R1.fastq.gz"), "wt") as f1, gzip.open(
        os.path.join(RAW, f"{name}_R2.fastq.gz"), "wt"
    ) as f2:
        for r1, r2 in pairs:
            f1.write("\n".join(r1) + "\n")
            f2.write("\n".join(r2) + "\n")


def main():
    rng = random.Random(SEED)
    for sample in ("S1_REP1", "S1_REP2", "C1_REP1", "C1_REP2"):
        pairs = read_pairs(sample)
        if not pairs:
            print(f"{sample}: empty, skipping")
            continue
        out = []
        n = len(pairs)
        # ~30% of pairs get duplicated at 2x/4x/8x/16x (>=4 count levels).
        # Duplicates get unique read ids (PCR duplicates are separate clusters).
        dup_idx = set(rng.sample(range(n), int(n * 0.3)))
        for i, p in enumerate(pairs):
            out.append(p)
            if i in dup_idx:
                mult = rng.choice([1, 3, 7, 15])  # total multiplicity 2/4/8/16
                for k in range(mult):
                    (h1, s1, q1n, q1s), (h2, s2, q2n, q2s) = p
                    # tag goes BEFORE the /1 /2 suffix — bwa mem requires
                    # mates to share the exact read name after stripping
                    # the trailing /1 or /2 (live: 'paired reads have
                    # different names: "C1_REP2_7/1_d3", "C1_REP2_7/2_d3"').
                    tag = f"_d{k}"
                    out.append(((h1[:-2] + tag + h1[-2:], s1, q1n, q1s), (h2[:-2] + tag + h2[-2:], s2, q2n, q2s)))
        write_pairs(sample, out)
        print(f"{sample}: {n} -> {len(out)} reads")
    print("duplicates added")


if __name__ == "__main__":
    main()
