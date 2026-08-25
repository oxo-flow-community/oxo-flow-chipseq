# oxo-flow-chipseq — ChIP-seq peak calling, QC and differential analysis

> ★ Verified · ⇄ Official port of [`nf-core/chipseq`](https://github.com/nf-core/chipseq) @ `2.1.0` — same tools, same versions, same commands. Part of the [oxo-flow-community catalog](https://oxo-flow-community.github.io/).

[![CI](https://github.com/oxo-flow-community/oxo-flow-chipseq/actions/workflows/ci.yml/badge.svg)](https://github.com/oxo-flow-community/oxo-flow-chipseq/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

ChIP-seq analysis pipeline: raw-read QC (FastQC), adapter trimming (Trim
Galore), BWA-MEM (default) or STAR alignment, library merging, Picard
mark-duplicates, BAMTools filtering with blacklist and orphan-read removal,
library complexity (preseq, phantompeakqualtools SPP), bigWig tracks, deepTools
QC plots, MACS3 peak calling with input controls in broad (default) or narrow
mode, HOMER peak annotation, FRiP scoring, consensus peaks across replicates
(MACS3 merge + featureCounts quantification + DESeq2 QC), an IGV session and a
MultiQC report. Runs paired-end samples; `aligner` and `narrow_peak` select the
alignment/peak branches.

## Installation

### 1. Install oxo-flow

Requires **oxo-flow >= 0.12.0**. The recommended route is the release binary:

```bash
curl -fL -o oxo-flow.tar.gz https://github.com/Traitome/oxo-flow/releases/latest/download/oxo-flow-latest-x86_64-unknown-linux-gnu.tar.gz
tar xzf oxo-flow.tar.gz && sudo mv oxo-flow /usr/local/bin/
```

Alternatively via conda: `conda install -c bioconda oxo-flow-cli` (note the
bioconda package may lag behind releases; binaries for other platforms are on
the [releases page](https://github.com/Traitome/oxo-flow/releases)).

### 2. Get this workflow

```bash
git clone https://github.com/oxo-flow-community/oxo-flow-chipseq.git
```

### 3. Requirements

**Reference data you must provide** — set under `[config]` in `main.oxoflow`
(paths are relative to the run directory):

- genome FASTA and its FASTA index (`fasta`, `fai`)
- annotation GTF (`gtf`)
- gene-body regions in BED format (`gene_bed` — upstream derives this from
  the GTF; this port takes it pre-built, or generates it with
  `make_gene_bed = true`)
- chromosome sizes file (`chrom_sizes`)
- blacklist regions in BED format (`blacklist`)
- BWA index directory (`bwa_index`) containing the BWA index files (`*.amb`,
  `*.ann`, `*.bwt`, `*.pac`, `*.sa`) for the reference FASTA — for
  `aligner = "star"`, a STAR index directory (`star_index`) instead, or build
  one with `make_star_index = true` (then set `star_index = "results/star/index"`)
- raw paired-end FASTQ reads (`raw_dir`, named `raw/{pair_id}_R{1,2}.fastq.gz`)
  with sample metadata declared in `[[pairs]]` and `ip_ids`

**Compute** — the largest rules (`bwa_mem`, `star_align`, `trimgalore`,
`star_genomegenerate`) request 12 threads and 72 GB of memory; most other
rules request 6 threads and 36 GB (see `modules/*.oxoflow`
`[rules.resources]`). A default single-rule run peaks around 72 GB — scale
the numbers down for small machines.

**Tool delivery** — containers with pinned images. Every rule pins its
biocontainers (or nf-core / Seqera Wave) image and the workflow runs with the
**docker** backend, so you need Docker (or Singularity) installed.

## Usage

```bash
oxo-flow validate main.oxoflow
oxo-flow dry-run main.oxoflow       # prints the 186-instance plan (153 run, 33 gated off)
oxo-flow run main.oxoflow           # executes locally (docker backend)
oxo-flow debug main.oxoflow         # show expanded commands
```

Branch selection (each key toggles only its branch; the default plan is
unchanged when all are off):

| Key | Effect |
|---|---|
| `aligner = "star"` | STAR alignment instead of BWA; requires `star_index` (or `make_star_index = true` + `star_index = "results/star/index"`) |
| `narrow_peak = true` | whole peak chain in narrowPeak mode (MACS3 args, `narrow_peak/` dirs, consensus merge columns, IGV/MultiQC outputs) |
| `make_gene_bed = true` | derive `results/genome/gene.bed` from the GTF (then point `gene_bed` at it) |
| `make_blacklist_regions = true` | derive `results/genome/chrom.sizes.include_regions.bed` from blacklist + chrom sizes (then point `blacklist` at it) |

The workflow runs with the **docker** backend (all tools pinned to the
upstream container images). Replace the `test/fixtures` input paths in
`main.oxoflow` `[config]` with your own reads and references, and edit
`[[pairs]]`/`ip_ids` to match your samplesheet. Results are written to
`results/`.

## Source

Ported from **[nf-core/chipseq](https://github.com/nf-core/chipseq)** at tag
`2.1.0` (commit `76e2382b6d443db4dc2396e6831d1243256d80b0`), licensed MIT.
Created 2026-08-15; this workflow may lag behind upstream releases. See
[NOTICE.md](NOTICE.md) for attribution details and the upstream license note.

## Fidelity

Ported with upstream defaults: `aligner=bwa`, paired-end, `narrow_peak=false`
(broad peaks), `with_umi=false`. One row per upstream process; steps not
ported are listed with reasons.

| Upstream process | oxo-flow rule | Tool (version) | Notes |
|---|---|---|---|
| FASTQ_FASTQC_UMITOOLS_TRIMGALORE → FASTQC | `fastqc` | fastqc 0.12.1 | identical command (`--quiet --threads --memory`, 10GB cap); UMITOOLS_EXTRACT branch not ported (`with_umi=false` by default) |
| FASTQ_FASTQC_UMITOOLS_TRIMGALORE → TRIMGALORE | `trimgalore` | trim-galore 0.6.7 | identical (`--fastqc --cores n-4 --paired --gzip`, conditional `--nextseq/--clip_r1/--clip_r2/--three_prime_clip_*`) |
| BWA_MEM | `bwa_mem` | bwa 0.7.17, samtools 1.17 | identical (`-M -R '@RG...'`, secondary filter `-F 0x0100`, `-t` cores, `sort -T`); index lookup via `find` over `config.bwa_index`, same as upstream; gated on `aligner='bwa'` |
| STAR_GENOMEGENERATE (local) | `star_genomegenerate` | star 2.6.1d | gated (`aligner='star' && make_star_index`); identical command (`--genomeSAindexNbases` awk over the `.fai`); index written to `results/star/index/` |
| STAR_ALIGN (local) | `star_align` | star 2.6.1d | gated (`aligner='star'`); identical command (defaults + `--outSAMtype BAM Unsorted`, RG line with ID/SM + conditional CN); `Aligned.out.bam` renamed to `.Lb.bam` so the shared downstream chain is reused; logs published to `results/bwa/library/log/` for MultiQC auto-discovery |
| BAM_SORT_STATS_SAMTOOLS (library, SAMTOOLS_SORT + SAMTOOLS_INDEX + BAM_STATS_SAMTOOLS) | `sort_align` + `index_align` + `stats_align`/`flagstat_align`/`idxstats_align` | samtools 1.17 | identical commands (`samtools cat \| samtools sort` pipeline, `samtools index -@`, stats/flagstat/idxstats) |
| PICARD_MERGESAMFILES_LIBRARY | `mergesamfiles` | picard 3.2.0 | upstream symlink branch for single-library samples replicated exactly (`ln -s`); the multi-library MergeSamFiles branch is off the default path |
| BAM_MARKDUPLICATES_PICARD (PICARD_MARKDUPLICATES + SAMTOOLS_INDEX + BAM_STATS_SAMTOOLS) | `markduplicates` + `index_markdup` + `stats_markdup`/`flagstat_markdup`/`idxstats_markdup` | picard 3.2.0, samtools 1.17 | identical (`--ASSUME_SORTED true --REMOVE_DUPLICATES false --VALIDATION_STRINGENCY LENIENT --TMP_DIR tmp`, `XMX = memory*1024*8/10` heap) |
| BAM_FILTER_BAMTOOLS → BAMTOOLS_FILTER | `bamtools_filter` | samtools 1.17, bamtools 2.5.2 | identical (`-F 0x004 -F 0x0008 -f 0x001`, conditional `-F 0x0400`/`-q 1` on `keep_dups`/`keep_multi_map`, `-L blacklist`, `assets/bamtools_filter_pe.json`) |
| BAM_FILTER_BAMTOOLS → SAMTOOLS_SORT (name sort) | `sort_name` | samtools 1.17 | identical (`samtools cat \| samtools sort -n`, prefix `.mLb.flT.name_sorted`) |
| BAM_REMOVE_ORPHANS | `bam_remove_orphans` | python 3.8 | identical (`bampe_rm_orphan.py ... --only_fr_pairs`) |
| BAM_FILTER_BAMTOOLS → BAM_SORT_STATS_SAMTOOLS | `sort_filter` + `index_filter` + `stats_filter`/`flagstat_filter`/`idxstats_filter` | samtools 1.17 | identical commands, prefix `.mLb.clN.sorted` |
| PRESEQ_LCEXTRAP | `preseq` | preseq 3.2.0 | identical (`lc_extrap -verbose -bam -seed 1 -pe`, command log to stderr) |
| PICARD_COLLECTMULTIPLEMETRICS | `picard_collectmultiplemetrics` | picard 3.2.0 | identical (`-Xmx` heap, `--VALIDATION_STRINGENCY LENIENT --TMP_DIR tmp`, reference, `mv *.CollectMultipleMetrics.*`) |
| PHANTOM_PEAK_QUALTOOLS | `phantompeakqualtools` | r-base 3.5.1, phantompeakqualtools 1.2.2 | identical (`RUN_SPP=which run_spp.R`, `Rscript --max-ppsize=500000 -e "library(caTools); source(..)" -c= -savp= -savd= -out= -p=threads`) |
| MULTIQC_CUSTOM_PHANTOM_PEAK_QUALTOOLS | `multiqc_custom_phantompeakqualtools` | r-base 3.5.1 | identical (cross.correlation RData table, `$9`/`$10` NSC/RSC awk, header concat) |
| BEDTOOLS_GENOMECOV | `bedtools_genomecov` | bedtools 2.30.0 | identical (`-bg -scale 1e6/mapped -pc`, sort, scale-factor file) |
| UCSC_BEDGRAPHTOBIGWIG | `ucsc_bedgraphtobigwig` | ucsc-bedgraphtobigwig 445 | identical |
| DEEPTOOLS_COMPUTEMATRIX (scale-regions) | `deeptools_computematrix` | deeptools 3.5.5 | identical (regionBodyLength 1000, ±3000, `--missingDataAsZero --skipZeros --smartLabels`) |
| DEEPTOOLS_PLOTPROFILE / DEEPTOOLS_PLOTHEATMAP | `deeptools_plotprofile` / `deeptools_plotheatmap` | deeptools 3.5.5 | identical |
| DEEPTOOLS_PLOTFINGERPRINT | `deeptools_plotfingerprint` | deeptools 3.5.5 | identical (`--skipZeros --numberOfSamples 500000 --labels ip control`, paired bamfiles); control-only samples skipped like upstream |
| KHMER_UNIQUEKMERS | `khmer` | khmer 3.0.0a3 | identical (`unique-kmers.py -k read_length -R report`, `grep ^number`); gated on `macs_gsize` being empty, same as upstream |
| MACS2_CALLPEAK (nf-core) | `macs3_callpeak` | macs3 3.0.1 | identical flags (`--keep-dup all --broad --broad-cutoff 0.1`, conditional `--bdg --SPMR` on `save_macs_pileup`, gsize from khmer or config); treatment/control pairs map the upstream `ch_ip`/`ch_ip_control_bam` join — control-only samples are skipped via `optional = true`; gated on `narrow_peak=false`; pileup `.bdg` files moved to the macs3 dir when `save_macs_pileup=true` |
| FRIP_SCORE | `frip_score` | bedtools 2.30.0, samtools 1.17 | identical (intersectBed `-c -f 0.20`, flagstat `mapped (` non-primary fraction); gated on `narrow_peak=false` |
| MULTIQC_CUSTOM_PEAKS | `multiqc_custom_peaks` | bash, awk | identical (`wc -l` peak count, FRiP header concat); peak-count for control-only samples skipped like upstream; gated on `narrow_peak=false` |
| HOMER_ANNOTATEPEAKS | `homer_annotatepeaks` | homer 4.11 | identical (`-gid -gtf -cpu`); gated on `narrow_peak=false` |
| PLOT_MACS3_QC | `plot_macs3_qc` | r-base 3.5.1, macs3 3.0.1 | identical (`-i` comma paths, `-s` paths minus `_peaks.broadPeak`, `-o qc -p macs3_peak`); gated on `narrow_peak=false` |
| PLOT_HOMER_ANNOTATEPEAKS | `plot_homer_annotatepeaks` | r-base 3.5.1, homer 4.11 | identical (comma paths, summary + MQC header concat); gated on `narrow_peak=false` |
| BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2 → MACS3_CONSENSUS (local) | `macs3_consensus` | bedtools 2.30.0, macs3 3.0.1 | identical (mergeBed collapse, `macs3_merged_expand.py --min_replicates`, awk BED/SAF conversion, `plot_peak_intersect.r`, antibody.txt); gated on `narrow_peak=false` |
| BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2 → HOMER_ANNOTATEPEAKS (consensus) | `homer_annotate_consensus` | homer 4.11 | identical; gated on `narrow_peak=false` |
| ANNOTATE_BOOLEAN_PEAKS (local) | `annotate_boolean_peaks` | ubuntu 20.04 | identical (`cut -f2-`, sorted paste); gated on `narrow_peak=false` |
| SUBREAD_FEATURECOUNTS | `subread_featurecounts` | subread 2.0.1 | identical (`-F SAF -O --fracOverlap 0.2 -p -s 0`, counts IP-sample BAMs only); gated on `narrow_peak=false` |
| DESEQ2_QC (local) | `deseq2_qc` | mulled DESeq2 1.38.0 | identical (`--id_col 1 --sample_suffix .mLb.clN.sorted.bam --count_col 7 --vst TRUE`, header sed `_1` suffixes, mv to deseq2/); gated on `narrow_peak=false` |
| IGV (local) | `igv` | python 3.8.3 | identical (bigWig/Peak/bed find, consensus guard, antibody.txt, `igv_files_to_session.py ... --path_prefix '../../'`, genome.fa publish); gated on `narrow_peak=false` |
| MULTIQC (local) | `multiqc` | multiqc 1.23 | upstream mechanism replicated: `multiqc_config.yml` staged in cwd, `multiqc -f .`; config `path_filters`/`report_section_order` adapted to this port's `results/` layout; report, `multiqc_data/` and `multiqc_plots/` moved into `results/multiqc/broad_peak/`; gated on `narrow_peak=false` |
| MACS2_CALLPEAK / FRIP_SCORE / MULTIQC_CUSTOM_PEAKS / HOMER_ANNOTATEPEAKS / PLOT_MACS3_QC / PLOT_HOMER_ANNOTATEPEAKS / MACS3_CONSENSUS / HOMER_ANNOTATEPEAKS (consensus) / ANNOTATE_BOOLEAN_PEAKS / SUBREAD_FEATURECOUNTS / DESEQ2_QC / IGV / MULTIQC — narrow_peak mode | `*_narrow` rules (13) | same tools | identical commands with the upstream narrow_peak layout: MACS3 without `--broad/--broad-cutoff` (plus `*_summits.bed`), `narrowPeak` suffixes everywhere, `narrow_peak/` dirs (macs3, consensus, igv, multiqc), consensus merge of columns 2-10 (`collapse` x9) and `--is_narrow_peak`; all gated on `narrow_peak=true` |
| GTF2BED | `gtf2bed` | perl 5.26.2 | gated (`make_gene_bed`); runs the upstream `bin/gtf2bed` script verbatim; output fixed at `results/genome/gene.bed` (upstream names it after the GTF basename) |
| GENOME_BLACKLIST_REGIONS | `blacklist_regions` | bedtools 2.30.0 | gated (`make_blacklist_regions`); identical `sortBed \| complementBed` pipeline producing `chrom.sizes.include_regions.bed` |
| prepare_genome (gunzip/untar, custom getchromsizes, BWA/Bowtie2/Chromap index builders) / samplesheet_check | — | — | **not ported** — Nextflow plumbing / reference derivation; the port takes pre-built references (or builds gene.bed / include-regions / STAR index with the gated rules above) |
| UMI handling (UMITOOLS_EXTRACT, umi_extract) | — | — | **not ported** — `with_umi=false` is hardcoded in chipseq.nf at 2.1.0 (dead branch; the parameter does not exist) |
| save_align_intermeds / save_mapped / save_tracks outputs | — | — | intermediate BAM/bedGraph publish branches are `false` by default upstream; the port always behaves as `save_align_intermeds=true`; `save_macs_pileup` IS ported (conditional pileup publication in both macs3 rules) |
| DUMP_SOFTWARE_VERSIONS / pipeline summary + software versions sections of MultiQC | — | — | **not ported** — Nextflow metadata plumbing (paramsSummaryMap/softwareVersionsToYAML); `multiqc_data/` and `multiqc_plots/` ARE published |
| Multi-antibody consensus (`consensus_cluster` grouping) | — | — | single antibody (`config.antibody`) per run; upstream multi-antibody grouping is out of scope |
| Bowtie2 / Chromap aligner alternatives | — | — | **not ported** — only the BWA (default) and STAR branches are ported |

### Known divergences

- **Sample metadata**: nf-core reads a samplesheet (`--input`); oxo-flow uses
  `[[pairs]]` in `main.oxoflow` (pair_id, experiment, control). `ip_ids`
  (samples that receive peak calling) must be kept in sync with `[[pairs]]`.
  Upstream runs MACS3/FRiP/plotFingerprint only for samples that have a
  control; the port mirrors this exactly: per-pair rules whose `{control}`
  input is empty for control-only samples are skipped at run time
  (`optional = true`).
- **Reference inputs**: upstream derives references from `--genome`/iGenomes
  (GTF2BED for gene body regions, blacklist check); this port consumes
  pre-built files (`fasta`, `fai`, `gtf`, `gene_bed`, `chrom_sizes`,
  `blacklist`, `bwa_index` prefix).
- **MultiQC config**: `path_filters` and `report_section_order` were adapted
  to the port's `results/` layout and the single-antibody assumption;
  module order and custom-content sections are otherwise identical.
- **STAR output tree**: STAR alignment results stay under `results/bwa/`
  (bwa_mem/star_align both produce `{pair_id}.Lb.bam`, so the whole downstream
  chain is shared). Upstream also routes STAR BAMs through the same
  `results/bwa/` subworkflows — the difference is only that upstream STAR
  alignment writes `{id}.Aligned.out.bam`; the port renames it to `.Lb.bam`.
  STAR logs are published to `results/bwa/library/log/` (upstream feeds them
  to MultiQC only, from the workdir).
- **STAR index generation**: `star_genomegenerate` computes
  `--genomeSAindexNbases` from the provided `.fai` (upstream re-derives the
  index with `samtools faidx` inside its workdir); this port consumes the
  supplied index so the reference files are never touched. Index outputs go
  to a fixed `results/star/index/` directory — set
  `star_index = "results/star/index"` with `make_star_index = true`.
- **GTF2BED output name**: `gtf2bed` writes `results/genome/gene.bed`
  (upstream names it `{gtf.baseName}.bed`); point `gene_bed` at the generated
  file when using `make_gene_bed = true`.
- **No `--narrow-cutoff` at 2.1.0**: the upstream `narrow_peak` parameter has
  no cutoff argument in this release (it only switches file names, MACS3
  broad args and the consensus merge columns) — the port matches that
  behaviour exactly.
- **known limitation**: with `skip_consensus_peaks = true` the IGV/MultiQC
  consensus inputs (`consensus_peaks.bed`, featureCounts summary) are absent;
  keep the default (`false`) or set `skip_igv`/`skip_multiqc` together.

## Test

```bash
bash test/run.sh
```

Runs `validate`, `lint` and `dry-run` (plus a debug-expansion check) against
`main.oxoflow` and exits non-zero on any failure. CI runs it on every push.

## License

This workflow is licensed under the Apache License, Version 2.0 — see
[LICENSE](LICENSE). The port is derived from nf-core/chipseq (MIT); the
upstream license text is preserved verbatim in
[LICENSE.upstream](LICENSE.upstream).
