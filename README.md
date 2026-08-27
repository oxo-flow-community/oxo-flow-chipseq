# oxo-flow-chipseq — ChIP-seq peak calling, QC and differential analysis

> ★ Verified · ⇄ Official port of [`nf-core/chipseq`](https://github.com/nf-core/chipseq) @ `2.1.0` — same tools, same versions, same commands. Part of the [oxo-flow-community catalog](https://oxo-flow-community.github.io/).

[![CI](https://github.com/oxo-flow-community/oxo-flow-chipseq/actions/workflows/ci.yml/badge.svg)](https://github.com/oxo-flow-community/oxo-flow-chipseq/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

ChIP-seq analysis pipeline: raw-read QC (FastQC), adapter trimming (Trim
Galore), BWA-MEM (default), Bowtie2, Chromap or STAR alignment, library
merging, Picard mark-duplicates, BAMTools filtering with blacklist and
orphan-read removal,
library complexity (preseq, phantompeakqualtools SPP), bigWig tracks, deepTools
QC plots, MACS3 peak calling with input controls in broad (default) or narrow
mode, HOMER peak annotation, FRiP scoring, consensus peaks across replicates
(MACS3 merge + featureCounts quantification + DESeq2 QC), an IGV session and a
MultiQC report. Runs paired-end samples; `aligner` and `narrow_peak` select the
alignment/peak branches, and an optional multi-antibody mode runs the
consensus chain once per antibody (see
[Multi-antibody runs](#multi-antibody-runs)).

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

- genome FASTA (`fasta`; the FASTA index `fai` is only needed for
  `make_star_index`)
- annotation GTF (`gtf`)
- gene-body regions in BED format (`gene_bed` — upstream derives this from
  the GTF; this port takes it pre-built, or generates it with
  `make_gene_bed = true`)
- chromosome sizes file (`chrom_sizes`, or derive it from the FASTA with
  `make_chrom_sizes = true`)
- blacklist regions in BED format (`blacklist`)
- aligner index — for the default `aligner = "bwa"` a BWA index directory
  (`bwa_index`) containing the BWA index files (`*.amb`, `*.ann`, `*.bwt`,
  `*.pac`, `*.sa`) for the reference FASTA; for `aligner = "bowtie2"` a
  Bowtie2 index directory (`bowtie2_index`), for `aligner = "chromap"` a
  Chromap index file (`chromap_index`, upstream default `genome.index`), and
  for `aligner = "star"` a STAR index directory (`star_index`). Each index
  can be built from the FASTA with the matching gated rule
  (`make_bwa_index` / `make_bowtie2_index` / `make_chromap_index` /
  `make_star_index`), then point the config key at the generated
  `results/...` path
- raw paired-end FASTQ reads (`raw_dir`, named `raw/{pair_id}_R{1,2}.fastq.gz`)
  with sample metadata declared in `[[pairs]]` and `ip_ids`

**Compute** — the largest rules (`bwa_mem`, `bowtie2_align`, `star_align`,
`trimgalore`, `star_genomegenerate`) request 12 threads and 72 GB of memory,
`bowtie2_index_build` 12 threads and 36 GB; most other rules request 6
threads and 36 GB (`chromap_align`, `chromap_index_build` included — see
`modules/*.oxoflow` `[rules.resources]`). A default single-rule run peaks
around 72 GB — scale the numbers down for small machines.

**Tool delivery** — containers with pinned images. Every rule pins its
biocontainers (or nf-core / Seqera Wave) image and the workflow runs with the
**docker** backend, so you need Docker (or Singularity) installed.

## Usage

```bash
oxo-flow validate main.oxoflow
oxo-flow dry-run main.oxoflow       # prints the 198-instance plan (153 run, 45 gated off)
oxo-flow run main.oxoflow           # executes locally (docker backend)
oxo-flow debug main.oxoflow         # show expanded commands
```

Branch selection (each key toggles only its branch; the default plan is
unchanged when all are off):

| Key | Effect |
|---|---|
| `aligner = "bowtie2"` | Bowtie2 alignment instead of BWA; requires `bowtie2_index` (pre-built Bowtie2 index directory) |
| `aligner = "chromap"` | Chromap alignment instead of BWA; requires `chromap_index` (pre-built Chromap index file) |
| `aligner = "star"` | STAR alignment instead of BWA; requires `star_index` (or `make_star_index = true` + `star_index = "results/star/index"`) |
| `narrow_peak = true` | whole peak chain in narrowPeak mode (MACS3 args, `narrow_peak/` dirs, consensus merge columns, IGV/MultiQC outputs) |
| `make_gene_bed = true` | derive `results/genome/gene.bed` from the GTF (then point `gene_bed` at it) |
| `make_blacklist_regions = true` | derive `results/genome/chrom.sizes.include_regions.bed` from blacklist + chrom sizes (then point `blacklist` at it) |
| `make_chrom_sizes = true` | derive `results/genome/chrom.sizes` + `results/genome/genome.fa.fai` from the FASTA (then point `chrom_sizes`/`fai` at them) |
| `make_bwa_index = true` | build the BWA index into `results/bwa/index` (requires `aligner = "bwa"`; then point `bwa_index` at it) |
| `make_bowtie2_index = true` | build the Bowtie2 index into `results/bowtie2/index` (requires `aligner = "bowtie2"`; then point `bowtie2_index` at it) |
| `make_chromap_index = true` | build the Chromap index into `results/chromap/index/genome.index` (requires `aligner = "chromap"`; then point `chromap_index` at it) |

The workflow runs with the **docker** backend (all tools pinned to the
upstream container images). Replace the `test/fixtures` input paths in
`main.oxoflow` `[config]` with your own reads and references, and edit
`[[pairs]]`/`ip_ids` (plus the `metadata_file` sample sheet for
multi-antibody runs) to match your samplesheet. Results are written to
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
| BOWTIE2_ALIGN | `bowtie2_align` | bowtie2 2.5.2, samtools 1.18 | gated (`aligner='bowtie2'`); identical command (`find`-located index with `.bt2l` fallback, paired-end `-1/-2`, `--threads`, RG flags `--rg-id/--rg` with `SM` minus `_T\d+` + conditional `CN`, `sort_bam=false`/`save_unaligned=false` so `samtools view` emits the BAM directly); stderr teed to `{pair_id}.Lb.bowtie2.log`, published under `results/bwa/library/log/` for MultiQC auto-discovery (upstream drops the log — no MultiQC bowtie2 section there) |
| CHROMAP_CHROMAP | `chromap_align` | chromap 0.2.6, samtools 1.20 | gated (`aligner='chromap'`); identical command (`-l 2000 --low-mem --SAM -t -x -r -1/-2`, then `samtools addreplacerg -r '@RG...'` + `samtools view -bh`); barcodes/whitelist/chr-order inputs are empty upstream (`[]`), so no such flags |
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
| CUSTOM_GETCHROMSIZES | `getchromsizes` | samtools 1.20 | gated (`make_chrom_sizes`); identical script (`samtools faidx` + `cut -f 1,2`); the fasta is symlinked into `results/genome/` so the user's reference files are never touched; outputs fixed at `results/genome/chrom.sizes` + `results/genome/genome.fa.fai` (upstream names them `{fasta}.sizes` / `{fasta}.fai`) |
| BWA_INDEX | `bwa_index_build` | bwa 0.7.18 | gated (`aligner='bwa' && make_bwa_index`); identical command (`bwa index -p {prefix} {fasta}`); prefix fixed at `results/bwa/index/genome` (upstream names the files after the fasta basename — `bwa_mem` locates the index by `find`, so the prefix is transparent) |
| BOWTIE2_BUILD | `bowtie2_index_build` | bowtie2 2.5.2 | gated (`aligner='bowtie2' && make_bowtie2_index`); identical command (`bowtie2-build --threads`); index base fixed at `results/bowtie2/index/genome` (upstream names the files after the fasta basename — `bowtie2_align` locates the index by `find`, so the base is transparent) |
| CHROMAP_INDEX | `chromap_index_build` | chromap 0.2.6 | gated (`aligner='chromap' && make_chromap_index`); identical command (`chromap -i -t -r -o`); output fixed at `results/chromap/index/genome.index` (upstream names it `{fasta.baseName}.index`) |
| prepare_genome (gunzip/untar, GFFREAD) / samplesheet_check | — | — | **not ported** — compressed-reference gunzip/untar convenience, GFFREAD (GFF3 -> GTF) and samplesheet validation/staging (the port consumes pre-built plain reference files; the samplesheet analogue is `[[pairs]]`). The reference-derivation steps that make sense on plain files ARE ported as gated rules — getchromsizes and the BWA/Bowtie2/Chromap/STAR index builders (rows above) |
| UMI handling (UMITOOLS_EXTRACT, umi_extract) | — | — | **not ported** — `with_umi=false` is hardcoded in chipseq.nf at 2.1.0 (dead branch; the parameter does not exist) |
| save_reference / save_trimmed / save_unaligned outputs | — | — | publish branches that are `false` by default upstream; the port always behaves as `save_align_intermeds=true` (intermediates are kept). Upstream 2.1.0 has no `save_mapped` / `save_tracks` params. `save_macs_pileup` IS ported (conditional pileup publication in both macs3 rules) |
| DUMP_SOFTWARE_VERSIONS / pipeline summary + software versions sections of MultiQC | — | — | **not ported** — Nextflow metadata plumbing (paramsSummaryMap/softwareVersionsToYAML); `multiqc_data/` and `multiqc_plots/` ARE published |
| Multi-antibody consensus (`consensus_cluster` grouping) | `macs3_consensus_multi` / `homer_annotate_consensus_multi` / `annotate_boolean_peaks_multi` / `subread_featurecounts_multi` / `deseq2_qc_multi` / `igv_multi` / `multiqc_multi` (+ `*_narrow_multi` variants) | same tools | upstream groups the consensus chain by `meta.antibody` (groupTuple `by: antibody`); the port does the same with the engine's metadata binding — `[workflow] metadata_file` (TSV: sample + antibody columns) + `input_groups` `group_by = "meta.antibody"` runs the consensus chain once per distinct antibody with per-antibody inputs, and the IGV/MultiQC rules collect every antibody. Gated on `config.multi_antibody` (default `false` — the single-antibody `config.antibody` path is byte-identical); see [Multi-antibody runs](#multi-antibody-runs) |

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
  pre-built files (`fasta`, `gtf`, `gene_bed`, `chrom_sizes`, `blacklist`,
  index prefixes). The derivation steps that make sense on plain reference
  files are ported as gated rules — `make_gene_bed`, `make_blacklist_regions`,
  `make_chrom_sizes`, `make_bwa_index`, `make_bowtie2_index`,
  `make_chromap_index`, `make_star_index` (see the Fidelity table).
- **Fixed output names for derived references**: the gated reference rules
  write fixed names — `results/genome/gene.bed`,
  `results/genome/chrom.sizes`, `results/genome/genome.fa.fai`,
  `results/bwa/index/genome.*`, `results/bowtie2/index/genome.*.bt2`,
  `results/chromap/index/genome.index`, `results/star/index/` — where
  upstream names outputs after the FASTA/GTF basename. The consuming rules
  locate indexes by `find`, so the prefix is transparent; point the config
  keys at the generated paths when the rules are enabled.
- **Multi-library merge not ported**: upstream groups libraries by
  `meta.id` minus `_T\d+` and merges variable-size library sets before
  markduplicates. The merge itself is expressible in oxo-flow
  (`expand_inputs` glob over the libraries), but every downstream rule keys
  on `pair_id`; a second per-sample grouping dimension would restructure the
  whole pipeline, so only the single-library default path (the upstream
  symlink shortcut, replicated exactly) is ported.
- **MultiQC config**: `path_filters` and `report_section_order` were adapted
  to the port's `results/` layout and the `antibody.txt` convention; the
  `multiqc_multi` variant scans the per-antibody `consensus/{antibody}/`
  trees instead of relying on antibody-specific inputs. Module order and
  custom-content sections are otherwise identical.
- **Aligner output tree**: every aligner's results stay under `results/bwa/`
  (bwa_mem/star_align/bowtie2_align/chromap_align all produce
  `results/bwa/library/{pair_id}.Lb.bam`, so the whole downstream chain is
  shared — upstream writes `results/{aligner}/library/` instead). Upstream
  also routes STAR BAMs through the same `results/bwa/` subworkflows — the
  difference is only that upstream STAR alignment writes
  `{id}.Aligned.out.bam`; the port renames it to `.Lb.bam`. STAR logs are
  published to `results/bwa/library/log/` (upstream feeds them to MultiQC
  only, from the workdir); the bowtie2 log (`{pair_id}.Lb.bowtie2.log`) is
  published to the same directory for MultiQC auto-discovery (upstream does
  not feed it to MultiQC at all).
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

## Multi-antibody runs

Upstream groups the consensus chain (consensus peaks, annotation,
featureCounts quantification, DESeq2 QC) by the samplesheet's `antibody`
column. The port replicates that with the engine's metadata binding:

- `[workflow] metadata_file` — a TSV whose first column is the sample id
  (matching `[[pairs]]` `pair_id`) plus an `antibody` column; samples
  without an antibody (controls) are excluded from the groups. The shipped
  `test/fixtures/samples.tsv` carries two antibody replicate sets
  (S1* = H3K4me3, S2* = H3K27ac) with their controls (C1*/C2*).
- `config.multi_antibody = true` — runs `macs3_consensus_multi`,
  `homer_annotate_consensus_multi`, `annotate_boolean_peaks_multi`,
  `subread_featurecounts_multi` and `deseq2_qc_multi` once per distinct
  antibody (each with its own BAMs and peaks), plus `igv_multi` /
  `multiqc_multi` collecting every antibody. With the default `false` the
  single-antibody `config.antibody` path runs, byte-identical to a
  workflow without a metadata table.

The easiest way to use it is the shipped profile, which flips the toggle
and replaces the four-sample fixture lists with the full eight-sample set
(`profiles/multi_antibody.toml`; the merge mode comes from
`profile_mode = "override"` in `main.oxoflow` — see
[the oxo-flow profile docs](https://oxo-flow.readthedocs.io/en/latest/commands/run/#profiles)):

```bash
oxo-flow run main.oxoflow --profile multi_antibody
```

For your own data, keep `metadata_file`'s first column in sync with
`[[pairs]]`, set `multi_antibody = true` (plus your full `pair_ids` /
`ip_ids` lists), and declare the second antibody's pairs with
`when = "config.multi_antibody"` so they only fan out in multi mode.
When the metadata table has fewer than two distinct antibodies with
replicates, the multi rules fall back to the same per-sample behaviour
as the single-antibody path.

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
