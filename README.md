# MGERT2

[![Unit Tests](https://github.com/andrewgull/MGERT2/actions/workflows/units.yml/badge.svg)](https://github.com/andrewgull/MGERT2/actions/workflows/units.yml)
[![Code Coverage](https://github.com/andrewgull/MGERT2/actions/workflows/test-and-coverage.yml/badge.svg)](https://github.com/andrewgull/MGERT2/actions/workflows/test-and-coverage.yml)
[![Snakemake Dry Run](https://github.com/andrewgull/MGERT2/actions/workflows/snakemake-dry-run.yml/badge.svg)](https://github.com/andrewgull/MGERT2/actions/workflows/snakemake-dry-run.yml)
[![Python Code Style](https://github.com/andrewgull/MGERT2/actions/workflows/pycodestyle.yml/badge.svg)](https://github.com/andrewgull/MGERT2/actions/workflows/pycodestyle.yml)
![CodeRabbit Pull Request Reviews](https://img.shields.io/coderabbit/prs/github/andrewgull/MGERT2?utm_source=oss&utm_medium=github&utm_campaign=andrewgull%2FMGERT2&labelColor=171717&color=FF570A&link=https%3A%2F%2Fcoderabbit.ai&label=CodeRabbit+Reviews)
![GitHub commit date](https://img.shields.io/github/last-commit/andrewgull/MGERT2?label=last%20commit&color=orange&style=flat)

## Introduction

This is a rework of [MGERT](https://github.com/andrewgull/MGERT) —
a big-ass 6-year-old Python script slain by dependency hell.

---

## User Guide

### What the pipeline does

MGERT2 takes one or more genome assemblies and a target TE family name, then:

1. Builds a repeat database per genome and runs **RepeatModeler2** to
   discover TE families.
2. Collects consensus sequences matching the target TE name.
3. Runs **RepeatMasker** to locate all copies in each genome and plots
   a Kimura-distance repeat landscape.
4. Extracts the genomic TE copy sequences as FASTA.
5. Calls ORFs, filters them, and scores intrinsic coding potential.
6. Searches for protein domains via **RPS-BLAST** against a user-supplied
   PSSM database; annotates each hit with a domain type (e.g. RT, EN).
7. Classifies each ORF by confidence level and produces a summary report.
8. Pools ORF sequences across all genomes, aligns them with **MAFFT**,
   and infers a maximum-likelihood phylogenetic tree with **IQ-TREE2**
   — one tree per TE family.

### Installation

Install [Pixi](https://pixi.prefix.dev/latest/installation/), then
clone the repository and run:

```bash
pixi install
```

This creates a fully reproducible environment with all tools and packages.

### Input files

| Path | Description |
| --- | --- |
| Any genome assembly path | Gzip-compressed FASTA (`.fasta.gz`, `.fa.gz`) |
| `data/pfam/*.smp` | PSSM files compiled into the RPS-BLAST domain database |
| `data/domains.csv` | TSV mapping each `.smp` filename to a domain type label |

`data/domains.csv` is a two-column tab-separated file with no header:

```text
pfam00078.smp   RT
cd00304.smp     EN
```

Genome files can live anywhere — their paths are listed directly in
`config.yaml`. Sample names are derived automatically from the filename
by stripping compression and FASTA suffixes:
`/data/X.laevis.fasta.gz` → sample name `X.laevis`.

### Configuration

Edit `config/config.yaml` before running.

#### Genomes and TE name

```yaml
# List one or more genome assemblies; paths can be absolute or relative.
# Compression (.gz, .bz2, .xz) and FASTA (.fasta, .fa, .fna, .ffn, .fas)
# suffixes are stripped to derive the sample name.
genomes:
  - "data/organism1.fasta.gz"
  - "data/organism2.fna.gz"

te_name: "Penelope"  # TE family name to search for in RepeatModeler output
```

#### RepeatMasker

```yaml
repeatmasker:
  engine: "rmblast"    # or "crossmatch", "abblast"
  extra: " -s -nocut"  # additional RepeatMasker flags
```

#### ORF analysis

```yaml
orf_analysis:
  min_orf_aa: 100            # minimum ORF length (amino acids)
  require_start_codon: true  # false = include 5'-truncated ORFs
  max_orfs_per_te: 3         # keep at most N ORFs per TE copy
  require_stop_codon: true
```

#### Domain search

```yaml
orf_analysis:
  rpsblast:
    pssm_dir: "data/pfam"           # directory of *.smp PSSM files
    domains_csv: "data/domains.csv" # .smp filename → domain type mapping
    evalue: 1e-5
    max_target_seqs: 10
    min_query_coverage: 0.35
    min_bitscore: 50
```

#### Domain filter

Controls which ORFs and TEs are flagged as domain-supported in the
classification output.

```yaml
orf_analysis:
  domain_filter:
    # Domain types (from domains.csv) that must be present.
    # Empty list = accept any domain hit.
    required_domains: []
    # How many of required_domains an ORF must have to pass the per-ORF filter.
    min_domain_types_per_orf: 1
    # true = required domains may be spread across different ORFs of one TE.
    allow_split_across_orfs: false
    # Which ORFs count toward TE-level domain coverage when split is on.
    # "all" = every ORF; "passing" = only ORFs that individually pass.
    te_coverage_scope: "all"
```

Common scenarios:

- **Any domain hit is enough** — leave `required_domains: []`
- **Each ORF must have RT** — `required_domains: ["RT"]`
- **TE needs RT and EN, possibly in separate ORFs** —
  `required_domains: ["RT", "EN"]`, `allow_split_across_orfs: true`
- **Same, but only count ORFs that pass per-ORF filter** —
  add `te_coverage_scope: "passing"`

The filter adds two columns to the classification output:

- `passes_domain_filter` — per-ORF boolean
- `te_domain_complete` — per-TE boolean

#### Confidence classification

```yaml
orf_analysis:
  classification:
    high_min_aa: 300           # minimum length for high-confidence
    high_min_intrinsic: 0.60   # minimum intrinsic coding score
    putative_min_aa: 150
    putative_min_intrinsic: 0.50
```

Each ORF is labelled `high_confidence_coding`, `putative_coding`, or
`unlikely_coding`.

#### Phylogeny

```yaml
phylogeny:
  model: "MFP"    # ModelFinder Plus (auto); or a fixed model e.g. "LG+G4"
  bootstrap: 1000 # ultrafast bootstrap replicates (>=1000 for publication)
  seed: 12345     # fixed seed for reproducible trees
```

ORF sequences from all genomes are pooled per TE family, aligned with
MAFFT, and passed to IQ-TREE2. One tree is produced per `te_name`.

### Running

```bash
# Validate configuration without running any jobs
pixi run dry-run

# Run full analysis using conda envs
pixi run full-conda

# or using containers
pixi run full-container
```

### Outputs

| Path | Contents |
| --- | --- |
| `results/plots/` | Kimura-distance repeat landscape plots (per sample) |
| `results/extracted_te/` | Extracted TE copy sequences in FASTA format |
| `results/orf/*_orf_classification.tsv` | Per-ORF confidence and domain flags |
| `results/reports/*_orf_summary.tsv` | Counts by confidence class |
| `results/reports/*_orf_summary.html` | HTML version of the summary |
| `results/phylogeny/{te_name}.treefile` | Best ML tree in Newick format |
| `results/phylogeny/{te_name}.contree` | Consensus tree with bootstrap values |
| `results/phylogeny/{te_name}.iqtree` | Full IQ-TREE2 analysis report |

---
