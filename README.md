# MGERT2

[![Unit Tests](https://github.com/andrewgull/MGERT2/actions/workflows/units.yml/badge.svg)](https://github.com/andrewgull/MGERT2/actions/workflows/units.yml)
[![Code Coverage](https://github.com/andrewgull/MGERT2/actions/workflows/test-and-coverage.yml/badge.svg)](https://github.com/andrewgull/MGERT2/actions/workflows/test-and-coverage.yml)
[![Integration test](https://github.com/andrewgull/MGERT2/actions/workflows/integration.yml/badge.svg)](https://github.com/andrewgull/MGERT2/actions/workflows/integration.yml)
[![Snakemake Dry Run](https://github.com/andrewgull/MGERT2/actions/workflows/snakemake-dry-run.yml/badge.svg)](https://github.com/andrewgull/MGERT2/actions/workflows/snakemake-dry-run.yml)
[![Snakefile Formatting Check](https://github.com/andrewgull/MGERT2/actions/workflows/snakefmt.yml/badge.svg)](https://github.com/andrewgull/MGERT2/actions/workflows/snakefmt.yml)
[![Python Code Style](https://github.com/andrewgull/MGERT2/actions/workflows/pycodestyle.yml/badge.svg)](https://github.com/andrewgull/MGERT2/actions/workflows/pycodestyle.yml)
[![Markdown Style](https://github.com/andrewgull/MGERT2/actions/workflows/markdown-lint.yml/badge.svg)](https://github.com/andrewgull/MGERT2/actions/workflows/markdown-lint.yml)
![CodeRabbit Pull Request Reviews](https://img.shields.io/coderabbit/prs/github/andrewgull/MGERT2?utm_source=oss&utm_medium=github&utm_campaign=andrewgull%2FMGERT2&labelColor=171717&color=FF570A&link=https%3A%2F%2Fcoderabbit.ai&label=CodeRabbit+Reviews)
![GitHub commit date](https://img.shields.io/github/last-commit/andrewgull/MGERT2?label=last%20commit&color=orange&style=flat)

## Introduction

This is a rework of [MGERT](https://github.com/andrewgull/MGERT) -
aka a big-ass 6-year-old Python script which was slain by "dependency hell".

---

## User Guide

### What the pipeline does

MGERT2 takes one or more genome assemblies and a target TE family name, then:

1. Builds a repeat database and runs **RepeatModeler2** to find TE families.
2. Collects consensus sequences matching the target TE name.
3. Runs **RepeatMasker** to locate all copies in the genome.
4. Extracts the genomic TE sequences as FASTA.
5. Calls ORFs, scores coding potential, and searches for protein domains
   via **RPS-BLAST**.
6. Classifies each ORF by confidence and produces a summary report.

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
| `data/{sample}.fasta.gz` | Gzip-compressed genome assembly |
| `data/pfam/*.smp` | PSSM files for the RPS-BLAST domain database |
| `data/domains.csv` | TSV mapping each `.smp` file to a domain type |

`data/domains.csv` is a two-column tab-separated file with no header:

```text
pfam00078.smp   RT
cd00304.smp     EN
```

### Configuration

Edit `config/config.yaml` before running.

#### Required settings

```yaml
samples:
  - "my_genome"   # basename of data/my_genome.fasta.gz

te_name: "Penelope"  # TE family name to find in RepeatModeler output
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

#### Domain filter

Controls which ORFs and TEs are flagged as domain-supported.

```yaml
orf_analysis:
  rpsblast:
    pssm_dir: "data/pfam"          # directory of *.smp PSSM files
    domains_csv: "data/domains.csv" # .smp → domain type mapping
    evalue: 1e-5
    min_query_coverage: 0.35
    min_bitscore: 50

  domain_filter:
    # Domain types (from domains.csv) that must be present.
    # Empty list = accept any domain hit.
    required_domains: []
    # How many of required_domains each ORF must match.
    min_domain_types_per_orf: 1
    # true = RT and EN may be in different ORFs of the same TE.
    allow_split_across_orfs: false
    # Which ORFs count toward TE-level coverage when split is on.
    # "all" = every ORF; "passing" = only ORFs that pass per-ORF filter.
    te_coverage_scope: "all"
```

Common scenarios:

- **Any domain hit is enough** — leave `required_domains: []`
- **Each ORF must have RT** — `required_domains: ["RT"]`
- **TE needs RT and EN, possibly in separate ORFs** —
  `required_domains: ["RT","EN"]`, `allow_split_across_orfs: true`
- **Same, but count only ORFs that individually pass** —
  add `te_coverage_scope: "passing"`

The filter adds two columns to the classification output:

- `passes_domain_filter` — per-ORF boolean
- `te_domain_complete` — per-TE boolean

#### Confidence classification

```yaml
orf_analysis:
  classification:
    high_min_aa: 300           # minimum length for high-confidence
    high_min_intrinsic: 0.60   # minimum intrinsic score
    putative_min_aa: 150
    putative_min_intrinsic: 0.50
```

Each ORF is labelled `high_confidence_coding`, `putative_coding`, or
`unlikely_coding`.

### Running

```bash
# Validate configuration without running jobs
pixi run dry-run

# Run the full pipeline (4 cores)
pixi run full-run
```

### Outputs

| Path | Contents |
| --- | --- |
| `results/plots/` | Kimura-distance repeat landscape plots |
| `results/extracted_te/` | Extracted TE copy sequences (FASTA) |
| `results/orf/*_orf_classification.tsv` | Per-ORF confidence and domain flags |
| `results/reports/*_orf_summary.tsv` | Counts by confidence class |
| `results/reports/*_orf_summary.html` | HTML version of the summary |

---
