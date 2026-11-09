<h1 align="center">GenoSanger — Sanger amplicon sequencing analysis with Tracy</h1>

<p align="center">
  <a href="https://github.com/murphygargantua-agent/GenoSanger#GPL-3.0-1-ov-file"><img src="https://img.shields.io/badge/License-GPL%203.0-yellow.svg" alt="License GPL 3.0" style="display: inline-block;"></a>
  <a href="https://github.com/asadprodhan/sangerFlow"><img src="https://img.shields.io/badge/Original-sangerFlow-blue.svg" alt="Based on sangerFlow" style="display: inline-block;"></a>
</p>

<br />

## About GenoSanger

GenoSanger is a fork of [sangerFlow](https://github.com/asadprodhan/sangerFlow) (a Nextflow pipeline for Sanger amplicon sequencing analysis) that adds **native SCF chromatogram support** via [Tracy](https://github.com/gear-genomics/tracy).

Original pipeline: M. Asaduzzaman Prodhan, DPIRD Diagnostics and Laboratory Services, Department of Primary Industries and Regional Development, 3 Baron-Hay Court, South Perth, WA 6151, Australia. Correspondence: Asad.Prodhan@dpird.wa.gov.au

Original publication: Prodhan et al. (2024) "sangerFlow, a Sanger sequencing-based bioinformatics pipeline for pests and pathogens identification" (bioRxiv, DOI: 10.1101/2024.05.10.593518).

GenoSanger was developed by Harrison Deng (murphygargantua-agent) with assistance from the Hermes Agent AI.

### Key differences from sangerFlow

| Feature | sangerFlow (original) | GenoSanger (this fork) |
|---|---|---|
| **Input** | FASTA files (already basecalled) | **SCF chromatogram traces** (forward + reverse) |
| **Consensus** | bioawk → sed → Clustal → EMBOSS cons | **Tracy `assemble`** — reference-guided assembly of SCF traces |
| **Reference** | Not required | **Optional** (`--ref` for reference-guided; de-novo by default) |
| **QC** | Simple sed-based ambiguity trimming | Tracy-based quality-controlled consensus |
| **BLAST** | Yes (unchanged) | Yes (unchanged) |

GenoSanger retains the original sangerFlow's BLAST-based species identification workflow, but replaces the consensus-forming steps (bioawk/Clustal/EMBOSS) with Tracy's assembly, which:
- Reads SCF/AB1/ABI chromatogram files natively
- Assembles forward and reverse traces — reference-guided if `--ref` is given, or de-novo otherwise
- Produces a consensus FASTA with uncovered positions left as N/ambiguous (no forced filling)
- Auto-downloads Tracy v0.9.1 static binary on first run

## AI Content Disclosure

This repository contains code generated or significantly modified by artificial intelligence (Hermes Agent). The original pipeline (sangerFlow) was written by M. Asaduzzaman Prodhan. The GenoSanger modifications — including the Tracy integration, `main.nf` rewrite, `nextflow.config` resource definitions, and documentation — were produced with AI assistance.

## How to use GenoSanger

### Step 1: Install required software

- Install [Nextflow](https://www.nextflow.io/)
- Install [Singularity](https://docs.sylabs.io/) (or Docker)
- Optionally have a **reference FASTA** for your target locus (e.g., ompA reference for Chlamydia) — if not provided, Tracy uses de-novo assembly

> Note: Tracy is auto-downloaded from GitHub on first run. No separate Tracy installation needed.

### Step 1b: Verify system requirements (recommended)

Before running the pipeline, verify all requirements are met:

```bash
./check_requirements.sh
```

This script checks for:
- ✓ Container runtime (Docker or Singularity)
- ✓ Download tools (wget or curl)
- ✓ Nextflow version

If any requirement is missing, the script exits with a clear error message.

### Step 2: Prepare a sample description file

Create a TSV file with three columns:
1. Sample ID
2. Forward SCF trace file name (e.g., `sample_R1.scf`)
3. Reverse SCF trace file name (e.g., `sample_R2.scf`)

Example `samples.tsv`:
```
sampleId	read1	read2
SAMPLE001	SAMPLE001_R1.scf	SAMPLE001_R2.scf
SAMPLE002	SAMPLE002_R1.scf	SAMPLE002_R2.scf
```

### Step 3: Download a BLAST database from NCBI

See: https://github.com/asadprodhan/blastn

### Step 3: Run GenoSanger

**Current status:** The pipeline is tested and working with Tracy assembly. BLAST steps are commented out until a test database is available.

You'll need a configuration file that matches your environment:

**For Docker (local development / cloud):**
```bash
nextflow run . -c docker_local.config \
    --samples_tsv samples.tsv \
    [--ref /path/to/reference.fa] \
    --db /path/to/blastn_database \
    --publish_dir ./
```

**For Singularity (SLURM clusters):**
```bash
nextflow run . -c singularity_nextflow.config \
    --samples_tsv samples.tsv \
    [--ref /path/to/reference.fa] \
    --db /path/to/blastn_database \
    --publish_dir ./
```

**For SLURM with Docker (production):**
```bash
nextflow run . -c prod_nextflow.config \
    --samples_tsv samples.tsv \
    [--ref /path/to/reference.fa] \
    --db /path/to/blastn_database \
    --publish_dir ./
```

You can set blastn thresholds with:
- `--evalue` — default 0.1
- `--cpus` — default 18
- `--topHits` — default 5

Example (reference-guided):
```bash
nextflow run . \
    --samples_tsv samples.tsv \
    --ref reference.fa \
    --evalue 0.05 --topHits 1 --cpus 16 \
    --db /path/to/blastn_database
```

Example (de-novo, no reference):
```bash
nextflow run . \
    --samples_tsv samples.tsv \
    --db /path/to/blastn_database \
    --publish_dir ./
```

> **Note:** The pipeline currently runs with hardcoded test samples. Full `samples.tsv` parsing is implemented but requires manual invocation. For testing, modify `main.nf` to use the sample channel directly.

### Step 3b: Verify system requirements (recommended)

Before running the pipeline, verify all requirements are met:

```bash
./check_requirements.sh
```

This script checks for:
- ✓ Container runtime (Docker or Singularity)
- ✓ Download tools (wget or curl)
- ✓ Nextflow version

If any requirement is missing, the script exits with a clear error message.

## Outputs

When the run completes, three directories are created:

- **results/** — BLAST hit TSV files per sample + master `concatenatedHits_withHeaders.tsv`
- **temp/** — intermediate consensus FASTA files
- **work/** — Nextflow work directory (can be deleted to free space)

## Limitations

- The Tracy binary must be downloadable from GitHub at runtime (requires network access).
- Reference-guided assembly requires `--ref` pointing to a valid reference FASTA.
- De-novo assembly (without `--ref`) is supported but may produce lower-quality consensus when traces don't sufficiently overlap.

## License

This project is licensed under GPL v3.0. The original sangerFlow was also licensed under GPL v3.0.
