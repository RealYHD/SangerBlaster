# GenoSanger — Sanger amplicon sequencing analysis with Tracy

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
- Tracy is installed via Bioconda (recommended)

## AI Content Disclosure

This repository contains code generated or significantly modified by artificial intelligence (Hermes Agent). The original pipeline (sangerFlow) was written by M. Asaduzzaman Prodhan. The GenoSanger modifications — including the Tracy integration, `main.nf` rewrite, config file definitions, and documentation — were produced with AI assistance.

## How to use GenoSanger

### Step 1: Install required software

**Option A: Use Conda/Mamba with Bioconda (recommended for SLURM servers)**

```bash
# Install Mamba (faster conda) if not already installed
conda install -y mamba -c conda-forge

# Create and activate environment with Tracy and BLAST
mamba create -y -n genosanger -c conda-forge -c bioconda tracy blast
source activate genosanger

# Verify installation
tracy --version  # Should show v0.9.1+ or newer
```

> **Important for SLURM servers:** You must have `conda` or `mamba` available in your environment PATH when running the pipeline. Nextflow will automatically create and use a Conda environment containing Tracy via the `conda 'bioconda::tracy'` directive in `main.nf`. If you encounter "tracy: command not found" errors, ensure conda is installed and accessible before running the pipeline.

### Step 1b: Verify system requirements (recommended)

Before running the pipeline, verify all requirements are met:

```bash
./check_requirements.sh
```

This script checks for:
- ✓ Conda/Mamba with `tracy` package installed from Bioconda
- ✓ Nextflow version
- ✓ BLAST database availability (if specified)

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

### Step 2: Download a BLAST database from NCBI

See: https://github.com/asadprodhan/blastn

**Quick test option:** Use the provided sample BLAST database for initial testing (created in this repo's `test_db/` directory):
```bash
# Build sample BLAST database (already done in test_db/)
# If you need to rebuild:
makeblastdb -in test_db/ompA_test.fasta -dbtype nucl -out test_db/ompA_test
```

### Step 3: Run GenoSanger (Quick Start Example)

**Prerequisites:**
- SCF/AB1 chromatogram files (forward and reverse pairs)
- A BLAST database (download from NCBI or use `test_db/`)
- Conda environment with `tracy` and `blast` installed from Bioconda

**Example using test data in this repository:**

1. Create a `samples.tsv` file:
```bash
cat > samples.tsv << 'EOF'
sampleId	read1	read2
LHMD002	./SangerSequencingChromatograms/LHMD002-07_F1_1.ab1	./SangerSequencingChromatograms/LHMD002-07_R1_1.ab1
EOF
```

2. Run the pipeline (Conda method tested E2E on local Linux):
```bash
nextflow run main.nf -c slurm_conda.config \
    --samples_tsv samples.tsv \
    --publish_dir ./test_output \
    --db /path/to/blastn_database/base_name
```

> **Note:** The `--db` parameter should point to the **base name** of your BLAST database (the prefix of the `.ndb`, `.nhr`, etc. files). For example, if your BLAST database files are located at `/path/to/db/ompA_test.ndb`, use `--db /path/to/db/ompA_test`.

3. Check results:
```bash
ls test_output/temp/  # Consensus FASTA files
# LHMD002.consensus.fasta
```

## Resource Definitions

The pipeline includes resource definitions (CPU, memory, time) for optimized SLURM execution. Use `slurm_conda.config` for the appropriate configuration.

## Outputs

When the run completes, three directories are created:

- **results/** — BLAST hit TSV files per sample + master `concatenatedHits_withHeaders.tsv`
- **temp/** — intermediate consensus FASTA files
- **work/** — Nextflow work directory (can be deleted to free space)

## Limitations

- Reference-guided assembly requires `--ref` pointing to a valid reference FASTA.
- De-novo assembly (without `--ref`) is supported but may produce lower-quality consensus when traces don't sufficiently overlap.

## License

This project is licensed under GPL v3.0. The original sangerFlow was also licensed under GPL v3.0.
