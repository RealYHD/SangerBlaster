<p align="center">
  <a href="https://github.com/murphygargantua-agent/SangerBlaster#GPL-3.0-1-ov-file"><img src="https://img.shields.io/badge/License-GPL%203.0-yellow.svg" alt="License GPL 3.0" style="display: inline-block;"></a>
  <a href="https://github.com/asadprodhan/sangerFlow"><img src="https://img.shields.io/badge/Inspired%20by-sangerFlow-blue.svg" alt="Inspired by sangerFlow" style="display: inline-block;"></a>
</p>

<br />

## About SangerBlaster

SangerBlaster is inspired by [sangerFlow](https://github.com/asadprodhan/sangerFlow) (a Nextflow pipeline for Sanger amplicon sequencing analysis). It adds **native SCF chromatogram support** via [Tracy](https://github.com/gear-genomics/tracy) and replaces the consensus-forming steps with reference-guided or de-novo assembly of SCF traces.

Original pipeline: M. Asaduzzaman Prodhan, DPIRD Diagnostics and Laboratory Services, Department of Primary Industries and Regional Development, 3 Baron-Hay Court, South Perth, WA 6151, Australia. Correspondence: Asad.Prodhan@dpird.wa.gov.au

Original publication: Prodhan et al. (2024) "sangerFlow, a Sanger sequencing-based bioinformatics pipeline for pests and pathogens identification" (bioRxiv, DOI: 10.1101/2024.05.10.593518).

SangerBlaster was developed by Harrison Deng (murphygargantua-agent) with assistance from the Hermes Agent AI.

## AI Content Disclosure

This repository contains code generated or significantly modified by artificial intelligence (Hermes Agent). The original pipeline (sangerFlow) was written by M. Asaduzzaman Prodhan. The SangerBlaster modifications — including the Tracy integration, `main.nf` rewrite, config file definitions, and documentation — were produced with AI assistance.

## How to use SangerBlaster

### Step 1: Install required software

**Option A: Use Conda/Mamba with Bioconda (recommended for SLURM servers)**

```bash
# Install Mamba (faster conda) if not already installed
conda install -y mamba -c conda-forge

# Create and activate environment with Tracy and BLAST
mamba create -y -n sangerblaster -c conda-forge -c bioconda tracy blast
source activate sangerblaster

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

### Step 3: Download a BLAST database from NCBI

See: https://github.com/asadprodhan/blastn

**Quick test option:** Use the provided sample BLAST database for initial testing (the FASTA source is in this repo's `data/` directory):
```bash
# Build sample BLAST database (already done in data/):
makeblastdb -in data/ompA_test.fasta -dbtype nucl -out data/ompA_test
```

### Step 4: Run SangerBlaster (Quick Start Example)

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
nextflow run main.nf -c slurm_conda.config.example \
    --samples_tsv samples.tsv \
    --publish_dir ./test_output \
    --db /path/to/blastn_database/base_name
```

> **⚠️ CRITICAL:** The `--db` parameter is **required** for the BLAST workflow. Without it, the pipeline will exit immediately with an error.  
> The `--db` parameter must point to an **absolute path** to the base name of your BLAST database (the prefix of the index files), because blastn runs inside its own per-sample work directory and cannot resolve relative paths. For example:
> - If your files are `/path/to/db/ompA_test.ndb`, use `--db /path/to/db/ompA_test`
> - If your files are `/path/to/db/blastdb.nhd`, use `--db /path/to/db/blastdb`

3. Check results:
```bash
ls test_output/temp/  # Consensus FASTA files
# LHMD002.consensus.fasta
```

## Resource Definitions

The pipeline includes resource definitions (CPU, memory, time) for optimized SLURM execution. Use `slurm_conda.config.example` for the appropriate configuration.

## Automated Testing

This repository includes a GitHub Actions workflow (`.github/workflows/test_pipeline.yml`) that runs an end-to-end test on every push to `main`:

1. Sets up Miniconda with Bioconda channels
2. Creates a Conda environment with `tracy` and `blast`
3. Builds a test BLAST database from `data/ompA_test.fasta`
4. Runs the full pipeline with sample SCF data (LHMD002, LHMD004, LHMD014)
5. Verifies all output files are generated correctly
6. Tests duplicate sample IDs in `samples.tsv`:
   - Duplicates without `--dedupe_samples` → fail fast
   - Duplicates with identical reads + `--dedupe_samples` → collapse to one run
   - Duplicates with different reads + `--dedupe_samples` → fail fast
7. Tests a single-sample run (the original CI baseline)

This ensures the pipeline works correctly in CI before deployment.

## Container Fallback (Optional)

If Conda is not available on your system, a Singularity container can be used instead. The script `build_sangerblaster_container.sh` builds a Singularity image with Tracy pre-installed:

```bash
./build_sangerblaster_container.sh
```

However, the Conda method is strongly recommended because:
- It's faster to set up and update
- No container runtime (Singularity/Docker) required on SLURM clusters
- Bioconda provides automatic dependency management

See `docker_local.config.example` for a Docker-based configuration if needed.

## Outputs

When the run completes, three directories are created:

- **results/** — BLAST hit TSV files per sample + master `concatenatedHits_withHeaders.tsv`
  - The concatenated hits file has a `qaccver` column that carries the **sample ID** (the consensus FASTA header is renamed to the sample ID by `tracyAssemble`), so every hit can be traced back to the sample it came from.
- **temp/** — intermediate consensus FASTA files
- **work/** — Nextflow work directory (can be deleted to free space)

### Graceful failure handling

The pipeline keeps running even when individual samples fail:

- **Tracy assembly** uses `errorStrategy 'ignore'`, so a failing sample only exits with Tracy's failure code (recorded in `work/<hash>/.exitcode`) and the rest of the run continues. Nextflow logs a `NOTE:` line for each ignored sample, which is normal.
- **BLAST** also uses `errorStrategy 'ignore'` with optional outputs, so a failed consensus simply produces no hits downstream instead of failing the whole run. `blastnTopHits` behaves the same way: if a consensus generates zero hits, no topHits file is emitted and that sample is skipped by the master concatenation.

### Pipeline-level hard failures

These conditions stop the run immediately — they are validated in the workflow before any process runs:

- **Missing BLAST database** (`--db` not found) — validation stops the run with a clear error rather than failing mid-run.
- **Duplicate sample IDs** in `samples.tsv` — each `sampleId` must be unique, otherwise every downstream output file name collides with itself. Pass `--dedupe_samples` to allow duplicates when all forward and reverse reads are identical (only one copy is processed); it still fails if the reads differ across duplicates.
- **No data at all** — if *no* sample produces BLAST hits (every assembly failed), `blastHits_concat` exits with status 1, because a run that produced nothing is an error, not a result. A single-sample or mixed run (some samples fail, some succeed) is fine: the hits from the successful samples are concatenated as usual.

## Limitations

- Reference-guided assembly requires `--ref` pointing to a valid reference FASTA.
- De-novo assembly (without `--ref`) is supported but may produce lower-quality consensus when traces don't sufficiently overlap.

## Differences from sangerFlow

| Feature | sangerFlow (original) | SangerBlaster (this pipeline) |
||---|---|---|
|| **Input** | FASTA files (already basecalled) | **SCF chromatogram traces** (forward + reverse) |
|| **Consensus** | bioawk → sed → Clustal → EMBOSS cons | **Tracy `assemble`** — reference-guided assembly of SCF traces |
|| **Reference** | Not required | **Optional** (`--ref` for reference-guided; de-novo by default) |
|| **QC** | Simple sed-based ambiguity trimming | Tracy-based quality-controlled consensus |
|| **BLAST** | Yes (unchanged) | Yes (unchanged) |

SangerBlaster retains the original sangerFlow's BLAST-based species identification workflow, but replaces the consensus-forming steps (bioawk/Clustal/EMBOSS) with Tracy's assembly, which:
- Reads SCF/AB1/ABI chromatogram files natively
- Assembles forward and reverse traces — reference-guided if `--ref` is given, or de-novo otherwise
- Produces a consensus FASTA with uncovered positions left as N/ambiguous (no forced filling)
- Tracy is installed via Bioconda (recommended)

## License

This project is licensed under GPL v3.0. The original sangerFlow was also licensed under GPL v3.0.
