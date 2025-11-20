⚙️ Environment Setup Guide

This document describes how to install all dependencies required to run the
FluHub: Vaccine Design and Evolutionary Analysis Pipeline.

The recommended setup uses a pinned Conda environment (`env/environment.yml`) to ensure that anyone can reproduce the results on both local systems and HPC clusters. A pip-only installation is also possible but is less strictly reproducible.

---

## 1. System Requirements

The pipeline has been tested on:

- Linux: Ubuntu 22.04, NeSI Aoraki HPC
- macOS: Apple Silicon (M2), Ventura / Sonoma

Minimum recommended specifications:

- Python: 3.11
- RAM: 8 GB
- CPU cores: 4

---

## 2. Core Command-Line Tools

The following external tools are required. Versions shown here match the pinned Conda environment in `env/environment.yml`.

Tool | Version used | Purpose
---- | ------------ | -------
MAFFT | 7.526 | Multiple sequence alignment per lineage and year
IQ-TREE | 3.0.1 | Maximum-likelihood phylogenetic tree inference
CD-HIT | 4.8.1 | Clustering of sequences at 99% identity
SeqKit | 2.10.1 | Sequence filtering, trimming, and ambiguity removal
Python | 3.11 | Core scripting environment

These tools are installed automatically if you use the Conda environment described below.

---

## 3. Recommended: Conda Environment (Reproducible)

This is the preferred way to set up the environment, especially on Aoraki and other Linux systems.

### 3.1 Create the environment

From the repository root:

```bash
conda env create -f env/environment.yml -n flu-vaccine-pipeline
conda activate flu-vaccine-pipeline
```

This will install:

- mafft=7.526
- iqtree=3.0.1
- cd-hit=4.8.1
- seqkit=2.10.1
- python=3.11
- and the Python packages listed below (with pinned versions where appropriate).

### 3.2 Verify installation

```bash
mafft --version
iqtree2 -h | head -n 1
cd-hit -h | head -n 1
seqkit version
python --version
```

If these commands run and show the expected or very similar versions, the environment is correctly set up.

---

## 4. Python Packages (Pinned Versions)

Within the Conda environment, the pipeline uses the following Python libraries:

Package | Version | Purpose
------- | ------- | -------
biopython | 1.86 | FASTA parsing, sequence manipulation, alignment I/O
pandas | 2.3.3 | Data wrangling, summary tables
numpy | 2.3.5 | Array operations, numerical computations
matplotlib | 3.10.8 | Plot generation for distances and entropy
seaborn | 0.13.2 | Statistical visualisation
scipy | 1.16.3 | Shannon entropy and numerical utilities
tqdm | (unpinned) | Progress bars for long-running computations

These versions match the current Conda environment (`conda list`) and are pinned in `env/environment.yml` (except `tqdm`, which is not version-pinned but not critical for numerical reproducibility).

To confirm inside the environment:

```bash
conda list biopython pandas numpy matplotlib seaborn scipy tqdm
```

---

## 5. Alternative: Pip-Only Installation (Not Recommended)

If you choose not to use Conda, you can install the Python dependencies with `pip` in a virtual environment. This is less reproducible across systems and assumes you have installed MAFFT, IQ-TREE, CD-HIT, and SeqKit separately via `apt`, `brew`, modules, etc.

### 5.1 Create and activate a virtual environment

```bash
python3.11 -m venv venv
source venv/bin/activate
```

### 5.2 Install Python packages

```bash
pip install \
  biopython==1.86 \
  pandas==2.3.3 \
  numpy==2.3.5 \
  matplotlib==3.10.8 \
  seaborn==0.13.2 \
  scipy==1.16.3 \
  tqdm
```

### 5.3 Verify Python package import

```bash
python -c "import Bio, pandas, numpy, matplotlib, seaborn, scipy, tqdm; print('✅ Python packages loaded successfully!')"
```

> Note: In pip-only mode, you are responsible for installing and version-matching `mafft`, `iqtree2`, `cd-hit`, and `seqkit` yourself.

---

## 6. Quick Environment “Smoke Test”

Once the environment is set up (Conda or pip), you can do a quick sanity check.

From the repository root:

```bash
# Activate the environment
conda activate flu-vaccine-pipeline  # or: source venv/bin/activate

# Python import test
python -c "import Bio, pandas, numpy, matplotlib, seaborn, scipy, tqdm; print('✅ Python environment OK')"

# Command-line tools
mafft --version
iqtree2 -h | head -n 1
cd-hit -h | head -n 1
seqkit version
```

If all of these succeed without errors, you are ready to run the pipeline.

---

## 7. Running the Pipeline (Full Workflow Overview)

The complete pipeline consists of multiple steps, implemented as shell scripts in `scripts/` and SLURM submission scripts in `scripts/slurm/`.

### 7.1 Serial (non-SLURM) example for one lineage (H1N1)

From the repository root:

```bash
# 1) Filter raw regional GISAID FASTA files
bash scripts/01_filter_sequences.sh H1N1

# 2) Cluster sequences at 99% identity (CD-HIT)
bash scripts/02_cluster_sequences.sh H1N1

# 3) Split clustered sequences by year (e.g. 2009–2025)
bash scripts/02b_split_by_year.sh H1N1 2009 2025

# 4) Align sequences per year (MAFFT)
#    The last argument is the number of threads, e.g. 4
bash scripts/03_align_sequences.sh H1N1 2009 2025 4

# 5) Per-year consensus / medoid / ancestral design
bash scripts/04_per_year_analysis.sh H1N1 2009 2025 4

# 6) Distance calculations to vaccine / medoid / ancestral designs
bash scripts/05_calculate_distances.sh H1N1 2009 2025
```

You would typically repeat the same pattern for `H3N2` and `VicB`.

### 7.2 SLURM / HPC usage examples (Aoraki)

On NeSI Aoraki, the same steps can be run via SLURM, using the submission scripts under `scripts/slurm/`. For example:

```bash
# Per-year analysis for a given lineage
sbatch --export=LINEAGE=H1N1 scripts/slurm/04_per_year_analysis.slurm

# Distance calculations for the same lineage
sbatch --export=LINEAGE=H1N1 scripts/slurm/05_calculate_distances.slurm
```

Refer to the comments inside each `scripts/slurm/*.slurm` file for the relevant SBATCH options (partition, time, CPUs, memory).

---

## 8. Where to Go Next

For a detailed description of each step, input/output layout, and expected files, see:

- `docs/WORKFLOW.md` – step-by-step description of the pipeline
- `docs/DATA_ACQUISITION.md` – how to obtain and organise raw HA sequences
- `README.md` – quick-start overview and high-level project description

This guide focuses solely on software and environment setup to ensure that anyone can reproduce the analysis with the same toolchain and library versions.


