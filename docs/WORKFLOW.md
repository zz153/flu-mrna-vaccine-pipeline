# Pipeline Workflow Guide

This document provides step-by-step instructions for running the complete flu vaccine design pipeline.

---

## Prerequisites

1. **Environment Setup:**
```bash
   conda env create -f env/environment.yml
   conda activate flu-vaccine-pipeline
```

2. **Data Acquisition:**
   - Follow instructions in [DATA_ACQUISITION.md](DATA_ACQUISITION.md)
   - Ensure raw FASTA files are in `data/raw/[Lineage]/`

---

## Pipeline Overview

The pipeline consists of the following steps:

1. **Data Cleaning & Quality Control**
2. **Sequence Alignment** (MAFFT)
3. **Phylogenetic Tree Building** (IQ-TREE)
4. **Vaccine Design** (Consensus, Medoid, ASR, COBRA)
5. **Distance Analysis** (p-distance and ML distance)
6. **Diversity Analysis** (Shannon entropy, sequence logos)
7. **Statistical Analysis & Visualization**

---

## Step-by-Step Execution

### Step 1: Clean and Filter Sequences

**Purpose:** Remove low-quality sequences, duplicates, and cluster at 99% identity

**For each lineage (H1N1, H3N2, VicB):**
```bash
# Example for H1N1
LINEAGE="H1N1"
YEAR_START=2009
YEAR_END=2025

# Combine regional files (if needed)
cat data/raw/${LINEAGE}/*.fasta > data/raw/${LINEAGE}/${LINEAGE}_all_regions.fasta

# Filter sequences (≥550 aa, no ambiguous characters, years 2009-2025)
seqkit seq -m 550 data/raw/${LINEAGE}/${LINEAGE}_all_regions.fasta | \
  seqkit grep -s -r -p "^[ACDEFGHIKLMNPQRSTVWY]+$" > \
  data/processed/${LINEAGE}/${LINEAGE}_filtered.fasta

# Remove duplicates
seqkit rmdup -s data/processed/${LINEAGE}/${LINEAGE}_filtered.fasta > \
  data/processed/${LINEAGE}/${LINEAGE}_unique.fasta

# Cluster at 99% identity using CD-HIT
cd-hit -i data/processed/${LINEAGE}/${LINEAGE}_unique.fasta \
  -o data/processed/${LINEAGE}/${LINEAGE}_cdhit99.fasta \
  -c 0.99 -n 5 -M 8000 -T 4
```

**Output:** `data/processed/[Lineage]/[Lineage]_cdhit99.fasta`

---

### Step 2: Multiple Sequence Alignment

**Purpose:** Align sequences using MAFFT
```bash
# Run alignment for each lineage
bash scripts/modules/run_yearly_alignments.sh H1N1 2009 2025

# Or manually:
mafft --auto --thread 4 \
  data/processed/${LINEAGE}/${LINEAGE}_cdhit99.fasta > \
  results/alignments/${LINEAGE}_aligned.fasta
```

**Output:** `results/alignments/[Lineage]_aligned.fasta`

---

### Step 3: Build Phylogenetic Trees

**Purpose:** Construct maximum-likelihood trees with IQ-TREE
```bash
# Submit SLURM job (on HPC)
sbatch scripts/modules/run_yearly_asr_generic.slurm H1N1 2009 2025

# Or run locally:
iqtree2 -s results/alignments/${LINEAGE}_aligned.fasta \
  -m LG+G4 -bb 1000 -nt AUTO \
  --prefix results/trees/${LINEAGE}
```

**Output:** 
- `results/trees/[Lineage].treefile`
- `results/trees/[Lineage].iqtree` (log file)

---

### Step 4: Create Vaccine Designs

**Purpose:** Generate candidate vaccine sequences using multiple strategies

#### 4a. Consensus Sequence
```bash
bash scripts/modules/run_yearly_consensus.sh H1N1 2009 2025
```

#### 4b. Medoid Sequence
```bash
sbatch scripts/modules/run_yearly_medoid.slurm H1N1 2009 2025
```

#### 4c. Ancestral Sequence Reconstruction (ASR)
```bash
sbatch scripts/modules/run_yearly_asr_generic.slurm H1N1 2009 2025
```

#### 4d. COBRA Design
```bash
sbatch scripts/modules/run_yearly_cobra.slurm H1N1 2009 2025
```

**Output:** `results/designs/[Lineage]_[Year]_[DesignType].fasta`

---

### Step 5: Calculate Evolutionary Distances

**Purpose:** Compute pairwise distances between vaccine designs and circulating strains
```bash
sbatch scripts/modules/run_yearly_distances.slurm H1N1 2009 2025
```

**Output:**
- `results/distances/[Lineage]_[Year]_distance_matrix.csv`
- `results/distances/[Lineage]_distance_summary.csv`

---

### Step 6: Diversity Analysis

**Purpose:** Analyze sequence variability and conservation patterns
```bash
# Shannon entropy calculation
python scripts/utils/calculate_entropy.py \
  --alignment results/alignments/${LINEAGE}_aligned.fasta \
  --output results/entropy/${LINEAGE}_entropy.csv

# Generate sequence logos
python scripts/utils/generate_sequence_logo.py \
  --alignment results/alignments/${LINEAGE}_aligned.fasta \
  --output results/figures/${LINEAGE}_logo.png
```

**Output:**
- `results/entropy/[Lineage]_entropy.csv`
- `results/figures/[Lineage]_logo.png`

---

### Step 7: Statistical Analysis & Visualization

**Purpose:** Generate summary statistics and publication-quality figures
```bash
# Comprehensive distance analysis
python scripts/analyze_distances.py \
  --lineage H1N1 \
  --output results/figures/

# Generate heatmaps
python scripts/plot_distances.py \
  --input results/distances/H1N1_distance_summary.csv \
  --output results/figures/H1N1_heatmap.png
```

**Output:**
- `results/figures/[Lineage]_heatmap.png`
- `results/tables/[Lineage]_summary_statistics.csv`

---

## Running the Complete Pipeline

### Option 1: Run All Steps Automatically (Coming Soon)
```bash
bash scripts/run_full_pipeline.sh --lineage H1N1 --start-year 2009 --end-year 2025
```

### Option 2: Run Steps Individually

Follow Steps 1-7 above sequentially for each lineage.

---

## Expected Runtime

On NeSI Aoraki HPC (using SLURM):
- **Data cleaning:** ~5-10 minutes per lineage
- **Alignment:** ~30-60 minutes per lineage
- **Tree building:** ~2-4 hours per lineage
- **Vaccine designs:** ~1-2 hours per lineage
- **Distance analysis:** ~30 minutes per lineage
- **Total:** ~6-10 hours per lineage (parallelizable)

---

## Troubleshooting

### Common Issues

1. **Memory errors during tree building:**
   - Increase memory allocation in SLURM script
   - Use `-mem 32G` or higher

2. **CD-HIT clustering fails:**
   - Reduce `-M` parameter (memory)
   - Split input file into smaller chunks

3. **MAFFT alignment too slow:**
   - Use `--auto` instead of `--maxiterate`
   - Consider `mafft-linsi` for very large datasets

---

## Citation

If you use this pipeline, please cite:
- MAFFT: Katoh & Standley (2013)
- IQ-TREE: Nguyen et al. (2015)
- CD-HIT: Fu et al. (2012)
- GISAID: Shu & McCauley (2017)
