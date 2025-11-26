# Pipeline Workflow Guide

This document describes the complete flu vaccine design pipeline workflow.

---

## Quick Start
```bash
# Clone repository
git clone https://github.com/zz153/flu-vaccine-design-distance-analysis.git
cd flu-vaccine-design-distance-analysis

# Create conda environment
conda env create -f env/environment.yml
conda activate flu-vaccine-pipeline

# Run complete pipeline
sbatch scripts/slurm/00_run_full_pipeline.slurm
```

**That's it!** Data is included, no additional downloads required.

---

## Pipeline Overview

The pipeline processes three influenza lineages:
- **H1N1** (clustered at 99% identity)
- **H3N2** (clustered at 99% identity)  
- **VicB** (unclustered - low diversity)

### Processing Steps

1. **Filter sequences** - Quality control, length filtering
2. **Split by year** - Organize sequences 2009-2025
3. **Cluster sequences** - CD-HIT at 99% (H1N1, H3N2 only)
4. **Per-year analysis** - For each year:
   - Align sequences (MAFFT)
   - Build phylogenetic tree (IQ-TREE)
   - Generate 4 vaccine designs:
     - Consensus sequence
     - Medoid (most representative)
     - Ancestral (phylogenetic reconstruction)
     - COBRA (2-round clustering)
5. **Calculate distances** - ML distances (LG+G4 model)
6. **Visualize results** - Generate publication figures

---

## Directory Structure
```
flu-vaccine-design-distance-analysis/
├── data/
│   ├── raw/                    # Input sequences (included)
│   │   ├── H1N1/
│   │   ├── H3N2/
│   │   └── VicB/
│   └── processed/              # Filtered/clustered sequences
├── scripts/
│   ├── 01_filter_sequences.sh
│   ├── 02_split_by_year.sh
│   ├── 03_cluster_sequences.sh
│   ├── 04_per_year_analysis.sh
│   ├── 05_calculate_distances.sh
│   ├── 06_visualize_distances.sh
│   └── slurm/
│       ├── 00_run_full_pipeline.slurm  # Master script
│       ├── 01_filter_sequences.slurm
│       ├── 02_split_by_year.slurm
│       ├── 03_cluster_sequences.slurm
│       ├── 04_per_year_analysis.slurm
│       ├── 05_calculate_distances.slurm
│       └── 06_visualize_distances.slurm
└── results/
    ├── per_year_clustered/
    │   ├── H1N1/
    │   │   ├── alignments/
    │   │   ├── designs/
    │   │   ├── distances/
    │   │   ├── trees/
    │   │   ├── figures/        # Publication-ready outputs
    │   │   └── logs/
    │   └── H3N2/
    └── per_year_unclustered/
        └── VicB/
```

---

## Running the Pipeline

### Method 1: Complete Pipeline (Recommended)

Submit the master SLURM script that runs all steps with proper dependencies:
```bash
sbatch scripts/slurm/00_run_full_pipeline.slurm
```

**What happens:**
1. Master script submits 18 jobs
2. Jobs run in dependency order:
   - Filter → Split → Cluster → Per-year → Distances → Visualize
3. Jobs run in parallel where possible (3 lineages simultaneously)
4. Total runtime: 2-4 hours on HPC

**Monitor progress:**
```bash
# Check job status
squeue -u $USER

# Watch real-time (updates every 5 seconds)
watch -n 5 squeue -u $USER

# Check logs
tail -f logs/master_pipeline_*.out
```

---

### Method 2: Individual Steps

Run each step manually for a specific lineage:

#### Step 1: Filter Sequences
```bash
# Option A: Run with bash
bash scripts/01_filter_sequences.sh H1N1

# Option B: Submit as SLURM job
sbatch --export=LINEAGE=H1N1 scripts/slurm/01_filter_sequences.slurm
```

**Output:** `data/processed/H1N1/H1N1_filtered.fasta`

---

#### Step 2: Split by Year
```bash
# Run with bash
bash scripts/02_split_by_year.sh H1N1 2009 2025
```

**Output:** `data/processed/H1N1/per_year/H1N1_2009.fasta`, etc.

---

#### Step 3: Cluster Sequences (H1N1, H3N2 only)
```bash
# Option A: Run with bash
bash scripts/03_cluster_sequences.sh H1N1

# Option B: Submit as SLURM job
sbatch --export=LINEAGE=H1N1 scripts/slurm/03_cluster_sequences.slurm
```

**Clusters at 99% identity using CD-HIT**

**Output:** `data/processed/H1N1/clustered/H1N1_2009_clustered.fasta`

---

#### Step 4: Per-Year Analysis
```bash
# Option A: Run with bash
bash scripts/04_per_year_analysis.sh H1N1 2009 2025 clustered

# Option B: Submit as SLURM job
sbatch --export=LINEAGE=H1N1,START_YEAR=2009,END_YEAR=2025,USE_UNCLUSTERED=false \
       scripts/slurm/04_per_year_analysis.slurm
```

**For each year, creates:**
- Aligned sequences
- Phylogenetic tree
- 4 vaccine designs (consensus, medoid, ancestral, COBRA)
- Tree with designs included

**Output:** `results/per_year_clustered/H1N1/`

---

#### Step 5: Calculate Distances
```bash
# Option A: Run with bash
bash scripts/05_calculate_distances.sh H1N1 2009 2025 clustered

# Option B: Submit as SLURM job
sbatch --export=LINEAGE=H1N1,SOURCE_TYPE=clustered \
       scripts/slurm/05_calculate_distances.slurm
```

**For each design and year:**
- Computes ML distance (LG+G4 model) between design and each circulating strain
- Generates summary statistics

**Output:** `results/per_year_clustered/H1N1/distances/`

---

#### Step 6: Visualize Results
```bash
# Run with bash
bash scripts/06_visualize_distances.sh H1N1 clustered
```

**Creates:**
- `H1N1_yearly_comparison.png` - Distance trends over time
- `H1N1_overall_comparison.png` - Mean distances by design type
- `H1N1_summary_statistics.csv` - Numerical summary

**Output:** `results/per_year_clustered/H1N1/figures/`

---

## Expected Runtime

On Aoraki HPC (8 CPUs, 16GB RAM):

| Step | H1N1 | H3N2 | VicB | Notes |
|------|------|------|------|-------|
| Filter | 2 min | 2 min | 2 min | Sequential |
| Split | 1 min | 1 min | 1 min | Sequential |
| Cluster | 5 min | 5 min | N/A | Parallel |
| Per-year | 50 min | 67 min | 4-6 hrs | **Longest step** |
| Distances | 10 min | 12 min | 15 min | Sequential |
| Visualize | 2 min | 2 min | 2 min | Sequential |
| **Total** | **~1.2 hrs** | **~1.5 hrs** | **~5-7 hrs** | Parallel |

**Overall pipeline:** 2-4 hours (depends on VicB large datasets)

---

## Pipeline Details

### Vaccine Design Strategies

1. **Consensus Sequence**
   - Most common amino acid at each position
   - Simple, interpretable
   - May not represent any real strain

2. **Medoid Sequence**
   - Real strain closest to all others
   - Guaranteed viable sequence
   - Represents "average" strain

3. **Ancestral Sequence**
   - Reconstructed using phylogenetic tree (IQ-TREE ASR)
   - Theoretically optimal for past evolution
   - May include extinct variants

4. **COBRA (Computationally Optimized Broadly Reactive Antigen)**
   - 2-round clustering approach
   - Combines diversity from multiple clusters
   - Designed for broad coverage

### Distance Metrics

**ML Distance (LG+G4):**
- Maximum likelihood distance
- Le-Gascuel amino acid substitution model
- Gamma distribution for rate variation
- Accounts for multiple substitutions at same site
- **More accurate than simple p-distance**

Calculated between each vaccine design and all circulating strains for that year.

### Quality Control

- Minimum sequence length: 550 amino acids
- Remove sequences with ambiguous characters
- Remove duplicate sequences
- Cluster at 99% identity (H1N1, H3N2)
- Skip years with <3 sequences

---

## Output Files

### Per Lineage (9 files total):

**Publication Figures (3):**
- `*_yearly_comparison.png` - Time series of distances
- `*_overall_comparison.png` - Design comparison
- `*_summary_statistics.csv` - Numerical data

**Intermediate Files (tracked in results/, excluded from Git):**
- Alignments (MAFFT output)
- Trees (IQ-TREE output)
- Designs (FASTA sequences)
- Distance matrices (CSV files)

---

## Troubleshooting

### Job Fails with Memory Error

**Solution:** Increase memory in SLURM script
```bash
#SBATCH --mem=32G  # Default is 16G
```

### Per-Year Analysis Stuck on Large Year

**Symptom:** Job runs for >2 hours on single year

**Explanation:** Years with 500+ sequences take much longer
- 2023-2024 VicB can have 600+ sequences
- IQ-TREE with 1000 bootstraps is slow
- Expected behavior, not a bug

**Solution:** Be patient or reduce bootstrap replicates in script

### Visualization Creates Old Heatmap

**Solution:** You have old scripts - pull latest:
```bash
git pull origin main
```

Current version creates 3 plots (no heatmap).

### Distance Files Empty

**Symptom:** CSV files exist but have no data

**Check:** ML distance files from IQ-TREE
```bash
ls results/per_year_clustered/H1N1/trees/*_with_designs.mldist
```

If missing, per-year analysis didn't complete properly.

---

## Advanced Usage

### Run Single Year
```bash
# Process only 2023 for H1N1
bash scripts/04_per_year_analysis.sh H1N1 2023 2023 clustered
```

### Run Specific Design Strategy
Edit `scripts/04_per_year_analysis.sh` and comment out unwanted designs.

### Change Clustering Threshold
Edit `scripts/03_cluster_sequences.sh`:
```bash
cd-hit -c 0.99  # Change to 0.95 for 95% identity
```

### Use Different Substitution Model
Edit `scripts/04_per_year_analysis.sh`:
```bash
iqtree -m LG+G4  # Change to JTT+I+G4 or other
```

---

## Citation

If you use this pipeline, please cite:

**Software:**
```
Rana, Z. (2025). Influenza Vaccine Design and Distance Analysis Pipeline.
GitHub: https://github.com/zz153/flu-vaccine-design-distance-analysis
```

**Tools:**
- MAFFT: Katoh & Standley (2013) *Mol Biol Evol* 30(4):772-780
- IQ-TREE: Nguyen et al. (2015) *Mol Biol Evol* 32(1):268-274
- CD-HIT: Fu et al. (2012) *Bioinformatics* 28(23):3150-3152

---

## Support

For questions or issues:
1. Check this WORKFLOW.md
2. Review [REPRODUCIBILITY.md](REPRODUCIBILITY.md)
3. Open an issue on GitHub
