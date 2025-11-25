# Influenza Vaccine Design and Distance Analysis Pipeline

A fully reproducible computational pipeline for designing and evaluating influenza vaccine candidates using evolutionary analysis, phylogenetic methods, and distance-based metrics.

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/)
[![Conda](https://img.shields.io/badge/conda-environment-green.svg)](env/environment.yml)

---

## 🎯 Overview

FluHub generates and evaluates four distinct vaccine design strategies for influenza hemagglutinin (HA) sequences:

1. **Consensus** - Computationally derived sequence with most common amino acids
2. **Medoid** - Actual circulating strain with minimum distance to all others
3. **Ancestral (ASR)** - Phylogenetically reconstructed ancestral sequence
4. **COBRA** - Computationally Optimized Broadly Reactive Antigen (2-round clustering)

Each design is evaluated against circulating strains using both simple p-distance and maximum likelihood evolutionary distances (LG+G4 model).

---

## 📊 Key Findings

Analysis of H1N1 and H3N2 (2009-2025) reveals:

**H1N1 (Lower evolutionary rate):**
- **Consensus performs best**: 1.42% mean distance to circulating strains
- **Medoid**: 1.79% mean distance
- **Ancestral**: 1.80% mean distance  
- **COBRA**: 2.99% mean distance (higher variability in recent years)

**H3N2 (Higher evolutionary rate - ~10× faster than H1N1):**
- **Consensus performs best**: 13.20% mean distance
- **COBRA**: 13.69% mean distance (better suited for high diversity)
- **Medoid**: 13.84% mean distance
- **Ancestral**: 18.31% mean distance (performs poorly in rapidly evolving lineages)

---

## 🧬 Sequence Data

### Lineages Analyzed
### Lineages Analyzed
- **A/H1N1** - 16 years (2009-2025, excluding 2015) - *Clustered at 99% identity*
- **A/H3N2** - 17 years (2009-2025) - *Clustered at 99% identity*
- **B/Victoria (VicB)** - 17 years (2009-2025) - *Unclustered (low diversity)*
### Geographic Coverage
- Europe
- USA
- Oceania
- Asia

### Data Source
Sequences obtained from GISAID (requires authenticated access).

---

## 🔬 Pipeline Workflow
```
┌─────────────────────────────────────────────────────────────┐
│  Script 01: Filter & Clean Sequences                        │
│  • Remove sequences < 550 aa                                │
│  • Remove ambiguous amino acids                              │
│  • Retain only 2009-2025                                    │
└─────────────────┬───────────────────────────────────────────┘
                  │
┌─────────────────▼───────────────────────────────────────────┐
│  Script 02: Cluster Sequences (CD-HIT 99%)                  │
│  • Reduce redundancy while preserving diversity             │
└─────────────────┬───────────────────────────────────────────┘
                  │
┌─────────────────▼───────────────────────────────────────────┐
│  Script 02b: Split by Year                                  │
│  • Create per-year datasets                                 │
└─────────────────┬───────────────────────────────────────────┘
                  │
┌─────────────────▼───────────────────────────────────────────┐
│  Script 03: Combined Alignment (MAFFT)                      │
│  • Generate lineage-wide reference alignment                │
└─────────────────┬───────────────────────────────────────────┘
                  │
┌─────────────────▼───────────────────────────────────────────┐
│  Script 04: Per-Year Analysis & Vaccine Design              │
│  • Align sequences (MAFFT)                                  │
│  • Build phylogenetic tree (IQ-TREE, LG+G4)                │
│  • Generate 4 vaccine designs:                              │
│    1. Consensus (most common amino acids)                   │
│    2. Medoid (minimum distance strain)                      │
│    3. Ancestral (ASR from tree root)                        │
│    4. COBRA (2-round CD-HIT clustering: 95% → 90%)          │
│  • Re-align COBRA to match gap structure                    │
│  • Rebuild tree with all designs for ML distance extraction │
└─────────────────┬───────────────────────────────────────────┘
                  │
┌─────────────────▼───────────────────────────────────────────┐
│  Script 05: Calculate Evolutionary Distances                │
│  • P-distance (observed differences) for all 4 designs      │
│  • ML distance (LG+G4 model) from trees with designs        │
│  • Generate per-year and summary statistics                 │
└─────────────────┬───────────────────────────────────────────┘
                  │
┌─────────────────▼───────────────────────────────────────────┐
│  Script 06: Visualize Performance                           │
│  • Year-by-year comparison plots                            │
│  • Clean heatmaps (no cluttered numbers)                    │
│  • Overall performance comparison                           │
│  • Summary statistics tables                                │
└─────────────────────────────────────────────────────────────┘
```

---

## 🚀 Quick Start

### 1. Clone the Repository
```bash
git clone https://github.com/zz153/flu-vaccine-design-distance-analysis.git
cd flu-vaccine-design-distance-analysis
```

### 2. Set Up Environment
```bash
# Create conda environment
conda env create -f env/environment.yml

# Activate environment
conda activate flu-vaccine-pipeline

# Verify installation
bash scripts/run_full_pipeline.sh check
```

### 3. Add Your Data

Place GISAID HA sequences in `data/raw/[LINEAGE]/`:
```
data/raw/
├── H1N1/
│   ├── Asia_H1N1_gisaid_epiflu_sequence.fasta
│   ├── Europe_H1N1_gisaid_epiflu_sequence.fasta
│   ├── USA_H1N1_HA_gisaid_epiflu_sequence.fasta
│   └── Oceania_H1N1_HA_gisaid_epiflu_sequence.fasta
├── H3N2/ (same structure)
└── VicB/ (same structure)
```
---

## 🎮 Execution Modes

This pipeline supports **two execution modes**:

| Mode | Best For | Time* | Requirements |
|------|----------|-------|--------------|
| **Bash (Local)** | Testing, debugging, no HPC | 6-12h | Linux/Mac, 8+ cores, 32GB RAM |
| **SLURM (HPC)** | Production, parallel processing | 2-4h | HPC cluster with SLURM |

*Complete 3-lineage analysis. Both modes produce **identical, reproducible results**.

---

### 🖥️ Option A: Local Execution (Bash)

**Best for:** Testing, development, single-lineage analysis, or systems without SLURM access

#### Complete Pipeline (All Lineages)
```bash
bash scripts/run_full_pipeline.sh local
```

#### Single Lineage (Example: H1N1)
```bash
bash scripts/01_filter_sequences.sh H1N1
bash scripts/02_cluster_sequences.sh H1N1
bash scripts/02b_split_by_year.sh H1N1
bash scripts/03_align_sequences.sh H1N1
bash scripts/04_per_year_analysis.sh H1N1 2009 2025 4
bash scripts/05_calculate_distances.sh H1N1 2009 2025 clustered
bash scripts/06_visualize_distances.sh H1N1 clustered
```

#### Quick Test (Single Year)
```bash
# Test pipeline on 2023 data only
bash scripts/04_per_year_analysis.sh H1N1 2023 2023 4
bash scripts/05_calculate_distances.sh H1N1 2023 2023 clustered
bash scripts/06_visualize_distances.sh H1N1 clustered
```

**Pros:** ✅ Immediate execution • ✅ Real-time output • ✅ Easy debugging  
**Cons:** ⏳ Sequential (slower) • 🔒 Single lineage at a time

---

### ⚡ Option B: HPC Execution (SLURM)

**Best for:** Production runs, complete analysis, parallel processing of all lineages

#### Complete Pipeline (Recommended)
```bash
# Submit master orchestrator - handles everything automatically
sbatch scripts/slurm/00_run_full_pipeline.slurm

# Monitor all jobs
squeue -u $USER

# Check master log
tail -f logs/pipeline_master_*.out
```

**What this does:**
- ✅ Processes H1N1, H3N2, and VicB in parallel
- ✅ Automatically sets up job dependencies
- ✅ Manages resource allocation
- ✅ Ensures correct execution order

#### Individual Lineage
```bash
# Example: H1N1 complete pipeline
sbatch --export=LINEAGE=H1N1 scripts/slurm/01_filter_sequences.slurm
sbatch --export=LINEAGE=H1N1 scripts/slurm/02_cluster_sequences.slurm
sbatch --export=LINEAGE=H1N1 scripts/slurm/02b_split_by_year.slurm
sbatch --export=LINEAGE=H1N1 scripts/slurm/03_align_sequences.slurm
sbatch --export=LINEAGE=H1N1 scripts/slurm/04_per_year_analysis.slurm
sbatch --export=LINEAGE=H1N1 scripts/slurm/05_calculate_distances.slurm
sbatch --export=LINEAGE=H1N1 scripts/slurm/06_visualize_distances.slurm
```

**Pros:** ⚡ Fast (parallel) • 🔄 Background processing • 📊 All lineages simultaneously  
**Cons:** ⏰ Queue wait time • 🖥️ Requires HPC access

---

### 🔀 Hybrid Approach (Recommended for Development)
```bash
# 1. Test locally first
bash scripts/04_per_year_analysis.sh H1N1 2023 2023 4
bash scripts/05_calculate_distances.sh H1N1 2023 2023 clustered

# 2. If successful, run full analysis on HPC
sbatch scripts/slurm/00_run_full_pipeline.slurm
```

---

## 📁 Repository Structure
```
flu-vaccine-design-distance-analysis/
├── data/
│   ├── raw/                    # Raw GISAID sequences (user-provided)
│   └── processed/              # Cleaned and clustered sequences
│       └── {LINEAGE}/
│           ├── {LINEAGE}_clean.fasta
│           ├── {LINEAGE}_clustered.fasta
│           └── per_year/       # Split by year
├── results/
│   ├── per_year_clustered/     # H1N1, H3N2 results
│   │   └── {LINEAGE}/
│   │       ├── alignments/     # MAFFT alignments
│   │       ├── trees/          # IQ-TREE phylogenies & ML distances
│   │       ├── designs/        # 4 vaccine designs (× 2 versions each)
│   │       ├── distances/      # P-distance & ML distance CSVs
│   │       └── figures/        # Publication-quality plots
│   └── per_year_unclustered/   # VicB results
├── scripts/
│   ├── 01_filter_sequences.sh
│   ├── 02_cluster_sequences.sh
│   ├── 02b_split_by_year.sh
│   ├── 03_align_sequences.sh
│   ├── 04_per_year_analysis.sh
│   ├── 05_calculate_distances.sh
│   ├── 06_visualize_distances.sh
│   ├── run_full_pipeline.sh        # Master bash script
│   └── slurm/                      # SLURM wrappers for HPC
│       ├── 00_run_full_pipeline.slurm  # Master orchestrator
│       ├── 01_filter_sequences.slurm
│       ├── 02_cluster_sequences.slurm
│       ├── 02b_split_by_year.slurm
│       ├── 03_align_sequences.slurm
│       ├── 04_per_year_analysis.slurm
│       ├── 05_calculate_distances.slurm
│       └── 06_visualize_distances.slurm
├── env/
│   └── environment.yml         # Conda environment specification
├── README.md                   # This file
├── REPRODUCIBILITY.md          # Detailed reproduction guide
└── .gitignore
```

---

## 🔧 Dependencies

All dependencies are specified in `env/environment.yml`:

- **Python 3.11** with Biopython, NumPy, Pandas, Matplotlib, Seaborn
- **MAFFT 7.526** - Multiple sequence alignment
- **IQ-TREE 3.0.1** - Maximum likelihood phylogenetic inference
- **CD-HIT 4.8.1** - Sequence clustering
- **SeqKit 2.10.1** - Sequence manipulation

---

## 📈 Output Files

### Per-Year Results

For each year and lineage:

| File | Description |
|------|-------------|
| `{LINEAGE}_{YEAR}_aligned.fasta` | MAFFT alignment of circulating strains |
| `{LINEAGE}_{YEAR}.treefile` | IQ-TREE phylogenetic tree |
| `{LINEAGE}_{YEAR}_with_designs.treefile` | Tree including all 4 designs |
| `{LINEAGE}_{YEAR}_consensus.fasta` | Consensus vaccine design (ungapped) |
| `{LINEAGE}_{YEAR}_medoid.fasta` | Medoid vaccine design |
| `{LINEAGE}_{YEAR}_ancestral.fasta` | Ancestral vaccine design (ASR) |
| `{LINEAGE}_{YEAR}_cobra.fasta` | COBRA vaccine design |
| `{LINEAGE}_{YEAR}_{design}_aligned.fasta` | Aligned version (for distances) |
| `{LINEAGE}_{YEAR}_{design}_distances.csv` | Distances to all strains |

### Summary Files

| File | Description |
|------|-------------|
| `distance_summary_{design}.csv` | Mean/median/SD distances per year |
| `{LINEAGE}_yearly_comparison.png` | Year-by-year distance plots |
| `{LINEAGE}_distance_heatmap.png` | Heatmap of all designs × years |
| `{LINEAGE}_overall_comparison.png` | 4-panel summary comparison |
| `{LINEAGE}_summary_statistics.csv` | Overall performance metrics |

---

## 🎨 Visualization Examples

### Year-by-Year Comparison
Shows mean p-distance and ML distance for all 4 designs across years.

### Distance Heatmap
Clean visualization of mean p-distances (no cluttered numbers).

### Overall Comparison
4-panel figure showing:
- Average performance across all years
- Median performance
- Strain counts per year
- Distance variability

---

## 🧪 Design Strategies Explained

### 1. Consensus Sequence
- **Method**: Most common amino acid at each position
- **Pros**: Computationally efficient, represents population average
- **Cons**: May not correspond to any real strain
- **Performance**: **Best** (H1N1: 1.42%, H3N2: 13.20%)

### 2. Medoid Sequence
- **Method**: Actual strain with minimum sum of distances to all others
- **Pros**: Real sequence, exists in nature
- **Cons**: May not capture all diversity
- **Performance**: Good (H1N1: 1.79%, H3N2: 13.84%)

### 3. Ancestral Sequence (ASR)
- **Method**: Maximum likelihood reconstruction of tree root sequence
- **Pros**: Phylogenetically informed, evolutionary perspective
- **Cons**: Hypothetical sequence from the past
- **Performance**: Variable (H1N1: 1.80%, H3N2: 18.31%)

### 4. COBRA (Computationally Optimized Broadly Reactive Antigen)
- **Method**: 2-round CD-HIT clustering (95% → 90% identity)
  1. Cluster all strains at 95% identity
  2. Align cluster representatives
  3. Create consensus
  4. Cluster again at 90% identity
  5. Final consensus of representatives
  6. Re-align to match original alignment gap structure
- **Pros**: Designed for broad coverage across diversity
- **Cons**: Complex, may over-optimize for past diversity
- **Performance**: Variable (H1N1: 2.99%, H3N2: 13.69%)

---

## 📊 Distance Metrics

### P-Distance (Simple)
```
p = (number of differences) / (number of aligned positions)
```
- Ignores gap-gap positions
- Ignores positions where either sequence has a gap
- Fast, interpretable

### ML Distance (Evolutionary)
```
Using LG+G4 model (IQ-TREE)
```
- Le-Gascuel amino acid substitution matrix
- Gamma rate heterogeneity (4 categories)
- Accounts for multiple substitutions at same site
- More accurate for evolutionary distances

---

## 🔬 Reproducibility

This pipeline is designed for complete reproducibility:

✅ **Version-controlled scripts** - All analysis code in Git  
✅ **Conda environment** - Exact software versions specified  
✅ **Parameterized workflows** - No hard-coded paths  
✅ **Dual execution modes** - Both bash and SLURM produce identical results  
✅ **Documented methods** - Clear README and comprehensive [REPRODUCIBILITY.md](REPRODUCIBILITY.md)

### To Reproduce Results:

1. Clone this repository
2. Create conda environment from `env/environment.yml`
3. Obtain HA sequences from GISAID (requires account)
4. Run scripts 01-06 in sequence (bash or SLURM)
5. Results should match published figures and statistics

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for detailed step-by-step instructions.

---

## 📝 Citation

If you use this pipeline in your research, please cite:
```
[Your citation information here]
```

---

## 👥 Authors

**Zohaib Rana**  
Postdoctoral Fellow  
Department of Biochemistry  
University of Otago  
RNA & Cancer Therapeutics

---

## 📄 License

[Specify your license here - e.g., MIT, GPL, etc.]

---

## 🤝 Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request with clear description of changes

---

## 📧 Contact

For questions or issues:
- Open an issue on GitHub
- Email: zohaib.rana@otago.ac.nz

---

## 🙏 Acknowledgments

- GISAID for sequence data access
- University of Otago HPC facility (Aoraki)
- Open-source bioinformatics community

---

**Last Updated**: November 2025  
**Pipeline Version**: v1.0  
**Execution Modes**: Bash (local) & SLURM (HPC) - both fully supported
