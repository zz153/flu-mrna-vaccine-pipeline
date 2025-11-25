# Reproducibility Guide

Complete instructions for reproducing the influenza vaccine design analysis.

---

## 🎯 System Requirements

- **OS:** Linux or macOS (tested on Ubuntu 24.04)
- **RAM:** 32GB minimum (64GB recommended for VicB)
- **Storage:** 100GB free space
- **CPU:** 8+ cores recommended
- **Time:** 8-12 hours for complete analysis

---

## 📥 Step 1: Clone Repository
```bash
git clone https://github.com/zz153/flu-mrna-vaccine-pipeline.git
cd flu-mrna-vaccine-pipeline
```

---

## 🔧 Step 2: Set Up Environment
```bash
# Create conda environment
conda env create -f env/environment.yml

# Activate environment
conda activate flu-vaccine-pipeline

# Verify installation
bash scripts/run_full_pipeline.sh check
```

**Expected output:**
```
[✓] Conda environment 'flu-vaccine-pipeline' found
[✓] MAFFT installed
[✓] IQ-TREE installed
[✓] CD-HIT installed
[✓] All prerequisites satisfied!
```

---

## 📊 Step 3: Obtain Data

### Download from GISAID

1. **Register** at https://gisaid.org (free)
2. **Follow** instructions in `docs/DATA_ACQUISITION.md`
3. **Download** HA protein sequences for H1N1, H3N2, VicB (2009-2025)
4. **Place** in `data/raw/`:
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

## 🚀 Step 4: Run Pipeline

### **Recommended: Automatic Execution**
```bash
bash scripts/run_full_pipeline.sh local
```

Processes all 3 lineages sequentially (~8-12 hours)

### **Alternative: SLURM (HPC)**
```bash
sbatch scripts/slurm/00_run_full_pipeline.slurm
squeue -u $USER  # Monitor
```

Parallel execution (~2-4 hours)

### **Alternative: Single Lineage**
```bash
bash scripts/01_filter_sequences.sh H1N1
bash scripts/02_cluster_sequences.sh H1N1
bash scripts/02b_split_by_year.sh H1N1
bash scripts/03_align_sequences.sh H1N1
bash scripts/04_per_year_analysis.sh H1N1 2009 2025 4
bash scripts/05_calculate_distances.sh H1N1 2009 2025 clustered
bash scripts/06_visualize_distances.sh H1N1 clustered
```

---

## ✅ Step 5: Verify Results

### File Counts
```bash
# H1N1: 16 years
find results/per_year_clustered/H1N1/designs -name "*.fasta" | wc -l
# Expected: 128

# H3N2: 17 years
find results/per_year_clustered/H3N2/designs -name "*.fasta" | wc -l
# Expected: 136

# VicB: 17 years
find results/per_year_unclustered/VicB/designs -name "*.fasta" | wc -l
# Expected: 136
```

### Summary Statistics
```bash
cat results/per_year_clustered/H1N1/figures/H1N1_summary_statistics.csv
```

**Expected:**
```
Design,Mean_Distance,Median_Distance,Min_Distance,Max_Distance,Avg_SD,Years_Analyzed
Consensus,0.0142,0.0106,0.0102,0.0247,0.0162,16
Medoid,0.0179,0.0150,0.0106,0.0263,0.0172,16
Ancestral,0.0180,0.0150,0.0118,0.0275,0.0169,16
Cobra,0.0299,0.0221,0.0106,0.1092,0.0166,16
```
```bash
cat results/per_year_clustered/H3N2/figures/H3N2_summary_statistics.csv
```

**Expected:**
```
Design,Mean_Distance,Median_Distance,Min_Distance,Max_Distance,Avg_SD,Years_Analyzed
Consensus,0.1320,0.1295,0.0819,0.2039,0.0202,17
Cobra,0.1369,0.1331,0.0819,0.2229,0.0203,17
Medoid,0.1384,0.1357,0.0823,0.2207,0.0210,17
Ancestral,0.1831,0.1793,0.1151,0.2690,0.0207,17
```
```bash
cat results/per_year_unclustered/VicB/figures/VicB_summary_statistics.csv
```

**Expected:**
```
Design,Mean_Distance,Median_Distance,Min_Distance,Max_Distance,Avg_SD,Years_Analyzed
Consensus,0.0050,0.0034,0.0036,0.0067,0.0029,16
Medoid,0.0056,0.0051,0.0036,0.0076,0.0030,16
Cobra,0.0085,0.0077,0.0036,0.0191,0.0029,16
Ancestral,0.0086,0.0069,0.0053,0.0153,0.0033,16
```

---

## 📊 Reproducibility Levels

### Statistical Reproducibility (Same Data)
- Vaccine sequences: 100% identical
- P-distances: 100% identical
- ML distances: ≥99.9% identical
- Variation: ±0.01%

### Biological Reproducibility (New GISAID Data)
- Trends: Consistent
- Exact numbers: May vary ±0.5%
- Conclusions: Robust

---

## 🐛 Troubleshooting

**Command not found:**
```bash
conda activate flu-vaccine-pipeline
```

**Out of memory:**
```bash
bash scripts/04_per_year_analysis.sh H1N1 2009 2025 2  # Use 2 threads
```

**Different results:**
- New GISAID data (updated daily)
- Contact authors for exact processed data

---

## 📞 Support

- GitHub Issues: https://github.com/zz153/flu-mrna-vaccine-pipeline/issues
- Email: zohaib.rana@otago.ac.nz

---

## 📚 Citation

Software used:
- MAFFT: Katoh & Standley (2013)
- IQ-TREE: Nguyen et al. (2015)
- CD-HIT: Fu et al. (2012)

---

**Last Updated:** November 23, 2025  
**Pipeline Version:** v1.0

## Pre-Generated Figures

The repository includes pre-generated publication-ready figures for immediate viewing:
- `results/per_year_clustered/H1N1/figures/` (H1N1 lineage)
- `results/per_year_clustered/H3N2/figures/` (H3N2 lineage)  
- `results/per_year_unclustered/VicB/figures/` (VicB lineage)

**To reproduce:** Simply run the pipeline - all outputs are automatically overwritten with fresh results. No cleanup required.
```bash
sbatch scripts/slurm/00_run_full_pipeline.slurm
```
