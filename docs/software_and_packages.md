# 🧰 Software and Packages Used

This document lists all software, command-line tools, and Python dependencies used in the **FluHub: Vaccine Design and Evolutionary Analysis Pipeline**.  
The environment ensures reproducibility across systems (local or HPC).

---

## ⚙️ System Environment

- **Platform:** Linux (tested on Ubuntu 22.04 and NeSI Aoraki HPC)
- **Workflow Manager:** SLURM 23+
- **Shell:** GNU Bash 5+
- **Package Management:** Conda 24+ (via `environment.yml`)

---

## 🔹 Core Bioinformatics Software

| Tool | Version | Description |
|------|----------|-------------|
| **MAFFT** | ≥7.520 | Multiple sequence alignment per lineage and year |
| **IQ-TREE 2** | ≥2.3 | Maximum-likelihood phylogenetic tree inference and evolutionary distance estimation |
| **CD-HIT** | ≥4.8 | Clustering of sequences at 99% identity to remove redundancy |
| **SeqKit** | ≥2.8 | Sequence manipulation, length filtering, and removal of ambiguous residues |
| **Python** | ≥3.11 | Core scripting and analysis environment |


---

## 🔹 Python Packages

| Package | Purpose |
|----------|----------|
| **biopython** | FASTA parsing, sequence manipulation, and alignment handling |
| **pandas** | Data wrangling and summary tables |
| **numpy** | Matrix operations, distance, and entropy calculations |
| **matplotlib** | Plot generation for distances and entropy |
| **seaborn** | Aesthetic statistical visualizations |
| **logomaker** | Sequence logo generation from alignments |
| **scipy** | Shannon entropy and mathematical utilities |
| **tqdm** | Progress bars for long computations |
| **argparse** | CLI argument parsing in Python scripts |

---

## 🔹 HPC and Automation Tools

| Component | Description |
|------------|-------------|
| **SLURM** | Job submission and resource allocation across multiple nodes |
| **Bash scripts** (`*.sh`) | Automate per-year and per-lineage analyses |
| **SLURM scripts** (`*.slurm`) | Schedule and parallelize MAFFT, IQ-TREE, and COBRA jobs |

---

## 🔹 Environment Setup

To reproduce the full environment:

```bash
conda env create -f environment.yml
conda activate flu_design
