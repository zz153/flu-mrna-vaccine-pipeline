# 🧰 Software and Packages Used

This repository lists all software, command-line tools, and Python dependencies used in the **FluHub: Vaccine Design and Evolutionary Analysis Pipeline**.  
The environment ensures full reproducibility across systems (local or HPC) via the [`environment.yml`](./environment.yml) file.

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
| **MAFFT** | v7.505 | Multiple sequence alignment per lineage and year |
| **IQ-TREE 2** | v2.2.6 | Maximum-likelihood phylogenetic tree inference and evolutionary distance estimation |
| **CD-HIT** | v4.7 | Clustering of sequences at 99% identity to remove redundancy |
| **SeqKit** | v2.10.1 | Sequence manipulation, length filtering, and removal of ambiguous residues |
| **Python** | v3.11.13 | Core scripting and analysis environment |

---

## 🔹 Python Packages

| Package | Version | Purpose |
|----------|----------|----------|
| **biopython** | 1.83 | FASTA parsing, sequence manipulation, and alignment handling |
| **pandas** | 2.3.3 | Data wrangling and summary tables |
| **numpy** | 2.0.1 | Matrix operations, distance, and entropy calculations |
| **matplotlib** | 3.10.7 | Plot generation for distances and entropy |
| **seaborn** | 0.13.2 | Aesthetic statistical visualizations |
| **logomaker** | 0.8.7 | Sequence logo generation from alignments |
| **scipy** | 1.16.2 | Shannon entropy and mathematical utilities |
| **tqdm** | 4.67.1 | Progress bars for long computations |
| **argparse** | Built-in | CLI argument parsing in Python scripts |

---

## 🔹 HPC and Automation Tools

| Component | Description |
|------------|-------------|
| **SLURM** | Job submission and resource allocation across multiple nodes |
| **Bash scripts (`*.sh`)** | Automate per-year and per-lineage analyses |
| **SLURM scripts (`*.slurm`)** | Schedule and parallelize MAFFT, IQ-TREE, and COBRA jobs |

---

## 🔹 Environment Setup

Reproduce the full Conda environment:

```bash
conda env create -f environment.yml
conda activate flu_design

