# ⚙️ Environment Setup Guide

This document describes how to install all dependencies required to run the  
**FluHub: Vaccine Design and Evolutionary Analysis Pipeline**.  
It ensures a reproducible setup on both local systems and HPC clusters, without requiring Conda environments.

---

## 🧩 1. System Requirements

The pipeline has been tested on:

- **Linux (Ubuntu 22.04, NeSI Aoraki HPC)**
- **macOS (Apple Silicon M2, Ventura and Sonoma)**

Minimum recommended specifications:
- Python ≥ 3.11  
- 8 GB RAM  
- 4 CPU cores  

---

## 🧰 2. Required Software and Versions

| Tool | Version Used | Purpose |
|------|---------------|----------|
| **MAFFT** | v7.505 (2022/Apr/10) | Multiple sequence alignment per lineage and year |
| **IQ-TREE 2** | v2.2.6 (COVID-edition, Dec 2023) | Maximum-likelihood phylogenetic tree inference |
| **CD-HIT** | v4.7 (built Jul 2018) | Clustering of sequences at 99 % identity |
| **SeqKit** | v2.10.1 (2025) | Sequence filtering, trimming, and ambiguity removal |
| **Python** | v3.11.13 | Core scripting environment |

### 🧱 Installation

#### 🐧 On Linux (Aoraki / Ubuntu)
Most tools are available through **modules** or **Bioconda**:
bash
module load mafft iqtree cd-hit
conda install -c bioconda seqkit
🍎 On macOS (Apple Silicon)
Use Homebrew for native ARM builds:

bash
Copy code
brew install mafft iqtree cd-hit seqkit
Verify:

bash
Copy code
mafft --version
iqtree2 --version
cd-hit -h | head -n 1
seqkit version
## 🐍 3. Python Packages and Versions
The following Python libraries were used in analyses:

| **Package** |	**Version** |	**Purpose** |
|-------------|-------------|-------------|
| biopython |	1.83	| FASTA parsing, sequence manipulation, alignment handling |
| pandas	| 2.3.3	| Data wrangling and summary tables |
| numpy	| 2.0.1	| Matrix operations and entropy calculations |
| matplotlib	| 3.10.7	| Plot generation for distances and entropy |
| seaborn	| 0.13.2	| Statistical visualization |
| scipy	| 1.16.2	| Shannon entropy and numerical utilities |
| tqdm	| 4.67.1	| Progress bars for computations |
| logomaker	| 0.8.7 | Sequence logo generation from alignments |

Installation
Install all packages at once:

bash
Copy code
pip install biopython==1.83 pandas==2.3.3 numpy==2.0.1 \
matplotlib==3.10.7 seaborn==0.13.2 scipy==1.16.2 tqdm==4.67.1 logomaker==0.8.7
Verify installation:

bash
Copy code
python -c "import Bio, pandas, numpy, matplotlib, seaborn, scipy, tqdm, logomaker; print('✅ Python packages loaded successfully!')"


