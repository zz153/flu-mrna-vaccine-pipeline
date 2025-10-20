# FluHub: Vaccine Design and Evolutionary Analysis Pipeline

## 🧫 Sequence Acquisition and Pre-Processing

### 1️⃣ Obtain Sequences
Download HA (hemagglutinin) protein sequences for the **main influenza lineages**:
- **A/H1N1**
- **A/H3N2**
- **B/Victoria (VicB)**  

from public databases (e.g., GISAID, NCBI, or Influenza Research Database), covering the regions:
**Europe, USA, Oceania, and Asia**,  
for the years **2009–2025**.

Organize them by lineage and year (e.g. `H1N1_2025.fasta`, `H3N2_2018.fasta`, etc.).

---

### 2️⃣ Clean Up Sequences

Perform stringent filtering to ensure data quality and consistency across years and lineages:

- Remove sequences shorter than 550 amino acids
- Remove duplicate sequences
- Cluster sequences at 99% identity using `cd-hit`
- Remove ambiguous amino acids (non-ACDEFGHIKLMNPQRSTVWY characters)
- Retain only years 2009–2025

### 3️⃣ Vaccine Design & Evolutionary Analyses

For each lineage (H1N1, H3N2, VicB) and year (2009–2025, plus combined dataset):

**Design sequences**

Consensus

Medoid

Ancestral (ASR)

COBRA design

**Compute p-distances**
Calculate pairwise genetic distances between:

Circulating strains

Designed vaccine candidates

Generate **per-year** and **combined heatmaps** and summary tables of distances.

### 4️⃣ Sequence Diversity Analyses

To visualize evolutionary variability and antigenic drift:

**Shannon Entropy Plots:**
Compute per-position Shannon entropy to quantify amino acid variability across alignments.

**Sequence Logos:**
Generate sequence logo plots showing amino acid frequency and conservation patterns for each year or region.

These analyses help identify conserved vs variable regions within HA and evaluate design robustness against circulating diversity.

### 5️⃣ Output Summary

Each lineage produces:

**Output File	Description**
*_aligned.fasta	Cleaned and aligned sequence data
*_treefile	IQ-TREE phylogenetic tree
*_asr.fasta	Ancestral sequence reconstruction
*_designs.fasta	Consensus, medoid, COBRA, and ASR design sequences
*_distance_analysis_*.png	Heatmaps of mean/median evolutionary and p-distance
*_entropy_plot.png	Shannon entropy profile across positions
*_sequence_logo.png	Sequence logo showing conservation and diversity
*_summary.csv	Yearly and combined distance/entropy summaries
### 🧩 Summary

In short, the workflow:

- Downloads sequences (H1N1, H3N2, VicB; 2009–2025; global regions)

- Cleans and clusters them

- Designs vaccine candidates (consensus, medoid, ASR, COBRA)

- Computes p-distances between design and circulating strains

- Generates entropy and sequence logo plots to visualize sequence diversity

This provides a **comprehensive evolutionary** and **antigenic landscape** for influenza vaccine design.

