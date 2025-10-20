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

- Remove sequences **shorter than 550 amino acids**
- Remove **duplicate sequences**
- **Cluster sequences at 99% identity** using `cd-hit`
- **Remove ambiguous amino acids** (non-ACDEFGHIKLMNPQRSTVWY characters)
- **Retain only years 2009–2025**

