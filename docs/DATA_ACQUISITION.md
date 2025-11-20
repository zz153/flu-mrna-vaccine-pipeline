# Data Acquisition Guide

## Overview

This pipeline requires hemagglutinin (HA) protein sequences for three influenza lineages:
- **A/H1N1pdm09**
- **A/H3N2**
- **B/Victoria**

**Time period:** 2009-2025  
**Geographic regions:** Europe, USA, Oceania, Asia (downloaded separately)

---

## Step 1: Access GISAID Database

1. Create a free account at [GISAID](https://www.gisaid.org/)
2. Navigate to **EpiFlu** database
3. Agree to the Data Access Agreement

---

## Step 2: Download Sequences by Region

### Download Strategy

For EACH lineage (H1N1, H3N2, VicB), download sequences SEPARATELY for each region.

**Why region-specific downloads?**
- Allows regional diversity analysis
- Better quality control per region
- Easier to track geographic patterns

---

### Search Criteria (Repeat for Each Region)

**For H1N1pdm09:**
- **Type:** A
- **Subtype:** H1N1pdm09
- **Segment:** HA (Hemagglutinin)
- **Host:** Human
- **Collection Date:** 2009-01-01 to 2025-12-31
- **Location:** [Select one: Europe / North America:USA / Oceania / Asia]
- **Complete sequences only:** Yes
- **High coverage only:** Yes
- **Protein:** Check (to get amino acid sequences)

**Download Format:**
- Format: FASTA or XLS (Excel file with sequences and metadata)
- File naming: `[Region]_[Lineage]_gisaid_epiflu_isolates.xls`

**Repeat for H3N2 and VicB**

---

## Step 3: Expected Downloaded Files

After downloading, you should have files like:
```
Asia_H1N1_gisaid_epiflu_isolates.xls
Europe_H1N1_gisaid_epiflu_isolates.xls
USA_H1N1_gisaid_epiflu_isolates.xls
Oceania_H1N1_gisaid_epiflu_isolates.xls

Asia_H3N2_gisaid_epiflu_isolates.xls
Europe_H3N2_gisaid_epiflu_isolates.xls
USA_H3N2_gisaid_epiflu_isolates.xls
Oceania_H3N2_gisaid_epiflu_isolates.xls

Asia_VicB_gisaid_epiflu_isolates.xls
Europe_VicB_gisaid_epiflu_isolates.xls
USA_VicB_gisaid_epiflu_isolates.xls
Oceania_VicB_gisaid_epiflu_isolates.xls
```

---

## Step 4: Extract FASTA Sequences

If you downloaded XLS files with metadata, extract the sequences:

**Option 1: Manual (in Excel)**
1. Open the XLS file
2. Find the column with HA sequences
3. Export as FASTA format

**Option 2: Using Python** (we'll provide a script)
```bash
python scripts/utils/extract_sequences_from_xls.py --input Asia_H1N1_gisaid_epiflu_isolates.xls --output data/raw/H1N1/Asia_H1N1.fasta
```

---

## Step 5: Organize Data

Place extracted FASTA files in subdirectories by lineage:
```
data/raw/H1N1/
  ├── Asia_H1N1.fasta
  ├── Europe_H1N1.fasta
  ├── USA_H1N1.fasta
  └── Oceania_H1N1.fasta

data/raw/H3N2/
  ├── Asia_H3N2.fasta
  ├── Europe_H3N2.fasta
  ├── USA_H3N2.fasta
  └── Oceania_H3N2.fasta

data/raw/VicB/
  ├── Asia_VicB.fasta
  ├── Europe_VicB.fasta
  ├── USA_VicB.fasta
  └── Oceania_VicB.fasta
```

---

## Step 6: Combine Regional Files (Optional)

The pipeline can work with either:
- **Individual regional files** (for region-specific analysis)
- **Combined files per lineage** (for pan-temporal analysis)

To combine all regions:
```bash
cat data/raw/H1N1/*.fasta > data/raw/H1N1/H1N1_all_regions.fasta
```

---

## Data Privacy

⚠️ **IMPORTANT:** GISAID sequences are subject to terms of use:
- Do not redistribute raw sequence files publicly
- Acknowledge GISAID and submitting laboratories in publications
- Raw data is `.gitignore`d and not included in this repository

---

## Expected Data Volume

Approximate sequences per lineage per region (as of 2025):
- **H1N1:** 10,000-25,000 per region
- **H3N2:** 20,000-40,000 per region  
- **VicB:** 7,000-15,000 per region

After quality filtering (≥550 aa) and CD-HIT clustering (99% identity), expect ~1,000-5,000 sequences per lineage per year.
