# Data Acquisition Guide

## Overview

This pipeline uses **hemagglutinin (HA) protein sequences** for three influenza lineages:

- A/H1N1pdm09
- A/H3N2
- B/Victoria (VicB)

**Time period:** 2009–2025  
**Geographic regions:** Europe, USA, Oceania, Asia (downloaded separately)

The data are downloaded **directly as protein FASTA files** from the GISAID EpiFlu database.

---

## Step 1: Access GISAID EpiFlu

1. Create a free account at **GISAID** (https://gisaid.org)  
2. Log in and navigate to the **EpiFlu™** database  
3. Read and accept the **Data Access Agreement (DAA)**

You must comply with GISAID’s terms of use.

---

## Step 2: Set Search Criteria (per lineage & region)

For **each lineage** (H1N1, H3N2, VicB) and **each region** (Asia, Europe, USA, Oceania), run a separate query.

### Example: A/H1N1pdm09 (HA protein)

In EpiFlu search:

- **Type:** A  
- **Subtype:** H1N1pdm09  
- **Segment:** HA (Hemagglutinin)  
- **Host:** Human  
- **Collection date:** `2009-01-01` to `2025-12-31`  
- **Location:** One of:
  - Asia  
  - Europe  
  - North America → USA  
  - Oceania  
- **Sequence completeness / quality:**
  - Complete sequences only: **Yes**  
  - High coverage only: **Yes** (or equivalent quality filter in the interface)
- **Sequence type / output:**
  - Protein / AA sequence view (HA protein)
  - Export as **FASTA**

Repeat analogous filters for:

- **A/H3N2** (Type A, Subtype H3N2, Segment HA)  
- **B/Victoria** (Type B, Lineage Victoria, Segment HA)

---

## Step 3: Download as FASTA

For each (lineage, region) combination:

1. Apply the filters above.
2. Use the **Download** function to export **protein FASTA**.
3. Save with a clear name, e.g.:

For H1N1:

- `Asia_H1N1_gisaid_epiflu_sequence.fasta`
- `Europe_H1N1_gisaid_epiflu_sequence.fasta`
- `USA_H1N1_HA_gisaid_epiflu_sequence.fasta`
- `Oceania_H1N1_HA_gisaid_epiflu_sequence.fasta`

For H3N2:

- `Asia_H3N2_HA_gisaid_epiflu_sequence.fasta`
- `Europe_H3N2_HA_gisaid_epiflu_sequence.fasta`
- `USA_H3N2_HA_gisaid_epiflu_sequence.fasta`
- `Oceania_H3N2_HA_gisaid_epiflu_sequence.fasta`

For VicB:

- `Asia_VicB_gisaid_epiflu_sequence.fasta`
- `Europe_VicB_gisaid_epiflu_sequence.fasta`
- `USA_VicB_gisaid_epiflu_sequence.fasta`
- `Oceania_VicB_gisaid_epiflu_sequence.fasta`

> The exact filenames do not matter as long as they match what your pipeline expects.
> In this repository, we use the naming pattern above.

---

## Step 4: Organize Files in `data/raw/`

After downloading, place files into lineage-specific subdirectories:

```text
data/raw/H1N1/
  ├── Asia_H1N1_gisaid_epiflu_sequence.fasta
  ├── Europe_H1N1_gisaid_epiflu_sequence.fasta
  ├── Oceania_H1N1_HA_gisaid_epiflu_sequence.fasta
  └── USA_H1N1_HA_gisaid_epiflu_sequence.fasta

data/raw/H3N2/
  ├── Asia_H3N2_HA_gisaid_epiflu_sequence.fasta
  ├── Europe_H3N2_HA_gisaid_epiflu_sequence.fasta
  ├── Oceania_H3N2_HA_gisaid_epiflu_sequence.fasta
  └── USA_H3N2_HA_gisaid_epiflu_sequence.fasta

data/raw/VicB/
  ├── Asia_VicB_gisaid_epiflu_sequence.fasta
  ├── Europe_VicB_gisaid_epiflu_sequence.fasta
  ├── Oceania_VicB_gisaid_epiflu_sequence.fasta
  └── USA_VicB_gisaid_epiflu_sequence.fasta

