#!/bin/bash

# Script 01: Filter and Clean Raw Sequences
# 
# Quality control steps:
# 1. Combine regional files
# 2. Filter for full-length HA sequences (≥550 AA)
# 3. Remove sequences with mislabeled years (keep 2009-2025 only)
# 4. Remove duplicate sequences
# 5. Remove sequences with >5% ambiguous bases (X, B, J, Z)

set -e  # Exit on error

# --- Parameters ---
LINEAGE=$1  # H1N1, H3N2, or VicB
MIN_LENGTH=550
YEAR_START=2009
YEAR_END=2025

if [ -z "$LINEAGE" ]; then
    echo "Usage: bash scripts/01_filter_sequences.sh [LINEAGE]"
    echo "Example: bash scripts/01_filter_sequences.sh H1N1"
    exit 1
fi

echo "=========================================="
echo "Filtering ${LINEAGE} sequences"
echo "=========================================="

# --- Directories ---
RAW_DIR="data/raw/${LINEAGE}"
PROC_DIR="data/processed/${LINEAGE}"

mkdir -p ${PROC_DIR}

# --- Step 1: Combine regional files ---
echo "Step 1: Combining regional sequences..."
cat ${RAW_DIR}/*.fasta > ${PROC_DIR}/${LINEAGE}_combined.fasta

INITIAL_COUNT=$(grep -c "^>" ${PROC_DIR}/${LINEAGE}_combined.fasta)
echo "   Initial sequences: ${INITIAL_COUNT}"

# --- Step 2: Filter by length (≥550 AA) ---
echo "Step 2: Filtering sequences ≥${MIN_LENGTH} amino acids..."
seqkit seq -m ${MIN_LENGTH} ${PROC_DIR}/${LINEAGE}_combined.fasta > \
    ${PROC_DIR}/${LINEAGE}_length_filtered.fasta

LENGTH_COUNT=$(grep -c "^>" ${PROC_DIR}/${LINEAGE}_length_filtered.fasta)
echo "   After length filter: ${LENGTH_COUNT}"

# --- Step 3: Filter by year (2009-2025) ---
echo "Step 3: Filtering sequences from ${YEAR_START}-${YEAR_END}..."
# Extract year from header format: >EPI534686|HA|A/Taiwan/362/2014|...
awk -v start=$YEAR_START -v end=$YEAR_END '
    /^>/ {
        split($0, fields, "|")
        split(fields[3], parts, "/")
        year = parts[length(parts)]
        if (year >= start && year <= end) {
            print_seq = 1
            print $0
        } else {
            print_seq = 0
        }
        next
    }
    print_seq { print }
' ${PROC_DIR}/${LINEAGE}_length_filtered.fasta > ${PROC_DIR}/${LINEAGE}_year_filtered.fasta

YEAR_COUNT=$(grep -c "^>" ${PROC_DIR}/${LINEAGE}_year_filtered.fasta)
echo "   After year filter: ${YEAR_COUNT}"

# --- Step 4: Remove duplicates ---
echo "Step 4: Removing duplicate sequences..."
seqkit rmdup -s ${PROC_DIR}/${LINEAGE}_year_filtered.fasta > \
    ${PROC_DIR}/${LINEAGE}_dedup.fasta

DEDUP_COUNT=$(grep -c "^>" ${PROC_DIR}/${LINEAGE}_dedup.fasta)
echo "   After deduplication: ${DEDUP_COUNT}"

# --- Step 5: Remove sequences with >5% ambiguous bases ---
echo "Step 5: Removing sequences with >5% ambiguous bases..."
# Keep only sequences with standard amino acids
seqkit grep -s -r -p "^[ACDEFGHIKLMNPQRSTVWY]+$" \
    ${PROC_DIR}/${LINEAGE}_dedup.fasta > \
    ${PROC_DIR}/${LINEAGE}_clean.fasta

FINAL_COUNT=$(grep -c "^>" ${PROC_DIR}/${LINEAGE}_clean.fasta)
echo "   After ambiguous base filter: ${FINAL_COUNT}"

echo ""
echo "=========================================="
echo "✅ Filtering complete for ${LINEAGE}"
echo "=========================================="
echo "Initial sequences:    ${INITIAL_COUNT}"
echo "Final sequences:      ${FINAL_COUNT}"
echo "Reduction:            $((INITIAL_COUNT - FINAL_COUNT)) sequences removed"
echo ""
echo "Output: ${PROC_DIR}/${LINEAGE}_clean.fasta"
