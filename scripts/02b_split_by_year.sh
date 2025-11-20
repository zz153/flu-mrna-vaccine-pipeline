#!/bin/bash

# Script 02b: Split Sequences by Year
#
# Extracts sequences for each year (2009-2025) separately
# Creates individual FASTA files for per-year analysis

set -e

# --- Parameters ---
LINEAGE=$1
YEAR_START=2009
YEAR_END=2025

if [ -z "$LINEAGE" ]; then
    echo "Usage: bash scripts/02b_split_by_year.sh [LINEAGE]"
    echo "Example: bash scripts/02b_split_by_year.sh H1N1"
    exit 1
fi

echo "=========================================="
echo "Splitting ${LINEAGE} sequences by year"
echo "=========================================="

# --- Directories ---
INPUT="data/processed/${LINEAGE}/${LINEAGE}_cdhit99.fasta"
OUTPUT_DIR="data/processed/${LINEAGE}/per_year"

mkdir -p ${OUTPUT_DIR}

# Check input
if [ ! -f "$INPUT" ]; then
    echo "Error: Input file not found: $INPUT"
    echo "Run clustering first: bash scripts/02_cluster_sequences.sh ${LINEAGE}"
    exit 1
fi

# --- Split by year ---
echo "Extracting sequences for each year..."

for YEAR in $(seq ${YEAR_START} ${YEAR_END}); do
    OUTPUT_FILE="${OUTPUT_DIR}/${LINEAGE}_${YEAR}.fasta"
    
    # Extract sequences from this year
    # Header format: >EPI534686|HA|A/Taiwan/362/2014|...
    awk -v year="$YEAR" '
        /^>/ {
            split($0, fields, "|")
            split(fields[3], parts, "/")
            seq_year = parts[length(parts)]
            if (seq_year == year) {
                print_seq = 1
                print $0
            } else {
                print_seq = 0
            }
            next
        }
        print_seq { print }
    ' ${INPUT} > ${OUTPUT_FILE}
    
    # Count sequences
    if [ -s ${OUTPUT_FILE} ]; then
        COUNT=$(grep -c "^>" ${OUTPUT_FILE})
        echo "  ${YEAR}: ${COUNT} sequences"
    else
        echo "  ${YEAR}: 0 sequences (file removed)"
        rm ${OUTPUT_FILE}
    fi
done

echo ""
echo "=========================================="
echo "✅ Year splitting complete for ${LINEAGE}"
echo "=========================================="
echo "Output directory: ${OUTPUT_DIR}"
echo ""
echo "Files created:"
ls -1 ${OUTPUT_DIR}/${LINEAGE}_*.fasta | wc -l
