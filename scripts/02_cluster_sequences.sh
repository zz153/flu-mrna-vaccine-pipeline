#!/bin/bash

# Script 02: Cluster sequences at 99% identity using CD-HIT
#
# Purpose: Reduce near-identical sequences to representative sequences
# This removes sequences that are >99% similar, keeping only diverse representatives

set -e

# --- Parameters ---
LINEAGE=$1  # H1N1, H3N2, or VicB
IDENTITY=0.99
THREADS=4

if [ -z "$LINEAGE" ]; then
    echo "Usage: bash scripts/02_cluster_sequences.sh [LINEAGE]"
    echo "Example: bash scripts/02_cluster_sequences.sh H1N1"
    exit 1
fi

echo "=========================================="
echo "Clustering ${LINEAGE} at ${IDENTITY} identity"
echo "=========================================="

# --- Directories ---
PROC_DIR="data/processed/${LINEAGE}"
INPUT="${PROC_DIR}/${LINEAGE}_clean.fasta"
OUTPUT="${PROC_DIR}/${LINEAGE}_cdhit99.fasta"

# --- Run CD-HIT ---
cd-hit -i ${INPUT} \
    -o ${OUTPUT} \
    -c ${IDENTITY} \
    -n 5 \
    -M 8000 \
    -T ${THREADS} \
    -d 0

# --- Count sequences ---
BEFORE=$(grep -c "^>" ${INPUT})
AFTER=$(grep -c "^>" ${OUTPUT})
REMOVED=$((BEFORE - AFTER))

echo ""
echo "=========================================="
echo "✅ Clustering complete for ${LINEAGE}"
echo "=========================================="
echo "Before clustering:  ${BEFORE}"
echo "After clustering:   ${AFTER}"
echo "Sequences removed:  ${REMOVED}"
echo ""
echo "Output: ${OUTPUT}"
