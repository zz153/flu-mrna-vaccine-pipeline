#!/bin/bash

# Script 03: Multiple Sequence Alignment with MAFFT
#
# Aligns sequences for both clustered and unclustered tracks
# Uses MAFFT with automatic algorithm selection

set -e

# --- Parameters ---
LINEAGE=$1      # H1N1, H3N2, or VicB
TRACK=$2        # clustered or unclustered
THREADS=4

if [ -z "$LINEAGE" ] || [ -z "$TRACK" ]; then
    echo "Usage: bash scripts/03_align_sequences.sh [LINEAGE] [TRACK]"
    echo "Example: bash scripts/03_align_sequences.sh H1N1 clustered"
    echo "         bash scripts/03_align_sequences.sh H1N1 unclustered"
    echo ""
    echo "TRACK options: clustered, unclustered"
    exit 1
fi

echo "=========================================="
echo "Aligning ${LINEAGE} (${TRACK} track)"
echo "=========================================="

# --- Set input/output based on track ---
PROC_DIR="data/processed/${LINEAGE}"
RESULTS_DIR="results/${TRACK}"

if [ "$TRACK" == "clustered" ]; then
    INPUT="${PROC_DIR}/${LINEAGE}_cdhit99.fasta"
elif [ "$TRACK" == "unclustered" ]; then
    INPUT="${PROC_DIR}/${LINEAGE}_clean.fasta"
else
    echo "Error: TRACK must be 'clustered' or 'unclustered'"
    exit 1
fi

OUTPUT="${RESULTS_DIR}/alignments/${LINEAGE}_aligned.fasta"

# Check input exists
if [ ! -f "$INPUT" ]; then
    echo "Error: Input file not found: $INPUT"
    exit 1
fi

mkdir -p ${RESULTS_DIR}/alignments

# Count sequences
SEQ_COUNT=$(grep -c "^>" ${INPUT})
echo "Input sequences: ${SEQ_COUNT}"

# --- Run MAFFT ---
echo "Running MAFFT alignment..."
echo "This may take several minutes..."

mafft --auto --thread ${THREADS} ${INPUT} > ${OUTPUT}

# Verify output
if [ -f "$OUTPUT" ]; then
    ALIGNED_COUNT=$(grep -c "^>" ${OUTPUT})
    echo ""
    echo "=========================================="
    echo "✅ Alignment complete for ${LINEAGE} (${TRACK})"
    echo "=========================================="
    echo "Sequences aligned: ${ALIGNED_COUNT}"
    echo "Output: ${OUTPUT}"
else
    echo "❌ Error: Alignment failed!"
    exit 1
fi
