#!/bin/bash -e
#SBATCH --job-name=yearly_medoid
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --output=logs/yearly_medoid.%j.out
#SBATCH --error=logs/yearly_medoid.%j.err

# =========================================================
#  Generic Yearly MEDOID Calculator
#  Usage:
#     sbatch run_yearly_medoid.slurm H1N1 h1n1
# =========================================================

# ---- Environment ----
ENV=/home/ranzo85p/conda_envs/bio
PY=$ENV/bin/python

PREFIX=${1:-H1N1}     # e.g. H1N1 or H3N2
OUT_DIR=${2:-h1n1}    # base folder with yearly subdirectories

mkdir -p logs

echo "=== Tool versions ==="
$PY --version
$PY -c "import Bio; print('Biopython', Bio.__version__)"
echo "=== Starting MEDOID for all years in ${OUT_DIR} ==="

# =========================================================
#  Loop through each year
# =========================================================
for year_dir in ${OUT_DIR}/*; do
    year=$(basename "$year_dir")
    aln="${year_dir}/${PREFIX}_${year}_aligned.fasta"
    designdir="${year_dir}/designs"
    logdir="${year_dir}/logs"

    mkdir -p "$designdir" "$logdir"

    echo "---------------------------------------------"
    echo "[CHECK] Year: ${year}"
    echo "[CHECK] Alignment: ${aln}"

    # --- verify alignment ---
    if [[ ! -s "$aln" ]]; then
        echo "[SKIP] Missing or empty alignment for ${year}"
        continue
    fi

    seqs=$(grep -c "^>" "$aln" || echo 0)
    if [[ $seqs -lt 2 ]]; then
        echo "[SKIP] Only ${seqs} sequences (need ≥ 2) for ${year}"
        continue
    fi

    # --- Calculate MEDOID ---
    echo "[INFO] Calculating MEDOID for ${year} (${seqs} seqs) ..."
    $PY <<EOF
from collections import Counter
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Read alignment
aln = AlignIO.read("${aln}", "fasta")
seq_ids = [r.id for r in aln]
seqs = [str(r.seq) for r in aln]

if not seqs:
    import sys
    sys.stderr.write("[ERROR] No sequences in alignment for ${year}\n")
    sys.exit(1)

# Hamming distance (ignoring gaps on either side)
def hamming_ungapped(a, b):
    return sum((x != y) and (x != '-') and (y != '-') for x, y in zip(a, b))

# ===== MEDOID =====
# Find sequence that minimizes sum of distances to all others
scores = [sum(hamming_ungapped(s, t) for t in seqs) for s in seqs]
medoid_idx = scores.index(min(scores))
medoid_seq = seqs[medoid_idx].replace('-', 'X')  # ungap
medoid_id = seq_ids[medoid_idx]

medoid_record = SeqRecord(
    Seq(medoid_seq),
    id="${PREFIX}_${year}_MEDOID",
    description=f"(from {medoid_id})"
)

# Write output
SeqIO.write([medoid_record], "${designdir}/${PREFIX}_${year}_MEDOID.fasta", "fasta")

print(f"[OK] MEDOID written to ${designdir}/${PREFIX}_${year}_MEDOID.fasta")
EOF

    if [[ $? -eq 0 ]]; then
        echo "[DONE] ${year} MEDOID complete ✅"
    else
        echo "[ERROR] Failed to calculate MEDOID for ${year}"
    fi
done

echo "=== All yearly MEDOID calculations complete 🚀 ==="
