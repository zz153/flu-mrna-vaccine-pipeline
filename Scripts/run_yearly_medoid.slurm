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
