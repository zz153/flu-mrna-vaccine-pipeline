#!/bin/bash -e
#SBATCH --job-name=yearly_asr
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --output=logs/yearly_asr.%j.out
#SBATCH --error=logs/yearly_asr.%j.err

# =========================================================
#  Generic Yearly ASR Runner (IQ-TREE 2 + Biopython)
#  Usage:
#     sbatch run_yearly_asr_generic.slurm H1N1 h1n1
# =========================================================

# ---- Environment ----
ENV=/home/ranzo85p/conda_envs/bio
PY=$ENV/bin/python
IQT=$ENV/bin/iqtree2

PREFIX=${1:-H1N1}     # e.g. H1N1 or H3N2
OUT_DIR=${2:-h1n1}    # base folder with yearly subdirectories

mkdir -p logs

echo "=== Tool versions ==="
$PY --version
$PY -c "import Bio; print('Biopython', Bio.__version__)"
$IQT --version
echo "=== Starting ASR for all years in ${OUT_DIR} ==="

# =========================================================
#  Loop through each year
# =========================================================
for year_dir in ${OUT_DIR}/*; do
    year=$(basename "$year_dir")
    aln="${year_dir}/${PREFIX}_${year}_aligned.fasta"
    asrdir="${year_dir}/designs/asr"
    designdir="${year_dir}/designs"
    logdir="${year_dir}/logs"

    # Create directories FIRST before any tee/redirect operations
    mkdir -p "$asrdir" "$designdir" "$logdir"

    echo "---------------------------------------------"
    echo "[CHECK] Year: ${year}"
    echo "[CHECK] Alignment: ${aln}"

    # --- verify alignment ---
    if [[ ! -s "$aln" ]]; then
        echo "[SKIP] Missing or empty alignment for ${year}"
        continue
    fi

    seqs=$(grep -c "^>" "$aln" || echo 0)
    if [[ $seqs -lt 4 ]]; then
        echo "[SKIP] Only ${seqs} sequences (need ≥ 4) for ${year}"
        continue
    fi

    # --- Run ASR ---
    echo "[INFO] Running ASR for ${year} (${seqs} seqs) ..."
    aln_abs=$(realpath "$aln")
    log_abs=$(realpath "${logdir}/iqtree_asr_${year}.log")
    cd "$asrdir"
    $IQT -s "$aln_abs" -m LG+G -asr -nt 16 -pre iqtree_asr_${year} --redo \
        > >(tee "$log_abs") 2>&1
    cd - >/dev/null

    STATE="${asrdir}/iqtree_asr_${year}.state"
    if [[ ! -f "$STATE" ]]; then
        echo "[ERROR] ASR failed for ${year} (no .state file)"
        continue
    fi

    # --- Extract root sequence ---
    echo "[INFO] Extracting root sequence for ${year}"
    $PY <<EOF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

root_seq, root_node = [], None
with open("${STATE}") as f:
    for line in f:
        if line.startswith("#") or line.startswith("Node\t"): continue
        parts = line.strip().split("\t")
        if len(parts) < 3: continue
        node, site, state = parts[0], parts[1], parts[2]
        if root_node is None: root_node = node
        if node == root_node:
            root_seq.append("X" if state == "-" else state)

if root_seq:
    out_fa = "${designdir}/${PREFIX}_${year}_ASR.fasta"
    SeqIO.write([SeqRecord(Seq("".join(root_seq)), id="${PREFIX}_${year}_ASR", description="")],
                out_fa, "fasta")
    print(f"[OK] Root sequence written to {out_fa}")
else:
    import sys; sys.stderr.write(f"[ERROR] Could not reconstruct root sequence for ${year}\n")
    sys.exit(1)
EOF

    echo "[DONE] ${year} ASR + root extraction complete ✅"
done

echo "=== All yearly ASR tasks complete 🚀 ==="

