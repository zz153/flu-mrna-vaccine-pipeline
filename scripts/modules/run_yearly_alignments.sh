# Multiple Alignments
# ============================

# --- Parameters ---
ENV_PATH="/home/$USER/conda_envs/bio"
INPUT_DIR=""
OUT_BASE="h1n1"
PREFIX="H1N1"
MIN_SEQS=6
CPUS=${SLURM_CPUS_PER_TASK:-8}

# --- Parse Arguments ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --env) ENV_PATH="$2"; shift 2;;
    --input-dir) INPUT_DIR="$2"; shift 2;;
    --out) OUT_BASE="$2"; shift 2;;
    --prefix) PREFIX="$2"; shift 2;;
    --min-seqs) MIN_SEQS="$2"; shift 2;;
    --cpus) CPUS="$2"; shift 2;;
    *) echo "[ERROR] Unknown argument: $1"; exit 1;;
  esac
done

if [[ -z "${INPUT_DIR}" ]]; then
  echo "[ERROR] You must specify --input-dir"
  exit 1
fi

# --- Tool Paths ---
MAFFT="${ENV_PATH}/bin/mafft"

# --- Function ---
log(){ echo "[$(date +'%F %T')] $*"; }

mkdir -p logs

# --- Run MAFFT for each yearly fasta ---
for fasta in ${INPUT_DIR}/${PREFIX}_*.fasta; do
  year=$(basename "$fasta" .fasta | awk -F '_' '{print $2}')
  outdir="${OUT_BASE}/${year}"
  mkdir -p "${outdir}/logs"

  log "Aligning ${PREFIX}_${year} ..."
  nseq=$(grep -c '^>' "$fasta" || echo 0)
  if (( nseq < MIN_SEQS )); then
    log "[SKIP] ${year} has only ${nseq} sequences (<${MIN_SEQS})"
    continue
  fi

  "$MAFFT" --thread ${CPUS} --auto "$fasta" \
    > "${outdir}/${PREFIX}_${year}_aligned.fasta" \
    2> "${outdir}/logs/mafft.log"

  log "Alignment complete → ${outdir}/${PREFIX}_${year}_aligned.fasta"
done

log "=== ALL ALIGNMENTS DONE 🚀 ==="
