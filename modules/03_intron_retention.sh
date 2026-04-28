#!/bin/bash
# ==============================================================================
# NeOcean Module 03: Intron Retention Detection (STAR + iREAD)
# Output: IR ORF FASTA for peptide screening
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../workflow/common.sh"

# ------------------------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------------------------
SAMPLE_DIR="${NEOCEAN_SAMPLE_DIR:-}"
LOG_DIR="${NEOCEAN_LOG_DIR:-${SAMPLE_DIR}/logs}"
STATUS_DIR="${NEOCEAN_STATUS_DIR:-${SAMPLE_DIR}/status}"

while [[ $# -gt 0 ]]; do
    case $1 in
        --sample-dir) SAMPLE_DIR="$2"; shift 2 ;;
        --log-dir) LOG_DIR="$2"; shift 2 ;;
        --status-dir) STATUS_DIR="$2"; shift 2 ;;
        *) shift ;;
    esac
done

[[ -z "$SAMPLE_DIR" ]] && { log_error "SAMPLE_DIR not set"; exit 1; }

# ------------------------------------------------------------------------------
# Resolve paths
# ------------------------------------------------------------------------------
SAMPLE_NAME="${NEOCEAN_SAMPLE_NAME}"
FQ1="${NEOCEAN_RNA_FQ1}"
FQ2="${NEOCEAN_RNA_FQ2}"
STAR_INDEX="${NEOCEAN_STAR_INDEX}"
GENOME_FASTA="${NEOCEAN_GENOME_FASTA}"
IREAD_DIR="${NEOCEAN_IREAD_DIR}"
IREAD_BED="${NEOCEAN_IREAD_BED}"
STAR_THREADS="${NEOCEAN_STAR_THREADS:-20}"
CONDA_DENOVO="${NEOCEAN_CONDA_DENOVO:-de_nove}"

IR_DIR="${SAMPLE_DIR}/intron_retention"
IR_STAR_DIR="${IR_DIR}/STAR"
IR_IREAD_DIR="${IR_DIR}/iREAD"
mkdir -p "$IR_STAR_DIR" "$IR_IREAD_DIR"

# ------------------------------------------------------------------------------
# Pre-flight checks
# ------------------------------------------------------------------------------
log_info "[03_intron_retention] Starting intron retention detection for sample: ${SAMPLE_NAME}"

check_tools STAR bedtools awk || exit 1
check_file "$FQ1" "RNA-seq R1" || exit 1
check_file "$FQ2" "RNA-seq R2" || exit 1
check_dir "$STAR_INDEX" "STAR index" || exit 1
check_file "$GENOME_FASTA" "Reference genome" || exit 1

if [[ -n "$IREAD_BED" ]] && [[ -f "$IREAD_BED" ]]; then
    log_info "[03_intron_retention] iREAD BED file: ${IREAD_BED}"
else
    log_warn "[03_intron_retention] iREAD BED file not found. IR detection may fail."
fi

# ------------------------------------------------------------------------------
# Step 1: STAR alignment for IR
# ------------------------------------------------------------------------------
IR_BAM="${IR_STAR_DIR}/${SAMPLE_NAME}Aligned.sortedByCoord.out.bam"
IR_LOG_FINAL="${IR_STAR_DIR}/${SAMPLE_NAME}Log.final.out"

if [[ -f "$IR_BAM" ]] && [[ -f "$IR_LOG_FINAL" ]]; then
    log_info "[03_intron_retention] IR BAM exists. Skipping STAR alignment."
else
    log_info "[03_intron_retention] Running STAR for intron retention..."
    
    STAR \
        --twopassMode Basic \
        --quantMode GeneCounts \
        --runThreadN "$STAR_THREADS" \
        --genomeDir "$STAR_INDEX" \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "${IR_STAR_DIR}/${SAMPLE_NAME}" \
        --readFilesCommand gunzip -c \
        --readFilesIn "$FQ1" "$FQ2" \
        > "${LOG_DIR}/intron_retention_star.log" 2>&1
    
    log_info "[03_intron_retention] STAR alignment completed."
fi

# ------------------------------------------------------------------------------
# Step 2: Extract uniquely mapped reads count
# ------------------------------------------------------------------------------
UNIQUE_READS_FILE="${IR_STAR_DIR}/unique_reads_count.txt"
if [[ -f "$IR_LOG_FINAL" ]]; then
    grep "Uniquely mapped reads number" "$IR_LOG_FINAL" | awk '{print $NF/2}' > "$UNIQUE_READS_FILE"
    log_info "[03_intron_retention] Unique reads count: $(cat "$UNIQUE_READS_FILE")"
else
    log_error "[03_intron_retention] STAR log final out not found."
    exit 1
fi

# ------------------------------------------------------------------------------
# Step 3: Run iREAD
# ------------------------------------------------------------------------------
IR_RESULTS="${IR_IREAD_DIR}/results"
mkdir -p "$IR_RESULTS"

IREAD_OUTPUT="${IR_RESULTS}/${SAMPLE_NAME}Aligned.sortedByCoord.out.ir.txt"

if [[ -f "$IREAD_OUTPUT" ]]; then
    log_info "[03_intron_retention] iREAD output exists. Skipping."
else
    if [[ -n "$IREAD_DIR" ]] && [[ -f "${IREAD_DIR}/iread.py" ]]; then
        log_info "[03_intron_retention] Running iREAD..."
        
        python "${IREAD_DIR}/iread.py" \
            "$IR_BAM" \
            "$IREAD_BED" \
            -o "$IR_RESULTS" \
            -t "$(cat "$UNIQUE_READS_FILE")" \
            > "${LOG_DIR}/intron_retention_iread.log" 2>&1
        
        log_info "[03_intron_retention] iREAD completed."
    else
        log_warn "[03_intron_retention] iREAD not available. Skipping."
    fi
fi

# ------------------------------------------------------------------------------
# Step 4: Extract retained introns and get sequences
# ------------------------------------------------------------------------------
IR_BED="${IR_RESULTS}/${SAMPLE_NAME}_yes.bed"
IR_FA="${IR_IREAD_DIR}/${SAMPLE_NAME}_intron_retention.fa"

if [[ -f "$IREAD_OUTPUT" ]]; then
    log_info "[03_intron_retention] Extracting retained introns..."
    
    # Extract 'yes' calls and convert to BED format
    # Note: chromosome name handling (with or without 'chr' prefix)
    awk '$7 == "yes" {gsub("-", "\t", $1); print "chr"$1}' "$IREAD_OUTPUT" > "$IR_BED"
    
    if [[ -s "$IR_BED" ]]; then
        bedtools getfasta \
            -fi "$GENOME_FASTA" \
            -bed "$IR_BED" \
            -fo "$IR_FA" \
            -name \
            > "${LOG_DIR}/intron_retention_bedtools.log" 2>&1
        
        log_info "[03_intron_retention] Retained intron sequences: ${IR_FA}"
    else
        log_warn "[03_intron_retention] No retained introns detected."
        touch "$IR_FA"
    fi
else
    log_warn "[03_intron_retention] iREAD output not found. Creating placeholder."
    touch "$IR_FA"
fi

# ------------------------------------------------------------------------------
# Step 5: ORF prediction
# ------------------------------------------------------------------------------
GETORF_OUT="${IR_IREAD_DIR}/${SAMPLE_NAME}_IR.getorf.fa"

if [[ -s "$IR_FA" ]]; then
    log_info "[03_intron_retention] Running getorf..."
    conda_run "$CONDA_DENOVO" getorf \
        -sequence "$IR_FA" \
        -outseq "$GETORF_OUT" \
        -find 1 \
        -minsize 24 \
        > "${LOG_DIR}/intron_retention_getorf.log" 2>&1
    log_info "[03_intron_retention] ORF prediction completed: ${GETORF_OUT}"
else
    log_warn "[03_intron_retention] Empty IR FASTA. Skipping getorf."
    touch "$GETORF_OUT"
fi

log_info "[03_intron_retention] Module completed."
exit 0
