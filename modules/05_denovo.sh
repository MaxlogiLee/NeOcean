#!/bin/bash
# ==============================================================================
# NeOcean Module 05: De Novo Transcript Assembly
# Steps: StringTie -> Cuffmerge -> gffread -> getorf
# Output: De novo ORF FASTA for peptide screening
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
BAM_SORTED="${SAMPLE_DIR}/star/${SAMPLE_NAME}.sorted.2passAligned.out.bam"
REF_GTF="${NEOCEAN_REF_GTF}"
GENOME_FASTA="${NEOCEAN_GENOME_FASTA}"
STRINGTIE_THREADS="${NEOCEAN_STRINGTIE_THREADS:-20}"
CONDA_CUFFLINKS="${NEOCEAN_CONDA_CUFFLINKS:-cufflinks_env}"
CONDA_DENOVO="${NEOCEAN_CONDA_DENOVO:-de_nove}"

DENOVO_DIR="${SAMPLE_DIR}/denovo"
mkdir -p "$DENOVO_DIR"

# ------------------------------------------------------------------------------
# Pre-flight checks
# ------------------------------------------------------------------------------
log_info "[05_denovo] Starting de novo assembly for sample: ${SAMPLE_NAME}"

check_tools stringtie || exit 1
check_file "$BAM_SORTED" "Sorted BAM" || exit 1
check_file "$REF_GTF" "Reference GTF" || exit 1
check_file "$GENOME_FASTA" "Reference genome" || exit 1

# ------------------------------------------------------------------------------
# Step 1: StringTie assembly (without reference guidance for novel transcripts)
# ------------------------------------------------------------------------------
PASS_GTF="${DENOVO_DIR}/${SAMPLE_NAME}_2pass.gtf"

if [[ -f "$PASS_GTF" ]]; then
    log_info "[05_denovo] StringTie GTF exists. Skipping."
else
    log_info "[05_denovo] Running StringTie de novo assembly..."
    stringtie "$BAM_SORTED" \
        -o "$PASS_GTF" \
        -m 100 \
        -c 1 \
        -p "$STRINGTIE_THREADS" \
        > "${LOG_DIR}/denovo_stringtie.log" 2>&1
    log_info "[05_denovo] StringTie assembly completed."
fi

# ------------------------------------------------------------------------------
# Step 2: Cuffmerge with reference annotation
# ------------------------------------------------------------------------------
GTF_LIST="${DENOVO_DIR}/${SAMPLE_NAME}_cuffmergegtf.list"
MERGED_DIR="${DENOVO_DIR}/merged_asm"
MERGED_GTF="${DENOVO_DIR}/${SAMPLE_NAME}_de_novo_merged_candidates.gtf"

echo "$PASS_GTF" > "$GTF_LIST"

if [[ -f "$MERGED_GTF" ]]; then
    log_info "[05_denovo] Cuffmerge output exists. Skipping."
else
    log_info "[05_denovo] Running Cuffmerge..."
    
    conda_run "$CONDA_CUFFLINKS" cuffmerge \
        -o "$MERGED_DIR" \
        -g "$REF_GTF" \
        "$GTF_LIST" \
        > "${LOG_DIR}/denovo_cuffmerge.log" 2>&1
    
    mv "${MERGED_DIR}/merged.gtf" "$MERGED_GTF"
    log_info "[05_denovo] Cuffmerge completed."
fi

# ------------------------------------------------------------------------------
# Step 3: Extract transcript sequences
# ------------------------------------------------------------------------------
TRANSCRIPT_FA="${DENOVO_DIR}/${SAMPLE_NAME}_transcript_de_novo.fa"

if [[ -f "$TRANSCRIPT_FA" ]]; then
    log_info "[05_denovo] Transcript FASTA exists. Skipping."
else
    log_info "[05_denovo] Extracting transcript sequences with gffread..."
    
    conda_run "$CONDA_CUFFLINKS" gffread \
        "$MERGED_GTF" \
        -g "$GENOME_FASTA" \
        -w "$TRANSCRIPT_FA" \
        > "${LOG_DIR}/denovo_gffread.log" 2>&1
    
    log_info "[05_denovo] Transcript sequences extracted."
fi

# ------------------------------------------------------------------------------
# Step 4: ORF prediction
# ------------------------------------------------------------------------------
GETORF_OUT="${DENOVO_DIR}/${SAMPLE_NAME}_transcript_de_novo_getorf.fa"

if [[ -f "$GETORF_OUT" ]]; then
    log_info "[05_denovo] ORF file exists. Skipping."
else
    log_info "[05_denovo] Running getorf..."
    conda_run "$CONDA_DENOVO" getorf \
        -find 1 \
        -minsize 24 \
        -sequence "$TRANSCRIPT_FA" \
        -outseq "$GETORF_OUT" \
        > "${LOG_DIR}/denovo_getorf.log" 2>&1
    log_info "[05_denovo] ORF prediction completed: ${GETORF_OUT}"
fi

log_info "[05_denovo] Module completed."
exit 0
