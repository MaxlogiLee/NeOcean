#!/bin/bash
# ==============================================================================
# NeOcean Module 02: Fusion Gene Detection (SOAPfuse)
# Output: Fusion ORF FASTA for peptide screening
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
SOAPFUSE_PROG="${NEOCEAN_SOAPFUSE_PROG}"
SOAPFUSE_CONFIG="${NEOCEAN_SOAPFUSE_CONFIG}"
CONDA_DENOVO="${NEOCEAN_CONDA_DENOVO:-de_nove}"

FUSION_DIR="${SAMPLE_DIR}/fusion"
mkdir -p "$FUSION_DIR"

# ------------------------------------------------------------------------------
# Pre-flight checks
# ------------------------------------------------------------------------------
log_info "[02_fusion] Starting fusion detection for sample: ${SAMPLE_NAME}"

check_tools perl samtools || exit 1
check_file "$FQ1" "RNA-seq R1" || exit 1
check_file "$FQ2" "RNA-seq R2" || exit 1

if [[ -n "$SOAPFUSE_PROG" ]] && [[ -d "$SOAPFUSE_PROG" ]]; then
    log_info "[02_fusion] SOAPfuse program directory: ${SOAPFUSE_PROG}"
else
    log_warn "[02_fusion] SOAPfuse program directory not configured. Fusion detection may be limited."
fi

# ------------------------------------------------------------------------------
# Step 1: Prepare SOAPfuse input structure
# ------------------------------------------------------------------------------
FUSION_DATA_DIR="${FUSION_DIR}/data/${SAMPLE_NAME}/Lib-${SAMPLE_NAME}"
mkdir -p "$FUSION_DATA_DIR"

# Create symlinks to clean data
ln -sf "$FQ1" "${FUSION_DATA_DIR}/Run-${SAMPLE_NAME}_1.fq.gz"
ln -sf "$FQ2" "${FUSION_DATA_DIR}/Run-${SAMPLE_NAME}_2.fq.gz"

# Create sample list
LIST_FILE="${FUSION_DIR}/${SAMPLE_NAME}.soapfuse.list"
printf '%s\t%s\t%s\t%s\n' \
    "$SAMPLE_NAME" \
    "Lib-${SAMPLE_NAME}" \
    "Run-${SAMPLE_NAME}" \
    "150" \
    > "$LIST_FILE"

# ------------------------------------------------------------------------------
# Step 2: Run SOAPfuse
# ------------------------------------------------------------------------------
SOAPFUSE_DONE="${FUSION_DIR}/soapfuse.done"

if [[ -f "$SOAPFUSE_DONE" ]]; then
    log_info "[02_fusion] SOAPfuse output exists. Skipping."
else
    if [[ -n "$SOAPFUSE_PROG" ]] && [[ -f "${SOAPFUSE_PROG}/SOAPfuse-RUN.pl" ]]; then
        log_info "[02_fusion] Running SOAPfuse..."
        
        export PERL5LIB="${SOAPFUSE_PROG}/source/bin/perl_module/:${PERL5LIB:-}"
        
        perl "${SOAPFUSE_PROG}/SOAPfuse-RUN.pl" \
            -c "$SOAPFUSE_CONFIG" \
            -fd "${FUSION_DIR}/data/" \
            -l "$LIST_FILE" \
            -o "$FUSION_DIR" \
            > "${LOG_DIR}/fusion_soapfuse.log" 2>&1
        
        touch "$SOAPFUSE_DONE"
        log_info "[02_fusion] SOAPfuse completed."
    else
        log_warn "[02_fusion] SOAPfuse not available. Skipping fusion detection."
        log_warn "[02_fusion] If fusion results exist from previous run, peptide extraction will still proceed."
    fi
fi

# ------------------------------------------------------------------------------
# Step 3: Extract fusion peptide sequences
# ------------------------------------------------------------------------------
PEPTIDE_CHAIN="${FUSION_DIR}/final_fusion_genes/${SAMPLE_NAME}/analysis/For_peptides_analysis/${SAMPLE_NAME}.trans.fusion.peptide.chain"
FUSION_FA="${FUSION_DIR}/final_fusion_genes/${SAMPLE_NAME}/analysis/For_peptides_analysis/${SAMPLE_NAME}_fusion.fa"
FUSION_FA_CLEAN="${FUSION_DIR}/final_fusion_genes/${SAMPLE_NAME}/analysis/For_peptides_analysis/${SAMPLE_NAME}_fusion_clean.fa"

if [[ -f "$PEPTIDE_CHAIN" ]]; then
    log_info "[02_fusion] Extracting fusion peptide sequences..."
    
    while IFS= read -r line; do
        col1=$(echo "$line" | awk '{print $1}')
        col2=$(echo "$line" | awk '{print $2}')
        echo ">$col1" >> "$FUSION_FA"
        echo "$col2" >> "$FUSION_FA"
    done < "$PEPTIDE_CHAIN"
    
    # Remove pipe characters
    sed 's/|//g' "$FUSION_FA" > "$FUSION_FA_CLEAN"
    
    log_info "[02_fusion] Fusion FASTA generated: ${FUSION_FA_CLEAN}"
else
    log_warn "[02_fusion] Peptide chain file not found: ${PEPTIDE_CHAIN}"
    log_warn "[02_fusion] Creating empty placeholder for downstream compatibility."
    touch "$FUSION_FA_CLEAN"
fi

# ------------------------------------------------------------------------------
# Step 4: ORF prediction
# ------------------------------------------------------------------------------
GETORF_OUT="${FUSION_DIR}/final_fusion_genes/${SAMPLE_NAME}/analysis/For_peptides_analysis/${SAMPLE_NAME}_fusion.getorf.fa"

if [[ -s "$FUSION_FA_CLEAN" ]]; then
    log_info "[02_fusion] Running getorf on fusion sequences..."
    conda_run "$CONDA_DENOVO" getorf \
        -sequence "$FUSION_FA_CLEAN" \
        -outseq "$GETORF_OUT" \
        -find 1 \
        -minsize 24 \
        > "${LOG_DIR}/fusion_getorf.log" 2>&1
    log_info "[02_fusion] ORF prediction completed: ${GETORF_OUT}"
else
    log_warn "[02_fusion] Empty fusion FASTA. Skipping getorf."
    touch "$GETORF_OUT"
fi

# ------------------------------------------------------------------------------
# Verify outputs
# ------------------------------------------------------------------------------
if [[ ! -f "$GETORF_OUT" ]]; then
    log_warn "[02_fusion] Final ORF file not produced. Creating empty placeholder."
    touch "$GETORF_OUT"
fi

log_info "[02_fusion] Module completed."
exit 0
