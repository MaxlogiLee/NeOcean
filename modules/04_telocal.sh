#!/bin/bash
# ==============================================================================
# NeOcean Module 04: Transposable Element Expression (TElocal)
# Output: TE ORF FASTA for peptide screening
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
REF_GTF_TELOCAL="${NEOCEAN_REF_GTF_TELOCAL}"
TE_INDEX="${NEOCEAN_TE_INDEX}"
TE_LOCATIONS="${NEOCEAN_TE_LOCATIONS}"
GENOME_FASTA="${NEOCEAN_GENOME_FASTA}"
CONDA_TELOCAL="${NEOCEAN_CONDA_TELOCAL:-telocal-xly}"
CONDA_DENOVO="${NEOCEAN_CONDA_DENOVO:-de_nove}"

TELOCAL_DIR="${SAMPLE_DIR}/telocal"
mkdir -p "$TELOCAL_DIR"

# ------------------------------------------------------------------------------
# Pre-flight checks
# ------------------------------------------------------------------------------
log_info "[04_telocal] Starting TElocal analysis for sample: ${SAMPLE_NAME}"

check_tools bedtools awk grep || exit 1
check_file "$BAM_SORTED" "Sorted BAM" || exit 1
check_file "$REF_GTF_TELOCAL" "TElocal GTF" || exit 1
check_file "$TE_INDEX" "TE index" || exit 1
check_file "$GENOME_FASTA" "Reference genome" || exit 1

if [[ -n "$TE_LOCATIONS" ]] && [[ -f "$TE_LOCATIONS" ]]; then
    log_info "[04_telocal] TE locations file: ${TE_LOCATIONS}"
else
    log_warn "[04_telocal] TE locations file not found. Coordinate extraction may fail."
fi

# ------------------------------------------------------------------------------
# Step 1: Run TElocal
# ------------------------------------------------------------------------------
CNT_TABLE="${TELOCAL_DIR}/${SAMPLE_NAME}.cntTable"

if [[ -f "$CNT_TABLE" ]]; then
    log_info "[04_telocal] TElocal count table exists. Skipping."
else
    log_info "[04_telocal] Running TElocal..."
    
    if command -v TElocal >/dev/null 2>&1; then
        TElocal \
            -b "$BAM_SORTED" \
            --GTF "$REF_GTF_TELOCAL" \
            --TE "$TE_INDEX" \
            --sortByPos \
            --project "${TELOCAL_DIR}/${SAMPLE_NAME}" \
            > "${LOG_DIR}/telocal_quant.log" 2>&1
        
        log_info "[04_telocal] TElocal quantification completed."
    else
        log_warn "[04_telocal] TElocal command not found in PATH."
        log_warn "[04_telocal] If using conda, ensure environment '${CONDA_TELOCAL}' is configured."
        
        # Try conda fallback
        if command -v conda >/dev/null 2>&1; then
            log_info "[04_telocal] Attempting conda run with ${CONDA_TELOCAL}..."
            conda_run "$CONDA_TELOCAL" TElocal \
                -b "$BAM_SORTED" \
                --GTF "$REF_GTF_TELOCAL" \
                --TE "$TE_INDEX" \
                --sortByPos \
                --project "${TELOCAL_DIR}/${SAMPLE_NAME}" \
                > "${LOG_DIR}/telocal_quant.log" 2>&1
            log_info "[04_telocal] TElocal quantification completed (conda)."
        else
            log_error "[04_telocal] TElocal not available. Cannot proceed."
            exit 1
        fi
    fi
fi

# ------------------------------------------------------------------------------
# Step 2: Extract TE counts
# ------------------------------------------------------------------------------
TE_TXT="${TELOCAL_DIR}/${SAMPLE_NAME}_TE.txt"
GENES_TXT="${TELOCAL_DIR}/${SAMPLE_NAME}_genes.txt"

if [[ -f "$CNT_TABLE" ]]; then
    log_info "[04_telocal] Extracting TE and gene counts..."
    
    # Extract TE part (contains ":" in ID)
    awk 'NR==1{h=$0} NR>1 && /:/{print (f++?"":h) $0}' "$CNT_TABLE" > "$TE_TXT"
    
    # Extract gene part (no ":" in ID)
    awk 'NR==1{h=$0} NR>1 && !/:/ {print (g++?"":h) $0}' "$CNT_TABLE" > "$GENES_TXT"
    
    log_info "[04_telocal] TE entries: $(wc -l < "$TE_TXT"), Gene entries: $(wc -l < "$GENES_TXT")"
else
    log_warn "[04_telocal] Count table not found. Creating placeholders."
    touch "$TE_TXT" "$GENES_TXT"
fi

# ------------------------------------------------------------------------------
# Step 3: Filter and extract TE locations
# ------------------------------------------------------------------------------
TE_FILTERED="${TELOCAL_DIR}/${SAMPLE_NAME}_filtered_10_TE.txt"
TE_NAMES="${TELOCAL_DIR}/${SAMPLE_NAME}_TE_name.txt"
TE_LOCATIONS_OUT="${TELOCAL_DIR}/${SAMPLE_NAME}_TE_name_location.txt"
TE_BED="${TELOCAL_DIR}/${SAMPLE_NAME}_TE_name_location.bed"
TE_BED_NEW="${TELOCAL_DIR}/${SAMPLE_NAME}_TE_name_location_new.bed"
TE_FA_RAW="${TELOCAL_DIR}/${SAMPLE_NAME}_TE_name_location_new.fa"
TE_FA="${TELOCAL_DIR}/${SAMPLE_NAME}_TE.fa"

if [[ -s "$TE_TXT" ]]; then
    log_info "[04_telocal] Filtering TE counts > 10..."
    awk '$2 > 10' "$TE_TXT" > "$TE_FILTERED"
    awk -F: '{print $1}' "$TE_FILTERED" > "$TE_NAMES"
    
    if [[ -n "$TE_LOCATIONS" ]] && [[ -f "$TE_LOCATIONS" ]] && [[ -s "$TE_NAMES" ]]; then
        log_info "[04_telocal] Mapping TE names to genomic locations..."
        
        # Extract IDs (skip header)
        tail -n +2 "$TE_NAMES" > "${TELOCAL_DIR}/${SAMPLE_NAME}_TE_name.list"
        
        # Batch grep
        grep -w -F -f "${TELOCAL_DIR}/${SAMPLE_NAME}_TE_name.list" "$TE_LOCATIONS" > "$TE_LOCATIONS_OUT"
        
        # Convert to BED
        awk 'BEGIN {OFS="\t"} {split($2, a, "[:-]"); print a[1], a[2], a[3], $1, 0, a[4]}' "$TE_LOCATIONS_OUT" > "$TE_BED"
        
        # Fix empty strand
        awk 'BEGIN {OFS="\t"} $6 == "" {$6="-"} {print}' "$TE_BED" > "$TE_BED_NEW"
        
        # Extract sequences
        bedtools getfasta \
            -fi "$GENOME_FASTA" \
            -bed "$TE_BED_NEW" \
            -fo "$TE_FA_RAW" \
            > "${LOG_DIR}/telocal_bedtools.log" 2>&1
        
        # Rename headers to TE names
        awk 'NR==FNR {header[FNR]=$4; next} /^>/ {sub(/>.*/, ">" header[++count])} 1' \
            "$TE_BED_NEW" "$TE_FA_RAW" > "$TE_FA"
        
        log_info "[04_telocal] TE sequences extracted: ${TE_FA}"
    else
        log_warn "[04_telocal] Cannot map TE locations. Creating placeholder."
        touch "$TE_FA"
    fi
else
    log_warn "[04_telocal] Empty TE count table. Creating placeholder."
    touch "$TE_FA"
fi

# ------------------------------------------------------------------------------
# Step 4: ORF prediction
# ------------------------------------------------------------------------------
GETORF_OUT="${TELOCAL_DIR}/${SAMPLE_NAME}_TElocal.getorf.fa"

if [[ -s "$TE_FA" ]]; then
    log_info "[04_telocal] Running getorf..."
    conda_run "$CONDA_DENOVO" getorf \
        -find 1 \
        -minsize 24 \
        -sequence "$TE_FA" \
        -outseq "$GETORF_OUT" \
        > "${LOG_DIR}/telocal_getorf.log" 2>&1
    log_info "[04_telocal] ORF prediction completed: ${GETORF_OUT}"
else
    log_warn "[04_telocal] Empty TE FASTA. Skipping getorf."
    touch "$GETORF_OUT"
fi

log_info "[04_telocal] Module completed."
exit 0
