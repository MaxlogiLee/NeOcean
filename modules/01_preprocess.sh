#!/bin/bash
# ==============================================================================
# NeOcean Module 01: RNA-seq Preprocessing
# Steps: STAR 1-pass -> SJ index -> STAR 2-pass -> Sort/Index -> StringTie
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../workflow/common.sh"

# ------------------------------------------------------------------------------
# Parse arguments (also accept env vars from main.sh)
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
# Resolve paths from environment
# ------------------------------------------------------------------------------
SAMPLE_NAME="${NEOCEAN_SAMPLE_NAME}"
FQ1="${NEOCEAN_RNA_FQ1}"
FQ2="${NEOCEAN_RNA_FQ2}"
STAR_INDEX="${NEOCEAN_STAR_INDEX}"
GENOME_FASTA="${NEOCEAN_GENOME_FASTA}"
STAR_THREADS="${NEOCEAN_STAR_THREADS:-20}"
STRINGTIE_THREADS="${NEOCEAN_STRINGTIE_THREADS:-20}"

STAR_BAM_DIR="${SAMPLE_DIR}/star"
mkdir -p "$STAR_BAM_DIR"

# ------------------------------------------------------------------------------
# Pre-flight checks
# ------------------------------------------------------------------------------
log_info "[01_preprocess] Starting preprocessing for sample: ${SAMPLE_NAME}"

check_tools STAR samtools stringtie || exit 1
check_file "$FQ1" "RNA-seq R1" || exit 1
check_file "$FQ2" "RNA-seq R2" || exit 1
check_dir "$STAR_INDEX" "STAR index" || exit 1
check_file "$GENOME_FASTA" "Reference genome" || exit 1

# ------------------------------------------------------------------------------
# Step 1: STAR 1-pass (generate splice junctions)
# ------------------------------------------------------------------------------
SJ_TAB="${STAR_BAM_DIR}/${SAMPLE_NAME}.SJ.out.tab"
EXIT_1PASS="${STAR_BAM_DIR}/${SAMPLE_NAME}.1pass.exit"

if [[ -f "$SJ_TAB" ]] && [[ -f "$EXIT_1PASS" ]] && [[ "$(cat "$EXIT_1PASS")" == "0" ]]; then
    log_info "[01_preprocess] 1-pass SJ.out.tab exists and succeeded. Skipping."
else
    log_info "[01_preprocess] Running STAR 1-pass..."
    echo 1 > "$EXIT_1PASS"
    
    STAR \
        --genomeDir "$STAR_INDEX" \
        --readFilesType Fastx \
        --readFilesIn "$FQ1" "$FQ2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${STAR_BAM_DIR}/${SAMPLE_NAME}." \
        --outSAMtype None \
        --outSAMmode None \
        --runThreadN "$STAR_THREADS" \
        --outFilterMultimapScoreRange 1 \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 10 \
        --alignIntronMax 500000 \
        --alignMatesGapMax 1000000 \
        --sjdbScore 2 \
        --alignSJDBoverhangMin 1 \
        --genomeLoad NoSharedMemory \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterScoreMinOverLread 0.33 \
        --sjdbOverhang 100 \
        --outSAMattributes NH HI NM MD AS nM XS \
        --outSAMattrIHstart 0 \
        > "${STAR_BAM_DIR}/${SAMPLE_NAME}.1pass.log" 2>&1
    
    if [[ $? -eq 0 ]]; then
        log_info "[01_preprocess] STAR 1-pass completed successfully."
        echo 0 > "$EXIT_1PASS"
    else
        log_error "[01_preprocess] STAR 1-pass failed."
        echo 1 > "$EXIT_1PASS"
        exit 1
    fi
fi

# ------------------------------------------------------------------------------
# Step 2: Build sample-specific genome index from SJ
# ------------------------------------------------------------------------------
GENOMEDIR="${STAR_BAM_DIR}/genomedir"
SA_INDEX="${GENOMEDIR}/SAindex"

if [[ -f "$SA_INDEX" ]]; then
    log_info "[01_preprocess] Sample-specific genome index exists. Skipping."
else
    log_info "[01_preprocess] Building genome index from 1-pass SJ..."
    mkdir -p "$GENOMEDIR"
    
    if [[ ! -s "$SJ_TAB" ]]; then
        log_error "[01_preprocess] SJ.out.tab is missing or empty after 1-pass."
        exit 1
    fi
    
    STAR \
        --runMode genomeGenerate \
        --genomeDir "$GENOMEDIR" \
        --genomeFastaFiles "$GENOME_FASTA" \
        --sjdbOverhang 100 \
        --runThreadN "$STAR_THREADS" \
        --limitSjdbInsertNsj 7000000 \
        --sjdbFileChrStartEnd "$SJ_TAB" \
        > "${STAR_BAM_DIR}/${SAMPLE_NAME}.genomeGenerate.log" 2>&1
    
    if [[ $? -ne 0 ]]; then
        log_error "[01_preprocess] Genome index generation failed."
        exit 1
    fi
    log_info "[01_preprocess] Genome index built successfully."
fi

# ------------------------------------------------------------------------------
# Step 3: STAR 2-pass alignment
# ------------------------------------------------------------------------------
BAM_2PASS_UNSORTED="${STAR_BAM_DIR}/${SAMPLE_NAME}.unsorted.2passAligned.out.bam"
EXIT_2PASS="${STAR_BAM_DIR}/${SAMPLE_NAME}.2pass.exit"

if [[ -f "$BAM_2PASS_UNSORTED" ]] && [[ -f "$EXIT_2PASS" ]] && [[ "$(cat "$EXIT_2PASS")" == "0" ]]; then
    log_info "[01_preprocess] 2-pass BAM exists and succeeded. Skipping."
else
    log_info "[01_preprocess] Running STAR 2-pass..."
    echo 1 > "$EXIT_2PASS"
    
    STAR \
        --runThreadN "$STAR_THREADS" \
        --genomeDir "$GENOMEDIR" \
        --readFilesType Fastx \
        --readFilesIn "$FQ1" "$FQ2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${STAR_BAM_DIR}/${SAMPLE_NAME}.unsorted.2pass" \
        --outFilterMultimapScoreRange 1 \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 10 \
        --alignIntronMax 500000 \
        --alignMatesGapMax 1000000 \
        --sjdbScore 2 \
        --alignSJDBoverhangMin 1 \
        --genomeLoad NoSharedMemory \
        --limitBAMsortRAM 98765432100 \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterScoreMinOverLread 0.33 \
        --sjdbOverhang 100 \
        --outSAMunmapped Within \
        --outSAMtype BAM Unsorted \
        --outSAMattributes NH HI NM MD AS nM XS \
        --outSAMattrIHstart 0 \
        --limitSjdbInsertNsj 7000000 \
        --outFilterType BySJout \
        > "${STAR_BAM_DIR}/${SAMPLE_NAME}.2pass.log" 2>&1
    
    if [[ $? -eq 0 ]]; then
        log_info "[01_preprocess] STAR 2-pass completed successfully."
        echo 0 > "$EXIT_2PASS"
    else
        log_error "[01_preprocess] STAR 2-pass failed."
        echo 1 > "$EXIT_2PASS"
        exit 1
    fi
fi

# ------------------------------------------------------------------------------
# Step 4: Sort, index BAM; run StringTie
# ------------------------------------------------------------------------------
BAM_SORTED="${STAR_BAM_DIR}/${SAMPLE_NAME}.sorted.2passAligned.out.bam"
BAM_INDEX="${BAM_SORTED}.bai"
STRINGTIE_GTF="${STAR_BAM_DIR}/${SAMPLE_NAME}.stringtie.gtf"

if [[ -f "$STRINGTIE_GTF" ]] && [[ -f "$BAM_INDEX" ]]; then
    log_info "[01_preprocess] Sorted BAM and StringTie GTF exist. Skipping."
else
    log_info "[01_preprocess] Sorting and indexing BAM..."
    samtools sort -@ "$STAR_THREADS" -o "$BAM_SORTED" "$BAM_2PASS_UNSORTED"
    samtools index "$BAM_SORTED"
    
    log_info "[01_preprocess] Running StringTie..."
    stringtie "$BAM_SORTED" \
        -o "$STRINGTIE_GTF" \
        --rf \
        -m 100 \
        -c 1 \
        -p "$STRINGTIE_THREADS" \
        -l "$SAMPLE_NAME" \
        -v > "${STAR_BAM_DIR}/${SAMPLE_NAME}.stringtie.log" 2>&1
    
    # Cleanup unsorted BAM to save space
    if [[ -f "$STRINGTIE_GTF" ]]; then
        rm -f "$BAM_2PASS_UNSORTED"
        rm -rf "${STAR_BAM_DIR}"/*_STARtmp
        log_info "[01_preprocess] Cleanup completed."
    else
        log_error "[01_preprocess] StringTie failed to produce GTF."
        exit 1
    fi
fi

# ------------------------------------------------------------------------------
# Verify outputs
# ------------------------------------------------------------------------------
OUTPUT_FILES=(
    "$BAM_SORTED"
    "$BAM_INDEX"
    "$STRINGTIE_GTF"
    "$SJ_TAB"
)

check_files "${OUTPUT_FILES[@]}" || exit 1

log_info "[01_preprocess] All preprocessing outputs verified."
exit 0
