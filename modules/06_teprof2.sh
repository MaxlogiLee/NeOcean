#!/bin/bash
# ==============================================================================
# NeOcean Module 06: TEProf2 - Transposable Element Profiling 2
# Full 12-step pipeline for TE-derived transcript identification
# Output: TEProf2 ORF FASTA for peptide screening
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
STAR_BAM_DIR="${SAMPLE_DIR}/star"
STRINGTIE_GTF="${STAR_BAM_DIR}/${SAMPLE_NAME}.stringtie.gtf"
BAM_SORTED="${STAR_BAM_DIR}/${SAMPLE_NAME}.sorted.2passAligned.out.bam"
BAM_INDEX="${BAM_SORTED}.bai"
REF_GTF="${NEOCEAN_REF_GTF}"

TEPROF2_BIN="${NEOCEAN_TEPROF2_BIN}"
TEPROF2_DOC="${NEOCEAN_TEPROF2_DOC}"
TEPROF2_ARGS="${NEOCEAN_TEPROF2_ARGS}"
CONDA_TEPROF2="${NEOCEAN_CONDA_TEPROF2:-teprof2}"
CONDA_CUFFLINKS="${NEOCEAN_CONDA_CUFFLINKS:-cufflinks_env}"
CONDA_DENOVO="${NEOCEAN_CONDA_DENOVO:-de_nove}"

TEPROF2_DIR="${SAMPLE_DIR}/teprof2"
mkdir -p "$TEPROF2_DIR"
mkdir -p "${TEPROF2_DIR}/filterreadstats"

# ------------------------------------------------------------------------------
# Pre-flight checks
# ------------------------------------------------------------------------------
log_info "[06_teprof2] Starting TEProf2 analysis for sample: ${SAMPLE_NAME}"

check_tools parallel samtools bedtools || exit 1
check_file "$STRINGTIE_GTF" "StringTie GTF" || exit 1
check_file "$BAM_SORTED" "Sorted BAM" || exit 1
check_file "$BAM_INDEX" "BAM index" || exit 1
check_file "$REF_GTF" "Reference GTF" || exit 1

if [[ -n "$TEPROF2_BIN" ]] && [[ -d "$TEPROF2_BIN" ]]; then
    log_info "[06_teprof2] TEProf2 bin directory: ${TEPROF2_BIN}"
else
    log_warn "[06_teprof2] TEProf2 bin directory not configured."
fi

# ------------------------------------------------------------------------------
# Setup input symlinks
# ------------------------------------------------------------------------------
ln -sf "$STRINGTIE_GTF" "${TEPROF2_DIR}/${SAMPLE_NAME}.stringtie.gtf"
ln -sf "$BAM_SORTED" "${TEPROF2_DIR}/${SAMPLE_NAME}.stringtie.bam"
ln -sf "$BAM_INDEX" "${TEPROF2_DIR}/${SAMPLE_NAME}.stringtie.bam.bai"

cd "$TEPROF2_DIR" || exit 1

# ------------------------------------------------------------------------------
# Step 1: Annotate GTF
# ------------------------------------------------------------------------------
STEP1_OUT="${SAMPLE_NAME}.stringtie.gtf_annotated_filtered_test_all"

if [[ -f "$STEP1_OUT" ]]; then
    log_info "[06_teprof2] Step 1 annotation exists. Skipping."
else
    log_info "[06_teprof2] Step 1: Annotating GTF..."
    conda_run "$CONDA_TEPROF2" python2 \
        "${TEPROF2_BIN}/rmskhg38_annotate_gtf_update_test_tpm.py" \
        "${SAMPLE_NAME}.stringtie.gtf" \
        "$TEPROF2_ARGS" \
        > "${LOG_DIR}/teprof2_step1.log" 2>&1
    log_info "[06_teprof2] Step 1 completed."
fi

# ------------------------------------------------------------------------------
# Step 2: Process annotation TPM
# ------------------------------------------------------------------------------
STEP2_OUT="${STEP1_OUT}_c"

if [[ -f "$STEP2_OUT" ]]; then
    log_info "[06_teprof2] Step 2 TPM processing exists. Skipping."
else
    log_info "[06_teprof2] Step 2: Processing annotation TPM..."
    conda_run "$CONDA_TEPROF2" python \
        "${TEPROF2_BIN}/annotationtpmprocess.py" \
        "$STEP1_OUT" \
        > "${LOG_DIR}/teprof2_step2.log" 2>&1
    log_info "[06_teprof2] Step 2 completed."
fi

# ------------------------------------------------------------------------------
# Step 3: Aggregate processed annotations
# ------------------------------------------------------------------------------
STEP3_OUT="filter_combined_candidates.tsv"

if [[ -f "$STEP3_OUT" ]]; then
    log_info "[06_teprof2] Step 3 aggregation exists. Skipping."
else
    log_info "[06_teprof2] Step 3: Aggregating annotations..."
    conda_run "$CONDA_TEPROF2" \
        "${TEPROF2_BIN}/aggregateProcessedAnnotation.R" \
        -a "$TEPROF2_ARGS" \
        -e T \
        -d "$TEPROF2_DIR" \
        > "${LOG_DIR}/teprof2_step3.log" 2>&1
    log_info "[06_teprof2] Step 3 completed."
fi

# ------------------------------------------------------------------------------
# Step 4: Calculate read information
# ------------------------------------------------------------------------------
STEP4_CMDS="filterreadcommands.txt"
STEP4_STATS="filter_read_stats.txt"

if [[ -f "$STEP4_STATS" ]]; then
    log_info "[06_teprof2] Step 4 read stats exist. Skipping."
else
    log_info "[06_teprof2] Step 4: Calculating read information..."
    conda_run "$CONDA_TEPROF2" python \
        "${TEPROF2_BIN}/commandsmax_speed.py" \
        "$STEP3_OUT" \
        ./ \
        > "${LOG_DIR}/teprof2_step4.log" 2>&1
    
    if [[ -f "$STEP4_CMDS" ]]; then
        conda_run "$CONDA_TEPROF2" parallel -j 10 < "$STEP4_CMDS" \
            > "${LOG_DIR}/teprof2_step4_parallel.log" 2>&1
        
        find ./filterreadstats/ -name "*.stats" -type f -maxdepth 1 -print0 | \
            xargs -0 -n128 -P1 grep e > ./resultgrep_filterreadstatsdone.txt
        cat ./resultgrep_filterreadstatsdone.txt | sed 's/\:/\t/g' > "$STEP4_STATS"
        log_info "[06_teprof2] Step 4 completed."
    else
        log_warn "[06_teprof2] Step 4: filterreadcommands.txt not generated."
    fi
fi

# ------------------------------------------------------------------------------
# Step 5: Filter candidates
# ------------------------------------------------------------------------------
STEP5_OUT="candidate_transcripts.gff3"

if [[ -f "$STEP5_OUT" ]]; then
    log_info "[06_teprof2] Step 5 filtering exists. Skipping."
else
    log_info "[06_teprof2] Step 5: Filtering candidates..."
    conda_run "$CONDA_TEPROF2" \
        "${TEPROF2_DOC}/filterReadCandidates.R" \
        > "${LOG_DIR}/teprof2_step5.log" 2>&1
    rm -f Step4.RData
    log_info "[06_teprof2] Step 5 completed."
fi

# ------------------------------------------------------------------------------
# Step 6: Merge with reference GTF
# ------------------------------------------------------------------------------
STEP6_MERGED="reference_merged_candidates.gtf"

if [[ -f "$STEP6_MERGED" ]]; then
    log_info "[06_teprof2] Step 6 merge exists. Skipping."
else
    log_info "[06_teprof2] Step 6: Merging with reference GTF..."
    
    conda_run "$CONDA_TEPROF2" gffread \
        -E candidate_transcripts.gff3 \
        -T -o candidate_transcripts.gtf \
        > "${LOG_DIR}/teprof2_step6a.log" 2>&1
    
    echo candidate_transcripts.gtf > cuffmergegtf.list
    
    conda_run "$CONDA_CUFFLINKS" cuffmerge \
        -o ./merged_asm_full \
        -g "$REF_GTF" \
        cuffmergegtf.list \
        > "${LOG_DIR}/teprof2_step6b.log" 2>&1
    
    mv ./merged_asm_full/merged.gtf "$STEP6_MERGED"
    conda_run "$CONDA_TEPROF2" gffread \
        -E "$STEP6_MERGED" \
        -o- > reference_merged_candidates.gff3 \
        > "${LOG_DIR}/teprof2_step6c.log" 2>&1
    
    log_info "[06_teprof2] Step 6 completed."
fi

# ------------------------------------------------------------------------------
# Step 7: Annotate merged GTF
# ------------------------------------------------------------------------------
STEP7_OUT="reference_merged_candidates.gff3_annotated"

if [[ -f "$STEP7_OUT" ]]; then
    log_info "[06_teprof2] Step 7 re-annotation exists. Skipping."
else
    log_info "[06_teprof2] Step 7: Annotating merged GTF..."
    conda_run "$CONDA_TEPROF2" python2.7 \
        "${TEPROF2_BIN}/rmskhg38_annotate_gtf_update_test_tpm_cuff.py" \
        reference_merged_candidates.gff3 \
        "$TEPROF2_ARGS" \
        > "${LOG_DIR}/teprof2_step7.log" 2>&1
    log_info "[06_teprof2] Step 7 completed."
fi

# ------------------------------------------------------------------------------
# Step 8: Transcript-level quantification
# ------------------------------------------------------------------------------
STEP8_CMDS="quantificationCommands2.txt"

if [[ -f "$STEP8_CMDS" ]]; then
    log_info "[06_teprof2] Step 8 quantification commands exist. Skipping generation."
else
    log_info "[06_teprof2] Step 8: Generating quantification commands..."
    find ./ -maxdepth 1 -name "*bam" | while read file; do
        xbase=${file##*/}
        echo "samtools view -q 255 -h \"$file\" | stringtie - -o \"${xbase%.*}.gtf\" -e -b \"${xbase%.*}_stats\" -p 2 -m 100 -c 1 -G reference_merged_candidates.gtf" >> "$STEP8_CMDS"
    done
fi

if [[ -f "$STEP8_CMDS" ]]; then
    conda_run "$CONDA_TEPROF2" parallel -j 20 < "$STEP8_CMDS" \
        > "${LOG_DIR}/teprof2_step8.log" 2>&1
    log_info "[06_teprof2] Step 8 completed."
fi

# ------------------------------------------------------------------------------
# Step 9: Process expression output
# ------------------------------------------------------------------------------
log_info "[06_teprof2] Step 9: Processing expression output..."

conda_run "$CONDA_TEPROF2" \
    "${TEPROF2_DOC}/mergeAnnotationProcess.R" \
    > "${LOG_DIR}/teprof2_step9.log" 2>&1

# 9.1 Intron coverage
find . -maxdepth 2 -name "*i_data.ctab" > ctab_i.txt
if [[ -s ctab_i.txt ]]; then
    cat ctab_i.txt | while read ID; do
        fileid=$(echo "$ID" | awk -F "/" '{print $2}')
        cat <(printf 'chr\tstrand\tstart\tend\t'${fileid/_stats/}'\n') \
            <(grep -F -f candidate_introns.txt "$ID" | awk -F'\t' '{ print $2"\t"$3"\t"$4"\t"$5"\t"$6 }') \
            > "${ID}_cand"
    done
    
    cat <(find . -maxdepth 2 -name "*i_data.ctab_cand" | head -1 | while read file; do cat "$file" | awk '{print $1"\t"$2"\t"$3"\t"$4}'; done;) > table_i_all
    find . -maxdepth 2 -name "*i_data.ctab_cand" | while read file; do
        paste -d'\t' <(cat table_i_all) <(cat "$file" | awk '{print $5}') > table_i_all_temp
        mv -f table_i_all_temp table_i_all
    done
fi

# 9.2 Transcript-level expression
ls ./*stats/t_data.ctab > ctablist.txt 2>/dev/null || true
if [[ -s ctablist.txt ]]; then
    cat ctablist.txt | while read file; do
        echo "${TEPROF2_BIN}/stringtieExpressionFrac.py $file" >> stringtieExpressionFracCommands.txt
    done
    
    conda_run "$CONDA_TEPROF2" parallel -j 20 < stringtieExpressionFracCommands.txt \
        > "${LOG_DIR}/teprof2_step9b.log" 2>&1
    
    # Aggregate stats
    ls ./*stats/t_data.ctab_frac_tot > ctab_frac_tot_files.txt 2>/dev/null || true
    ls ./*stats/t_data.ctab_tpm > ctab_tpm_files.txt 2>/dev/null || true
    
    if [[ -s ctab_frac_tot_files.txt ]]; then
        cat <(echo "TranscriptID") <(find . -name "*ctab_frac_tot" | head -1 | while read file; do sort "$file" | awk '{print $1}'; done;) > table_frac_tot
        cat ctab_frac_tot_files.txt | while read file; do
            fileid=$(echo "$file" | awk -F "/" '{print $2}')
            paste -d'\t' <(cat table_frac_tot) <(cat <(echo ${fileid/_stats/}) <(sort "$file" | awk '{print $2}')) > table_frac_tot_temp
            mv table_frac_tot_temp table_frac_tot
        done
    fi
    
    if [[ -s ctab_tpm_files.txt ]]; then
        cat <(echo "TranscriptID") <(find . -name "*ctab_tpm" | head -1 | while read file; do sort "$file" | awk '{print $1}'; done;) > table_tpm
        cat ctab_tpm_files.txt | while read file; do
            fileid=$(echo "$file" | awk -F "/" '{print $2}')
            paste -d'\t' <(cat table_tpm) <(cat <(echo ${fileid/_stats/}) <(sort "$file" | awk '{print $2}')) > table_tpm_temp
            mv table_tpm_temp table_tpm
        done
    fi
    
    # Filter for candidates
    if [[ -f candidate_names.txt ]]; then
        cat <(head -1 table_frac_tot) <(grep -Ff candidate_names.txt table_frac_tot) > table_frac_tot_cand
        cat <(head -1 table_tpm) <(grep -Ff candidate_names.txt table_tpm) > table_tpm_cand
    fi
fi

log_info "[06_teprof2] Step 9 completed."

# ------------------------------------------------------------------------------
# Step 10: Final statistics
# ------------------------------------------------------------------------------
STEP10_DONE="step10.done"

if [[ -f "$STEP10_DONE" ]]; then
    log_info "[06_teprof2] Step 10 statistics exist. Skipping."
else
    log_info "[06_teprof2] Step 10: Generating final statistics..."
    conda_run "$CONDA_TEPROF2" \
        "${TEPROF2_DOC}/finalStatisticsOutput_new.R" \
        -e T \
        -a "$TEPROF2_ARGS" \
        > "${LOG_DIR}/teprof2_step10.log" 2>&1
    rm -f Step10.RData
    touch "$STEP10_DONE"
    log_info "[06_teprof2] Step 10 completed."
fi

# ------------------------------------------------------------------------------
# Step 11: Translate transcripts (Kozak method)
# ------------------------------------------------------------------------------
CANDIDATES_FA="candidates.fa"

if [[ -f "$CANDIDATES_FA" ]]; then
    log_info "[06_teprof2] Step 11 translation exists. Skipping."
else
    log_info "[06_teprof2] Step 11: Translating transcripts with Kozak method..."
    conda_run "$CONDA_TEPROF2" \
        "${TEPROF2_DOC}/translationPart1.R" \
        -g 'BSgenome.Hsapiens.UCSC.hg38' \
        > "${LOG_DIR}/teprof2_step11.log" 2>&1
    log_info "[06_teprof2] Step 11 completed."
fi

# ------------------------------------------------------------------------------
# Step 12: ORF prediction
# ------------------------------------------------------------------------------
GETORF_OUT="${TEPROF2_DIR}/${SAMPLE_NAME}_TEprof2.getorf.fa"

if [[ -f "$GETORF_OUT" ]]; then
    log_info "[06_teprof2] Step 12 ORF exists. Skipping."
else
    log_info "[06_teprof2] Step 12: Predicting ORFs..."
    conda_run "$CONDA_DENOVO" getorf \
        -find 1 \
        -minsize 24 \
        -sequence "$CANDIDATES_FA" \
        -outseq "$GETORF_OUT" \
        > "${LOG_DIR}/teprof2_step12.log" 2>&1
    log_info "[06_teprof2] Step 12 completed: ${GETORF_OUT}"
fi

cd - > /dev/null || true

log_info "[06_teprof2] All 12 steps completed."
exit 0
