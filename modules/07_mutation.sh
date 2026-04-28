#!/bin/bash
# ==============================================================================
# NeOcean Module 07: Somatic Mutation -> Neoantigen Peptides
# Steps: HLA typing (seq2HLA) -> Mutect2 -> VEP -> pVACseq -> mutation.fa
# Supports: paired tumor-normal & tumor-only (PON) modes
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
# Resolve paths from environment
# ------------------------------------------------------------------------------
SAMPLE_NAME="${NEOCEAN_SAMPLE_NAME}"
SEQUENCING_TYPE="${NEOCEAN_SEQUENCING_TYPE:-WGS}"
TUMOR_BAM="${NEOCEAN_TUMOR_BAM:-}"
NORMAL_BAM="${NEOCEAN_NORMAL_BAM:-}"
PON_VCF="${NEOCEAN_PON_VCF:-}"
HLA_FILE="${NEOCEAN_HLA_FILE:-}"
HLA_ALLELES="${NEOCEAN_HLA_ALLELES:-}"
RNA_FQ1="${NEOCEAN_RNA_FQ1}"
RNA_FQ2="${NEOCEAN_RNA_FQ2}"
GENOME_FASTA="${NEOCEAN_GENOME_FASTA}"
THREADS="${NEOCEAN_THREADS:-40}"

# Mutation-specific resources
SEQ2HLA="${NEOCEAN_MUTATION_SEQ2HLA:-}"
KNOWN_SITES="${NEOCEAN_MUTATION_DBSNP:-}"
GNOMAD="${NEOCEAN_MUTATION_GNOMAD:-}"
VEP_CACHE="${NEOCEAN_MUTATION_VEP_CACHE:-}"
VEP_PLUGINS="${NEOCEAN_MUTATION_VEP_PLUGINS:-}"
IEDB_DIR="${NEOCEAN_MUTATION_IEDB_DIR:-}"
CONDA_SEQ2HLA="${NEOCEAN_CONDA_SEQ2HLA:-seq2hla}"

MUTATION_DIR="${SAMPLE_DIR}/mutation"
mkdir -p "${MUTATION_DIR}"/{hla,vcf,vep,pvacseq,fa}

# Global variable to pass HLA alleles between steps
HLA_STRING=""

# ------------------------------------------------------------------------------
# Pre-flight checks
# ------------------------------------------------------------------------------
log_info "[07_mutation] Starting mutation analysis for sample: ${SAMPLE_NAME}"

check_tools gatk bcftools samtools tabix || exit 1
check_file "$GENOME_FASTA" "Reference genome" || exit 1

# ------------------------------------------------------------------------------
# Step 1: HLA Typing (seq2HLA from tumor RNA-seq)
# ------------------------------------------------------------------------------
run_hla_typing() {
    log_info "[07_mutation] Step 1: HLA typing via seq2HLA..."

    if [[ ! -f "$RNA_FQ1" || ! -f "$RNA_FQ2" ]]; then
        log_warn "[07_mutation] RNA-seq FASTQ not found. HLA typing skipped."
        return 1
    fi

    if [[ -z "$SEQ2HLA" || ! -f "$SEQ2HLA" ]]; then
        log_warn "[07_mutation] seq2HLA script not configured. Skipping HLA typing."
        return 1
    fi

    local hla_out="${MUTATION_DIR}/hla/${SAMPLE_NAME}"
    local hla_result="${hla_out}-ClassI.HLAgenotype4digits"

    if [[ -f "$hla_result" ]]; then
        log_info "[07_mutation] seq2HLA output exists. Parsing alleles..."
    else
        log_info "[07_mutation] Running seq2HLA..."
        conda_run "$CONDA_SEQ2HLA" python2 "$SEQ2HLA" \
            -1 "$RNA_FQ1" \
            -2 "$RNA_FQ2" \
            -p "$THREADS" \
            -r "$hla_out" \
            > "${LOG_DIR}/mutation_hla.log" 2>&1
    fi

    if [[ ! -f "$hla_result" ]]; then
        log_error "[07_mutation] seq2HLA output not found: ${hla_result}"
        return 1
    fi

    # Parse Class I HLA alleles (4-digit resolution), filter NA
    HLA_STRING=$(awk -F'\t' 'NR>1 && NF>=3 && $2!="NA" && $3!="NA" {
        print "HLA-"$1"*"gensub(/:/, ":", "g", $2);
        print "HLA-"$1"*"gensub(/:/, ":", "g", $3)
    }' "$hla_result" | sort -u | paste -sd "," -)

    if [[ -z "$HLA_STRING" ]]; then
        log_error "[07_mutation] Failed to parse HLA alleles from seq2HLA output."
        return 1
    fi

    log_info "[07_mutation] HLA type (seq2HLA): ${HLA_STRING}"
}

# Try seq2HLA first; if it fails, fall back to config-provided HLA
run_hla_typing || {
    if [[ -n "$HLA_ALLELES" ]]; then
        HLA_STRING="$HLA_ALLELES"
        log_info "[07_mutation] Using HLA alleles from config: ${HLA_STRING}"
    elif [[ -n "$HLA_FILE" && -f "$HLA_FILE" ]]; then
        HLA_STRING=$(tr '\n' ',' < "$HLA_FILE" | sed 's/,$//')
        log_info "[07_mutation] Using HLA alleles from file: ${HLA_STRING}"
    else
        log_warn "[07_mutation] HLA typing failed and no fallback HLA provided."
        log_warn "[07_mutation] pVACseq will fail unless HLA alleles are supplied."
    fi
}

# ------------------------------------------------------------------------------
# Step 2: DNA Preprocessing (optional; skip if BAM already exists)
# ------------------------------------------------------------------------------
# NOTE: DNA preprocessing is disabled by default. If tumor/normal DNA FASTQs
# are available and BAMs do not exist, uncomment the lines below and define
# TUMOR_DNA_R1, TUMOR_DNA_R2, NORMAL_DNA_R1, NORMAL_DNA_R2 in the config.
# run_dna_preprocessing() {
#     local sample_id=$1; local r1=$2; local r2=$3
#     local out_bam="${SAMPLE_DIR}/preprocessed/${sample_id}.final.bam"
#     if [[ -f "$out_bam" ]]; then
#         log_info "[07_mutation] ${sample_id} BAM exists, skipping DNA preprocessing."
#         return 0
#     fi
#     log_info "[07_mutation] Preprocessing DNA for ${sample_id}..."
#     mkdir -p "${SAMPLE_DIR}/preprocessed"
#     local rg="@RG\tID:${sample_id}\tSM:${sample_id}\tLB:WGS\tPL:ILLUMINA\tPU:unit1"
#     bwa mem -t "$THREADS" -R "$rg" "$GENOME_FASTA" "$r1" "$r2" | \
#         samtools view -@ "$THREADS" -Sb - | \
#         samtools sort -@ "$THREADS" -o "${SAMPLE_DIR}/preprocessed/${sample_id}.sorted.bam" -
#     samtools index "${SAMPLE_DIR}/preprocessed/${sample_id}.sorted.bam"
#     gatk MarkDuplicates \
#         -I "${SAMPLE_DIR}/preprocessed/${sample_id}.sorted.bam" \
#         -O "${SAMPLE_DIR}/preprocessed/${sample_id}.marked.bam" \
#         -M "${SAMPLE_DIR}/preprocessed/${sample_id}.markdup.metrics" \
#         --REMOVE_DUPLICATES false
#     samtools index "${SAMPLE_DIR}/preprocessed/${sample_id}.marked.bam"
#     gatk BaseRecalibrator -R "$GENOME_FASTA" \
#         -I "${SAMPLE_DIR}/preprocessed/${sample_id}.marked.bam" \
#         --known-sites "$KNOWN_SITES" \
#         -O "${SAMPLE_DIR}/preprocessed/${sample_id}.recal.table"
#     gatk ApplyBQSR -R "$GENOME_FASTA" \
#         -I "${SAMPLE_DIR}/preprocessed/${sample_id}.marked.bam" \
#         --bqsr-recal-file "${SAMPLE_DIR}/preprocessed/${sample_id}.recal.table" \
#         -O "$out_bam"
#     samtools index "$out_bam"
#     rm -f "${SAMPLE_DIR}/preprocessed/${sample_id}.sorted.bam"
#     log_info "[07_mutation] ${sample_id} DNA preprocessing completed."
# }

# ------------------------------------------------------------------------------
# Step 3: Somatic Mutation Detection (Mutect2)
# ------------------------------------------------------------------------------
run_mutect2() {
    log_info "[07_mutation] Step 3: Somatic mutation detection (Mutect2)..."

    local vcf_raw="${MUTATION_DIR}/vcf/somatic.unfiltered.vcf.gz"
    local vcf_filtered="${MUTATION_DIR}/vcf/somatic.filtered.vcf.gz"
    local vcf_passed="${MUTATION_DIR}/vcf/somatic.PASS.vcf.gz"
    local f1r2_stats="${MUTATION_DIR}/vcf/f1r2.tar.gz"
    local rom_model="${MUTATION_DIR}/vcf/rom.tar.gz"

    if [[ -f "$vcf_passed" ]]; then
        log_info "[07_mutation] Mutect2 PASS VCF exists. Skipping."
        return 0
    fi

    if [[ ! -f "$TUMOR_BAM" ]]; then
        log_error "[07_mutation] Tumor BAM missing: ${TUMOR_BAM}"
        exit 1
    fi

    # Mutect2: paired vs tumor-only
    if [[ -n "${NORMAL_BAM}" && -f "${NORMAL_BAM}" ]]; then
        log_info "[07_mutation] Mutect2: tumor-normal paired mode"
        gatk Mutect2 \
            -R "$GENOME_FASTA" \
            -I "$TUMOR_BAM" \
            -I "$NORMAL_BAM" \
            --normal "$(basename "$NORMAL_BAM" .final.bam)" \
            -O "$vcf_raw" \
            --germline-resource "$GNOMAD" \
            --f1r2-tar-gz "$f1r2_stats" \
            > "${LOG_DIR}/mutation_mutect2.log" 2>&1
    else
        log_info "[07_mutation] Mutect2: tumor-only mode"
        local pon_arg=""
        [[ -n "$PON_VCF" && -f "$PON_VCF" ]] && pon_arg="--panel-of-normals ${PON_VCF}"

        gatk Mutect2 \
            -R "$GENOME_FASTA" \
            -I "$TUMOR_BAM" \
            -O "$vcf_raw" \
            --germline-resource "$GNOMAD" \
            ${pon_arg:+${pon_arg}} \
            --f1r2-tar-gz "$f1r2_stats" \
            > "${LOG_DIR}/mutation_mutect2.log" 2>&1
    fi

    # Filter
    gatk LearnReadOrientationModel -I "$f1r2_stats" -O "$rom_model" \
        > "${LOG_DIR}/mutation_rom.log" 2>&1

    gatk FilterMutectCalls \
        -V "$vcf_raw" -R "$GENOME_FASTA" -O "$vcf_filtered" \
        --ob-priors "$rom_model" \
        > "${LOG_DIR}/mutation_filter.log" 2>&1

    # Extract PASS
    bcftools view -f PASS -O z -o "$vcf_passed" "$vcf_filtered" \
        > "${LOG_DIR}/mutation_pass.log" 2>&1
    bcftools index "$vcf_passed"

    local n_pass
    n_pass=$(bcftools view -H "$vcf_passed" | wc -l)
    log_info "[07_mutation] PASS variants: ${n_pass}"
}

# ------------------------------------------------------------------------------
# Step 4: VEP Annotation
# ------------------------------------------------------------------------------
run_vep() {
    log_info "[07_mutation] Step 4: VEP annotation..."

    local input_vcf="${MUTATION_DIR}/vcf/somatic.PASS.vcf.gz"
    local output_vcf="${MUTATION_DIR}/vep/${SAMPLE_NAME}.vep.ann.vcf.gz"

    if [[ -f "$output_vcf" ]]; then
        log_info "[07_mutation] VEP output exists. Skipping."
        return 0
    fi

    if [[ ! -f "$input_vcf" ]]; then
        log_error "[07_mutation] Input VCF missing: ${input_vcf}"
        exit 1
    fi

    check_tool vep || exit 1

    vep -i "$input_vcf" \
        -o "$output_vcf" \
        --vcf \
        --compress_output bgzip \
        --cache \
        --assembly GRCh38 \
        --dir_cache "$VEP_CACHE" \
        --fasta "$GENOME_FASTA" \
        --pick \
        --transcript_version \
        --symbol \
        --terms SO \
        --plugin Wildtype,"${VEP_PLUGINS}" \
        --plugin Downstream,"${VEP_PLUGINS}" \
        --fork "$THREADS" \
        > "${LOG_DIR}/mutation_vep.log" 2>&1

    tabix -p vcf "$output_vcf"
    log_info "[07_mutation] VEP completed: ${output_vcf}"
}

# ------------------------------------------------------------------------------
# Step 5: pVACseq Neoantigen Prediction
# ------------------------------------------------------------------------------
run_pvacseq() {
    log_info "[07_mutation] Step 5: pVACseq prediction..."

    if [[ -z "${HLA_STRING}" ]]; then
        log_error "[07_mutation] HLA_STRING undefined. HLA typing may have failed."
        exit 1
    fi

    local vep_vcf="${MUTATION_DIR}/vep/${SAMPLE_NAME}.vep.ann.vcf.gz"
    local pvac_out="${MUTATION_DIR}/pvacseq"

    if [[ ! -f "$vep_vcf" ]]; then
        log_error "[07_mutation] VEP VCF missing: ${vep_vcf}"
        exit 1
    fi

    check_tool pvacseq || exit 1

    pvacseq run \
        "$vep_vcf" \
        "$SAMPLE_NAME" \
        "$HLA_STRING" \
        NetMHCpan MHCflurry \
        "$pvac_out" \
        -e1 8,9,10,11 \
        --pass-only \
        --binding-threshold 500 \
        --top-score-metric median \
        --tdna-vaf 0.05 \
        --normal-vaf 0.02 \
        --tdna-cov 10 \
        --normal-cov 5 \
        --maximum-transcript-support-level 5 \
        --netmhc-stab \
        --run-reference-proteome-similarity \
        --iedb-install-directory "$IEDB_DIR" \
        --keep-tmp-files \
        > "${LOG_DIR}/mutation_pvacseq.log" 2>&1

    log_info "[07_mutation] pVACseq completed."
}

# ------------------------------------------------------------------------------
# Step 6: Extract mutation.fa for downstream Peptide Screening
# ------------------------------------------------------------------------------
extract_mutation_fa() {
    log_info "[07_mutation] Step 6: Extracting mutation.fa..."

    local filtered_tsv="${MUTATION_DIR}/pvacseq/MHC_Class_I/${SAMPLE_NAME}/${SAMPLE_NAME}.filtered.tsv"
    local mutation_fa="${MUTATION_DIR}/fa/${SAMPLE_NAME}_mutation.fa"
    local mutation_final="${SAMPLE_DIR}/mutation/mutation.fa"

    if [[ ! -f "$filtered_tsv" ]]; then
        log_error "[07_mutation] pVACseq filtered.tsv not found: ${filtered_tsv}"
        exit 1
    fi

    # Extract MT Epitope Seq (column 13), de-duplicate, output FASTA
    # Header: gene_chrom_pos_ref_alt for traceability
    awk -F'\t' 'NR>1 && $24=="True" {
        seq = $13;
        if (seq != "" && !seen[seq]++) {
            print ">"$1"_"$2"_"$3"_"$4"_"$5;
            print seq
        }
    }' "$filtered_tsv" > "$mutation_fa"

    # Symlink to peptide screening expected location
    mkdir -p "$(dirname "$mutation_final")"
    ln -sf "$mutation_fa" "$mutation_final"

    local count
    count=$(grep -c "^>" "$mutation_fa" || echo 0)
    log_info "[07_mutation] mutation.fa generated: ${count} peptides."
}

# ------------------------------------------------------------------------------
# Main execution
# ------------------------------------------------------------------------------

run_mutect2
run_vep
run_pvacseq
extract_mutation_fa

log_info "[07_mutation] All steps completed successfully."
exit 0
