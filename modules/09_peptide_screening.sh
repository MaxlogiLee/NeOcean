#!/bin/bash
# ==============================================================================
# NeOcean Module 09: Peptide Screening & Mass Spectrometry Integration
# Integrates MS data with predicted neoantigen sequences from all upstream modules.
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../workflow/common.sh"

# ------------------------------------------------------------------------------
# Tool script directories (bundled with the pipeline code)
# ------------------------------------------------------------------------------
TOOLS_R_DIR="${SCRIPT_DIR}/../tools/r"
TOOLS_PY_DIR="${SCRIPT_DIR}/../tools/python"

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
GREP_CODE_DIR="${NEOCEAN_GREP_CODE_DIR}"
BN_FASTA="${NEOCEAN_BN_FASTA}"
UNIPROT_FASTA="${NEOCEAN_UNIPROT_FASTA}"
NCBI_FASTA="${NEOCEAN_NCBI_FASTA}"

GREP_DIR="${SAMPLE_DIR}/grep"
GREP_DATA="${GREP_DIR}/data"
GREP_RES="${GREP_DIR}/res"
GREP_RES_SAMPLE="${GREP_RES}/${SAMPLE_NAME}"
mkdir -p "${GREP_DIR}/code" "${GREP_DATA}/MSdata" "${GREP_DATA}/Decoy" \
    "${GREP_RES}/peptide" "${GREP_RES_SAMPLE}"

# Upstream module outputs
FUSION_ORF="${SAMPLE_DIR}/fusion/final_fusion_genes/${SAMPLE_NAME}/analysis/For_peptides_analysis/${SAMPLE_NAME}_fusion.getorf.fa"
IR_ORF="${SAMPLE_DIR}/intron_retention/iREAD/${SAMPLE_NAME}_IR.getorf.fa"
TELOCAL_ORF="${SAMPLE_DIR}/telocal/${SAMPLE_NAME}_TElocal.getorf.fa"
DENOVO_ORF="${SAMPLE_DIR}/denovo/${SAMPLE_NAME}_transcript_de_novo_getorf.fa"
TEPROF2_ORF="${SAMPLE_DIR}/teprof2/${SAMPLE_NAME}_TEprof2.getorf.fa"
MUTATION_FA="${SAMPLE_DIR}/mutation/mutation.fa"

MSDATA_XLSX="${GREP_DATA}/MSdata/${SAMPLE_NAME}.xlsx"

# ------------------------------------------------------------------------------
# Helper: Link upstream ORF with fallback
# ------------------------------------------------------------------------------
link_orf() {
    local src="$1"
    local dst_dir="$2"
    local dst_name="$3"
    
    mkdir -p "$dst_dir"
    if [[ -f "$src" ]] && [[ -s "$src" ]]; then
        ln -sf "$src" "${dst_dir}/${dst_name}"
        log_info "[09_peptide_screening] Linked: ${src} -> ${dst_dir}/${dst_name}"
    else
        log_warn "[09_peptide_screening] Missing ORF: ${src}"
        log_warn "[09_peptide_screening] Creating empty placeholder: ${dst_dir}/${dst_name}"
        touch "${dst_dir}/${dst_name}"
    fi
}

# ------------------------------------------------------------------------------
# Pre-flight checks
# ------------------------------------------------------------------------------
log_info "[09_peptide_screening] Starting peptide screening for sample: ${SAMPLE_NAME}"

if [[ ! -f "$MSDATA_XLSX" ]]; then
    log_warn "[09_peptide_screening] MS data not found: ${MSDATA_XLSX}"
    log_warn "[09_peptide_screening] Please place mass spectrometry results at the above path."
fi

# ------------------------------------------------------------------------------
# Step 0: Prepare directory structure and link upstream ORFs
# ------------------------------------------------------------------------------
log_info "[09_peptide_screening] Step 0: Preparing directories and linking ORFs..."

mkdir -p "${GREP_DATA}/1.mutation" "${GREP_DATA}/2.Fusion" "${GREP_DATA}/3.IR" \
    "${GREP_DATA}/4.TE_local" "${GREP_DATA}/5.TEprof2" "${GREP_DATA}/6.Denovo" \
    "${GREP_DATA}/BN" "${GREP_DATA}/NCBI" "${GREP_DATA}/Uniport"

link_orf "$MUTATION_FA" "${GREP_DATA}/1.mutation" "mutation.fa"
link_orf "$FUSION_ORF" "${GREP_DATA}/2.Fusion" "${SAMPLE_NAME}_fusion.getorf.fa"
link_orf "$IR_ORF" "${GREP_DATA}/3.IR" "${SAMPLE_NAME}_IR.getorf.fa"
link_orf "$TELOCAL_ORF" "${GREP_DATA}/4.TE_local" "${SAMPLE_NAME}_TElocal.getorf.fa"
link_orf "$TEPROF2_ORF" "${GREP_DATA}/5.TEprof2" "${SAMPLE_NAME}_TEprof2.getorf.fa"
link_orf "$DENOVO_ORF" "${GREP_DATA}/6.Denovo" "${SAMPLE_NAME}_de_novo.getorf.fa"

# Link reference databases
if [[ -n "$BN_FASTA" ]] && [[ -f "$BN_FASTA" ]]; then
    ln -sf "$BN_FASTA" "${GREP_DATA}/BN/BN_processed.fa"
else
    log_warn "[09_peptide_screening] BN database not found."
fi

if [[ -n "$UNIPROT_FASTA" ]] && [[ -f "$UNIPROT_FASTA" ]]; then
    ln -sf "$UNIPROT_FASTA" "${GREP_DATA}/Uniport/uniprot_protein.fa"
else
    log_warn "[09_peptide_screening] UniProt database not found."
fi

if [[ -n "$NCBI_FASTA" ]] && [[ -f "$NCBI_FASTA" ]]; then
    ln -sf "$NCBI_FASTA" "${GREP_DATA}/NCBI/ncbi_protein.fa"
else
    log_warn "[09_peptide_screening] NCBI database not found."
fi

cd "$GREP_DIR" || exit 1

# ------------------------------------------------------------------------------
# Step 1: Process mass spectrometry results
# ------------------------------------------------------------------------------
log_info "[09_peptide_screening] Step 1: Processing MS results..."

if [[ -f "$MSDATA_XLSX" ]]; then
    if [[ -f "${TOOLS_R_DIR}/ms_results.R" ]]; then
        # Determine I/L transform mode from sample config (default TRUE)
        DO_IL="TRUE"
        if yaml_exists "$NEOCEAN_SAMPLE_CONFIG" "sample.ms.no_il_transform" 2>/dev/null; then
            DO_IL="FALSE"
        fi
        
        Rscript "${TOOLS_R_DIR}/ms_results.R" \
            "$SAMPLE_NAME" \
            "$MSDATA_XLSX" \
            "$GREP_DIR" \
            "$DO_IL" \
            > "${LOG_DIR}/peptide_ms_results.log" 2>&1
        log_info "[09_peptide_screening] MS results processed."
    else
        log_warn "[09_peptide_screening] ms_results.R not found at ${TOOLS_R_DIR}"
    fi
else
    log_warn "[09_peptide_screening] Skipping MS processing (no MS data)."
fi

# ------------------------------------------------------------------------------
# Step 2: Parallel peptide grep against all databases
# ------------------------------------------------------------------------------
log_info "[09_peptide_screening] Step 2: Parallel peptide grep..."

declare -A DB_MAP=(
    ["BN"]="${GREP_DATA}/BN/BN_processed.fa"
    ["NCBI"]="${GREP_DATA}/NCBI/ncbi_protein.fa"
    ["Uniport"]="${GREP_DATA}/Uniport/uniprot_protein.fa"
    ["Mutation"]="${GREP_DATA}/1.mutation/mutation.fa"
    ["Fusion"]="${GREP_DATA}/2.Fusion/${SAMPLE_NAME}_fusion.getorf.fa"
    ["IR"]="${GREP_DATA}/3.IR/${SAMPLE_NAME}_IR.getorf.fa"
    ["TElocal"]="${GREP_DATA}/4.TE_local/${SAMPLE_NAME}_TElocal.getorf.fa"
    ["TEprof2"]="${GREP_DATA}/5.TEprof2/${SAMPLE_NAME}_TEprof2.getorf.fa"
    ["DN"]="${GREP_DATA}/6.Denovo/${SAMPLE_NAME}_de_novo.getorf.fa"
)

if [[ -f "${TOOLS_R_DIR}/peptide_grep.R" ]]; then
    for DB_NAME in BN NCBI Uniport Mutation Fusion IR TElocal TEprof2 DN; do
        DB_PATH="${DB_MAP[$DB_NAME]}"
        
        if [[ ! -f "$DB_PATH" ]]; then
            log_warn "[09_peptide_screening] Skipping ${DB_NAME} (not found: ${DB_PATH})"
            continue
        fi
        
        log_info "[09_peptide_screening]   -> Grep against ${DB_NAME}..."
        
        nohup Rscript "${TOOLS_R_DIR}/peptide_grep.R" \
            "$SAMPLE_NAME" \
            "${GREP_RES}/peptide/${SAMPLE_NAME}_peptide_list.txt" \
            "$DB_NAME" \
            "$DB_PATH" \
            FALSE \
            "$GREP_RES_SAMPLE" \
            > "${GREP_RES}/peptide/${DB_NAME}_grep.log" 2>&1 &
    done
    
    wait
    log_info "[09_peptide_screening] Step 2 completed."
else
    log_warn "[09_peptide_screening] peptide_grep.R not found. Skipping Step 2."
fi

# ------------------------------------------------------------------------------
# Step 3: Coordinate mapping (Denovo, TEprof2, IR, TElocal)
# ------------------------------------------------------------------------------
log_info "[09_peptide_screening] Step 3: Coordinate mapping..."

# 3a. Denovo coordinate mapping
if [[ -f "${TOOLS_R_DIR}/gtf2refbed.R" ]] && [[ -f "${TOOLS_R_DIR}/dn_teprof2_loc.R" ]]; then
    if [[ -f "${SAMPLE_DIR}/denovo/${SAMPLE_NAME}_de_novo_merged_candidates.gtf" ]]; then
        log_info "[09_peptide_screening] Mapping Denovo coordinates..."
        
        Rscript "${TOOLS_R_DIR}/gtf2refbed.R" \
            "${SAMPLE_DIR}/denovo/${SAMPLE_NAME}_de_novo_merged_candidates.gtf" \
            "${GREP_DATA}/6.Denovo/${SAMPLE_NAME}_de_novo_merged_candidates.refbed" \
            > "${LOG_DIR}/peptide_dn_refbed.log" 2>&1
        
        Rscript "${TOOLS_R_DIR}/dn_teprof2_loc.R" \
            "${GREP_RES_SAMPLE}/DN_grep.txt" \
            "${GREP_DATA}/6.Denovo/${SAMPLE_NAME}_de_novo_merged_candidates.refbed" \
            "${GREP_RES_SAMPLE}/DN2loc_grep.txt" \
            FALSE \
            > "${LOG_DIR}/peptide_dn_loc.log" 2>&1
    fi
else
    log_warn "[09_peptide_screening] Denovo mapping scripts not found."
fi

# 3b. TEprof2 coordinate mapping
if [[ -f "${TOOLS_R_DIR}/gtf2refbed.R" ]] && [[ -f "${TOOLS_R_DIR}/dn_teprof2_loc.R" ]] \
    && [[ -f "${TOOLS_R_DIR}/teprof2_refbed.R" ]]; then
    if [[ -f "${SAMPLE_DIR}/teprof2/${SAMPLE_NAME}.stringtie.gtf" ]]; then
        log_info "[09_peptide_screening] Mapping TEprof2 coordinates..."
        
        mkdir -p "${GREP_DATA}/5.TEprof2/gtf" "${GREP_DATA}/5.TEprof2/refbed"
        
        grep "exon" "${SAMPLE_DIR}/teprof2/${SAMPLE_NAME}.stringtie.gtf" \
            > "${GREP_DATA}/5.TEprof2/gtf/${SAMPLE_NAME}_stringtie_exon.gtf"
        
        Rscript "${TOOLS_R_DIR}/gtf2refbed.R" \
            "${GREP_DATA}/5.TEprof2/gtf/${SAMPLE_NAME}_stringtie_exon.gtf" \
            "${GREP_DATA}/5.TEprof2/refbed/${SAMPLE_NAME}_stringtie.refbed" \
            > "${LOG_DIR}/peptide_teprof2_refbed.log" 2>&1
        
        # Generate TEProf2-specific refbed linking ORF headers to transcript IDs
        if [[ -f "${SAMPLE_DIR}/teprof2/candidates.fa" ]]; then
            # Find Step10 CSV (TEProf2 statistics output)
            STEP10_CSV=$(find "${SAMPLE_DIR}/teprof2" -maxdepth 1 -name "*All_TE-derived_Alternative_Isoforms_Statistics.csv" | head -1)
            if [[ -n "$STEP10_CSV" ]]; then
                Rscript "${TOOLS_R_DIR}/teprof2_refbed.R" \
                    "$SAMPLE_NAME" \
                    "$STEP10_CSV" \
                    "${GREP_DATA}/5.TEprof2/refbed/${SAMPLE_NAME}_stringtie.refbed" \
                    "${SAMPLE_DIR}/teprof2/candidates.fa" \
                    "${GREP_DATA}/5.TEprof2/refbed/${SAMPLE_NAME}_TEprof2.refbed" \
                    > "${LOG_DIR}/peptide_teprof2_custom_refbed.log" 2>&1
                
                TEPROF2_REFBED="${GREP_DATA}/5.TEprof2/refbed/${SAMPLE_NAME}_TEprof2.refbed"
            else
                log_warn "[09_peptide_screening] TEProf2 Step10 CSV not found. Using generic refbed."
                TEPROF2_REFBED="${GREP_DATA}/5.TEprof2/refbed/${SAMPLE_NAME}_stringtie.refbed"
            fi
        else
            TEPROF2_REFBED="${GREP_DATA}/5.TEprof2/refbed/${SAMPLE_NAME}_stringtie.refbed"
        fi
        
        Rscript "${TOOLS_R_DIR}/dn_teprof2_loc.R" \
            "${GREP_RES_SAMPLE}/TEprof2_grep.txt" \
            "$TEPROF2_REFBED" \
            "${GREP_RES_SAMPLE}/TEprof22loc_grep.txt" \
            FALSE \
            > "${LOG_DIR}/peptide_teprof2_loc.log" 2>&1
    fi
else
    log_warn "[09_peptide_screening] TEprof2 mapping scripts not found."
fi

# 3c. IR & TElocal coordinate mapping
if [[ -f "${TOOLS_R_DIR}/ir_telocal_loc.R" ]]; then
    # IR coordinate mapping
    if [[ -f "${GREP_RES_SAMPLE}/IR_grep.txt" ]]; then
        log_info "[09_peptide_screening] Mapping IR coordinates..."
        Rscript "${TOOLS_R_DIR}/ir_telocal_loc.R" \
            "${GREP_RES_SAMPLE}/IR_grep.txt" \
            "${GREP_RES_SAMPLE}/IR2loc_grep.txt" \
            FALSE \
            > "${LOG_DIR}/peptide_ir_loc.log" 2>&1
    fi
    
    # TElocal coordinate mapping
    if [[ -f "${GREP_RES_SAMPLE}/TElocal_grep.txt" ]]; then
        log_info "[09_peptide_screening] Mapping TElocal coordinates..."
        Rscript "${TOOLS_R_DIR}/ir_telocal_loc.R" \
            "${GREP_RES_SAMPLE}/TElocal_grep.txt" \
            "${GREP_RES_SAMPLE}/TElocal2loc_grep.txt" \
            FALSE \
            > "${LOG_DIR}/peptide_telocal_loc.log" 2>&1
    fi
else
    log_warn "[09_peptide_screening] IR/TElocal mapping script not found."
fi

# ------------------------------------------------------------------------------
# Step 4: Integrate MS and grep results + FDR calculation
# ------------------------------------------------------------------------------
log_info "[09_peptide_screening] Step 4: Integrating MS/grep results and FDR..."

# 4a. Integrate results
if [[ -f "${TOOLS_R_DIR}/integrate_result.R" ]]; then
    if [[ -f "${GREP_RES_SAMPLE}.xlsx" ]]; then
        Rscript "${TOOLS_R_DIR}/integrate_result.R" \
            "${GREP_RES_SAMPLE}.xlsx" \
            "${GREP_RES}/sp_file.txt" \
            "${GREP_RES}/database_file.txt" \
            > "${LOG_DIR}/peptide_integrate.log" 2>&1
        log_info "[09_peptide_screening] MS/grep integration completed."
    else
        log_warn "[09_peptide_screening] MS result xlsx not found for integration."
    fi
else
    log_warn "[09_peptide_screening] integrate_result.R not found."
fi

# 4b. Generate decoy
if [[ -f "${TOOLS_PY_DIR}/reverse_fasta.py" ]]; then
    log_info "[09_peptide_screening] Generating decoy sequences..."
    python "${TOOLS_PY_DIR}/reverse_fasta.py" \
        -i "$UNIPROT_FASTA" \
        -o "${GREP_DATA}/Decoy/Decoy_uniprotkb.fasta" \
        > "${LOG_DIR}/peptide_decoy.log" 2>&1
else
    log_warn "[09_peptide_screening] reverse_fasta.py not found."
fi

# 4c. Merge decoy FASTA
DECOY_INPUT="${GREP_DATA}/Decoy/Decoy_uniprotkb.fasta"
DECOY_MERGED="${GREP_DATA}/Decoy/Decoy_merge_uniprotkb.fasta"

if [[ -f "$DECOY_INPUT" ]]; then
    if command -v seqkit >/dev/null 2>&1; then
        seqkit seq -w 0 "$DECOY_INPUT" > "$DECOY_MERGED"
    else
        log_warn "[09_peptide_screening] seqkit not found. Decoy may be multi-line."
        cp "$DECOY_INPUT" "$DECOY_MERGED"
    fi
fi

# 4d. Target-decoy FDR calculation
FDR_SCRIPT="${TOOLS_R_DIR}/target_decoy_fdr.R"

if [[ -f "$FDR_SCRIPT" ]]; then
    if [[ -f "${GREP_RES}/Uniport_cut_peptidelist.txt" ]] && [[ -f "$UNIPROT_FASTA" ]] && [[ -f "$DECOY_MERGED" ]]; then
        log_info "[09_peptide_screening] Calculating FDR..."
        Rscript "$FDR_SCRIPT" \
            "${GREP_RES}/Uniport_cut_peptidelist.txt" \
            "$GREP_RES" \
            0.05 \
            "$UNIPROT_FASTA" \
            "$DECOY_MERGED" \
            > "${LOG_DIR}/peptide_fdr.log" 2>&1
        log_info "[09_peptide_screening] FDR calculation completed."
    else
        log_warn "[09_peptide_screening] Missing inputs for FDR calculation."
    fi
else
    log_warn "[09_peptide_screening] target_decoy_fdr.R not found."
fi

# ------------------------------------------------------------------------------
# Step 5: Re-search unmatched peptides
# ------------------------------------------------------------------------------
log_info "[09_peptide_screening] Step 5: Re-searching unmatched peptides..."

if [[ -f "${TOOLS_R_DIR}/peptide_grep.R" ]] && [[ -f "${GREP_RES}/target_decoy_signif_sequences.txt" ]]; then
    nohup Rscript "${TOOLS_R_DIR}/peptide_grep.R" \
        target_decoy \
        "${GREP_RES}/target_decoy_signif_sequences.txt" \
        Uniport_cut \
        "$UNIPROT_FASTA" \
        FALSE \
        "$GREP_RES" \
        > "${LOG_DIR}/peptide_research.log" 2>&1
    log_info "[09_peptide_screening] Step 5 completed."
else
    log_warn "[09_peptide_screening] Skipping Step 5 (missing inputs or script)."
fi

# ------------------------------------------------------------------------------
# Step 6: Final cis/trans classification
# ------------------------------------------------------------------------------
log_info "[09_peptide_screening] Step 6: Final cis/trans classification..."

if [[ -f "${TOOLS_R_DIR}/integrate_cis_trans.R" ]]; then
    nohup Rscript "${TOOLS_R_DIR}/integrate_cis_trans.R" \
        "${GREP_RES}/Uniport_cut_grep.txt" \
        "${GREP_RES}/Uniport_cut_peptide.txt" \
        "${GREP_RES}/target_decoy_result.xls" \
        "${GREP_RES}/0_unmatched_MS_result.xlsx" \
        "${GREP_RES}/sp_file.txt" \
        > "${LOG_DIR}/peptide_cis_trans.log" 2>&1
    log_info "[09_peptide_screening] Step 6 completed."
else
    log_warn "[09_peptide_screening] integrate_cis_trans.R not found. Skipping Step 6."
fi

cd - > /dev/null || true

log_info "[09_peptide_screening] Module completed."
exit 0
