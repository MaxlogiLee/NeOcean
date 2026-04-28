#!/bin/bash
# ==============================================================================
# NeOcean Pipeline - Main Controller
# Version: 2.0 (Standardized)
# Usage:
#   bash workflow/main.sh --global-config config/config.yaml \
#                         --sample-config config/samples/MCJ.yaml \
#                         [--modules preprocess,fusion,...] \
#                         [--resume] [--dry-run]
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ------------------------------------------------------------------------------
# Parse Arguments
# ------------------------------------------------------------------------------
GLOBAL_CONFIG=""
SAMPLE_CONFIG=""
MODULES=""
RESUME=false
DRY_RUN=false

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Required:
  -g, --global-config FILE    Global configuration YAML
  -s, --sample-config FILE    Sample configuration YAML

Optional:
  -m, --modules LIST          Comma-separated module list (default: all)
  -r, --resume                Skip modules with existing checkpoints
  -n, --dry-run               Print execution plan without running
  -h, --help                  Show this help

Available modules:
  preprocess, fusion, intron_retention, telocal,
  denovo, teprof2, mutation, peptide_screening
EOF
    exit 0
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -g|--global-config) GLOBAL_CONFIG="$2"; shift 2 ;;
        -s|--sample-config) SAMPLE_CONFIG="$2"; shift 2 ;;
        -m|--modules) MODULES="$2"; shift 2 ;;
        -r|--resume) RESUME=true; shift ;;
        -n|--dry-run) DRY_RUN=true; shift ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

if [[ -z "$GLOBAL_CONFIG" ]] || [[ -z "$SAMPLE_CONFIG" ]]; then
    echo "ERROR: --global-config and --sample-config are required." >&2
    usage
fi

if [[ ! -f "$GLOBAL_CONFIG" ]]; then
    echo "ERROR: Global config not found: $GLOBAL_CONFIG" >&2
    exit 1
fi

if [[ ! -f "$SAMPLE_CONFIG" ]]; then
    echo "ERROR: Sample config not found: $SAMPLE_CONFIG" >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# Load Common Library
# ------------------------------------------------------------------------------
source "${SCRIPT_DIR}/common.sh"
print_banner

# ------------------------------------------------------------------------------
# Resolve Config Values
# ------------------------------------------------------------------------------
log_info "Loading configurations..."
log_info "  Global:  $GLOBAL_CONFIG"
log_info "  Sample:  $SAMPLE_CONFIG"

# Helper to resolve a config key (tries sample config first, then global)
resolve_config() {
    local key="$1"
    local value=""
    value=$(parse_yaml "$SAMPLE_CONFIG" "$key" 2>/dev/null || true)
    if [[ -z "$value" ]]; then
        value=$(parse_yaml "$GLOBAL_CONFIG" "$key" 2>/dev/null || true)
    fi
    echo "$value"
}

# Sample info
SAMPLE_NAME=$(resolve_config "sample.name")
WORK_DIR=$(resolve_config "sample.work_dir")
SEQUENCING_TYPE=$(resolve_config "sample.sequencing_type")

# Global paths
NEOCEAN_ROOT=$(resolve_config "neocean.root_dir")
RESOURCE_DIR=$(resolve_config "neocean.resource_dir")

# Reference paths
GENOME_FASTA=$(resolve_config "reference.genome_fasta")
STAR_INDEX=$(resolve_config "reference.star_index")
REF_GTF=$(resolve_config "reference.gtf")
REF_GTF_TELOCAL=$(resolve_config "reference.gtf_telocal")
TE_INDEX=$(resolve_config "reference.te_index")
TE_LOCATIONS=$(resolve_config "reference.te_locations")
IREAD_BED=$(resolve_config "reference.iread_bed")

# Resource paths
TEPROF2_BIN=$(resolve_config "resources.teprof2.bin_dir")
TEPROF2_DOC=$(resolve_config "resources.teprof2.doc_dir")
TEPROF2_ARGS=$(resolve_config "resources.teprof2.arguments")
SOAPFUSE_PROG=$(resolve_config "resources.soapfuse.program_dir")
SOAPFUSE_CONFIG=$(resolve_config "resources.soapfuse.config")
IREAD_DIR=$(resolve_config "resources.iread.program_dir")
GREP_CODE_DIR=$(resolve_config "resources.grep.code_dir")
BN_FASTA=$(resolve_config "resources.grep.bn_fasta")
UNIPROT_FASTA=$(resolve_config "resources.grep.uniprot_fasta")
NCBI_FASTA=$(resolve_config "resources.grep.ncbi_fasta")

# Compute
THREADS=$(resolve_config "compute.threads")
PARALLEL_JOBS=$(resolve_config "compute.parallel_jobs")
STAR_THREADS=$(resolve_config "compute.star_threads")
STRINGTIE_THREADS=$(resolve_config "compute.stringtie_threads")

# Conda envs
CONDA_TEPROF2=$(resolve_config "conda_envs.teprof2")
CONDA_CUFFLINKS=$(resolve_config "conda_envs.cufflinks")
CONDA_DENOVO=$(resolve_config "conda_envs.denovo")
CONDA_TELOCAL=$(resolve_config "conda_envs.telocal")

# Sample data paths
RNA_FQ1=$(resolve_config "sample.rna.fq1")
RNA_FQ2=$(resolve_config "sample.rna.fq2")
TUMOR_BAM=$(resolve_config "sample.dna.tumor_bam")
NORMAL_BAM=$(resolve_config "sample.dna.normal_bam")
HLA_FILE=$(resolve_config "sample.hla.hla_file")
HLA_ALLELES_RAW=$(parse_yaml "$SAMPLE_CONFIG" "sample.hla.alleles" 2>/dev/null || true)
HLA_ALLELES=$(echo "$HLA_ALLELES_RAW" | paste -sd "," -)
PON_VCF=$(resolve_config "sample.pon.pon_vcf")

# Mutation module resources
MUTATION_SEQ2HLA=$(resolve_config "resources.mutation.seq2hla")
MUTATION_DBSNP=$(resolve_config "resources.mutation.dbsnp")
MUTATION_GNOMAD=$(resolve_config "resources.mutation.gnomad")
MUTATION_VEP_CACHE=$(resolve_config "resources.mutation.vep_cache")
MUTATION_VEP_PLUGINS=$(resolve_config "resources.mutation.vep_plugins")
MUTATION_IEDB_DIR=$(resolve_config "resources.mutation.iedb_dir")
CONDA_SEQ2HLA=$(resolve_config "conda_envs.seq2hla")

# ------------------------------------------------------------------------------
# Validate Essential Configs
# ------------------------------------------------------------------------------
log_info "Validating configuration..."

[[ -z "$SAMPLE_NAME" ]] && { log_error "sample.name is required"; exit 1; }
[[ -z "$WORK_DIR" ]] && { log_error "sample.work_dir is required"; exit 1; }

# Expand any env vars in paths
WORK_DIR=$(eval echo "$WORK_DIR")
GENOME_FASTA=$(eval echo "$GENOME_FASTA")
STAR_INDEX=$(eval echo "$STAR_INDEX")
REF_GTF=$(eval echo "$REF_GTF")

# Check reference files
check_file "$GENOME_FASTA" "Reference genome" || exit 1
check_dir "$STAR_INDEX" "STAR index" || exit 1
check_file "$REF_GTF" "Reference GTF" || exit 1

# Check RNA-seq data
check_file "$RNA_FQ1" "RNA-seq R1" || exit 1
check_file "$RNA_FQ2" "RNA-seq R2" || exit 1

# Check optional DNA data
dna_available=false
if [[ -n "$TUMOR_BAM" ]] && [[ -f "$TUMOR_BAM" ]]; then
    dna_available=true
    log_info "DNA data available for mutation module."
else
    log_warn "DNA data not found. Mutation module will be skipped unless external VCF provided."
fi

# ------------------------------------------------------------------------------
# Initialize Sample Directory Structure
# ------------------------------------------------------------------------------
SAMPLE_DIR="${WORK_DIR}/${SAMPLE_NAME}"
LOG_DIR="${SAMPLE_DIR}/logs"
STATUS_DIR="${SAMPLE_DIR}/status"

log_info "Initializing sample workspace: $SAMPLE_DIR"
init_sample_dirs "$SAMPLE_DIR"

# Global log
MAIN_LOG="${LOG_DIR}/neocean_${SAMPLE_NAME}_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$MAIN_LOG") 2>&1

log_info "=========================================================================="
log_info " Sample:     $SAMPLE_NAME"
log_info " Work Dir:   $SAMPLE_DIR"
log_info " Sequencing: $SEQUENCING_TYPE"
log_info " Threads:    $THREADS"
log_info " Resume:     $RESUME"
log_info " Dry-run:    $DRY_RUN"
log_info "=========================================================================="

# ------------------------------------------------------------------------------
# Export Environment Variables for Sub-scripts
# All modules receive parameters via env vars to minimize argument passing.
# ------------------------------------------------------------------------------
export NEOCEAN_SAMPLE_NAME="$SAMPLE_NAME"
export NEOCEAN_SAMPLE_DIR="$SAMPLE_DIR"
export NEOCEAN_RESOURCE_DIR="$RESOURCE_DIR"
export NEOCEAN_GLOBAL_CONFIG="$GLOBAL_CONFIG"
export NEOCEAN_SAMPLE_CONFIG="$SAMPLE_CONFIG"
export NEOCEAN_SEQUENCING_TYPE="$SEQUENCING_TYPE"

export NEOCEAN_GENOME_FASTA="$GENOME_FASTA"
export NEOCEAN_STAR_INDEX="$STAR_INDEX"
export NEOCEAN_REF_GTF="$REF_GTF"
export NEOCEAN_REF_GTF_TELOCAL="$REF_GTF_TELOCAL"
export NEOCEAN_TE_INDEX="$TE_INDEX"
export NEOCEAN_TE_LOCATIONS="$TE_LOCATIONS"
export NEOCEAN_IREAD_BED="$IREAD_BED"

export NEOCEAN_RNA_FQ1="$RNA_FQ1"
export NEOCEAN_RNA_FQ2="$RNA_FQ2"
export NEOCEAN_TUMOR_BAM="$TUMOR_BAM"
export NEOCEAN_NORMAL_BAM="$NORMAL_BAM"
export NEOCEAN_HLA_FILE="$HLA_FILE"
export NEOCEAN_HLA_ALLELES="$HLA_ALLELES"
export NEOCEAN_PON_VCF="$PON_VCF"

export NEOCEAN_MUTATION_SEQ2HLA="$MUTATION_SEQ2HLA"
export NEOCEAN_MUTATION_DBSNP="$MUTATION_DBSNP"
export NEOCEAN_MUTATION_GNOMAD="$MUTATION_GNOMAD"
export NEOCEAN_MUTATION_VEP_CACHE="$MUTATION_VEP_CACHE"
export NEOCEAN_MUTATION_VEP_PLUGINS="$MUTATION_VEP_PLUGINS"
export NEOCEAN_MUTATION_IEDB_DIR="$MUTATION_IEDB_DIR"
export NEOCEAN_CONDA_SEQ2HLA="$CONDA_SEQ2HLA"

export NEOCEAN_THREADS="$THREADS"
export NEOCEAN_PARALLEL_JOBS="$PARALLEL_JOBS"
export NEOCEAN_STAR_THREADS="$STAR_THREADS"
export NEOCEAN_STRINGTIE_THREADS="$STRINGTIE_THREADS"

export NEOCEAN_TEPROF2_BIN="$TEPROF2_BIN"
export NEOCEAN_TEPROF2_DOC="$TEPROF2_DOC"
export NEOCEAN_TEPROF2_ARGS="$TEPROF2_ARGS"
export NEOCEAN_SOAPFUSE_PROG="$SOAPFUSE_PROG"
export NEOCEAN_SOAPFUSE_CONFIG="$SOAPFUSE_CONFIG"
export NEOCEAN_IREAD_DIR="$IREAD_DIR"
export NEOCEAN_GREP_CODE_DIR="$GREP_CODE_DIR"
export NEOCEAN_BN_FASTA="$BN_FASTA"
export NEOCEAN_UNIPROT_FASTA="$UNIPROT_FASTA"
export NEOCEAN_NCBI_FASTA="$NCBI_FASTA"

export NEOCEAN_CONDA_TEPROF2="$CONDA_TEPROF2"
export NEOCEAN_CONDA_CUFFLINKS="$CONDA_CUFFLINKS"
export NEOCEAN_CONDA_DENOVO="$CONDA_DENOVO"
export NEOCEAN_CONDA_TELOCAL="$CONDA_TELOCAL"

export NEOCEAN_LOG_DIR="$LOG_DIR"
export NEOCEAN_STATUS_DIR="$STATUS_DIR"
export NEOCEAN_DRY_RUN="$DRY_RUN"

# ------------------------------------------------------------------------------
# Determine Module Execution Order
# ------------------------------------------------------------------------------
ALL_MODULES=("preprocess" "fusion" "intron_retention" "telocal" "denovo" "teprof2" "mutation" "peptide_screening")

if [[ -n "$MODULES" ]]; then
    IFS=',' read -ra RUN_MODULES <<< "$MODULES"
else
    RUN_MODULES=("${ALL_MODULES[@]}")
fi

# Module dependencies (downstream requires upstream)
# preprocess -> all others
# fusion, intron_retention, telocal, denovo, teprof2, mutation -> peptide_screening

declare -A MODULE_DEPS=(
    ["peptide_screening"]="preprocess fusion intron_retention telocal denovo teprof2 mutation"
    ["mutation"]="preprocess"
    ["teprof2"]="preprocess"
    ["denovo"]="preprocess"
    ["telocal"]="preprocess"
    ["intron_retention"]="preprocess"
    ["fusion"]="preprocess"
)

# ------------------------------------------------------------------------------
# Pre-flight Checks
# ------------------------------------------------------------------------------
log_info "Pre-flight dependency check..."

# Check module scripts exist
MODULES_DIR="${PROJECT_ROOT}/modules"
for mod in "${RUN_MODULES[@]}"; do
    mod_script="${MODULES_DIR}/$(printf '%02d' $((10#$(grep -n "$mod" <<< "${ALL_MODULES[*]}" | cut -d: -f1))))_${mod}.sh"
    # Actually use fixed naming
    case $mod in
        preprocess) mod_script="${MODULES_DIR}/01_preprocess.sh" ;;
        fusion) mod_script="${MODULES_DIR}/02_fusion.sh" ;;
        intron_retention) mod_script="${MODULES_DIR}/03_intron_retention.sh" ;;
        telocal) mod_script="${MODULES_DIR}/04_telocal.sh" ;;
        denovo) mod_script="${MODULES_DIR}/05_denovo.sh" ;;
        teprof2) mod_script="${MODULES_DIR}/06_teprof2.sh" ;;
        mutation) mod_script="${MODULES_DIR}/07_mutation.sh" ;;
        build_bn) mod_script="${MODULES_DIR}/08_build_bn.sh" ;;
        peptide_screening) mod_script="${MODULES_DIR}/09_peptide_screening.sh" ;;
    esac
    
    if [[ ! -f "$mod_script" ]]; then
        log_error "Module script not found: $mod_script"
        exit 1
    fi
    export "NEOCEAN_MOD_${mod^^}_SCRIPT=$mod_script"
done

log_info "All module scripts found."

# ------------------------------------------------------------------------------
# Dry-run Mode
# ------------------------------------------------------------------------------
if [[ "$DRY_RUN" == "true" ]]; then
    log_info "DRY-RUN MODE: Execution plan for sample '${SAMPLE_NAME}'"
    echo ""
    echo "  Module Execution Order:"
    for mod in "${RUN_MODULES[@]}"; do
        local ckpt="$(get_checkpoint_file "$STATUS_DIR" "$mod")"
        local status="RUN"
        [[ -f "$ckpt" ]] && status="SKIP (checkpoint found)"
        printf "    %-20s %s\n" "$mod" "$status"
    done
    echo ""
    log_info "Dry-run complete. No jobs were executed."
    exit 0
fi

# ------------------------------------------------------------------------------
# Execute Modules
# ------------------------------------------------------------------------------
TOTAL_MODULES=${#RUN_MODULES[@]}
CURRENT=0

for mod in "${RUN_MODULES[@]}"; do
    CURRENT=$((CURRENT + 1))
    log_info "[${CURRENT}/${TOTAL_MODULES}] Processing module: ${mod}"
    
    # Skip if checkpoint exists and resume is true
    if [[ "$RESUME" == "true" ]] && is_module_done "$STATUS_DIR" "$mod"; then
        log_info "  -> Skipping (checkpoint found)"
        continue
    fi
    
    # Check dependencies
    if [[ -n "${MODULE_DEPS[$mod]:-}" ]]; then
        for dep in ${MODULE_DEPS[$mod]}; do
            if [[ ! -f "$(get_checkpoint_file "$STATUS_DIR" "$dep")" ]]; then
                log_error "Dependency not satisfied: ${mod} requires ${dep}"
                exit 1
            fi
        done
    fi
    
    # Special skip logic
    if [[ "$mod" == "mutation" ]] && [[ "$dna_available" != "true" ]]; then
        # Check if external VCF provided
        external_vcf=$(resolve_config "sample.dna.vcf")
        if [[ -z "$external_vcf" ]] || [[ ! -f "$external_vcf" ]]; then
            log_warn "No DNA data or external VCF for mutation module. Skipping."
            touch "$(get_checkpoint_file "$STATUS_DIR" "mutation")"
            continue
        fi
    fi
    
    # Run module
    mod_script_var="NEOCEAN_MOD_${mod^^}_SCRIPT"
    mod_script="${!mod_script_var}"
    
    if ! run_module "$mod" "$mod_script" "$SAMPLE_DIR" "$LOG_DIR" "$STATUS_DIR"; then
        log_error "Pipeline halted due to failure in module: ${mod}"
        exit 1
    fi
done

# ------------------------------------------------------------------------------
# Completion
# ------------------------------------------------------------------------------
log_info "=========================================================================="
log_info " NeOcean Pipeline Completed Successfully"
log_info " Sample:      ${SAMPLE_NAME}"
log_info " Output:      ${SAMPLE_DIR}"
log_info " Main Log:    ${MAIN_LOG}"
log_info " Status Dir:  ${STATUS_DIR}"
log_info "=========================================================================="

exit 0
