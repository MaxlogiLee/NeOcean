#!/bin/bash
# ==============================================================================
# NeOcean Module 08: Build Background Normal (BN) Protein Database
# Placeholder module. Should be run once per project, not per sample.
# Output: BN_processed.fa for peptide screening
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
        --project-dir) PROJECT_DIR="$2"; shift 2 ;;
        *) shift ;;
    esac
done

[[ -z "$SAMPLE_DIR" ]] && { log_error "SAMPLE_DIR not set"; exit 1; }

# ------------------------------------------------------------------------------
# Resolve paths
# ------------------------------------------------------------------------------
SAMPLE_NAME="${NEOCEAN_SAMPLE_NAME}"
RESOURCE_DIR="${NEOCEAN_RESOURCE_DIR}"
UNIPROT_FASTA="${NEOCEAN_UNIPROT_FASTA}"
BN_FASTA="${NEOCEAN_BN_FASTA}"

# BN database is typically project-wide, stored in resources
BN_DIR="${RESOURCE_DIR}/Grep/BN"
mkdir -p "$BN_DIR"

# ------------------------------------------------------------------------------
# Pre-flight checks
# ------------------------------------------------------------------------------
log_info "[08_build_bn] Starting BN database build"

# If BN already exists and is non-empty, skip
if [[ -n "$BN_FASTA" ]] && [[ -s "$BN_FASTA" ]]; then
    log_info "[08_build_bn] BN database already exists: ${BN_FASTA}"
    log_info "[08_build_bn] Module completed (skipped)."
    exit 0
fi

# Check if UniProt is available as source
if [[ -n "$UNIPROT_FASTA" ]] && [[ -f "$UNIPROT_FASTA" ]]; then
    log_info "[08_build_bn] UniProt source found: ${UNIPROT_FASTA}"
else
    log_warn "[08_build_bn] UniProt FASTA not found. Cannot build BN database."
    log_warn "[08_build_bn] Please provide UniProt human proteome FASTA."
    exit 0
fi

# ------------------------------------------------------------------------------
# Build BN database (placeholder)
# ------------------------------------------------------------------------------
# TODO: Implement BN database construction logic
# 
# Typical workflow:
# 1. Download/format UniProt human reference proteome
# 2. (Optional) Remove known tumor-specific isoforms
# 3. (Optional) Filter by HLA if project-specific
# 4. Format to single-line FASTA (seqkit seq -w 0)
# 5. Save to BN_processed.fa

BN_OUTPUT="${BN_DIR}/BN_processed.fa"

log_info "[08_build_bn] BN construction placeholder."
log_info "[08_build_bn] Expected steps:"
log_info "[08_build_bn]   1. Format UniProt reference"
log_info "[08_build_bn]   2. Remove tumor-specific entries (optional)"
log_info "[08_build_bn]   3. Standardize FASTA format (single-line sequences)"
log_info "[08_build_bn]   4. Output to ${BN_OUTPUT}"

# For now, if UniProt exists, just reformat it as a placeholder
if command -v seqkit >/dev/null 2>&1; then
    seqkit seq -w 0 "$UNIPROT_FASTA" > "$BN_OUTPUT"
    log_info "[08_build_bn] Created BN database from UniProt (reformatted)."
else
    log_warn "[08_build_bn] seqkit not found. Cannot reformat."
    cp "$UNIPROT_FASTA" "$BN_OUTPUT"
fi

# Update config reference
log_info "[08_build_bn] BN database location: ${BN_OUTPUT}"
log_info "[08_build_bn] Please update global config: resources.grep.bn_fasta"

log_info "[08_build_bn] Module completed (placeholder)."
exit 0
