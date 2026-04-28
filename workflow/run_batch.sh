#!/bin/bash
# ==============================================================================
# NeOcean Batch Runner
# Execute NeOcean for multiple samples sequentially or in parallel.
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Required:
  -g, --global-config FILE    Global configuration YAML
  -l, --sample-list FILE      Sample list (one sample config path per line)
  -d, --sample-dir DIR        Directory containing sample YAML files

Optional:
  -m, --modules LIST          Comma-separated module list
  -r, --resume                Resume mode
  -n, --dry-run               Dry-run mode
  -p, --parallel N            Run N samples in parallel (default: 1)
  -h, --help                  Show this help

Examples:
  # Run all samples in a directory
  bash workflow/run_batch.sh -g config/config.yaml -d config/samples/

  # Run specific samples from a list
  bash workflow/run_batch.sh -g config/config.yaml -l sample_list.txt -m preprocess,peptide_screening

  # Dry-run all samples
  bash workflow/run_batch.sh -g config/config.yaml -d config/samples/ -n
EOF
    exit 0
}

GLOBAL_CONFIG=""
SAMPLE_LIST=""
SAMPLE_DIR=""
MODULES=""
RESUME=""
DRY_RUN=""
PARALLEL=1

while [[ $# -gt 0 ]]; do
    case $1 in
        -g|--global-config) GLOBAL_CONFIG="$2"; shift 2 ;;
        -l|--sample-list) SAMPLE_LIST="$2"; shift 2 ;;
        -d|--sample-dir) SAMPLE_DIR="$2"; shift 2 ;;
        -m|--modules) MODULES="$2"; shift 2 ;;
        -r|--resume) RESUME="--resume"; shift ;;
        -n|--dry-run) DRY_RUN="--dry-run"; shift ;;
        -p|--parallel) PARALLEL="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

if [[ -z "$GLOBAL_CONFIG" ]]; then
    echo "ERROR: --global-config is required." >&2
    exit 1
fi

# Build sample list
SAMPLES=()
if [[ -n "$SAMPLE_LIST" ]] && [[ -f "$SAMPLE_LIST" ]]; then
    while IFS= read -r line; do
        [[ -n "$line" ]] && SAMPLES+=("$line")
    done < "$SAMPLE_LIST"
elif [[ -n "$SAMPLE_DIR" ]] && [[ -d "$SAMPLE_DIR" ]]; then
    for f in "$SAMPLE_DIR"/*.yaml; do
        [[ -f "$f" ]] && SAMPLES+=("$f")
    done
else
    echo "ERROR: Either --sample-list or --sample-dir is required." >&2
    exit 1
fi

TOTAL=${#SAMPLES[@]}
echo "=========================================================================="
echo " NeOcean Batch Runner"
echo " Total samples: ${TOTAL}"
echo " Parallel: ${PARALLEL}"
echo " Global config: ${GLOBAL_CONFIG}"
echo "=========================================================================="

# Build main.sh options
OPTS=""
[[ -n "$MODULES" ]] && OPTS="${OPTS} --modules ${MODULES}"
[[ -n "$RESUME" ]] && OPTS="${OPTS} ${RESUME}"
[[ -n "$DRY_RUN" ]] && OPTS="${OPTS} ${DRY_RUN}"

# Run samples
CURRENT=0
for sample_config in "${SAMPLES[@]}"; do
    CURRENT=$((CURRENT + 1))
    echo ""
    echo "[${CURRENT}/${TOTAL}] Processing: $(basename "$sample_config")"
    
    if [[ "$PARALLEL" -gt 1 ]]; then
        bash "${SCRIPT_DIR}/main.sh" \
            --global-config "$GLOBAL_CONFIG" \
            --sample-config "$sample_config" \
            $OPTS &
        
        # Control parallelism
        if (( CURRENT % PARALLEL == 0 )) && (( CURRENT < TOTAL )); then
            echo "Waiting for batch to complete..."
            wait
        fi
    else
        bash "${SCRIPT_DIR}/main.sh" \
            --global-config "$GLOBAL_CONFIG" \
            --sample-config "$sample_config" \
            $OPTS
    fi
done

# Wait for remaining background jobs
if [[ "$PARALLEL" -gt 1 ]]; then
    wait
fi

echo ""
echo "=========================================================================="
echo " Batch processing completed."
echo "=========================================================================="
