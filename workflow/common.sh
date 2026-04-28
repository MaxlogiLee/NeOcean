#!/bin/bash
# ==============================================================================
# NeOcean Pipeline - Common Function Library
# Version: 2.0 (Standardized)
# ==============================================================================
set -euo pipefail

# ------------------------------------------------------------------------------
# Logging
# ------------------------------------------------------------------------------
NEOCEAN_LOG_LEVEL="${NEOCEAN_LOG_LEVEL:-INFO}"

_log() {
    local level="$1"
    local msg="$2"
    local ts
    ts=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${ts}] [${level}] ${msg}"
}

log_debug() { [[ "${NEOCEAN_LOG_LEVEL}" == "DEBUG" ]] && _log "DEBUG" "$1" || true; }
log_info()  { _log "INFO" "$1"; }
log_warn()  { _log "WARN" "$1" >&2; }
log_error() { _log "ERROR" "$1" >&2; }

# ------------------------------------------------------------------------------
# YAML Parsing (requires python3 + PyYAML)
# ------------------------------------------------------------------------------
parse_yaml() {
    local yaml_file="$1"
    local key="$2"
    python3 -c "
import yaml, sys
try:
    with open('$yaml_file', 'r') as f:
        data = yaml.safe_load(f)
    keys = '$key'.split('.')
    for k in keys:
        if data is None:
            break
        data = data.get(k, None)
    if data is None:
        sys.exit(1)
    if isinstance(data, bool):
        print('true' if data else 'false')
    elif isinstance(data, list):
        print('\n'.join(str(x) for x in data))
    else:
        print(data)
except Exception as e:
    sys.exit(1)
" 2>/dev/null
}

yaml_exists() {
    local yaml_file="$1"
    local key="$2"
    parse_yaml "$yaml_file" "$key" >/dev/null 2>&1
}

# ------------------------------------------------------------------------------
# File/Directory Checks
# ------------------------------------------------------------------------------
check_file() {
    local file="$1"
    local label="${2:-File}"
    if [[ ! -f "$file" ]]; then
        log_error "${label} not found: $file"
        return 1
    fi
    log_debug "${label} OK: $file"
}

check_dir() {
    local dir="$1"
    local label="${2:-Directory}"
    if [[ ! -d "$dir" ]]; then
        log_error "${label} not found: $dir"
        return 1
    fi
    log_debug "${label} OK: $dir"
}

check_files() {
    local failed=0
    for f in "$@"; do
        if [[ ! -f "$f" ]]; then
            log_error "Missing file: $f"
            failed=1
        fi
    done
    return $failed
}

# ------------------------------------------------------------------------------
# Tool Checks
# ------------------------------------------------------------------------------
check_tool() {
    local tool="$1"
    if ! command -v "$tool" >/dev/null 2>&1; then
        log_error "Required tool not found in PATH: $tool"
        return 1
    fi
    log_debug "Tool OK: $tool"
}

check_tools() {
    local failed=0
    for t in "$@"; do
        check_tool "$t" || failed=1
    done
    return $failed
}

# ------------------------------------------------------------------------------
# Conda Wrapper (preserves original environment names)
# ------------------------------------------------------------------------------
conda_run() {
    local env_name="$1"
    shift
    if command -v conda >/dev/null 2>&1; then
        conda run -n "$env_name" --no-capture-output "$@"
    else
        log_error "conda not found. Cannot activate environment: $env_name"
        return 1
    fi
}

# ------------------------------------------------------------------------------
# Checkpoint System (Resume support)
# ------------------------------------------------------------------------------
get_checkpoint_file() {
    local status_dir="$1"
    local module_name="$2"
    echo "${status_dir}/${module_name}.done"
}

is_module_done() {
    local status_dir="$1"
    local module_name="$2"
    local ckpt
    ckpt=$(get_checkpoint_file "$status_dir" "$module_name")
    [[ -f "$ckpt" ]]
}

mark_module_done() {
    local status_dir="$1"
    local module_name="$2"
    local ckpt
    ckpt=$(get_checkpoint_file "$status_dir" "$module_name")
    mkdir -p "$status_dir"
    touch "$ckpt"
    log_info "Checkpoint created: ${module_name}.done"
}

clear_checkpoint() {
    local status_dir="$1"
    local module_name="$2"
    local ckpt
    ckpt=$(get_checkpoint_file "$status_dir" "$module_name")
    [[ -f "$ckpt" ]] && rm -f "$ckpt" && log_info "Checkpoint cleared: ${module_name}.done"
}

# ------------------------------------------------------------------------------
# Module Runner
# ------------------------------------------------------------------------------
run_module() {
    local module_name="$1"
    local module_script="$2"
    local sample_dir="$3"
    local log_dir="$4"
    local status_dir="$5"
    shift 5
    # remaining args passed to module script

    local ckpt
    ckpt=$(get_checkpoint_file "$status_dir" "$module_name")

    if [[ -f "$ckpt" ]]; then
        log_info "Module '${module_name}' already completed (checkpoint found). Skipping."
        return 0
    fi

    log_info "============================================================================"
    log_info "Starting module: ${module_name}"
    log_info "Sample directory: ${sample_dir}"
    log_info "Log file: ${log_dir}/${module_name}.log"
    log_info "============================================================================"

    mkdir -p "$log_dir" "$status_dir"

    if [[ ! -f "$module_script" ]]; then
        log_error "Module script not found: $module_script"
        return 1
    fi

    local start_ts end_ts elapsed
    start_ts=$(date +%s)

    set +e
    bash "$module_script" \
        --sample-dir "$sample_dir" \
        --log-dir "$log_dir" \
        --status-dir "$status_dir" \
        "$@" \
        > "${log_dir}/${module_name}.log" 2>&1
    local exit_code=$?
    set -e

    end_ts=$(date +%s)
    elapsed=$((end_ts - start_ts))

    if [[ $exit_code -eq 0 ]]; then
        mark_module_done "$status_dir" "$module_name"
        log_info "Module '${module_name}' completed successfully (${elapsed}s)."
        return 0
    else
        log_error "Module '${module_name}' failed with exit code ${exit_code} (${elapsed}s)."
        log_error "Check log: ${log_dir}/${module_name}.log"
        return 1
    fi
}

# ------------------------------------------------------------------------------
# Directory Setup for a Sample
# ------------------------------------------------------------------------------
init_sample_dirs() {
    local sample_dir="$1"
    mkdir -p "${sample_dir}/logs"
    mkdir -p "${sample_dir}/status"
    mkdir -p "${sample_dir}/star"
    mkdir -p "${sample_dir}/fusion"
    mkdir -p "${sample_dir}/intron_retention"
    mkdir -p "${sample_dir}/telocal"
    mkdir -p "${sample_dir}/denovo"
    mkdir -p "${sample_dir}/teprof2"
    mkdir -p "${sample_dir}/mutation"
    mkdir -p "${sample_dir}/grep/data"
    mkdir -p "${sample_dir}/grep/res"
}

# ------------------------------------------------------------------------------
# Print pipeline banner
# ------------------------------------------------------------------------------
print_banner() {
    cat <<'EOF'
==========================================================================
  NeOcean: Tumor Neoantigen Discovery Pipeline v2.0
==========================================================================
EOF
}
