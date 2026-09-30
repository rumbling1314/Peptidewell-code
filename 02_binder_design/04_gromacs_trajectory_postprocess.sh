#!/bin/bash
# -----------------------------------------------------------------------------
# GROMACS 分子动力学轨迹后处理
#
# 需在包含 md.tpr 和 md.xtc 的目录中运行，或通过 TARGET_DIR 指定目录。
# 处理步骤:
#   1. trjconv -pbc nojump
#   2. trjconv -center -pbc mol -ur compact
#   3. trjconv -fit rot+trans
#   4. 导出末帧/指定时刻 md_fit.gro
# -----------------------------------------------------------------------------

set -eo pipefail

# 可通过环境变量覆盖，例如:
#   TARGET_DIR=/path/to/MD/001/output bash 04_gromacs_trajectory_postprocess.sh  # 001 为示例编号，具体以实际项目为准
WORK_DIR="${TARGET_DIR:-$(pwd)}"
cd "$WORK_DIR" || { echo "错误: 无法进入 $WORK_DIR"; exit 1; }

DRY_RUN="${DRY_RUN:-0}"
FORCE="${FORCE:-0}"
ALLOW_LOGIN_NODE="${ALLOW_LOGIN_NODE:-0}"
DUMP_TIME_PS="${DUMP_TIME_PS:-50000}"
INDEX_FILE="$WORK_DIR/index.ndx"

has_tpr="NO"
has_xtc="NO"
nojump_done="NA"
center_done="NA"
fit_done="NA"
fit_gro_done="NA"
status="FAILED"
GMX_SOURCE="UNSET"
declare -a GMX_CMD=()

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$WORK_DIR/postprocess.log" >&2
}

append_command_output() {
    "$@" >>"$WORK_DIR/postprocess.log" 2>&1
}

emit_summary() {
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$WORK_DIR" "$has_tpr" "$has_xtc" "$nojump_done" "$center_done" "$fit_done" "$fit_gro_done" "$status"
}

finish_with_status() {
    status="$1"
    log "Final status: $status"
    emit_summary
    exit 0
}

trap 'rc=$?; if [[ $rc != 0 ]]; then log "Unexpected shell error with exit code $rc"; status="FAILED_SHELL"; emit_summary; fi' EXIT

resolve_gmx() {
    if [[ -n "${GMX:-}" ]]; then
        if command -v "$GMX" >/dev/null 2>&1; then
            GMX_CMD=("$GMX")
            GMX_SOURCE="env:GMX"
            return 0
        fi
        if [[ -x "$GMX" ]]; then
            GMX_CMD=("$GMX")
            GMX_SOURCE="env:GMX"
            return 0
        fi
    fi

    if [[ -f /opt/gromacsgpu/bin/GMXRC ]]; then
        set +u
        # shellcheck disable=SC1091
        source /opt/gromacsgpu/bin/GMXRC
        set -u
        if command -v gmx >/dev/null 2>&1; then
            GMX_CMD=("$(command -v gmx)")
            GMX_SOURCE="/opt/gromacsgpu/bin/GMXRC"
            return 0
        fi
    fi

    if [[ -x /opt/gromacs_amd/bin/gmx ]]; then
        GMX_CMD=(/opt/gromacs_amd/bin/gmx)
        GMX_SOURCE="/opt/gromacs_amd/bin/gmx"
        return 0
    fi

    if [[ -f /opt/gromacs_amd/bin/GMXRC ]]; then
        set +u
        # shellcheck disable=SC1091
        source /opt/gromacs_amd/bin/GMXRC
        set -u
        if command -v gmx >/dev/null 2>&1; then
            GMX_CMD=("$(command -v gmx)")
            GMX_SOURCE="/opt/gromacs_amd/bin/GMXRC"
            return 0
        fi
    fi

    if command -v mamba >/dev/null 2>&1; then
        for env_name in ${GMX_ENV_CANDIDATES:-gromacs base}; do
            if mamba run -n "$env_name" gmx --version >>"$WORK_DIR/postprocess.log" 2>&1; then
                GMX_CMD=(mamba run -n "$env_name" gmx)
                GMX_SOURCE="mamba:$env_name"
                return 0
            fi
        done
    fi

    if command -v conda >/dev/null 2>&1; then
        for env_name in ${GMX_ENV_CANDIDATES:-gromacs base}; do
            if conda run -n "$env_name" gmx --version >>"$WORK_DIR/postprocess.log" 2>&1; then
                GMX_CMD=(conda run -n "$env_name" gmx)
                GMX_SOURCE="conda:$env_name"
                return 0
            fi
        done
    fi

    if command -v gmx >/dev/null 2>&1; then
        GMX_CMD=("$(command -v gmx)")
        GMX_SOURCE="PATH"
        return 0
    fi

    return 1
}

declare -A GROUP_IDS=()

capture_groups() {
    local idx=0
    local line=""
    local group_name=""
    while IFS= read -r line; do
        if [[ "$line" =~ ^\[[[:space:]]*(.+)[[:space:]]*\]$ ]]; then
            group_name="${BASH_REMATCH[1]}"
            group_name="$(printf '%s' "$group_name" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"
            GROUP_IDS["$group_name"]="$idx"
            ((idx+=1))
        fi
    done < "$INDEX_FILE"
}

require_group() {
    local group_name="$1"
    if [[ -n "${GROUP_IDS[$group_name]:-}" ]]; then
        printf '%s' "${GROUP_IDS[$group_name]}"
        return 0
    fi
    return 1
}

choose_fit_group() {
    local candidate=""
    for candidate in Backbone Protein C-alpha; do
        if [[ -n "${GROUP_IDS[$candidate]:-}" ]]; then
            printf '%s' "${GROUP_IDS[$candidate]}"
            return 0
        fi
    done
    return 1
}

run_gmx_with_input() {
    local input_text="$1"
    shift
    log "Running: ${GMX_CMD[*]} $*"
    printf '%b' "$input_text" | "${GMX_CMD[@]}" "$@" >>"$WORK_DIR/postprocess.log" 2>&1
}

if [[ -f md.tpr ]]; then
    has_tpr="YES"
fi
if [[ -f md.xtc ]]; then
    has_xtc="YES"
fi

log "System directory: $WORK_DIR"
log "DRY_RUN=$DRY_RUN FORCE=$FORCE ALLOW_LOGIN_NODE=$ALLOW_LOGIN_NODE"

if [[ "$has_tpr" != "YES" || "$has_xtc" != "YES" ]]; then
    log "Missing required inputs: md.tpr or md.xtc"
    finish_with_status "MISSING_FILES"
fi

if [[ -f md_fit.xtc && "$FORCE" != "1" ]]; then
    nojump_done="YES"
    center_done="YES"
    fit_done="YES"
    if [[ -f md_fit.gro ]]; then
        fit_gro_done="YES"
    else
        fit_gro_done="NO"
    fi
    log "Existing md_fit.xtc detected and FORCE!=1; skipping"
    finish_with_status "SKIPPED_ALREADY_DONE"
fi

if [[ "$DRY_RUN" == "1" ]]; then
    log "DRY_RUN=1; no GROMACS commands will be executed"
    finish_with_status "DRY_RUN"
fi

if [[ -z "${SLURM_JOB_ID:-}" && "$ALLOW_LOGIN_NODE" != "1" ]]; then
    log "Refusing to run trajectory processing outside a Slurm allocation"
    finish_with_status "REFUSED_LOGIN_NODE"
fi

if ! resolve_gmx; then
    log "Unable to locate a usable GROMACS command"
    finish_with_status "FAILED_ENV"
fi

log "Using GROMACS from: $GMX_SOURCE"
if ! append_command_output "${GMX_CMD[@]}" --version; then
    log "gmx --version failed for source $GMX_SOURCE"
    finish_with_status "FAILED_ENV"
fi

if ! run_gmx_with_input 'q\n' make_ndx -f md.tpr -o "$INDEX_FILE"; then
    log "make_ndx failed"
    finish_with_status "FAILED_MAKE_NDX"
fi

capture_groups

system_group="$(require_group System || true)"
protein_group="$(require_group Protein || true)"
fit_group="$(choose_fit_group || true)"

if [[ -z "$system_group" || -z "$protein_group" || -z "$fit_group" ]]; then
    log "Required index groups not found. System=$system_group Protein=$protein_group Fit=$fit_group"
    finish_with_status "FAILED_GROUPS"
fi

log "Resolved groups: System=$system_group Protein=$protein_group Fit=$fit_group"

if [[ ! -f md_nojump.xtc || "$FORCE" == "1" ]]; then
    if ! run_gmx_with_input "$system_group\n" trjconv -s md.tpr -f md.xtc -n "$INDEX_FILE" -o md_nojump.xtc -pbc nojump; then
        log "Nojump step failed"
        finish_with_status "FAILED_NOJUMP"
    fi
else
    log "Skipping nojump; output already exists"
fi
nojump_done="YES"

if [[ ! -f md_center.xtc || "$FORCE" == "1" ]]; then
    if ! run_gmx_with_input "$protein_group\n$system_group\n" trjconv -s md.tpr -f md_nojump.xtc -n "$INDEX_FILE" -o md_center.xtc -center -pbc mol -ur compact; then
        log "Center step failed"
        finish_with_status "FAILED_CENTER"
    fi
else
    log "Skipping center; output already exists"
fi
center_done="YES"

if [[ ! -f md_fit.xtc || "$FORCE" == "1" ]]; then
    if ! run_gmx_with_input "$fit_group\n$system_group\n" trjconv -s md.tpr -f md_center.xtc -n "$INDEX_FILE" -o md_fit.xtc -fit rot+trans; then
        log "Fit step failed"
        finish_with_status "FAILED_FIT"
    fi
else
    log "Skipping fit; output already exists"
fi
fit_done="YES"

if [[ ! -f md_fit.gro || "$FORCE" == "1" ]]; then
    if ! run_gmx_with_input "$system_group\n" trjconv -s md.tpr -f md_fit.xtc -n "$INDEX_FILE" -o md_fit.gro -dump "$DUMP_TIME_PS"; then
        log "Dump step failed"
        finish_with_status "FAILED_DUMP"
    fi
else
    log "Skipping fit.gro export; output already exists"
fi
fit_gro_done="YES"

log "Trajectory post-processing completed successfully"
finish_with_status "SUCCESS"
