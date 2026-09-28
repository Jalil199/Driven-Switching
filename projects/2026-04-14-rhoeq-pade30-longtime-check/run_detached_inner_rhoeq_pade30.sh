#!/usr/bin/env bash

set -u -o pipefail

repo_root="$1"
script_dir="$2"
julia_bin="$3"
results_root="$4"
log_path="$5"
status_path="$6"
progress_path="$7"

cd "${repo_root}" || exit 1
export JULIA_NUM_THREADS=1
export FERMILINK_RESULTS_ROOT="${results_root}"
export FERMILINK_PROGRESS_PATH="${progress_path}"
export FERMILINK_PROGRESS_EVERY_STEPS="${FERMILINK_PROGRESS_EVERY_STEPS:-200}"

printf '{"status":"running","launched_at_utc":"%s","results_root":"%s","log_path":"%s","progress_path":"%s"}\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${results_root}" "${log_path}" "${progress_path}" > "${status_path}"

set +e
"${julia_bin}" "${script_dir}/run_rhoeq_pade30_longtime_check.jl" >> "${log_path}" 2>&1
exit_code=$?

printf '{"status":"finished","finished_at_utc":"%s","exit_code":%s,"results_root":"%s","log_path":"%s","progress_path":"%s"}\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${exit_code}" "${results_root}" "${log_path}" "${progress_path}" > "${status_path}"
exit "${exit_code}"
