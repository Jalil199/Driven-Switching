#!/usr/bin/env bash

set -u -o pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
stamp="${1:-$(date -u +%Y%m%dT%H%M%SZ)}"
julia_bin="$(command -v julia)"

results_root="${script_dir}/results/detached_${stamp}"
log_path="${script_dir}/run_${stamp}.log"
pid_path="${script_dir}/run_${stamp}.pid"
status_path="${script_dir}/run_${stamp}.status.json"
progress_path="${script_dir}/run_${stamp}.progress.json"
unit_name="fermilink-rhoeq-pade30-${stamp}"

mkdir -p "${results_root}"
printf '%s\n' "${results_root}" > "${script_dir}/latest_results_root.txt"
printf '%s\n' "${log_path}" > "${script_dir}/latest_log_path.txt"
printf '%s\n' "${status_path}" > "${script_dir}/latest_status_path.txt"
printf '%s\n' "${progress_path}" > "${script_dir}/latest_progress_path.txt"
printf '%s\n' "${unit_name}" > "${script_dir}/latest_unit_name.txt"

if command -v systemd-run >/dev/null 2>&1 && systemctl --user is-active default.target >/dev/null 2>&1; then
    systemd-run --user --unit="${unit_name}" --same-dir --collect --quiet \
        "${script_dir}/run_detached_inner_rhoeq_pade30.sh" \
        "${repo_root}" \
        "${script_dir}" \
        "${julia_bin}" \
        "${results_root}" \
        "${log_path}" \
        "${status_path}" \
        "${progress_path}"
    printf '%s\n' "${unit_name}" > "${pid_path}"
    printf '{\"status\":\"launched\",\"launcher\":\"systemd-run\",\"unit_name\":\"%s\",\"pid_path\":\"%s\",\"status_path\":\"%s\",\"results_root\":\"%s\",\"log_path\":\"%s\",\"progress_path\":\"%s\"}\n' \
        "${unit_name}" "${pid_path}" "${status_path}" "${results_root}" "${log_path}" "${progress_path}"
else
    nohup "${script_dir}/run_detached_inner_rhoeq_pade30.sh" \
        "${repo_root}" \
        "${script_dir}" \
        "${julia_bin}" \
        "${results_root}" \
        "${log_path}" \
        "${status_path}" \
        "${progress_path}" \
        </dev/null >/dev/null 2>&1 &

    launcher_pid=$!
    printf '%s\n' "${launcher_pid}" > "${pid_path}"
    printf '{\"status\":\"launched\",\"launcher\":\"nohup\",\"launcher_pid\":%s,\"pid_path\":\"%s\",\"status_path\":\"%s\",\"results_root\":\"%s\",\"log_path\":\"%s\",\"progress_path\":\"%s\"}\n' \
        "${launcher_pid}" "${pid_path}" "${status_path}" "${results_root}" "${log_path}" "${progress_path}"
fi
