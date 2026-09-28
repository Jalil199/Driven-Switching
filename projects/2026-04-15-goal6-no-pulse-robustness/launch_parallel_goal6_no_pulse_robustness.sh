#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
timestamp="$(date -u +%Y%m%dT%H%M%SZ)"
log_path="${script_dir}/launch_${timestamp}.log"
pid_path="${script_dir}/launch_${timestamp}.pid"

nohup python3 "${script_dir}/run_goal6_no_pulse_robustness_parallel.py" >"${log_path}" 2>&1 &
runner_pid=$!
printf '%s\n' "${runner_pid}" >"${pid_path}"

cat <<EOF
launched=true
timestamp_utc=${timestamp}
pid=${runner_pid}
log_path=${log_path}
pid_path=${pid_path}
status_path=${script_dir}/run_status.json
EOF
