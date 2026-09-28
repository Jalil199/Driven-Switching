from __future__ import annotations

import csv
import json
import os
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parent
RESULTS_ROOT = PROJECT_DIR / "results" / "brancha_fixed_llg_pade_soc_sweep"
STATUS_PATH = PROJECT_DIR / "run_status.json"
MANIFEST_JSON_PATH = PROJECT_DIR / "cases_manifest.json"
CASE_MANIFEST_TSV_PATH = RESULTS_ROOT / "case_manifest.tsv"
SCAN_CONFIG_PATH = RESULTS_ROOT / "scan_config.json"
SCAN_SUMMARY_PATH = RESULTS_ROOT / "scan_summary.tsv"
LAUNCHER_LOG = PROJECT_DIR / "parallel_launcher.log"
CASE_DRIVER = PROJECT_DIR / "run_brancha_fixed_llg_pade_soc_case.jl"
ANALYZER = PROJECT_DIR / "analyze_brancha_rhoeq_pade_soc.py"
JULIA = Path("/home/jalil/.local/bin/julia")
MAX_PARALLEL = min(4, max(1, (os.cpu_count() or 1) // 2))


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sanitized_token(value: float) -> str:
    return str(value).replace("-", "m").replace(".", "p").replace("+", "")


def build_cases() -> list[dict]:
    shared = {
        "branch_reference_run_id": "goal6_A_main_0003",
        "initial_configuration": "FM",
        "Hs_update_mode": "nonequilibrium_only",
        "Nx": 2,
        "Ny": 2,
        "dt": 0.1,
        "t_end": 1000.0,
        "pulse_start": 1000.0,
        "pulse_ramp": 5.0,
        "pulse_width": 1000.0,
        "Bx": 0.0,
        "left_bias_scale": 0.0,
        "right_bias_scale": 1.0,
        "bias_amplitude": 0.0,
        "j_sd": 0.3,
        "j_ex": -0.02,
        "j_ani": 0.03,
        "damping": 0.3,
        "N_λ1": 49,
    }
    cases = []
    for gamma_so in (0.9, 0.3, 0.1, 0.05):
        for n_lambda2 in (20, 30, 40):
            run_id = f"brancha_fixedllg_gso{sanitized_token(gamma_so)}_pade{n_lambda2}"
            case = dict(shared)
            case.update(
                {
                    "run_id": run_id,
                    "gamma_so": gamma_so,
                    "N_λ2": n_lambda2,
                    "results_dir": str(RESULTS_ROOT / run_id),
                    "log_path": str(PROJECT_DIR / f"{run_id}.log"),
                    "progress_path": str(RESULTS_ROOT / run_id / "progress.json"),
                }
            )
            cases.append(case)
    return cases


def write_log(message: str) -> None:
    with LAUNCHER_LOG.open("a") as fh:
        fh.write(f"[{now_utc()}] {message}\n")


def write_case_manifest(cases: list[dict]) -> None:
    headers = [
        "run_id",
        "branch_reference_run_id",
        "gamma_so",
        "N_λ2",
        "dt",
        "t_end",
        "pulse_start",
        "pulse_width",
        "j_sd",
        "j_ex",
        "j_ani",
        "damping",
        "results_dir",
        "log_path",
        "progress_path",
    ]
    with CASE_MANIFEST_TSV_PATH.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        for case in cases:
            writer.writerow({key: case[key] for key in headers})


def write_scan_config(cases: list[dict]) -> None:
    config = {
        "scan_id": RESULTS_ROOT.name,
        "results_root": str(RESULTS_ROOT),
        "requested_case_count": len(cases),
        "parallel_execution": True,
        "parallel_mode": "independent_serial_processes",
        "max_parallel_cases": MAX_PARALLEL,
        "reference_case": "goal6_A_main_0003",
        "objective": "test pre-pulse sigma_delta decay versus rhoeq Pade pole count and reduced gamma_so in a Branch-A-like fixed-LLG no-pulse setup",
        "timestamp_utc": now_utc(),
    }
    SCAN_CONFIG_PATH.write_text(json.dumps(config, indent=2) + "\n")


def build_command(case: dict) -> list[str]:
    args = {
        "run_id": case["run_id"],
        "initial_configuration": case["initial_configuration"],
        "Hs_update_mode": case["Hs_update_mode"],
        "Nx": case["Nx"],
        "Ny": case["Ny"],
        "dt": case["dt"],
        "t_end": case["t_end"],
        "pulse_start": case["pulse_start"],
        "pulse_ramp": case["pulse_ramp"],
        "pulse_width": case["pulse_width"],
        "Bx": case["Bx"],
        "left_bias_scale": case["left_bias_scale"],
        "right_bias_scale": case["right_bias_scale"],
        "bias_amplitude": case["bias_amplitude"],
        "γso": case["gamma_so"],
        "j_sd": case["j_sd"],
        "j_ex": case["j_ex"],
        "j_ani": case["j_ani"],
        "damping": case["damping"],
        "N_λ1": case["N_λ1"],
        "N_λ2": case["N_λ2"],
        "results_root": str(RESULTS_ROOT),
    }
    cmd = [str(JULIA), str(CASE_DRIVER)]
    cmd.extend(f"{key}={value}" for key, value in args.items())
    return cmd


def load_progress(case: dict) -> dict | None:
    progress_path = Path(case["progress_path"])
    if not progress_path.exists():
        return None
    try:
        return json.loads(progress_path.read_text())
    except json.JSONDecodeError:
        return {"status": "progress_unreadable"}


def load_summary(case: dict) -> dict | None:
    summary_path = RESULTS_ROOT / case["run_id"] / "summary.json"
    if not summary_path.exists():
        return None
    try:
        return json.loads(summary_path.read_text())
    except json.JSONDecodeError:
        return {"status": "summary_unreadable"}


def status_snapshot(cases: list[dict], running: dict[str, subprocess.Popen], completed: int, note: str) -> None:
    case_states = []
    for case in cases:
        state = {
            "run_id": case["run_id"],
            "gamma_so": case["gamma_so"],
            "rhoeq_pade_poles": case["N_λ2"],
            "results_dir": case["results_dir"],
            "log_path": case["log_path"],
            "progress_path": case["progress_path"],
            "status": "pending",
        }
        summary = load_summary(case)
        progress = load_progress(case)
        proc = running.get(case["run_id"])
        if summary is not None:
            state["status"] = summary.get("status", "done")
            state["final_sigma_delta_absmax"] = summary.get("final_sigma_delta_absmax")
            state["late_sigma_delta_absmax_mean"] = summary.get("late_sigma_delta_absmax_mean")
            state["t_final"] = summary.get("t_final")
        elif proc is not None and proc.poll() is None:
            state["status"] = "running"
            state["pid"] = proc.pid
            if progress is not None:
                state["fraction_complete"] = progress.get("fraction_complete")
                state["t_now"] = progress.get("t_now")
                state["wall_elapsed_s"] = progress.get("wall_elapsed_s")
                state["progress_status"] = progress.get("status")
        elif proc is not None:
            state["status"] = "exited_without_summary"
            state["returncode"] = proc.returncode
        case_states.append(state)

    status = {
        "status": "running" if completed < len(cases) else "finished",
        "completed_cases": completed,
        "total_cases": len(cases),
        "max_parallel_cases": MAX_PARALLEL,
        "results_root": str(RESULTS_ROOT),
        "timestamp_utc": now_utc(),
        "note": note,
        "case_states": case_states,
    }
    STATUS_PATH.write_text(json.dumps(status, indent=2) + "\n")


def write_scan_summary(cases: list[dict]) -> None:
    headers = [
        "run_id",
        "gamma_so",
        "rhoeq_pade_poles",
        "status",
        "t_final",
        "max_sigma_delta_abs",
        "final_sigma_delta_absmax",
        "late_sigma_delta_absmax_mean",
        "late10_sigma_delta_absmax_slope",
        "final_to_peak_sigma_delta_ratio",
        "results_dir",
    ]
    with SCAN_SUMMARY_PATH.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        for case in cases:
            summary = load_summary(case)
            if summary is None:
                continue
            writer.writerow(
                {
                    "run_id": case["run_id"],
                    "gamma_so": case["gamma_so"],
                    "rhoeq_pade_poles": case["N_λ2"],
                    "status": summary.get("status", "unknown"),
                    "t_final": summary.get("t_final"),
                    "max_sigma_delta_abs": summary.get("max_sigma_delta_abs"),
                    "final_sigma_delta_absmax": summary.get("final_sigma_delta_absmax"),
                    "late_sigma_delta_absmax_mean": summary.get("late_sigma_delta_absmax_mean"),
                    "late10_sigma_delta_absmax_slope": summary.get("late10_sigma_delta_absmax_slope"),
                    "final_to_peak_sigma_delta_ratio": summary.get("final_to_peak_sigma_delta_ratio"),
                    "results_dir": case["results_dir"],
                }
            )


def main() -> None:
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    cases = build_cases()
    MANIFEST_JSON_PATH.write_text(json.dumps(cases, indent=2) + "\n")
    write_case_manifest(cases)
    write_scan_config(cases)
    running: dict[str, subprocess.Popen] = {}
    pending = cases.copy()
    completed = 0
    status_snapshot(cases, running, completed, "launcher_initialized")

    while completed < len(cases):
        while pending and sum(proc.poll() is None for proc in running.values()) < MAX_PARALLEL:
            case = pending.pop(0)
            Path(case["results_dir"]).mkdir(parents=True, exist_ok=True)
            log_path = Path(case["log_path"])
            with log_path.open("w") as log_fh:
                proc = subprocess.Popen(
                    build_command(case),
                    cwd=PROJECT_DIR,
                    stdout=log_fh,
                    stderr=subprocess.STDOUT,
                    env={**os.environ, "JULIA_NUM_THREADS": "1"},
                )
            running[case["run_id"]] = proc
            write_log(f"launched {case['run_id']} pid={proc.pid}")

        completed = 0
        for case in cases:
            summary = load_summary(case)
            proc = running.get(case["run_id"])
            if summary is not None:
                completed += 1
            elif proc is not None and proc.poll() is not None:
                completed += 1
        status_snapshot(cases, running, completed, "parallel_batch_running")
        if completed >= len(cases):
            break
        time.sleep(30)

    for run_id, proc in running.items():
        returncode = proc.wait()
        write_log(f"completed {run_id} returncode={returncode}")

    write_scan_summary(cases)
    subprocess.run(["python3", str(ANALYZER)], check=True, cwd=PROJECT_DIR)
    status_snapshot(cases, running, len(cases), "all_cases_finished_and_analyzed")
    write_log("analysis completed")


if __name__ == "__main__":
    main()
