from __future__ import annotations

import csv
import json
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parent
RESULTS_ROOT = PROJECT_DIR / "results" / "matched_no_pulse_controls"
STATUS_PATH = PROJECT_DIR / "run_status.json"
MANIFEST_PATH = PROJECT_DIR / "cases_manifest.json"
SCAN_SUMMARY_PATH = RESULTS_ROOT / "scan_summary.tsv"
CASE_MANIFEST_TSV_PATH = RESULTS_ROOT / "case_manifest.tsv"
SCAN_CONFIG_PATH = RESULTS_ROOT / "scan_config.json"
LAUNCHER_LOG = PROJECT_DIR / "parallel_launcher.log"
CASE_DRIVER = PROJECT_DIR.parent / "2026-04-08-tdnegf-spin-baseline" / "run_switching_case.jl"
JULIA = Path("/home/jalil/.local/bin/julia")


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sanitized_token(value: float) -> str:
    return str(value).replace("-", "m").replace(".", "p").replace("+", "")


def mode_tag(mode: str) -> str:
    return "hsfull" if mode == "full_spin_density" else "hsneq"


def control_run_id(spec: dict) -> str:
    return "_".join(
        [
            spec["label"],
            "nopulse",
            f"gso{sanitized_token(spec['gamma_so'])}",
            f"jsd{sanitized_token(spec['j_sd'])}",
            f"jex{sanitized_token(spec['j_ex'])}",
            f"jani{sanitized_token(spec['j_ani'])}",
            f"damp{sanitized_token(spec['damping'])}",
            f"twidth{sanitized_token(spec['pulse_width'])}",
            mode_tag(spec["Hs_update_mode"]),
        ]
    )


def load_cases() -> list[dict]:
    cases = json.loads(MANIFEST_PATH.read_text())
    for spec in cases:
        spec["run_id"] = control_run_id(spec)
        spec["results_dir"] = str(RESULTS_ROOT / spec["run_id"])
    return cases


def write_log(message: str) -> None:
    with LAUNCHER_LOG.open("a") as fh:
        fh.write(f"[{now_utc()}] {message}\n")


def write_case_manifest(cases: list[dict]) -> None:
    headers = [
        "label",
        "source_run_id",
        "run_id",
        "Hs_update_mode",
        "bias_amplitude",
        "gamma_so",
        "j_sd",
        "j_ex",
        "j_ani",
        "damping",
        "pulse_start",
        "pulse_ramp",
        "pulse_width",
        "t_end",
        "results_dir",
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
        "control_type": "goal6_branchb_hsfull_family_neighbor_no_pulse_followup",
        "parallel_execution": True,
        "parallel_mode": "independent_serial_processes",
        "source_project": str(PROJECT_DIR.parent / "2026-04-10-goal6-retained-switch"),
        "timestamp_utc": now_utc(),
    }
    SCAN_CONFIG_PATH.write_text(json.dumps(config, indent=2) + "\n")


def status_snapshot(cases: list[dict], running: dict[str, subprocess.Popen], completed: int, note: str) -> dict:
    case_states = []
    for case in cases:
        summary_path = RESULTS_ROOT / case["run_id"] / "summary.json"
        state = {
            "label": case["label"],
            "source_run_id": case["source_run_id"],
            "run_id": case["run_id"],
            "results_dir": case["results_dir"],
            "log_path": str(PROJECT_DIR / f"{case['run_id']}.log"),
            "status": "pending",
        }
        proc = running.get(case["run_id"])
        if summary_path.exists():
            try:
                summary = json.loads(summary_path.read_text())
                state["status"] = summary.get("status", "done")
                state["avg_Mz_final"] = summary.get("avg_Mz_final")
                state["avg_Mz_min"] = summary.get("avg_Mz_min")
                state["t_final"] = summary.get("t_final")
            except json.JSONDecodeError:
                state["status"] = "summary_unreadable"
        elif proc is not None and proc.poll() is None:
            state["status"] = "running"
            state["pid"] = proc.pid
        elif proc is not None:
            state["status"] = "exited_without_summary"
            state["returncode"] = proc.returncode
        case_states.append(state)

    status = {
        "status": "running" if completed < len(cases) else "finished",
        "completed_cases": completed,
        "total_cases": len(cases),
        "results_root": str(RESULTS_ROOT),
        "timestamp_utc": now_utc(),
        "note": note,
        "case_states": case_states,
    }
    STATUS_PATH.write_text(json.dumps(status, indent=2) + "\n")
    return status


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
        "results_root": str(RESULTS_ROOT),
    }
    cmd = [str(JULIA), str(CASE_DRIVER)]
    cmd.extend(f"{key}={value}" for key, value in args.items())
    return cmd


def read_summary(case: dict) -> dict:
    return json.loads((RESULTS_ROOT / case["run_id"] / "summary.json").read_text())


def write_scan_summary(cases: list[dict]) -> None:
    headers = [
        "label",
        "source_run_id",
        "run_id",
        "status",
        "t_final",
        "avg_Mz_final",
        "avg_Mz_min",
        "max_sigma_delta_abs",
        "pre_pulse_sigma_delta_absmax",
        "results_dir",
    ]
    with SCAN_SUMMARY_PATH.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        for case in cases:
            summary = read_summary(case)
            writer.writerow(
                {
                    "label": case["label"],
                    "source_run_id": case["source_run_id"],
                    "run_id": case["run_id"],
                    "status": summary.get("status", "unknown"),
                    "t_final": summary.get("t_final"),
                    "avg_Mz_final": summary.get("avg_Mz_final"),
                    "avg_Mz_min": summary.get("avg_Mz_min"),
                    "max_sigma_delta_abs": summary.get("max_sigma_delta_abs"),
                    "pre_pulse_sigma_delta_absmax": summary.get("pre_pulse_sigma_delta_absmax"),
                    "results_dir": case["results_dir"],
                }
            )


def main() -> None:
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    cases = load_cases()
    write_case_manifest(cases)
    write_scan_config(cases)
    running: dict[str, subprocess.Popen] = {}

    for case in cases:
        log_path = PROJECT_DIR / f"{case['run_id']}.log"
        with log_path.open("w") as log_fh:
            proc = subprocess.Popen(
                build_command(case),
                cwd=PROJECT_DIR,
                stdout=log_fh,
                stderr=subprocess.STDOUT,
            )
        running[case["run_id"]] = proc
        write_log(f"launched {case['run_id']} pid={proc.pid}")

    completed = 0
    status_snapshot(cases, running, completed, "parallel_batch_launched")

    while completed < len(cases):
        completed = 0
        for run_id, proc in list(running.items()):
            if proc.poll() is not None:
                completed += 1
            else:
                summary_path = RESULTS_ROOT / run_id / "summary.json"
                if summary_path.exists():
                    completed += 1
        status_snapshot(cases, running, completed, "parallel_batch_running")
        if completed >= len(cases):
            break
        time.sleep(60)

    for run_id, proc in running.items():
        returncode = proc.wait()
        write_log(f"completed {run_id} returncode={returncode}")

    write_scan_summary(cases)
    subprocess.run(["python3", str(PROJECT_DIR / "analyze_branchb_no_pulse_followup.py")], check=True, cwd=PROJECT_DIR)
    status_snapshot(cases, running, len(cases), "all_cases_finished_and_analyzed")
    write_log("analysis completed")


if __name__ == "__main__":
    main()
