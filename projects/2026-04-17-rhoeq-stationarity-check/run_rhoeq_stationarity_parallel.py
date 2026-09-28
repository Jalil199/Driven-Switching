from __future__ import annotations

import csv
import json
import os
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parent
RESULTS_ROOT = PROJECT_DIR / "results" / "stationarity_cases"
STATUS_PATH = PROJECT_DIR / "run_status.json"
MANIFEST_JSON_PATH = PROJECT_DIR / "cases_manifest.json"
CASE_MANIFEST_TSV_PATH = RESULTS_ROOT / "case_manifest.tsv"
SCAN_CONFIG_PATH = RESULTS_ROOT / "scan_config.json"
SCAN_SUMMARY_PATH = RESULTS_ROOT / "scan_summary.tsv"
SUMMARY_MD = PROJECT_DIR / "summary.md"
LAUNCHER_LOG = PROJECT_DIR / "parallel_launcher.log"
CASE_DRIVER = PROJECT_DIR / "run_rhoeq_stationarity_case.jl"
JULIA = Path("/home/jalil/.local/bin/julia")
MAX_PARALLEL = min(3, max(1, (os.cpu_count() or 1) // 4))


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
        "t_end": 20.0,
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
    cases: list[dict] = []
    for gamma_so in (0.9, 0.1, 0.05):
        for n_lambda2 in (20, 40):
            run_id = f"rhoeq_stationarity_gso{sanitized_token(gamma_so)}_pade{n_lambda2}"
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
        "objective": "test whether rho_eq(H_ab) with zero auxiliaries is stationary for the block TDNEGF RHS at fixed Branch-A-like H_ab and zero bias",
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
    path = Path(case["progress_path"])
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError:
        return {"status": "progress_unreadable"}


def load_summary(case: dict) -> dict | None:
    path = RESULTS_ROOT / case["run_id"] / "summary.json"
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
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
            state["rhs_total_absmax"] = summary.get("rhs_total_absmax")
            state["rhs_psi_absmax"] = summary.get("rhs_psi_absmax")
            state["rhoeq_guess_sigma_delta_absmax_final"] = summary.get("rhoeq_guess_sigma_delta_absmax_final")
            state["zero_init_sigma_delta_absmax_final"] = summary.get("zero_init_sigma_delta_absmax_final")
        elif proc is not None and proc.poll() is None:
            state["status"] = "running"
            state["pid"] = proc.pid
            if progress is not None:
                state["phase"] = progress.get("phase")
                state["fraction_complete"] = progress.get("fraction_complete")
                state["t_now"] = progress.get("t_now")
                state["wall_elapsed_s"] = progress.get("wall_elapsed_s")
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


def write_scan_summary(cases: list[dict]) -> list[dict]:
    headers = [
        "run_id",
        "gamma_so",
        "rhoeq_pade_poles",
        "status",
        "rhs_rho_absmax",
        "rhs_psi_absmax",
        "rhs_omega11_absmax",
        "rhs_omega12_absmax",
        "rhs_omega21_absmax",
        "rhs_total_absmax",
        "rhoeq_guess_sigma_delta_absmax_peak",
        "rhoeq_guess_sigma_delta_absmax_final",
        "rhoeq_guess_rho_minus_rhoeq_absmax_final",
        "zero_init_sigma_delta_absmax_final",
        "zero_init_rho_minus_rhoeq_absmax_final",
        "results_dir",
    ]
    rows: list[dict] = []
    with SCAN_SUMMARY_PATH.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        for case in cases:
            summary = load_summary(case)
            if summary is None:
                continue
            row = {
                "run_id": case["run_id"],
                "gamma_so": case["gamma_so"],
                "rhoeq_pade_poles": case["N_λ2"],
                "status": summary.get("status", "unknown"),
                "rhs_rho_absmax": summary.get("rhs_rho_absmax"),
                "rhs_psi_absmax": summary.get("rhs_psi_absmax"),
                "rhs_omega11_absmax": summary.get("rhs_omega11_absmax"),
                "rhs_omega12_absmax": summary.get("rhs_omega12_absmax"),
                "rhs_omega21_absmax": summary.get("rhs_omega21_absmax"),
                "rhs_total_absmax": summary.get("rhs_total_absmax"),
                "rhoeq_guess_sigma_delta_absmax_peak": summary.get("rhoeq_guess_sigma_delta_absmax_peak"),
                "rhoeq_guess_sigma_delta_absmax_final": summary.get("rhoeq_guess_sigma_delta_absmax_final"),
                "rhoeq_guess_rho_minus_rhoeq_absmax_final": summary.get("rhoeq_guess_rho_minus_rhoeq_absmax_final"),
                "zero_init_sigma_delta_absmax_final": summary.get("zero_init_sigma_delta_absmax_final"),
                "zero_init_rho_minus_rhoeq_absmax_final": summary.get("zero_init_rho_minus_rhoeq_absmax_final"),
                "results_dir": case["results_dir"],
            }
            writer.writerow(row)
            rows.append(row)
    return rows


def write_summary_md(rows: list[dict]) -> None:
    rows_ok = [row for row in rows if row["status"] == "done"]
    status = "done" if len(rows_ok) == len(rows) and rows else "blocked"
    with SUMMARY_MD.open("w") as fh:
        fh.write("- prompt: intentemos a ver si funciona\n")
        fh.write("- procedures: se armó un diagnóstico directo sobre el caso de referencia `goal6_A_main_0003` con espines fijos y bias cero. Para cada combinación `gamma_so ∈ {0.9, 0.1, 0.05}` y `N_λ2 ∈ {20, 40}` se construyó `u_guess = [rho_eq(H_ab), 0, 0]`, se evaluó el RHS bloqueado `eom_tdnegf_blocks!` en `t=0`, y luego se propagaron en paralelo dos relajaciones cortas (`t_end=200`): una desde `u_guess` y otra desde el arranque estándar `u0=0`.\n")
        fh.write(f"- generated data: {SCAN_SUMMARY_PATH}\n")
        if rows_ok:
            strongest_rhs = max(rows_ok, key=lambda row: float(row["rhs_total_absmax"]))
            smallest_rhs = min(rows_ok, key=lambda row: float(row["rhs_total_absmax"]))
            strongest_drift = max(rows_ok, key=lambda row: float(row["rhoeq_guess_sigma_delta_absmax_peak"]))
            fh.write(
                "- generated data: "
                f"menor `rhs_total_absmax` en `{smallest_rhs['run_id']}` = {float(smallest_rhs['rhs_total_absmax']):.6e}; "
                f"mayor `rhs_total_absmax` en `{strongest_rhs['run_id']}` = {float(strongest_rhs['rhs_total_absmax']):.6e}; "
                f"mayor deriva desde `u_guess` en `{strongest_drift['run_id']}` con `rhoeq_guess_sigma_delta_absmax_peak={float(strongest_drift['rhoeq_guess_sigma_delta_absmax_peak']):.6e}`.\n"
            )
            fh.write("- debugging clue: si `u_guess=[rho_eq,0,0]` fuera un estacionario consistente del sistema extendido, deberían salir `rhs_total_absmax≈0` y `rhoeq_guess_sigma_delta_absmax_peak≈0`; estos tests permiten separar esa hipótesis de la mera relajación desde `u0=0`.\n")
        fh.write(f"- status: {status}\n")


def main() -> None:
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    cases = build_cases()
    MANIFEST_JSON_PATH.write_text(json.dumps(cases, indent=2) + "\n")
    write_case_manifest(cases)
    write_scan_config(cases)

    running: dict[str, subprocess.Popen] = {}
    pending = cases.copy()
    completed = 0
    write_log(f"Launching {len(cases)} stationarity cases with max_parallel={MAX_PARALLEL}")

    while completed < len(cases):
        while pending and len(running) < MAX_PARALLEL:
            case = pending.pop(0)
            log_path = Path(case["log_path"])
            log_path.parent.mkdir(parents=True, exist_ok=True)
            env = os.environ.copy()
            env["JULIA_NUM_THREADS"] = "1"
            env["FERMILINK_PROGRESS_EVERY_STEPS"] = "250"
            with log_path.open("w") as log_fh:
                proc = subprocess.Popen(
                    build_command(case),
                    cwd=PROJECT_DIR.parents[1],
                    stdout=log_fh,
                    stderr=subprocess.STDOUT,
                    env=env,
                )
            running[case["run_id"]] = proc
            write_log(f"started {case['run_id']} pid={proc.pid}")

        finished_ids: list[str] = []
        for run_id, proc in running.items():
            if proc.poll() is None:
                continue
            finished_ids.append(run_id)
            completed += 1
            write_log(f"finished {run_id} returncode={proc.returncode}")
        for run_id in finished_ids:
            running.pop(run_id, None)

        status_snapshot(cases, running, completed, "stationarity batch active")
        if completed < len(cases):
            time.sleep(5.0)

    rows = write_scan_summary(cases)
    write_summary_md(rows)
    status_snapshot(cases, running, completed, "stationarity batch finished")
    write_log("all stationarity cases completed")


if __name__ == "__main__":
    main()
