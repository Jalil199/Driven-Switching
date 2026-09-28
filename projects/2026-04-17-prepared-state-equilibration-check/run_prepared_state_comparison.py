from __future__ import annotations

import csv
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parent
RESULTS_ROOT = PROJECT_DIR / "results" / "prepared_state_comparison"
SUMMARY_TSV = RESULTS_ROOT / "comparison_summary.tsv"
SUMMARY_JSON = RESULTS_ROOT / "comparison_summary.json"
STATUS_JSON = PROJECT_DIR / "run_status.json"
CASE_DRIVER = PROJECT_DIR.parent / "2026-04-08-tdnegf-spin-baseline" / "run_switching_case.jl"
JULIA = Path("/home/jalil/.local/bin/julia")


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def build_cases() -> list[dict]:
    shared = {
        "initial_configuration": "FM",
        "Hs_update_mode": "nonequilibrium_only",
        "Nx": 2,
        "Ny": 2,
        "dt": 0.1,
        "t_end": 200.0,
        "Bx": 0.0,
        "γso": 0.9,
        "j_sd": 0.3,
        "j_ex": -0.02,
        "j_ani": 0.03,
        "damping": 0.3,
        "results_root": str(RESULTS_ROOT),
    }
    return [
        {
            **shared,
            "run_id": "baseline_nopulse_staticref",
            "analysis_reference_mode": "static_rhoeq",
            "equilibration_time": 0.0,
            "bias_amplitude": 0.0,
            "left_bias_scale": 0.0,
            "right_bias_scale": 0.0,
            "pulse_start": 200.0,
            "pulse_ramp": 1.0,
            "pulse_width": 1.0,
        },
        {
            **shared,
            "run_id": "prepared_nopulse_preparedref",
            "analysis_reference_mode": "prepared_state",
            "equilibration_time": 200.0,
            "bias_amplitude": 0.0,
            "left_bias_scale": 0.0,
            "right_bias_scale": 0.0,
            "pulse_start": 200.0,
            "pulse_ramp": 1.0,
            "pulse_width": 1.0,
        },
        {
            **shared,
            "run_id": "baseline_pulse_staticref",
            "analysis_reference_mode": "static_rhoeq",
            "equilibration_time": 0.0,
            "bias_amplitude": 1.5,
            "left_bias_scale": 0.0,
            "right_bias_scale": 1.0,
            "pulse_start": 100.0,
            "pulse_ramp": 5.0,
            "pulse_width": 100.0,
        },
        {
            **shared,
            "run_id": "prepared_pulse_preparedref",
            "analysis_reference_mode": "prepared_state",
            "equilibration_time": 200.0,
            "bias_amplitude": 1.5,
            "left_bias_scale": 0.0,
            "right_bias_scale": 1.0,
            "pulse_start": 100.0,
            "pulse_ramp": 5.0,
            "pulse_width": 100.0,
        },
    ]


def write_status(status: str, note: str, cases: list[dict]) -> None:
    STATUS_JSON.write_text(
        json.dumps(
            {
                "status": status,
                "timestamp_utc": now_utc(),
                "note": note,
                "results_root": str(RESULTS_ROOT),
                "cases": [
                    {
                        "run_id": case["run_id"],
                        "analysis_reference_mode": case["analysis_reference_mode"],
                        "equilibration_time": case["equilibration_time"],
                    }
                    for case in cases
                ],
            },
            indent=2,
        )
        + "\n"
    )


def run_case(case: dict) -> dict:
    cmd = [str(JULIA), str(CASE_DRIVER)]
    cmd.extend(f"{key}={value}" for key, value in case.items())
    subprocess.run(cmd, check=True, cwd=PROJECT_DIR.parents[1])
    summary_path = RESULTS_ROOT / case["run_id"] / "summary.json"
    return json.loads(summary_path.read_text())


def main() -> None:
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    cases = build_cases()
    write_status("running", "Launching sequential prepared-state comparison cases.", cases)

    summaries = []
    for case in cases:
        summaries.append(run_case(case))

    headers = [
        "run_id",
        "status",
        "analysis_reference_mode",
        "equilibration_time",
        "bias_amplitude",
        "pulse_start",
        "pre_pulse_sigma_delta_absmax",
        "max_sigma_delta_abs",
        "avg_Mz_initial",
        "avg_Mz_final",
        "results_dir",
    ]
    with SUMMARY_TSV.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        for summary in summaries:
            writer.writerow({key: summary.get(key) for key in headers})

    baseline_nopulse = next(s for s in summaries if s["run_id"] == "baseline_nopulse_staticref")
    prepared_nopulse = next(s for s in summaries if s["run_id"] == "prepared_nopulse_preparedref")
    baseline_pulse = next(s for s in summaries if s["run_id"] == "baseline_pulse_staticref")
    prepared_pulse = next(s for s in summaries if s["run_id"] == "prepared_pulse_preparedref")

    comparison = {
        "timestamp_utc": now_utc(),
        "results_root": str(RESULTS_ROOT),
        "baseline_nopulse_pre_pulse_sigma_delta_absmax": baseline_nopulse["pre_pulse_sigma_delta_absmax"],
        "prepared_nopulse_pre_pulse_sigma_delta_absmax": prepared_nopulse["pre_pulse_sigma_delta_absmax"],
        "baseline_pulse_pre_pulse_sigma_delta_absmax": baseline_pulse["pre_pulse_sigma_delta_absmax"],
        "prepared_pulse_pre_pulse_sigma_delta_absmax": prepared_pulse["pre_pulse_sigma_delta_absmax"],
        "nopulse_reduction_factor": baseline_nopulse["pre_pulse_sigma_delta_absmax"] / prepared_nopulse["pre_pulse_sigma_delta_absmax"],
        "pulse_reduction_factor": baseline_pulse["pre_pulse_sigma_delta_absmax"] / prepared_pulse["pre_pulse_sigma_delta_absmax"],
    }
    SUMMARY_JSON.write_text(json.dumps(comparison, indent=2) + "\n")
    write_status("finished", "Prepared-state comparison completed.", cases)


if __name__ == "__main__":
    main()
