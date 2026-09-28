from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


PROJECT_DIR = Path(__file__).resolve().parent
RESULTS_ROOT = PROJECT_DIR / "results" / "matched_no_pulse_controls"
GOAL6_ROOT = (
    PROJECT_DIR.parent
    / "2026-04-10-goal6-retained-switch"
    / "results"
    / "queued_goal6_20260410T191626Z_first117"
)
ANALYSIS_DIR = RESULTS_ROOT / "analysis"

CASE_MAP = {
    "goal6_A_main_0003": {
        "pulse_run_id": "goal6_A_main_0003_pa1p5_gso0p9_jsd0p3_jexm0p02_jani0p03_damp0p3_twidth1000p0_hsneq",
        "control_run_id": "goal6_A_main_0003_nopulse_gso0p9_jsd0p3_jexm0p02_jani0p03_damp0p3_twidth1000p0_hsneq",
    },
    "goal6_A_main_0034": {
        "pulse_run_id": "goal6_A_main_0034_pa1p5_gso0p9_jsd0p3_jexm0p02_jani0p01_damp0p7_twidth2000p0_hsneq",
        "control_run_id": "goal6_A_main_0034_nopulse_gso0p9_jsd0p3_jexm0p02_jani0p01_damp0p7_twidth2000p0_hsneq",
    },
}


def read_tsv(path: Path) -> dict[str, np.ndarray]:
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        columns: dict[str, list[float]] = {}
        for row in reader:
            for key, value in row.items():
                columns.setdefault(key, []).append(float(value))
    return {key: np.asarray(values, dtype=float) for key, values in columns.items()}


def sigma_delta_absmax(data: dict[str, np.ndarray]) -> np.ndarray:
    cols = [name for name in data if name.startswith("sigma_delta_")]
    stack = np.vstack([np.abs(data[name]) for name in cols])
    return stack.max(axis=0)


def summarize_case(label: str, pulse_data: dict[str, np.ndarray], control_data: dict[str, np.ndarray]) -> dict[str, float | str]:
    pulse_sigma = sigma_delta_absmax(pulse_data)
    control_sigma = sigma_delta_absmax(control_data)
    excess_sigma = pulse_sigma - control_sigma
    pulse_avgmz = pulse_data["avg_Mz"]
    control_avgmz = control_data["avg_Mz"]
    excess_avgmz = pulse_avgmz - control_avgmz
    pulse_start = 1000.0
    pre_mask = pulse_data["t"] < pulse_start
    post_mask = ~pre_mask

    return {
        "label": label,
        "pulse_run_id": CASE_MAP[label]["pulse_run_id"],
        "control_run_id": CASE_MAP[label]["control_run_id"],
        "pulse_sigma_delta_absmax_peak": float(np.max(pulse_sigma)),
        "control_sigma_delta_absmax_peak": float(np.max(control_sigma)),
        "excess_sigma_delta_absmax_peak": float(np.max(np.abs(excess_sigma))),
        "pulse_sigma_delta_absmax_final": float(pulse_sigma[-1]),
        "control_sigma_delta_absmax_final": float(control_sigma[-1]),
        "excess_sigma_delta_absmax_final": float(excess_sigma[-1]),
        "pre_pulse_excess_sigma_delta_absmax_peak": float(np.max(np.abs(excess_sigma[pre_mask]))),
        "post_pulse_excess_sigma_delta_absmax_peak": float(np.max(np.abs(excess_sigma[post_mask]))),
        "pulse_avg_Mz_final": float(pulse_avgmz[-1]),
        "control_avg_Mz_final": float(control_avgmz[-1]),
        "excess_avg_Mz_final": float(excess_avgmz[-1]),
        "pulse_avg_Mz_min": float(np.min(pulse_avgmz)),
        "control_avg_Mz_min": float(np.min(control_avgmz)),
        "excess_avg_Mz_min": float(np.min(excess_avgmz)),
    }


def write_summary(rows: list[dict[str, float | str]]) -> None:
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    json_path = ANALYSIS_DIR / "pulse_minus_no_pulse_summary.json"
    json_path.write_text(json.dumps(rows, indent=2))

    tsv_path = ANALYSIS_DIR / "pulse_minus_no_pulse_summary.tsv"
    headers = list(rows[0].keys())
    with tsv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def make_plot(series: dict[str, tuple[dict[str, np.ndarray], dict[str, np.ndarray]]]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex="col")
    for col, label in enumerate(["goal6_A_main_0003", "goal6_A_main_0034"]):
        pulse_data, control_data = series[label]
        t = pulse_data["t"]
        pulse_sigma = sigma_delta_absmax(pulse_data)
        control_sigma = sigma_delta_absmax(control_data)
        excess_sigma = pulse_sigma - control_sigma

        axes[0, col].plot(t, pulse_sigma, label="pulse", lw=1.5)
        axes[0, col].plot(t, control_sigma, label="no-pulse", lw=1.2)
        axes[0, col].plot(t, excess_sigma, label="pulse - no-pulse", lw=1.2, ls="--")
        axes[0, col].axvline(1000.0, color="k", lw=0.8, ls=":")
        axes[0, col].set_title(label)
        axes[0, col].set_ylabel("max |sigma_delta|")
        axes[0, col].legend(frameon=False, fontsize=8)

        axes[1, col].plot(t, pulse_data["avg_Mz"], label="pulse", lw=1.5)
        axes[1, col].plot(t, control_data["avg_Mz"], label="no-pulse", lw=1.2)
        axes[1, col].plot(t, pulse_data["avg_Mz"] - control_data["avg_Mz"], label="pulse - no-pulse", lw=1.2, ls="--")
        axes[1, col].axvline(1000.0, color="k", lw=0.8, ls=":")
        axes[1, col].set_xlabel("t")
        axes[1, col].set_ylabel("avg Mz")

    fig.tight_layout()
    fig.savefig(ANALYSIS_DIR / "goal6_pulse_vs_no_pulse_controls.png", dpi=180)


def make_mz_plot(series: dict[str, tuple[dict[str, np.ndarray], dict[str, np.ndarray]]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    for ax, label in zip(axes, ["goal6_A_main_0003", "goal6_A_main_0034"], strict=True):
        pulse_data, control_data = series[label]
        t = pulse_data["t"]

        ax.plot(t, pulse_data["avg_Mz"], label="pulse", lw=1.8)
        ax.plot(t, control_data["avg_Mz"], label="no-pulse", lw=1.5)
        ax.axvline(1000.0, color="k", lw=0.8, ls=":")
        ax.axhline(0.0, color="0.6", lw=0.8, ls="--")
        ax.set_title(label)
        ax.set_xlabel("t")
        ax.set_ylabel("avg Mz")
        ax.set_ylim(-1.05, 1.05)
        ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    fig.savefig(ANALYSIS_DIR / "goal6_avg_Mz_vs_time_pulse_vs_no_pulse.png", dpi=180)


def main() -> None:
    rows = []
    series = {}
    for label, run_ids in CASE_MAP.items():
        pulse_data = read_tsv(GOAL6_ROOT / run_ids["pulse_run_id"] / "observables.tsv")
        control_data = read_tsv(RESULTS_ROOT / run_ids["control_run_id"] / "observables.tsv")
        rows.append(summarize_case(label, pulse_data, control_data))
        series[label] = (pulse_data, control_data)

    write_summary(rows)
    make_plot(series)
    make_mz_plot(series)


if __name__ == "__main__":
    main()
