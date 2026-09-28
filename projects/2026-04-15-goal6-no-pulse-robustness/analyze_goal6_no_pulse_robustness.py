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
MANIFEST_PATH = PROJECT_DIR / "cases_manifest.json"
ANALYSIS_DIR = RESULTS_ROOT / "analysis"


def read_tsv(path: Path) -> dict[str, np.ndarray]:
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols: dict[str, list[float]] = {}
        for row in reader:
            for key, value in row.items():
                cols.setdefault(key, []).append(float(value))
    return {key: np.asarray(values, dtype=float) for key, values in cols.items()}


def sigma_delta_absmax(data: dict[str, np.ndarray]) -> np.ndarray:
    cols = [name for name in data if name.startswith("sigma_delta_")]
    return np.vstack([np.abs(data[name]) for name in cols]).max(axis=0)


def control_run_id(spec: dict) -> str:
    def tok(value: float) -> str:
        return str(value).replace("-", "m").replace(".", "p").replace("+", "")

    mode = "hsfull" if spec["Hs_update_mode"] == "full_spin_density" else "hsneq"
    return "_".join(
        [
            spec["label"],
            "nopulse",
            f"gso{tok(spec['gamma_so'])}",
            f"jsd{tok(spec['j_sd'])}",
            f"jex{tok(spec['j_ex'])}",
            f"jani{tok(spec['j_ani'])}",
            f"damp{tok(spec['damping'])}",
            f"twidth{tok(spec['pulse_width'])}",
            mode,
        ]
    )


def first_cross_zero_time(values: np.ndarray, t: np.ndarray) -> float | None:
    mask = values <= 0.0
    if not np.any(mask):
        return None
    return float(t[int(np.argmax(mask))])


def summarize_case(spec: dict, pulse_data: dict[str, np.ndarray], control_data: dict[str, np.ndarray]) -> dict[str, float | str | None]:
    pulse_sigma = sigma_delta_absmax(pulse_data)
    control_sigma = sigma_delta_absmax(control_data)
    excess_sigma = pulse_sigma - control_sigma
    pulse_avgmz = pulse_data["avg_Mz"]
    control_avgmz = control_data["avg_Mz"]
    pulse_start = float(spec["pulse_start"])
    pre_mask = pulse_data["t"] < pulse_start
    post_mask = ~pre_mask

    return {
        "label": spec["label"],
        "pulse_run_id": spec["source_run_id"],
        "control_run_id": control_run_id(spec),
        "pulse_avg_Mz_final": float(pulse_avgmz[-1]),
        "control_avg_Mz_final": float(control_avgmz[-1]),
        "pulse_avg_Mz_min": float(np.min(pulse_avgmz)),
        "control_avg_Mz_min": float(np.min(control_avgmz)),
        "control_first_cross_zero_time": first_cross_zero_time(control_avgmz, control_data["t"]),
        "pulse_first_cross_zero_time": first_cross_zero_time(pulse_avgmz, pulse_data["t"]),
        "pulse_sigma_delta_absmax_peak": float(np.max(pulse_sigma)),
        "control_sigma_delta_absmax_peak": float(np.max(control_sigma)),
        "pre_pulse_excess_sigma_delta_absmax_peak": float(np.max(np.abs(excess_sigma[pre_mask]))),
        "post_pulse_excess_sigma_delta_absmax_peak": float(np.max(np.abs(excess_sigma[post_mask]))),
        "excess_sigma_delta_absmax_final": float(excess_sigma[-1]),
    }


def write_outputs(rows: list[dict]) -> None:
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    (ANALYSIS_DIR / "pulse_minus_no_pulse_summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    headers = list(rows[0].keys())
    with (ANALYSIS_DIR / "pulse_minus_no_pulse_summary.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def make_mz_plot(specs: list[dict], series: dict[str, tuple[dict[str, np.ndarray], dict[str, np.ndarray]]]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=False, sharey=True)
    axes_flat = axes.flatten()
    for ax, spec in zip(axes_flat, specs, strict=True):
        pulse_data, control_data = series[spec["label"]]
        ax.plot(pulse_data["t"], pulse_data["avg_Mz"], lw=1.8, label="pulse")
        ax.plot(control_data["t"], control_data["avg_Mz"], lw=1.5, label="no-pulse")
        ax.axvline(float(spec["pulse_start"]), color="k", lw=0.8, ls=":")
        ax.axhline(0.0, color="0.6", lw=0.8, ls="--")
        ax.set_title(spec["label"])
        ax.set_xlabel("t")
        ax.set_ylabel("avg Mz")
        ax.set_ylim(-1.05, 1.05)
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(ANALYSIS_DIR / "goal6_avg_Mz_vs_time_pulse_vs_no_pulse_robustness.png", dpi=180)


def make_summary_note(rows: list[dict]) -> None:
    lines = [
        "# Goal 6 No-Pulse Robustness Expansion",
        "",
        "## Control outcomes",
    ]
    for row in rows:
        crossing = "none" if row["control_first_cross_zero_time"] is None else f"{row['control_first_cross_zero_time']:.1f}"
        lines.append(
            f"- {row['label']}: control `avg_Mz_final={row['control_avg_Mz_final']:.6f}`, `avg_Mz_min={row['control_avg_Mz_min']:.6f}`, `first_cross_zero_time={crossing}`, "
            f"pulse/control raw `max|sigma_delta|={row['pulse_sigma_delta_absmax_peak']:.6f}/{row['control_sigma_delta_absmax_peak']:.6f}`, "
            f"post-pulse excess sigma peak `={row['post_pulse_excess_sigma_delta_absmax_peak']:.6f}`."
        )
    (ANALYSIS_DIR / "robustness_review.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    specs = json.loads(MANIFEST_PATH.read_text())
    rows = []
    series = {}
    for spec in specs:
        control_id = control_run_id(spec)
        pulse_data = read_tsv(GOAL6_ROOT / spec["source_run_id"] / "observables.tsv")
        control_data = read_tsv(RESULTS_ROOT / control_id / "observables.tsv")
        rows.append(summarize_case(spec, pulse_data, control_data))
        series[spec["label"]] = (pulse_data, control_data)
    write_outputs(rows)
    make_mz_plot(specs, series)
    make_summary_note(rows)


if __name__ == "__main__":
    main()
