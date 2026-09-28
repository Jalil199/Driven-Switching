from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parent
REPO_ROOT = PROJECT_DIR.parent.parent
PHASE2_SUMMARY = REPO_ROOT / "projects/2026-04-15-goal6-phase2-clean-family-mechanism/results/phase2_mechanism_clean_families/analysis/phase2_mechanism_summary.tsv"
GOAL6_ROOT = REPO_ROOT / "projects/2026-04-10-goal6-retained-switch/results/queued_goal6_20260410T191626Z_first117"
OUT_TABLE = PROJECT_DIR / "brancha_sigma_delta_summary.tsv"
OUT_FIG = PROJECT_DIR / "brancha_sigma_delta_overview.png"


def load_branch_a_cases() -> pd.DataFrame:
    df = pd.read_csv(PHASE2_SUMMARY, sep="\t")
    return df[df["label"].str.startswith("goal6_A_")].copy()


def site_component_columns(prefix: str, axis: str) -> list[str]:
    return [f"{prefix}_{axis}_site{i}" for i in range(1, 5)]


def build_case_metrics(row: pd.Series) -> tuple[pd.DataFrame, dict[str, float | str]]:
    run_id = str(row["pulse_run_id"])
    obs_path = GOAL6_ROOT / run_id / "observables.tsv"
    obs = pd.read_csv(obs_path, sep="\t")
    pulse_start = float(row["pulse_start"])
    pulse_end = pulse_start + float(row["pulse_width"]) + 5.0

    avg = pd.DataFrame({"t": obs["t"].to_numpy(float)})
    for axis in "xyz":
        avg[f"sigma_delta_{axis}_avg"] = obs[site_component_columns("sigma_delta", axis)].mean(axis=1)
    avg["sigma_delta_absmax"] = obs[
        [col for axis in "xyz" for col in site_component_columns("sigma_delta", axis)]
    ].abs().max(axis=1)
    avg["avg_Mz"] = obs["avg_Mz"].to_numpy(float)

    pre_mask = avg["t"] < pulse_start
    pulse_mask = (avg["t"] >= pulse_start) & (avg["t"] <= pulse_end)
    post_mask = avg["t"] > pulse_end

    metrics: dict[str, float | str] = {
        "label": str(row["label"]),
        "family": str(row["family"]),
        "pulse_run_id": run_id,
        "pulse_start": pulse_start,
        "pulse_end": pulse_end,
        "avg_Mz_final": float(row["pulse_avg_Mz_final"]),
        "sigma_delta_absmax_peak": float(avg["sigma_delta_absmax"].max()),
        "sigma_delta_absmax_peak_time": float(avg.loc[avg["sigma_delta_absmax"].idxmax(), "t"]),
        "sigma_delta_absmax_final": float(avg["sigma_delta_absmax"].iloc[-1]),
    }
    for axis in "xyz":
        series = avg[f"sigma_delta_{axis}_avg"]
        metrics[f"sigma_delta_{axis}_avg_pre_mean"] = float(series[pre_mask].mean())
        metrics[f"sigma_delta_{axis}_avg_pre_absmax"] = float(series[pre_mask].abs().max())
        metrics[f"sigma_delta_{axis}_avg_pulse_absmax"] = float(series[pulse_mask].abs().max())
        metrics[f"sigma_delta_{axis}_avg_post_absmax"] = float(series[post_mask].abs().max())
        metrics[f"sigma_delta_{axis}_avg_final"] = float(series.iloc[-1])
    return avg, metrics


def make_figure(case_traces: list[tuple[pd.Series, pd.DataFrame]]) -> None:
    ncases = len(case_traces)
    fig, axes = plt.subplots(ncases, 1, figsize=(11, 2.6 * ncases), sharex=True, constrained_layout=True)
    if ncases == 1:
        axes = [axes]

    colors = {"x": "#d62728", "y": "#1f77b4", "z": "#2ca02c", "abs": "#222222"}

    for ax, (row, trace) in zip(axes, case_traces):
        pulse_start = float(row["pulse_start"])
        pulse_end = pulse_start + float(row["pulse_width"]) + 5.0
        ax.axvspan(pulse_start, pulse_end, color="#f2d388", alpha=0.35, lw=0)
        for axis in "xyz":
            ax.plot(trace["t"], trace[f"sigma_delta_{axis}_avg"], lw=1.6, color=colors[axis], label=rf"$\langle \sigma_{{\Delta,{axis}}}\rangle$")
        ax.plot(trace["t"], trace["sigma_delta_absmax"], lw=1.1, color=colors["abs"], ls="--", label=r"$|\sigma_\Delta|_{\max}$")
        ax.axhline(0.0, color="#666666", lw=0.8)
        ax.set_ylabel(r"$\sigma-\sigma_{\mathrm{eq}}$")
        ax.set_title(
            f"{row['label']} | width={float(row['pulse_width']):.0f} | final avg_Mz={float(row['pulse_avg_Mz_final']):.3f}",
            fontsize=10,
        )
        ax.grid(alpha=0.2, lw=0.4)

    axes[0].legend(loc="upper right", ncol=4, fontsize=8, frameon=False)
    axes[-1].set_xlabel("t")
    fig.suptitle("Branch A successful runs: site-averaged $\\sigma(t)-\\sigma_{eq}(t)$ and per-time absolute maximum", fontsize=13)
    fig.savefig(OUT_FIG, dpi=220)


def main() -> None:
    cases = load_branch_a_cases()
    case_traces: list[tuple[pd.Series, pd.DataFrame]] = []
    metrics = []
    for _, row in cases.iterrows():
        trace, case_metrics = build_case_metrics(row)
        case_traces.append((row, trace))
        metrics.append(case_metrics)

    pd.DataFrame(metrics).to_csv(OUT_TABLE, sep="\t", index=False)
    make_figure(case_traces)


if __name__ == "__main__":
    main()
