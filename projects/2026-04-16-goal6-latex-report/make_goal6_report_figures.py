from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parent
FIG_DIR = PROJECT_DIR / "figures"

GOAL5_ANALYSIS_DIR = (
    PROJECT_DIR.parent
    / "2026-04-08-tdnegf-spin-baseline"
    / "results"
    / "queued_parameter_scan_20260410T001920Z_goal5"
    / "analysis"
)
GOAL5_SELECTED = GOAL5_ANALYSIS_DIR / "goal5_selected_pulse_window_candidates.tsv"
GOAL5_CONTAMINATED_RUN = (
    PROJECT_DIR.parent
    / "2026-04-08-tdnegf-spin-baseline"
    / "results"
    / "queued_parameter_scan_20260410T001920Z_goal5"
    / "scan_0111_pa2p0_gso0p7_jsd0p3_jexm0p01_bx0p0_jani0p0_damp0p3_tstart1000p0_tramp5p0_twidth1500p0_lbs0p0_rbs1p0_hsfull"
)
GOAL5_SUMMARY = GOAL5_ANALYSIS_DIR / "summary.md"

GOAL6_PROJECT = PROJECT_DIR.parent / "2026-04-10-goal6-retained-switch"
GOAL6_ANALYSIS_DIR = (
    GOAL6_PROJECT / "results" / "queued_goal6_20260410T191626Z_first117" / "analysis"
)
GOAL6_IN_PULSE = GOAL6_ANALYSIS_DIR / "goal6_in_pulse_ranked_runs.tsv"
GOAL6_RETENTION = GOAL6_ANALYSIS_DIR / "goal6_retention_ranked_runs.tsv"
GOAL6_BRANCH_SUMMARY = GOAL6_ANALYSIS_DIR / "goal6_branch_summary.tsv"

PHASE2_PROJECT = PROJECT_DIR.parent / "2026-04-15-goal6-phase2-clean-family-mechanism"
PHASE2_MANIFEST = PHASE2_PROJECT / "cases_manifest.json"
PHASE2_SUMMARY = (
    PHASE2_PROJECT
    / "results"
    / "phase2_mechanism_clean_families"
    / "analysis"
    / "phase2_mechanism_summary.tsv"
)

BRANCHB_PROJECT = PROJECT_DIR.parent / "2026-04-15-goal6-branchb-no-pulse-followup"
BRANCHB_SCAN_SUMMARY = (
    BRANCHB_PROJECT
    / "results"
    / "matched_no_pulse_controls"
    / "scan_summary.tsv"
)

FIXED_LLG_T1000 = (
    PROJECT_DIR.parent
    / "2026-04-13-fixed-llg-no-pulse-relaxation"
    / "results"
    / "fixed_llg_no_pulse_20260413T201501"
    / "fixedllg_gso0p05"
    / "observables.tsv"
)
RHOEQ_T2000 = (
    PROJECT_DIR.parent
    / "2026-04-13-rhoeq-longtime-relaxation"
    / "results"
    / "rhoeq_longtime_20260413T205653"
    / "pade20_gso0p05_t2000"
    / "observables.tsv"
)

SIGMA_DELTA_RE = re.compile(r"^sigma_delta_([xyz])_site(\d+)$")
COMPONENTS = ("x", "y", "z")
COMPONENT_COLORS = {"x": "#2563eb", "y": "#f59e0b", "z": "#059669"}
MZ_COLORS = {
    "clean": "#0f766e",
    "delayed": "#c2410c",
    "contaminated": "#b91c1c",
    "branch_b": "#7c3aed",
    "pulse": "#b91c1c",
    "control": "#475569",
}
PULSE_FILL = "#dbeafe"
PULSE_EDGE = "#1d4ed8"
FIG_BG = "#fcfcfb"
GRID_COLOR = "#d4d4d4"


plt.rcParams.update(
    {
        "figure.facecolor": FIG_BG,
        "axes.facecolor": "#ffffff",
        "axes.edgecolor": "#94a3b8",
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.titlesize": 11,
        "axes.labelsize": 11,
        "xtick.color": "#334155",
        "ytick.color": "#334155",
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "grid.color": GRID_COLOR,
        "grid.alpha": 0.35,
        "grid.linewidth": 0.6,
        "legend.frameon": False,
        "savefig.facecolor": FIG_BG,
        "savefig.bbox": "tight",
    }
)


def ensure_dir() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)


def read_observables(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def label_from_run_id(run_id: str) -> str:
    return run_id.replace("_pa", "\npa", 1).replace("_hs", "\nhs", 1)


def sigma_component_absmax(df: pd.DataFrame) -> dict[str, np.ndarray]:
    traces: dict[str, np.ndarray] = {}
    for component in COMPONENTS:
        cols = [col for col in df.columns if SIGMA_DELTA_RE.match(col) and SIGMA_DELTA_RE.match(col).group(1) == component]
        values = df[cols].to_numpy(dtype=float)
        traces[component] = np.max(np.abs(values), axis=1)
    return traces


def sigma_delta_global_absmax(df: pd.DataFrame) -> np.ndarray:
    cols = [col for col in df.columns if SIGMA_DELTA_RE.match(col)]
    values = df[cols].to_numpy(dtype=float)
    return np.max(np.abs(values), axis=1)


def add_pulse_window(ax: plt.Axes, pulse_start: float, pulse_width: float) -> None:
    pulse_end = pulse_start + pulse_width
    ax.axvspan(pulse_start, pulse_end, color=PULSE_FILL, alpha=0.8, lw=0, zorder=0)
    ax.axvline(pulse_start, color=PULSE_EDGE, lw=1.1, ls=(0, (3, 2)), alpha=0.9)
    ax.axvline(pulse_end, color=PULSE_EDGE, lw=1.0, ls=(0, (3, 2)), alpha=0.8)
    ymax = ax.get_ylim()[1]
    ax.text(
        0.5 * (pulse_start + pulse_end),
        ymax - 0.08 * (ax.get_ylim()[1] - ax.get_ylim()[0]),
        "pulse on",
        ha="center",
        va="top",
        fontsize=8.5,
        color=PULSE_EDGE,
    )


def style_axis(ax: plt.Axes, *, ylabel: str, xlabel: str = "time", ylim: tuple[float, float] | None = None) -> None:
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True)


def add_case_box(ax: plt.Axes, text: str) -> None:
    ax.text(
        0.015,
        0.04,
        text,
        transform=ax.transAxes,
        fontsize=8.5,
        va="bottom",
        ha="left",
        bbox={"boxstyle": "round,pad=0.28", "facecolor": "#ffffff", "alpha": 0.9, "edgecolor": "#cbd5e1"},
    )


def write_goal5_representative_dynamics() -> None:
    selected = pd.read_csv(GOAL5_SELECTED, sep="\t")
    clean_in_pulse = selected.iloc[0]
    delayed_retention = selected.iloc[1]
    contaminated = {
        "run_id": GOAL5_CONTAMINATED_RUN.name,
        "results_dir": str(GOAL5_CONTAMINATED_RUN),
        "pulse_start": 1000.0,
        "pulse_width": 1500.0,
        "title": "Goal 5: strongest final negative state,\nbut already negative before the pulse",
        "color": "tab:red",
    }
    cases = [
        {
            "run_id": clean_in_pulse["run_id"],
            "results_dir": clean_in_pulse["results_dir"],
            "pulse_start": float(clean_in_pulse["pulse_start"]),
            "pulse_width": float(clean_in_pulse["pulse_width"]),
            "title": "Clean pulse-triggered inversion without retention",
            "color": MZ_COLORS["clean"],
        },
        {
            "run_id": delayed_retention["run_id"],
            "results_dir": delayed_retention["results_dir"],
            "pulse_start": float(delayed_retention["pulse_start"]),
            "pulse_width": float(delayed_retention["pulse_width"]),
            "title": "Moderate in-pulse response followed by delayed fall",
            "color": MZ_COLORS["delayed"],
        },
        contaminated,
    ]
    contaminated["title"] = "Most negative final state, but contaminated before pulse"
    contaminated["color"] = MZ_COLORS["contaminated"]

    fig, axes = plt.subplots(3, 1, figsize=(10.5, 10.5), sharex=False, sharey=True)
    for ax, case in zip(axes, cases):
        df = read_observables(Path(case["results_dir"]) / "observables.tsv")
        t = df["t"].to_numpy(dtype=float)
        mz = df["avg_Mz"].to_numpy(dtype=float)
        ax.set_ylim(-1.1, 1.05)
        add_pulse_window(ax, case["pulse_start"], case["pulse_width"])
        ax.plot(t, mz, lw=2.4, color=case["color"], solid_capstyle="round")
        ax.axhline(0.0, color="#334155", lw=1.0, ls="--", alpha=0.8)
        ax.axhline(-0.8, color="#94a3b8", lw=0.9, ls=(0, (4, 3)))
        ax.set_title(case["title"], fontsize=11)
        style_axis(ax, ylabel=r"$\langle M_z \rangle$", xlabel="time", ylim=(-1.1, 1.05))
        add_case_box(ax, label_from_run_id(case["run_id"]))
    fig.suptitle("Three distinct time-domain outcomes that shaped the switching campaign", fontsize=15, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig(FIG_DIR / "goal5_representative_dynamics.png", dpi=220)
    plt.close(fig)


def write_goal6_representative_dynamics() -> None:
    in_pulse = pd.read_csv(GOAL6_IN_PULSE, sep="\t")
    retention = pd.read_csv(GOAL6_RETENTION, sep="\t")

    selected_ids = [
        "goal6_A_main_0003_pa1p5_gso0p9_jsd0p3_jexm0p02_jani0p03_damp0p3_twidth1000p0_hsneq",
        "goal6_A_main_0034_pa1p5_gso0p9_jsd0p3_jexm0p02_jani0p01_damp0p7_twidth2000p0_hsneq",
        "goal6_B_main_0042_pa2p0_gso0p9_jsd0p3_jexm0p01_jani0p0_damp0p5_twidth3000p0_hsfull",
    ]

    rows = []
    for run_id in selected_ids:
        row = in_pulse.loc[in_pulse["run_id"] == run_id]
        if row.empty:
            row = retention.loc[retention["run_id"] == run_id]
        rows.append(row.iloc[0])

    titles = [
        "Branch A: strong in-pulse switching with retention",
        "Branch A: delayed but clean retained switching",
        "Branch B: promising pulsed trace, unresolved baseline",
    ]
    colors = [MZ_COLORS["clean"], "#2563eb", MZ_COLORS["branch_b"]]

    fig, axes = plt.subplots(3, 1, figsize=(10.5, 10.5), sharex=False, sharey=True)
    for ax, row, title, color in zip(axes, rows, titles, colors):
        df = read_observables(Path(row["results_dir"]) / "observables.tsv")
        t = df["t"].to_numpy(dtype=float)
        mz = df["avg_Mz"].to_numpy(dtype=float)
        ax.set_ylim(-1.1, 1.05)
        add_pulse_window(ax, 1000.0, float(row["pulse_width"]))
        ax.plot(t, mz, lw=2.4, color=color, solid_capstyle="round")
        ax.axhline(0.0, color="#334155", lw=1.0, ls="--", alpha=0.8)
        ax.axhline(-0.8, color="#94a3b8", lw=0.9, ls=(0, (4, 3)))
        ax.set_title(title, fontsize=11)
        style_axis(ax, ylabel=r"$\langle M_z \rangle$", xlabel="time", ylim=(-1.1, 1.05))
        add_case_box(ax, label_from_run_id(str(row["run_id"])))
    fig.suptitle("Representative routes toward retained switching in the focused follow-up", fontsize=15, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig(FIG_DIR / "goal6_representative_dynamics.png", dpi=220)
    plt.close(fig)


def load_phase2_cases() -> list[dict]:
    return json.loads(PHASE2_MANIFEST.read_text())


def write_clean_family_avgmz_grid(cases: list[dict]) -> None:
    fig, axes = plt.subplots(3, 2, figsize=(12.5, 12), sharex=False, sharey=True)
    axes = axes.flatten()

    for ax, case in zip(axes, cases):
        pulse_df = read_observables(Path(case["pulse_results_dir"]) / "observables.tsv")
        control_df = read_observables(Path(case["control_results_dir"]) / "observables.tsv")
        t = pulse_df["t"].to_numpy(dtype=float)
        pulse_mz = pulse_df["avg_Mz"].to_numpy(dtype=float)
        control_mz = control_df["avg_Mz"].to_numpy(dtype=float)
        ax.set_ylim(-1.1, 1.05)
        add_pulse_window(ax, float(case["pulse_start"]), float(case["pulse_width"]))
        ax.plot(t, pulse_mz, lw=2.3, label="pulse", color=MZ_COLORS["pulse"])
        ax.plot(t, control_mz, lw=2.0, label="no pulse", color=MZ_COLORS["control"])
        ax.axhline(0.0, color="#94a3b8", lw=0.9, ls="--")
        ax.set_title(case["label"].replace("goal6_", ""), fontsize=10)
        style_axis(ax, ylabel=r"$\langle M_z \rangle$", xlabel="time", ylim=(-1.1, 1.05))

    axes[-1].axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    axes[-1].legend(handles, labels, loc="center", fontsize=10)
    fig.suptitle("Matched clean families: switching appears only in the pulsed twin", fontsize=15, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig(FIG_DIR / "clean_family_avgMz_grid.png", dpi=220)
    plt.close(fig)


def write_clean_family_excess_grid(cases: list[dict]) -> None:
    fig, axes = plt.subplots(3, 2, figsize=(12.5, 12), sharex=False, sharey=True)
    axes = axes.flatten()

    for ax, case in zip(axes, cases):
        pulse_df = read_observables(Path(case["pulse_results_dir"]) / "observables.tsv")
        control_df = read_observables(Path(case["control_results_dir"]) / "observables.tsv")
        t = pulse_df["t"].to_numpy(dtype=float)
        pulse_traces = sigma_component_absmax(pulse_df)
        control_traces = sigma_component_absmax(control_df)
        excess_traces = {comp: pulse_traces[comp] - control_traces[comp] for comp in COMPONENTS}

        for comp in COMPONENTS:
            ax.plot(
                t,
                excess_traces[comp],
                lw=2.3 if comp == "z" else 1.5,
                color=COMPONENT_COLORS[comp],
                label=fr"$\Delta \sigma_{{{comp}}}$",
            )
        ymax = max(float(np.max(np.abs(excess_traces[comp]))) for comp in COMPONENTS)
        ax.set_ylim(-0.03, max(0.21, ymax * 1.18))
        add_pulse_window(ax, float(case["pulse_start"]), float(case["pulse_width"]))
        ax.axhline(0.0, color="#94a3b8", lw=0.9, ls="--")
        ax.set_title(case["label"].replace("goal6_", ""), fontsize=10)
        style_axis(ax, ylabel=r"$\Delta \max_i |\sigma_{\alpha,i}|$", xlabel="time")

    axes[-1].axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    axes[-1].legend(handles, labels, loc="center", fontsize=10)
    fig.suptitle("Matched electronic excess: the pulse-specific signal turns on after onset", fontsize=15, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig(FIG_DIR / "clean_family_excess_grid.png", dpi=220)
    plt.close(fig)


def write_branchb_followup_grid() -> None:
    scan_df = pd.read_csv(BRANCHB_SCAN_SUMMARY, sep="\t")
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=False, sharey=True)
    axes = axes.flatten()

    for ax, (_, row) in zip(axes, scan_df.iterrows()):
        obs_path = Path(row["results_dir"]) / "observables.tsv"
        df = read_observables(obs_path)
        t = df["t"].to_numpy(dtype=float)
        avg_mz = df["avg_Mz"].to_numpy(dtype=float)
        ax.plot(t, avg_mz, lw=2.2, color=MZ_COLORS["branch_b"])
        ax.axhline(0.0, color="#334155", lw=1.0, ls="--")
        ax.set_title(str(row["label"]).replace("goal6_B_followup_", "B followup "), fontsize=10)
        style_axis(ax, ylabel=r"$\langle M_z \rangle$", xlabel="time", ylim=(-0.3, 1.05))

    fig.suptitle("Branch B no-pulse controls: baseline ambiguity is a family property", fontsize=15, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig(FIG_DIR / "branchb_nopulse_followup_avgMz_grid.png", dpi=220)
    plt.close(fig)


def write_rhoeq_relaxation_figure() -> None:
    df_1000 = read_observables(FIXED_LLG_T1000)
    df_2000 = read_observables(RHOEQ_T2000)

    trace_1000 = sigma_delta_global_absmax(df_1000)
    trace_2000 = sigma_delta_global_absmax(df_2000)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].plot(df_1000["t"], trace_1000, lw=2.2, color="#2563eb", label=r"fixed spins, $t_{\max}=1000$")
    axes[0].plot(df_2000["t"], trace_2000, lw=2.2, color="#ea580c", label=r"improved long-time run, $t_{\max}=2000$")
    axes[0].set_title("Global view")
    style_axis(axes[0], ylabel=r"$\max_{i,\alpha} |\sigma_{\Delta,\alpha i}(t)|$", xlabel="time")
    axes[0].legend(fontsize=9)

    late_mask_1000 = df_1000["t"].to_numpy(dtype=float) >= 600.0
    late_mask_2000 = df_2000["t"].to_numpy(dtype=float) >= 600.0
    axes[1].plot(df_1000.loc[late_mask_1000, "t"], trace_1000[late_mask_1000], lw=2.2, color="#2563eb")
    axes[1].plot(df_2000.loc[late_mask_2000, "t"], trace_2000[late_mask_2000], lw=2.2, color="#ea580c")
    axes[1].set_yscale("log")
    axes[1].set_title("Late-time decay")
    style_axis(axes[1], ylabel=r"$\max_{i,\alpha} |\sigma_{\Delta,\alpha i}(t)|$", xlabel="time")

    fig.suptitle("No-pulse electronic excess decays instead of mimicking a persistent driven response", fontsize=15, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig(FIG_DIR / "rhoeq_relaxation_tendency.png", dpi=220)
    plt.close(fig)


def write_report_metrics() -> None:
    phase2_rows = list(csv.DictReader(PHASE2_SUMMARY.open(), delimiter="\t"))
    branchb_rows = list(csv.DictReader(BRANCHB_SCAN_SUMMARY.open(), delimiter="\t"))
    goal6_rows = list(csv.DictReader(GOAL6_BRANCH_SUMMARY.open(), delimiter="\t"))
    goal5_selected = pd.read_csv(GOAL5_SELECTED, sep="\t")

    strongest = max(phase2_rows, key=lambda row: float(row["excess_sigma_delta_absmax_peak"]))
    metrics = {
        "goal5_selected_candidates": int(len(goal5_selected)),
        "goal5_strict_clean_in_pulse_candidates": int(goal5_selected["clean_pulse_window_strong_cross"].sum()),
        "goal6_clean_family_pair_count": len(phase2_rows),
        "goal6_branch_count": len(goal6_rows),
        "branchb_followup_case_count": len(branchb_rows),
        "strongest_clean_family_case": strongest["label"],
        "strongest_clean_family_peak": float(strongest["excess_sigma_delta_absmax_peak"]),
        "strongest_clean_family_peak_time": float(strongest["excess_sigma_delta_absmax_peak_time"]),
        "branchb_all_negative_excursion": all(float(row["avg_Mz_min"]) < 0.0 for row in branchb_rows),
        "branchb_any_negative_final": any(float(row["avg_Mz_final"]) < 0.0 for row in branchb_rows),
    }
    (PROJECT_DIR / "report_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")


def main() -> None:
    ensure_dir()
    cases = load_phase2_cases()
    write_goal5_representative_dynamics()
    write_goal6_representative_dynamics()
    write_clean_family_avgmz_grid(cases)
    write_clean_family_excess_grid(cases)
    write_branchb_followup_grid()
    write_rhoeq_relaxation_figure()
    write_report_metrics()


if __name__ == "__main__":
    main()
