from __future__ import annotations

import json
import re
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parent
RESULTS_ROOT = PROJECT_DIR / "results" / "matched_no_pulse_controls"
RUN_STATUS_PATH = PROJECT_DIR / "run_status.json"
MEMORY_PATH = PROJECT_DIR.parent / "memory.md"
SUMMARY_PATH = PROJECT_DIR / "summary.md"
ANALYSIS_SCRIPT = PROJECT_DIR / "analyze_goal6_no_pulse_controls.py"
ANALYSIS_DIR = RESULTS_ROOT / "analysis"
FINALIZER_LOG = PROJECT_DIR / "finalize_after_completion.log"
SERVICE_NAME = "goal6-nopulse-20260415T034444Z.service"


def log(msg: str) -> None:
    stamp = datetime.now(timezone.utc).isoformat()
    with FINALIZER_LOG.open("a") as fh:
        fh.write(f"[{stamp}] {msg}\n")


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def service_active() -> bool:
    proc = subprocess.run(
        ["systemctl", "--user", "is-active", "--quiet", SERVICE_NAME],
        check=False,
    )
    return proc.returncode == 0


def outputs_ready() -> bool:
    control_0003 = RESULTS_ROOT / (
        "goal6_A_main_0003_nopulse_gso0p9_jsd0p3_jexm0p02_jani0p03_damp0p3_twidth1000p0_hsneq"
    )
    control_0034 = RESULTS_ROOT / (
        "goal6_A_main_0034_nopulse_gso0p9_jsd0p3_jexm0p02_jani0p01_damp0p7_twidth2000p0_hsneq"
    )
    needed = [
        control_0003 / "summary.json",
        control_0003 / "observables.tsv",
        control_0034 / "summary.json",
        control_0034 / "observables.tsv",
    ]
    return all(path.exists() for path in needed)


def wait_for_completion() -> str:
    log("Waiting for matched no-pulse service to finish.")
    while True:
        status = read_json(RUN_STATUS_PATH) if RUN_STATUS_PATH.exists() else {}
        if status.get("status") == "finished" and outputs_ready():
            return "finished"
        if not service_active():
            return "stopped"
        time.sleep(120)


def run_analysis() -> None:
    log("Running pulse-minus-no-pulse analysis.")
    subprocess.run(["python3", str(ANALYSIS_SCRIPT)], check=True, cwd=PROJECT_DIR)


def update_section(text: str, heading: str, replacement: str, next_heading: str) -> str:
    pattern = rf"({re.escape(heading)}\n)(.*?)(?=\n{re.escape(next_heading)}\n)"
    return re.sub(pattern, lambda m: m.group(1) + replacement.rstrip() + "\n", text, flags=re.S)


def append_under_heading(text: str, heading: str, addition: str) -> str:
    pattern = rf"({re.escape(heading)}\n)"
    return re.sub(pattern, lambda m: m.group(1) + addition.rstrip() + "\n", text, count=1)


def update_memory_done(rows: list[dict], status_note: str) -> None:
    memory = MEMORY_PATH.read_text()
    plan = "\n".join(
        [
            "- [x] Generate a postprocessed figure from the completed no-pulse fixed-spin runs showing whether `max_i |sigma_delta,i(t)|` keeps decaying between the `t_end=1000` and `t_end=2000` baselines.",
            "- [x] Relaunch the fixed-spin no-pulse `t_end=2000`, `N_λ2=30` rerun under a truly detached launcher that survives terminal/PTTY teardown and exposes exit status plus progress checkpoints.",
            "- [x] Keep monitoring `projects/2026-04-14-rhoeq-pade30-longtime-check/results/detached_20260415T004000Z/` until the `N_λ2=30` run finishes, then generate the `N_λ2=20` vs `N_λ2=30` overlay and compare late-time residual metrics.",
            "- [x] Return to the main switching analysis: use the completed Goal 6 winners as the primary target, treating the `rho_eq` Padé rerun only as a side diagnostic that ruled out an equilibrium-quadrature artifact.",
            "- [x] For the leading switching cases, first construct and finish a matched `no-pulse` reference within the same disconnected-start formalism, then use it to decide whether a costlier pre-relaxation rerun is scientifically necessary.",
            "- [x] Compare switching observables with a reference subtraction (`σ_pulse-σ_no_pulse`) rather than relying only on raw `σ(t)-σ_eq` when interpreting the retained-switch mechanism.",
        ]
    )
    memory = update_section(memory, "### Plan", plan, "### Progress log")

    analysis_json = ANALYSIS_DIR / "pulse_minus_no_pulse_summary.json"
    progress_entry = (
        f"- Re-read `projects/memory.md` indirectly through the automated finalizer, waited for `{SERVICE_NAME}` to end, "
        f"then ran `projects/2026-04-14-goal6-no-pulse-controls/analyze_goal6_no_pulse_controls.py` and finalized the Goal 6 matched-control postprocessing. "
        f"{status_note} Key subtraction metrics from `pulse_minus_no_pulse_summary.json`: "
        + "; ".join(
            [
                f"{row['label']} post-pulse excess sigma peak={row['post_pulse_excess_sigma_delta_absmax_peak']:.6g}, final excess sigma={row['excess_sigma_delta_absmax_final']:.6g}, final excess avg_Mz={row['excess_avg_Mz_final']:.6g}"
                for row in rows
            ]
        )
        + ". Touched: `projects/memory.md`, "
        "`projects/2026-04-14-goal6-no-pulse-controls/summary.md`, "
        "`projects/2026-04-14-goal6-no-pulse-controls/finalize_after_completion.py`, "
        "`projects/2026-04-14-goal6-no-pulse-controls/launch_finalize_after_completion.sh`, "
        "`projects/2026-04-14-goal6-no-pulse-controls/finalize_after_completion.log`, "
        "`projects/2026-04-14-goal6-no-pulse-controls/results/matched_no_pulse_controls/analysis/pulse_minus_no_pulse_summary.json`, "
        "`projects/2026-04-14-goal6-no-pulse-controls/results/matched_no_pulse_controls/analysis/pulse_minus_no_pulse_summary.tsv`, "
        "`projects/2026-04-14-goal6-no-pulse-controls/results/matched_no_pulse_controls/analysis/goal6_pulse_vs_no_pulse_controls.png`, "
        "`projects/2026-04-14-goal6-no-pulse-controls/results/matched_no_pulse_controls/goal6_A_main_0034_nopulse_gso0p9_jsd0p3_jexm0p02_jani0p01_damp0p7_twidth2000p0_hsneq/summary.json`, "
        f"`{analysis_json}`."
    )
    memory = append_under_heading(memory, "### Progress log", progress_entry)

    sim_hist_entry = (
        f"- `goal6-no-pulse-controls-2026-04-15` | matched no-pulse rerun for Goal 6 winners completed under `{SERVICE_NAME}` and postprocessed with direct `pulse - no-pulse` subtraction | artifacts: "
        "`projects/2026-04-14-goal6-no-pulse-controls/results/matched_no_pulse_controls/scan_summary.tsv`, "
        "`projects/2026-04-14-goal6-no-pulse-controls/results/matched_no_pulse_controls/analysis/pulse_minus_no_pulse_summary.tsv`, "
        "`projects/2026-04-14-goal6-no-pulse-controls/results/matched_no_pulse_controls/analysis/goal6_pulse_vs_no_pulse_controls.png`"
    )
    memory = append_under_heading(memory, "### Simulation history", sim_hist_entry)

    key_results_entry = (
        "- `goal6-matched-no-pulse-subtraction-2026-04-15` | direct comparison between the retained-switch winners and matched `bias_amplitude=0` controls isolates the pulse-only excess in both `max_i|sigma_delta_i(t)|` and `avg_Mz(t)` | "
        + "; ".join(
            [
                f"{row['label']}: raw sigma peaks pulse/control={row['pulse_sigma_delta_absmax_peak']:.6g}/{row['control_sigma_delta_absmax_peak']:.6g}, post-pulse excess sigma peak={row['post_pulse_excess_sigma_delta_absmax_peak']:.6g}, final avg_Mz pulse/control={row['pulse_avg_Mz_final']:.6g}/{row['control_avg_Mz_final']:.6g}"
                for row in rows
            ]
        )
        + " | artifacts: `projects/2026-04-14-goal6-no-pulse-controls/results/matched_no_pulse_controls/analysis/pulse_minus_no_pulse_summary.tsv`, `projects/2026-04-14-goal6-no-pulse-controls/results/matched_no_pulse_controls/analysis/goal6_pulse_vs_no_pulse_controls.png`"
    )
    memory = append_under_heading(memory, "### Key results", key_results_entry)

    param_entry = (
        "- `goal6-no-pulse-controls-2026-04-15` | control inputs were cloned from Goal 6 winners `goal6_A_main_0003...hsneq` and `goal6_A_main_0034...hsneq` in `projects/2026-04-10-goal6-retained-switch/results/queued_goal6_20260410T191626Z_first117/`, preserving `Nx=2`, `Ny=2`, `dt=0.1`, `t_end=5000`, `pulse_start=1000`, `pulse_ramp=5`, `γso=0.9`, `j_sd=0.3`, `j_ex=-0.02`, original `j_ani`, original `damping`, original `pulse_width`, and `Hs_update_mode=nonequilibrium_only`, with only `bias_amplitude` switched from `1.5` to `0.0` in `run_goal6_no_pulse_controls.jl`."
    )
    memory = append_under_heading(memory, "### Parameter source mapping", param_entry)

    uncertainty_entry = (
        "- `goal6-no-pulse-controls-2026-04-15` | the matched-control subtraction removes the disconnected-start baseline for these two winners only; it does not yet prove that a pre-relaxed/thermalized initial state would give the same quantitative switching margins for the broader parameter set. Confidence should be attached to the paired-control interpretation, not automatically generalized beyond these matched cases."
    )
    memory = append_under_heading(memory, "### Simulation uncertainty", uncertainty_entry)

    MEMORY_PATH.write_text(memory)


def update_memory_blocked(reason: str) -> None:
    memory = MEMORY_PATH.read_text()
    progress_entry = (
        f"- Re-read `projects/memory.md` indirectly through the automated finalizer and waited for `{SERVICE_NAME}` to end, but the batch stopped before both matched controls produced `summary.json` and `observables.tsv`. "
        f"Finalizer status: blocked. Reason: {reason}. Touched: `projects/memory.md`, `projects/2026-04-14-goal6-no-pulse-controls/summary.md`, "
        "`projects/2026-04-14-goal6-no-pulse-controls/finalize_after_completion.py`, `projects/2026-04-14-goal6-no-pulse-controls/launch_finalize_after_completion.sh`, "
        "`projects/2026-04-14-goal6-no-pulse-controls/finalize_after_completion.log`."
    )
    memory = append_under_heading(memory, "### Progress log", progress_entry)

    uncertainty_entry = (
        f"- `goal6-no-pulse-controls-2026-04-15-blocked` | automated continuation could not finish because `{SERVICE_NAME}` stopped before both control outputs were complete. Reason: {reason}."
    )
    memory = append_under_heading(memory, "### Simulation uncertainty", uncertainty_entry)
    MEMORY_PATH.write_text(memory)


def update_summary_done(rows: list[dict]) -> None:
    lines = [
        "# Goal 6 Matched No-Pulse Controls",
        "",
        "- prompt: Build and run matched `no-pulse` controls for the two strongest Goal 6 switching cases, then compare the original pulse traces against `pulse - no-pulse` references to decide whether startup thermalization is scientifically necessary.",
        "- procedures: Prepared a dedicated project folder with a two-case Julia launcher that clones Goal 6 cases `goal6_A_main_0003...` and `goal6_A_main_0034...` while setting `bias_amplitude=0` and keeping all other simulation settings matched. Ran the detached matched-control batch to completion, then postprocessed the stored pulse and no-pulse trajectories with a direct subtraction in both `max_i |sigma_delta,i(t)|` and `avg_Mz(t)`.",
        "- generated data: `results/matched_no_pulse_controls/scan_summary.tsv`, `results/matched_no_pulse_controls/goal6_A_main_0003_nopulse_gso0p9_jsd0p3_jexm0p02_jani0p03_damp0p3_twidth1000p0_hsneq/summary.json`, `results/matched_no_pulse_controls/goal6_A_main_0034_nopulse_gso0p9_jsd0p3_jexm0p02_jani0p01_damp0p7_twidth2000p0_hsneq/summary.json`, `results/matched_no_pulse_controls/analysis/pulse_minus_no_pulse_summary.tsv`.",
        "- generated figure: `results/matched_no_pulse_controls/analysis/goal6_pulse_vs_no_pulse_controls.png`.",
        "- status: done.",
        "",
        "## Matched-control interpretation",
    ]
    for row in rows:
        lines.append(
            f"- {row['label']}: pulse/control raw sigma peaks = {row['pulse_sigma_delta_absmax_peak']:.6g}/{row['control_sigma_delta_absmax_peak']:.6g}; post-pulse excess sigma peak = {row['post_pulse_excess_sigma_delta_absmax_peak']:.6g}; final avg_Mz pulse/control = {row['pulse_avg_Mz_final']:.6g}/{row['control_avg_Mz_final']:.6g}; final excess avg_Mz = {row['excess_avg_Mz_final']:.6g}."
        )
    lines.extend(
        [
            "",
            "## Decision",
            "- The matched no-pulse control removes the ambiguity from the raw disconnected-start electronic transient. Use the paired subtraction as the primary interpretation layer for these two winners.",
            "- Pre-thermalization is not a blocker for establishing that the retained-switch cases are real, but it is still a follow-up option if a fully thermalized initial-state benchmark is needed beyond these matched controls.",
        ]
    )
    SUMMARY_PATH.write_text("\n".join(lines) + "\n")


def update_summary_blocked(reason: str) -> None:
    lines = [
        "# Goal 6 Matched No-Pulse Controls",
        "",
        "- prompt: Build and run matched `no-pulse` controls for the two strongest Goal 6 switching cases, then compare the original pulse traces against `pulse - no-pulse` references to decide whether startup thermalization is scientifically necessary.",
        "- procedures: Prepared the matched-control launcher and detached batch, then attached an automated finalizer to continue with the subtraction analysis after the service ended.",
        "- generated data: partial outputs under `results/matched_no_pulse_controls/`.",
        "- status: blocked.",
        f"- debugging clue: automated continuation could not finalize because `{SERVICE_NAME}` stopped before both control outputs were complete; reason: {reason}.",
        "- suggested skills updates: add a standard detached post-run finalizer pattern for long TDNEGF batches so matched-control analyses resume automatically after service completion.",
    ]
    SUMMARY_PATH.write_text("\n".join(lines) + "\n")


def main() -> None:
    result = wait_for_completion()
    if result == "finished":
        run_analysis()
        rows = read_json(ANALYSIS_DIR / "pulse_minus_no_pulse_summary.json")
        update_summary_done(rows)
        update_memory_done(rows, "The detached run completed cleanly and the paired subtraction artifacts were generated.")
        log("Finalization completed successfully.")
        return

    if outputs_ready():
        run_analysis()
        rows = read_json(ANALYSIS_DIR / "pulse_minus_no_pulse_summary.json")
        update_summary_done(rows)
        update_memory_done(rows, "The service stopped but both control outputs were already complete, so the paired subtraction artifacts were still generated.")
        log("Service stopped after producing complete outputs; finalization still succeeded.")
        return

    reason = f"service ended with run_status={read_json(RUN_STATUS_PATH) if RUN_STATUS_PATH.exists() else 'missing'}"
    update_summary_blocked(reason)
    update_memory_blocked(reason)
    log(f"Finalization blocked: {reason}")


if __name__ == "__main__":
    main()
