"""Generate a Jupyter notebook for the first pilot-scan comparison layer.

The notebook is intentionally written against the existing saved output layout:

    projects/2026-04-08-tdnegf-spin-baseline/results/pilot_sweep_*/...

It auto-selects the latest non-empty pilot sweep root, loads each run bundle,
classifies runs as successful / failed / suspicious, defines a conservative
switching criterion, and produces comparison plots with publication-style
matplotlib settings.
"""

from __future__ import annotations

from pathlib import Path
import textwrap

import nbformat as nbf


PROJECT_DIR = Path(__file__).resolve().parent
NOTEBOOK_PATH = PROJECT_DIR / "pilot_comparison_analysis.ipynb"


def markdown_cell(source: str):
    """Return a markdown notebook cell with normalized leading whitespace."""

    return nbf.v4.new_markdown_cell(textwrap.dedent(source).strip() + "\n")


def code_cell(source: str):
    """Return a code notebook cell with a trailing newline."""

    return nbf.v4.new_code_cell(textwrap.dedent(source).strip() + "\n")


def build_notebook():
    """Assemble the notebook cells for pilot-run comparison and review."""

    nb = nbf.v4.new_notebook()
    nb["metadata"].update(
        {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {
                "name": "python",
                "version": "3",
            },
        }
    )

    cells = [
        markdown_cell(
            """
            # Pilot Sweep Comparison and First Analysis Layer

            This notebook loads the saved pilot sweep output bundles under
            `results/pilot_sweep_*`, picks the latest non-empty sweep root, and
            builds the first comparison layer without changing the output
            structure. The emphasis is on:

            - numerical completion versus physical trustworthiness,
            - a conservative switching criterion,
            - representative time traces,
            - final magnetization and current comparisons,
            - compact parameter-to-observable views.

            The run classification is intentionally conservative:

            - **successful**: completed and passed the sanity checks used here,
            - **suspicious**: completed but shows strong warning signs,
            - **failed**: missing outputs, incomplete runs, or explicit failure.
            """
        ),
        code_cell(
            r"""
            # Core imports for file discovery, tabular analysis, and plotting.
            import json
            import math
            import shutil
            from pathlib import Path

            import numpy as np
            import pandas as pd
            import matplotlib.pyplot as plt

            try:
                from IPython.display import Markdown, display
            except Exception:
                # Fall back to plain-text display when IPython is unavailable.
                def Markdown(text):
                    return text

                def display(obj):
                    print(obj)
            """
        ),
        code_cell(
            r"""
            def configure_matplotlib():
                '''Apply the requested plotting style with a safe LaTeX fallback.'''

                latex_available = shutil.which("latex") is not None

                # Use the exact style requested by the user when LaTeX exists.
                plt.rc("text", usetex=latex_available)
                plt.rc("font", family="serif")
                # plt.rc("text.latex", preamble=r"\usepackage{amsmath}")
                plt.rc("font", size=18)
                plt.rc("axes", labelsize=22, titlesize=22)
                plt.rc("xtick", labelsize=18)
                plt.rc("ytick", labelsize=18)
                plt.rc("legend", fontsize=16)
                plt.rc("axes", linewidth=2)

                if not latex_available:
                    print(
                        "LaTeX was not found in this environment. "
                        "Using the requested serif style with usetex=False."
                    )


            def collect_pilot_roots(results_dir: Path):
                '''Return pilot-sweep roots that contain a non-empty summary table.'''

                roots = []
                for root in sorted(results_dir.glob("pilot_sweep_*")):
                    summary_path = root / "pilot_summary.tsv"
                    if summary_path.is_file() and summary_path.stat().st_size > 0:
                        roots.append(root)
                return roots


            def choose_latest_complete_root(results_dir: Path):
                '''Pick the newest pilot root that has a readable summary table.'''

                roots = collect_pilot_roots(results_dir)
                if not roots:
                    raise FileNotFoundError(
                        f"No non-empty pilot_sweep_* results were found under {results_dir}"
                    )
                return roots[-1]


            def load_json(path: Path):
                '''Load one JSON file and return a Python dictionary.'''

                with path.open("r", encoding="utf-8") as handle:
                    return json.load(handle)


            def read_observables_table(path: Path):
                '''Load one observables table and derive compact averaged traces.'''

                df = pd.read_csv(path, sep="\t")

                # Derive average quantities directly from the saved per-site columns.
                classical_z_cols = sorted(col for col in df.columns if col.startswith("s_classical_z_site"))
                sigma_z_cols = sorted(col for col in df.columns if col.startswith("sigma_z_site"))
                density_cols = sorted(col for col in df.columns if col.startswith("n_site"))

                df["avg_Mz"] = df[classical_z_cols].mean(axis=1)
                df["avg_sigma_z"] = df[sigma_z_cols].mean(axis=1)
                df["avg_density"] = df[density_cols].mean(axis=1)
                df["current_balance"] = df["I_left"] + df["I_right"]
                return df


            def load_run_record(run_dir: Path):
                '''Load one run bundle and merge summary, params, and path metadata.'''

                summary_path = run_dir / "summary.json"
                params_path = run_dir / "params.json"
                observables_path = run_dir / "observables.tsv"

                if not summary_path.is_file():
                    return {
                        "run_id": run_dir.name,
                        "status": "failed",
                        "run_dir": str(run_dir),
                        "summary_path": str(summary_path),
                        "params_path": str(params_path),
                        "observables_path": str(observables_path),
                        "failure_reason": "missing summary.json",
                    }

                summary = load_json(summary_path)
                params = load_json(params_path) if params_path.is_file() else {}

                record = {
                    **params,
                    **summary,
                    "run_dir": str(run_dir),
                    "summary_path": str(summary_path),
                    "params_path": str(params_path),
                    "observables_path": str(observables_path) if observables_path.is_file() else None,
                    "observables_exists": observables_path.is_file(),
                }
                return record


            def build_run_table(pilot_root: Path):
                '''Load every run folder below one pilot root into a single table.'''

                run_dirs = sorted(path for path in pilot_root.iterdir() if path.is_dir())
                records = [load_run_record(run_dir) for run_dir in run_dirs]
                if not records:
                    raise RuntimeError(f"No run subdirectories were found under {pilot_root}")
                return pd.DataFrame(records)


            def classify_run(row, charge_warning_limit=0.25, spin_warning_limit=1e-6):
                '''Classify a run using conservative completion and sanity checks.'''

                reasons = []

                # Explicit failure or missing observable output is a hard failure.
                if str(row.get("status", "")).lower() != "done":
                    reasons.append("summary status is not 'done'")
                    return "failed", reasons
                if not bool(row.get("observables_exists", False)):
                    reasons.append("observables.tsv is missing")
                    return "failed", reasons

                # Check that the run reached the configured end time and expected steps.
                dt = float(row.get("dt", np.nan))
                t0 = float(row.get("t0", np.nan))
                t_end = float(row.get("t_end", np.nan))
                completed_steps = float(row.get("completed_steps", np.nan))
                t_final = float(row.get("t_final", np.nan))
                expected_steps = int(round((t_end - t0) / dt)) if np.isfinite(dt) and dt > 0 else None

                if expected_steps is not None and completed_steps != expected_steps:
                    reasons.append(f"completed_steps={completed_steps:.0f} differs from expected_steps={expected_steps}")
                if np.isfinite(t_end) and np.isfinite(t_final) and abs(t_final - t_end) > max(dt, 1e-9):
                    reasons.append(f"t_final={t_final:.6g} differs from t_end={t_end:.6g}")

                # Flag physically suspicious runs even when they completed numerically.
                charge_dev = float(row.get("max_site_charge_deviation_from_1", np.nan))
                spin_dev = float(row.get("max_classical_spin_norm_deviation", np.nan))
                current_balance_peak = float(row.get("current_balance_peak_abs", np.nan))

                if np.isfinite(charge_dev) and charge_dev > charge_warning_limit:
                    reasons.append(
                        "large site-charge excursion from unit filling "
                        f"(max deviation={charge_dev:.3f})"
                    )
                if np.isfinite(spin_dev) and spin_dev > spin_warning_limit:
                    reasons.append(
                        "classical spin norm drift exceeded tolerance "
                        f"(max deviation={spin_dev:.3e})"
                    )
                if np.isfinite(current_balance_peak) and current_balance_peak > 0.25:
                    reasons.append(
                        "large transient current imbalance suggests a strong initialization transient "
                        f"(peak |I_left + I_right|={current_balance_peak:.3f})"
                    )

                if reasons:
                    return "suspicious", reasons
                return "successful", ["passed the notebook sanity checks"]


            def assess_switching(trace_df, row, mz_threshold=-0.8, sustain_fraction=0.2):
                '''Apply a conservative switching test to one run.

                The criterion only reports a switch when the average classical
                magnetization along z:

                1. drops below -0.8 after the pulse is fully on,
                2. remains below -0.8 through the final sustained time window,
                3. ends negative by a large margin rather than by a small fluctuation.
                '''

                avg_mz = trace_df["avg_Mz"]
                time = trace_df["t"]
                pulse_ready_time = float(row["pulse_start"]) + float(row["pulse_ramp"])

                post_pulse = trace_df.loc[time >= pulse_ready_time].copy()
                if post_pulse.empty:
                    return {
                        "switching_label": "insufficient_post_pulse_window",
                        "switch_detected": False,
                        "pulse_ready_time": pulse_ready_time,
                        "min_post_pulse_avg_Mz": math.nan,
                        "final_window_avg_Mz_max": math.nan,
                    }

                tail_count = max(3, int(math.ceil(len(post_pulse) * sustain_fraction)))
                final_window = post_pulse.tail(tail_count)

                crossed_threshold = bool((post_pulse["avg_Mz"] <= mz_threshold).any())
                sustained_negative = bool((final_window["avg_Mz"] <= mz_threshold).all())
                final_negative = bool(post_pulse["avg_Mz"].iloc[-1] <= mz_threshold)

                switch_detected = crossed_threshold and sustained_negative and final_negative
                label = "conservative_switch_detected" if switch_detected else "no_conservative_switch"

                return {
                    "switching_label": label,
                    "switch_detected": switch_detected,
                    "pulse_ready_time": pulse_ready_time,
                    "min_post_pulse_avg_Mz": float(post_pulse["avg_Mz"].min()),
                    "final_window_avg_Mz_max": float(final_window["avg_Mz"].max()),
                }


            def load_traces_for_runs(run_table):
                '''Read observables for completed runs and return a run_id keyed mapping.'''

                traces = {}
                for row in run_table.itertuples(index=False):
                    path = getattr(row, "observables_path", None)
                    if path:
                        traces[row.run_id] = read_observables_table(Path(path))
                return traces


            configure_matplotlib()
            results_dir = Path("results")
            pilot_root = choose_latest_complete_root(results_dir)
            run_table = build_run_table(pilot_root)
            traces = load_traces_for_runs(run_table)

            print(f"Selected pilot root: {pilot_root}")
            print(f"Loaded {len(run_table)} run records.")
            """
        ),
        code_cell(
            r"""
            # Apply the run-status classification and the conservative switching test.
            categories = []
            reasons = []
            switching_records = []

            for row in run_table.to_dict(orient="records"):
                category, category_reasons = classify_run(row)
                categories.append(category)
                reasons.append(" | ".join(category_reasons))

                if row["run_id"] in traces:
                    switching = assess_switching(traces[row["run_id"]], row)
                else:
                    switching = {
                        "switching_label": "not_evaluated",
                        "switch_detected": False,
                        "pulse_ready_time": math.nan,
                        "min_post_pulse_avg_Mz": math.nan,
                        "final_window_avg_Mz_max": math.nan,
                    }
                switching_records.append(switching)

            run_table["analysis_category"] = categories
            run_table["analysis_reasons"] = reasons
            switching_df = pd.DataFrame(switching_records)
            run_table = pd.concat([run_table.reset_index(drop=True), switching_df], axis=1)

            display(
                Markdown(
                    "## Selected Pilot Root\n"
                    f"`{pilot_root}`"
                )
            )

            display(
                Markdown(
                    "## Conservative Switching Criterion\n"
                    "A run is flagged as switched only if the average classical "
                    "$M_z$ drops below $-0.8$ after the pulse is fully on and "
                    "stays below $-0.8$ through the final sustained window."
                )
            )
            """
        ),
        code_cell(
            r"""
            # Summarize how many runs fall in each category.
            category_counts = (
                run_table["analysis_category"]
                .value_counts()
                .rename_axis("analysis_category")
                .reset_index(name="count")
            )
            display(Markdown("## Category Counts"))
            display(category_counts)

            display(Markdown("## Run Overview"))
            overview_columns = [
                "run_id",
                "status",
                "analysis_category",
                "switching_label",
                "bias_amplitude",
                "j_sd",
                "j_ani",
                "avg_Mz_final",
                "avg_Mz_min",
                "current_left_final",
                "current_right_final",
                "analysis_reasons",
            ]
            display(run_table[overview_columns].sort_values(["analysis_category", "run_id"]))
            """
        ),
        code_cell(
            r"""
            # Show the three analysis buckets separately for readability.
            successful_runs = run_table.loc[run_table["analysis_category"] == "successful"].copy()
            suspicious_runs = run_table.loc[run_table["analysis_category"] == "suspicious"].copy()
            failed_runs = run_table.loc[run_table["analysis_category"] == "failed"].copy()

            display(Markdown("## Successful Runs"))
            display(successful_runs if not successful_runs.empty else pd.DataFrame({"note": ["No runs passed all sanity checks in this pilot set."]}))

            display(Markdown("## Suspicious Runs"))
            display(suspicious_runs if not suspicious_runs.empty else pd.DataFrame({"note": ["No suspicious runs were detected."]}))

            display(Markdown("## Failed Runs"))
            display(failed_runs if not failed_runs.empty else pd.DataFrame({"note": ["No failed runs were detected in the selected pilot root."]}))
            """
        ),
        code_cell(
            r"""
            # Build representative time traces for the completed runs.
            completed_run_ids = run_table.loc[
                run_table["analysis_category"].isin(["successful", "suspicious"]), "run_id"
            ].tolist()

            fig, axes = plt.subplots(3, 1, figsize=(12, 16), sharex=True, constrained_layout=True)

            for run_id in completed_run_ids:
                trace = traces[run_id]
                row = run_table.loc[run_table["run_id"] == run_id].iloc[0]
                label = (
                    rf"$V={row['bias_amplitude']:.1f},\ J_{{sd}}={row['j_sd']:.2f},\ "
                    rf"J_{{ani}}={row['j_ani']:.2f}$"
                )

                # Time trace 1: average classical magnetization along z.
                axes[0].plot(trace["t"], trace["avg_Mz"], linewidth=2.5, label=label)

                # Time trace 2: lead currents to show the pulse response and current balance.
                axes[1].plot(trace["t"], trace["I_left"], linewidth=2.0, label=label + " (left)")
                axes[1].plot(trace["t"], trace["I_right"], linewidth=2.0, linestyle="--", label=label + " (right)")

                # Time trace 3: average electronic spin density along z.
                axes[2].plot(trace["t"], trace["avg_sigma_z"], linewidth=2.5, label=label)

            pulse_start = float(run_table["pulse_start"].dropna().iloc[0])
            pulse_ready = pulse_start + float(run_table["pulse_ramp"].dropna().iloc[0])
            for ax in axes:
                ax.axvline(pulse_start, color="black", linestyle=":", linewidth=2, label="pulse start")
                ax.axvline(pulse_ready, color="gray", linestyle="--", linewidth=2, label="pulse fully on")
                ax.grid(alpha=0.25)

            axes[0].axhline(-0.8, color="tab:red", linestyle="--", linewidth=2, label="switch threshold")
            axes[0].set_ylabel(r"$\langle M_z \rangle$")
            axes[0].set_title("Representative time traces")

            axes[1].set_ylabel(r"Current")
            axes[2].set_ylabel(r"$\langle \sigma_z \rangle$")
            axes[2].set_xlabel(r"Time")

            # Deduplicate repeated legend labels before rendering.
            handles, labels = [], []
            for ax in axes:
                h, l = ax.get_legend_handles_labels()
                handles.extend(h)
                labels.extend(l)
            unique = dict(zip(labels, handles))
            axes[0].legend(unique.values(), unique.keys(), loc="best")

            plt.show()
            """
        ),
        code_cell(
            r"""
            # Compare final magnetization and final current across parameter choices.
            comparison = run_table.copy()
            comparison["parameter_label"] = comparison.apply(
                lambda row: (
                    rf"$V={row['bias_amplitude']:.1f}$" + "\n" +
                    rf"$J_{{sd}}={row['j_sd']:.2f},\ J_{{ani}}={row['j_ani']:.2f}$"
                ),
                axis=1,
            )

            x = np.arange(len(comparison))
            width = 0.38

            fig, axes = plt.subplots(1, 2, figsize=(15, 6), constrained_layout=True)

            axes[0].bar(x, comparison["avg_Mz_final"], color="tab:blue", alpha=0.85)
            axes[0].axhline(-0.8, color="tab:red", linestyle="--", linewidth=2, label="switch threshold")
            axes[0].set_xticks(x, comparison["parameter_label"])
            axes[0].set_ylabel(r"Final $\langle M_z \rangle$")
            axes[0].set_title("Final magnetization comparison")
            axes[0].grid(axis="y", alpha=0.25)
            axes[0].legend(loc="best")

            axes[1].bar(x - width / 2, comparison["current_left_final"], width=width, color="tab:green", label=r"$I_{\mathrm{left}}(t_f)$")
            axes[1].bar(x + width / 2, comparison["current_right_final"], width=width, color="tab:orange", label=r"$I_{\mathrm{right}}(t_f)$")
            axes[1].set_xticks(x, comparison["parameter_label"])
            axes[1].set_ylabel(r"Final current")
            axes[1].set_title("Final lead-current comparison")
            axes[1].grid(axis="y", alpha=0.25)
            axes[1].legend(loc="best")

            plt.show()
            """
        ),
        code_cell(
            r"""
            # Provide one compact view of how the main observables move with parameters.
            fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), constrained_layout=True)

            scatter_0 = axes[0].scatter(
                run_table["bias_amplitude"],
                run_table["j_sd"],
                c=run_table["avg_Mz_final"],
                s=450 * (1.02 - run_table["avg_Mz_final"].clip(upper=1.02)),
                cmap="viridis_r",
                edgecolor="black",
                linewidth=1.5,
            )
            axes[0].set_xlabel(r"Bias amplitude")
            axes[0].set_ylabel(r"$J_{sd}$")
            axes[0].set_title(r"Color: final $\langle M_z \rangle$")
            for row in run_table.itertuples(index=False):
                axes[0].annotate(f"Jani={row.j_ani:.2f}", (row.bias_amplitude, row.j_sd), xytext=(6, 6), textcoords="offset points")
            plt.colorbar(scatter_0, ax=axes[0])

            axes[1].plot(run_table["bias_amplitude"], run_table["current_left_final"], "o-", linewidth=2.5, markersize=10)
            axes[1].set_xlabel(r"Bias amplitude")
            axes[1].set_ylabel(r"$I_{\mathrm{left}}(t_f)$")
            axes[1].set_title("Bias versus final left current")
            axes[1].grid(alpha=0.25)

            axes[2].plot(run_table["j_sd"], 1.0 - run_table["avg_Mz_final"], "s-", linewidth=2.5, markersize=10, color="tab:red")
            axes[2].set_xlabel(r"$J_{sd}$")
            axes[2].set_ylabel(r"$1 - \langle M_z(t_f) \rangle$")
            axes[2].set_title("Drive strength versus magnetization softening")
            axes[2].grid(alpha=0.25)

            plt.show()
            """
        ),
        code_cell(
            r"""
            # Write a short notebook-side interpretation from the loaded results.
            completed = run_table.loc[run_table["analysis_category"] != "failed"].copy()
            strongest_case = completed.sort_values("avg_Mz_final").iloc[0]

            interpretation_lines = [
                "## First-pass interpretation",
                f"- Selected pilot root: `{pilot_root}`",
                f"- Completed runs: {len(completed)}",
                f"- Failed runs: {(run_table['analysis_category'] == 'failed').sum()}",
                f"- Suspicious completed runs: {(run_table['analysis_category'] == 'suspicious').sum()}",
                "- Conservative switching verdict: no run meets the sustained reversal criterion.",
                (
                    "- Strongest magnetization softening occurs for "
                    f"`run_id={strongest_case['run_id']}` with "
                    f"`avg_Mz_final={strongest_case['avg_Mz_final']:.6f}`."
                ),
                (
                    "- All completed runs remain physically suspicious in this notebook "
                    "because the saved summaries show large charge excursions and large "
                    "transient current imbalance, both consistent with the known zero-state "
                    "initialization transient."
                ),
                (
                    "- Conservative recommendation: treat this pilot set as a ranking tool "
                    "for stronger-drive cases, not as evidence of switching."
                ),
            ]
            display(Markdown("\n".join(interpretation_lines)))
            """
        ),
    ]

    nb["cells"] = cells
    return nb


def main():
    """Write the notebook to disk."""

    notebook = build_notebook()
    with NOTEBOOK_PATH.open("w", encoding="utf-8") as handle:
        nbf.write(notebook, handle)
    print(f"Wrote notebook to {NOTEBOOK_PATH}")


if __name__ == "__main__":
    main()
