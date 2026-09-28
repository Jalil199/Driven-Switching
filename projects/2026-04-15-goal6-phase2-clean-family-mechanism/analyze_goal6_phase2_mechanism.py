from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


PROJECT_DIR = Path(__file__).resolve().parent
MANIFEST_PATH = PROJECT_DIR / "cases_manifest.json"
RESULTS_ROOT = PROJECT_DIR / "results" / "phase2_mechanism_clean_families"
ANALYSIS_DIR = RESULTS_ROOT / "analysis"
SUMMARY_JSON = ANALYSIS_DIR / "phase2_mechanism_summary.json"
SUMMARY_TSV = ANALYSIS_DIR / "phase2_mechanism_summary.tsv"
REVIEW_MD = ANALYSIS_DIR / "phase2_mechanism_review.md"
CASE_MANIFEST_TSV = RESULTS_ROOT / "case_manifest.tsv"
SCAN_CONFIG_JSON = RESULTS_ROOT / "scan_config.json"

SIGMA_RE = re.compile(r"^sigma_delta_([xyz])_site(\d+)$")
COMPONENTS = ("x", "y", "z")
DROP_THRESHOLD = 0.9


def read_tsv(path: Path) -> dict[str, np.ndarray]:
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        columns: dict[str, list[float]] = {}
        for row in reader:
            for key, value in row.items():
                columns.setdefault(key, []).append(float(value))
    return {key: np.asarray(values, dtype=float) for key, values in columns.items()}


def write_case_manifest(cases: list[dict]) -> None:
    headers = [
        "label",
        "family",
        "Hs_update_mode",
        "pulse_start",
        "pulse_width",
        "pulse_run_id",
        "pulse_results_dir",
        "control_run_id",
        "control_results_dir",
    ]
    with CASE_MANIFEST_TSV.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerows([{key: case[key] for key in headers} for case in cases])


def write_scan_config(cases: list[dict]) -> None:
    config = {
        "scan_id": RESULTS_ROOT.name,
        "analysis_type": "goal6_phase2_pulse_minus_no_pulse_mechanism_clean_families",
        "requested_case_count": len(cases),
        "drop_threshold_avg_Mz": DROP_THRESHOLD,
        "notes": [
            "Fase 2 only for clean families from Step 1 plus the two mandatory baseline pairs.",
            "No new TDNEGF simulations were launched; this is post-processing of completed pulse/control pairs.",
        ],
    }
    SCAN_CONFIG_JSON.write_text(json.dumps(config, indent=2) + "\n")


def first_crossing_time(t: np.ndarray, values: np.ndarray, threshold: float) -> float | None:
    hits = np.flatnonzero(values <= threshold)
    if hits.size == 0:
        return None
    return float(t[hits[0]])


def sigma_delta_cube(data: dict[str, np.ndarray]) -> tuple[list[int], dict[str, np.ndarray]]:
    site_ids: set[int] = set()
    component_map: dict[str, dict[int, np.ndarray]] = {component: {} for component in COMPONENTS}
    for key, values in data.items():
        match = SIGMA_RE.match(key)
        if not match:
            continue
        component, site_str = match.groups()
        site = int(site_str)
        site_ids.add(site)
        component_map[component][site] = values
    ordered_sites = sorted(site_ids)
    cubes = {
        component: np.vstack([component_map[component][site] for site in ordered_sites])
        for component in COMPONENTS
    }
    return ordered_sites, cubes


def component_absmax_traces(cubes: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    return {component: np.max(np.abs(cube), axis=0) for component, cube in cubes.items()}


def global_absmax_from_cubes(cubes: dict[str, np.ndarray]) -> np.ndarray:
    return np.max(np.vstack([np.abs(cube) for cube in cubes.values()]), axis=0)


def pulse_window_mask(t: np.ndarray, pulse_start: float, pulse_width: float) -> np.ndarray:
    window_end = min(float(t[-1]), pulse_start + pulse_width + 500.0)
    return (t >= max(0.0, pulse_start - 200.0)) & (t <= window_end)


def summarise_case(case: dict, pulse_data: dict[str, np.ndarray], control_data: dict[str, np.ndarray]) -> tuple[dict[str, object], dict[str, object]]:
    t = pulse_data["t"]
    pulse_start = float(case["pulse_start"])
    pulse_width = float(case["pulse_width"])
    ordered_sites, pulse_cubes = sigma_delta_cube(pulse_data)
    _, control_cubes = sigma_delta_cube(control_data)
    excess_cubes = {component: pulse_cubes[component] - control_cubes[component] for component in COMPONENTS}

    pulse_component_absmax = component_absmax_traces(pulse_cubes)
    control_component_absmax = component_absmax_traces(control_cubes)
    excess_component_absmax = component_absmax_traces(excess_cubes)

    pulse_global_absmax = global_absmax_from_cubes(pulse_cubes)
    control_global_absmax = global_absmax_from_cubes(control_cubes)
    excess_global_absmax = global_absmax_from_cubes(excess_cubes)

    pre_mask = t < pulse_start
    post_mask = ~pre_mask
    pulse_avg_mz = pulse_data["avg_Mz"]
    control_avg_mz = control_data["avg_Mz"]
    excess_avg_mz = pulse_avg_mz - control_avg_mz

    peak_index = int(np.argmax(excess_global_absmax))
    peak_component = max(COMPONENTS, key=lambda component: float(np.max(excess_component_absmax[component])))
    strong_drop_time = first_crossing_time(t, pulse_avg_mz, DROP_THRESHOLD)
    zero_cross_time = first_crossing_time(t, pulse_avg_mz, 0.0)

    summary_row: dict[str, object] = {
        "label": case["label"],
        "family": case["family"],
        "Hs_update_mode": case["Hs_update_mode"],
        "pulse_run_id": case["pulse_run_id"],
        "control_run_id": case["control_run_id"],
        "pulse_start": pulse_start,
        "pulse_width": pulse_width,
        "pulse_avg_Mz_final": float(pulse_avg_mz[-1]),
        "control_avg_Mz_final": float(control_avg_mz[-1]),
        "pulse_avg_Mz_min": float(np.min(pulse_avg_mz)),
        "control_avg_Mz_min": float(np.min(control_avg_mz)),
        "excess_avg_Mz_final": float(excess_avg_mz[-1]),
        "pulse_sigma_delta_absmax_peak": float(np.max(pulse_global_absmax)),
        "control_sigma_delta_absmax_peak": float(np.max(control_global_absmax)),
        "excess_sigma_delta_absmax_peak": float(excess_global_absmax[peak_index]),
        "excess_sigma_delta_absmax_peak_time": float(t[peak_index]),
        "peak_excess_component": peak_component,
        "pre_pulse_excess_sigma_delta_absmax_peak": float(np.max(excess_global_absmax[pre_mask])),
        "post_pulse_excess_sigma_delta_absmax_peak": float(np.max(excess_global_absmax[post_mask])),
        "first_strong_drop_time_avg_Mz_le_0p9": strong_drop_time,
        "first_zero_cross_time_avg_Mz": zero_cross_time,
    }

    for component in COMPONENTS:
        component_peak_index = int(np.argmax(excess_component_absmax[component]))
        summary_row[f"excess_sigma_delta_{component}_peak"] = float(excess_component_absmax[component][component_peak_index])
        summary_row[f"excess_sigma_delta_{component}_peak_time"] = float(t[component_peak_index])

    series = {
        "t": t,
        "ordered_sites": ordered_sites,
        "pulse_component_absmax": pulse_component_absmax,
        "control_component_absmax": control_component_absmax,
        "excess_component_absmax": excess_component_absmax,
        "pulse_global_absmax": pulse_global_absmax,
        "control_global_absmax": control_global_absmax,
        "excess_global_absmax": excess_global_absmax,
        "pulse_avg_mz": pulse_avg_mz,
        "control_avg_mz": control_avg_mz,
        "excess_avg_mz": excess_avg_mz,
        "excess_cubes": excess_cubes,
        "pulse_start": pulse_start,
        "pulse_width": pulse_width,
        "strong_drop_time": strong_drop_time,
        "zero_cross_time": zero_cross_time,
        "peak_excess_time": float(t[peak_index]),
    }
    return summary_row, series


def write_summary(rows: list[dict[str, object]]) -> None:
    SUMMARY_JSON.write_text(json.dumps(rows, indent=2) + "\n")
    headers = list(rows[0].keys())
    with SUMMARY_TSV.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def add_marker_lines(ax: plt.Axes, series: dict[str, object]) -> None:
    pulse_start = float(series["pulse_start"])
    ax.axvline(pulse_start, color="k", lw=0.9, ls=":", label="pulse_start")
    if series["strong_drop_time"] is not None:
        ax.axvline(float(series["strong_drop_time"]), color="0.35", lw=0.9, ls="--", label="Mz<=0.9")
    if series["zero_cross_time"] is not None:
        ax.axvline(float(series["zero_cross_time"]), color="0.55", lw=0.9, ls="-.", label="Mz=0")
    ax.axvline(float(series["peak_excess_time"]), color="crimson", lw=0.9, ls="--", label="peak excess")


def make_trace_figure(case: dict, series: dict[str, object]) -> None:
    t = np.asarray(series["t"])
    pulse_start = float(series["pulse_start"])
    zoom_end = min(float(t[-1]), pulse_start + float(case["pulse_width"]) + 600.0)
    zoom_mask = (t >= pulse_start - 200.0) & (t <= zoom_end)
    colors = {"x": "tab:blue", "y": "tab:orange", "z": "tab:green"}

    fig, axes = plt.subplots(2, 2, figsize=(14, 8), sharex="col")
    axes[0, 0].plot(t, np.asarray(series["pulse_avg_mz"]), lw=1.8, label="pulse")
    axes[0, 0].plot(t, np.asarray(series["control_avg_mz"]), lw=1.4, label="no-pulse")
    axes[0, 0].plot(t, np.asarray(series["excess_avg_mz"]), lw=1.2, ls="--", label="pulse - no-pulse")
    axes[0, 0].axhline(0.0, color="0.75", lw=0.8, ls="--")
    add_marker_lines(axes[0, 0], series)
    axes[0, 0].set_title(f"{case['label']} avg_Mz full")
    axes[0, 0].set_ylabel("avg_Mz")
    axes[0, 0].legend(frameon=False, fontsize=8, ncol=2)

    axes[0, 1].plot(t[zoom_mask], np.asarray(series["pulse_avg_mz"])[zoom_mask], lw=1.8, label="pulse")
    axes[0, 1].plot(t[zoom_mask], np.asarray(series["control_avg_mz"])[zoom_mask], lw=1.4, label="no-pulse")
    axes[0, 1].plot(t[zoom_mask], np.asarray(series["excess_avg_mz"])[zoom_mask], lw=1.2, ls="--", label="pulse - no-pulse")
    axes[0, 1].axhline(0.0, color="0.75", lw=0.8, ls="--")
    add_marker_lines(axes[0, 1], series)
    axes[0, 1].set_title("avg_Mz zoom")

    for component in COMPONENTS:
        color = colors[component]
        axes[1, 0].plot(t, np.asarray(series["pulse_component_absmax"][component]), color=color, lw=1.6, label=f"pulse {component}")
        axes[1, 0].plot(t, np.asarray(series["control_component_absmax"][component]), color=color, lw=1.1, alpha=0.55, label=f"no-pulse {component}")
        axes[1, 0].plot(t, np.asarray(series["excess_component_absmax"][component]), color=color, lw=1.1, ls="--", label=f"excess {component}")
        axes[1, 1].plot(t[zoom_mask], np.asarray(series["pulse_component_absmax"][component])[zoom_mask], color=color, lw=1.6, label=f"pulse {component}")
        axes[1, 1].plot(t[zoom_mask], np.asarray(series["control_component_absmax"][component])[zoom_mask], color=color, lw=1.1, alpha=0.55, label=f"no-pulse {component}")
        axes[1, 1].plot(t[zoom_mask], np.asarray(series["excess_component_absmax"][component])[zoom_mask], color=color, lw=1.1, ls="--", label=f"excess {component}")

    add_marker_lines(axes[1, 0], series)
    add_marker_lines(axes[1, 1], series)
    axes[1, 0].set_title("max_site |sigma_delta_component| full")
    axes[1, 1].set_title("max_site |sigma_delta_component| zoom")
    axes[1, 0].set_ylabel("magnitude")
    axes[1, 0].set_xlabel("t")
    axes[1, 1].set_xlabel("t")
    axes[1, 0].legend(frameon=False, fontsize=7, ncol=3)
    fig.tight_layout()
    fig.savefig(ANALYSIS_DIR / f"{case['label']}_mechanism_traces.png", dpi=180)
    plt.close(fig)


def make_heatmap_figure(case: dict, series: dict[str, object]) -> None:
    t = np.asarray(series["t"])
    mask = pulse_window_mask(t, float(series["pulse_start"]), float(series["pulse_width"]))
    ordered_sites = list(series["ordered_sites"])
    fig, axes = plt.subplots(3, 1, figsize=(13, 8), sharex=True)

    max_abs = 0.0
    for component in COMPONENTS:
        max_abs = max(max_abs, float(np.max(np.abs(np.asarray(series["excess_cubes"][component])[:, mask]))))
    vmax = max_abs if max_abs > 0.0 else 1.0

    for ax, component in zip(axes, COMPONENTS, strict=True):
        matrix = np.asarray(series["excess_cubes"][component])[:, mask]
        image = ax.imshow(
            matrix,
            aspect="auto",
            origin="lower",
            cmap="coolwarm",
            vmin=-vmax,
            vmax=vmax,
            extent=[float(t[mask][0]), float(t[mask][-1]), ordered_sites[0] - 0.5, ordered_sites[-1] + 0.5],
        )
        ax.axvline(float(series["pulse_start"]), color="k", lw=0.9, ls=":")
        ax.axvline(float(series["peak_excess_time"]), color="crimson", lw=0.9, ls="--")
        if series["strong_drop_time"] is not None:
            ax.axvline(float(series["strong_drop_time"]), color="0.25", lw=0.9, ls="--")
        ax.set_ylabel(f"site ({component})")
        ax.set_yticks(ordered_sites)
        fig.colorbar(image, ax=ax, pad=0.01, label="pulse - no-pulse")

    axes[-1].set_xlabel("t")
    fig.suptitle(f"{case['label']} excess sigma_delta by site/component")
    fig.tight_layout()
    fig.savefig(ANALYSIS_DIR / f"{case['label']}_excess_heatmaps.png", dpi=180)
    plt.close(fig)


def write_review(rows: list[dict[str, object]]) -> None:
    lines = [
        "# Goal 6 Phase 2 Mechanism Review",
        "",
        f"- Heuristic marker: `first_strong_drop_time_avg_Mz_le_0p9` is the first time the pulse twin reaches `avg_Mz <= {DROP_THRESHOLD}`.",
        "- Excess signal means site/component-resolved `sigma_delta_pulse - sigma_delta_no_pulse`.",
        "",
        "## Clean-family pulse-specific readout",
    ]
    for row in rows:
        zero_text = "none" if row["first_zero_cross_time_avg_Mz"] is None else f"{row['first_zero_cross_time_avg_Mz']:.1f}"
        drop_text = "none" if row["first_strong_drop_time_avg_Mz_le_0p9"] is None else f"{row['first_strong_drop_time_avg_Mz_le_0p9']:.1f}"
        lines.append(
            "- "
            f"{row['label']}: excess peak `{row['excess_sigma_delta_absmax_peak']:.6f}` at `t={row['excess_sigma_delta_absmax_peak_time']:.1f}` "
            f"(component `{row['peak_excess_component']}`), `pulse_start={row['pulse_start']:.1f}`, "
            f"`first strong drop={drop_text}`, `first zero crossing={zero_text}`, "
            f"`pulse avg_Mz_final={row['pulse_avg_Mz_final']:.6f}`, `control avg_Mz_final={row['control_avg_Mz_final']:.6f}`."
        )

    sorted_rows = sorted(rows, key=lambda row: float(row["excess_sigma_delta_absmax_peak"]), reverse=True)
    strongest = sorted_rows[0]
    earliest = min(rows, key=lambda row: math.inf if row["first_strong_drop_time_avg_Mz_le_0p9"] is None else float(row["first_strong_drop_time_avg_Mz_le_0p9"]))
    lines.extend(
        [
            "",
            "## Short interpretation",
            "- The raw startup `max|sigma_delta|` peak remains nearly identical in pulse and no-pulse twins, so mechanism attribution still requires subtraction rather than raw amplitude reading.",
            f"- The strongest clean-family excess occurs in `{strongest['label']}` with peak `{strongest['excess_sigma_delta_absmax_peak']:.6f}` at `t={strongest['excess_sigma_delta_absmax_peak_time']:.1f}`.",
            f"- The earliest strong magnetization drop among these clean families is `{earliest['label']}` at `t={earliest['first_strong_drop_time_avg_Mz_le_0p9']:.1f}` using the `avg_Mz<=0.9` marker.",
        ]
    )
    REVIEW_MD.write_text("\n".join(lines) + "\n")


def main() -> None:
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    cases = json.loads(MANIFEST_PATH.read_text())
    write_case_manifest(cases)
    write_scan_config(cases)

    rows: list[dict[str, object]] = []
    for case in cases:
        pulse_data = read_tsv(Path(case["pulse_results_dir"]) / "observables.tsv")
        control_data = read_tsv(Path(case["control_results_dir"]) / "observables.tsv")
        row, series = summarise_case(case, pulse_data, control_data)
        rows.append(row)
        make_trace_figure(case, series)
        make_heatmap_figure(case, series)

    write_summary(rows)
    write_review(rows)


if __name__ == "__main__":
    main()
