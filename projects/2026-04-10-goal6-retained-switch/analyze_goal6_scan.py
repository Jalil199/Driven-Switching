"""Post-process Goal 6 retained-switch follow-up results."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parent
MODE_NEQ = "nonequilibrium_only"
MODE_FULL = "full_spin_density"


@dataclass(frozen=True)
class Criteria:
    clean_pre_threshold: float = 0.0
    strong_pulse_min: float = -0.8
    large_pulse_drop: float = 1.2
    preferred_final: float = -0.5
    strong_final: float = -0.8


CRITERIA = Criteria()


def latest_scan_root() -> Path:
    roots = sorted((PROJECT_DIR / "results").glob("queued_goal6_*"))
    if not roots:
        raise FileNotFoundError("No queued_goal6_* result roots found")
    return roots[-1]


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def read_summary(scan_root: Path) -> pd.DataFrame:
    df = pd.read_csv(scan_root / "scan_summary.tsv", sep="\t")
    numeric_cols = [
        "pulse_amp", "gamma_so", "j_sd", "j_ex", "j_ani", "damping", "Bx",
        "left_bias_scale", "right_bias_scale", "pulse_start", "pulse_ramp",
        "pulse_width", "t_end", "t_final", "avg_Mz_final", "avg_Mz_min",
        "avg_Mz_max", "current_balance_peak_abs", "current_peak_abs",
        "max_site_charge_deviation_from_1", "max_classical_spin_norm_deviation",
        "max_sigma_noneq_abs", "max_sigma_eq_abs", "max_sigma_delta_abs",
    ]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def load_trace(row: pd.Series, usecols: list[str] | None = None) -> pd.DataFrame:
    return pd.read_csv(Path(str(row["results_dir"])) / "observables.tsv", sep="\t", usecols=usecols)


def first_crossing_time(times: np.ndarray, values: np.ndarray, threshold: float) -> float:
    idx = np.flatnonzero(values <= threshold)
    return float(times[idx[0]]) if idx.size else np.nan


def pulse_value(times: np.ndarray, row: pd.Series) -> np.ndarray:
    start = float(row["pulse_start"])
    ramp = float(row["pulse_ramp"])
    width = float(row["pulse_width"])
    amp = float(row["bias_amplitude"])
    tau = times - start
    out = np.zeros_like(times, dtype=float)
    rise = (tau >= 0.0) & (tau < ramp)
    plateau = (tau >= ramp) & (tau < width)
    fall = (tau >= width) & (tau < width + ramp)
    if ramp > 0.0:
        out[rise] = np.sin((math.pi / 2.0) * tau[rise] / ramp) ** 2 * amp
        out[fall] = np.cos((math.pi / 2.0) * (tau[fall] - width) / ramp) ** 2 * amp
    out[plateau] = amp
    return out


def empty_metrics() -> dict[str, object]:
    return {
        "trace_present": False,
        "pre_pulse_min_avg_Mz": np.nan,
        "pre_pulse_cross_zero": False,
        "clean_pre_pulse": False,
        "avg_Mz_at_pulse_start": np.nan,
        "avg_Mz_at_pulse_end": np.nan,
        "pulse_window_min_avg_Mz": np.nan,
        "pulse_window_time_min_avg_Mz": np.nan,
        "pulse_window_drop_from_start": np.nan,
        "pulse_window_rebound_to_end": np.nan,
        "post_pulse_min_avg_Mz": np.nan,
        "post_pulse_time_min_avg_Mz": np.nan,
        "final_rebound_from_global_min": np.nan,
        "first_cross_zero_time": np.nan,
        "first_cross_m0p8_time": np.nan,
        "cross_zero_during_pulse": False,
        "cross_m0p8_during_pulse": False,
        "preferred_final_retention": False,
        "strong_final_retention": False,
        "fallback_delayed_strong": False,
    }


def trajectory_metrics(row: pd.Series) -> dict[str, object]:
    obs_path = Path(str(row["results_dir"])) / "observables.tsv"
    if not obs_path.is_file():
        return empty_metrics()

    obs = pd.read_csv(obs_path, sep="\t", usecols=["t", "avg_Mz"])
    times = obs["t"].to_numpy(float)
    mz = obs["avg_Mz"].to_numpy(float)
    pulse_start = float(row["pulse_start"])
    pulse_end = pulse_start + float(row["pulse_width"]) + float(row["pulse_ramp"])
    pre = times < pulse_start
    pulse = (times >= pulse_start) & (times <= pulse_end)
    post = times >= pulse_end

    pre_vals = mz[pre]
    pulse_vals = mz[pulse]
    pulse_times = times[pulse]
    post_vals = mz[post]
    post_times = times[post]

    pre_min = float(np.nanmin(pre_vals)) if pre_vals.size else np.nan
    pre_cross_zero = bool(pre_vals.size and np.nanmin(pre_vals) <= 0.0)
    start_mz = float(pre_vals[-1]) if pre_vals.size else np.nan
    end_candidates = np.flatnonzero(times <= pulse_end)
    end_mz = float(mz[end_candidates[-1]]) if end_candidates.size else np.nan
    if pulse_vals.size:
        pulse_min_idx = int(np.nanargmin(pulse_vals))
        pulse_min = float(pulse_vals[pulse_min_idx])
        pulse_min_time = float(pulse_times[pulse_min_idx])
    else:
        pulse_min = np.nan
        pulse_min_time = np.nan
    if post_vals.size:
        post_min_idx = int(np.nanargmin(post_vals))
        post_min = float(post_vals[post_min_idx])
        post_min_time = float(post_times[post_min_idx])
    else:
        post_min = np.nan
        post_min_time = np.nan

    first_zero = first_crossing_time(times, mz, 0.0)
    first_m0p8 = first_crossing_time(times, mz, CRITERIA.strong_pulse_min)
    final_mz = float(row["avg_Mz_final"])
    global_min = float(np.nanmin(mz)) if mz.size else np.nan
    clean_pre = bool(np.isfinite(pre_min) and pre_min > CRITERIA.clean_pre_threshold and not pre_cross_zero)
    pulse_drop = start_mz - pulse_min if np.isfinite(start_mz) and np.isfinite(pulse_min) else np.nan
    pulse_rebound = end_mz - pulse_min if np.isfinite(end_mz) and np.isfinite(pulse_min) else np.nan

    return {
        "trace_present": True,
        "pre_pulse_min_avg_Mz": pre_min,
        "pre_pulse_cross_zero": pre_cross_zero,
        "clean_pre_pulse": clean_pre,
        "avg_Mz_at_pulse_start": start_mz,
        "avg_Mz_at_pulse_end": end_mz,
        "pulse_window_min_avg_Mz": pulse_min,
        "pulse_window_time_min_avg_Mz": pulse_min_time,
        "pulse_window_drop_from_start": pulse_drop,
        "pulse_window_rebound_to_end": pulse_rebound,
        "post_pulse_min_avg_Mz": post_min,
        "post_pulse_time_min_avg_Mz": post_min_time,
        "final_rebound_from_global_min": final_mz - global_min if np.isfinite(global_min) else np.nan,
        "first_cross_zero_time": first_zero,
        "first_cross_m0p8_time": first_m0p8,
        "cross_zero_during_pulse": bool(np.isfinite(first_zero) and pulse_start <= first_zero <= pulse_end),
        "cross_m0p8_during_pulse": bool(np.isfinite(first_m0p8) and pulse_start <= first_m0p8 <= pulse_end),
        "preferred_final_retention": bool(final_mz <= CRITERIA.preferred_final),
        "strong_final_retention": bool(final_mz <= CRITERIA.strong_final),
        "fallback_delayed_strong": bool(post_min <= CRITERIA.strong_pulse_min if np.isfinite(post_min) else False),
    }


def add_metrics(df: pd.DataFrame) -> pd.DataFrame:
    metrics = pd.DataFrame([trajectory_metrics(row) for _, row in df.iterrows()])
    out = pd.concat([df.reset_index(drop=True), metrics.reset_index(drop=True)], axis=1)
    out["large_pulse_drop"] = out["pulse_window_drop_from_start"] >= CRITERIA.large_pulse_drop
    out["strong_in_pulse"] = out["pulse_window_min_avg_Mz"] <= CRITERIA.strong_pulse_min
    out["in_pulse_rank_score"] = (
        5.0 * out["clean_pre_pulse"].astype(int)
        + 4.0 * out["strong_in_pulse"].astype(int)
        + 3.0 * out["large_pulse_drop"].astype(int)
        - out["pulse_window_min_avg_Mz"].fillna(1.0)
        + out["pulse_window_drop_from_start"].fillna(-2.0)
        - 0.25 * out["final_rebound_from_global_min"].clip(lower=0).fillna(4.0)
        + 0.5 * out["preferred_final_retention"].astype(int)
    )
    out["retention_rank_score"] = (
        5.0 * out["clean_pre_pulse"].astype(int)
        + 4.0 * out["strong_final_retention"].astype(int)
        + 3.0 * out["preferred_final_retention"].astype(int)
        + 1.5 * out["fallback_delayed_strong"].astype(int)
        - out["avg_Mz_final"].fillna(1.0)
        - 0.20 * out["final_rebound_from_global_min"].clip(lower=0).fillna(4.0)
        - 0.0001 * out["first_cross_zero_time"].fillna(1e7)
    )
    return out


def common_cols() -> list[str]:
    return [
        "run_id", "branch", "Hs_update_mode",
        "pulse_amp", "gamma_so", "j_sd", "j_ex", "j_ani", "damping",
        "pulse_width", "t_end",
        "pre_pulse_min_avg_Mz", "pre_pulse_cross_zero", "clean_pre_pulse",
        "avg_Mz_at_pulse_start", "pulse_window_min_avg_Mz",
        "pulse_window_time_min_avg_Mz", "pulse_window_drop_from_start",
        "pulse_window_rebound_to_end", "post_pulse_min_avg_Mz",
        "post_pulse_time_min_avg_Mz", "avg_Mz_final",
        "final_rebound_from_global_min", "first_cross_zero_time",
        "first_cross_m0p8_time", "large_pulse_drop", "strong_in_pulse",
        "preferred_final_retention", "strong_final_retention",
        "fallback_delayed_strong", "max_sigma_delta_abs", "max_sigma_noneq_abs",
        "max_sigma_eq_abs", "current_balance_peak_abs", "current_peak_abs",
        "max_site_charge_deviation_from_1", "results_dir",
    ]


def save_tables(df: pd.DataFrame, analysis_dir: Path) -> None:
    cols = [col for col in common_cols() if col in df.columns]
    in_pulse = df.sort_values(
        ["clean_pre_pulse", "strong_in_pulse", "large_pulse_drop", "in_pulse_rank_score"],
        ascending=[False, False, False, False],
    )
    in_pulse[cols].to_csv(analysis_dir / "goal6_in_pulse_ranked_runs.tsv", sep="\t", index=False)

    retention = df.sort_values(
        ["clean_pre_pulse", "strong_final_retention", "preferred_final_retention", "fallback_delayed_strong", "retention_rank_score"],
        ascending=[False, False, False, False, False],
    )
    retention[cols].to_csv(analysis_dir / "goal6_retention_ranked_runs.tsv", sep="\t", index=False)
    retention[retention["clean_pre_pulse"].astype(bool)].head(16)[cols].to_csv(
        analysis_dir / "goal6_best_retention_candidates.tsv", sep="\t", index=False
    )

    branch_summary = (
        df.groupby("branch", dropna=False)
        .agg(
            runs=("run_id", "count"),
            done=("trace_present", "sum"),
            clean_pre=("clean_pre_pulse", "sum"),
            strong_in_pulse=("strong_in_pulse", "sum"),
            large_pulse_drop=("large_pulse_drop", "sum"),
            preferred_final_retention=("preferred_final_retention", "sum"),
            strong_final_retention=("strong_final_retention", "sum"),
            min_pulse_window_avg_Mz=("pulse_window_min_avg_Mz", "min"),
            min_post_pulse_avg_Mz=("post_pulse_min_avg_Mz", "min"),
            min_avg_Mz_final=("avg_Mz_final", "min"),
            max_sigma_delta_abs=("max_sigma_delta_abs", "max"),
            max_current_balance_peak_abs=("current_balance_peak_abs", "max"),
        )
        .reset_index()
    )
    branch_summary.to_csv(analysis_dir / "goal6_branch_summary.tsv", sep="\t", index=False)

    cross = build_hs_mode_crosscheck(df)
    cross.to_csv(analysis_dir / "goal6_hs_mode_crosscheck.tsv", sep="\t", index=False)


def build_hs_mode_crosscheck(df: pd.DataFrame) -> pd.DataFrame:
    center = df[df["branch"].eq("C_center_hs_mode_crosscheck")].copy()
    if center.empty:
        return pd.DataFrame()
    key_cols = ["pulse_amp", "gamma_so", "j_sd", "j_ex", "j_ani", "damping", "pulse_width"]
    rows = []
    for key, group in center.groupby(key_cols, dropna=False):
        if set(group["Hs_update_mode"]) != {MODE_NEQ, MODE_FULL}:
            continue
        neq = group[group["Hs_update_mode"].eq(MODE_NEQ)].iloc[0]
        full = group[group["Hs_update_mode"].eq(MODE_FULL)].iloc[0]
        row = {col: val for col, val in zip(key_cols, key)}
        row.update({
            "run_id_neq": neq["run_id"],
            "run_id_full": full["run_id"],
            "pulse_window_min_neq": neq["pulse_window_min_avg_Mz"],
            "pulse_window_min_full": full["pulse_window_min_avg_Mz"],
            "delta_pulse_window_min_full_minus_neq": full["pulse_window_min_avg_Mz"] - neq["pulse_window_min_avg_Mz"],
            "post_pulse_min_neq": neq["post_pulse_min_avg_Mz"],
            "post_pulse_min_full": full["post_pulse_min_avg_Mz"],
            "delta_post_pulse_min_full_minus_neq": full["post_pulse_min_avg_Mz"] - neq["post_pulse_min_avg_Mz"],
            "avg_Mz_final_neq": neq["avg_Mz_final"],
            "avg_Mz_final_full": full["avg_Mz_final"],
            "delta_avg_Mz_final_full_minus_neq": full["avg_Mz_final"] - neq["avg_Mz_final"],
            "max_sigma_delta_abs_neq": neq["max_sigma_delta_abs"],
            "max_sigma_delta_abs_full": full["max_sigma_delta_abs"],
            "delta_sigma_delta_full_minus_neq": full["max_sigma_delta_abs"] - neq["max_sigma_delta_abs"],
            "current_balance_peak_abs_neq": neq["current_balance_peak_abs"],
            "current_balance_peak_abs_full": full["current_balance_peak_abs"],
            "delta_current_balance_full_minus_neq": full["current_balance_peak_abs"] - neq["current_balance_peak_abs"],
        })
        rows.append(row)
    return pd.DataFrame(rows)


def plot_trajectories(df: pd.DataFrame, path: Path, title: str, rows: pd.DataFrame) -> None:
    rows = rows.head(8)
    if rows.empty:
        return
    fig, axes = plt.subplots(4, 2, figsize=(14, 13), sharex=False, sharey=True, constrained_layout=True)
    axes_flat = axes.ravel()
    for ax, (_, row) in zip(axes_flat, rows.iterrows()):
        obs = load_trace(row, ["t", "avg_Mz", "neel_z"])
        t = obs["t"].to_numpy(float)
        pulse = pulse_value(t, row)
        pulse_scaled = (pulse / max(float(row["bias_amplitude"]), 1e-12)) * 2.0 - 1.0
        pulse_start = float(row["pulse_start"])
        pulse_end = pulse_start + float(row["pulse_width"]) + float(row["pulse_ramp"])
        ax.plot(t, obs["avg_Mz"], color="#174a7c", lw=1.5, label="avg_Mz")
        ax.plot(t, obs["neel_z"], color="#6c757d", lw=0.9, ls="--", alpha=0.7, label="neel_z")
        ax.plot(t, pulse_scaled, color="#b23a48", lw=0.8, alpha=0.5, label="pulse")
        ax.axvspan(pulse_start, pulse_end, color="#f4a261", alpha=0.12)
        ax.axhline(0.0, color="black", lw=0.7, alpha=0.35)
        ax.axhline(-0.8, color="#b23a48", lw=0.7, ls=":", alpha=0.7)
        ax.set_ylim(-1.05, 1.05)
        ax.grid(alpha=0.2)
        ax.set_title(
            f"{row['run_id']}\n"
            f"branch={row['branch']}, min_pulse={row['pulse_window_min_avg_Mz']:.3f}, "
            f"post_min={row['post_pulse_min_avg_Mz']:.3f}, final={row['avg_Mz_final']:.3f}"
        )
    for ax in axes_flat[len(rows):]:
        ax.axis("off")
    axes_flat[0].legend(loc="lower left", ncol=3)
    fig.suptitle(title, y=1.01)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def save_figures(df: pd.DataFrame, analysis_dir: Path) -> None:
    plot_trajectories(
        df,
        analysis_dir / "goal6_top_branch_a_trajectories.png",
        "Goal 6 top Branch A in-pulse trajectories",
        df[df["branch"].isin(["A_main_hsneq_stabilization", "A_high_damping_strip"])].sort_values(
            ["clean_pre_pulse", "strong_in_pulse", "large_pulse_drop", "in_pulse_rank_score"],
            ascending=[False, False, False, False],
        ),
    )
    plot_trajectories(
        df,
        analysis_dir / "goal6_top_branch_b_trajectories.png",
        "Goal 6 top Branch B retention trajectories",
        df[df["branch"].eq("B_main_hsfull_retention")].sort_values(
            ["clean_pre_pulse", "strong_final_retention", "preferred_final_retention", "fallback_delayed_strong", "retention_rank_score"],
            ascending=[False, False, False, False, False],
        ),
    )
    plot_trajectories(
        df,
        analysis_dir / "goal6_hs_mode_crosscheck_trajectories.png",
        "Goal 6 Hs update mode center crosschecks",
        df[df["branch"].eq("C_center_hs_mode_crosscheck")].sort_values(["pulse_amp", "Hs_update_mode"]),
    )

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    done = df[df["trace_present"].astype(bool)]
    for ax, value, title in [
        (axes[0], "pulse_window_min_avg_Mz", "Pulse-window minimum"),
        (axes[1], "post_pulse_min_avg_Mz", "Post-pulse minimum"),
    ]:
        sc = ax.scatter(
            done["pulse_width"],
            done["damping"],
            c=done[value],
            s=70,
            cmap="coolwarm_r",
            edgecolor="black",
            linewidth=0.35,
        )
        ax.set_xlabel("pulse_width")
        ax.set_ylabel("damping")
        ax.set_title(title)
        fig.colorbar(sc, ax=ax, label=value)
    fig.savefig(analysis_dir / "goal6_minimum_maps.png", bbox_inches="tight")
    plt.close(fig)


def write_summary(scan_root: Path, analysis_dir: Path, df: pd.DataFrame) -> None:
    done = df[df["trace_present"].astype(bool)]
    best_in = done.sort_values(
        ["clean_pre_pulse", "strong_in_pulse", "large_pulse_drop", "in_pulse_rank_score"],
        ascending=[False, False, False, False],
    ).head(1)
    best_ret = done.sort_values(
        ["clean_pre_pulse", "strong_final_retention", "preferred_final_retention", "fallback_delayed_strong", "retention_rank_score"],
        ascending=[False, False, False, False, False],
    ).head(1)
    lines = [
        "# Goal 6 Analysis Summary",
        "",
        f"- scan root: `{scan_root}`",
        f"- completed traces analyzed: {len(done)} / {len(df)}",
        f"- clean pre-pulse traces: {int(done['clean_pre_pulse'].sum()) if not done.empty else 0}",
        f"- strong in-pulse traces (`avg_Mz <= -0.8`): {int(done['strong_in_pulse'].sum()) if not done.empty else 0}",
        f"- preferred retained final traces (`avg_Mz_final <= -0.5`): {int(done['preferred_final_retention'].sum()) if not done.empty else 0}",
        f"- strong retained final traces (`avg_Mz_final <= -0.8`): {int(done['strong_final_retention'].sum()) if not done.empty else 0}",
    ]
    if not best_in.empty:
        row = best_in.iloc[0]
        lines.append(
            f"- best in-pulse candidate: `{row['run_id']}` with pulse min {row['pulse_window_min_avg_Mz']:.6g}, "
            f"drop {row['pulse_window_drop_from_start']:.6g}, final {row['avg_Mz_final']:.6g}"
        )
    if not best_ret.empty:
        row = best_ret.iloc[0]
        lines.append(
            f"- best retention candidate: `{row['run_id']}` with final {row['avg_Mz_final']:.6g}, "
            f"post-pulse min {row['post_pulse_min_avg_Mz']:.6g}, pulse min {row['pulse_window_min_avg_Mz']:.6g}"
        )
    (analysis_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str]) -> int:
    scan_root = Path(argv[1]).resolve() if len(argv) > 1 else latest_scan_root()
    analysis_dir = scan_root / "analysis"
    ensure_dir(analysis_dir)
    df = add_metrics(read_summary(scan_root))
    df.to_csv(analysis_dir / "goal6_all_runs_with_metrics.tsv", sep="\t", index=False)
    save_tables(df, analysis_dir)
    save_figures(df, analysis_dir)
    write_summary(scan_root, analysis_dir, df)
    print(analysis_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
