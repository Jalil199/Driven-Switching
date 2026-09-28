"""Post-process the completed Goal 5 queued switching scan.

The run-time ``full_reversal_detected`` flag is intentionally loose: it only
checks final sign change. This script adds trajectory-based criteria and writes
compact tables/figures for deciding the next scan.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parent
SCAN_ROOT = PROJECT_DIR / "results" / "queued_parameter_scan_20260410T001920Z_goal5"
SUMMARY_PATH = SCAN_ROOT / "scan_summary.tsv"
ANALYSIS_DIR = SCAN_ROOT / "analysis"

MODE_NEQ = "nonequilibrium_only"
MODE_FULL = "full_spin_density"


@dataclass(frozen=True)
class Criteria:
    final_sign_threshold: float = 0.0
    strong_final_threshold: float = -0.8
    tail_mean_threshold: float = -0.8
    tail_fraction_threshold: float = 0.8
    tail_fraction_window: float = 0.10
    clean_pre_pulse_min_threshold: float = 0.0
    large_pulse_drop_threshold: float = 0.8


CRITERIA = Criteria()


KEY_COLS = [
    "pulse_amp",
    "gamma_so",
    "j_sd",
    "j_ex",
    "j_ani",
    "damping",
    "Bx",
    "left_bias_scale",
    "right_bias_scale",
    "pulse_start",
    "pulse_ramp",
    "pulse_width",
    "t_end",
]


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def read_summary() -> pd.DataFrame:
    df = pd.read_csv(SUMMARY_PATH, sep="\t")
    for col in [
        "pulse_amp",
        "gamma_so",
        "j_sd",
        "j_ex",
        "j_ani",
        "damping",
        "Bx",
        "pulse_start",
        "pulse_ramp",
        "pulse_width",
        "t_end",
        "avg_Mz_final",
        "avg_Mz_min",
        "avg_Mz_max",
        "current_balance_peak_abs",
        "max_site_charge_deviation_from_1",
    ]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


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


def load_trace(row: pd.Series, usecols: list[str] | None = None) -> pd.DataFrame:
    obs_path = Path(str(row["results_dir"])) / "observables.tsv"
    return pd.read_csv(obs_path, sep="\t", usecols=usecols)


def first_crossing_time(times: np.ndarray, values: np.ndarray, threshold: float) -> float:
    idx = np.flatnonzero(values <= threshold)
    if idx.size == 0:
        return np.nan
    return float(times[idx[0]])


def trajectory_metrics(row: pd.Series) -> dict[str, float | bool | str]:
    obs_path = Path(str(row["results_dir"])) / "observables.tsv"
    if not obs_path.is_file():
        return {
            "trace_present": False,
            "tail_mean_avg_Mz": np.nan,
            "tail_min_avg_Mz": np.nan,
            "tail_max_avg_Mz": np.nan,
            "tail_fraction_below_m0p8": np.nan,
            "post_pulse_fraction_below_zero": np.nan,
            "post_pulse_fraction_below_m0p8": np.nan,
            "avg_Mz_initial": np.nan,
            "pre_pulse_min_avg_Mz": np.nan,
            "pre_pulse_max_avg_Mz": np.nan,
            "pre_pulse_drop_from_initial": np.nan,
            "pre_pulse_range_avg_Mz": np.nan,
            "clean_pre_pulse_positive": False,
            "avg_Mz_at_pulse_start": np.nan,
            "avg_Mz_at_pulse_end": np.nan,
            "post_pulse_min_avg_Mz": np.nan,
            "post_pulse_time_min_avg_Mz": np.nan,
            "pulse_to_final_delta_avg_Mz": np.nan,
            "pulse_window_min_avg_Mz": np.nan,
            "pulse_window_max_avg_Mz": np.nan,
            "pulse_window_mean_avg_Mz": np.nan,
            "pulse_window_time_min_avg_Mz": np.nan,
            "pulse_window_drop_from_start": np.nan,
            "pulse_window_start_to_end_delta": np.nan,
            "pulse_window_cross_zero": False,
            "pulse_window_cross_m0p8": False,
            "clean_pulse_window_sign_cross": False,
            "clean_pulse_window_strong_cross": False,
            "first_cross_zero_time": np.nan,
            "first_cross_m0p8_time": np.nan,
            "pre_pulse_cross_zero": False,
            "pre_pulse_cross_m0p8": False,
            "delayed_strong_cross_after_pulse": False,
            "clean_large_pulse_drop": False,
            "pulse_induced_final_sign_change": False,
            "pulse_induced_strong_final": False,
            "sustained_tail_reversal": False,
        }

    obs = pd.read_csv(obs_path, sep="\t", usecols=["t", "avg_Mz"])
    times = obs["t"].to_numpy(dtype=float)
    mz = obs["avg_Mz"].to_numpy(dtype=float)

    final_t = np.nanmax(times)
    tail_start = final_t - CRITERIA.tail_fraction_window * final_t
    tail_mask = times >= tail_start
    pulse_done = float(row["pulse_start"]) + float(row["pulse_width"]) + float(row["pulse_ramp"])
    pulse_mask = (times >= float(row["pulse_start"])) & (times <= pulse_done)
    post_mask = times >= pulse_done
    pre_pulse_mask = times < float(row["pulse_start"])
    pulse_start_idx = np.flatnonzero(pre_pulse_mask)
    pulse_end_idx = np.flatnonzero(times <= pulse_done)
    mz_at_pulse_start = float(mz[pulse_start_idx[-1]]) if pulse_start_idx.size else np.nan
    mz_at_pulse_end = float(mz[pulse_end_idx[-1]]) if pulse_end_idx.size else np.nan
    first_zero = first_crossing_time(times, mz, CRITERIA.final_sign_threshold)
    first_m0p8 = first_crossing_time(times, mz, CRITERIA.strong_final_threshold)
    final_mz = float(row["avg_Mz_final"])

    tail_vals = mz[tail_mask]
    post_vals = mz[post_mask]
    post_times = times[post_mask]
    pre_vals = mz[pre_pulse_mask]
    pulse_vals = mz[pulse_mask]
    pulse_times = times[pulse_mask]
    initial_mz = float(mz[0]) if mz.size else np.nan
    if pre_vals.size:
        pre_min = float(np.nanmin(pre_vals))
        pre_max = float(np.nanmax(pre_vals))
    else:
        pre_min = np.nan
        pre_max = np.nan
    pre_drop_from_initial = initial_mz - pre_min if np.isfinite(initial_mz) and np.isfinite(pre_min) else np.nan
    pre_range = pre_max - pre_min if np.isfinite(pre_max) and np.isfinite(pre_min) else np.nan
    tail_mean = float(np.nanmean(tail_vals)) if tail_vals.size else np.nan
    tail_frac = float(np.nanmean(tail_vals <= CRITERIA.strong_final_threshold)) if tail_vals.size else np.nan
    post_frac_zero = float(np.nanmean(post_vals <= CRITERIA.final_sign_threshold)) if post_vals.size else np.nan
    post_frac_strong = float(np.nanmean(post_vals <= CRITERIA.strong_final_threshold)) if post_vals.size else np.nan
    if post_vals.size:
        post_min_idx = int(np.nanargmin(post_vals))
        post_min = float(post_vals[post_min_idx])
        post_time_min = float(post_times[post_min_idx])
    else:
        post_min = np.nan
        post_time_min = np.nan
    if pulse_vals.size:
        pulse_min_idx = int(np.nanargmin(pulse_vals))
        pulse_min = float(pulse_vals[pulse_min_idx])
        pulse_time_min = float(pulse_times[pulse_min_idx])
        pulse_max = float(np.nanmax(pulse_vals))
        pulse_mean = float(np.nanmean(pulse_vals))
    else:
        pulse_min = np.nan
        pulse_time_min = np.nan
        pulse_max = np.nan
        pulse_mean = np.nan
    pulse_drop_from_start = mz_at_pulse_start - pulse_min if np.isfinite(mz_at_pulse_start) else np.nan
    pulse_start_to_end_delta = mz_at_pulse_end - mz_at_pulse_start if np.isfinite(mz_at_pulse_start) else np.nan
    pulse_cross_zero = bool(pulse_vals.size and np.nanmin(pulse_vals) <= CRITERIA.final_sign_threshold)
    pulse_cross_m0p8 = bool(pulse_vals.size and np.nanmin(pulse_vals) <= CRITERIA.strong_final_threshold)
    clean_pulse_sign_cross = bool(
        mz_at_pulse_start > CRITERIA.final_sign_threshold and pulse_cross_zero
    )
    clean_pulse_strong_cross = bool(
        mz_at_pulse_start > CRITERIA.final_sign_threshold and pulse_cross_m0p8
    )
    clean_pre_pulse_positive = bool(
        np.isfinite(pre_min)
        and pre_min > CRITERIA.clean_pre_pulse_min_threshold
        and mz_at_pulse_start > CRITERIA.final_sign_threshold
    )
    delayed_strong_cross_after_pulse = bool(
        np.isfinite(first_m0p8)
        and first_m0p8 > pulse_done
        and mz_at_pulse_start > CRITERIA.final_sign_threshold
        and not (np.isfinite(first_zero) and first_zero < float(row["pulse_start"]))
    )
    clean_large_pulse_drop = bool(
        clean_pre_pulse_positive
        and np.isfinite(pulse_drop_from_start)
        and pulse_drop_from_start >= CRITERIA.large_pulse_drop_threshold
    )

    return {
        "trace_present": True,
        "tail_mean_avg_Mz": tail_mean,
        "tail_min_avg_Mz": float(np.nanmin(tail_vals)) if tail_vals.size else np.nan,
        "tail_max_avg_Mz": float(np.nanmax(tail_vals)) if tail_vals.size else np.nan,
        "tail_fraction_below_m0p8": tail_frac,
        "post_pulse_fraction_below_zero": post_frac_zero,
        "post_pulse_fraction_below_m0p8": post_frac_strong,
        "avg_Mz_initial": initial_mz,
        "pre_pulse_min_avg_Mz": pre_min,
        "pre_pulse_max_avg_Mz": pre_max,
        "pre_pulse_drop_from_initial": pre_drop_from_initial,
        "pre_pulse_range_avg_Mz": pre_range,
        "clean_pre_pulse_positive": clean_pre_pulse_positive,
        "avg_Mz_at_pulse_start": mz_at_pulse_start,
        "avg_Mz_at_pulse_end": mz_at_pulse_end,
        "post_pulse_min_avg_Mz": post_min,
        "post_pulse_time_min_avg_Mz": post_time_min,
        "pulse_to_final_delta_avg_Mz": final_mz - mz_at_pulse_start,
        "pulse_window_min_avg_Mz": pulse_min,
        "pulse_window_max_avg_Mz": pulse_max,
        "pulse_window_mean_avg_Mz": pulse_mean,
        "pulse_window_time_min_avg_Mz": pulse_time_min,
        "pulse_window_drop_from_start": pulse_drop_from_start,
        "pulse_window_start_to_end_delta": pulse_start_to_end_delta,
        "pulse_window_cross_zero": pulse_cross_zero,
        "pulse_window_cross_m0p8": pulse_cross_m0p8,
        "clean_pulse_window_sign_cross": clean_pulse_sign_cross,
        "clean_pulse_window_strong_cross": clean_pulse_strong_cross,
        "first_cross_zero_time": first_zero,
        "first_cross_m0p8_time": first_m0p8,
        "pre_pulse_cross_zero": bool(np.isfinite(first_zero) and first_zero < float(row["pulse_start"])),
        "pre_pulse_cross_m0p8": bool(np.isfinite(first_m0p8) and first_m0p8 < float(row["pulse_start"])),
        "delayed_strong_cross_after_pulse": delayed_strong_cross_after_pulse,
        "clean_large_pulse_drop": clean_large_pulse_drop,
        "pulse_induced_final_sign_change": bool(
            mz_at_pulse_start > CRITERIA.final_sign_threshold
            and final_mz <= CRITERIA.final_sign_threshold
        ),
        "pulse_induced_strong_final": bool(
            mz_at_pulse_start > CRITERIA.final_sign_threshold
            and final_mz <= CRITERIA.strong_final_threshold
        ),
        "sustained_tail_reversal": bool(
            tail_mean <= CRITERIA.tail_mean_threshold
            and tail_frac >= CRITERIA.tail_fraction_threshold
        ),
    }


def add_analysis_columns(df: pd.DataFrame) -> pd.DataFrame:
    done = df["status"].eq("done") & df["observables_present"].astype(bool)
    metric_rows = []
    for _, row in df.iterrows():
        if bool(done.loc[row.name]):
            metric_rows.append(trajectory_metrics(row))
        else:
            metric_rows.append(trajectory_metrics(row))
    metrics = pd.DataFrame(metric_rows)
    out = pd.concat([df.reset_index(drop=True), metrics.reset_index(drop=True)], axis=1)
    out["final_sign_change"] = out["avg_Mz_final"] <= CRITERIA.final_sign_threshold
    out["strong_final_reversal"] = out["avg_Mz_final"] <= CRITERIA.strong_final_threshold
    out["deep_transient_reversal"] = out["avg_Mz_min"] <= CRITERIA.strong_final_threshold
    out["weak_final_near_zero"] = out["final_sign_change"] & ~out["strong_final_reversal"]
    out["run_rank_score"] = (
        -out["avg_Mz_final"].fillna(1.0)
        + 0.25 * -out["tail_mean_avg_Mz"].fillna(1.0)
        + 0.10 * -out["avg_Mz_min"].fillna(1.0)
    )
    out["pulse_window_rank_score"] = (
        out["pulse_window_drop_from_start"].fillna(-np.inf)
        + 0.25 * -out["pulse_window_min_avg_Mz"].fillna(1.0)
        + 0.10 * out["clean_pulse_window_sign_cross"].astype(int)
        + 0.20 * out["clean_pulse_window_strong_cross"].astype(int)
    )
    out["clean_pulse_response_rank_score"] = (
        out["pulse_window_drop_from_start"].fillna(-np.inf)
        + 0.50 * out["clean_pulse_window_strong_cross"].astype(int)
        + 0.35 * out["clean_pulse_window_sign_cross"].astype(int)
        + 0.25 * out["pulse_induced_final_sign_change"].astype(int)
        + 0.25 * out["delayed_strong_cross_after_pulse"].astype(int)
        + 0.15 * out["strong_final_reversal"].astype(int)
        - 0.30 * out["pre_pulse_drop_from_initial"].clip(lower=0).fillna(0.0)
    )
    return out


def pair_modes(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for key, group in df.groupby(KEY_COLS, dropna=False):
        if set(group["Hs_update_mode"]) != {MODE_NEQ, MODE_FULL}:
            continue
        neq = group[group["Hs_update_mode"].eq(MODE_NEQ)].iloc[0]
        full = group[group["Hs_update_mode"].eq(MODE_FULL)].iloc[0]
        row = {col: val for col, val in zip(KEY_COLS, key)}
        row.update(
            {
                "run_id_neq": neq["run_id"],
                "run_id_full": full["run_id"],
                "avg_Mz_final_neq": neq["avg_Mz_final"],
                "avg_Mz_final_full": full["avg_Mz_final"],
                "delta_final_full_minus_neq": full["avg_Mz_final"] - neq["avg_Mz_final"],
                "avg_Mz_min_neq": neq["avg_Mz_min"],
                "avg_Mz_min_full": full["avg_Mz_min"],
                "tail_mean_neq": neq["tail_mean_avg_Mz"],
                "tail_mean_full": full["tail_mean_avg_Mz"],
                "pulse_window_drop_neq": neq["pulse_window_drop_from_start"],
                "pulse_window_drop_full": full["pulse_window_drop_from_start"],
                "pulse_window_min_neq": neq["pulse_window_min_avg_Mz"],
                "pulse_window_min_full": full["pulse_window_min_avg_Mz"],
                "clean_large_pulse_drop_neq": bool(neq["clean_large_pulse_drop"]),
                "clean_large_pulse_drop_full": bool(full["clean_large_pulse_drop"]),
                "clean_pulse_window_sign_cross_neq": bool(neq["clean_pulse_window_sign_cross"]),
                "clean_pulse_window_sign_cross_full": bool(full["clean_pulse_window_sign_cross"]),
                "clean_pulse_window_strong_cross_neq": bool(neq["clean_pulse_window_strong_cross"]),
                "clean_pulse_window_strong_cross_full": bool(full["clean_pulse_window_strong_cross"]),
                "delayed_strong_cross_after_pulse_neq": bool(neq["delayed_strong_cross_after_pulse"]),
                "delayed_strong_cross_after_pulse_full": bool(full["delayed_strong_cross_after_pulse"]),
                "strong_final_neq": bool(neq["strong_final_reversal"]),
                "strong_final_full": bool(full["strong_final_reversal"]),
                "sustained_tail_neq": bool(neq["sustained_tail_reversal"]),
                "sustained_tail_full": bool(full["sustained_tail_reversal"]),
                "any_final_sign_change": bool(neq["final_sign_change"] or full["final_sign_change"]),
                "any_clean_pulse_window_sign_cross": bool(
                    neq["clean_pulse_window_sign_cross"] or full["clean_pulse_window_sign_cross"]
                ),
                "any_clean_pulse_window_strong_cross": bool(
                    neq["clean_pulse_window_strong_cross"] or full["clean_pulse_window_strong_cross"]
                ),
                "any_clean_large_pulse_drop": bool(neq["clean_large_pulse_drop"] or full["clean_large_pulse_drop"]),
                "any_delayed_strong_cross_after_pulse": bool(
                    neq["delayed_strong_cross_after_pulse"] or full["delayed_strong_cross_after_pulse"]
                ),
                "any_strong_final": bool(neq["strong_final_reversal"] or full["strong_final_reversal"]),
                "any_sustained_tail": bool(neq["sustained_tail_reversal"] or full["sustained_tail_reversal"]),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def save_tables(df: pd.DataFrame, paired: pd.DataFrame) -> None:
    cols = [
        "run_id",
        "Hs_update_mode",
        "pulse_amp",
        "gamma_so",
        "j_sd",
        "j_ex",
        "j_ani",
        "damping",
        "avg_Mz_final",
        "avg_Mz_min",
        "tail_mean_avg_Mz",
        "tail_fraction_below_m0p8",
        "post_pulse_fraction_below_zero",
        "post_pulse_fraction_below_m0p8",
        "avg_Mz_initial",
        "pre_pulse_min_avg_Mz",
        "pre_pulse_drop_from_initial",
        "pre_pulse_range_avg_Mz",
        "clean_pre_pulse_positive",
        "avg_Mz_at_pulse_start",
        "avg_Mz_at_pulse_end",
        "post_pulse_min_avg_Mz",
        "post_pulse_time_min_avg_Mz",
        "pulse_to_final_delta_avg_Mz",
        "pulse_window_min_avg_Mz",
        "pulse_window_mean_avg_Mz",
        "pulse_window_time_min_avg_Mz",
        "pulse_window_drop_from_start",
        "pulse_window_start_to_end_delta",
        "pulse_window_cross_zero",
        "pulse_window_cross_m0p8",
        "clean_pulse_window_sign_cross",
        "clean_pulse_window_strong_cross",
        "first_cross_zero_time",
        "first_cross_m0p8_time",
        "pre_pulse_cross_zero",
        "pre_pulse_cross_m0p8",
        "delayed_strong_cross_after_pulse",
        "clean_large_pulse_drop",
        "pulse_induced_final_sign_change",
        "pulse_induced_strong_final",
        "final_sign_change",
        "strong_final_reversal",
        "deep_transient_reversal",
        "sustained_tail_reversal",
        "clean_pulse_response_rank_score",
        "results_dir",
    ]
    ranked = df.sort_values(["strong_final_reversal", "avg_Mz_final", "tail_mean_avg_Mz"], ascending=[False, True, True])
    ranked[cols].to_csv(ANALYSIS_DIR / "goal5_ranked_runs.tsv", sep="\t", index=False)

    pulse_ranked = df.sort_values(
        [
            "clean_pre_pulse_positive",
            "clean_large_pulse_drop",
            "clean_pulse_window_strong_cross",
            "delayed_strong_cross_after_pulse",
            "pulse_induced_final_sign_change",
            "clean_pulse_window_sign_cross",
            "clean_pulse_response_rank_score",
            "post_pulse_min_avg_Mz",
            "pulse_window_min_avg_Mz",
        ],
        ascending=[False, False, False, False, False, False, False, True, True],
    )
    pulse_ranked[cols].to_csv(ANALYSIS_DIR / "goal5_pulse_window_ranked_runs.tsv", sep="\t", index=False)

    pair_cols = KEY_COLS + [
        "run_id_neq",
        "run_id_full",
        "avg_Mz_final_neq",
        "avg_Mz_final_full",
        "delta_final_full_minus_neq",
        "avg_Mz_min_neq",
        "avg_Mz_min_full",
        "tail_mean_neq",
        "tail_mean_full",
        "pulse_window_drop_neq",
        "pulse_window_drop_full",
        "pulse_window_min_neq",
        "pulse_window_min_full",
        "clean_large_pulse_drop_neq",
        "clean_large_pulse_drop_full",
        "clean_pulse_window_sign_cross_neq",
        "clean_pulse_window_sign_cross_full",
        "clean_pulse_window_strong_cross_neq",
        "clean_pulse_window_strong_cross_full",
        "delayed_strong_cross_after_pulse_neq",
        "delayed_strong_cross_after_pulse_full",
        "strong_final_neq",
        "strong_final_full",
        "sustained_tail_neq",
        "sustained_tail_full",
        "any_final_sign_change",
        "any_clean_pulse_window_sign_cross",
        "any_clean_pulse_window_strong_cross",
        "any_clean_large_pulse_drop",
        "any_delayed_strong_cross_after_pulse",
        "any_strong_final",
        "any_sustained_tail",
    ]
    paired.sort_values(["any_strong_final", "avg_Mz_final_full", "avg_Mz_final_neq"], ascending=[False, True, True])[
        pair_cols
    ].to_csv(ANALYSIS_DIR / "goal5_paired_hs_mode_comparison.tsv", sep="\t", index=False)

    by_jani = (
        df.groupby("j_ani", dropna=False)
        .agg(
            runs=("run_id", "count"),
            final_sign_changes=("final_sign_change", "sum"),
            clean_pulse_window_sign_crosses=("clean_pulse_window_sign_cross", "sum"),
            clean_pulse_window_strong_crosses=("clean_pulse_window_strong_cross", "sum"),
            strong_final_reversals=("strong_final_reversal", "sum"),
            sustained_tail_reversals=("sustained_tail_reversal", "sum"),
            min_avg_Mz_final=("avg_Mz_final", "min"),
            max_pulse_window_drop_from_start=("pulse_window_drop_from_start", "max"),
            min_pulse_window_avg_Mz=("pulse_window_min_avg_Mz", "min"),
            min_tail_mean_avg_Mz=("tail_mean_avg_Mz", "min"),
        )
        .reset_index()
    )
    by_jani.to_csv(ANALYSIS_DIR / "goal5_by_j_ani.tsv", sep="\t", index=False)

    selected = select_pulse_window_candidates(df)
    selected.to_csv(ANALYSIS_DIR / "goal5_selected_pulse_window_candidates.tsv", sep="\t", index=False)


def candidate_tier(row: pd.Series) -> str:
    if bool(row["clean_pulse_window_strong_cross"]):
        return "A_strong_clean_in_pulse"
    if bool(row["delayed_strong_cross_after_pulse"]):
        return "B_delayed_strong_after_pulse"
    if bool(row["clean_pulse_window_sign_cross"]) and bool(row["pulse_induced_final_sign_change"]):
        return "C_in_pulse_sign_and_negative_final"
    if bool(row["clean_pulse_window_sign_cross"]) and float(row["pulse_window_min_avg_Mz"]) <= -0.5:
        return "D_near_strong_clean_in_pulse"
    if bool(row["clean_pulse_window_sign_cross"]) and float(row["pulse_window_min_avg_Mz"]) <= -0.3:
        return "E_large_clean_in_pulse_drop"
    if bool(row["clean_large_pulse_drop"]):
        return "F_large_clean_drop_no_inversion"
    if bool(row["clean_pulse_window_sign_cross"]):
        return "G_clean_sign_cross_watchlist"
    return "H_not_selected"


def select_pulse_window_candidates(df: pd.DataFrame, limit: int = 16) -> pd.DataFrame:
    """Select clean pulse-response candidates for manual inspection.

    The main filter rejects trajectories whose magnetic order is already
    sign-changing before the pulse. Within that clean subset, ranking favors a
    large start-to-min drop during the pulse and then promotes any in-pulse,
    delayed strong, or final inversion evidence.
    """

    clean = df[
        df["clean_pre_pulse_positive"].astype(bool)
        & df["clean_large_pulse_drop"].astype(bool)
        & ~df["pre_pulse_cross_zero"].astype(bool)
    ].copy()
    clean = clean.sort_values(
        [
            "clean_pulse_window_strong_cross",
            "delayed_strong_cross_after_pulse",
            "pulse_induced_final_sign_change",
            "clean_pulse_window_sign_cross",
            "clean_pulse_response_rank_score",
            "post_pulse_min_avg_Mz",
            "pulse_window_min_avg_Mz",
        ],
        ascending=[False, False, False, False, False, True, True],
    )
    selected = clean.head(limit).copy()
    selected.insert(0, "candidate_rank", np.arange(1, len(selected) + 1))
    selected.insert(1, "candidate_tier", [candidate_tier(row) for _, row in selected.iterrows()])

    cols = [
        "candidate_rank",
        "candidate_tier",
        "run_id",
        "Hs_update_mode",
        "pulse_amp",
        "gamma_so",
        "j_sd",
        "j_ex",
        "j_ani",
        "damping",
        "Bx",
        "bias_amplitude",
        "left_bias_scale",
        "right_bias_scale",
        "pulse_start",
        "pulse_ramp",
        "pulse_width",
        "max_sigma_delta_abs",
        "avg_Mz_initial",
        "pre_pulse_min_avg_Mz",
        "pre_pulse_drop_from_initial",
        "pre_pulse_range_avg_Mz",
        "avg_Mz_at_pulse_start",
        "pulse_window_min_avg_Mz",
        "pulse_window_time_min_avg_Mz",
        "pulse_window_drop_from_start",
        "avg_Mz_at_pulse_end",
        "post_pulse_min_avg_Mz",
        "post_pulse_time_min_avg_Mz",
        "avg_Mz_final",
        "clean_large_pulse_drop",
        "clean_pulse_window_sign_cross",
        "clean_pulse_window_strong_cross",
        "delayed_strong_cross_after_pulse",
        "pulse_induced_final_sign_change",
        "strong_final_reversal",
        "clean_pulse_response_rank_score",
        "post_pulse_fraction_below_zero",
        "post_pulse_fraction_below_m0p8",
        "results_dir",
    ]
    return selected[cols]


def setup_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 13,
            "legend.fontsize": 9,
            "figure.dpi": 140,
            "savefig.dpi": 180,
        }
    )


def plot_top_trajectories(df: pd.DataFrame) -> None:
    top = df.sort_values(["avg_Mz_final", "tail_mean_avg_Mz", "avg_Mz_min"]).head(8)
    fig, axes = plt.subplots(4, 2, figsize=(14, 13), sharex=True, sharey=True, constrained_layout=True)
    axes_flat = axes.ravel()

    for ax, (_, row) in zip(axes_flat, top.iterrows()):
        obs = load_trace(row, ["t", "avg_Mz", "neel_z"])
        t = obs["t"].to_numpy(float)
        pulse = pulse_value(t, row)
        pulse_scaled = (pulse / max(float(row["bias_amplitude"]), 1e-12)) * 2.0 - 1.0

        ax.plot(t, obs["avg_Mz"], color="#174a7c", lw=1.8, label="avg_Mz")
        ax.plot(t, obs["neel_z"], color="#6c757d", lw=1.0, ls="--", alpha=0.75, label="neel_z")
        ax.plot(t, pulse_scaled, color="#b23a48", lw=0.9, alpha=0.55, label="pulse scaled")
        ax.axhline(0.0, color="black", lw=0.8, alpha=0.35)
        ax.axhline(-0.8, color="#b23a48", lw=0.8, ls=":", alpha=0.65)
        ax.set_ylim(-1.05, 1.05)
        ax.grid(alpha=0.2)
        ax.set_title(
            f"{row['run_id']}\n"
            f"pa={row['pulse_amp']}, gso={row['gamma_so']}, jsd={row['j_sd']}, "
            f"jex={row['j_ex']}, jani={row['j_ani']}, damp={row['damping']}, "
            f"{row['Hs_update_mode']}"
        )

    for ax in axes[-1, :]:
        ax.set_xlabel("time")
    for ax in axes[:, 0]:
        ax.set_ylabel("order parameter")
    axes_flat[0].legend(loc="lower left", ncol=3)
    fig.suptitle("Goal 5 strongest final avg_Mz trajectories", y=1.01)
    fig.savefig(ANALYSIS_DIR / "goal5_top_trajectories.png", bbox_inches="tight")
    plt.close(fig)


def plot_top_pulse_window_trajectories(df: pd.DataFrame) -> None:
    top = df.sort_values(
        [
            "clean_pre_pulse_positive",
            "clean_large_pulse_drop",
            "clean_pulse_window_strong_cross",
            "delayed_strong_cross_after_pulse",
            "pulse_induced_final_sign_change",
            "clean_pulse_window_sign_cross",
            "clean_pulse_response_rank_score",
            "post_pulse_min_avg_Mz",
            "pulse_window_min_avg_Mz",
        ],
        ascending=[False, False, False, False, False, False, False, True, True],
    ).head(8)
    fig, axes = plt.subplots(4, 2, figsize=(14, 13), sharex=True, sharey=True, constrained_layout=True)
    axes_flat = axes.ravel()

    for ax, (_, row) in zip(axes_flat, top.iterrows()):
        obs = load_trace(row, ["t", "avg_Mz", "neel_z"])
        t = obs["t"].to_numpy(float)
        pulse = pulse_value(t, row)
        pulse_scaled = (pulse / max(float(row["bias_amplitude"]), 1e-12)) * 2.0 - 1.0
        pulse_start = float(row["pulse_start"])
        pulse_end = pulse_start + float(row["pulse_width"]) + float(row["pulse_ramp"])

        ax.plot(t, obs["avg_Mz"], color="#174a7c", lw=1.8, label="avg_Mz")
        ax.plot(t, obs["neel_z"], color="#6c757d", lw=1.0, ls="--", alpha=0.75, label="neel_z")
        ax.plot(t, pulse_scaled, color="#b23a48", lw=0.9, alpha=0.55, label="pulse scaled")
        ax.axvspan(pulse_start, pulse_end, color="#f4a261", alpha=0.12)
        ax.axhline(0.0, color="black", lw=0.8, alpha=0.35)
        ax.axhline(-0.8, color="#b23a48", lw=0.8, ls=":", alpha=0.65)
        ax.scatter(
            [row["pulse_window_time_min_avg_Mz"]],
            [row["pulse_window_min_avg_Mz"]],
            color="#b23a48",
            s=28,
            zorder=4,
        )
        ax.set_ylim(-1.05, 1.05)
        ax.grid(alpha=0.2)
        ax.set_title(
            f"{row['run_id']}\n"
            f"drop={row['pulse_window_drop_from_start']:.3f}, "
            f"min_pulse={row['pulse_window_min_avg_Mz']:.3f}, "
            f"start={row['avg_Mz_at_pulse_start']:.3f}, final={row['avg_Mz_final']:.3f}, "
            f"{row['Hs_update_mode']}"
        )

    for ax in axes[-1, :]:
        ax.set_xlabel("time")
    for ax in axes[:, 0]:
        ax.set_ylabel("order parameter")
    axes_flat[0].legend(loc="lower left", ncol=3)
    fig.suptitle("Goal 5 clean large pulse responses", y=1.01)
    fig.savefig(ANALYSIS_DIR / "goal5_top_pulse_window_trajectories.png", bbox_inches="tight")
    plt.close(fig)


def plot_hs_mode_comparison(paired: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)

    sc = axes[0].scatter(
        paired["avg_Mz_final_neq"],
        paired["avg_Mz_final_full"],
        c=paired["j_ani"],
        s=60,
        cmap="viridis",
        edgecolor="black",
        linewidth=0.35,
    )
    axes[0].plot([-1.0, 1.0], [-1.0, 1.0], color="black", lw=1.0, ls="--", alpha=0.5)
    axes[0].axhline(-0.8, color="#b23a48", ls=":", lw=1.0)
    axes[0].axvline(-0.8, color="#b23a48", ls=":", lw=1.0)
    axes[0].set_xlabel("avg_Mz_final nonequilibrium_only")
    axes[0].set_ylabel("avg_Mz_final full_spin_density")
    axes[0].set_title("A/B final response")
    fig.colorbar(sc, ax=axes[0], label="j_ani")

    axes[1].hist(paired["delta_final_full_minus_neq"], bins=28, color="#2a9d8f", edgecolor="white")
    axes[1].axvline(0.0, color="black", lw=1.0)
    axes[1].set_xlabel("full - nonequilibrium avg_Mz_final")
    axes[1].set_ylabel("paired points")
    axes[1].set_title("Hs_update_mode shift")

    grouped = paired.groupby("j_ani").agg(
        any_sign=("any_final_sign_change", "sum"),
        any_strong=("any_strong_final", "sum"),
        any_sustained=("any_sustained_tail", "sum"),
    )
    x = np.arange(len(grouped.index))
    width = 0.25
    axes[2].bar(x - width, grouped["any_sign"], width, label="final sign")
    axes[2].bar(x, grouped["any_strong"], width, label="strong final")
    axes[2].bar(x + width, grouped["any_sustained"], width, label="sustained tail")
    axes[2].set_xticks(x)
    axes[2].set_xticklabels([str(v) for v in grouped.index])
    axes[2].set_xlabel("j_ani")
    axes[2].set_ylabel("paired points")
    axes[2].set_title("Reversal counts by anisotropy")
    axes[2].legend()
    axes[2].grid(axis="y", alpha=0.2)

    fig.savefig(ANALYSIS_DIR / "goal5_hs_mode_comparison.png", bbox_inches="tight")
    plt.close(fig)


def plot_parameter_maps(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)

    for ax, mode in zip(axes[0], [MODE_NEQ, MODE_FULL]):
        sub = df[df["Hs_update_mode"].eq(mode)]
        pivot = sub.pivot_table(index="j_ani", columns="pulse_amp", values="avg_Mz_final", aggfunc="min")
        im = ax.imshow(pivot.to_numpy(), aspect="auto", origin="lower", cmap="RdBu_r", vmin=-1.0, vmax=1.0)
        ax.set_title(f"Minimum final avg_Mz: {mode}")
        ax.set_xlabel("pulse_amp")
        ax.set_ylabel("j_ani")
        ax.set_xticks(np.arange(len(pivot.columns)))
        ax.set_xticklabels([str(v) for v in pivot.columns])
        ax.set_yticks(np.arange(len(pivot.index)))
        ax.set_yticklabels([str(v) for v in pivot.index])
        fig.colorbar(im, ax=ax, label="min avg_Mz_final")

    scatter = axes[1, 0].scatter(
        df["avg_Mz_min"],
        df["avg_Mz_final"],
        c=df["j_ani"],
        s=55,
        cmap="viridis",
        edgecolor="black",
        linewidth=0.25,
    )
    axes[1, 0].axhline(0.0, color="black", lw=0.8, alpha=0.4)
    axes[1, 0].axhline(-0.8, color="#b23a48", lw=0.9, ls=":")
    axes[1, 0].axvline(-0.8, color="#b23a48", lw=0.9, ls=":")
    axes[1, 0].set_xlabel("minimum avg_Mz during run")
    axes[1, 0].set_ylabel("final avg_Mz")
    axes[1, 0].set_title("Transient depth vs final state")
    axes[1, 0].grid(alpha=0.2)
    fig.colorbar(scatter, ax=axes[1, 0], label="j_ani")

    by = (
        df.groupby(["pulse_amp", "gamma_so", "Hs_update_mode"])
        .agg(best_final=("avg_Mz_final", "min"))
        .reset_index()
    )
    for mode, color in [(MODE_NEQ, "#174a7c"), (MODE_FULL, "#b23a48")]:
        sub = by[by["Hs_update_mode"].eq(mode)]
        labels = [f"pa={r.pulse_amp}\ngso={r.gamma_so}" for r in sub.itertuples()]
        axes[1, 1].scatter(labels, sub["best_final"], s=80, label=mode, color=color, alpha=0.8)
    axes[1, 1].axhline(0.0, color="black", lw=0.8, alpha=0.4)
    axes[1, 1].axhline(-0.8, color="#b23a48", lw=0.9, ls=":")
    axes[1, 1].set_ylabel("best final avg_Mz")
    axes[1, 1].set_title("Best final response by pulse/SOC")
    axes[1, 1].tick_params(axis="x", rotation=30)
    axes[1, 1].grid(axis="y", alpha=0.2)
    axes[1, 1].legend()

    fig.savefig(ANALYSIS_DIR / "goal5_parameter_maps.png", bbox_inches="tight")
    plt.close(fig)


def plot_pulse_window_maps(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)

    sc = axes[0].scatter(
        df["avg_Mz_at_pulse_start"],
        df["pulse_window_min_avg_Mz"],
        c=df["pulse_window_drop_from_start"],
        s=58,
        cmap="magma",
        edgecolor="black",
        linewidth=0.25,
    )
    axes[0].axhline(0.0, color="black", lw=0.8, alpha=0.4)
    axes[0].axhline(-0.8, color="#b23a48", lw=0.9, ls=":")
    axes[0].axvline(0.0, color="black", lw=0.8, alpha=0.4)
    axes[0].set_xlabel("avg_Mz at pulse_start")
    axes[0].set_ylabel("minimum avg_Mz during pulse")
    axes[0].set_title("Pulse-window depth")
    axes[0].grid(alpha=0.2)
    fig.colorbar(sc, ax=axes[0], label="drop from pulse_start")

    clean = df[df["avg_Mz_at_pulse_start"] > 0.0]
    pivot = clean.pivot_table(
        index="j_ani",
        columns="pulse_amp",
        values="pulse_window_drop_from_start",
        aggfunc="max",
    )
    im = axes[1].imshow(pivot.to_numpy(), aspect="auto", origin="lower", cmap="viridis")
    axes[1].set_title("Best clean pulse-window drop")
    axes[1].set_xlabel("pulse_amp")
    axes[1].set_ylabel("j_ani")
    axes[1].set_xticks(np.arange(len(pivot.columns)))
    axes[1].set_xticklabels([str(v) for v in pivot.columns])
    axes[1].set_yticks(np.arange(len(pivot.index)))
    axes[1].set_yticklabels([str(v) for v in pivot.index])
    fig.colorbar(im, ax=axes[1], label="max start-to-min drop")

    by = (
        clean.groupby(["pulse_amp", "gamma_so", "Hs_update_mode"])
        .agg(best_drop=("pulse_window_drop_from_start", "max"))
        .reset_index()
    )
    for mode, color in [(MODE_NEQ, "#174a7c"), (MODE_FULL, "#b23a48")]:
        sub = by[by["Hs_update_mode"].eq(mode)]
        labels = [f"pa={r.pulse_amp}\ngso={r.gamma_so}" for r in sub.itertuples()]
        axes[2].scatter(labels, sub["best_drop"], s=80, label=mode, color=color, alpha=0.8)
    axes[2].set_ylabel("largest clean pulse-window drop")
    axes[2].set_title("Pulse/SOC ranking")
    axes[2].tick_params(axis="x", rotation=30)
    axes[2].grid(axis="y", alpha=0.2)
    axes[2].legend()

    fig.savefig(ANALYSIS_DIR / "goal5_pulse_window_maps.png", bbox_inches="tight")
    plt.close(fig)


def write_markdown_summary(df: pd.DataFrame, paired: pd.DataFrame) -> None:
    top = df.sort_values(["avg_Mz_final", "tail_mean_avg_Mz", "avg_Mz_min"]).iloc[0]
    top_pulse = df.sort_values(
        [
            "clean_pre_pulse_positive",
            "clean_large_pulse_drop",
            "clean_pulse_window_strong_cross",
            "delayed_strong_cross_after_pulse",
            "pulse_induced_final_sign_change",
            "clean_pulse_window_sign_cross",
            "clean_pulse_response_rank_score",
            "post_pulse_min_avg_Mz",
            "pulse_window_min_avg_Mz",
        ],
        ascending=[False, False, False, False, False, False, False, True, True],
    ).iloc[0]
    sign_count = int(df["final_sign_change"].sum())
    strong_count = int(df["strong_final_reversal"].sum())
    sustained_count = int(df["sustained_tail_reversal"].sum())
    deep_transient_count = int(df["deep_transient_reversal"].sum())
    pre_pulse_zero_count = int(df["pre_pulse_cross_zero"].sum())
    pre_pulse_strong_count = int(df["pre_pulse_cross_m0p8"].sum())
    pulse_induced_sign_count = int(df["pulse_induced_final_sign_change"].sum())
    pulse_induced_strong_count = int(df["pulse_induced_strong_final"].sum())
    clean_large_drop_count = int(df["clean_large_pulse_drop"].sum())
    clean_pulse_sign_count = int(df["clean_pulse_window_sign_cross"].sum())
    clean_pulse_strong_count = int(df["clean_pulse_window_strong_cross"].sum())
    delayed_strong_count = int(df["delayed_strong_cross_after_pulse"].sum())
    paired_sign = int(paired["any_final_sign_change"].sum())
    paired_clean_large_drop = int(paired["any_clean_large_pulse_drop"].sum())
    paired_clean_pulse_sign = int(paired["any_clean_pulse_window_sign_cross"].sum())
    paired_clean_pulse_strong = int(paired["any_clean_pulse_window_strong_cross"].sum())
    paired_delayed_strong = int(paired["any_delayed_strong_cross_after_pulse"].sum())
    paired_strong = int(paired["any_strong_final"].sum())
    paired_sustained = int(paired["any_sustained_tail"].sum())
    full_minus_mean = float(paired["delta_final_full_minus_neq"].mean())
    full_better = int((paired["delta_final_full_minus_neq"] < 0.0).sum())
    selected = select_pulse_window_candidates(df)

    by_jani = (
        df.groupby("j_ani")
        .agg(
            final_sign=("final_sign_change", "sum"),
            clean_pulse_sign=("clean_pulse_window_sign_cross", "sum"),
            clean_pulse_strong=("clean_pulse_window_strong_cross", "sum"),
            clean_large_drop=("clean_large_pulse_drop", "sum"),
            delayed_strong=("delayed_strong_cross_after_pulse", "sum"),
            strong_final=("strong_final_reversal", "sum"),
            sustained_tail=("sustained_tail_reversal", "sum"),
            best_final=("avg_Mz_final", "min"),
            best_pulse_drop=("pulse_window_drop_from_start", "max"),
            best_pulse_min=("pulse_window_min_avg_Mz", "min"),
        )
        .reset_index()
    )

    lines = [
        "# Goal 5 Completed-Scan Analysis",
        "",
        "- prompt: Hagamos el analisis",
        "- procedures: Loaded the completed Goal 5 scan summary and all available avg_Mz trajectories; added final-sign, strong-final, transient-depth, sustained-tail, clean pre-pulse, pulse-window response, and delayed/post-pulse inversion criteria; generated ranked run tables, paired Hs_update_mode comparisons, and summary figures.",
        f"- generated data: {ANALYSIS_DIR / 'goal5_ranked_runs.tsv'}",
        f"- generated data: {ANALYSIS_DIR / 'goal5_pulse_window_ranked_runs.tsv'}",
        f"- generated data: {ANALYSIS_DIR / 'goal5_selected_pulse_window_candidates.tsv'}",
        f"- generated data: {ANALYSIS_DIR / 'goal5_paired_hs_mode_comparison.tsv'}",
        f"- generated data: {ANALYSIS_DIR / 'goal5_by_j_ani.tsv'}",
        f"- generated figure: {ANALYSIS_DIR / 'goal5_top_trajectories.png'}",
        f"- generated figure: {ANALYSIS_DIR / 'goal5_top_pulse_window_trajectories.png'}",
        f"- generated figure: {ANALYSIS_DIR / 'goal5_hs_mode_comparison.png'}",
        f"- generated figure: {ANALYSIS_DIR / 'goal5_parameter_maps.png'}",
        f"- generated figure: {ANALYSIS_DIR / 'goal5_pulse_window_maps.png'}",
        "",
        "## Main counts",
        "",
        f"- Branch runs analyzed: {len(df)}",
        f"- Paired A/B parameter points analyzed: {len(paired)}",
        f"- Final sign-changing branch runs: {sign_count}",
        f"- Strong final reversals, avg_Mz_final <= -0.8: {strong_count}",
        f"- Deep transient reversals, avg_Mz_min <= -0.8: {deep_transient_count}",
        f"- Sustained-tail reversals, tail_mean <= -0.8 and at least 80% of final 10% below -0.8: {sustained_count}",
        f"- Runs crossing avg_Mz=0 before pulse_start: {pre_pulse_zero_count}",
        f"- Runs crossing avg_Mz=-0.8 before pulse_start: {pre_pulse_strong_count}",
        f"- Pulse-induced final sign changes, positive at pulse_start and negative at final time: {pulse_induced_sign_count}",
        f"- Pulse-induced strong final reversals, positive at pulse_start and avg_Mz_final <= -0.8: {pulse_induced_strong_count}",
        f"- Clean large pulse drops, pre-pulse-positive with pulse-window drop >= {CRITERIA.large_pulse_drop_threshold}: {clean_large_drop_count}",
        f"- Clean pulse-window sign crosses, positive at pulse_start and crossing avg_Mz=0 during the pulse: {clean_pulse_sign_count}",
        f"- Clean pulse-window strong crosses, positive at pulse_start and crossing avg_Mz=-0.8 during the pulse: {clean_pulse_strong_count}",
        f"- Delayed strong crosses after pulse, clean-positive at pulse_start and first avg_Mz=-0.8 after pulse end: {delayed_strong_count}",
        f"- Parameter points with any final sign change: {paired_sign}",
        f"- Parameter points with any clean large pulse drop: {paired_clean_large_drop}",
        f"- Parameter points with any clean pulse-window sign cross: {paired_clean_pulse_sign}",
        f"- Parameter points with any clean pulse-window strong cross: {paired_clean_pulse_strong}",
        f"- Parameter points with any delayed strong post-pulse cross: {paired_delayed_strong}",
        f"- Parameter points with any strong final reversal: {paired_strong}",
        f"- Parameter points with any sustained-tail reversal: {paired_sustained}",
        "",
        "## Best case",
        "",
        f"- run_id: {top['run_id']}",
        f"- Hs_update_mode: {top['Hs_update_mode']}",
        f"- parameters: pulse_amp={top['pulse_amp']}, gamma_so={top['gamma_so']}, j_sd={top['j_sd']}, j_ex={top['j_ex']}, j_ani={top['j_ani']}, damping={top['damping']}",
        f"- avg_Mz_final: {top['avg_Mz_final']:.10f}",
        f"- avg_Mz_min: {top['avg_Mz_min']:.10f}",
        f"- avg_Mz_at_pulse_start: {top['avg_Mz_at_pulse_start']:.10f}",
        f"- avg_Mz_at_pulse_end: {top['avg_Mz_at_pulse_end']:.10f}",
        f"- tail_mean_avg_Mz: {top['tail_mean_avg_Mz']:.10f}",
        f"- tail_fraction_below_-0.8: {top['tail_fraction_below_m0p8']:.4f}",
        f"- first_cross_zero_time: {top['first_cross_zero_time']:.4f}",
        f"- first_cross_-0.8_time: {top['first_cross_m0p8_time']:.4f}",
        "",
        "## Best pulse-window case",
        "",
        f"- run_id: {top_pulse['run_id']}",
        f"- Hs_update_mode: {top_pulse['Hs_update_mode']}",
        f"- parameters: pulse_amp={top_pulse['pulse_amp']}, gamma_so={top_pulse['gamma_so']}, j_sd={top_pulse['j_sd']}, j_ex={top_pulse['j_ex']}, j_ani={top_pulse['j_ani']}, damping={top_pulse['damping']}",
        f"- avg_Mz_at_pulse_start: {top_pulse['avg_Mz_at_pulse_start']:.10f}",
        f"- pulse_window_min_avg_Mz: {top_pulse['pulse_window_min_avg_Mz']:.10f}",
        f"- pulse_window_time_min_avg_Mz: {top_pulse['pulse_window_time_min_avg_Mz']:.4f}",
        f"- pulse_window_drop_from_start: {top_pulse['pulse_window_drop_from_start']:.10f}",
        f"- avg_Mz_at_pulse_end: {top_pulse['avg_Mz_at_pulse_end']:.10f}",
        f"- avg_Mz_final: {top_pulse['avg_Mz_final']:.10f}",
        f"- clean_pulse_window_sign_cross: {top_pulse['clean_pulse_window_sign_cross']}",
        f"- clean_pulse_window_strong_cross: {top_pulse['clean_pulse_window_strong_cross']}",
        "",
        "## Selected pulse-window candidates",
        "",
        f"- Selection rule: require clean pre-pulse-positive dynamics with no pre-pulse zero crossing and pulse-window drop >= {CRITERIA.large_pulse_drop_threshold}; then rank by pulse-window drop, with priority for in-pulse strong/sign inversion, delayed strong post-pulse inversion, and negative final state.",
        f"- Selected watchlist size: {len(selected)}",
        f"- Strict strong clean in-pulse candidates, pulse_window_min_avg_Mz <= -0.8: {int(selected['clean_pulse_window_strong_cross'].sum())}",
        f"- Delayed strong post-pulse candidates in selected watchlist: {int(selected['delayed_strong_cross_after_pulse'].sum())}",
        f"- Negative-final candidates in selected watchlist: {int(selected['pulse_induced_final_sign_change'].sum())}",
        "",
        "## Hs_update_mode comparison",
        "",
        f"- Mean full_spin_density minus nonequilibrium_only final avg_Mz shift: {full_minus_mean:.6f}",
        f"- full_spin_density produced a lower final avg_Mz in {full_better}/{len(paired)} paired points.",
        "",
        "## By anisotropy",
        "",
    ]
    for row in by_jani.itertuples(index=False):
        lines.append(
            f"- j_ani={row.j_ani}: final_sign={int(row.final_sign)}, "
            f"clean_large_drop={int(row.clean_large_drop)}, "
            f"clean_pulse_sign={int(row.clean_pulse_sign)}, "
            f"clean_pulse_strong={int(row.clean_pulse_strong)}, "
            f"delayed_strong={int(row.delayed_strong)}, "
            f"strong_final={int(row.strong_final)}, sustained_tail={int(row.sustained_tail)}, "
            f"best_final={row.best_final:.6f}, best_pulse_drop={row.best_pulse_drop:.6f}, "
            f"best_pulse_min={row.best_pulse_min:.6f}"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- The original sign-change count overstates true switching: most sign-changing runs finish close to zero rather than near -1.",
            "- For choosing the next scan, the best criterion is now clean pre-pulse-positive behavior plus a large pulse-window drop, with in-pulse or final inversion used as a positive ranking signal rather than a mandatory filter.",
            "- The only strong final switching case is the full_spin_density branch at pulse_amp=2.0, gamma_so=0.7, j_sd=0.3, j_ex=-0.01, j_ani=0.0, damping=0.3.",
            "- That best case is not a clean pulse-induced switching event under this protocol: it has already crossed avg_Mz=0 and avg_Mz=-0.8 before pulse_start=1000, and it fails the strict sustained-tail criterion.",
            "- Low anisotropy is still the controlling feature. No j_ani >= 0.08 case even changes final sign, and no j_ani > 0.0 case reaches strong final reversal.",
            "- The next local scan should first control the pre-pulse relaxation problem, for example by shortening pulse_start or initializing from a relaxed reference, before widening pulse_amp around 1.8-2.3.",
            "- status: done",
            "",
        ]
    )
    (ANALYSIS_DIR / "summary.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    ensure_dir(ANALYSIS_DIR)
    setup_plot_style()
    summary = read_summary()
    analyzed = add_analysis_columns(summary)
    paired = pair_modes(analyzed)
    save_tables(analyzed, paired)
    plot_top_trajectories(analyzed)
    plot_top_pulse_window_trajectories(analyzed)
    plot_hs_mode_comparison(paired)
    plot_parameter_maps(analyzed)
    plot_pulse_window_maps(analyzed)
    write_markdown_summary(analyzed, paired)
    print(f"Wrote analysis to {ANALYSIS_DIR}")


if __name__ == "__main__":
    main()
