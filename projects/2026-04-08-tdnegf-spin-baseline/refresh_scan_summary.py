"""Refresh a queued-scan summary table from manifest and per-run outputs.

This is intentionally non-invasive: it does not touch the validated single-run
driver or the live queued controller. It merges `case_manifest.tsv` with any
available `summary.json` / `observables.tsv` files and writes one row per run.
"""

from __future__ import annotations

from pathlib import Path
import csv
import json
import sys


SUMMARY_HEADERS = [
    "run_id", "status", "initial_configuration",
    "Hs_update_mode",
    "Nx", "Ny",
    "pulse_amp", "bias_amplitude", "gamma_so", "γso",
    "j_sd", "j_ex", "j_ani", "damping", "Bx",
    "left_bias_scale", "right_bias_scale",
    "pulse_start", "pulse_ramp", "pulse_width", "t_end", "t_final",
    "completed_steps",
    "current_left_final", "current_right_final", "current_peak_abs", "current_balance_peak_abs",
    "avg_M_final", "avg_Mz_initial", "avg_Mz_final", "avg_Mz_min", "avg_Mz_max",
    "neel_final", "neel_z_peak_abs",
    "max_site_charge_deviation_from_1", "max_classical_spin_norm_deviation",
    "max_sigma_noneq_abs", "max_sigma_eq_abs", "max_sigma_delta_abs",
    "pre_pulse_sample_count", "pre_pulse_duration_realized",
    "pre_pulse_sigma_noneq_absmax", "pre_pulse_sigma_eq_absmax", "pre_pulse_sigma_delta_absmax",
    "full_reversal_detected", "partial_reversal_flag",
    "llg_drive_convention", "equilibrium_reference_saved", "nonequilibrium_minus_equilibrium_saved",
    "results_dir", "error", "timestamp_utc",
    "params_present", "summary_present", "observables_present",
]


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_json(path: Path) -> dict[str, object]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def build_row(manifest_row: dict[str, str], run_dir: Path) -> dict[str, object]:
    params_path = run_dir / "params.json"
    summary_path = run_dir / "summary.json"
    observables_path = run_dir / "observables.tsv"

    row: dict[str, object] = {key: "missing" for key in SUMMARY_HEADERS}
    row.update(manifest_row)
    row["bias_amplitude"] = manifest_row.get("pulse_amp", "missing")
    row["γso"] = manifest_row.get("gamma_so", "missing")
    row["results_dir"] = manifest_row.get("results_dir", str(run_dir))
    row["status"] = "pending"
    row["params_present"] = params_path.is_file()
    row["summary_present"] = summary_path.is_file()
    row["observables_present"] = observables_path.is_file()

    if params_path.is_file():
        row["status"] = "dispatched"
    if observables_path.is_file():
        row["status"] = "running_or_complete"
    if summary_path.is_file():
        summary = load_json(summary_path)
        for key, value in summary.items():
            if key in row:
                row[key] = value

    return row


def refresh_scan_summary(results_root: Path) -> Path:
    manifest_path = results_root / "case_manifest.tsv"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Missing manifest: {manifest_path}")

    manifest = load_tsv(manifest_path)
    rows = []
    for manifest_row in manifest:
        run_dir = Path(manifest_row["results_dir"])
        rows.append(build_row(manifest_row, run_dir))

    summary_path = results_root / "scan_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_HEADERS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return summary_path


def main(argv: list[str]) -> int:
    if len(argv) > 2:
        print("Usage: refresh_scan_summary.py [results_root]", file=sys.stderr)
        return 2

    if len(argv) == 2:
        results_root = Path(argv[1]).resolve()
    else:
        results_root = (Path(__file__).resolve().parent / "results").resolve()
        roots = sorted(path for path in results_root.glob("queued_parameter_scan_*") if path.is_dir())
        if not roots:
            raise FileNotFoundError(f"No queued_parameter_scan_* roots found under {results_root}")
        results_root = roots[-1]

    summary_path = refresh_scan_summary(results_root)
    print(summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
