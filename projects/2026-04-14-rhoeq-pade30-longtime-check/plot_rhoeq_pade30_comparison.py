from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
BASELINE = ROOT.parent / "2026-04-13-rhoeq-longtime-relaxation" / "results" / "rhoeq_longtime_20260413T205653" / "pade20_gso0p05_t2000" / "observables.tsv"


def load_curve(path: Path):
    times = []
    values = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        sigma_cols = [
            name
            for name in reader.fieldnames
            if name.startswith("sigma_delta_") and "_site" in name
        ]
        for row in reader:
            times.append(float(row["t"]))
            values.append(max(abs(float(row[name])) for name in sigma_cols))
    return times, values


def main():
    candidates = [path for path in (ROOT / "results").iterdir() if path.is_dir()]
    latest_root = max(candidates, key=lambda path: path.stat().st_mtime)
    pade30 = latest_root / "pade30_gso0p05_t2000" / "observables.tsv"
    out = latest_root / "rhoeq_pade20_vs_pade30_t2000.png"

    t20, y20 = load_curve(BASELINE)
    t30, y30 = load_curve(pade30)

    plt.figure(figsize=(9, 4.8))
    plt.plot(t20, y20, label="N_λ2 = 20", linewidth=2.0)
    plt.plot(t30, y30, label="N_λ2 = 30", linewidth=1.8, linestyle="--")
    plt.yscale("log")
    plt.xlabel("t")
    plt.ylabel(r"max$_i$ |sigma_delta,i(t)|")
    plt.title("Long-time no-pulse relaxation: Padé pole check")
    plt.grid(True, which="both", alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out, dpi=180)
    print(out)


if __name__ == "__main__":
    main()
