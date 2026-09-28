from pathlib import Path
import csv
import math

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent.parent.parent
SHORT_PATH = ROOT / "projects/2026-04-13-fixed-llg-no-pulse-relaxation/results/fixed_llg_no_pulse_20260413T201501/fixedllg_gso0p05/observables.tsv"
LONG_PATH = ROOT / "projects/2026-04-13-rhoeq-longtime-relaxation/results/rhoeq_longtime_20260413T205653/pade20_gso0p05_t2000/observables.tsv"
OUT_DIR = Path(__file__).resolve().parent / "results"
OUT_PATH = OUT_DIR / "rhoeq_relaxation_tendency.png"


def load_series(path: Path):
    with path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        times = []
        sigma_absmax = []
        for row in reader:
            t = float(row["t"])
            site_max = 0.0
            for site in range(1, 5):
                dx = float(row[f"sigma_delta_x_site{site}"])
                dy = float(row[f"sigma_delta_y_site{site}"])
                dz = float(row[f"sigma_delta_z_site{site}"])
                site_max = max(site_max, math.sqrt(dx * dx + dy * dy + dz * dz))
            times.append(t)
            sigma_absmax.append(site_max)
    return times, sigma_absmax


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    t_short, y_short = load_series(SHORT_PATH)
    t_long, y_long = load_series(LONG_PATH)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)

    ax = axes[0]
    ax.plot(t_short, y_short, label="hasta t=1000", lw=2.0, color="#c44e52")
    ax.plot(t_long, y_long, label="hasta t=2000", lw=2.0, color="#4c72b0")
    ax.set_xlabel("t")
    ax.set_ylabel(r"max$_i$ |$\sigma_{\Delta,i}(t)$|")
    ax.set_title("Relajacion sin pulso")
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False)

    ax = axes[1]
    ax.plot(t_long, y_long, lw=2.0, color="#4c72b0")
    ax.set_xlim(800, 2000)
    ax.set_yscale("log")
    ax.set_xlabel("t")
    ax.set_ylabel(r"max$_i$ |$\sigma_{\Delta,i}(t)$|")
    ax.set_title("Cola tardia, corrida larga")
    ax.grid(True, which="both", alpha=0.3)
    ax.axhline(y_long[-1], color="#dd8452", ls="--", lw=1.5, label=f"valor final = {y_long[-1]:.4f}")
    ax.legend(frameon=False)

    fig.suptitle("Decaimiento de sigma_delta hacia cero para gamma_so = 0.05, N_lambda2 = 20")
    fig.savefig(OUT_PATH, dpi=180)
    print(OUT_PATH)


if __name__ == "__main__":
    main()
