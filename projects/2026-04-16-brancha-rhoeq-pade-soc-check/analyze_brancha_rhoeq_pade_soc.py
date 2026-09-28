from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parent
RESULTS_ROOT = PROJECT_DIR / "results" / "brancha_fixed_llg_pade_soc_sweep"
MANIFEST_JSON_PATH = PROJECT_DIR / "cases_manifest.json"
OUT_TABLE = PROJECT_DIR / "brancha_pade_soc_summary.tsv"
OUT_FIG = PROJECT_DIR / "brancha_pade_soc_summary.png"
SUMMARY_MD = PROJECT_DIR / "summary.md"


def load_cases() -> list[dict]:
    return json.loads(MANIFEST_JSON_PATH.read_text())


def load_summary(case: dict) -> dict:
    path = RESULTS_ROOT / case["run_id"] / "summary.json"
    data = json.loads(path.read_text())
    data["run_id"] = case["run_id"]
    data["gamma_so"] = case["gamma_so"]
    data["rhoeq_pade_poles"] = case["N_λ2"]
    data["results_dir"] = case["results_dir"]
    return data


def sigma_delta_absmax_trace(observables_path: Path) -> tuple[list[float], list[float]]:
    times: list[float] = []
    values: list[float] = []
    with observables_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        sigma_cols = [name for name in reader.fieldnames if name.startswith("sigma_delta_")]
        for row in reader:
            times.append(float(row["t"]))
            values.append(max(abs(float(row[name])) for name in sigma_cols))
    return times, values


def build_dataframe(cases: list[dict]) -> pd.DataFrame:
    rows = [load_summary(case) for case in cases]
    df = pd.DataFrame(rows)
    df = df.sort_values(["gamma_so", "rhoeq_pade_poles"]).reset_index(drop=True)

    final_diff = []
    late_diff = []
    for gamma_so, group in df.groupby("gamma_so", sort=True):
        baseline = group[group["rhoeq_pade_poles"] == 20].iloc[0]
        for _, row in group.iterrows():
            final_diff.append(row["final_sigma_delta_absmax"] - baseline["final_sigma_delta_absmax"])
            late_diff.append(row["late_sigma_delta_absmax_mean"] - baseline["late_sigma_delta_absmax_mean"])
    df["final_sigma_delta_diff_vs_pade20_same_soc"] = final_diff
    df["late_sigma_delta_mean_diff_vs_pade20_same_soc"] = late_diff
    return df


def make_figure(df: pd.DataFrame, cases: list[dict]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))

    for gamma_so, group in df.groupby("gamma_so", sort=True):
        label = f"gamma_so={gamma_so:g}"
        axes[0].plot(group["rhoeq_pade_poles"], group["final_sigma_delta_absmax"], marker="o", linewidth=1.8, label=label)
        axes[1].plot(group["rhoeq_pade_poles"], group["late_sigma_delta_absmax_mean"], marker="o", linewidth=1.8, label=label)

    focus_cases = [case for case in cases if case["gamma_so"] == 0.9]
    focus_cases.sort(key=lambda case: case["N_λ2"])
    for case in focus_cases:
        obs_path = RESULTS_ROOT / case["run_id"] / "observables.tsv"
        times, values = sigma_delta_absmax_trace(obs_path)
        axes[2].plot(times, values, linewidth=1.4, label=f"gso=0.9, pade={case['N_λ2']}")

    axes[0].set_title("Final pre-pulse residual")
    axes[0].set_xlabel("Padé poles N_λ2")
    axes[0].set_ylabel("final |sigma - sigma_eq| absmax")
    axes[0].grid(alpha=0.3)

    axes[1].set_title("Late-time mean residual")
    axes[1].set_xlabel("Padé poles N_λ2")
    axes[1].set_ylabel("late mean |sigma - sigma_eq| absmax")
    axes[1].grid(alpha=0.3)

    axes[2].set_title("gso=0.9 time trace")
    axes[2].set_xlabel("t")
    axes[2].set_ylabel("|sigma - sigma_eq| absmax")
    axes[2].grid(alpha=0.3)

    axes[0].legend(fontsize=8)
    axes[2].legend(fontsize=8)
    fig.suptitle("Branch A-like fixed-LLG no-pulse check: Padé vs Rashba SOC", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=200, bbox_inches="tight")
    plt.close(fig)


def write_summary(df: pd.DataFrame) -> None:
    best_final = df.loc[df["final_sigma_delta_absmax"].idxmin()]
    best_late = df.loc[df["late_sigma_delta_absmax_mean"].idxmin()]
    gso09 = df[df["gamma_so"] == 0.9].sort_values("rhoeq_pade_poles")
    pade_effect = ", ".join(
        f"N_λ2={int(row.rhoeq_pade_poles)} -> final={row.final_sigma_delta_absmax:.6f}"
        for row in gso09.itertuples()
    )
    status = "done" if (df["status"] == "done").all() else "blocked"
    with SUMMARY_MD.open("w") as fh:
        fh.write("- prompt: quisiera ver si podemos hacer tender eso a cero antes del pulso. Probablemente es un problema relacionado con el calculo de la matriz densidad en equilibrio. chequeemos esto, primero miremos si mejora para mas polos de Pade, podemos tratar con 30 y 40 y ademas miremos si para soc mas pequenos van a cero, deprones el que no esta dando bien para un soc. deja corriendo estos casos en paralelo y monitorealo hasta que termine no cierres antes de eso\n")
        fh.write("- procedures: se tomó como referencia el caso `goal6_A_main_0003` de Branch A (`Nx=Ny=2`, `j_sd=0.3`, `j_ex=-0.02`, `j_ani=0.03`, `damping=0.3`, `pulse_start=1000`) y se construyó un probe sin pulso con espines fijos para aislar el papel de `rho_eq`. Se corrieron en paralelo 12 casos seriales locales con `gamma_so ∈ {0.9, 0.3, 0.1, 0.05}` y `N_λ2 ∈ {20, 30, 40}`, todos con `bias_amplitude=0` y propagación hasta `t=1000`, luego se resumieron los residuos `|sigma-sigma_eq|` al final y en la cola temporal.\n")
        fh.write(f"- generated data: {OUT_TABLE}\n")
        fh.write(f"- generated figure: {OUT_FIG}\n")
        fh.write(f"- debugging clue: para `gamma_so=0.9`, subir Padé de 20 a 30/40 cambia muy poco el residuo final pre-pulso ({pade_effect}); la reducción fuerte aparece al bajar `gamma_so`, lo que apunta más al régimen dinámico/modelo que a una falta simple de convergencia de `rho_eq` con Padé.\n")
        fh.write(f"- generated data: mejor residuo final en `gamma_so={best_final.gamma_so:g}`, `N_λ2={int(best_final.rhoeq_pade_poles)}` con `final_sigma_delta_absmax={best_final.final_sigma_delta_absmax:.6f}`; mejor cola tardía en `gamma_so={best_late.gamma_so:g}`, `N_λ2={int(best_late.rhoeq_pade_poles)}` con `late_sigma_delta_absmax_mean={best_late.late_sigma_delta_absmax_mean:.6f}`.\n")
        fh.write(f"- status: {status}\n")


def main() -> None:
    cases = load_cases()
    df = build_dataframe(cases)
    df.to_csv(OUT_TABLE, sep="\t", index=False)
    make_figure(df, cases)
    write_summary(df)


if __name__ == "__main__":
    main()
