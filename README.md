# Driven-Switching

Pulse-driven switching of classical magnetic moments in an open electronic
system, simulated with the TDNEGF+LLG stack (April 2026). The campaign was run
through the FermiLink/Codex agent workflow, so the repository is organised as a
sequence of dated experiments plus a single "unified memory" file.

## Physical system

- Small square lattice (`Nx`, `Ny` of order 2-3) with two metallic leads.
- Rashba spin-orbit coupling (`gamma_so`), s-d exchange `j_sd` between
  electrons and classical spins, exchange `j_ex` and anisotropy `j_ani` among
  the classical spins, Gilbert damping.
- A bias pulse on the right lead (`pulse_amp`, `pulse_start`, `pulse_width`)
  drives the electrons; the question is whether the resulting spin-transfer
  torque produces a clean, retained reversal of the average magnetisation
  `avg_Mz`.
- Two update conventions for the spin Hamiltonian were compared:
  `hsneq` (`Hs_update_mode=nonequilibrium_only`) and `hsfull`
  (`full_spin_density`).

Starting point: the TDNEGF example notebook
`examples/notebooks/02_spin_dynamics_and_electrons.ipynb`.

## Layout

| Path | Content |
|---|---|
| `goal1.md` … `goal6.md` | Task specifications for each phase of the campaign |
| `AGENTS.md` | FermiLink agent policy (package routing, memory rules) |
| `projects/memory.md` | Unified memory: file map, simulation history, key results, parameter provenance |
| `projects/<date>-<name>/` | One experiment each: run scripts (`.jl`, `.py`, `.sh`), `summary.md`, metrics `.json`, figures |
| `projects/2026-04-16-goal6-latex-report/report.tex` | Self-contained LaTeX report of the Goal 5/6 results |

Simulation outputs (`projects/*/results/`, about 48 GB of `observables.tsv`
and `params.json`) and logs are **not** in git. They live only on
`keldysh:~/jalil_codes/Projects_2026/Driven-Switching/`.

## Main campaigns and findings

- **Baseline and pilot (2026-04-08/09).** Reproduced the notebook example;
  no reversal in short runs.
- **Goal 5 (2026-04-10).** 160 parameter points x 2 `Hs_update_mode` branches
  (320 jobs, `t_end=4000`). Sign changes only at `j_ani=0`; strongest final
  reversal `avg_Mz_final=-0.907`.
- **Goal 6 (2026-04-10/13).** Focused follow-up, 117 traces. Best clean
  in-pulse retained switch: `goal6_A_main_0003…hsneq`
  (`pulse_window_min≈-0.92`, `avg_Mz_final≈-1.0`).
- **No-pulse controls (2026-04-14/15).** Matched `bias_amplitude=0` twins of
  the winners do not switch, confirming the pulse as the cause. Branch B
  (`hsfull`) families are baseline-ambiguous.
- **Equilibrium diagnostics (2026-04-13/18).** The static `rho_eq(H_ab)` is not
  a stationary point of the enlarged block TDNEGF equations; a coupled
  zero-bias pre-relaxation plus prepared-state reference reduces the pre-pulse
  residual `sigma-sigma_eq` from ~0.18 to ~0.018.

Each experiment's `summary.md` states its prompt, procedure, generated data,
results and status. `projects/memory.md` is the authoritative index.

## Requirements

Julia with the `TDNEGF` package (this repository was run against the
FermiLink-installed copy of `TDNEGF`), Python 3 with pandas/matplotlib for the
analysis scripts.
