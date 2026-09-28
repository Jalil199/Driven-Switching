# TDNEGF Spin-Electron Baseline Workflow Note

## Goal

Reproduce the existing notebook workflow in `examples/notebooks/02_spin_dynamics_and_electrons.ipynb` with minimal changes, then prepare a single-run driver suitable for later parameter sweeps.

## Notebook setup

- Electronic device geometry: square-lattice central region with `Nx=3`, `Ny=2`, `Nσ=2`, `N_orb=1`, so there are 6 lattice sites and 12 electronic basis states.
- Electronic Hamiltonian: `build_H_ab` constructs nearest-neighbor hopping on the square lattice with hopping `γ=1.0` and spin-orbit coupling `γso=0.5+0im`.
- Leads: two semi-infinite square-lattice leads are attached to the first and last `x` columns.
  - Left lead uses `build_ξ_an(...; xcol=1, y_coup=1:Ny)`.
  - Right lead uses `build_ξ_an(...; xcol=Nx, y_coup=1:Ny)`.
- Lead self-energy ingredients: `load_poles_square(49, 20)`, then `build_Σᴸ_nλ`, `build_Σᴳ_nλ`, and `build_χ_nλ` with `β=33.0` and `γ=1.0`.
- Classical-spin lattice: Sunny `System` on a `3x2x1` inhomogeneous lattice, initialized ferromagnetically because the notebook run uses `j_ex=-0.05`.
- Spin parameters used in the notebook run:
  - `j_sd=0.1`
  - `j_ex=-0.05`
  - `j_ani=0.1`
  - `Bx=0.0`
  - LLG integrator: `Langevin(Δt; damping=0.5, kT=0.0)`
- Pulse/bias protocol:
  - `Δ_left(t) = smooth_pulse(t-100; ti=10, Δ=1000) * 1.0`
  - `Δ_blocks = [+Δ_left, -Δ_left]`
  - For a run ending near `t=200`, this means the bias starts turning on at `t=100`, reaches its plateau near `t=110`, and does not turn off within the simulated window.

## Runtime workflow

1. Build the static electronic model and block-lead data.
2. Construct `ExperimentalBlockRHSParams(H_ab, blocks, Δ_blocks, p_model)` and a zero initial state.
3. Build the Sunny classical-spin system.
4. Propagate electrons with `ODEProblem(..., eom_tdnegf_blocks!, ...)` and `Vern7`.
5. After each electron step, advance the classical spins with Sunny `step!`.
6. Extract observables from `pointer_blocks(...)`.
7. Update the classical-spin effective field from the nonequilibrium electron spin density.
8. Update the electronic Hamiltonian through the `s-d` coupling to the current classical spin texture.
9. Update the lead bias `Δ_blocks`.

## Source-function mapping

- Model construction:
  - `build_H_ab` in `src/hamiltonians.jl`
  - `load_poles_square` in `src/poles.jl`
  - `build_Σᴸ_nλ`, `build_Σᴳ_nλ`, `build_χ_nλ`, `build_ξ_an` in `src/selfenergy.jl`
  - `ExperimentalBlockRHSParams` in `src/eom_tdnegf.jl`
- Time propagation:
  - `eom_tdnegf_blocks!` in `src/eom_tdnegf.jl`
  - DifferentialEquations `ODEProblem` and `Vern7` drive the electron ODE integration
- Lead coupling and self-energy data:
  - `SelfEnergyBlock` plus the `Σ/χ/ξ` builders above
  - Heterogeneous auxiliary layout handled by `pointer_blocks` in `src/selfenergy_aux_layout.jl`
- Classical-spin updates:
  - Notebook-local `update_H_s!` updates Sunny local fields using the nonequilibrium electron spin density
  - `update_H_e!` in `src/hamiltonians.jl` updates the electronic Hamiltonian from the classical spin texture
  - Sunny `step!` advances the classical spins
- Observable extraction:
  - `obs_n_i!`, `obs_σ_i!`, `obs_σ_i_eq!`, `obs_s_i!`, `obs_Ixα!` in `src/observables.jl`
  - `ρ_eq` in `src/stationary.jl` generates the equilibrium reference spin density used by the notebook

## Observables used by the notebook

- Site-resolved charge density: `obs.n_i`
- Site-resolved nonequilibrium electron spin density: `obs.σx_i`
- Site-resolved equilibrium electron spin density: `obs.σx_i_eq`
- Site-resolved classical spins: `obs.sx_i`
- Lead charge currents: `obs.Iα`
- Lead spin currents: `obs.Iαx`
- The notebook plots mainly:
  - classical-spin `z` components on selected sites,
  - lead charge current versus time,
  - lead spin current components versus time,
  - site-resolved electron `s^z`.
