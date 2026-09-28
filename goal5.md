# Goal 5

Refine the next switching scan so it explores a somewhat broader threshold-seeking control space, but still keeps the same simple all-points-in-parallel workflow as before.

## Added parameters

Keep the original Goal 5 controls:

- `pulse_amp`
- `gamma_so`
- `j_sd`
- `j_ex`
- `Bx`

Add five more parameters to the primary planning space:

- `j_ani`
- `damping`
- `pulse_start`
- `pulse_ramp`
- `pulse_width`

Also add one discrete modeling choice to compare in the next Goal 5 tests:

- `Hs_update_mode ∈ {nonequilibrium_only, full_spin_density}`

Interpret this as:

- `nonequilibrium_only`: update `H_s` with the nonequilibrium part only, `σ(t) - σ_eq`
- `full_spin_density`: update `H_s` with the full instantaneous spin density, `σ(t)`

So Goal 5 now carries:

- `10` primary switching/protocol controls:
  `pulse_amp`, `gamma_so`, `j_sd`, `j_ex`, `Bx`, `j_ani`, `damping`, `pulse_start`, `pulse_ramp`, `pulse_width`
- `1` discrete model branch:
  `Hs_update_mode ∈ {nonequilibrium_only, full_spin_density}`

For the next batch, do not try to sweep all 10 controls at once. Keep the pulse timing and bias partition configurable in the driver, but fixed in the first improved batch.

## Bias convention change

Use a right-lead-only pulse as the next default trial:

- `left_bias_scale = 0.0`
- `right_bias_scale = 1.0`

This replaces the old anti-symmetric convention:

- old: `Δ_left(t) = +Δ(t)`, `Δ_right(t) = -Δ(t)`
- new Goal 5 trial: `Δ_left(t) = 0`, `Δ_right(t) = +Δ(t)`

Keep this as the default improved protocol for now. If a promising corner appears later, the lead partition can be widened in a follow-up batch, but it is not part of the immediate expanded sweep.

## `H_s` update comparison

In addition to the lead-bias change, Goal 5 should explicitly compare both cases for the spin-Hamiltonian update:

- remove the equilibrium part before `update_H_s!`
- keep the equilibrium part in `update_H_s!`

Recommended interpretation for the next targeted tests:

- Case A: `Hs_update_mode = nonequilibrium_only` using `σx_delta = σx_i - σx_i_eq`
- Case B: `Hs_update_mode = full_spin_density` using `σx_full = σx_i`

This should be treated as a first-order branch in the next pilot batch, because changing the lead protocol and changing the magnetic feedback prescription at the same time may alter both torque strength and transient current balance.

For execution planning, keep the same global parallel launch pattern used before:

- the full selected batch of parameter points runs in parallel across the queue
- for each parameter point, launch both `Hs_update_mode` variants as its paired A/B jobs
- one job with `Hs_update_mode = nonequilibrium_only`
- one job with `Hs_update_mode = full_spin_density`
- allow at most `34` jobs running simultaneously
- if the full batch is larger than `34` jobs, keep the remaining jobs in the loop/queue and launch them as earlier jobs finish

So the intended structure is not "one serial scan, each with an internal comparison", but rather "all chosen scan points in parallel, and each scan point expands to a parallel A/B pair".

This keeps wall-clock time aligned with the previous workflow and preserves a clean A/B comparison under identical non-`Hs_update_mode` parameters.

## Recommended improved batch

Keep the queue logic from the earlier workflow:

- every selected parameter point runs in parallel
- every point is duplicated into the two `Hs_update_mode` branches
- total jobs = `(number of selected parameter points) x 2`

To keep Goal 5 simple but better than the earlier batch, use a 6-parameter sweep instead of reopening many staged strips immediately.

Recommended fixed settings:

- `t_end = 4000.0`
- `pulse_start = 1000.0`
- `pulse_ramp = 5.0`
- `pulse_width = 1500.0`
- `Bx = 0.0`
- `left_bias_scale = 0.0`
- `right_bias_scale = 1.0`

Recommended parallel sweep:

- `pulse_amp`: `1.5`, `2.0`
- `gamma_so`: `0.7`, `0.9`
- `j_sd`: `0.2`, `0.3`
- `j_ex`: `-0.02`, `-0.01`
- `j_ani`: `0.00`, `0.03`, `0.05`, `0.08`, `0.10`
- `damping`: `0.3`, `0.5`

This gives:

- `2 x 2 x 2 x 2 x 5 x 2 = 160` parameter points
- with the A/B `Hs_update_mode` split: `160 x 2 = 320` jobs
- execution policy: run up to `34` jobs at the same time and let the remaining `286` jobs continue in loop as slots free up

This keeps the same style as before:

- one compact batch
- all selected points belong to one parallel batch submission plan
- broader `j_ani` coverage without reopening a much more complicated staged campaign
- still focused on the knobs most likely to matter for switching

If this 64-point batch still shows no clear crossing, the next widening step should be small and local:

- either add one more `pulse_amp` value, for example `2.5`
- or open a short pulse-timing strip around the best corner
- but do not jump immediately to a full multi-stage campaign

## Suggested baseline for the next targeted batch

- `t_end = 4000.0`
- `pulse_start = 1000.0`
- `pulse_ramp = 5.0`
- `pulse_width = 1500.0`
- `j_ani = 0.05`
- `damping = 0.5`
- `Bx = 0.0`
- `left_bias_scale = 0.0`
- `right_bias_scale = 1.0`
- `Hs_update_mode ∈ {nonequilibrium_only, full_spin_density}`
