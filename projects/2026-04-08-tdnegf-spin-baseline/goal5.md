# Goal 5

## Current execution target

Goal 5 now uses a broader threshold-seeking batch while keeping the same queue style as the earlier scan:

- all selected parameter points belong to one global batch,
- each point expands into two parallel `Hs_update_mode` variants,
- up to `34` jobs may run simultaneously,
- and the remaining jobs backfill as slots open.

## Control set

Primary controls:

- `pulse_amp`
- `gamma_so`
- `j_sd`
- `j_ex`
- `Bx`
- `j_ani`
- `damping`
- `pulse_start`
- `pulse_ramp`
- `pulse_width`

Discrete branch:

- `Hs_update_mode ∈ {nonequilibrium_only, full_spin_density}`

Interpretation:

- `nonequilibrium_only`: update `H_s` with `σ(t) - σ_eq`
- `full_spin_density`: update `H_s` with `σ(t)`

## Default improved protocol

- `t_end = 4000.0`
- `pulse_start = 1000.0`
- `pulse_ramp = 5.0`
- `pulse_width = 1500.0`
- `Bx = 0.0`
- `left_bias_scale = 0.0`
- `right_bias_scale = 1.0`

## Immediate improved batch

- `pulse_amp`: `1.5`, `2.0`
- `gamma_so`: `0.7`, `0.9`
- `j_sd`: `0.2`, `0.3`
- `j_ex`: `-0.02`, `-0.01`
- `j_ani`: `0.00`, `0.03`, `0.05`, `0.08`, `0.10`
- `damping`: `0.3`, `0.5`

Counts:

- `160` parameter points
- `320` total jobs after the `Hs_update_mode` A/B split
- `34` maximum concurrent jobs

## Targeted baseline pair

- `j_ani = 0.05`
- `damping = 0.5`
- `Hs_update_mode ∈ {nonequilibrium_only, full_spin_density}`

## Next widening rule

If no clear switching threshold appears, widen only locally:

- add `pulse_amp = 2.5`, or
- open a short pulse-timing strip around the best corner.
