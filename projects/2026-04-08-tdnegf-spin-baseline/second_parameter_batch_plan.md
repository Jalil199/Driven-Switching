# Second Parameter Batch Proposal

## Goal 5 executable batch

The next Goal 5 scan should widen the threshold-seeking control space, but keep the same simple queue model as before:

- every selected parameter point belongs to one global parallel batch,
- every parameter point is duplicated into two `Hs_update_mode` jobs,
- the launcher may run up to `34` jobs simultaneously,
- and unfinished jobs simply wait in the queue until earlier jobs finish.

This keeps the old all-points-in-parallel workflow while making the magnetic-feedback prescription an explicit A/B comparison.

## Primary planning space

Goal 5 now carries ten primary controls plus one discrete branch:

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
- `Hs_update_mode ∈ {nonequilibrium_only, full_spin_density}`

Interpret the `Hs_update_mode` branches as:

- `nonequilibrium_only`: drive `update_H_s!` with `σ(t) - σ_eq`
- `full_spin_density`: drive `update_H_s!` with the full instantaneous `σ(t)`

Pulse timing and lead partition remain configurable in the driver, but fixed in the immediate improved batch below.

## Bias convention

The default Goal 5 protocol is now right-lead-only:

- `left_bias_scale = 0.0`
- `right_bias_scale = 1.0`

This replaces the earlier anti-symmetric trial. Lead partition should stay fixed in the first improved batch so the new bias convention and the new `H_s` update comparison can be interpreted cleanly.

## Recommended improved batch

Use the following fixed settings:

- `t_end = 4000.0`
- `pulse_start = 1000.0`
- `pulse_ramp = 5.0`
- `pulse_width = 1500.0`
- `Bx = 0.0`
- `left_bias_scale = 0.0`
- `right_bias_scale = 1.0`

Sweep these six controls in parallel:

- `pulse_amp`: `1.5`, `2.0`
- `gamma_so`: `0.7`, `0.9`
- `j_sd`: `0.2`, `0.3`
- `j_ex`: `-0.02`, `-0.01`
- `j_ani`: `0.00`, `0.03`, `0.05`, `0.08`, `0.10`
- `damping`: `0.3`, `0.5`

This gives:

- `2 x 2 x 2 x 2 x 5 x 2 = 160` parameter points
- `160 x 2 = 320` total jobs after the A/B `Hs_update_mode` split
- `34` maximum concurrent jobs

## Why this batch is the right next step

- It broadens the barrier and dissipation space without reopening a much larger staged campaign.
- It keeps focus on the knobs most likely to move the system toward switching.
- It compares the two magnetic-feedback prescriptions under identical non-`Hs_update_mode` parameters.
- It preserves the queue structure that has already been validated in this project.

## Suggested targeted baseline

For one paired pilot point under the new convention, use:

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

## Next widening rule

If this 160-point batch still shows no clear threshold crossing, widen only locally:

- add `pulse_amp = 2.5`, or
- open a short pulse-timing strip around the best corner,

but do not jump immediately to a full multi-stage campaign.
