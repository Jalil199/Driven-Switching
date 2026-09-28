# Goal 6

## Objective

Turn the two clean Goal 5 pulse-response mechanisms into a cleaner retained switch.

Goal 5 did not yet prove one robust pulse-induced switching trajectory. It found two useful mechanisms:

- `scan_0063...hsneq`: strong clean inversion during the pulse, but it rebounds to positive `avg_Mz`.
- `scan_0152...hsfull`: weaker in-pulse inversion, but better delayed/post-pulse retention.

Goal 6 should therefore be a focused local follow-up, not a broad global sweep.

## Success Criteria

Primary clean-switch target:

- pre-pulse state remains positive and stable:
  - `pre_pulse_min_avg_Mz > 0`
  - no pre-pulse crossing of `avg_Mz = 0`
- pulse-induced response is strong:
  - `pulse_window_drop_from_start >= 1.2`
  - preferred: `pulse_window_min_avg_Mz <= -0.8`
- retained or stabilized negative state:
  - preferred: `avg_Mz_final <= -0.5`
  - strong target: `avg_Mz_final <= -0.8`
  - fallback useful target: post-pulse minimum `<= -0.8` with slower rebound than Goal 5

Secondary diagnostics:

- monitor `sigma_delta` and full spin-density diagnostics to distinguish transient electronic torque from retained magnetic response.
- keep current-balance diagnostics in the summary tables, but do not rank primarily by current unless numerical instability appears.

## Fixed Protocol Unless Explicitly Varied

- `t_end = 5000.0`
- `pulse_start = 1000.0`
- `pulse_ramp = 5.0`
- `Bx = 0.0`
- `left_bias_scale = 0.0`
- `right_bias_scale = 1.0`

Use `t_end = 5000.0` instead of `4000.0` so delayed or retained switching has a longer post-pulse window.

For any case with `pulse_width = 3000`, use `t_end = 6000.0`. With `pulse_start = 1000`, that widest pulse ends near `t = 4000`, so `t_end = 6000` preserves a meaningful post-pulse relaxation/retention window.

## Branch A: Stabilize the Strong `hsneq` In-Pulse Inversion

Purpose:

- start from the strongest clean in-pulse Goal 5 event.
- test whether longer pulse action and stronger damping prevent the rebound.

Center:

- source candidate: `scan_0063...hsneq`
- `Hs_update_mode = nonequilibrium_only`
- `pulse_amp = 1.5`
- `gamma_so = 0.9`
- `j_sd = 0.3`
- `j_ex = -0.02`
- `j_ani = 0.03`
- `damping = 0.3`
- Goal 5 result:
  - `avg_Mz_at_pulse_start ≈ 0.990`
  - `pulse_window_min_avg_Mz ≈ -0.921`
  - `avg_Mz_final ≈ 0.99995`

Immediate local sweep:

- `pulse_width`: `1000`, `1500`, `2000`, `2500`, `3000`
- `damping`: `0.3`, `0.5`, `0.7`
- `j_ani`: `0.00`, `0.01`, `0.03`, `0.05`
- fixed:
  - `pulse_amp = 1.5`
  - `gamma_so = 0.9`
  - `j_sd = 0.3`
  - `j_ex = -0.02`
  - `Hs_update_mode = nonequilibrium_only`

Count:

- `5 x 3 x 4 = 60` jobs

High-damping retention check:

- run after or alongside the immediate Branch A sweep, but keep it marked separately from the main Cartesian block.
- use `damping = 1.0` only at the most informative Branch A points:
  - all `pulse_width` values at the center `j_ani = 0.03`, or
  - the best two `j_ani` values from the `damping <= 0.7` results if the first sweep is analyzed before launching this strip.
- purpose: test whether stronger damping suppresses rebound without erasing the in-pulse inversion.
- expected count if run as a simple center-`j_ani` strip: `5` jobs.

Optional second strip if Branch A still rebounds:

- take the best `pulse_width`, `damping`, and `j_ani` from the immediate local sweep.
- vary:
  - `pulse_amp`: `1.25`, `1.5`, `1.75`, `2.0`
  - `j_ex`: `-0.03`, `-0.02`, `-0.01`
- count: `4 x 3 = 12` jobs

Interpretation:

- If longer `pulse_width` plus higher `damping` gives negative final `avg_Mz`, then Goal 5 already found the correct torque mechanism and Goal 6 mainly needs retention tuning.
- If all cases still rebound strongly, the `nonequilibrium_only` branch may be a transient excitation mechanism rather than the best retained-switch mechanism.

## Branch B: Strengthen the Delayed `hsfull` Retention Candidate

Purpose:

- start from the Goal 5 case that retained a negative final state better.
- test whether the delayed reversal can be pulled into the pulse window or made more strongly retained.

Center:

- source candidate: `scan_0152...hsfull`
- `Hs_update_mode = full_spin_density`
- `pulse_amp = 2.0`
- `gamma_so = 0.9`
- `j_sd = 0.3`
- `j_ex = -0.01`
- `j_ani = 0.0`
- `damping = 0.5`
- Goal 5 result:
  - `avg_Mz_at_pulse_start ≈ 0.820`
  - `pulse_window_min_avg_Mz ≈ -0.134`
  - `post_pulse_min_avg_Mz ≈ -0.851`
  - `avg_Mz_final ≈ -0.508`

Immediate local sweep:

- `pulse_width`: `1500`, `2000`, `2500`, `3000`
- `pulse_amp`: `1.75`, `2.0`, `2.25`
- `damping`: `0.4`, `0.5`, `0.7`, `1.0`
- fixed:
  - `gamma_so = 0.9`
  - `j_sd = 0.3`
  - `j_ex = -0.01`
  - `j_ani = 0.0`
  - `Hs_update_mode = full_spin_density`

Count:

- `4 x 3 x 4 = 48` jobs

Do not include `pulse_amp = 2.5` in the first Branch B batch. First test whether `2.25` strengthens retention or moves the delayed minimum into the pulse window. Add `2.5` only if the `1.75-2.25` strip shows monotonic improvement without destabilizing diagnostics.

Optional second strip if Branch B improves but remains delayed:

- take the best `pulse_width`, `pulse_amp`, and `damping` from the immediate local sweep.
- vary:
  - `j_ex`: `-0.02`, `-0.01`, `0.0`
  - `j_ani`: `0.0`, `0.01`, `0.02`, `0.03`
- count: `3 x 4 = 12` jobs

Interpretation:

- If longer or stronger pulses move the minimum from post-pulse into the pulse window, this branch becomes the best retained-switch candidate.
- If it keeps switching only after the pulse, then Goal 6 should treat it as a delayed relaxation pathway and tune hold/turnoff rather than only amplitude.

## Branch C: Direct A/B Check at the Two Centers

Purpose:

- verify that the previous `Hs_update_mode` preference was not caused by only one parameter corner.

Run the two center points with both update modes:

- `scan_0063` center:
  - `Hs_update_mode = nonequilibrium_only`
  - `Hs_update_mode = full_spin_density`
- `scan_0152` center:
  - `Hs_update_mode = nonequilibrium_only`
  - `Hs_update_mode = full_spin_density`

Also repeat the best result from Branch A and Branch B with the opposite `Hs_update_mode`.

Count:

- `4` center checks
- plus `2` best-opposite-mode checks after Branch A/B finish

Interpretation:

- If the mode preference survives the A/B check, keep separate `hsneq` and `hsfull` branches.
- If one mode dominates both centers, collapse later Goal 6 follow-up into that mode.

## Recommended First Goal 6 Batch

Run one bounded batch with:

- Branch A immediate local sweep: `60` jobs
- Branch B immediate local sweep: `48` jobs
- Branch C center A/B checks: `4` jobs

Total first batch:

- `112` jobs

Recommended add-on if compute budget allows:

- Branch A high-damping center-`j_ani` strip: `5` jobs
- total with this add-on: `117` jobs

The `112`-job batch remains the clean first batch. The `5` high-damping jobs are useful but should be labeled separately so they do not blur the main Branch A sweep statistics.

Execution policy:

- use the same queued workflow as Goal 5.
- cap simultaneous jobs at `34`.
- do not launch the optional second strips until the first batch is analyzed.

## Ranking Table for Goal 6

Do not collapse all Goal 6 outcomes into one score. Produce two primary leaderboards plus one diagnostic crosscheck table.

Leaderboard 1: in-pulse switching

1. clean pre-pulse behavior.
2. in-pulse strong inversion, `pulse_window_min_avg_Mz <= -0.8`.
3. large in-pulse drop, `pulse_window_drop_from_start >= 1.2`.
4. lower rebound after the in-pulse minimum.
5. retained final inversion as a secondary tie-breaker.
6. acceptable numerical diagnostics.

Leaderboard 2: retention and delayed switching

1. clean pre-pulse behavior.
2. retained final inversion, especially `avg_Mz_final <= -0.5` and then `<= -0.8`.
3. delayed post-pulse inversion, `post_pulse_min_avg_Mz <= -0.8`, if final state is not yet retained.
4. earlier crossing time, preferring pulse-window crossing over late post-pulse crossing.
5. lower rebound after the minimum.
6. acceptable numerical diagnostics.

Diagnostic crosscheck:

- direct `Hs_update_mode` paired comparisons at the two centers and at the best Branch A/B cases.
- report both raw values and deltas for `pulse_window_min_avg_Mz`, `post_pulse_min_avg_Mz`, `avg_Mz_final`, `sigma_delta`, and current-balance metrics.

Recommended output tables:

- `goal6_in_pulse_ranked_runs.tsv`
- `goal6_retention_ranked_runs.tsv`
- `goal6_branch_summary.tsv`
- `goal6_best_retention_candidates.tsv`
- `goal6_hs_mode_crosscheck.tsv`

Recommended figures:

- top Branch A trajectories
- top Branch B trajectories
- direct `Hs_update_mode` crosscheck trajectories
- pulse-window and post-pulse minimum maps

## Discussion Points Before Launch

Settled choices for the next launcher draft:

- Use `t_end = 6000` for `pulse_width = 3000` cases.
- Include Branch A `damping = 1.0` only as a labeled high-damping strip, not as part of the main `60`-job Cartesian block.
- Do not include Branch B `pulse_amp = 2.5` in the first batch.
- Use separate in-pulse and retention leaderboards.

Remaining questions:

- Should we add a controlled pulse turnoff or hold parameter if the code supports it?
- Should Branch C also repeat the high-damping Branch A winner with the opposite `Hs_update_mode`, or wait until the base Branch A/B winners are known?

## Current Recommendation

Start with the `112`-job first Goal 6 batch, plus the optional `5`-job Branch A high-damping strip if compute budget is acceptable.

Do not widen to a global scan yet. The immediate question is whether the two Goal 5 mechanisms can be combined:

- Branch A asks whether the strong `hsneq` transient can be stabilized.
- Branch B asks whether the retained `hsfull` response can be strengthened and shifted earlier.
