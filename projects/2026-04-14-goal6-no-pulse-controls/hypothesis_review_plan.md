# Goal 6 Hypothesis Review Plan

## Current hypothesis

The raw electronic scale `max_i |sigma_delta,i(t)|` in the Goal 6 winners is not by itself a pulse-specific switching indicator, because a large part of that signal may already be present in the disconnected-start transient. The pulse should therefore be interpreted through the matched difference `pulse - no-pulse`, not through the raw peak alone.

## Evidence already available

- `goal6_A_main_0003` with pulse:
  - `avg_Mz_final = -0.9997414118745744`
  - `max_sigma_delta_abs = 0.17943600247824956`
- `goal6_A_main_0003` matched `no-pulse` control:
  - `avg_Mz_final = 0.990400993573977`
  - `avg_Mz_min = 0.9835380764267725`
  - `max_sigma_delta_abs = 0.17943600247824956`
- Interim interpretation:
  - the magnetization outcome changes dramatically with the pulse,
  - but the raw `sigma_delta` peak does not,
  - so the hypothesis is partially supported already: the raw peak is dominated by a contribution that is not sufficient to explain switching.

## Pending decisive check

- Wait for `goal6_A_main_0034_nopulse_gso0p9_jsd0p3_jexm0p02_jani0p01_damp0p7_twidth2000p0_hsneq` to finish.
- Then run `analyze_goal6_no_pulse_controls.py` and inspect:
  - `pulse_sigma_delta_absmax_peak`
  - `control_sigma_delta_absmax_peak`
  - `excess_sigma_delta_absmax_peak`
  - `pre_pulse_excess_sigma_delta_absmax_peak`
  - `post_pulse_excess_sigma_delta_absmax_peak`
  - `pulse_avg_Mz_final`
  - `control_avg_Mz_final`
  - `excess_avg_Mz_final`

## Decision rules for later review

- If both controls stay near `avg_Mz_final > 0.9` while both pulse runs end near `avg_Mz_final < -0.99`, then switching itself is robust without needing pre-thermalization to prove existence.
- If `control_sigma_delta_absmax_peak` stays close to the pulse peak in both matched pairs, treat raw `sigma_delta` amplitude as a startup-transient observable rather than a pulse-only signature.
- If `post_pulse_excess_sigma_delta_absmax_peak` is small compared with the raw pulse peak, do not use raw `sigma_delta` magnitude as the mechanism claim.
- If `post_pulse_excess_sigma_delta_absmax_peak` is clearly isolated after `t = 1000`, then the pulse still leaves a measurable electronic excess even if the raw peak is contaminated by startup dynamics.
- If the `0034` no-pulse control drifts negative or close to zero by itself, escalate pre-thermalization from "nice-to-have for attribution" to "needed before trusting the mechanism."

## Next review checklist

- [ ] Confirm `projects/2026-04-14-goal6-no-pulse-controls/run_status.json` reaches `completed_cases = 2`.
- [ ] Confirm `summary.json` and `observables.tsv` exist for the `0034` no-pulse control.
- [ ] Confirm the finalizer produced `results/matched_no_pulse_controls/analysis/pulse_minus_no_pulse_summary.json`.
- [ ] Read the paired metrics for `0003` and `0034`.
- [ ] Decide whether pre-thermalization is unnecessary, recommended for attribution only, or required before further mechanism claims.
