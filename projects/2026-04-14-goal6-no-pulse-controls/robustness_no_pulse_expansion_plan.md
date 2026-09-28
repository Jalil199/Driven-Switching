# Goal 6 Robustness Expansion Plan

## Objective

Test whether the current conclusion remains true beyond the two already matched pairs:

- pulse runs switch in `avg_Mz(t)`,
- matched `no-pulse` controls do not switch,
- raw `max_i |sigma_delta,i(t)|` is not a pulse-only signature.

The goal here is robustness, not mechanism. The main question is whether the magnetization conclusion survives across distinct Goal 6 switching families, not just `goal6_A_main_0003` and `goal6_A_main_0034`.

## Current baseline

Already confirmed:

- `goal6_A_main_0003` with pulse ends near `avg_Mz_final=-0.99974`; its matched `no-pulse` control ends near `+0.99040`.
- `goal6_A_main_0034` with pulse ends near `avg_Mz_final=-0.99987`; its matched `no-pulse` control ends near `+0.99710`.
- In both pairs, raw `max_sigma_delta_abs` is essentially unchanged by removing the pulse.

So the next useful test is breadth: do more winners still stay positive without pulse?

## Selection strategy

Do not spend the next batch on many near-duplicate widths from the same family. Prefer one representative from each scientifically distinct switching family:

1. strong in-pulse retained switching,
2. delayed retained switching,
3. different damping strip,
4. different `Hs_update_mode` branch.

This gives better falsification power per run than sampling many width clones inside one family.

## Recommended phases

### Phase 1: four new matched `no-pulse` controls

These are the highest-value additions because they extend the claim across distinct families rather than only adding more widths inside already-covered families.

1. `goal6_A_main_0010_pa1p5_gso0p9_jsd0p3_jexm0p02_jani0p01_damp0p7_twidth1000p0_hsneq`
Reason: same delayed-retention family as `0034`, but with shorter pulse width; checks whether the positive no-pulse result is family-wide or only specific to the 2000-width case.

2. `goal6_A_main_0059_pa1p5_gso0p9_jsd0p3_jexm0p02_jani0p03_damp0p7_twidth3000p0_hsneq`
Reason: distinct retained-switch family with `j_ani=0.03`, `damping=0.7`, and late switching; extends robustness away from the already-tested `damping=0.3` in-pulse winner.

3. `goal6_B_main_0042_pa2p0_gso0p9_jsd0p3_jexm0p01_jani0p0_damp0p5_twidth3000p0_hsfull`
Reason: strongest retained result in Branch B and the only top-ranked case using `Hs_update_mode=full_spin_density`; critical branch-level cross-check.

4. `goal6_A_high_damping_0001_pa1p5_gso0p9_jsd0p3_jexm0p02_jani0p03_damp1p0_twidth1000p0_hsneq`
Reason: probes the high-damping strip, which is physically distinct from the current `damping=0.3/0.7` baseline.

### Phase 2: only if Phase 1 is clean

If all four Phase 1 controls stay clearly positive, add width- or damping-neighbor controls to tighten the statement inside already-promising families:

1. `goal6_A_main_0011_pa1p5_gso0p9_jsd0p3_jexm0p02_jani0p03_damp0p7_twidth1000p0_hsneq`
2. `goal6_A_main_0035_pa1p5_gso0p9_jsd0p3_jexm0p02_jani0p03_damp0p7_twidth2000p0_hsneq`
3. `goal6_A_main_0030_pa1p5_gso0p9_jsd0p3_jexm0p02_jani0p01_damp0p5_twidth2000p0_hsneq`

These add intra-family confidence, but they are lower priority than Phase 1 because they are more redundant scientifically.

### Phase 3: only if something fails or looks ambiguous

Escalate only if at least one matched `no-pulse` control shows spontaneous drift toward reversal.

Escalation targets:

1. replicate the surprising family with one neighboring width,
2. inspect `avg_Mz(t)` and `sigma_delta(t)` timing around the spontaneous drift,
3. reconsider initialization or explicit pre-thermalization only for that failing family.

## Success and escalation criteria

Treat Phase 1 as supportive of robustness if every new matched control satisfies all of:

- `avg_Mz_final > +0.9`,
- `avg_Mz_min > 0`,
- no zero crossing in `avg_Mz(t)`,
- pulse twin still ends with strong negative retention.

Treat the robustness claim as weakened if any control shows one of:

- `avg_Mz_final < +0.5`,
- `avg_Mz_min < 0`,
- any sustained post-pulse drift toward negative `avg_Mz`,
- strong family dependence where only one parameter strip remains positive without pulse.

Escalate immediately if any control ends near or below zero. That would mean the present “pulse is required for switching” statement is not yet general.

## Metrics to compare for each new pair

Primary metrics:

- `avg_Mz_final`
- `avg_Mz_min`
- `first_cross_zero_time`
- whether `avg_Mz(t)` stays positive in the control

Secondary metrics:

- `max_sigma_delta_abs`
- `pre_pulse_sigma_delta_absmax`
- `post_pulse_excess_sigma_delta_absmax`

For the robustness question, magnetization metrics are decisive. `sigma_delta` remains secondary here.

## Minimal implementation path

Reuse the current project instead of creating a separate workflow.

1. Extend `run_goal6_no_pulse_controls.jl` so the case list is phase-based or manifest-driven rather than hard-coded to two entries.
2. Extend `analyze_goal6_no_pulse_controls.py` so `CASE_MAP` is generated from the manifest or a small JSON/TSV pairing file.
3. Write grouped outputs under the same root:
   - `results/matched_no_pulse_controls/scan_summary.tsv`
   - `results/matched_no_pulse_controls/analysis/pulse_minus_no_pulse_summary.tsv`
   - one updated `avg_Mz(t)` comparison figure
   - one ranked control-outcome table
4. Add a compact review note summarizing whether each family stayed positive without pulse.

## Review questions for later

1. Do all representative switching families still fail to switch without pulse?
2. Is the current conclusion valid only for Branch A `hsneq`, or also for Branch B `hsfull` and the high-damping strip?
3. Are there any families where the pulse seems unnecessary for magnetization reversal?
4. After Phase 1, is the next highest-value step still robustness expansion, or does one anomalous family force a mechanism/preparation study first?
