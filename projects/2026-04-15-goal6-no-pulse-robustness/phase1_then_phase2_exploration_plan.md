# Goal 6 Follow-Up Exploration Plan: Step 1 Then Step 2

## Scope

This plan keeps the next exploration split into two sequential questions:

1. robust generality of the `no-pulse` conclusion,
2. pulse-specific electronic mechanism after matched subtraction.

Do not mix them in one batch. Step 2 depends on Step 1 because the mechanism readout is more convincing only after the `no-pulse` baseline is shown to stay positive across multiple families.

## Current state

The active batch in `projects/2026-04-15-goal6-no-pulse-robustness/run_status.json` is already executing Step 1 Phase 1:

- `goal6_A_main_0010`
- `goal6_A_main_0059`
- `goal6_B_main_0042`
- `goal6_A_high_damping_0001`

All four are matched `bias_amplitude=0.0` controls for previously identified Goal 6 pulse-on switching winners.

## Step 1: robustness expansion

### Scientific question

Does the statement "the pulse is required for switching in `avg_Mz`" survive across distinct Goal 6 switching families rather than only the two already matched pairs (`0003`, `0034`)?

### Immediate execution plan

1. Let the active four-case batch finish.
2. Run the robustness postprocessor on all completed controls.
3. Build one compact review table across the six matched pairs available after completion:
   - existing pairs: `0003`, `0034`
   - active new pairs: `0010`, `0059`, `0042`, `A_high_damping_0001`
4. Decide whether robustness is now strong enough to stop, or whether a narrow Phase 2-of-Step-1 extension is needed.

### Decision metrics

For each matched `no-pulse` control:

- `avg_Mz_final`
- `avg_Mz_min`
- `first_cross_zero_time`
- whether `avg_Mz(t)` ever becomes negative

### Success gate to move on to Step 2

Advance directly to Step 2 if every new control satisfies:

- `avg_Mz_final > +0.9`
- `avg_Mz_min > 0`
- no zero crossing in `avg_Mz(t)`
- its pulse twin still shows strong retained negative final state

### Escalation inside Step 1

Do not move to Step 2 immediately if any control shows ambiguity. Instead do a narrow robustness extension first.

Trigger escalation if any control shows:

- `avg_Mz_final < +0.9`
- `avg_Mz_min <= 0`
- any zero crossing
- clear late drift toward zero or negative `avg_Mz`

### Step 1 Phase 2 only if needed

If the active four-case batch is mostly clean but one family remains too sparse or ambiguous, add only family-neighbor controls rather than opening a broad new scan. Priority neighbors:

1. `goal6_A_main_0011`
2. `goal6_A_main_0035`
3. `goal6_A_main_0030`

These are for tightening family-local robustness, not for exploratory mechanism work.

## Step 2: pulse-specific `sigma_delta` mechanism mapping

### Entry condition

Start Step 2 only after Step 1 is either:

- cleanly supportive across the representative families, or
- narrowed to a clearly defined anomalous family that we intentionally want to dissect.

### Scientific question

After subtracting the matched `no-pulse` baseline, what electronic response is actually induced by the pulse, and when does it appear relative to `pulse_start` and the later `avg_Mz` reversal?

### Dataset to analyze first

Use the matched pairs with the cleanest robustness evidence:

- mandatory baseline pairs: `0003`, `0034`
- plus every Step 1 family that passes the positivity gate

### Required outputs

For each matched pair:

1. `pulse - no-pulse` traces for every available `sigma_delta_*` component.
2. Site-resolved and component-resolved maps around the pulse window.
3. Zoomed curves around `pulse_start`.
4. Timing markers for:
   - `pulse_start`
   - peak excess `sigma_delta`
   - first strong drop in `avg_Mz`
   - zero crossing in `avg_Mz`, if present in the pulse twin

### Minimal postprocessing products

- one summary TSV with peak times and amplitudes of the excess signal
- one figure set for per-pair `sigma_delta` component traces
- one heatmap set showing excess signal by site/component/time
- one short interpretation note distinguishing startup transient from pulse-induced excess

### Step 2 decision target

The goal is not to re-prove switching. The goal is to identify whether the pulse leaves:

- a localized excess just after `pulse_start`,
- a delayed excess preceding magnetization reversal,
- or no strong unique electronic signature beyond what the startup transient already generates.

## Recommended order from here

1. Finish and review the active Step 1 robustness batch.
2. If all four new controls stay positive, write the six-pair robustness summary and stop Step 1 there.
3. Then launch Step 2 postprocessing on the clean matched pairs only.
4. Only after Step 2 decide whether an initialization or pre-thermalization stress test is scientifically necessary.

## Explicit non-goals for this turn

- No new pulse-on switching scan.
- No broader parameter exploration.
- No initialization variant yet.

Those are downstream tasks and should not be mixed into the current Step 1 to Step 2 sequence.
