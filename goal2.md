Good. The baseline workflow is now validated. Move one step forward only.

Next task:

1. Finalize a clean Julia single-run driver that executes one simulation from external parameters.
2. Create a separate Julia sweep launcher that defines multiple parameter sets and runs independent cases in parallel.
3. Use Julia `Distributed` for the parallel execution of independent runs.
4. Start with only a small pilot sweep over 2–3 key switching parameters.
5. Assign a unique `run_id` to every case.
6. Save all outputs in a dedicated results directory, with one subfolder per run.
7. Produce a global summary table with one row per run.
8. Reuse the validated baseline workflow and keep all code changes minimal and localized.

Implementation recommendation:

* Create one Julia file for a single run, for example `run_switching_case.jl`, it can be one for FM and AFM initial configuration.
* Create a separate Julia file for the pilot sweep, for example `run_pilot_sweep.jl`.
* Put the simulation logic in the single-run driver.
* Use the sweep file only to generate parameter combinations, assign `run_id`, dispatch jobs with `Distributed`, and collect summaries.
* Avoid duplicating physics code across files.


