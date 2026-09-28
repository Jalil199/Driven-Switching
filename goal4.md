Good. Keep the existing validated workflow and launch a queued parameter scan.

Use these parameters for the scan:

* pulse_amp = [0.5, 0.75, 1.0, 1.25]
* gamma_so = [0.3, 0.5, 0.7]
* j_sd = [0.05, 0.1, 0.2]
* j_ex = [-0.02, -0.05, -0.1, -0.2]
* Bx = [0.0, 1e-4]

Keep these fixed for now:

* j_ani = 0.1
* Nx = 2
* Ny = 2

Use these time settings:

* first check the baseline pre-pulse behavior of the nonequilibrium-minus-equilibrium electronic density / spin-density contribution that drives the LLG update
* verify carefully how the equilibrium density / spin density is computed in the reference example, how the nonequilibrium contribution is computed during time evolution, and confirm that the LLG is driven by the intended nonequilibrium-minus-equilibrium contribution following the same convention already used in the notebook
* if that check is reasonable, use `t_pulse_start = 1000.0`
* use `pulse_duration = 2000.0`
* use `t_total = 5000.0`

This should provide:

* a long pre-pulse settling window,
* a long pulse window,
* and a long post-pulse window to determine whether switching or partial reversal persists.

Also:

1. Keep the pulse delay fixed for the scan unless the baseline check clearly indicates it should be adjusted.
2. Use Julia `Distributed`.
3. Run the jobs with a queued workflow and a maximum of 30 threads / workers in parallel, so that when one run finishes, the next pending one starts automatically.
4. Keep the maximum parallel load capped at 30 and do not launch more than that concurrently.
5. Keep the simulation logic in the single-run driver and use the sweep launcher only to generate parameter combinations, assign `run_id`, dispatch jobs, and collect summaries.
6. Assign a unique `run_id` to every case.
7. Save one output folder per run.
8. Save a global summary table with one row per run.
9. Do not overwrite previous validated results.
10. Save all important observables and metadata that we will later need to analyze and plot in the notebook created earlier.
11. At minimum, save the parameter values, time grid, relevant charge/spin observables, magnetization dynamics, final-state summaries, run status, and any useful diagnostics for identifying switching, partial reversal, or suspicious runs.
12. Report clearly which files were added or modified.

