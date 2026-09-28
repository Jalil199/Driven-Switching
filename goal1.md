Your first task is NOT to launch a large parameter sweep yet.

Start by building a reliable baseline workflow from my existing TDNEGF example:
`examples/notebooks/02_spin_dynamics_and_electrons.ipynb`

Phase 1 — Understand the existing example

1. Inspect the notebook carefully and identify:

   * the exact model geometry,
   * the square-lattice setup,
   * the two leads,
   * the SOC terms,
   * the electron pulse definition,
   * the classical-spin dynamics parameters,
   * the observables currently computed and plotted.
2. Trace the notebook to the package source code and determine which functions are actually responsible for:

   * model construction,
   * time propagation,
   * self-energies / lead coupling,
   * classical-spin updates,
   * observable extraction.
3. Write a short structured note explaining the current workflow before changing anything.

Phase 2 — Reproduce one baseline run
4. Reproduce one run as close as possible to the notebook example.
5. Confirm that the simulation completes successfully.
6. Save the outputs in a clean results folder with:

* params.json
* summary.json
* observables file
* log.txt

7. Report whether the run appears physically sensible, including a brief comment on:

   * charge/current behavior,
   * spin dynamics behavior,
   * whether any apparent switching or partial reversal is seen,
   * any numerical instability or suspicious feature.

Phase 3 — Prepare a sweep-ready driver
8. Convert the notebook workflow into a single-run driver script that can accept external parameters.
9. Keep the code changes minimal and localized.
10. Do not perform a broad sweep yet.
11. Only after the baseline run is validated, propose a small pilot sweep over 2–3 key parameters.

Important constraints

* Reuse the current TDNEGF workflow already present in the repository.
* Do not refactor unrelated code.
* Do not invent new physics not already supported by the example unless clearly marked.
* Prefer reproducibility and clean outputs over ambitious exploration.
* At this stage, correctness is more important than breadth.

