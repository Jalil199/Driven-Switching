"""Generate a Julia Jupyter notebook for queued parameter-scan analysis."""

from __future__ import annotations

from pathlib import Path
import textwrap

import nbformat as nbf


PROJECT_DIR = Path(__file__).resolve().parent
NOTEBOOK_PATH = PROJECT_DIR / "queued_scan_analysis.ipynb"


def markdown_cell(source: str):
    return nbf.v4.new_markdown_cell(textwrap.dedent(source).strip() + "\n")


def code_cell(source: str):
    return nbf.v4.new_code_cell(textwrap.dedent(source).strip() + "\n")


def build_notebook():
    nb = nbf.v4.new_notebook()
    nb["metadata"].update(
        {
            "kernelspec": {
                "display_name": "Julia 1.12.4",
                "language": "julia",
                "name": "julia-1.12",
            },
            "language_info": {
                "name": "julia",
                "version": "1.12.4",
                "file_extension": ".jl",
                "mimetype": "application/julia",
            },
        }
    )

    cells = [
        markdown_cell(
            """
            # Queued Scan Analysis

            This notebook is intentionally written in **Julia** and uses
            **PyPlot** for all figures. It loads the latest
            `results/queued_parameter_scan_*` bundle that already contains
            `scan_summary.tsv`, then builds a compact review layer focused on:

            - completion and data-quality overview,
            - conservative ranking of the most responsive parameter points,
            - pre-pulse versus full-window diagnostics,
            - compact parameter-to-observable maps with the requested plot style,
            - time-domain views of `avg_Mz(t)` together with the applied pulse.
            """
        ),
        code_cell(
            r"""
            import Pkg

            ENV["MPLBACKEND"] = get(ENV, "MPLBACKEND", "Agg")

            required_packages = ["CSV", "DataFrames", "JSON3", "PyPlot"]
            missing_packages = [pkg for pkg in required_packages if Base.find_package(pkg) === nothing]

            if !isempty(missing_packages)
                println("Installing missing Julia packages: $(join(missing_packages, ", "))")
                Pkg.add(missing_packages)
            end
            """
        ),
        code_cell(
            r"""
            using CSV
            using DataFrames
            using JSON3
            using Markdown
            using Printf
            using PyPlot
            using Statistics
            """
        ),
        code_cell(
            r"""
            function configure_pyplot!()
                latex_available = Sys.which("latex") !== nothing

                PyPlot.rc("text", usetex=latex_available)
                PyPlot.rc("font", family="serif")
                PyPlot.rc("font", size=18)
                PyPlot.rc("axes", labelsize=22, titlesize=22, linewidth=2)
                PyPlot.rc("xtick", labelsize=18)
                PyPlot.rc("ytick", labelsize=18)
                PyPlot.rc("legend", fontsize=16)

                if !latex_available
                    println("LaTeX was not found in this environment. Using the requested serif style with usetex=False.")
                end

                return latex_available
            end

            function collect_ready_roots(results_dir::AbstractString)
                roots = sort([
                    path for path in readdir(results_dir; join=true)
                    if isdir(path) && startswith(basename(path), "queued_parameter_scan_")
                ])
                return [root for root in roots if isfile(joinpath(root, "scan_summary.tsv"))]
            end

            function pulse_envelope_value(t::Real; pulse_start::Real, pulse_ramp::Real, pulse_width::Real, bias_amplitude::Real)
                τ = Float64(t) - Float64(pulse_start)
                ti = Float64(pulse_ramp)
                Δ = Float64(pulse_width)

                if τ < 0.0
                    return 0.0
                elseif τ < ti
                    return sin((pi / 2) * τ / ti)^2 * Float64(bias_amplitude)
                elseif τ < Δ
                    return Float64(bias_amplitude)
                elseif τ < Δ + ti
                    return cos((pi / 2) * (τ - Δ) / ti)^2 * Float64(bias_amplitude)
                else
                    return 0.0
                end
            end

            pulse_trace(times, row) = [
                pulse_envelope_value(
                    t;
                    pulse_start=row.pulse_start,
                    pulse_ramp=row.pulse_ramp,
                    pulse_width=row.pulse_width,
                    bias_amplitude=row.bias_amplitude,
                ) for t in times
            ]

            function flatten_axes(axs)
                if axs isa AbstractArray
                    flattened = Any[]
                    for ax in axs
                        append!(flattened, flatten_axes(ax))
                    end
                    return flattened
                end

                try
                    collected = collect(axs)
                    isempty(collected) && return [axs]

                    flattened = Any[]
                    for ax in collected
                        append!(flattened, flatten_axes(ax))
                    end
                    return flattened
                catch
                    return [axs]
                end
            end

            function mean_over_columns(df::DataFrame, cols)
                isempty(cols) && return fill(NaN, nrow(df))
                normalized_cols = Symbol.(cols)
                columns = [Float64.(df[!, col]) for col in normalized_cols]
                return [mean(values) for values in zip(columns...)]
            end

            function sum_over_columns(df::DataFrame, cols)
                isempty(cols) && return fill(NaN, nrow(df))
                normalized_cols = Symbol.(cols)
                columns = [Float64.(df[!, col]) for col in normalized_cols]
                return [sum(values) for values in zip(columns...)]
            end

            function site_index(name)
                match_obj = match(r"site(\d+)$", String(name))
                match_obj === nothing && return typemax(Int)
                return parse(Int, match_obj.captures[1])
            end

            function site_columns(df::DataFrame, prefix::AbstractString)
                cols = filter(name -> startswith(String(name), prefix), names(df))
                return sort(cols; by=site_index)
            end

            function columns_to_matrix(df::DataFrame, cols)
                isempty(cols) && return Array{Float64}(undef, nrow(df), 0)
                normalized_cols = Symbol.(cols)
                matrix = hcat([Float64.(df[!, col]) for col in normalized_cols]...)
                return matrix
            end

            function late_prepulse_mask(times, pulse_start; fraction=0.2)
                pre_mask = times .< pulse_start
                pre_idx = findall(pre_mask)
                isempty(pre_idx) && return falses(length(times))
                n_keep = max(1, ceil(Int, fraction * length(pre_idx)))
                keep_idx = pre_idx[(end - n_keep + 1):end]
                mask = falses(length(times))
                mask[keep_idx] .= true
                return mask
            end

            function masked_mean(values, mask)
                idx = findall(mask)
                isempty(idx) && return NaN
                selected = [Float64(values[i]) for i in idx if isfinite(Float64(values[i]))]
                isempty(selected) && return NaN
                return mean(selected)
            end

            function columnwise_masked_mean(matrix, mask)
                size(matrix, 2) == 0 && return Float64[]
                idx = findall(mask)
                isempty(idx) && return fill(NaN, size(matrix, 2))
                means = Float64[]
                for col in Base.axes(matrix, 2)
                    vals = [Float64(matrix[i, col]) for i in idx if isfinite(Float64(matrix[i, col]))]
                    push!(means, isempty(vals) ? NaN : mean(vals))
                end
                return means
            end

            function load_observables(run_row)
                obs_path = joinpath(String(run_row.results_dir), "observables.tsv")
                isfile(obs_path) || error("Missing observables.tsv for $(run_row.run_id): $(obs_path)")
                return DataFrame(CSV.File(obs_path; delim='\t', silencewarnings=true))
            end

            latex_available = configure_pyplot!()
            """
        ),
        code_cell(
            r"""
            project_dir = pwd()
            results_dir = joinpath(project_dir, "results")
            ready_roots = collect_ready_roots(results_dir)

            isempty(ready_roots) && error("No queued_parameter_scan_* root with scan_summary.tsv was found.")

            scan_root = ready_roots[end]
            summary_path = joinpath(scan_root, "scan_summary.tsv")
            manifest_path = joinpath(scan_root, "case_manifest.tsv")
            config_path = joinpath(scan_root, "scan_config.json")

            summary_df = DataFrame(CSV.File(summary_path; delim='\t', silencewarnings=true))
            manifest_df = DataFrame(CSV.File(manifest_path; delim='\t', silencewarnings=true))
            config = JSON3.read(read(config_path, String))

            println("Using scan root: $(scan_root)")
            println("Rows in summary: $(nrow(summary_df)) / requested: $(config.requested_case_count)")
            """
        ),
        code_cell(
            r"""
            function bool_label(full_reversal, partial_reversal)
                if full_reversal === true
                    return "full_reversal"
                elseif partial_reversal === true
                    return "partial_reversal"
                else
                    return "no_reversal_or_pending"
                end
            end

            function suspicious_flag(row)
                charge_bad =
                    hasproperty(row, :max_site_charge_deviation_from_1) &&
                    !ismissing(row.max_site_charge_deviation_from_1) &&
                    row.max_site_charge_deviation_from_1 > 0.25
                balance_bad =
                    hasproperty(row, :current_balance_peak_abs) &&
                    !ismissing(row.current_balance_peak_abs) &&
                    row.current_balance_peak_abs > 0.25
                return charge_bad || balance_bad
            end

            function safe_fraction(num, den)
                if ismissing(num) || ismissing(den) || den == 0
                    return NaN
                end
                return Float64(num) / Float64(den)
            end

            function safe_value(x, default)
                return ismissing(x) ? default : Float64(x)
            end

            summary_df.suspicious = [suspicious_flag(row) for row in eachrow(summary_df)]
            summary_df.switching_label = [
                bool_label(row.full_reversal_detected, row.partial_reversal_flag) for row in eachrow(summary_df)
            ]
            summary_df.pre_pulse_fraction_of_peak = [
                safe_fraction(row.pre_pulse_sigma_delta_absmax, row.max_sigma_delta_abs) for row in eachrow(summary_df)
            ]
            summary_df.pulse_driven_score = [
                (1.0 - safe_value(row.avg_Mz_min, 1.0)) * 5.0 +
                safe_value(row.neel_z_peak_abs, 0.0) * 25.0 +
                safe_value(row.max_sigma_delta_abs, 0.0) * 0.05 -
                (isnan(row.pre_pulse_fraction_of_peak) ? 0.0 : row.pre_pulse_fraction_of_peak)
                for row in eachrow(summary_df)
            ]

            completed_df = filter(:status => ==("done"), summary_df)
            failed_df = filter(:status => ==("failed"), summary_df)
            pending_df = filter(row -> row.status ∉ ("done", "failed"), summary_df)

            status_counts = combine(groupby(summary_df, :status), nrow => :count)
            sort!(status_counts, :count, rev=true)

            display(Markdown.parse("## Status counts"))
            status_counts
            """
        ),
        code_cell(
            r"""
            candidate_cols = [
                :run_id, :pulse_amp, :gamma_so, :j_sd, :j_ex, :Bx,
                :avg_Mz_min, :avg_Mz_final, :neel_z_peak_abs,
                :max_sigma_delta_abs, :pre_pulse_sigma_delta_absmax,
                :pre_pulse_fraction_of_peak, :pulse_driven_score, :suspicious,
            ]

            ranked_completed = sort(completed_df, [:pulse_driven_score, :avg_Mz_min], rev=[true, false])

            display(Markdown.parse("## Highest-ranked completed candidates"))
            first(select(ranked_completed, candidate_cols), min(10, nrow(ranked_completed)))

            suspicious_cols = [
                :run_id, :status, :pulse_amp, :gamma_so, :j_sd, :j_ex, :Bx,
                :max_site_charge_deviation_from_1, :current_balance_peak_abs, :suspicious,
            ]

            display(Markdown.parse("## Runs flagged as suspicious"))
            first(select(filter(:suspicious => identity, summary_df), suspicious_cols), 20)
            """
        ),
        code_cell(
            r"""
            trace_candidates = filter(row -> isfile(joinpath(String(row.results_dir), "observables.tsv")), ranked_completed)
            n_trace_cases = min(4, nrow(trace_candidates))
            trace_cases = trace_candidates[1:n_trace_cases, :]

            trace_case_cols = [
                :run_id, :pulse_amp, :gamma_so, :j_sd, :j_ex, :Bx,
                :avg_Mz_min, :avg_Mz_final, :max_sigma_delta_abs, :pulse_driven_score,
            ]

            display(Markdown.parse("## Representative completed runs for time traces"))
            select(trace_cases, trace_case_cols)
            """
        ),
        code_cell(
            r"""
            function finite_extrema(values)
                finite_vals = [Float64(v) for v in values if !ismissing(v) && isfinite(Float64(v))]
                isempty(finite_vals) && return (0.0, 1.0)
                return (minimum(finite_vals), maximum(finite_vals))
            end

            function build_heatmap_matrix(df::DataFrame, rowcol::Symbol, colcol::Symbol, valcol::Symbol)
                rows = sort(unique(collect(skipmissing(df[!, rowcol]))))
                cols = sort(unique(collect(skipmissing(df[!, colcol]))))
                mat = fill(NaN, length(rows), length(cols))
                row_index = Dict(v => i for (i, v) in enumerate(rows))
                col_index = Dict(v => j for (j, v) in enumerate(cols))

                for row in eachrow(df)
                    value = row[valcol]
                    if !(ismissing(value) || ismissing(row[rowcol]) || ismissing(row[colcol]))
                        mat[row_index[row[rowcol]], col_index[row[colcol]]] = Float64(value)
                    end
                end

                return rows, cols, mat
            end

            isempty(completed_df) && error("No completed runs are available yet.")

            fig, axes = subplots(1, 3, figsize=(22, 6), constrained_layout=true)

            axes[1].bar(string.(status_counts.status), status_counts.count, color="#1f4e79", alpha=0.9)
            axes[1].set_title("Run status counts")
            axes[1].tick_params(axis="x", rotation=30)
            axes[1].grid(axis="y", alpha=0.25)

            scatter_1 = axes[2].scatter(
                completed_df.pulse_amp,
                completed_df.avg_Mz_final,
                c=completed_df.max_sigma_delta_abs,
                cmap="cividis",
                s=85,
                edgecolors="black",
                linewidths=0.4,
            )
            axes[2].set_title("Completed runs: final average \$M_z\$")
            axes[2].set_xlabel("pulse amplitude")
            axes[2].set_ylabel("final ⟨M_z⟩")
            axes[2].grid(alpha=0.2)
            colorbar(scatter_1, ax=axes[2], label="max |σΔ|")

            scatter_2 = axes[3].scatter(
                completed_df.pre_pulse_sigma_delta_absmax,
                completed_df.max_sigma_delta_abs,
                c=completed_df.avg_Mz_final,
                cmap="RdBu_r",
                s=85,
                edgecolors="black",
                linewidths=0.4,
            )
            xmin, xmax = finite_extrema(completed_df.pre_pulse_sigma_delta_absmax)
            axes[3].plot([xmin, xmax], [xmin, xmax], linestyle="--", color="black", linewidth=1.5)
            axes[3].set_title("Pre-pulse vs full-window σΔ")
            axes[3].set_xlabel("pre-pulse max |σΔ|")
            axes[3].set_ylabel("full-window max |σΔ|")
            axes[3].grid(alpha=0.2)
            colorbar(scatter_2, ax=axes[3], label="final ⟨M_z⟩")

            fig
            """
        ),
        code_cell(
            r"""
            isempty(trace_cases) && error("No completed runs with observables.tsv are available for time-trace plots.")

            ncols = 2
            nrows = ceil(Int, n_trace_cases / ncols)
            fig, axgrid = subplots(nrows, ncols, figsize=(20, 4.8 * nrows), constrained_layout=true, sharex=true)
            axes = flatten_axes(axgrid)

            for (ax, row) in zip(axes, eachrow(trace_cases))
                obs_df = load_observables(row)
                pulse_vals = pulse_trace(obs_df.t, row)

                ax.plot(obs_df.t, obs_df.avg_Mz, color="#0b3c5d", linewidth=2.4, label=L"\langle M_z \rangle")
                ax.axvspan(row.pulse_start, row.pulse_start + row.pulse_width, color="#f4d35e", alpha=0.18)
                ax.axvline(row.pulse_start, color="#c1121f", linestyle="--", linewidth=1.2)
                ax.axvline(row.pulse_start + row.pulse_width, color="#c1121f", linestyle="--", linewidth=1.2)
                ax.set_title(
                    @sprintf(
                        "%s | pa=%.2f, γso=%.2f, Jsd=%.2f, Jex=%.2f, Bx=%.4f",
                        String(row.run_id), row.pulse_amp, row.gamma_so, row.j_sd, row.j_ex, row.Bx
                    )
                )
                ax.set_xlabel("time")
                ax.set_ylabel(L"\langle M_z \rangle")
                ax.set_ylim(-1.05, 1.05)
                ax.grid(alpha=0.2)

                ax2 = ax.twinx()
                ax2.plot(obs_df.t, pulse_vals, color="#c1121f", linewidth=2.0, alpha=0.85, label="Δ_left(t)")
                ax2.set_ylabel("left-lead pulse Δ(t)")
                pulse_scale = max(abs(Float64(row.bias_amplitude)), 1e-8)
                ax2.set_ylim(-0.05 * pulse_scale, 1.15 * pulse_scale)
            end

            if n_trace_cases < length(axes)
                for ax in axes[(n_trace_cases + 1):length(axes)]
                    ax.remove()
                end
            end

            fig.suptitle("Time evolution of ⟨Mz⟩ with the applied left-lead pulse", y=1.02, fontsize=20)
            fig
            """
        ),
        code_cell(
            r"""
            top_case = trace_cases[1, :]
            top_obs = load_observables(top_case)
            top_pulse = pulse_trace(top_obs.t, top_case)

            sigma_z_cols = site_columns(top_obs, "sigma_z_site")
            sigma_eq_z_cols = site_columns(top_obs, "sigma_eq_z_site")
            sigma_delta_z_cols = site_columns(top_obs, "sigma_delta_z_site")
            local_sigma_z_matrix = columns_to_matrix(top_obs, sigma_z_cols)
            local_sigma_delta_z_matrix = columns_to_matrix(top_obs, sigma_delta_z_cols)
            total_sigma_z = sum_over_columns(top_obs, sigma_z_cols)
            total_sigma_eq_z = sum_over_columns(top_obs, sigma_eq_z_cols)
            total_sigma_delta_z = sum_over_columns(top_obs, sigma_delta_z_cols)
            mean_sigma_delta_z = mean_over_columns(top_obs, sigma_delta_z_cols)
            late_pre_mask = late_prepulse_mask(top_obs.t, top_case.pulse_start; fraction=0.2)
            late_pre_baseline_total_sigma_delta_z = masked_mean(total_sigma_delta_z, late_pre_mask)
            total_sigma_delta_z_relaxed = total_sigma_delta_z .- late_pre_baseline_total_sigma_delta_z
            local_sigma_delta_z_baseline = columnwise_masked_mean(local_sigma_delta_z_matrix, late_pre_mask)
            local_sigma_delta_z_relaxed = size(local_sigma_delta_z_matrix, 2) == 0 ? local_sigma_delta_z_matrix : (local_sigma_delta_z_matrix .- permutedims(local_sigma_delta_z_baseline))

            fig, axes = subplots(3, 1, figsize=(16, 12), constrained_layout=true, sharex=true)

            axes[1].plot(top_obs.t, top_obs.avg_Mz, color="#0b3c5d", linewidth=2.5)
            axes[1].plot(top_obs.t, top_obs.neel_z, color="#6c757d", linewidth=1.8, linestyle="--")
            axes[1].set_ylabel("spin order")
            axes[1].set_title("Best pulse-driven candidate: spin response over time")
            axes[1].legend([L"\langle M_z \rangle", L"N_z"], loc="best")
            axes[1].grid(alpha=0.2)

            axes[2].plot(top_obs.t, top_pulse, color="#c1121f", linewidth=2.2)
            axes[2].set_ylabel("Δ_left(t)")
            axes[2].set_title("Applied anti-symmetric pulse profile")
            axes[2].grid(alpha=0.2)

            axes[3].plot(top_obs.t, top_obs.I_left, color="#1d3557", linewidth=2.0)
            axes[3].plot(top_obs.t, top_obs.I_right, color="#e07a5f", linewidth=2.0)
            if any(isfinite, mean_sigma_delta_z)
                axes[3].plot(top_obs.t, mean_sigma_delta_z, color="#2a9d8f", linewidth=1.6, linestyle=":")
                axes[3].legend(["I_left", "I_right", "mean σΔ_z"], loc="best")
            else
                axes[3].legend(["I_left", "I_right"], loc="best")
            end
            axes[3].set_xlabel("time")
            axes[3].set_ylabel("transport / spin")
            axes[3].set_title("Currents and mean nonequilibrium spin density")
            axes[3].grid(alpha=0.2)

            for ax in axes
                ax.axvline(top_case.pulse_start, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
                ax.axvline(top_case.pulse_start + top_case.pulse_width, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
            end

            fig
            """
        ),
        markdown_cell(
            """
            ## Electronic Spin-Density Diagnostics

            Extra diagnostics for the strongest completed pulse-driven case:

            - total electronic spin density along `x`, `y`, and `z`,
            - raw nonequilibrium-minus-equilibrium traces and corresponding
              pulse-induced views after subtracting the late-pre-pulse baseline,
            - site-resolved local spin-density maps to reveal spatial structure
              in all three Cartesian components.

            The raw `σΔ = σ - σ_eq` trace can remain nonzero before the pulse here
            because the coupled dynamics already relaxes from the initial condition
            before `pulse_start`. The extra baseline-corrected views subtract the
            mean over the late pre-pulse window to isolate pulse-induced changes.
            """
        ),
        code_cell(
            r"""
            components = ["x", "y", "z"]
            fig, axes = subplots(length(components), 4, figsize=(26, 15), constrained_layout=true, sharex=false)

            for (row_idx, comp) in enumerate(components)
                sigma_cols = site_columns(top_obs, "sigma_$(comp)_site")
                sigma_eq_cols = site_columns(top_obs, "sigma_eq_$(comp)_site")
                sigma_delta_cols = site_columns(top_obs, "sigma_delta_$(comp)_site")

                local_sigma_matrix = columns_to_matrix(top_obs, sigma_cols)
                local_sigma_delta_matrix = columns_to_matrix(top_obs, sigma_delta_cols)
                total_sigma = sum_over_columns(top_obs, sigma_cols)
                total_sigma_eq = sum_over_columns(top_obs, sigma_eq_cols)
                total_sigma_delta = sum_over_columns(top_obs, sigma_delta_cols)

                late_pre_baseline_total_sigma_delta = masked_mean(total_sigma_delta, late_pre_mask)
                total_sigma_delta_relaxed = total_sigma_delta .- late_pre_baseline_total_sigma_delta
                local_sigma_delta_baseline = columnwise_masked_mean(local_sigma_delta_matrix, late_pre_mask)
                local_sigma_delta_relaxed = size(local_sigma_delta_matrix, 2) == 0 ? local_sigma_delta_matrix : (local_sigma_delta_matrix .- permutedims(local_sigma_delta_baseline))

                total_ax = axes[row_idx, 1]
                pulse_ax = axes[row_idx, 2]
                raw_map_ax = axes[row_idx, 3]
                corrected_map_ax = axes[row_idx, 4]

                sigma_handles = Any[]
                sigma_labels = String[]

                push!(sigma_handles, total_ax.plot(top_obs.t, total_sigma, color="#1d3557", linewidth=2.2)[1])
                push!(sigma_labels, "total σ_$(comp)")
                if any(isfinite, total_sigma_eq)
                    push!(sigma_handles, total_ax.plot(top_obs.t, total_sigma_eq, color="#6c757d", linewidth=1.6, linestyle="--")[1])
                    push!(sigma_labels, "total σ_eq,$(comp)")
                end
                if any(isfinite, total_sigma_delta)
                    push!(sigma_handles, total_ax.plot(top_obs.t, total_sigma_delta, color="#2a9d8f", linewidth=1.8, linestyle=":")[1])
                    push!(sigma_labels, "total σΔ_$(comp)")
                end
                total_ax.axvline(top_case.pulse_start, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
                total_ax.axvline(top_case.pulse_start + top_case.pulse_width, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
                total_ax.set_title("Raw total spin density along $(comp)")
                total_ax.set_xlabel("time")
                total_ax.set_ylabel("Σ σ_$(comp)")
                total_ax.grid(alpha=0.2)
                total_ax.legend(sigma_handles, sigma_labels, loc="best")

                pulse_ax.plot(top_obs.t, total_sigma_delta_relaxed, color="#2a9d8f", linewidth=2.0)
                pulse_ax.axhline(0.0, color="#6c757d", linestyle=":", linewidth=1.2)
                pulse_ax.axvline(top_case.pulse_start, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
                pulse_ax.axvline(top_case.pulse_start + top_case.pulse_width, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
                pulse_ax.set_title("Pulse-induced total σΔ_$(comp)")
                pulse_ax.set_xlabel("time")
                pulse_ax.set_ylabel("Σ σΔ_$(comp) - late-pre mean")
                pulse_ax.grid(alpha=0.2)

                heatmap_source = size(local_sigma_delta_matrix, 2) > 0 ? local_sigma_delta_matrix : local_sigma_matrix
                heatmap_title = size(local_sigma_delta_matrix, 2) > 0 ? "Raw local σΔ_$(comp) map" : "Raw local σ_$(comp) map"
                heatmap_label = size(local_sigma_delta_matrix, 2) > 0 ? "σΔ_$(comp)" : "σ_$(comp)"
                local_cols = size(local_sigma_delta_matrix, 2) > 0 ? sigma_delta_cols : sigma_cols
                site_ticks = collect(1:size(heatmap_source, 2))
                site_labels = [string(site_index(col)) for col in local_cols]

                if size(heatmap_source, 2) > 0
                    im = raw_map_ax.imshow(
                        permutedims(heatmap_source);
                        aspect="auto",
                        origin="lower",
                        extent=[minimum(top_obs.t), maximum(top_obs.t), 0.5, size(heatmap_source, 2) + 0.5],
                        cmap="coolwarm",
                    )
                    raw_map_ax.axvline(top_case.pulse_start, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
                    raw_map_ax.axvline(top_case.pulse_start + top_case.pulse_width, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
                    raw_map_ax.set_title(heatmap_title)
                    raw_map_ax.set_xlabel("time")
                    raw_map_ax.set_ylabel("site index")
                    raw_map_ax.set_yticks(site_ticks)
                    raw_map_ax.set_yticklabels(site_labels)
                    colorbar(im, ax=raw_map_ax, label=heatmap_label)
                else
                    raw_map_ax.text(0.5, 0.5, "No local columns found", ha="center", va="center", transform=raw_map_ax.transAxes)
                    raw_map_ax.set_axis_off()
                end

                relaxed_cols = sigma_delta_cols
                relaxed_ticks = collect(1:size(local_sigma_delta_relaxed, 2))
                relaxed_labels = [string(site_index(col)) for col in relaxed_cols]

                if size(local_sigma_delta_relaxed, 2) > 0
                    im_relaxed = corrected_map_ax.imshow(
                        permutedims(local_sigma_delta_relaxed);
                        aspect="auto",
                        origin="lower",
                        extent=[minimum(top_obs.t), maximum(top_obs.t), 0.5, size(local_sigma_delta_relaxed, 2) + 0.5],
                        cmap="coolwarm",
                    )
                    corrected_map_ax.axvline(top_case.pulse_start, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
                    corrected_map_ax.axvline(top_case.pulse_start + top_case.pulse_width, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
                    corrected_map_ax.set_title("Pulse-induced local σΔ_$(comp) map")
                    corrected_map_ax.set_xlabel("time")
                    corrected_map_ax.set_ylabel("site index")
                    corrected_map_ax.set_yticks(relaxed_ticks)
                    corrected_map_ax.set_yticklabels(relaxed_labels)
                    colorbar(im_relaxed, ax=corrected_map_ax, label="σΔ_$(comp) - late-pre mean")
                else
                    corrected_map_ax.text(0.5, 0.5, "No local σΔ columns found", ha="center", va="center", transform=corrected_map_ax.transAxes)
                    corrected_map_ax.set_axis_off()
                end
            end

            fig
            """
        ),
        code_cell(
            r"""
            zoom_mask = (top_obs.t .>= top_case.pulse_start - 40.0) .& (top_obs.t .<= top_case.pulse_start + 120.0)
            zoom_df = top_obs[zoom_mask, :]
            zoom_pulse = top_pulse[zoom_mask]

            fig, axes = subplots(1, 2, figsize=(18, 5.5), constrained_layout=true)

            axes[1].plot(zoom_df.t, zoom_df.avg_Mz, color="#0b3c5d", linewidth=2.4)
            axes[1].plot(zoom_df.t, zoom_df.neel_z, color="#6c757d", linewidth=1.8, linestyle="--")
            axes[1].set_title("Onset window: magnetic response")
            axes[1].set_xlabel("time")
            axes[1].set_ylabel("spin order")
            axes[1].grid(alpha=0.2)
            axes[1].legend([L"\langle M_z \rangle", L"N_z"], loc="best")

            axes[2].plot(zoom_df.t, zoom_pulse, color="#c1121f", linewidth=2.2)
            axes[2].plot(zoom_df.t, zoom_df.I_left, color="#1d3557", linewidth=1.8)
            axes[2].set_title("Onset window: pulse and left current")
            axes[2].set_xlabel("time")
            axes[2].set_ylabel("drive / current")
            axes[2].grid(alpha=0.2)
            axes[2].legend(["Δ_left(t)", "I_left"], loc="best")

            fig
            """
        ),
        code_cell(
            r"""
            fig, axes = subplots(1, 2, figsize=(18, 6), constrained_layout=true)

            gamma_vals_1, pulse_vals_1, final_mat = build_heatmap_matrix(completed_df, :gamma_so, :pulse_amp, :avg_Mz_final)
            im1 = axes[1].imshow(final_mat; aspect="auto", origin="lower", cmap="viridis")
            axes[1].set_title("avg_Mz_final heatmap")
            axes[1].set_xlabel("pulse_amp")
            axes[1].set_ylabel("gamma_so")
            axes[1].set_xticks(collect(0:length(pulse_vals_1)-1))
            axes[1].set_xticklabels(string.(pulse_vals_1))
            axes[1].set_yticks(collect(0:length(gamma_vals_1)-1))
            axes[1].set_yticklabels(string.(gamma_vals_1))
            colorbar(im1, ax=axes[1], label="avg_Mz_final")

            gamma_vals_2, jsd_vals_2, score_mat = build_heatmap_matrix(completed_df, :gamma_so, :j_sd, :pulse_driven_score)
            im2 = axes[2].imshow(score_mat; aspect="auto", origin="lower", cmap="magma")
            axes[2].set_title("pulse_driven_score heatmap")
            axes[2].set_xlabel("j_sd")
            axes[2].set_ylabel("gamma_so")
            axes[2].set_xticks(collect(0:length(jsd_vals_2)-1))
            axes[2].set_xticklabels(string.(jsd_vals_2))
            axes[2].set_yticks(collect(0:length(gamma_vals_2)-1))
            axes[2].set_yticklabels(string.(gamma_vals_2))
            colorbar(im2, ax=axes[2], label="pulse_driven_score")

            fig
            """
        ),
        markdown_cell(
            """
            ## Notes

            - This notebook assumes the Julia kernel has `CSV`, `DataFrames`,
              `JSON3`, and `PyPlot` available.
            - The figures are all generated through `PyPlot`, not through Python
              cells or matplotlib imports.
            - The pulse overlay uses the same anti-symmetric lead protocol as the
              driver: the plotted curve is `Δ_left(t)`, while the right lead gets
              `-Δ_left(t)`.
            - The ranking remains intentionally conservative because several of
              the strongest completed cases still have large pre-pulse
              `σΔ` contributions.
            """
        ),
    ]

    nb["cells"] = cells
    NOTEBOOK_PATH.write_text(nbf.writes(nb), encoding="utf-8")


if __name__ == "__main__":
    build_notebook()
