using Distributed
using Dates
using DelimitedFiles

tdnegf_project_root() = abspath(joinpath(@__DIR__, "..", "..", "..", "TDNEGF"))

project_root = tdnegf_project_root()
case_driver_path = joinpath(@__DIR__, "..", "2026-04-08-tdnegf-spin-baseline", "run_switching_case.jl")

const MAX_PARALLEL_LOAD = something(tryparse(Int, get(ENV, "TDNEGF_MAX_WORKERS", "34")), 34)
const DEFAULT_WORKER_COUNT = min(MAX_PARALLEL_LOAD, max(Sys.CPU_THREADS - 1, 1))
const INCLUDE_HIGH_DAMPING_STRIP = lowercase(get(ENV, "TDNEGF_GOAL6_INCLUDE_HIGH_DAMPING_STRIP", "true")) in ("1", "true", "yes")
const MANIFEST_ONLY = lowercase(get(ENV, "TDNEGF_MANIFEST_ONLY", "false")) in ("1", "true", "yes")

include(case_driver_path)

function sanitized_token(x::Real)
    s = replace(string(x), "-" => "m", "." => "p")
    return replace(s, "+" => "")
end

function mode_tag(mode::AbstractString)
    return mode == "full_spin_density" ? "hsfull" : "hsneq"
end

function t_end_for_width(pulse_width::Real)
    return pulse_width == 3000 ? 6000.0 : 5000.0
end

function make_case(;
    results_root::AbstractString,
    branch::AbstractString,
    branch_index::Integer,
    Hs_update_mode::AbstractString,
    pulse_amp::Real,
    gamma_so::Real,
    j_sd::Real,
    j_ex::Real,
    j_ani::Real,
    damping::Real,
    pulse_width::Real,
)
    run_id = join([
        "goal6",
        String(branch),
        lpad(string(branch_index), 4, '0'),
        "pa" * sanitized_token(pulse_amp),
        "gso" * sanitized_token(gamma_so),
        "jsd" * sanitized_token(j_sd),
        "jex" * sanitized_token(j_ex),
        "jani" * sanitized_token(j_ani),
        "damp" * sanitized_token(damping),
        "twidth" * sanitized_token(pulse_width),
        mode_tag(Hs_update_mode),
    ], "_")
    return RunConfig(;
        run_id=run_id,
        initial_configuration="FM",
        results_root=String(results_root),
        Nx=2,
        Ny=2,
        dt=0.1,
        t_end=t_end_for_width(pulse_width),
        pulse_start=1000.0,
        pulse_ramp=5.0,
        pulse_width=Float64(pulse_width),
        Bx=0.0,
        left_bias_scale=0.0,
        right_bias_scale=1.0,
        Hs_update_mode=String(Hs_update_mode),
        bias_amplitude=Float64(pulse_amp),
        γso=Float64(gamma_so),
        j_sd=Float64(j_sd),
        j_ex=Float64(j_ex),
        j_ani=Float64(j_ani),
        damping=Float64(damping),
    )
end

function branch_label(cfg::RunConfig)
    startswith(cfg.run_id, "goal6_A_high_damping") && return "A_high_damping_strip"
    startswith(cfg.run_id, "goal6_A_main") && return "A_main_hsneq_stabilization"
    startswith(cfg.run_id, "goal6_B_main") && return "B_main_hsfull_retention"
    startswith(cfg.run_id, "goal6_C_center") && return "C_center_hs_mode_crosscheck"
    return "unknown"
end

function build_scan_cases(results_root::AbstractString)
    cases = RunConfig[]

    idx = 0
    for pulse_width in [1000.0, 1500.0, 2000.0, 2500.0, 3000.0]
        for damping in [0.3, 0.5, 0.7]
            for j_ani in [0.0, 0.01, 0.03, 0.05]
                idx += 1
                push!(cases, make_case(;
                    results_root,
                    branch="A_main",
                    branch_index=idx,
                    Hs_update_mode="nonequilibrium_only",
                    pulse_amp=1.5,
                    gamma_so=0.9,
                    j_sd=0.3,
                    j_ex=-0.02,
                    j_ani,
                    damping,
                    pulse_width,
                ))
            end
        end
    end

    if INCLUDE_HIGH_DAMPING_STRIP
        for (idx_high, pulse_width) in enumerate([1000.0, 1500.0, 2000.0, 2500.0, 3000.0])
            push!(cases, make_case(;
                results_root,
                branch="A_high_damping",
                branch_index=idx_high,
                Hs_update_mode="nonequilibrium_only",
                pulse_amp=1.5,
                gamma_so=0.9,
                j_sd=0.3,
                j_ex=-0.02,
                j_ani=0.03,
                damping=1.0,
                pulse_width,
            ))
        end
    end

    idx = 0
    for pulse_width in [1500.0, 2000.0, 2500.0, 3000.0]
        for pulse_amp in [1.75, 2.0, 2.25]
            for damping in [0.4, 0.5, 0.7, 1.0]
                idx += 1
                push!(cases, make_case(;
                    results_root,
                    branch="B_main",
                    branch_index=idx,
                    Hs_update_mode="full_spin_density",
                    pulse_amp,
                    gamma_so=0.9,
                    j_sd=0.3,
                    j_ex=-0.01,
                    j_ani=0.0,
                    damping,
                    pulse_width,
                ))
            end
        end
    end

    center_specs = [
        ("C_center_A", 1, 1.5, 0.9, 0.3, -0.02, 0.03, 0.3, 1500.0, "nonequilibrium_only"),
        ("C_center_A", 2, 1.5, 0.9, 0.3, -0.02, 0.03, 0.3, 1500.0, "full_spin_density"),
        ("C_center_B", 3, 2.0, 0.9, 0.3, -0.01, 0.0, 0.5, 1500.0, "nonequilibrium_only"),
        ("C_center_B", 4, 2.0, 0.9, 0.3, -0.01, 0.0, 0.5, 1500.0, "full_spin_density"),
    ]
    for (branch, idx_center, pulse_amp, gamma_so, j_sd, j_ex, j_ani, damping, pulse_width, hs_mode) in center_specs
        push!(cases, make_case(;
            results_root,
            branch,
            branch_index=idx_center,
            Hs_update_mode=hs_mode,
            pulse_amp,
            gamma_so,
            j_sd,
            j_ex,
            j_ani,
            damping,
            pulse_width,
        ))
    end

    return cases
end

function ensure_workers(ncases::Int)
    return min(nworkers(), MAX_PARALLEL_LOAD, ncases)
end

function summary_value(summary::Dict{String,Any}, key::String)
    value = get(summary, key, missing)
    value === missing && return "missing"
    value isa AbstractVector && return join(string.(value), ",")
    return string(value)
end

function write_table(path::AbstractString, headers::Vector{String}, rows::Vector{Dict{String,Any}})
    open(path, "w") do io
        println(io, join(headers, '\t'))
        for row in rows
            values = [summary_value(row, key) for key in headers]
            writedlm(io, permutedims(values), '\t')
        end
    end
end

function manifest_row(cfg::RunConfig)
    return Dict(
        "run_id" => cfg.run_id,
        "branch" => branch_label(cfg),
        "initial_configuration" => cfg.initial_configuration,
        "Nx" => cfg.Nx,
        "Ny" => cfg.Ny,
        "Hs_update_mode" => cfg.Hs_update_mode,
        "j_ani" => cfg.j_ani,
        "damping" => cfg.damping,
        "Bx" => cfg.Bx,
        "pulse_amp" => cfg.bias_amplitude,
        "gamma_so" => cfg.γso,
        "j_sd" => cfg.j_sd,
        "j_ex" => cfg.j_ex,
        "left_bias_scale" => cfg.left_bias_scale,
        "right_bias_scale" => cfg.right_bias_scale,
        "pulse_start" => cfg.pulse_start,
        "pulse_ramp" => cfg.pulse_ramp,
        "pulse_width" => cfg.pulse_width,
        "t_end" => cfg.t_end,
        "results_dir" => results_dir(cfg),
    )
end

function write_scan_config(path::AbstractString, results_root::AbstractString, worker_cap::Int, ncases::Int)
    pairs = Pair{String,Any}[
        "scan_id" => basename(String(results_root)),
        "results_root" => String(results_root),
        "requested_case_count" => ncases,
        "branch_A_main_count" => 60,
        "branch_A_high_damping_strip_count" => INCLUDE_HIGH_DAMPING_STRIP ? 5 : 0,
        "branch_B_main_count" => 48,
        "branch_C_center_crosscheck_count" => 4,
        "max_parallel_load" => MAX_PARALLEL_LOAD,
        "actual_worker_count" => worker_cap,
        "fixed_t_end_default" => 5000.0,
        "fixed_t_end_for_pulse_width_3000" => 6000.0,
        "fixed_pulse_start" => 1000.0,
        "fixed_pulse_ramp" => 5.0,
        "fixed_Bx" => 0.0,
        "fixed_left_bias_scale" => 0.0,
        "fixed_right_bias_scale" => 1.0,
        "queue_policy" => "goal6_focused_local_followup_with_labeled_branches",
        "optional_second_strips_launched" => false,
        "timestamp_utc" => string(Dates.now(Dates.UTC)),
    ]
    write_json_object(path, pairs)
end

@everywhere begin
    using Pkg
    Pkg.activate($project_root)
end

if !MANIFEST_ONLY
    nprocs() == 1 && addprocs(DEFAULT_WORKER_COUNT; exeflags="--project=$(project_root)")
    @everywhere include($case_driver_path)
end

@everywhere function goal6_branch_label_from_run_id(run_id::AbstractString)
    startswith(run_id, "goal6_A_high_damping") && return "A_high_damping_strip"
    startswith(run_id, "goal6_A_main") && return "A_main_hsneq_stabilization"
    startswith(run_id, "goal6_B_main") && return "B_main_hsfull_retention"
    startswith(run_id, "goal6_C_center") && return "C_center_hs_mode_crosscheck"
    return "unknown"
end

@everywhere function safe_run(cfg::RunConfig)
    try
        summary = run_switching_case_and_write(cfg)
        summary["branch"] = goal6_branch_label_from_run_id(cfg.run_id)
        return summary
    catch err
        return Dict(
            "run_id" => cfg.run_id,
            "branch" => goal6_branch_label_from_run_id(cfg.run_id),
            "initial_configuration" => cfg.initial_configuration,
            "results_dir" => results_dir(cfg),
            "status" => "failed",
            "Hs_update_mode" => cfg.Hs_update_mode,
            "j_sd" => cfg.j_sd,
            "j_ex" => cfg.j_ex,
            "j_ani" => cfg.j_ani,
            "damping" => cfg.damping,
            "Bx" => cfg.Bx,
            "pulse_amp" => cfg.bias_amplitude,
            "bias_amplitude" => cfg.bias_amplitude,
            "gamma_so" => cfg.γso,
            "γso" => cfg.γso,
            "left_bias_scale" => cfg.left_bias_scale,
            "right_bias_scale" => cfg.right_bias_scale,
            "pulse_start" => cfg.pulse_start,
            "pulse_ramp" => cfg.pulse_ramp,
            "pulse_width" => cfg.pulse_width,
            "t_end" => cfg.t_end,
            "error" => sprint(showerror, err),
            "timestamp_utc" => string(Dates.now(Dates.UTC)),
        )
    end
end

function main()
    scan_id = get(ENV, "TDNEGF_SCAN_ID", "queued_goal6_" * Dates.format(Dates.now(Dates.UTC), "yyyymmddTHHMMSSZ"))
    results_root = joinpath(@__DIR__, "results", scan_id)
    mkpath(results_root)

    cases = build_scan_cases(results_root)
    worker_cap = MANIFEST_ONLY ? 0 : ensure_workers(length(cases))

    manifest_headers = [
        "run_id", "branch", "initial_configuration", "Nx", "Ny", "j_ani", "damping", "Bx",
        "Hs_update_mode",
        "pulse_amp", "gamma_so", "j_sd", "j_ex",
        "left_bias_scale", "right_bias_scale",
        "pulse_start", "pulse_ramp", "pulse_width", "t_end", "results_dir",
    ]
    write_table(joinpath(results_root, "case_manifest.tsv"), manifest_headers, manifest_row.(cases))
    write_scan_config(joinpath(results_root, "scan_config.json"), results_root, worker_cap, length(cases))

    if MANIFEST_ONLY
        println("Goal 6 manifest written. Results root: ", results_root)
        println("Case count: ", length(cases))
        return nothing
    end

    summaries = pmap(safe_run, cases)
    summary_headers = [
        "run_id", "branch", "status", "initial_configuration",
        "Hs_update_mode",
        "pulse_amp", "bias_amplitude", "gamma_so", "γso",
        "j_sd", "j_ex", "j_ani", "damping", "Bx",
        "left_bias_scale", "right_bias_scale",
        "pulse_start", "pulse_ramp", "pulse_width", "t_end", "t_final",
        "current_left_final", "current_right_final", "current_peak_abs", "current_balance_peak_abs",
        "avg_Mz_initial", "avg_Mz_final", "avg_Mz_min", "avg_Mz_max", "neel_z_peak_abs",
        "max_site_charge_deviation_from_1", "max_classical_spin_norm_deviation",
        "max_sigma_noneq_abs", "max_sigma_eq_abs", "max_sigma_delta_abs",
        "pre_pulse_sample_count", "pre_pulse_duration_realized",
        "pre_pulse_sigma_noneq_absmax", "pre_pulse_sigma_eq_absmax", "pre_pulse_sigma_delta_absmax",
        "full_reversal_detected", "partial_reversal_flag",
        "llg_drive_convention", "equilibrium_reference_saved", "nonequilibrium_minus_equilibrium_saved",
        "results_dir", "error",
    ]
    summary_path = joinpath(results_root, "scan_summary.tsv")
    write_table(summary_path, summary_headers, summaries)

    println("Goal 6 queued scan complete. Results root: ", results_root)
    println("Case manifest: ", joinpath(results_root, "case_manifest.tsv"))
    println("Global summary table: ", summary_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
