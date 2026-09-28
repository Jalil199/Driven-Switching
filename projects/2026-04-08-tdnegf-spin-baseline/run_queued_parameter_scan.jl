using Distributed
using Dates
using DelimitedFiles

tdnegf_project_root() = abspath(joinpath(@__DIR__, "..", "..", "..", "TDNEGF"))

project_root = tdnegf_project_root()
case_driver_path = joinpath(@__DIR__, "run_switching_case.jl")

const MAX_PARALLEL_LOAD = something(tryparse(Int, get(ENV, "TDNEGF_MAX_WORKERS", "34")), 34)
const DEFAULT_RESULTS_PREFIX = "queued_parameter_scan_"
const DEFAULT_WORKER_COUNT = min(MAX_PARALLEL_LOAD, max(Sys.CPU_THREADS - 1, 1))

function sanitized_token(x::Real)
    s = replace(string(x), "-" => "m", "." => "p")
    return replace(s, "+" => "")
end

function requested_scan_values()
    return (
        pulse_amp=[1.5, 2.0],
        gamma_so=[0.7, 0.9],
        j_sd=[0.2, 0.3],
        j_ex=[-0.02, -0.01],
        Bx=[0.0],
        j_ani=[0.0, 0.03, 0.05, 0.08, 0.10],
        damping=[0.3, 0.5],
        pulse_start=[1000.0],
        pulse_ramp=[5.0],
        pulse_width=[1500.0],
        left_bias_scale=[0.0],
        right_bias_scale=[1.0],
        Hs_update_mode=["nonequilibrium_only", "full_spin_density"],
    )
end

function build_scan_cases(results_root::AbstractString)
    values = requested_scan_values()
    shared = (
        initial_configuration="FM",
        results_root=String(results_root),
        Nx=2,
        Ny=2,
        dt=0.1,
        t_end=4000.0,
        pulse_start=1000.0,
        pulse_ramp=5.0,
        pulse_width=1500.0,
        j_ani=0.05,
        damping=0.5,
        Bx=0.0,
        left_bias_scale=0.0,
        right_bias_scale=1.0,
    )

    cases = RunConfig[]
    point_index = 0
    for pulse_amp in values.pulse_amp
        for gamma_so in values.gamma_so
            for j_sd in values.j_sd
                for j_ex in values.j_ex
                    for Bx in values.Bx
                        for j_ani in values.j_ani
                            for damping in values.damping
                                for pulse_start in values.pulse_start
                                    for pulse_ramp in values.pulse_ramp
                                        for pulse_width in values.pulse_width
                                            for left_bias_scale in values.left_bias_scale
                                                for right_bias_scale in values.right_bias_scale
                                                    point_index += 1
                                                    point_tag = join([
                                                        "scan",
                                                        lpad(string(point_index), 4, '0'),
                                                        "pa" * sanitized_token(pulse_amp),
                                                        "gso" * sanitized_token(gamma_so),
                                                        "jsd" * sanitized_token(j_sd),
                                                        "jex" * sanitized_token(j_ex),
                                                        "bx" * sanitized_token(Bx),
                                                        "jani" * sanitized_token(j_ani),
                                                        "damp" * sanitized_token(damping),
                                                        "tstart" * sanitized_token(pulse_start),
                                                        "tramp" * sanitized_token(pulse_ramp),
                                                        "twidth" * sanitized_token(pulse_width),
                                                        "lbs" * sanitized_token(left_bias_scale),
                                                        "rbs" * sanitized_token(right_bias_scale),
                                                    ], "_")
                                                    for hs_update_mode in values.Hs_update_mode
                                                        mode_tag = hs_update_mode == "full_spin_density" ? "hsfull" : "hsneq"
                                                        run_id = point_tag * "_" * mode_tag
                                                        push!(cases, RunConfig(;
                                                            shared...,
                                                            run_id=run_id,
                                                            Hs_update_mode=hs_update_mode,
                                                            bias_amplitude=pulse_amp,
                                                            γso=gamma_so,
                                                            j_sd=j_sd,
                                                            j_ex=j_ex,
                                                            Bx=Bx,
                                                            j_ani=j_ani,
                                                            damping=damping,
                                                            pulse_start=pulse_start,
                                                            pulse_ramp=pulse_ramp,
                                                            pulse_width=pulse_width,
                                                            left_bias_scale=left_bias_scale,
                                                            right_bias_scale=right_bias_scale,
                                                        ))
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return cases
end

function ensure_workers(ncases::Int)
    return min(nworkers(), MAX_PARALLEL_LOAD, ncases)
end

nprocs() == 1 && addprocs(DEFAULT_WORKER_COUNT; exeflags="--project=$(project_root)")

@everywhere begin
    using Pkg
    Pkg.activate($project_root)
end
@everywhere include($case_driver_path)

@everywhere function safe_run(cfg::RunConfig)
    try
        return run_switching_case_and_write(cfg)
    catch err
        return Dict(
            "run_id" => cfg.run_id,
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
    parameter_point_count = ncases ÷ length(requested_scan_values().Hs_update_mode)
    pairs = Pair{String,Any}[
        "scan_id" => basename(String(results_root)),
        "results_root" => String(results_root),
        "requested_parameter_point_count" => parameter_point_count,
        "requested_case_count" => ncases,
        "max_parallel_load" => MAX_PARALLEL_LOAD,
        "actual_worker_count" => worker_cap,
        "pulse_amp_values" => requested_scan_values().pulse_amp,
        "gamma_so_values" => requested_scan_values().gamma_so,
        "j_sd_values" => requested_scan_values().j_sd,
        "j_ex_values" => requested_scan_values().j_ex,
        "Bx_values" => requested_scan_values().Bx,
        "j_ani_values" => requested_scan_values().j_ani,
        "damping_values" => requested_scan_values().damping,
        "pulse_start_values" => requested_scan_values().pulse_start,
        "pulse_ramp_values" => requested_scan_values().pulse_ramp,
        "pulse_width_values" => requested_scan_values().pulse_width,
        "left_bias_scale_values" => requested_scan_values().left_bias_scale,
        "right_bias_scale_values" => requested_scan_values().right_bias_scale,
        "Hs_update_mode_values" => requested_scan_values().Hs_update_mode,
        "fixed_Nx" => 2,
        "fixed_Ny" => 2,
        "fixed_t_end" => 4000.0,
        "llg_drive_convention" => "paired_ab_sigma_delta_vs_sigma_full",
        "equilibrium_reference_saved" => true,
        "nonequilibrium_minus_equilibrium_saved" => true,
        "queue_policy" => "all_parameter_points_parallel_with_paired_Hs_update_mode_backfill",
        "timestamp_utc" => string(Dates.now(Dates.UTC)),
    ]
    write_json_object(path, pairs)
end

function main()
    scan_id = get(ENV, "TDNEGF_SCAN_ID", DEFAULT_RESULTS_PREFIX * Dates.format(Dates.now(Dates.UTC), "yyyymmddTHHMMSS"))
    results_root = joinpath(@__DIR__, "results", scan_id)
    mkpath(results_root)

    cases = build_scan_cases(results_root)
    worker_cap = ensure_workers(length(cases))

    manifest_headers = [
        "run_id", "initial_configuration", "Nx", "Ny", "j_ani", "damping", "Bx",
        "Hs_update_mode",
        "pulse_amp", "gamma_so", "j_sd", "j_ex",
        "left_bias_scale", "right_bias_scale",
        "pulse_start", "pulse_ramp", "pulse_width", "t_end", "results_dir",
    ]
    write_table(joinpath(results_root, "case_manifest.tsv"), manifest_headers, manifest_row.(cases))
    write_scan_config(joinpath(results_root, "scan_config.json"), results_root, worker_cap, length(cases))

    summaries = pmap(safe_run, cases)
    summary_headers = [
        "run_id", "status", "initial_configuration",
        "Hs_update_mode",
        "pulse_amp", "bias_amplitude", "gamma_so", "γso",
        "j_sd", "j_ex", "j_ani", "damping", "Bx",
        "left_bias_scale", "right_bias_scale",
        "pulse_start", "pulse_ramp", "pulse_width", "t_final",
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

    println("Queued parameter scan complete. Results root: ", results_root)
    println("Case manifest: ", joinpath(results_root, "case_manifest.tsv"))
    println("Global summary table: ", summary_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
