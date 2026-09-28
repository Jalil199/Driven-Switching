using Distributed
using Dates
using DelimitedFiles

tdnegf_project_root() = abspath(joinpath(@__DIR__, "..", "..", "..", "TDNEGF"))

project_root = tdnegf_project_root()
case_driver_path = joinpath(@__DIR__, "run_switching_case.jl")

worker_count = min(3, max(Sys.CPU_THREADS - 1, 1))
nprocs() == 1 && addprocs(worker_count; exeflags="--project=$(project_root)")

@everywhere begin
    using Pkg
    Pkg.activate($project_root)
end
@everywhere include($case_driver_path)

function build_pilot_cases(results_root::AbstractString)
    shared = (
        initial_configuration="FM",
        results_root=String(results_root),
        dt=0.1,
        t_end=20.0,
        pulse_start=5.0,
        pulse_ramp=2.0,
        pulse_width=100.0,
    )
    return RunConfig[
        RunConfig(; shared..., run_id="pilot_fm_bias1p0_jsd0p10_jani0p10", bias_amplitude=1.0, j_sd=0.10, j_ani=0.10),
        RunConfig(; shared..., run_id="pilot_fm_bias1p2_jsd0p10_jani0p10", bias_amplitude=1.2, j_sd=0.10, j_ani=0.10),
        RunConfig(; shared..., run_id="pilot_fm_bias1p2_jsd0p12_jani0p08", bias_amplitude=1.2, j_sd=0.12, j_ani=0.08),
    ]
end

@everywhere function safe_run(cfg::RunConfig)
    try
        return run_switching_case_and_write(cfg)
    catch err
        return Dict(
            "run_id" => cfg.run_id,
            "initial_configuration" => cfg.initial_configuration,
            "results_dir" => results_dir(cfg),
            "status" => "failed",
            "j_sd" => cfg.j_sd,
            "j_ex" => cfg.j_ex,
            "j_ani" => cfg.j_ani,
            "bias_amplitude" => cfg.bias_amplitude,
            "pulse_start" => cfg.pulse_start,
            "pulse_ramp" => cfg.pulse_ramp,
            "pulse_width" => cfg.pulse_width,
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

function write_global_summary(path::AbstractString, summaries::Vector{Dict{String,Any}})
    headers = [
        "run_id",
        "status",
        "initial_configuration",
        "j_sd",
        "j_ex",
        "j_ani",
        "bias_amplitude",
        "pulse_start",
        "pulse_ramp",
        "pulse_width",
        "t_final",
        "current_left_final",
        "current_right_final",
        "avg_Mz_final",
        "avg_Mz_min",
        "neel_z_peak_abs",
        "full_reversal_detected",
        "partial_reversal_flag",
        "results_dir",
        "error",
    ]
    open(path, "w") do io
        println(io, join(headers, '\t'))
        for summary in summaries
            row = [summary_value(summary, key) for key in headers]
            writedlm(io, permutedims(row), '\t')
        end
    end
end

function main()
    sweep_id = "pilot_sweep_" * Dates.format(Dates.now(Dates.UTC), "yyyymmddTHHMMSS")
    results_root = joinpath(@__DIR__, "results", sweep_id)
    mkpath(results_root)

    cases = build_pilot_cases(results_root)
    summaries = pmap(safe_run, cases)
    summary_path = joinpath(results_root, "pilot_summary.tsv")
    write_global_summary(summary_path, summaries)

    println("Pilot sweep complete. Results root: ", results_root)
    println("Global summary table: ", summary_path)
    for summary in summaries
        println(summary["run_id"], " => ", summary["status"])
    end
end

main()
