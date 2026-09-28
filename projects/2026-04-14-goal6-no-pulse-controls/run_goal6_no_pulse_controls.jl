using Dates
using DelimitedFiles

case_driver_path = joinpath(@__DIR__, "..", "2026-04-08-tdnegf-spin-baseline", "run_switching_case.jl")
include(case_driver_path)

const PROJECT_ROOT = @__DIR__
const RESULTS_ROOT = joinpath(PROJECT_ROOT, "results", "matched_no_pulse_controls")
const STATUS_PATH = joinpath(PROJECT_ROOT, "run_status.json")

function sanitized_token(x::Real)
    s = replace(string(x), "-" => "m", "." => "p")
    return replace(s, "+" => "")
end

function mode_tag(mode::AbstractString)
    return mode == "full_spin_density" ? "hsfull" : "hsneq"
end

function run_id_for_case(label::AbstractString, cfg::RunConfig)
    return join([
        label,
        "nopulse",
        "gso" * sanitized_token(cfg.γso),
        "jsd" * sanitized_token(cfg.j_sd),
        "jex" * sanitized_token(cfg.j_ex),
        "jani" * sanitized_token(cfg.j_ani),
        "damp" * sanitized_token(cfg.damping),
        "twidth" * sanitized_token(cfg.pulse_width),
        mode_tag(cfg.Hs_update_mode),
    ], "_")
end

function matched_no_pulse_cases()
    base = [
        (
            label="goal6_A_main_0003",
            cfg=RunConfig(
                initial_configuration="FM",
                Hs_update_mode="nonequilibrium_only",
                Nx=2,
                Ny=2,
                dt=0.1,
                t_end=5000.0,
                pulse_start=1000.0,
                pulse_ramp=5.0,
                pulse_width=1000.0,
                Bx=0.0,
                left_bias_scale=0.0,
                right_bias_scale=1.0,
                bias_amplitude=0.0,
                γso=0.9,
                j_sd=0.3,
                j_ex=-0.02,
                j_ani=0.03,
                damping=0.3,
                results_root=RESULTS_ROOT,
            ),
            source_run_id="goal6_A_main_0003_pa1p5_gso0p9_jsd0p3_jexm0p02_jani0p03_damp0p3_twidth1000p0_hsneq",
        ),
        (
            label="goal6_A_main_0034",
            cfg=RunConfig(
                initial_configuration="FM",
                Hs_update_mode="nonequilibrium_only",
                Nx=2,
                Ny=2,
                dt=0.1,
                t_end=5000.0,
                pulse_start=1000.0,
                pulse_ramp=5.0,
                pulse_width=2000.0,
                Bx=0.0,
                left_bias_scale=0.0,
                right_bias_scale=1.0,
                bias_amplitude=0.0,
                γso=0.9,
                j_sd=0.3,
                j_ex=-0.02,
                j_ani=0.01,
                damping=0.7,
                results_root=RESULTS_ROOT,
            ),
            source_run_id="goal6_A_main_0034_pa1p5_gso0p9_jsd0p3_jexm0p02_jani0p01_damp0p7_twidth2000p0_hsneq",
        ),
    ]

    return [
        (
            label=spec.label,
            source_run_id=spec.source_run_id,
            cfg=(spec.cfg.run_id = run_id_for_case(spec.label, spec.cfg); spec.cfg),
        ) for spec in base
    ]
end

function write_manifest(path::AbstractString, cases)
    headers = [
        "label",
        "source_run_id",
        "run_id",
        "Hs_update_mode",
        "bias_amplitude",
        "gamma_so",
        "j_sd",
        "j_ex",
        "j_ani",
        "damping",
        "pulse_start",
        "pulse_ramp",
        "pulse_width",
        "t_end",
        "results_dir",
    ]
    open(path, "w") do io
        println(io, join(headers, '\t'))
        for case in cases
            cfg = case.cfg
            row = [
                case.label,
                case.source_run_id,
                cfg.run_id,
                cfg.Hs_update_mode,
                string(cfg.bias_amplitude),
                string(cfg.γso),
                string(cfg.j_sd),
                string(cfg.j_ex),
                string(cfg.j_ani),
                string(cfg.damping),
                string(cfg.pulse_start),
                string(cfg.pulse_ramp),
                string(cfg.pulse_width),
                string(cfg.t_end),
                results_dir(cfg),
            ]
            writedlm(io, permutedims(row), '\t')
        end
    end
end

function write_scan_config(path::AbstractString, cases)
    pairs = Pair{String,Any}[
        "scan_id" => basename(RESULTS_ROOT),
        "results_root" => RESULTS_ROOT,
        "requested_case_count" => length(cases),
        "control_type" => "matched_goal6_no_pulse_reference",
        "bias_amplitude_fixed" => 0.0,
        "matched_pulse_metadata_retained" => true,
        "source_project" => joinpath(PROJECT_ROOT, "..", "2026-04-10-goal6-retained-switch"),
        "timestamp_utc" => string(Dates.now(Dates.UTC)),
    ]
    write_json_object(path, pairs)
end

function scan_summary_row(case, summary)
    return Dict(
        "label" => case.label,
        "source_run_id" => case.source_run_id,
        "run_id" => case.cfg.run_id,
        "status" => get(summary, "status", "unknown"),
        "t_final" => get(summary, "t_final", missing),
        "avg_Mz_final" => get(summary, "avg_Mz_final", missing),
        "avg_Mz_min" => get(summary, "avg_Mz_min", missing),
        "max_sigma_delta_abs" => get(summary, "max_sigma_delta_abs", missing),
        "pre_pulse_sigma_delta_absmax" => get(summary, "pre_pulse_sigma_delta_absmax", missing),
        "results_dir" => results_dir(case.cfg),
    )
end

function write_scan_summary(path::AbstractString, rows)
    headers = [
        "label",
        "source_run_id",
        "run_id",
        "status",
        "t_final",
        "avg_Mz_final",
        "avg_Mz_min",
        "max_sigma_delta_abs",
        "pre_pulse_sigma_delta_absmax",
        "results_dir",
    ]
    open(path, "w") do io
        println(io, join(headers, '\t'))
        for row in rows
            values = [get(row, key, missing) for key in headers]
            writedlm(io, permutedims(string.(values)), '\t')
        end
    end
end

function write_status(;
    status::AbstractString,
    current_case::Union{Nothing,AbstractString}=nothing,
    completed_cases::Int=0,
    total_cases::Int=0,
    note::Union{Nothing,AbstractString}=nothing,
)
    pairs = Pair{String,Any}[
        "status" => String(status),
        "current_case" => something(current_case, "none"),
        "completed_cases" => completed_cases,
        "total_cases" => total_cases,
        "results_root" => RESULTS_ROOT,
        "timestamp_utc" => string(Dates.now(Dates.UTC)),
    ]
    note === nothing || push!(pairs, "note" => String(note))
    write_json_object(STATUS_PATH, pairs)
end

function main()
    mkpath(RESULTS_ROOT)
    cases = matched_no_pulse_cases()
    write_manifest(joinpath(RESULTS_ROOT, "case_manifest.tsv"), cases)
    write_scan_config(joinpath(RESULTS_ROOT, "scan_config.json"), cases)
    write_status(status="running", completed_cases=0, total_cases=length(cases), note="launcher_initialized")

    rows = Dict{String,Any}[]
    for (idx, case) in enumerate(cases)
        write_status(status="running", current_case=case.cfg.run_id, completed_cases=idx - 1, total_cases=length(cases))
        println("Running matched no-pulse control: ", case.cfg.run_id)
        summary = run_switching_case_and_write(case.cfg)
        push!(rows, scan_summary_row(case, summary))
        write_status(status="running", current_case=case.cfg.run_id, completed_cases=idx, total_cases=length(cases), note="case_finished")
    end

    write_scan_summary(joinpath(RESULTS_ROOT, "scan_summary.tsv"), rows)
    write_status(status="finished", completed_cases=length(cases), total_cases=length(cases), note="all_cases_finished")
end

main()
