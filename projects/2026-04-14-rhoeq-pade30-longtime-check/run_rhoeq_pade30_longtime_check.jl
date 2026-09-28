include(joinpath(@__DIR__, "..", "2026-04-13-rhoeq-longtime-relaxation", "run_rhoeq_longtime_relaxation.jl"))

using Dates

const BASELINE_OBS = joinpath(
    @__DIR__,
    "..",
    "2026-04-13-rhoeq-longtime-relaxation",
    "results",
    "rhoeq_longtime_20260413T205653",
    "pade20_gso0p05_t2000",
    "observables.tsv",
)

const BASELINE_SUMMARY = joinpath(
    @__DIR__,
    "..",
    "2026-04-13-rhoeq-longtime-relaxation",
    "results",
    "rhoeq_longtime_20260413T205653",
    "pade20_gso0p05_t2000",
    "summary.json",
)

function baseline_case(results_root::AbstractString)
    cfg = RunConfig(
        run_id="pade20_gso0p05_t2000_reference",
        initial_configuration="FM",
        Hs_update_mode="nonequilibrium_only",
        results_root=String(results_root),
        Nx=2,
        Ny=2,
        dt=0.1,
        t_end=2000.0,
        bias_amplitude=0.0,
        left_bias_scale=0.0,
        right_bias_scale=0.0,
        pulse_start=2000.0,
        pulse_ramp=1.0,
        pulse_width=1.0,
        j_sd=0.3,
        j_ex=-0.02,
        j_ani=0.03,
        damping=0.3,
        γso=0.05,
        N_λ1=49,
        N_λ2=20,
    )
    metrics = sigma_delta_metrics_from_file(BASELINE_OBS)
    static = static_rhoeq(cfg)
    summary = init_summary(cfg, metrics)
    summary["status"] = "done"
    summary["results_dir"] = dirname(BASELINE_OBS)
    summary["baseline_source_summary"] = BASELINE_SUMMARY
    enrich_summary!(summary, static.ρ, static.σeq, nothing)
    return (; cfg, summary, ρ=static.ρ, σeq=static.σeq)
end

function pade30_case(results_root::AbstractString)
    return RunConfig(
        run_id="pade30_gso0p05_t2000",
        initial_configuration="FM",
        Hs_update_mode="nonequilibrium_only",
        results_root=String(results_root),
        Nx=2,
        Ny=2,
        dt=0.1,
        t_end=2000.0,
        bias_amplitude=0.0,
        left_bias_scale=0.0,
        right_bias_scale=0.0,
        pulse_start=2000.0,
        pulse_ramp=1.0,
        pulse_width=1.0,
        j_sd=0.3,
        j_ex=-0.02,
        j_ani=0.03,
        damping=0.3,
        γso=0.05,
        N_λ1=49,
        N_λ2=30,
    )
end

function pade30_failure_summary(cfg::RunConfig, err)
    return Dict{String,Any}(
        "run_id" => cfg.run_id,
        "status" => "failed",
        "results_dir" => results_dir(cfg),
        "gamma_so" => cfg.γso,
        "j_sd" => cfg.j_sd,
        "j_ex" => cfg.j_ex,
        "j_ani" => cfg.j_ani,
        "rhoeq_pade_poles" => cfg.N_λ2,
        "probe_mode" => "fixed_llg_no_pulse_rhoeq_longtime",
        "static_spin_texture" => true,
        "lead_bias_forced_zero" => true,
        "error" => sprint(showerror, err),
        "timestamp_utc" => string(Dates.now(Dates.UTC)),
    )
end

function write_summary_md(path::AbstractString, prompt::AbstractString, summaries::Vector{Dict{String,Any}}, results_root::AbstractString)
    status = all(get(s, "status", "") == "done" for s in summaries) ? "done" : "blocked"
    open(path, "w") do io
        println(io, "- prompt: ", prompt)
        println(io, "- procedures: reused the completed `t_end=2000` fixed-spin no-pulse baseline for `N_λ2=20`, then ran the same disconnected-start TDNEGF setup with only the equilibrium Padé pole count changed to `N_λ2=30`. Wrote a two-row comparison table and a curve overlay plot for `max_i |sigma_delta,i(t)|` so the late-time agreement can be checked directly.")
        println(io, "- generated data: ", joinpath(results_root, "scan_summary.tsv"), " ; ", joinpath(results_root, "pade30_gso0p05_t2000", "observables.tsv"), " ; ", joinpath(results_root, "rhoeq_pade20_vs_pade30_t2000.png"))
        println(io, "- status: ", status)
    end
end

function main()
    run_stamp = Dates.format(Dates.now(Dates.UTC), "yyyymmddTHHMMSS")
    results_root = get(ENV, "FERMILINK_RESULTS_ROOT", joinpath(@__DIR__, "results", "rhoeq_pade30_longtime_" * run_stamp))
    mkpath(results_root)

    base = baseline_case(results_root)
    summaries = Dict{String,Any}[base.summary]

    cfg30 = pade30_case(results_root)
    mkpath(results_dir(cfg30))
    write_params_json(joinpath(results_dir(cfg30), "params.json"), cfg30)
    try
        result30 = run_fixed_llg_static_rhoeq_case(cfg30)
        write_observables(joinpath(results_dir(cfg30), "observables.tsv"), result30.obs, cfg30)
        metrics30 = sigma_delta_metrics_from_result(result30)
        σeq30 = sigma_eq_from_rho(result30.ρ, result30.p_model)
        summary30 = init_summary(cfg30, metrics30)
        enrich_summary!(summary30, result30.ρ, σeq30, (; ρ=base.ρ, σeq=base.σeq, summary=base.summary))
        write_json_dict(joinpath(results_dir(cfg30), "summary.json"), summary30)
        push!(summaries, summary30)
    catch err
        summary30 = pade30_failure_summary(cfg30, err)
        write_json_dict(joinpath(results_dir(cfg30), "summary.json"), summary30)
        push!(summaries, summary30)
        summary_tsv = joinpath(results_root, "scan_summary.tsv")
        write_scan_summary(summary_tsv, summaries)
        write_summary_md(
            joinpath(@__DIR__, "summary.md"),
            "pues propon una forma de que no muera y muestrame los resultados cuando esten",
            summaries,
            results_root,
        )
        rethrow()
    end

    summary_tsv = joinpath(results_root, "scan_summary.tsv")
    write_scan_summary(summary_tsv, summaries)
    write_summary_md(
        joinpath(@__DIR__, "summary.md"),
        "pues propon una forma de que no muera y muestrame los resultados cuando esten",
        summaries,
        results_root,
    )

    println("rho_eq pade30 long-time check complete. Results root: ", results_root)
    println("Summary table: ", summary_tsv)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
