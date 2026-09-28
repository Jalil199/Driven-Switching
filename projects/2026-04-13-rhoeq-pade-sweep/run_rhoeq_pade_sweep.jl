include(joinpath(@__DIR__, "..", "2026-04-13-fixed-llg-no-pulse-relaxation", "run_fixed_llg_no_pulse_scan.jl"))

using Dates
using DelimitedFiles

function sigma_eq_from_rho(ρ::Matrix{ComplexF64}, p_model)
    obs = ObservablesTDNEGF(p_model; N_tmax=1, N_leads=2)
    obs.idx = 1
    obs_σ_i_eq!(ρ, p_model, obs)
    return copy(obs.σx_i_eq[:, :, 1])
end

function rhoeq_metrics(cfg::RunConfig)
    pstate = init_params_blocks(cfg)
    p_model = pstate.p_model
    site_ranges_couple = [get_sub(i, p_model.N_loc) for i in 1:p_model.N_sites]
    sys = init_params_heisenberg_llg(cfg)
    update_H_e!(p_model, site_ranges_couple, sys.dipoles[:, :, 1, 1], cfg.j_sd)
    ρ = ρ_eq(0.0, cfg.β, p_model.H_ab, cfg.N_λ2, cfg.Nx, cfg.Ny, cfg.Nσ, cfg.N_orb)
    σeq = sigma_eq_from_rho(ρ, p_model)
    return (; ρ, σeq, p_model)
end

function build_cases(results_root::AbstractString)
    shared = (
        run_id="pending",
        initial_configuration="FM",
        Hs_update_mode="nonequilibrium_only",
        results_root=String(results_root),
        Nx=2,
        Ny=2,
        dt=0.1,
        t_end=1000.0,
        bias_amplitude=0.0,
        left_bias_scale=0.0,
        right_bias_scale=0.0,
        pulse_start=1000.0,
        pulse_ramp=1.0,
        pulse_width=1.0,
        j_sd=0.3,
        j_ex=-0.02,
        j_ani=0.03,
        damping=0.3,
        γso=0.05,
        N_λ1=49,
    )
    return RunConfig[
        RunConfig(; shared..., run_id="pade20_gso0p05", N_λ2=20),
        RunConfig(; shared..., run_id="pade30_gso0p05", N_λ2=30),
        RunConfig(; shared..., run_id="pade40_gso0p05", N_λ2=40),
    ]
end

function enrich_summary!(summary::Dict{String,Any}, cfg::RunConfig, baseline_metrics, metrics)
    summary["probe_mode"] = "fixed_llg_no_pulse_rhoeq_pade_sweep"
    summary["rhoeq_pade_poles"] = cfg.N_λ2
    summary["rhoeq_sigma_eq_absmax_static"] = maximum(abs, metrics.σeq)
    summary["rhoeq_trace_real"] = real(tr(metrics.ρ))
    summary["rhoeq_trace_imag"] = imag(tr(metrics.ρ))
    if baseline_metrics !== nothing
        summary["rhoeq_maxabs_diff_vs_pade20"] = maximum(abs, metrics.ρ .- baseline_metrics.ρ)
        summary["sigma_eq_maxabs_diff_vs_pade20"] = maximum(abs, metrics.σeq .- baseline_metrics.σeq)
    else
        summary["rhoeq_maxabs_diff_vs_pade20"] = 0.0
        summary["sigma_eq_maxabs_diff_vs_pade20"] = 0.0
    end
    return summary
end

function write_scan_summary(path::AbstractString, summaries::Vector{Dict{String,Any}})
    headers = [
        "run_id",
        "status",
        "rhoeq_pade_poles",
        "gamma_so",
        "rhoeq_sigma_eq_absmax_static",
        "sigma_eq_maxabs_diff_vs_pade20",
        "rhoeq_maxabs_diff_vs_pade20",
        "max_sigma_delta_abs",
        "final_sigma_delta_absmax",
        "late_sigma_delta_absmax_mean",
        "late_sigma_delta_absmax_min",
        "late_sigma_delta_absmax_max",
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

function write_summary_md(path::AbstractString, prompt::AbstractString, summaries::Vector{Dict{String,Any}}, results_root::AbstractString)
    status = all(get(s, "status", "") == "done" for s in summaries) ? "done" : "blocked"
    open(path, "w") do io
        println(io, "- prompt: ", prompt)
        println(io, "- procedures: created a fixed-spin, zero-bias no-pulse sweep that keeps the same `gamma_so=0.05`, `j_sd=0.3`, `j_ex=-0.02`, `j_ani=0.03`, `damping=0.3`, `Nx=Ny=2`, and `t_end=1000`, while varying only the equilibrium Padé-pole count `N_λ2 ∈ {20,30,40}`. For each case, computed the static `rho_eq`, measured its change relative to the `N_λ2=20` baseline, then propagated the electrons and recorded the late-time `sigma_delta` residual.")
        println(io, "- generated data: ", joinpath(results_root, "scan_summary.tsv"))
        println(io, "- status: ", status)
    end
end

function main()
    run_stamp = Dates.format(Dates.now(Dates.UTC), "yyyymmddTHHMMSS")
    results_root = joinpath(@__DIR__, "results", "rhoeq_pade_sweep_" * run_stamp)
    mkpath(results_root)

    cases = build_cases(results_root)
    summaries = Dict{String,Any}[]

    baseline_metrics = nothing
    for cfg in cases
        println("Running ", cfg.run_id, " ...")
        metrics = rhoeq_metrics(cfg)
        baseline_metrics === nothing && (baseline_metrics = metrics)
        summary = run_fixed_llg_no_pulse_case_and_write(cfg)
        enrich_summary!(summary, cfg, baseline_metrics, metrics)
        write_json_object(joinpath(results_dir(cfg), "summary.json"), summary)
        push!(summaries, summary)
    end

    summary_tsv = joinpath(results_root, "scan_summary.tsv")
    write_scan_summary(summary_tsv, summaries)
    write_summary_md(
        joinpath(@__DIR__, "summary.md"),
        "hay alguna forma de mejorar el rho_eq? talvez incluyendo mas polos en su calculo, y asi volvemos a comparar?",
        summaries,
        results_root,
    )

    println("rho_eq Padé sweep complete. Results root: ", results_root)
    println("Summary table: ", summary_tsv)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
