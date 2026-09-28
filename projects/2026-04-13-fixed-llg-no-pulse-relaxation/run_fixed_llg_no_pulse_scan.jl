include(joinpath(@__DIR__, "..", "2026-04-08-tdnegf-spin-baseline", "run_switching_case.jl"))

using Dates
using DelimitedFiles

@inline function tail_start_index(n::Int, fraction::Float64=0.1)
    return max(1, floor(Int, (1.0 - fraction) * n) + 1)
end

function summary_value(summary::Dict{String,Any}, key::String)
    value = get(summary, key, missing)
    value === missing && return "missing"
    value isa AbstractVector && return join(string.(value), ",")
    return string(value)
end

function run_fixed_llg_no_pulse_case(cfg::RunConfig)
    cfg.Hs_update_mode = validated_hs_update_mode(cfg.Hs_update_mode)
    pstate = init_params_blocks(cfg)
    p_model = pstate.p_model
    p_blocks = pstate.p_blocks
    u = copy(pstate.u0)
    sys = init_params_heisenberg_llg(cfg)

    site_ranges_couple = [get_sub(i, p_model.N_loc) for i in 1:p_model.N_sites]

    # Freeze the spin texture but include the static s-d background from the initial state.
    update_H_e!(p_model, site_ranges_couple, sys.dipoles[:, :, 1, 1], cfg.j_sd)
    p_blocks.Δ_blocks[1] = 0.0 + 0.0im
    p_blocks.Δ_blocks[2] = 0.0 + 0.0im

    prob = ODEProblem(eom_tdnegf_blocks!, u, (cfg.t0, cfg.t_end), p_blocks)
    intg = init(
        prob, Vern7();
        dt=cfg.dt, save_everystep=false, adaptive=true, dense=false,
        reltol=cfg.reltol, abstol=cfg.abstol,
    )

    n_steps = Int(round((cfg.t_end - cfg.t0) / cfg.dt))
    obs = ObservablesTDNEGF(p_model; N_tmax=n_steps, N_leads=length(p_blocks.blocks))

    for step_idx in 1:n_steps
        obs.idx = step_idx
        DifferentialEquations.step!(intg, cfg.dt, true)

        dv = pointer_blocks(intg.u, p_blocks.dims_ρ_ab, p_blocks.aux_layout)
        ρ = ρ_eq(0.0, cfg.β, p_model.H_ab, p_model.N_λ2, p_model.Nx, p_model.Ny, p_model.Nσ, p_model.N_orb)

        obs.t[step_idx] = intg.t
        obs_n_i!(dv, p_model, obs)
        obs_σ_i!(dv, p_model, obs)
        obs_Ixα!(dv, p_blocks, obs)
        obs_s_i!(sys.dipoles[:, :, 1, 1], p_model, obs)
        obs_σ_i_eq!(ρ, p_model, obs)
    end

    return (; obs, final_state=intg.u, sys, p_model, p_blocks)
end

function summarize_fixed_llg_case(result, cfg::RunConfig)
    base = summarize_run(result, cfg)
    obs = result.obs
    sigma_delta = obs.σx_i .- obs.σx_i_eq
    sigma_delta_absmax_t = [maximum(abs, @view sigma_delta[:, :, it]) for it in axes(sigma_delta, 3)]
    tail_idx = tail_start_index(length(sigma_delta_absmax_t))
    tail = @view sigma_delta_absmax_t[tail_idx:end]
    push!(base, "probe_mode" => "no_pulse_fixed_llg_static_sd")
    push!(base, "static_spin_texture" => true)
    push!(base, "lead_bias_forced_zero" => true)
    push!(base, "final_sigma_delta_absmax" => sigma_delta_absmax_t[end])
    push!(base, "late_sigma_delta_absmax_mean" => sum(tail) / length(tail))
    push!(base, "late_sigma_delta_absmax_min" => minimum(tail))
    push!(base, "late_sigma_delta_absmax_max" => maximum(tail))
    return base
end

function run_fixed_llg_no_pulse_case_and_write(cfg::RunConfig)
    case_dir = results_dir(cfg)
    mkpath(case_dir)
    write_params_json(joinpath(case_dir, "params.json"), cfg)
    try
        result = run_fixed_llg_no_pulse_case(cfg)
        write_observables(joinpath(case_dir, "observables.tsv"), result.obs, cfg)
        summary = summarize_fixed_llg_case(result, cfg)
        write_json_object(joinpath(case_dir, "summary.json"), summary)
        return summary_dict(summary)
    catch err
        summary = failure_summary(cfg, err)
        write_json_object(joinpath(case_dir, "summary.json"), summary)
        rethrow()
    end
end

function build_cases(results_root::AbstractString)
    shared = (
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
    )
    return RunConfig[
        RunConfig(; shared..., run_id="fixedllg_gso0p05", γso=0.05),
        RunConfig(; shared..., run_id="fixedllg_gso0p10", γso=0.10),
        RunConfig(; shared..., run_id="fixedllg_gso0p30", γso=0.30),
        RunConfig(; shared..., run_id="fixedllg_gso0p70", γso=0.70),
    ]
end

function write_scan_summary(path::AbstractString, summaries::Vector{Dict{String,Any}})
    headers = [
        "run_id",
        "status",
        "gamma_so",
        "j_sd",
        "j_ex",
        "j_ani",
        "t_final",
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
        println(io, "- procedures: created a fixed-spin, zero-bias relaxation probe based on `run_switching_case.jl`, applied the static s-d Hamiltonian from the initial FM texture once, propagated electrons only up to `t=1000`, and wrote per-case `observables.tsv` plus `summary.json` for four Rashba values.")
        println(io, "- generated data: ", joinpath(results_root, "scan_summary.tsv"))
        println(io, "- status: ", status)
    end
end

function main()
    run_stamp = Dates.format(Dates.now(Dates.UTC), "yyyymmddTHHMMSS")
    results_root = joinpath(@__DIR__, "results", "fixed_llg_no_pulse_" * run_stamp)
    mkpath(results_root)

    cases = build_cases(results_root)
    summaries = Dict{String,Any}[]
    for cfg in cases
        println("Running ", cfg.run_id, " ...")
        push!(summaries, run_fixed_llg_no_pulse_case_and_write(cfg))
    end

    summary_tsv = joinpath(results_root, "scan_summary.tsv")
    write_scan_summary(summary_tsv, summaries)
    write_summary_md(
        joinpath(@__DIR__, "summary.md"),
        "podriamos correr unos de prueba sin el pulso y con el llg fijo hasta 1000 para ver si la densidad se relaja a cero? si no estoy mal esto pasa cuando el parametro de rashba es pequeno",
        summaries,
        results_root,
    )

    println("Fixed-LLG no-pulse probe complete. Results root: ", results_root)
    println("Summary table: ", summary_tsv)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
