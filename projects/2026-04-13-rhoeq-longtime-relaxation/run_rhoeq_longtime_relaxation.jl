include(joinpath(@__DIR__, "..", "2026-04-13-fixed-llg-no-pulse-relaxation", "run_fixed_llg_no_pulse_scan.jl"))

using Dates
using DelimitedFiles

function env_int(name::AbstractString, default::Int)
    raw = get(ENV, name, "")
    isempty(raw) && return default
    value = tryparse(Int, raw)
    value === nothing && return default
    return value
end

function write_progress_json(path::AbstractString, cfg::RunConfig, step_idx::Int, n_steps::Int, t_now::Float64, wall_start::Float64, status::AbstractString)
    progress = Pair{String,Any}[
        "run_id" => cfg.run_id,
        "status" => String(status),
        "step_idx" => step_idx,
        "n_steps" => n_steps,
        "fraction_complete" => n_steps == 0 ? 0.0 : step_idx / n_steps,
        "t_now" => t_now,
        "t_end" => cfg.t_end,
        "wall_elapsed_s" => time() - wall_start,
        "timestamp_utc" => string(Dates.now(Dates.UTC)),
        "results_dir" => results_dir(cfg),
        "pid" => getpid(),
    ]
    write_json_object(path, progress)
end

function sigma_eq_from_rho(ρ::Matrix{ComplexF64}, p_model)
    obs = ObservablesTDNEGF(p_model; N_tmax=1, N_leads=2)
    obs.idx = 1
    obs_σ_i_eq!(ρ, p_model, obs)
    return copy(obs.σx_i_eq[:, :, 1])
end

function static_rhoeq(cfg::RunConfig)
    pstate = init_params_blocks(cfg)
    p_model = pstate.p_model
    sys = init_params_heisenberg_llg(cfg)
    site_ranges_couple = [get_sub(i, p_model.N_loc) for i in 1:p_model.N_sites]
    update_H_e!(p_model, site_ranges_couple, sys.dipoles[:, :, 1, 1], cfg.j_sd)
    ρ = ρ_eq(0.0, cfg.β, p_model.H_ab, cfg.N_λ2, cfg.Nx, cfg.Ny, cfg.Nσ, cfg.N_orb)
    σeq = sigma_eq_from_rho(ρ, p_model)
    return (; ρ, σeq)
end

function linear_tail_slope(ts::AbstractVector{<:Real}, ys::AbstractVector{<:Real})
    n = length(ts)
    n <= 1 && return 0.0
    tmean = sum(ts) / n
    ymean = sum(ys) / n
    denom = sum((t - tmean)^2 for t in ts)
    denom == 0 && return 0.0
    numer = sum((ts[i] - tmean) * (ys[i] - ymean) for i in eachindex(ts))
    return numer / denom
end

function sigma_delta_metrics(ts::AbstractVector{<:Real}, sigma_delta_absmax_t::AbstractVector{<:Real})
    tail10_idx = tail_start_index(length(sigma_delta_absmax_t), 0.1)
    tail5_idx = tail_start_index(length(sigma_delta_absmax_t), 0.05)
    tail10 = @view sigma_delta_absmax_t[tail10_idx:end]
    tail5 = @view sigma_delta_absmax_t[tail5_idx:end]
    prev5_start = max(1, tail5_idx - length(tail5))
    prev5_stop = max(prev5_start, tail5_idx - 1)
    prev5 = @view sigma_delta_absmax_t[prev5_start:prev5_stop]
    tail_ts = @view ts[tail10_idx:end]
    return Dict{String,Any}(
        "t_final" => ts[end],
        "max_sigma_delta_abs" => maximum(sigma_delta_absmax_t),
        "final_sigma_delta_absmax" => sigma_delta_absmax_t[end],
        "late_sigma_delta_absmax_mean" => sum(tail10) / length(tail10),
        "late_sigma_delta_absmax_min" => minimum(tail10),
        "late_sigma_delta_absmax_max" => maximum(tail10),
        "late5_sigma_delta_absmax_mean" => sum(tail5) / length(tail5),
        "prev5_sigma_delta_absmax_mean" => sum(prev5) / length(prev5),
        "late10_sigma_delta_absmax_slope" => linear_tail_slope(tail_ts, tail10),
        "tail10_first_sigma_delta_absmax" => tail10[1],
        "tail10_last_sigma_delta_absmax" => tail10[end],
        "tail10_drop_sigma_delta_absmax" => tail10[1] - tail10[end],
    )
end

function sigma_delta_metrics_from_result(result)
    obs = result.obs
    sigma_delta = obs.σx_i .- obs.σx_i_eq
    sigma_delta_absmax_t = [maximum(abs, @view sigma_delta[:, :, it]) for it in axes(sigma_delta, 3)]
    return sigma_delta_metrics(obs.t, sigma_delta_absmax_t)
end

function sigma_delta_metrics_from_file(path::AbstractString)
    open(path, "r") do io
        header = split(chomp(readline(io)), '\t')
        sigma_cols = findall(name -> startswith(name, "sigma_delta_"), header)
        t_col = findfirst(==("t"), header)
        t_col === nothing && error("Missing time column in observables file: $path")
        isempty(sigma_cols) && error("Missing sigma_delta columns in observables file: $path")

        ts = Float64[]
        sigma_delta_absmax_t = Float64[]
        for line in eachline(io)
            isempty(line) && continue
            fields = split(line, '\t')
            push!(ts, parse(Float64, fields[t_col]))
            row_absmax = 0.0
            for idx in sigma_cols
                row_absmax = max(row_absmax, abs(parse(Float64, fields[idx])))
            end
            push!(sigma_delta_absmax_t, row_absmax)
        end
        isempty(ts) && error("No data rows in observables file: $path")
        return sigma_delta_metrics(ts, sigma_delta_absmax_t)
    end
end

function init_summary(cfg::RunConfig, metrics::Dict{String,Any})
    summary = Dict{String,Any}(
        "run_id" => cfg.run_id,
        "status" => "done",
        "results_dir" => results_dir(cfg),
        "gamma_so" => cfg.γso,
        "j_sd" => cfg.j_sd,
        "j_ex" => cfg.j_ex,
        "j_ani" => cfg.j_ani,
        "rhoeq_pade_poles" => cfg.N_λ2,
        "probe_mode" => "fixed_llg_no_pulse_rhoeq_longtime",
        "static_spin_texture" => true,
        "lead_bias_forced_zero" => true,
    )
    merge!(summary, metrics)
    summary["final_to_peak_sigma_delta_ratio"] = summary["final_sigma_delta_absmax"] / summary["max_sigma_delta_abs"]
    return summary
end

function run_fixed_llg_static_rhoeq_case(cfg::RunConfig)
    cfg.Hs_update_mode = validated_hs_update_mode(cfg.Hs_update_mode)
    pstate = init_params_blocks(cfg)
    p_model = pstate.p_model
    p_blocks = pstate.p_blocks
    u = copy(pstate.u0)
    sys = init_params_heisenberg_llg(cfg)

    site_ranges_couple = [get_sub(i, p_model.N_loc) for i in 1:p_model.N_sites]

    update_H_e!(p_model, site_ranges_couple, sys.dipoles[:, :, 1, 1], cfg.j_sd)
    p_blocks.Δ_blocks[1] = 0.0 + 0.0im
    p_blocks.Δ_blocks[2] = 0.0 + 0.0im
    ρ = ρ_eq(0.0, cfg.β, p_model.H_ab, cfg.N_λ2, cfg.Nx, cfg.Ny, cfg.Nσ, cfg.N_orb)

    prob = ODEProblem(eom_tdnegf_blocks!, u, (cfg.t0, cfg.t_end), p_blocks)
    intg = init(
        prob, Vern7();
        dt=cfg.dt, save_everystep=false, adaptive=true, dense=false,
        reltol=cfg.reltol, abstol=cfg.abstol,
    )

    n_steps = Int(round((cfg.t_end - cfg.t0) / cfg.dt))
    obs = ObservablesTDNEGF(p_model; N_tmax=n_steps, N_leads=length(p_blocks.blocks))
    progress_path = get(ENV, "FERMILINK_PROGRESS_PATH", "")
    progress_every = max(1, env_int("FERMILINK_PROGRESS_EVERY_STEPS", 250))
    wall_start = time()
    if !isempty(progress_path)
        write_progress_json(progress_path, cfg, 0, n_steps, cfg.t0, wall_start, "running")
    end

    for step_idx in 1:n_steps
        obs.idx = step_idx
        DifferentialEquations.step!(intg, cfg.dt, true)

        dv = pointer_blocks(intg.u, p_blocks.dims_ρ_ab, p_blocks.aux_layout)
        obs.t[step_idx] = intg.t
        obs_n_i!(dv, p_model, obs)
        obs_σ_i!(dv, p_model, obs)
        obs_Ixα!(dv, p_blocks, obs)
        obs_s_i!(sys.dipoles[:, :, 1, 1], p_model, obs)
        obs_σ_i_eq!(ρ, p_model, obs)

        if !isempty(progress_path) && (step_idx == n_steps || step_idx % progress_every == 0)
            write_progress_json(progress_path, cfg, step_idx, n_steps, intg.t, wall_start, step_idx == n_steps ? "postprocessing" : "running")
        end
    end

    if !isempty(progress_path)
        write_progress_json(progress_path, cfg, n_steps, n_steps, intg.t, wall_start, "completed")
    end

    return (; obs, final_state=intg.u, sys, p_model, p_blocks, ρ)
end

function enrich_summary!(summary::Dict{String,Any}, ρ::Matrix{ComplexF64}, σeq::Matrix{Float64}, baseline=nothing)
    summary["rhoeq_sigma_eq_absmax_static"] = maximum(abs, σeq)
    if baseline !== nothing
        summary["rhoeq_maxabs_diff_vs_baseline"] = maximum(abs, ρ .- baseline.ρ)
        summary["sigma_eq_maxabs_diff_vs_baseline"] = maximum(abs, σeq .- baseline.σeq)
        summary["final_sigma_delta_diff_vs_baseline"] = summary["final_sigma_delta_absmax"] - baseline.summary["final_sigma_delta_absmax"]
        summary["late10_sigma_delta_mean_diff_vs_baseline"] = summary["late_sigma_delta_absmax_mean"] - baseline.summary["late_sigma_delta_absmax_mean"]
    else
        summary["rhoeq_maxabs_diff_vs_baseline"] = 0.0
        summary["sigma_eq_maxabs_diff_vs_baseline"] = 0.0
        summary["final_sigma_delta_diff_vs_baseline"] = 0.0
        summary["late10_sigma_delta_mean_diff_vs_baseline"] = 0.0
    end
    return summary
end

function write_scan_summary(path::AbstractString, summaries::Vector{Dict{String,Any}})
    headers = [
        "run_id",
        "status",
        "rhoeq_pade_poles",
        "gamma_so",
        "t_final",
        "rhoeq_sigma_eq_absmax_static",
        "rhoeq_maxabs_diff_vs_baseline",
        "sigma_eq_maxabs_diff_vs_baseline",
        "max_sigma_delta_abs",
        "final_sigma_delta_absmax",
        "late_sigma_delta_absmax_mean",
        "late5_sigma_delta_absmax_mean",
        "prev5_sigma_delta_absmax_mean",
        "late10_sigma_delta_absmax_slope",
        "tail10_first_sigma_delta_absmax",
        "tail10_last_sigma_delta_absmax",
        "tail10_drop_sigma_delta_absmax",
        "final_to_peak_sigma_delta_ratio",
        "final_sigma_delta_diff_vs_baseline",
        "late10_sigma_delta_mean_diff_vs_baseline",
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

write_json_dict(path::AbstractString, summary::Dict{String,Any}) = write_json_object(path, collect(summary))

function write_summary_md(path::AbstractString, prompt::AbstractString, summaries::Vector{Dict{String,Any}}, results_root::AbstractString)
    status = all(get(s, "status", "") == "done" for s in summaries) ? "done" : "blocked"
    open(path, "w") do io
        println(io, "- prompt: ", prompt)
        println(io, "- procedures: created a new fixed-spin, zero-bias no-pulse relaxation comparison for the same small-Rashba setup (`gamma_so=0.05`, `Nx=Ny=2`, `j_sd=0.3`, `j_ex=-0.02`, `j_ani=0.03`, `damping=0.3`) but extended the propagation to `t_end=2000` and compared the baseline `N_λ2=20` against an increased equilibrium pole count `N_λ2=60`. The driver computes the static `rho_eq` once per run, propagates the partitioned TDNEGF dynamics from `u0=0`, and records late-time residual metrics for `sigma_delta`.")
        println(io, "- generated data: ", joinpath(results_root, "scan_summary.tsv"))
        println(io, "- status: ", status)
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
    )
    return RunConfig[
        RunConfig(; shared..., run_id="pade20_gso0p05_t2000", N_λ2=20),
        RunConfig(; shared..., run_id="pade60_gso0p05_t2000", N_λ2=60),
    ]
end

function main()
    run_stamp = Dates.format(Dates.now(Dates.UTC), "yyyymmddTHHMMSS")
    results_root = get(ENV, "FERMILINK_RESULTS_ROOT", joinpath(@__DIR__, "results", "rhoeq_longtime_" * run_stamp))
    mkpath(results_root)

    cases = build_cases(results_root)
    summaries = Dict{String,Any}[]
    baseline = nothing

    for cfg in cases
        println("Running ", cfg.run_id, " ...")
        mkpath(results_dir(cfg))
        obs_path = joinpath(results_dir(cfg), "observables.tsv")
        if isfile(obs_path)
            metrics = sigma_delta_metrics_from_file(obs_path)
            static = static_rhoeq(cfg)
            ρ = static.ρ
            σeq = static.σeq
        else
            result = run_fixed_llg_static_rhoeq_case(cfg)
            write_observables(obs_path, result.obs, cfg)
            metrics = sigma_delta_metrics_from_result(result)
            ρ = result.ρ
            σeq = sigma_eq_from_rho(result.ρ, result.p_model)
        end
        summary = init_summary(cfg, metrics)
        if baseline === nothing
            baseline = (; ρ=ρ, σeq=σeq, summary=summary)
        end
        enrich_summary!(summary, ρ, σeq, baseline)
        write_params_json(joinpath(results_dir(cfg), "params.json"), cfg)
        write_json_dict(joinpath(results_dir(cfg), "summary.json"), summary)
        push!(summaries, summary)
    end

    summary_tsv = joinpath(results_root, "scan_summary.tsv")
    write_scan_summary(summary_tsv, summaries)
    write_summary_md(
        joinpath(@__DIR__, "summary.md"),
        "hagamos otra comparacion hasta 2000, con una mejor cantidad de polos para ver si la diferencia tiende a cero",
        summaries,
        results_root,
    )

    println("rho_eq long-time relaxation comparison complete. Results root: ", results_root)
    println("Summary table: ", summary_tsv)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
