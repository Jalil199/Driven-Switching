include(joinpath(@__DIR__, "..", "2026-04-08-tdnegf-spin-baseline", "run_switching_case.jl"))

using Dates
using LinearAlgebra

@inline function tail_start_index(n::Int, fraction::Float64=0.1)
    return max(1, floor(Int, (1.0 - fraction) * n) + 1)
end

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
        "gamma_so" => cfg.γso,
        "rhoeq_pade_poles" => cfg.N_λ2,
        "wall_elapsed_s" => time() - wall_start,
        "timestamp_utc" => string(Dates.now(Dates.UTC)),
        "results_dir" => results_dir(cfg),
        "pid" => getpid(),
    ]
    write_json_object(path, progress)
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

function sigma_eq_from_rho(ρ::Matrix{ComplexF64}, p_model)
    obs = ObservablesTDNEGF(p_model; N_tmax=1, N_leads=2)
    obs.idx = 1
    obs_σ_i_eq!(ρ, p_model, obs)
    return copy(obs.σx_i_eq[:, :, 1])
end

function sigma_delta_metrics(obs)
    sigma_delta = obs.σx_i .- obs.σx_i_eq
    sigma_delta_absmax_t = [maximum(abs, @view sigma_delta[:, :, it]) for it in axes(sigma_delta, 3)]
    tail10_idx = tail_start_index(length(sigma_delta_absmax_t), 0.1)
    tail10 = @view sigma_delta_absmax_t[tail10_idx:end]
    tail_ts = @view obs.t[tail10_idx:end]
    return Dict{String,Any}(
        "t_final" => obs.t[end],
        "max_sigma_delta_abs" => maximum(sigma_delta_absmax_t),
        "final_sigma_delta_absmax" => sigma_delta_absmax_t[end],
        "late_sigma_delta_absmax_mean" => sum(tail10) / length(tail10),
        "late_sigma_delta_absmax_min" => minimum(tail10),
        "late_sigma_delta_absmax_max" => maximum(tail10),
        "late10_sigma_delta_absmax_slope" => linear_tail_slope(tail_ts, tail10),
        "tail10_first_sigma_delta_absmax" => tail10[1],
        "tail10_last_sigma_delta_absmax" => tail10[end],
        "tail10_drop_sigma_delta_absmax" => tail10[1] - tail10[end],
        "pre_pulse_sigma_delta_absmax" => sigma_delta_absmax_t[end],
        "final_to_peak_sigma_delta_ratio" => sigma_delta_absmax_t[end] / maximum(sigma_delta_absmax_t),
    )
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
    ρ = ρ_eq(0.0, cfg.β, p_model.H_ab, p_model.N_λ2, p_model.Nx, p_model.Ny, p_model.Nσ, p_model.N_orb)

    prob = ODEProblem(eom_tdnegf_blocks!, u, (cfg.t0, cfg.t_end), p_blocks)
    intg = init(
        prob, Vern7();
        dt=cfg.dt, save_everystep=false, adaptive=true, dense=false,
        reltol=cfg.reltol, abstol=cfg.abstol,
    )

    n_steps = Int(round((cfg.t_end - cfg.t0) / cfg.dt))
    obs = ObservablesTDNEGF(p_model; N_tmax=n_steps, N_leads=length(p_blocks.blocks))
    progress_path = get(ENV, "FERMILINK_PROGRESS_PATH", "")
    progress_every = max(1, env_int("FERMILINK_PROGRESS_EVERY_STEPS", 500))
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

    return (; obs, ρ, p_model)
end

function build_summary(cfg::RunConfig, result)
    metrics = sigma_delta_metrics(result.obs)
    σeq = sigma_eq_from_rho(result.ρ, result.p_model)
    summary = Pair{String,Any}[
        "run_id" => cfg.run_id,
        "status" => "done",
        "results_dir" => results_dir(cfg),
        "branch_reference_run_id" => "goal6_A_main_0003",
        "probe_mode" => "brancha_fixed_llg_no_pulse_rhoeq_pade_soc_check",
        "static_spin_texture" => true,
        "lead_bias_forced_zero" => true,
        "gamma_so" => cfg.γso,
        "rhoeq_pade_poles" => cfg.N_λ2,
        "j_sd" => cfg.j_sd,
        "j_ex" => cfg.j_ex,
        "j_ani" => cfg.j_ani,
        "damping" => cfg.damping,
        "pulse_start" => cfg.pulse_start,
        "pulse_width" => cfg.pulse_width,
        "rhoeq_trace_real" => real(tr(result.ρ)),
        "rhoeq_trace_imag" => imag(tr(result.ρ)),
        "rhoeq_sigma_eq_absmax_static" => maximum(abs, σeq),
    ]
    append!(summary, collect(metrics))
    return summary
end

function main(args)
    cfg = parse_cli_args(RunConfig(), args)
    case_dir = results_dir(cfg)
    mkpath(case_dir)
    progress_path = joinpath(case_dir, "progress.json")
    ENV["FERMILINK_PROGRESS_PATH"] = progress_path
    write_params_json(joinpath(case_dir, "params.json"), cfg)
    try
        result = run_fixed_llg_static_rhoeq_case(cfg)
        write_observables(joinpath(case_dir, "observables.tsv"), result.obs, cfg)
        summary = build_summary(cfg, result)
        write_json_object(joinpath(case_dir, "summary.json"), summary)
        println("Run complete. Results written to: ", case_dir)
    catch err
        write_json_object(progress_path, Pair{String,Any}[
            "run_id" => cfg.run_id,
            "status" => "failed",
            "timestamp_utc" => string(Dates.now(Dates.UTC)),
            "error" => sprint(showerror, err),
            "results_dir" => case_dir,
        ])
        summary = failure_summary(cfg, err)
        write_json_object(joinpath(case_dir, "summary.json"), summary)
        rethrow()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
