include(joinpath(@__DIR__, "..", "2026-04-08-tdnegf-spin-baseline", "run_switching_case.jl"))

using Dates
using DifferentialEquations
using LinearAlgebra

function env_int(name::AbstractString, default::Int)
    raw = get(ENV, name, "")
    isempty(raw) && return default
    value = tryparse(Int, raw)
    value === nothing && return default
    return value
end

function write_progress_json(path::AbstractString, cfg::RunConfig, phase::AbstractString, step_idx::Int, n_steps::Int, t_now::Float64, wall_start::Float64, status::AbstractString)
    progress = Pair{String,Any}[
        "run_id" => cfg.run_id,
        "status" => String(status),
        "phase" => String(phase),
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

function sigma_eq_from_rho(ρ::AbstractMatrix{ComplexF64}, p_model)
    obs = ObservablesTDNEGF(p_model; N_tmax=1, N_leads=2)
    obs.idx = 1
    obs_σ_i_eq!(ρ, p_model, obs)
    return copy(obs.σx_i_eq[:, :, 1])
end

function sigma_from_rho(ρ::AbstractMatrix{ComplexF64}, p_model)
    ptr = (; ρ_ab = ρ)
    obs = ObservablesTDNEGF(p_model; N_tmax=1, N_leads=2)
    obs.idx = 1
    obs_σ_i!(ptr, p_model, obs)
    return copy(obs.σx_i[:, :, 1])
end

function sigma_delta_absmax_from_rho(ρ::AbstractMatrix{ComplexF64}, ρeq::AbstractMatrix{ComplexF64}, p_model)
    σ = sigma_from_rho(ρ, p_model)
    σeq = sigma_eq_from_rho(ρeq, p_model)
    return maximum(abs, σ .- σeq)
end

function build_rhoeq_guess_state(pstate, cfg::RunConfig)
    ρ = ρ_eq(0.0, cfg.β, pstate.p_model.H_ab, pstate.p_model.N_λ2, pstate.p_model.Nx, pstate.p_model.Ny, pstate.p_model.Nσ, pstate.p_model.N_orb)
    u = copy(pstate.u0)
    ptr = pointer_blocks(u, pstate.p_blocks.dims_ρ_ab, pstate.p_blocks.aux_layout)
    ptr.ρ_ab .= ρ
    return u, ρ
end

function rhs_sector_metrics(du::Vector{ComplexF64}, p_blocks)
    dptr = pointer_blocks(du, p_blocks.dims_ρ_ab, p_blocks.aux_layout)
    ρ_metrics = Dict(
        "rhs_rho_absmax" => maximum(abs, dptr.ρ_ab),
        "rhs_rho_fro" => norm(dptr.ρ_ab),
        "rhs_rho_trace_real" => real(tr(dptr.ρ_ab)),
        "rhs_rho_trace_imag" => imag(tr(dptr.ρ_ab)),
    )

    psi_max = 0.0
    psi_fro_sq = 0.0
    omega11_max = 0.0
    omega12_max = 0.0
    omega21_max = 0.0
    omega_fro_sq = 0.0

    for bptr in dptr.blocks
        psi_max = max(psi_max, maximum(abs, bptr.Ψ_anλ))
        psi_fro_sq += norm(bptr.Ψ_anλ)^2
    end
    for pairptr in dptr.Ω_pairs
        omega11_max = max(omega11_max, maximum(abs, pairptr.Ω11))
        omega12_max = max(omega12_max, maximum(abs, pairptr.Ω12))
        omega21_max = max(omega21_max, maximum(abs, pairptr.Ω21))
        omega_fro_sq += norm(pairptr.Ω11)^2 + norm(pairptr.Ω12)^2 + norm(pairptr.Ω21)^2
    end

    return merge(
        ρ_metrics,
        Dict(
            "rhs_psi_absmax" => psi_max,
            "rhs_psi_fro" => sqrt(psi_fro_sq),
            "rhs_omega11_absmax" => omega11_max,
            "rhs_omega12_absmax" => omega12_max,
            "rhs_omega21_absmax" => omega21_max,
            "rhs_omega_fro" => sqrt(omega_fro_sq),
            "rhs_total_absmax" => maximum(abs, du),
            "rhs_total_fro" => norm(du),
        ),
    )
end

function propagate_fixed_case(cfg::RunConfig, pstate, u_init::Vector{ComplexF64}, ρeq::Matrix{ComplexF64}, tag::AbstractString; progress_path::AbstractString="", wall_start::Float64=time())
    prob = ODEProblem(eom_tdnegf_blocks!, copy(u_init), (cfg.t0, cfg.t_end), pstate.p_blocks)
    intg = init(
        prob, Vern7();
        dt=cfg.dt, save_everystep=false, adaptive=false, dense=false,
        reltol=cfg.reltol, abstol=cfg.abstol,
    )

    n_steps = Int(round((cfg.t_end - cfg.t0) / cfg.dt))
    sigma_delta_trace = zeros(Float64, n_steps)
    rho_trace = zeros(Float64, n_steps)
    t_trace = zeros(Float64, n_steps)
    progress_every = max(1, env_int("FERMILINK_PROGRESS_EVERY_STEPS", 250))

    if !isempty(progress_path)
        write_progress_json(progress_path, cfg, tag, 0, n_steps, cfg.t0, wall_start, "running")
    end

    for step_idx in 1:n_steps
        DifferentialEquations.step!(intg, cfg.dt, true)
        ptr = pointer_blocks(intg.u, pstate.p_blocks.dims_ρ_ab, pstate.p_blocks.aux_layout)
        t_trace[step_idx] = intg.t
        sigma_delta_trace[step_idx] = sigma_delta_absmax_from_rho(ptr.ρ_ab, ρeq, pstate.p_model)
        rho_trace[step_idx] = maximum(abs, ptr.ρ_ab .- ρeq)

        if !isempty(progress_path) && (step_idx == n_steps || step_idx % progress_every == 0)
            status = tag == "zero_init_relax" && step_idx == n_steps ? "postprocessing" : "running"
            write_progress_json(progress_path, cfg, tag, step_idx, n_steps, intg.t, wall_start, status)
        end
    end

    ptr_final = pointer_blocks(intg.u, pstate.p_blocks.dims_ρ_ab, pstate.p_blocks.aux_layout)
    return Dict{String,Any}(
        "tag" => String(tag),
        "t_final" => intg.t,
        "sigma_delta_absmax_initial" => sigma_delta_trace[1],
        "sigma_delta_absmax_peak" => maximum(sigma_delta_trace),
        "sigma_delta_absmax_final" => sigma_delta_trace[end],
        "sigma_delta_absmax_mean" => sum(sigma_delta_trace) / length(sigma_delta_trace),
        "rho_minus_rhoeq_absmax_initial" => rho_trace[1],
        "rho_minus_rhoeq_absmax_peak" => maximum(rho_trace),
        "rho_minus_rhoeq_absmax_final" => rho_trace[end],
        "rho_minus_rhoeq_fro_final" => norm(ptr_final.ρ_ab .- ρeq),
        "rho_trace_real_final" => real(tr(ptr_final.ρ_ab)),
        "rho_trace_imag_final" => imag(tr(ptr_final.ρ_ab)),
        "final_state_path" => "",
    )
end

function build_summary(cfg::RunConfig, rhs_metrics::AbstractDict{String}, rhoeq_guess_metrics::AbstractDict{String}, zero_init_metrics::AbstractDict{String}, ρeq::Matrix{ComplexF64}, p_model)
    σeq = sigma_eq_from_rho(ρeq, p_model)
    summary = Pair{String,Any}[
        "run_id" => cfg.run_id,
        "status" => "done",
        "results_dir" => results_dir(cfg),
        "probe_mode" => "rhoeq_stationarity_check",
        "branch_reference_run_id" => "goal6_A_main_0003",
        "initial_configuration" => cfg.initial_configuration,
        "Hs_update_mode" => cfg.Hs_update_mode,
        "lead_bias_forced_zero" => true,
        "fixed_spin_texture" => true,
        "gamma_so" => cfg.γso,
        "rhoeq_pade_poles" => cfg.N_λ2,
        "Nx" => cfg.Nx,
        "Ny" => cfg.Ny,
        "j_sd" => cfg.j_sd,
        "j_ex" => cfg.j_ex,
        "j_ani" => cfg.j_ani,
        "damping" => cfg.damping,
        "beta" => cfg.β,
        "dt" => cfg.dt,
        "t_end" => cfg.t_end,
        "rhoeq_trace_real" => real(tr(ρeq)),
        "rhoeq_trace_imag" => imag(tr(ρeq)),
        "rhoeq_sigma_eq_absmax_static" => maximum(abs, σeq),
    ]
    append!(summary, collect(rhs_metrics))
    for (prefix, metrics) in (("rhoeq_guess", rhoeq_guess_metrics), ("zero_init", zero_init_metrics))
        for (key, value) in metrics
            key == "tag" && continue
            push!(summary, "$(prefix)_$(key)" => value)
        end
    end
    push!(summary, "rhoeq_guess_minus_zero_init_final_sigma_delta_absmax" => rhoeq_guess_metrics["sigma_delta_absmax_final"] - zero_init_metrics["sigma_delta_absmax_final"])
    push!(summary, "rhoeq_guess_minus_zero_init_final_rho_absmax" => rhoeq_guess_metrics["rho_minus_rhoeq_absmax_final"] - zero_init_metrics["rho_minus_rhoeq_absmax_final"])
    return summary
end

function main(args)
    cfg = parse_cli_args(RunConfig(; results_root=joinpath(@__DIR__, "results", "stationarity_cases")), args)
    case_dir = results_dir(cfg)
    mkpath(case_dir)
    progress_path = joinpath(case_dir, "progress.json")
    write_params_json(joinpath(case_dir, "params.json"), cfg)

    try
        wall_start = time()
        pstate = init_params_blocks(cfg)
        sys = init_params_heisenberg_llg(cfg)
        site_ranges_couple = [get_sub(i, pstate.p_model.N_loc) for i in 1:pstate.p_model.N_sites]
        update_H_e!(pstate.p_model, site_ranges_couple, sys.dipoles[:, :, 1, 1], cfg.j_sd)
        pstate.p_blocks.Δ_blocks[1] = 0.0 + 0.0im
        pstate.p_blocks.Δ_blocks[2] = 0.0 + 0.0im

        u_guess, ρeq = build_rhoeq_guess_state(pstate, cfg)
        du_guess = similar(u_guess)
        eom_tdnegf_blocks!(du_guess, u_guess, pstate.p_blocks, cfg.t0)
        rhs_metrics = rhs_sector_metrics(du_guess, pstate.p_blocks)

        rhoeq_guess_metrics = propagate_fixed_case(cfg, pstate, u_guess, ρeq, "rhoeq_guess_relax"; progress_path=progress_path, wall_start=wall_start)
        zero_init_metrics = propagate_fixed_case(cfg, pstate, pstate.u0, ρeq, "zero_init_relax"; progress_path=progress_path, wall_start=wall_start)

        summary = build_summary(cfg, rhs_metrics, rhoeq_guess_metrics, zero_init_metrics, ρeq, pstate.p_model)
        write_json_object(joinpath(case_dir, "summary.json"), summary)
        write_progress_json(progress_path, cfg, "done", Int(round((cfg.t_end - cfg.t0) / cfg.dt)), Int(round((cfg.t_end - cfg.t0) / cfg.dt)), cfg.t_end, wall_start, "completed")
        println("Stationarity check complete. Results written to: ", case_dir)
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
