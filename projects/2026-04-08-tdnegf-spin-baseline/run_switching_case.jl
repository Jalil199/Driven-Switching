tdnegf_project_root() = abspath(joinpath(@__DIR__, "..", "..", "..", "TDNEGF"))

using Pkg
Pkg.activate(tdnegf_project_root())

using TDNEGF
using DifferentialEquations
using LinearAlgebra
using Sunny
using DelimitedFiles
using Dates

default_run_id() = "run_" * Dates.format(Dates.now(Dates.UTC), "yyyymmddTHHMMSS")

Base.@kwdef mutable struct RunConfig
    run_id::String = default_run_id()
    initial_configuration::String = "FM"
    Hs_update_mode::String = "nonequilibrium_only"
    analysis_reference_mode::String = "static_rhoeq"
    Nx::Int = 3
    Ny::Int = 2
    Nσ::Int = 2
    N_orb::Int = 1
    γ::Float64 = 1.0
    γso::Float64 = 0.5
    Nα::Int = 2
    N_λ1::Int = 49
    N_λ2::Int = 20
    β::Float64 = 33.0
    dt::Float64 = 0.1
    t0::Float64 = 0.0
    t_end::Float64 = 200.0
    j_sd::Float64 = 0.1
    j_ex::Float64 = -0.05
    j_ani::Float64 = 0.1
    Bx::Float64 = 0.0
    a0::Float64 = 1.0
    damping::Float64 = 0.5
    kT::Float64 = 0.0
    bias_amplitude::Float64 = 1.0
    left_bias_scale::Float64 = 1.0
    right_bias_scale::Float64 = -1.0
    pulse_start::Float64 = 100.0
    pulse_ramp::Float64 = 10.0
    pulse_width::Float64 = 1000.0
    equilibration_time::Float64 = 0.0
    reltol::Float64 = 1e-6
    abstol::Float64 = 1e-8
    results_root::String = joinpath(@__DIR__, "results", "switching_runs")
end

results_dir(cfg::RunConfig) = joinpath(cfg.results_root, cfg.run_id)

function validated_hs_update_mode(mode::AbstractString)
    mode in ("nonequilibrium_only", "full_spin_density") || error(
        "Unknown Hs_update_mode=$(mode). Supported values: nonequilibrium_only, full_spin_density",
    )
    return String(mode)
end

function validated_analysis_reference_mode(mode::AbstractString)
    mode in ("static_rhoeq", "prepared_state") || error(
        "Unknown analysis_reference_mode=$(mode). Supported values: static_rhoeq, prepared_state",
    )
    return String(mode)
end

@inline linear_index(i::Int, j::Int, Ny::Int) = (i - 1) * Ny + j

@inline function ij_from_linear(l::Int, Ny::Int)
    i = (l - 1) ÷ Ny + 1
    j = (l - 1) % Ny + 1
    return i, j
end

function smooth_pulse(t::Float64; ti::Float64=10.0, Δ::Float64=100.0)
    t < 0.0 && return 0.0
    t < ti && return sin((pi / 2) * t / ti)^2
    t < Δ && return 1.0
    t < Δ + ti && return cos((pi / 2) * (t - Δ) / ti)^2
    return 0.0
end

function init_params_blocks(cfg::RunConfig)
    Rλ, zλ = load_poles_square(cfg.N_λ1, cfg.N_λ2)
    p_model = ModelParamsTDNEGF(
        Nx=cfg.Nx, Ny=cfg.Ny, Nσ=cfg.Nσ, N_orb=cfg.N_orb,
        Nα=cfg.Nα, N_λ1=cfg.N_λ1, N_λ2=cfg.N_λ2,
    )

    H_ab = build_H_ab(
        Nx=p_model.Nx, Ny=p_model.Ny, Nσ=p_model.Nσ, N_orb=p_model.N_orb,
        γ=cfg.γ, γso=cfg.γso + 0.0im,
    )
    ΣL_nλ = build_Σᴸ_nλ(Rλ, zλ, p_model.Ny, p_model.Nσ, p_model.N_orb, p_model.N_λ1, p_model.N_λ2; β=cfg.β, γ=1.0)
    ΣG_nλ = build_Σᴳ_nλ(Rλ, zλ, p_model.Ny, p_model.Nσ, p_model.N_orb, p_model.N_λ1, p_model.N_λ2; β=cfg.β, γ=1.0)
    χ_nλ = build_χ_nλ(zλ, p_model.Ny, p_model.Nσ, p_model.N_orb, p_model.N_λ1, p_model.N_λ2; β=cfg.β, γ=1.0)
    ξ_anR = build_ξ_an(p_model.Nx, p_model.Ny, p_model.Nσ, p_model.N_orb; xcol=p_model.Nx, y_coup=1:p_model.Ny)
    ξ_anL = build_ξ_an(p_model.Nx, p_model.Ny, p_model.Nσ, p_model.N_orb; xcol=1, y_coup=1:p_model.Ny)

    left_block = SelfEnergyBlock(:left, p_model.Nc, p_model.N_λ1, p_model.N_λ2, ΣL_nλ, ΣG_nλ, χ_nλ, ξ_anL)
    right_block = SelfEnergyBlock(:right, p_model.Nc, p_model.N_λ1, p_model.N_λ2, ΣL_nλ, ΣG_nλ, χ_nλ, ξ_anR)
    blocks = [left_block, right_block]

    p_model.H_ab .= H_ab
    p_model.H0_ab .= H_ab
    p_blocks = ExperimentalBlockRHSParams(p_model.H_ab, blocks, ComplexF64[0.0 + 0.0im, 0.0 + 0.0im], p_model)
    u0 = zeros(ComplexF64, p_blocks.dims_ρ_ab[1]^2 + p_blocks.aux_layout.total_size)
    return (p_model=p_model, p_blocks=p_blocks, u0=u0)
end

function init_params_heisenberg_llg(cfg::RunConfig)
    latvecs = lattice_vectors(cfg.a0, cfg.a0, 4 * cfg.a0, 90, 90, 90)
    positions = [[0.5, 0.5, 0.0]]
    cryst = Crystal(latvecs, positions)
    moments = [1 => Moment(s=1.0, g=1.0)]
    sys = System(cryst, moments, :dipole; dims=(cfg.Nx, cfg.Ny, 1))
    set_field!(sys, [cfg.Bx, 0.0, 0.0])
    set_exchange!(sys, cfg.j_ex, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, cfg.j_ex, Bond(1, 1, [0, 1, 0]))
    set_onsite_coupling!(sys, S -> -cfg.j_ani * S[3]^2, 1)
    sys = to_inhomogeneous(sys)
    remove_periodicity!(sys, (true, true, true))

    init_mode = uppercase(cfg.initial_configuration)
    if init_mode == "AFM"
        for I in eachindex(sys.dipoles)
            i, j = ij_from_linear(I, cfg.Ny)
            sys.dipoles[I] = Sunny.SVector(0.0, 0.0, (-1)^(i + j))
        end
    elseif init_mode == "FM"
        for I in eachindex(sys.dipoles)
            sys.dipoles[I] = Sunny.SVector(0.0, 0.0, 1.0)
        end
    else
        error("Unknown initial_configuration=$(cfg.initial_configuration). Supported values: FM, AFM")
    end
    return sys
end

function update_H_s!(sys::System{0}, σx_i::Array{Float64,2}, j_sd::Float64, Nx::Int, Ny::Int; B0::Vector{Float64}=[0.0, 0.0, 0.0])
    @inbounds for i in 1:Nx, j in 1:Ny
        idx = linear_index(i, j, Ny)
        sden = -σx_i[idx, :] * j_sd
        set_field_at!(sys, sden .+ B0, (i, j, 1, 1))
    end
    return nothing
end

function component_absmax(A)
    return maximum(abs, A)
end

function pre_pulse_mask(t::AbstractVector{<:Real}, pulse_start::Real)
    return [ti < pulse_start for ti in t]
end

function masked_absmax(A, mask)
    any(mask) || return 0.0
    return maximum(abs, @view A[:, :, mask])
end

function avg_M(dip::Array{Float64,3}, Nx::Int, Ny::Int)
    size_t = size(dip, 3)
    size_sys = Nx * Ny
    M_t = zeros(Float64, 3, size_t)
    for t_idx in 1:size_t, i in 1:size_sys, x in 1:3
        M_t[x, t_idx] += dip[i, x, t_idx]
    end
    return M_t / size_sys
end

function avg_N(dip::Array{Float64,3}, Nx::Int, Ny::Int)
    size_t = size(dip, 3)
    size_sys = Nx * Ny
    N_t = zeros(Float64, 3, size_t)
    for t_idx in 1:size_t, i in 1:size_sys
        i_x, j_y = ij_from_linear(i, Ny)
        for x in 1:3
            N_t[x, t_idx] += dip[i, x, t_idx] * (-1)^(i_x + j_y)
        end
    end
    return N_t / size_sys
end

function current_static_rhoeq(cfg::RunConfig, p_model)
    return ρ_eq(0.0, cfg.β, p_model.H_ab, p_model.N_λ2, p_model.Nx, p_model.Ny, p_model.Nσ, p_model.N_orb)
end

function current_analysis_reference(cfg::RunConfig, p_model, prepared_ρ_ref)
    if cfg.analysis_reference_mode == "prepared_state"
        prepared_ρ_ref === nothing && error("analysis_reference_mode=prepared_state requires a prepared reference density")
        return prepared_ρ_ref
    end
    return current_static_rhoeq(cfg, p_model)
end

function llg_drive_density(cfg::RunConfig, dv, p_model, prepared_ρ_ref; force_static_ref::Bool=false)
    obs = ObservablesTDNEGF(p_model; N_tmax=1, N_leads=length((1, 2)))
    obs.idx = 1
    obs_σ_i!(dv, p_model, obs)
    σx_full = copy(obs.σx_i[:, :, 1])
    σx_ref = begin
        ρref = force_static_ref ? current_static_rhoeq(cfg, p_model) : current_analysis_reference(cfg, p_model, prepared_ρ_ref)
        ref_obs = ObservablesTDNEGF(p_model; N_tmax=1, N_leads=length((1, 2)))
        ref_obs.idx = 1
        obs_σ_i_eq!(ρref, p_model, ref_obs)
        copy(ref_obs.σx_i_eq[:, :, 1])
    end
    return cfg.Hs_update_mode == "full_spin_density" ? σx_full : (σx_full .- σx_ref)
end

function prepare_equilibrated_state!(cfg::RunConfig, p_model, p_blocks, sys, u)
    cfg.equilibration_time <= 0.0 && return nothing

    site_ranges_couple = [get_sub(i, p_model.N_loc) for i in 1:p_model.N_sites]
    p_blocks.Δ_blocks[1] = 0.0 + 0.0im
    p_blocks.Δ_blocks[2] = 0.0 + 0.0im
    update_H_e!(p_model, site_ranges_couple, sys.dipoles[:, :, 1, 1], cfg.j_sd)

    prob = ODEProblem(eom_tdnegf_blocks!, u, (0.0, cfg.equilibration_time), p_blocks)
    intg = init(
        prob, Vern7();
        dt=cfg.dt, save_everystep=false, adaptive=true, dense=false,
        reltol=cfg.reltol, abstol=cfg.abstol,
    )
    integ = Langevin(cfg.dt; damping=cfg.damping, kT=cfg.kT)

    n_steps = Int(round(cfg.equilibration_time / cfg.dt))
    for _ in 1:n_steps
        DifferentialEquations.step!(intg, cfg.dt, true)
        Sunny.step!(sys, integ)

        dv = pointer_blocks(intg.u, p_blocks.dims_ρ_ab, p_blocks.aux_layout)
        hs_drive = llg_drive_density(cfg, dv, p_model, nothing; force_static_ref=true)
        update_H_s!(sys, hs_drive, cfg.j_sd, p_model.Nx, p_model.Ny)
        update_H_e!(p_model, site_ranges_couple, sys.dipoles[:, :, 1, 1], cfg.j_sd)
        p_blocks.Δ_blocks[1] = 0.0 + 0.0im
        p_blocks.Δ_blocks[2] = 0.0 + 0.0im
    end

    return copy(intg.u)
end

function run_switching_case(cfg::RunConfig)
    cfg.Hs_update_mode = validated_hs_update_mode(cfg.Hs_update_mode)
    cfg.analysis_reference_mode = validated_analysis_reference_mode(cfg.analysis_reference_mode)
    cfg.analysis_reference_mode == "prepared_state" && cfg.equilibration_time <= 0.0 && error(
        "analysis_reference_mode=prepared_state requires equilibration_time > 0",
    )
    pstate = init_params_blocks(cfg)
    p_model = pstate.p_model
    p_blocks = pstate.p_blocks
    u = copy(pstate.u0)
    sys = init_params_heisenberg_llg(cfg)

    prepared_state = prepare_equilibrated_state!(cfg, p_model, p_blocks, sys, u)
    if prepared_state !== nothing
        u .= prepared_state
    end
    prepared_ptr = pointer_blocks(u, p_blocks.dims_ρ_ab, p_blocks.aux_layout)
    prepared_ρ_ref = copy(prepared_ptr.ρ_ab)

    prob = ODEProblem(eom_tdnegf_blocks!, u, (cfg.t0, cfg.t_end), p_blocks)
    intg = init(
        prob, Vern7();
        dt=cfg.dt, save_everystep=false, adaptive=true, dense=false,
        reltol=cfg.reltol, abstol=cfg.abstol,
    )
    integ = Langevin(cfg.dt; damping=cfg.damping, kT=cfg.kT)

    N_steps = Int(round((cfg.t_end - cfg.t0) / cfg.dt))
    obs = ObservablesTDNEGF(p_model; N_tmax=N_steps, N_leads=length(p_blocks.blocks))
    site_ranges_couple = [get_sub(i, p_model.N_loc) for i in 1:p_model.N_sites]

    for step_idx in 1:N_steps
        obs.idx = step_idx
        DifferentialEquations.step!(intg, cfg.dt, true)
        Sunny.step!(sys, integ)

        dv = pointer_blocks(intg.u, p_blocks.dims_ρ_ab, p_blocks.aux_layout)
        ρ = current_analysis_reference(cfg, p_model, prepared_ρ_ref)

        obs.t[step_idx] = intg.t
        obs_n_i!(dv, p_model, obs)
        obs_σ_i!(dv, p_model, obs)
        obs_Ixα!(dv, p_blocks, obs)
        obs_s_i!(sys.dipoles[:, :, 1, 1], p_model, obs)
        obs_σ_i_eq!(ρ, p_model, obs)

        σx_full = obs.σx_i[:, :, step_idx]
        σx_delta = σx_full .- obs.σx_i_eq[:, :, step_idx]
        hs_drive = cfg.Hs_update_mode == "full_spin_density" ? σx_full : σx_delta
        update_H_s!(sys, hs_drive, cfg.j_sd, p_model.Nx, p_model.Ny)
        update_H_e!(p_model, site_ranges_couple, sys.dipoles[:, :, 1, 1], cfg.j_sd)
        lead_bias = smooth_pulse(intg.t - cfg.pulse_start; ti=cfg.pulse_ramp, Δ=cfg.pulse_width) * cfg.bias_amplitude + 0.0im
        p_blocks.Δ_blocks[1] = cfg.left_bias_scale * lead_bias
        p_blocks.Δ_blocks[2] = cfg.right_bias_scale * lead_bias
    end

    return (; obs, final_state=intg.u, sys, p_model, p_blocks, prepared_ρ_ref)
end

json_escape(s::AbstractString) = replace(replace(replace(replace(String(s), "\\" => "\\\\"), "\"" => "\\\""), "\n" => "\\n"), "\t" => "\\t")

function json_value(x)
    if x isa AbstractString
        return "\"" * json_escape(x) * "\""
    elseif x isa Bool
        return x ? "true" : "false"
    elseif x isa Integer
        return string(x)
    elseif x isa AbstractFloat
        return isfinite(x) ? string(x) : "null"
    elseif x isa AbstractVector
        return "[" * join(json_value.(collect(x)), ", ") * "]"
    else
        return "\"" * json_escape(string(x)) * "\""
    end
end

function write_json_object(path::AbstractString, pairs::Vector{Pair{String,Any}})
    open(path, "w") do io
        println(io, "{")
        for (idx, kv) in enumerate(pairs)
            suffix = idx == length(pairs) ? "" : ","
            println(io, "  ", json_value(kv.first), ": ", json_value(kv.second), suffix)
        end
        println(io, "}")
    end
end

function write_params_json(path::AbstractString, cfg::RunConfig)
    pairs = Pair{String,Any}[
        "run_id" => cfg.run_id,
        "initial_configuration" => cfg.initial_configuration,
        "Hs_update_mode" => cfg.Hs_update_mode,
        "analysis_reference_mode" => cfg.analysis_reference_mode,
        "Nx" => cfg.Nx,
        "Ny" => cfg.Ny,
        "Nσ" => cfg.Nσ,
        "N_orb" => cfg.N_orb,
        "γ" => cfg.γ,
        "γso" => cfg.γso,
        "Nα" => cfg.Nα,
        "N_λ1" => cfg.N_λ1,
        "N_λ2" => cfg.N_λ2,
        "β" => cfg.β,
        "dt" => cfg.dt,
        "t0" => cfg.t0,
        "t_end" => cfg.t_end,
        "j_sd" => cfg.j_sd,
        "j_ex" => cfg.j_ex,
        "j_ani" => cfg.j_ani,
        "Bx" => cfg.Bx,
        "a0" => cfg.a0,
        "damping" => cfg.damping,
        "kT" => cfg.kT,
        "bias_amplitude" => cfg.bias_amplitude,
        "left_bias_scale" => cfg.left_bias_scale,
        "right_bias_scale" => cfg.right_bias_scale,
        "pulse_start" => cfg.pulse_start,
        "pulse_ramp" => cfg.pulse_ramp,
        "pulse_width" => cfg.pulse_width,
        "equilibration_time" => cfg.equilibration_time,
        "reltol" => cfg.reltol,
        "abstol" => cfg.abstol,
        "results_root" => cfg.results_root,
        "results_dir" => results_dir(cfg),
    ]
    write_json_object(path, pairs)
end

function write_observables(path::AbstractString, obs, cfg::RunConfig)
    headers = String[
        "t",
        "I_left", "I_right",
        "I_left_sx", "I_left_sy", "I_left_sz",
        "I_right_sx", "I_right_sy", "I_right_sz",
        "avg_Mx", "avg_My", "avg_Mz",
        "neel_x", "neel_y", "neel_z",
    ]
    for site in 1:(cfg.Nx * cfg.Ny)
        push!(headers, "n_site$(site)")
        push!(headers, "sigma_x_site$(site)")
        push!(headers, "sigma_y_site$(site)")
        push!(headers, "sigma_z_site$(site)")
        push!(headers, "sigma_eq_x_site$(site)")
        push!(headers, "sigma_eq_y_site$(site)")
        push!(headers, "sigma_eq_z_site$(site)")
        push!(headers, "sigma_delta_x_site$(site)")
        push!(headers, "sigma_delta_y_site$(site)")
        push!(headers, "sigma_delta_z_site$(site)")
        push!(headers, "s_classical_x_site$(site)")
        push!(headers, "s_classical_y_site$(site)")
        push!(headers, "s_classical_z_site$(site)")
    end

    M = avg_M(obs.sx_i, cfg.Nx, cfg.Ny)
    N = avg_N(obs.sx_i, cfg.Nx, cfg.Ny)
    data = Matrix{Float64}(undef, length(obs.t), length(headers))
    for it in eachindex(obs.t)
        col = 1
        data[it, col] = obs.t[it]; col += 1
        data[it, col] = obs.Iα[1, it]; col += 1
        data[it, col] = obs.Iα[2, it]; col += 1
        data[it, col] = obs.Iαx[1, 1, it]; col += 1
        data[it, col] = obs.Iαx[1, 2, it]; col += 1
        data[it, col] = obs.Iαx[1, 3, it]; col += 1
        data[it, col] = obs.Iαx[2, 1, it]; col += 1
        data[it, col] = obs.Iαx[2, 2, it]; col += 1
        data[it, col] = obs.Iαx[2, 3, it]; col += 1
        data[it, col] = M[1, it]; col += 1
        data[it, col] = M[2, it]; col += 1
        data[it, col] = M[3, it]; col += 1
        data[it, col] = N[1, it]; col += 1
        data[it, col] = N[2, it]; col += 1
        data[it, col] = N[3, it]; col += 1
        for site in 1:(cfg.Nx * cfg.Ny)
            data[it, col] = obs.n_i[site, it]; col += 1
            data[it, col] = obs.σx_i[site, 1, it]; col += 1
            data[it, col] = obs.σx_i[site, 2, it]; col += 1
            data[it, col] = obs.σx_i[site, 3, it]; col += 1
            data[it, col] = obs.σx_i_eq[site, 1, it]; col += 1
            data[it, col] = obs.σx_i_eq[site, 2, it]; col += 1
            data[it, col] = obs.σx_i_eq[site, 3, it]; col += 1
            data[it, col] = obs.σx_i[site, 1, it] - obs.σx_i_eq[site, 1, it]; col += 1
            data[it, col] = obs.σx_i[site, 2, it] - obs.σx_i_eq[site, 2, it]; col += 1
            data[it, col] = obs.σx_i[site, 3, it] - obs.σx_i_eq[site, 3, it]; col += 1
            data[it, col] = obs.sx_i[site, 1, it]; col += 1
            data[it, col] = obs.sx_i[site, 2, it]; col += 1
            data[it, col] = obs.sx_i[site, 3, it]; col += 1
        end
    end

    open(path, "w") do io
        println(io, join(headers, '\t'))
        writedlm(io, data, '\t')
    end
end

function summarize_run(result, cfg::RunConfig)
    obs = result.obs
    M = avg_M(obs.sx_i, cfg.Nx, cfg.Ny)
    N = avg_N(obs.sx_i, cfg.Nx, cfg.Ny)
    current_total = obs.Iα[1, :] .+ obs.Iα[2, :]
    σ_delta = obs.σx_i .- obs.σx_i_eq
    avgMz = M[3, :]
    neelz = N[3, :]
    current_peak = maximum(abs.(vec(obs.Iα)))
    current_balance_peak = maximum(abs.(current_total))
    avgMz_initial = avgMz[1]
    avgMz_final = avgMz[end]
    avgMz_min = minimum(avgMz)
    avgMz_max = maximum(avgMz)
    neelz_peak = maximum(abs.(neelz))
    reversed = sign(avgMz_initial) != 0 && sign(avgMz_final) != sign(avgMz_initial)
    partial_reversal = avgMz_initial != 0.0 && avgMz_final / avgMz_initial < 0.5
    pre_mask = pre_pulse_mask(obs.t, cfg.pulse_start)
    pre_pulse_samples = count(pre_mask)

    return Pair{String,Any}[
        "run_id" => cfg.run_id,
        "initial_configuration" => cfg.initial_configuration,
        "Hs_update_mode" => cfg.Hs_update_mode,
        "analysis_reference_mode" => cfg.analysis_reference_mode,
        "results_dir" => results_dir(cfg),
        "status" => "done",
        "completed_steps" => length(obs.t),
        "t_final" => obs.t[end],
        "j_sd" => cfg.j_sd,
        "j_ex" => cfg.j_ex,
        "j_ani" => cfg.j_ani,
        "Bx" => cfg.Bx,
        "pulse_amp" => cfg.bias_amplitude,
        "bias_amplitude" => cfg.bias_amplitude,
        "left_bias_scale" => cfg.left_bias_scale,
        "right_bias_scale" => cfg.right_bias_scale,
        "gamma_so" => cfg.γso,
        "γso" => cfg.γso,
        "pulse_start" => cfg.pulse_start,
        "pulse_ramp" => cfg.pulse_ramp,
        "pulse_width" => cfg.pulse_width,
        "equilibration_time" => cfg.equilibration_time,
        "current_peak_abs" => current_peak,
        "current_left_final" => obs.Iα[1, end],
        "current_right_final" => obs.Iα[2, end],
        "current_balance_peak_abs" => current_balance_peak,
        "avg_M_final" => vec(M[:, end]),
        "avg_Mz_initial" => avgMz_initial,
        "avg_Mz_final" => avgMz_final,
        "avg_Mz_min" => avgMz_min,
        "avg_Mz_max" => avgMz_max,
        "neel_final" => vec(N[:, end]),
        "neel_z_peak_abs" => neelz_peak,
        "full_reversal_detected" => reversed,
        "partial_reversal_flag" => partial_reversal,
        "max_classical_spin_norm_deviation" => maximum(abs.(sqrt.(sum(obs.sx_i .^ 2; dims=2))[:, 1, :] .- 1.0)),
        "max_site_charge_deviation_from_1" => maximum(abs.(obs.n_i .- 1.0)),
        "max_sigma_noneq_abs" => component_absmax(obs.σx_i),
        "max_sigma_eq_abs" => component_absmax(obs.σx_i_eq),
        "max_sigma_delta_abs" => component_absmax(σ_delta),
        "pre_pulse_sample_count" => pre_pulse_samples,
        "pre_pulse_duration_realized" => pre_pulse_samples == 0 ? 0.0 : obs.t[findlast(pre_mask)] - obs.t[1],
        "pre_pulse_sigma_noneq_absmax" => masked_absmax(obs.σx_i, pre_mask),
        "pre_pulse_sigma_eq_absmax" => masked_absmax(obs.σx_i_eq, pre_mask),
        "pre_pulse_sigma_delta_absmax" => masked_absmax(σ_delta, pre_mask),
        "llg_drive_convention" => cfg.Hs_update_mode == "full_spin_density" ? "sigma_full" : "sigma_delta_noneq_minus_eq",
        "analysis_reference_saved" => cfg.analysis_reference_mode,
        "prepared_state_used" => cfg.equilibration_time > 0.0,
        "equilibrium_reference_saved" => true,
        "nonequilibrium_minus_equilibrium_saved" => true,
        "lead_bias_protocol" => cfg.left_bias_scale == 0.0 && cfg.right_bias_scale == 1.0 ? "right_only" : "scaled_two_lead",
        "timestamp_utc" => string(Dates.now(Dates.UTC)),
    ]
end

function failure_summary(cfg::RunConfig, err)
    return Pair{String,Any}[
        "run_id" => cfg.run_id,
        "initial_configuration" => cfg.initial_configuration,
        "Hs_update_mode" => cfg.Hs_update_mode,
        "analysis_reference_mode" => cfg.analysis_reference_mode,
        "results_dir" => results_dir(cfg),
        "status" => "failed",
        "j_sd" => cfg.j_sd,
        "j_ex" => cfg.j_ex,
        "j_ani" => cfg.j_ani,
        "bias_amplitude" => cfg.bias_amplitude,
        "left_bias_scale" => cfg.left_bias_scale,
        "right_bias_scale" => cfg.right_bias_scale,
        "pulse_start" => cfg.pulse_start,
        "pulse_ramp" => cfg.pulse_ramp,
        "pulse_width" => cfg.pulse_width,
        "equilibration_time" => cfg.equilibration_time,
        "error" => sprint(showerror, err),
        "timestamp_utc" => string(Dates.now(Dates.UTC)),
    ]
end

summary_dict(pairs::Vector{Pair{String,Any}}) = Dict(k => v for (k, v) in pairs)

function run_switching_case_and_write(cfg::RunConfig)
    case_dir = results_dir(cfg)
    mkpath(case_dir)
    write_params_json(joinpath(case_dir, "params.json"), cfg)
    try
        result = run_switching_case(cfg)
        write_observables(joinpath(case_dir, "observables.tsv"), result.obs, cfg)
        summary = summarize_run(result, cfg)
        write_json_object(joinpath(case_dir, "summary.json"), summary)
        return summary_dict(summary)
    catch err
        summary = failure_summary(cfg, err)
        write_json_object(joinpath(case_dir, "summary.json"), summary)
        rethrow()
    end
end

function parse_cli_args(cfg::RunConfig, args::Vector{String})
    for arg in args
        occursin('=', arg) || error("Expected key=value argument, got: $arg")
        key, value = split(arg, '='; limit=2)
        if key in ("run_id", "initial_configuration", "results_root", "Hs_update_mode", "analysis_reference_mode")
            setfield!(cfg, Symbol(key), String(value))
        elseif key in ("Nx", "Ny", "Nσ", "N_orb", "Nα", "N_λ1", "N_λ2")
            setfield!(cfg, Symbol(key), parse(Int, value))
        elseif key in ("γ", "γso", "β", "dt", "t0", "t_end", "j_sd", "j_ex", "j_ani", "Bx", "a0", "damping", "kT", "bias_amplitude", "left_bias_scale", "right_bias_scale", "pulse_start", "pulse_ramp", "pulse_width", "equilibration_time", "reltol", "abstol")
            setfield!(cfg, Symbol(key), parse(Float64, value))
        else
            error("Unknown parameter: $key")
        end
    end
    return cfg
end

function main(args)
    cfg = parse_cli_args(RunConfig(), args)
    summary = run_switching_case_and_write(cfg)
    println("Run complete. Results written to: ", summary["results_dir"])
    println("Final time: ", get(summary, "t_final", missing))
    println("Final average Mz: ", get(summary, "avg_Mz_final", missing))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
