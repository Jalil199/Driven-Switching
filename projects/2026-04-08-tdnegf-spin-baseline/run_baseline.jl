using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using TDNEGF
using DifferentialEquations
using LinearAlgebra
using Sunny
using DelimitedFiles
using Dates

Base.@kwdef mutable struct RunConfig
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
    pulse_start::Float64 = 100.0
    pulse_ramp::Float64 = 10.0
    pulse_width::Float64 = 1000.0
    reltol::Float64 = 1e-6
    abstol::Float64 = 1e-8
    results_dir::String = joinpath(@__DIR__, "results", "baseline_notebook_like")
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

    if cfg.j_ex >= 0.0
        for I in eachindex(sys.dipoles)
            i, j = ij_from_linear(I, cfg.Ny)
            sys.dipoles[I] = Sunny.SVector(0.0, 0.0, (-1)^(i + j))
        end
    else
        for I in eachindex(sys.dipoles)
            sys.dipoles[I] = Sunny.SVector(0.0, 0.0, 1.0)
        end
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

function run_spin_electron_baseline(cfg::RunConfig)
    pstate = init_params_blocks(cfg)
    p_model = pstate.p_model
    p_blocks = pstate.p_blocks
    u = copy(pstate.u0)
    sys = init_params_heisenberg_llg(cfg)

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
        ρ = ρ_eq(0.0, cfg.β, p_model.H_ab, p_model.N_λ2, p_model.Nx, p_model.Ny, p_model.Nσ, p_model.N_orb)

        obs.t[step_idx] = intg.t
        obs_n_i!(dv, p_model, obs)
        obs_σ_i!(dv, p_model, obs)
        obs_Ixα!(dv, p_blocks, obs)
        obs_s_i!(sys.dipoles[:, :, 1, 1], p_model, obs)
        obs_σ_i_eq!(ρ, p_model, obs)

        update_H_s!(sys, obs.σx_i[:, :, step_idx], cfg.j_sd, p_model.Nx, p_model.Ny)
        update_H_e!(p_model, site_ranges_couple, sys.dipoles[:, :, 1, 1], cfg.j_sd)
        Δ_left = smooth_pulse(intg.t - cfg.pulse_start; ti=cfg.pulse_ramp, Δ=cfg.pulse_width) * cfg.bias_amplitude + 0.0im
        p_blocks.Δ_blocks[1] = +Δ_left
        p_blocks.Δ_blocks[2] = -Δ_left
    end

    return (; obs, final_state=intg.u, sys, p_model, p_blocks)
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
        "pulse_start" => cfg.pulse_start,
        "pulse_ramp" => cfg.pulse_ramp,
        "pulse_width" => cfg.pulse_width,
        "reltol" => cfg.reltol,
        "abstol" => cfg.abstol,
        "results_dir" => cfg.results_dir,
    ]
    write_json_object(path, pairs)
end

function write_observables(path::AbstractString, obs, cfg::RunConfig)
    headers = String["t", "I_left", "I_right", "I_left_sx", "I_left_sy", "I_left_sz", "I_right_sx", "I_right_sy", "I_right_sz"]
    for site in 1:(cfg.Nx * cfg.Ny)
        push!(headers, "n_site$(site)")
        push!(headers, "sigma_z_site$(site)")
        push!(headers, "s_classical_z_site$(site)")
    end

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
        for site in 1:(cfg.Nx * cfg.Ny)
            data[it, col] = obs.n_i[site, it]; col += 1
            data[it, col] = obs.σx_i[site, 3, it]; col += 1
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

    return Pair{String,Any}[
        "status" => "done",
        "completed_steps" => length(obs.t),
        "t_final" => obs.t[end],
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
        "timestamp_utc" => string(Dates.now(Dates.UTC)),
    ]
end

function parse_cli_args(cfg::RunConfig, args::Vector{String})
    for arg in args
        occursin('=', arg) || error("Expected key=value argument, got: $arg")
        key, value = split(arg, '='; limit=2)
        if key == "results_dir"
            setfield!(cfg, :results_dir, String(value))
        elseif key in ("Nx", "Ny", "Nσ", "N_orb", "Nα", "N_λ1", "N_λ2")
            setfield!(cfg, Symbol(key), parse(Int, value))
        elseif key in ("γ", "γso", "β", "dt", "t0", "t_end", "j_sd", "j_ex", "j_ani", "Bx", "a0", "damping", "kT", "bias_amplitude", "pulse_start", "pulse_ramp", "pulse_width", "reltol", "abstol")
            setfield!(cfg, Symbol(key), parse(Float64, value))
        else
            error("Unknown parameter: $key")
        end
    end
    return cfg
end

function main(args)
    cfg = parse_cli_args(RunConfig(), args)
    mkpath(cfg.results_dir)
    write_params_json(joinpath(cfg.results_dir, "params.json"), cfg)
    result = run_spin_electron_baseline(cfg)
    write_observables(joinpath(cfg.results_dir, "observables.tsv"), result.obs, cfg)
    write_json_object(joinpath(cfg.results_dir, "summary.json"), summarize_run(result, cfg))
    println("Run complete. Results written to: ", cfg.results_dir)
    println("Final time: ", result.obs.t[end])
    println("Final average classical magnetization: ", vec(avg_M(result.obs.sx_i, cfg.Nx, cfg.Ny)[:, end]))
end

main(ARGS)
