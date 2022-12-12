function run_sbayesian(; verbose = true, kwargs...)

    ψ0 = dm(normalize(g + e))
    tf = 1.0
    dt = 1e-3
    Ω = 2π * 1.0
    Γ = 0.2
    γ = 0.5
    η = 0.8

    H = Ω/2 * σy
    J = [(σz, Γ/2), (σz, (1-η)γ/2)]
    C = [(σz, γ, η)]

    if verbose
        println("testing sbayesian...")
    end

    @time sol = sbayesian((0.0, tf), ψ0, H, J, C; dt=dt)
    @time sol = sbayesian((0.0, tf), ψ0, H, J, C; dt=dt)
    @time sol = sbayesian((0.0, tf), ψ0, H, J, C; dt=dt)
end

function run_ssbayesian(; verbose = true, kwargs...)

    ψ0 = dm(normalize(g + e))
    tf = 1.0
    dt = 1e-3
    Ω = 2π * 1.0
    Γ = 0.2
    γ = 0.5
    η = 0.8

    H = Ω/2 * σy
    J = [(σz, Γ/2), (σz, (1-η)γ/2)]
    C = [(σz, γ, η)]

    if verbose
        println("testing ssbayesian...")
    end

    @time sol = ssbayesian((0.0, tf), ψ0, H, J, C; dt=dt)
    @time sol = ssbayesian((0.0, tf), ψ0, H, J, C; dt=dt)
    @time sol = ssbayesian((0.0, tf), ψ0, H, J, C; dt=dt)
end

function test_sbayesian(; kwargs...)
    run_sbayesian(; kwargs...)
end

function test_ssbayesian(; kwargs...)
    run_ssbayesian(; kwargs...)
end