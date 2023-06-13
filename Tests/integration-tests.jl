function returns_solution(; solve=bayesian, kwargs...)

    Ω = 2π * 1.0
    Γ = 2.0
    η = 0.5

    ψ0 = normalize(g + e)
    tf = 10.0
    dt = 1e-3

    H = Ω/2 * σy
    J = [(σz, (1 - η) * Γ/2)]
    C = [(σz, Γ, η)]

    sol = solve((0.0, tf), ψ0, H, J, C; dt = dt)

    # return plot(blochtimeseries, ts, exps...)
    return isa(sol, Solution)
end


function test_integration(; solve=bayesian, kwargs...)

    @testset "returns solution ($(solve))" begin
        returns_solution(; solve=solve, kwargs...)   
    end
end