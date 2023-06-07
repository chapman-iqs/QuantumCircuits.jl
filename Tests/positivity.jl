ispositive(ψ::Ket; kwargs...) = ispositive(dm(ψ); kwargs...)
ispositive(ρs::Vector; kwargs...) = prod(ispositive.(ρs; kwargs...))
function ispositive(ρ::Operator; tolerance = 1e-6)
    return prod(eigvals(ρ.data) .<= 1 + tolerance)
end

function positive_trajectory(; solve=bayesian, tolerance=1e-6, kwargs...)

    Ω = 2π * 1.0
    Γ = 2.0
    η = 1.0

    ψ0 = normalize(g + e)
    tf = 10.0
    dt = 1e-3

    H = Ω/2 * σy
    J = η < 1.0 ? [(σz, (1 - η) * Γ/2)] : []
    C = [(σz, Γ, η)]

    sol = solve((0.0, tf), ψ0, H, J, C; dt = dt)
    sol = Solution(sol, qbasis)

    # return plot(blochtimeseries, ts, exps...)
    return ispositive(sol.ρ; tolerance=tolerance)
end

function test_positive_trajectory(; tolerance=1e-6, kwargs...)
    @testset "positive trajectory" begin
        @test positive_trajectory(; tolerance=tolerance, kwargs...)
    end
end

