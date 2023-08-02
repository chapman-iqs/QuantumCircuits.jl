ispositive(ψ::Ket; kwargs...) = ispositive(dm(ψ); kwargs...)
ispositive(ρs::Vector; kwargs...) = prod(ispositive.(ρs; kwargs...))
# should it be necessary to take real part of eigenvalues? rouchon has complex eigvals
function ispositive(ρ::Operator; tolerance = 1e-6)
    evals = eigvals(ρ.data)
    isreal = prod(abs.(imag.(evals)) .<= tolerance)
    lessthanone = prod(real.(evals) .<= 1 + tolerance)
    greaterthanzero = prod(real.(evals) .> -1 * tolerance)
    return isreal * lessthanone * greaterthanzero
end

function positive_trajectory(; makeplots=false, testset_id="", test_id="positive_trajectory", plotpath="test_result_plots", solve=bayesian, tolerance=1e-6, kwargs...)

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
    pass = ispositive(sol.ρ; tolerance=tolerance)
    pass_string = pass ? "passed" : "failed"

    make_test_plot(makeplots, sol, solve, plotpath, testset_id, string(pass_string, "_", test_id))

    # return plot(blochtimeseries, ts, exps...)
    return pass
end

function maintains_purity(; makeplots=false, testset_id="", test_id="maintains_purity", plotpath="test_result_plots", solve=bayesian, tolerance=1e-6, kwargs...)

    Ω = 2π * 1.0
    Γ = 2.0
    η = 1.0

    ψ0 = dm(normalize(g + e))
    tf = 10.0
    dt = 1e-3

    H = Ω/2 * σy
    J = []
    C = [(σz, Γ, η)]

    sol = solve((0.0, tf), ψ0, H, J, C; dt = dt)
    ps = purity.(sol.ρ)
    pass = prod(abs.(1 .- ps) .< tolerance)
    pass_string = pass ? "passed" : "failed"

    make_test_plot(makeplots, sol, solve, plotpath, testset_id, string(pass_string, "_", test_id))

    # return plot(blochtimeseries, ts, exps...)
    return pass
end

function test_positivity_purity(; tolerance=1e-6, kwargs...)

    @testset "positive trajectory" begin
        @test positive_trajectory(; tolerance=tolerance, kwargs...)
    end
    @testset "maintains purity" begin
        @test maintains_purity(; tolerance=tolerance, kwargs...)
    end
end



