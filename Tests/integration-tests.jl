function returns_solution(; solve=bayesian, makeplots=false, testset_id="", test_id="returns_solution", plotpath="test_result_plots", kwargs...)

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
    pass = isa(sol, Solution)
    pass_string = pass ? "passed" : "failed"
    
    make_test_plot(makeplots, sol, solve, plotpath, testset_id, string(pass_string, "_", test_id))

    # return plot(blochtimeseries, ts, exps...)
    return pass
end


function test_integration(; solve=bayesian, testset_id="integration", kwargs...)

    @testset "returns solution ($(solve))" begin
        Random.seed!()
       @test returns_solution(; solve=solve, testset_id=testset_id, kwargs...)   
    end
end