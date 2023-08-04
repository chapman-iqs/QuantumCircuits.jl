function ensemble_avg(Ω, Γ, η; N=300, solve=bayesian, ψ0 = randstate(SpinBasis(1//2)), tolerance=1/sqrt(N), plotpath="test_result_plots", test_id="ensemble_avg", testset_id="", makeplots=true, kwargs...)
    ψ0 = dm(normalize(g + e)) # initialize in the +x state
    tf = 2.5
    dt = 1e-3
    N = 300

    # set up Kraus operators
    H = Ω/2 * σy
    J = η < 1.0 ? [(σz, (1 - η) * Γ/2)] : [] # here, we use Γ/2 because of the relation between bloch coordinate dephasing and density matrix dephasing for a single qubit
    C(η) = [(σz, Γ, η)]

    # solve an ensemble of trajectories, and find the corresponding η = 0 solution
    ens = Ensemble(ensemble(solve, (0.0, tf), ψ0, H, J, C(η); dt=dt, N = N), qbasis, :qbasis);
    # sol0 = Solution(solve((0.0, tf), ψ0, H, [], C(0.0); dt=dt), qbasis, :qbasis); # this may not be the correct reference for the current rouchon code
    sol0 = Solution(solve((0.0, tf), ψ0, H, [(σz, Γ/2)], []; dt=dt), qbasis, :qbasis); # this may not be the correct reference for the current rouchon code
    solavg = Solution(ens.t, ensemble_average(ens.sols))

    atd = average_trace_distance(sol0.ρ, solavg.ρ)
    pass = atd < tolerance
    pass_string = pass ? "passed" : "failed"

    # plot the ensemble of trajectories and  η = 0 trajectory along with the master equation solution
    if makeplots
        path = mkpath(joinpath(plotpath, string(solve), testset_id, test_id))
        plot(blochensemble, ens, sol0; legend=:outerright, size=(1000,500), bottommargin=5mm, title="d = $atd")
        savefig(joinpath(path, "$(pass_string)_Ω_$(round(Ω/2π, digits=3))_MHz_Γ_$(Γ)_η_$(η).png"))
    end

    return pass

end

average_trace_distance(sol1::Solution, sol2::Solution) = average_trace_distance(sol1.ρ, sol2.ρ)
function average_trace_distance(ρ1s, ρ2s)
    trace_distances = map(zip(dm.(ρ1s), dm.(ρ2s))) do (ρ1, ρ2)
                        tracedistance(ρ1, ρ2)
                      end
    return mean(trace_distances)
end

function test_ensemble(; solve=bayesian, testset_id="ensemble", kwargs...)

     @testset "ensemble average" begin
        # @test ensemble_avg(2π, 0.5, 1.0; solve=solve, testset_id=testset_id, kwargs...)
        @test ensemble_avg(2π, 0.5, 0.9; solve=solve, testset_id=testset_id, kwargs...)
        @test ensemble_avg(2π, 0.5, 0.5; solve=solve, testset_id=testset_id, kwargs...)
        @test ensemble_avg(2π, 0.5, 0.1; solve=solve, testset_id=testset_id, kwargs...)
    end
end

