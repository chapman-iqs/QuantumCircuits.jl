"""
Eventually, want to fill this out with tests including

-   unit tests: for individual functions (e.g. meas) check that bayesian
    reproduces the "calculated-by-hand" answer

-   analytical checks: systems for which a master equation solution exists.
    Check that bayesian and rouchon reproduce analytical solution within
    some error

-   ensemble average: check that stochastic trajectories properly
    average to the η = 0 (deterministic) solution

-   cross-checks: check that bayesian and rouchon reproduce each other's
    results to within some acceptable error by feeding them the same
    record. Save the record as data and check expected results for each
    method across versions of QuantumCircuits.jl.

-   timing tests: benchmark performance and track across versions


"""
module Tests

using ..QuantumCircuits
using ..QuantumCircuits.SingleQubitOperators
using Test
using Plots, Measures, Random, Distributions
import LinearAlgebra: eigvals
include("../plots/single_qubit_plots.jl")

export runtests

include("timedelay.jl")
export test_timedelay

include("lindblad.jl")
export test_lindblad, lindblad_const, lindblad_timedep

# include("superoperators.jl")
# export test_superops, test_sbayesian, test_ssbayesian

include("positivity.jl")
export test_positive_trajectory, positive_trajectory

include("single-time-step.jl")
export test_single_timestep, bayesian_update, lindblad_Γ2_decay, ham_update

include("integration-tests.jl")
export test_integration, returns_solution

include("record.jl")
export test_readout, record_distribution

include("bayesian.jl")
export test_bayesian, bayesian_dt_robustness

# include("run.jl")
# include("unit.jl")
# include("analytical.jl")



function runtests(; functions = [
                                test_integration,
                                 test_positive_trajectory,
                                 test_single_timestep,
                                 test_lindblad, 
                                 test_timedelay,
                                 test_readout
                                 ], 
                    solvers = [bayesian, rouchon],
                    makeplots = false,
                    kwargs...)

    for solve in solvers
        println("\n\n**************** Running tests for ", solve, " solver: **************** \n")
        for f in functions
            f(; makeplots=makeplots, solve=solve, kwargs...)
        end
    end
end


function make_test_plot(makeplots::Bool, sol, solve, plotpath::String, testset_id::String, test_id::String)
    path = mkpath(joinpath(plotpath, string(solve), testset_id))
    sol = Solution(sol, qbasis)
    plot(blochtimeseries, sol; title=string(solve))
    savefig(joinpath(path, test_id))
end


end
