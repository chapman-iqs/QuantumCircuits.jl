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

export runtests, test_timedelay, test_lindblad, test_superops, test_sbayesian, test_ssbayesian


include("timedelay.jl")
include("lindblad.jl")
include("superoperators.jl")
include("run.jl")
# include("unit.jl")
# include("analytical.jl")


function runtests(; functions = [test_superops, test_lindblad, test_timedelay, test_sbayesian], kwargs...)
    for f in functions
        f(; kwargs...)
    end
end


end
