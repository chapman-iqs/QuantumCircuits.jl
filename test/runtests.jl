include("../src/QuantumCircuits.jl")

using .QuantumCircuits
using .QuantumCircuits.Tests


# runtests(; functions = [test_ssbayesian, test_sbayesian], verbose=true, n = 10)
runtests(;  
            functions=[ 
                        test_integration, 
                        test_single_timestep, 
                        test_lindblad, 
                        test_timedelay,
                        test_readout,
                        test_bayesian,
                        test_positivity_purity,
                        test_ensemble,
                        ],
            solvers = [
                        bayesian,
                        rouchon
                        ],
            verbose=true, 
            n = 3, 
            makeplots=true,
            plotpath="test_result_plots")
