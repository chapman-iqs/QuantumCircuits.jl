include("../src/QuantumCircuits.jl")

using .QuantumCircuits
using .QuantumCircuits.Tests


# runtests(; functions = [test_ssbayesian, test_sbayesian], verbose=true, n = 10)
runtests(;  verbose=true, 
            n = 3, 
            makeplots=true,
            plotpath="test_result_plots")
