include("../src/QuantumCircuits.jl")

using .QuantumCircuits
using .QuantumCircuits.Tests

runtests(; functions = [test_ssbayesian, test_sbayesian], verbose=true, n = 10)