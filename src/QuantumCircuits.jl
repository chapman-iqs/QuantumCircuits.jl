module QuantumCircuits

using Reexport
@reexport using QuantumOpticsBase
using Distributions, Distributed, ProgressMeter

include("types.jl")
export QOp, State, States, Solution, Ensemble
export Timescale, Times, Rate, Efficiency, Record, Records, Readout, Readouts

include("utils.jl")
export fidelity, expectations, purity, purities, average_purity, ensemble_average, coarse_grain, subselect, δ

include("readout.jl")
include("bayesian.jl")
include("rouchon.jl")
include("ensemble.jl")
export bayesian, rouchon, ensemble

include("../utilities/SingleQubitOperators.jl")
include("../utilities/TwoQubitOperators.jl")


end # module
