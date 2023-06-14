module QuantumCircuits

using Reexport
@reexport using QuantumOpticsBase
import QuantumOpticsBase: dm
using Distributions, Distributed, ProgressMeter

include("types.jl")
export QOp, State, States, Solution, Ensemble
export Timescale, Times, Rate, Efficiency, Record, Records, Readout, Readouts

include("utils.jl")
export fidelity, expectations, purity, purities, average_purity, ensemble_average, coarse_grain, subselect, δ

include("readout.jl")
include("gausskraus.jl")
include("measlindham.jl")

include("feedback.jl")
export Feedback, ForwardEstimation, NoFeedback, HamiltonianFE, HamLindFE, nofeedback, hamiltonianfe, hamlindfe
include("bayesian.jl")
include("rouchon.jl")
include("ensemble.jl")
export bayesian, rouchon, ensemble

include("superevolution.jl")
export sbayesian, ssbayesian, forwardtrajectory, forwardbayesian


include("../utilities/SingleQubitOperators.jl")
include("../utilities/TwoQubitOperators.jl")
include("../Tests/Tests.jl")


end # module
