module QubitPlots

using Suppressor: @suppress_err
@suppress_err using Plots

using LaTeXStrings
using Statistics
using QuantumCircuits

export bloch_plots, qubit_plot
export single_qubit_plots, bell_plot
export blochsphere, blochtimeseries, blochprojections

include("../operators/single_qubit_operators.jl")
include("../operators/two_qubit_operators.jl")
# include("../operators/three_qubit_operators")

include("single_qubit_plots.jl")
include("two_qubit_plots.jl")
include("recipes.jl")

end # module QubitPlots
