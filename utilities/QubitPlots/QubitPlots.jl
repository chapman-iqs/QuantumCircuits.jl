module QubitPlots

"---- Export symbols ----"

export blochwireframe, blochsphere
export BlochEnsemble, blochensemble
export SeriesHistogram, serieshistogram
export BlochTimeSeries, blochtimeseries
export qubit_plot
export blochprojections
export bell_plot


"---- Dependencies ----"

using Suppressor: @suppress_err
@suppress_err using Plots
import Makie, GLMakie, CairoMakie

using LaTeXStrings
using DataFrames
using Statistics
using Measures
using QuantumCircuits

include("../SingleQubitOperators.jl")
include("../TwoQubitOperators.jl")
include("../NotebookHelpers/NotebookHelpers.jl")

using .SingleQubitOperators
using .TwoQubitOperators
using .NotebookHelpers


"---- Definitions of functions and constants ----"

include("blochsphere.jl")
include("ensembleplots.jl")
include("single_qubit_plots.jl")
include("two_qubit_plots.jl")
include("heatplots.jl")

end # module QubitPlots
