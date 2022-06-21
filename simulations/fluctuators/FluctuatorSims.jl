module FluctuatorSims

export makefluctuators
export fluctuator_test, fluctuator_test_defaultpars, T2simulation, T2simulation_defaultpars
export plot_fluctuators, PSD, plotPSD
export writedata, exportpath

using Random, Distributions
using LsqFit
using ProgressMeter
using Dates, DataFrames, CSV
using OrderedCollections
using QuantumCircuits

# change this to wherever you want to save data
const datapath = "/Users/sachagreenfield/Desktop/Physics/Research/2021-Excitation-feedback/data/"
const defaultfolder = "fluctuator_sims"

include("../../utilities_new/operators/single_qubit_operators.jl")

include("fluctuators.jl")
include("simulations.jl")
include("analysis.jl")
include("writedata.jl")


end # module FluctuatorSims
