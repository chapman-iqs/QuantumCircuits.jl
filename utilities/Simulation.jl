module Simulation

export Simulation

using Reexport
using QuantumCircuits
using CSV
@reexport using Parameters
@reexport using DataFrames
import Dates

"""
	mutable struct Simulation(pars, data, path, description)

pars                    -- Parameter struct containing simulation parameters values as fields
data::Dict{DataFrame}   -- Dictionary of DataFrames containing data, where keys are filenames
path::String            -- path to data
description::String     -- optional description of the simulation

Wrapper for simulation data and parameters to create consistency in importing and exporting data.

"""
@with_kw mutable struct Simulation
    pars
    data::Dict{DataFrame} = Dict()
    path::String = "/Users/sachagreenfield/Desktop/Physics/Research/data/unnamed_data"
    description::String = ""
end
Simulation(pars) = Simulation(pars = pars)

function exportpath(path)
	ep = string(path, "/", Dates.now(), "/")
	replace(ep, ":" => "-")
end

function write(sim::Simulation)

	ep = exportpath(sim.path)
	mkpath(ep)
	println(string("Writing to path ," ep))

	writepars(sim.pars)
	writedescription(sim.description)
	writedata(sim.data)

end

function writepars(pars::T) where {T}

	d = Dict()

	for n in fieldnames(T)
		d[string(n)] = getfield(pars, n)
	end

	CSV.write(string(ep, "parameters.csv"), DataFrame(d))
end

function writedescription(des::String)
	CSV.write(string(ep, "description.txt"), des)
end

function writedata(data::Dict{DataFrame})
	for (k, v) in data
		filename = string(ep, k, ".csv")
		CSV.write(string(ep, filename), v)
	end
end





# function writedata(pars::OrderedDict, dataframes::Vector{DataFrame}, datanames::Vector{String},
# 					path_to_data = datapath, foldername = defaultfolder)
#
# 	ep = exportpath(path_to_data, foldername)
# 	mkpath(ep)
# 	println(string("Writing to path ", ep))
#
# 	CSV.write(string(ep, "parameters.csv"), pars)
#
# 	for (data, name) in zip(dataframes, datanames)
# 		CSV.write(string(ep, name, ".csv"), data)
# 	end
# end
#
# function writedata(pars::OrderedDict, dataframe::DataFrame, dataname::String,
# 					path_to_data = datapath, foldername = defaultfolder)
#
# 	ep = exportpath(path_to_data, foldername)
# 	mkpath(ep)
# 	println(string("Writing to path ", ep))
#
# 	CSV.write(string(ep, "parameters.csv"), pars)
# 	CSV.write(string(ep, dataname, ".csv"), dataframe)
# end





end # module Simulation
