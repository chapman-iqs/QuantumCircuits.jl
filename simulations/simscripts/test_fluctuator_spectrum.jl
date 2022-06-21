const description =
"""
Written: June 20, 2022

Description: Generate plots of fluctuator spectrum (power spectral density) for a specific method of
fluctuator simulation. Goal is to verify 1/f noise spectrum. Specifically, trying to optimize Nf (number of
fluctuators) to use the minimum number of fluctuators required to get a 1/f noise spectrum, since the
simulations are time intensive.

/Users/sachagreenfield/Desktop/GitHub/QuantumCircuits.jl/simulations/simscripts/test_fluctuator_spectrum.jl

"""

using Pkg
Pkg.activate(".")

using QuantumCircuits
using Suppressor: @suppress_err

include(string(pwd(), "/simulations/fluctuators/FluctuatorSims.jl"))
include(string(pwd(), "/utilities_new/helpers/helpers.jl"))
include(string(pwd(), "/utilities_new/operators/single_qubit_operators.jl"))
include(string(pwd(), "/utilities_new/plotting/QubitPlots.jl"))

using .FluctuatorSims
using .QubitPlots
using OrderedCollections
using CSV, DataFrames
@suppress_err using Plots 	# suppressing precompilation errors of the form, "The call to compile
							# cache failed to create a usable precompiled cache file"

const path = "/Users/sachagreenfield/Desktop/Physics/Research/2021-Excitation-feedback/data/"
const folder = "fluctuator_sims/test_fluctuator_spectrum"
ep = exportpath(path, folder)
mkpath(ep)
println(string("Writing to path ", ep))
open(string(ep, "description.txt"), "w") do f
    write(f, description)
end

pars = OrderedDict("N" => 100, "tf" => 100.0, "dt" => 20e-3)
CSV.write(string(ep, "parameters.csv"), pars)

for N in [10, 20, 50, 100, 500, 1000]
	p = pars
	p["N"] = N
	(ts, χtotal_traj) = fluctuator_test(p)
	(freqs, psd) = PSD(ts, χtotal_traj)

	df = DataFrame([collect(ts), χtotal_traj], ["ts", "χtotal"])
	CSV.write(string(ep, "χtotal_N", N, ".csv"), df)

	pf = plot_fluctuators(ts, χtotal_traj);
	savefig(string(ep, "/fluctuator_timeseries_N", N, ".png"))

	ps = plotPSD(freqs, psd; scale=20e3, logscale=true);
	savefig(string(ep, "/fluctuator_PSD_N", N, ".png"))
end
