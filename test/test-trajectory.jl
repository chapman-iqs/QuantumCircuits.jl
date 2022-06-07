N = 100						# number of trajectories
writeblochs = true		# export bloch coordinates to CSV?
loadrecords = false
loadpath = string("test/purestate_test_data/", "2022-01-15T14-52-32.301-dm", "/")

cd("/Users/sachagreenfield/Desktop/GitHub/QC-notebooks")

using QuantumCircuits

using Random
using Statistics
using Distributions
using DataFrames
using CSV
using Dates
using ProgressMeter

exportpath = let
				ep = string("test/purestate_test_data/", Dates.now(), "/")
				replace(ep, ":" => "-") end
mkpath(exportpath)

if loadrecords
	R1 = DataFrame(CSV.File(string(loadpath, "r1.csv")))
	N = size(R1)[2] # number of records, corresponding to number of simulations
	records = [R1[!,i] for i in 1:N]
	println("records loaded")
end


include("../utilities.jl")

"Qubit Hilbert space operators -------------------------------------------------------------------------------------"

# Basis
q = SpinBasis(1//2)

# Operators, using convention that |-z> is ground state
σx = sigmax(q)
σy = sigmay(q)
σz = sigmaz(q)
σp = sigmap(q)
σm = sigmam(q)
id = identityoperator(q)

g = spindown(q)
e = spinup(q)

"""
System parameters -------------------------------------------------------------------------------------------------
all times given in μs
"""

# initial state
# (x0,y0,z0) = (0.0, 0.0, 1.0)
# ρ0 = DenseOperator(0.5 * (id + x0*σx + y0*σy + z0*σz))
ρ0 = normalize(g + e)


"measurement parameters"

Γ = 1.0						# measurement rate
τ = 1/(2Γ) 				# measurement time
η = 1.0						# quantum efficiency
φ = 0.0						# measurement angle
# td = 200 * 1e-3 	# time delay for feedback

"decay parameters"
T1 = 40 	# energy decay time
Γ1 = 1/(2T1)# energy decay rate
T2 = 60 	# environmental dephasing time
Γ2 = 1/T2 	# environemntal dephasing rate


"simulation timescales"

T = (0.0, 10.0) 	# simulation duration
dt = 1e-3  		# integration time-step



"Kraus operators -------------------------------------------------------------"

ΩR = π/2  # Rabi rate in rad MHz

H = ΩR * σy
J = []
# J = [(σz, ((1 - η) * Γ)), (σm, Γ1), (σz, Γ2)]
C = [(σz, Γ, η)]

"Bayesian simulation ---------------------------------------------------------"

@time trajs = @showprogress map(1:N) do m
					sol = loadrecords ? bayesian(T, ρ0, H, J, C; records=[records[m]], dt=dt) :
										bayesian(T, ρ0, H, J, C; dt=dt)

					return traj(sol)
			end


"Write data --------------------------------------------------------------------"

if writeblochs

	t = collect(T[1]:dt:T[2])   # list of times

	dfx = DataFrame([tr.x for tr in trajs], [string("x", i) for i in 1:N])
	dfy = DataFrame([tr.y for tr in trajs], [string("y", i) for i in 1:N])
	dfz = DataFrame([tr.z for tr in trajs], [string("z", i) for i in 1:N])
	dft = DataFrame([t], ["t"])

	CSV.write(string(exportpath, "x.csv"), dfx)
	CSV.write(string(exportpath, "y.csv"), dfy)
	CSV.write(string(exportpath, "z.csv"), dfz)
	CSV.write(string(exportpath, "t.csv"), dft)

	dfr1 = DataFrame([tr.r[1] for tr in trajs], [string("r1", i) for i in 1:N])
	CSV.write(string(exportpath, "r1.csv"), dfr1)

	# if dualquad
	# 	dfr2 = DataFrame([tr.r[2] for tr in trajs], [string("r2", i) for i in 1:N])
	# 	CSV.write(string(exportpath, "r2.csv"), dfr2)
	# end
end


"Parameter dictionary  --------------------------------------------------------------------"

pars = Dict("N" => N,
			# "x0" => x0,
			# "y0" => y0,
			# "z0" => z0,
			"Γ" => Γ,
			"τ" => τ,
			"η" => η,
			"φ" => φ,
			"dt" => dt,
			"ΩR" => ΩR
			# "T1" => T1,
			# "T2" => T2
			)

CSV.write(string(exportpath, "parameters.csv"), pars)

println(string("Data output to ", exportpath))
