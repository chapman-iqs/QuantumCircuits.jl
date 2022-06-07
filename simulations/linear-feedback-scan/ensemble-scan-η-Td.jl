
N = 1000

cd("/Users/sachagreenfield/Desktop/GitHub/QuantumCircuits.jl")


using QuantumCircuits

using Random
using Statistics
using Distributions
using DataFrames
using CSV
using Dates
using ProgressMeter

exportpath = let
				ep = string("data/", Dates.now(), "/")
				replace(ep, ":" => "-") end
mkpath(exportpath)


include("../utilities.jl")
ideal = false
dualquad = false


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

ground = spindown(q)
excited = spinup(q)

"""
System parameters -------------------------------------------------------------------------------------------------
all times given in μs
"""

# initial state
(x0,y0,z0) = (0., 0.3, 0.91)
ρ0 = DenseOperator(0.5*(id + x0*σx + y0*σy + z0*σz))


"measurement parameters"

τm = 2			# time
Γm = 1/(2τm) 	# rate
# ηs = [0.3, 0.4]
ηs = range(0.05, 1, step=0.05)		# efficiency
φ = 0			# angle
# tds = range(0, 500, step=250) .* 1e-3  # time delay for feedback
tds = range(0, 500, step=50) .* 1e-3


"decay times and rates"

T1 = 40 	# energy decay time
Γ1 = 1/(2T1)# energy decay rate
T2 = 60 	# environmental dephasing time
Γ2 = 1/T2 	# environemntal dephasing rate


"simulation timescales"

T = (0, 8τm) # simulation duration
dt = 1e-3  # integration time-step

#

"Loop over η, td -------------------------------------------------------------"
ηl, tdl, fl, vl = [], [], [], []
v(ser) = var(ser[5 * Int64(floor(length(ser)/6)):end]) # calculates the variance of the last 1/6 of a series,
														# to test for convergence

ϕ = π 		# fixes plane of oscillations
σϕ = cos(ϕ)*σx + sin(ϕ)*σy


for η in ηs

	for td in tds

		println(string("Running η = ", η, ", td = ", td))

		"feedback target parameters"

		θs = 3π/10 		# target angle on Bloch sphere
		Rs = 1/sqrt(1/η + 2τm / T2)	# radius of target

		vec = (sin(θs) * cos(π/2), sin(θs) * sin(π/2), cos(θs))

		"feedback drive parameters"

		Δ0 = ideal ? -sin(2θs)/(4τm) : -sin(2θs)/(4τm*Rs^2)
		Δ1 = ideal ? sin(θs)/τm : sin(θs)/(τm*Rs)
		ΔS = 0.


		"Kraus operators -------------------------------------------------------------"

		H(t, r) = dualquad ?
					(Δ0 + Δ1*r[1]) * σϕ/2 + ΔS * r[2] * σz/2 :
					(Δ0 + Δ1*r[1]) * σϕ/2
		J = ideal ?
				[(σz, ((1-η)*Γm))] :
				[(σz, ((1-η)*Γm + Γ2)), (σm, Γ1)]

		C = dualquad ?
			[(exp(im * φ) * σz, τm, η/2), (exp(im * (φ + π/2)) * σz, τm, η/2)] :
			[(exp(im * φ) * σz, τm, η)]

		"Bayesian simulation ---------------------------------------------------------"

		# x, y, z = [], [], []
		# @showprogress for i in 1:N
		# 	sim = traj(bayesian(T, ρ0, H, J, C; dt=dt, td=td))
		# 	push!(x, last(sim.x))
		# 	push!(y, last(sim.y))
		# 	push!(z, last(sim.z))
		# end
		#
		# x, y, z = mean.([x, y, z])
		# r = sqrt(x^2 + y^2 + z^2)
		# f = (1 + r)/2

		@time trajs = @showprogress map(1:N) do m
					return sim = traj(bayesian(T, ρ0, H, J, C; dt=dt, td=td))
				end

		xm = mean([tr.x for tr in trajs])
		ym = mean([tr.y for tr in trajs])
		zm = mean([tr.z for tr in trajs])
		vm = v(xm) + v(ym) + v(zm)

		x, y, z = [last(a) for a in [xm, ym, zm]]
		r = sqrt(x^2 + y^2 + z^2)
		f = (1 + r)/2

		push!(ηl, η)
		push!(tdl, td)
		push!(fl, f)
		push!(vl, vm)

		df = DataFrame([ηl, tdl, fl, vl], ["η", "td", "F", "var"])
		CSV.write(string(exportpath, "scan-η-td.csv"), df)

		println()

	end

end



"Write data --------------------------------------------------------------------"

df = DataFrame([ηl, tdl, fl, vl], ["η", "td", "F", "var"])
# df = DataFrame([ηl, tdl, fl], ["η", "td", "F"])
CSV.write(string(exportpath, "scan-η-td.csv"), df)



"Parameter dictionary"

pars = Dict("N" => N,
			"x0" => x0,
			"y0" => y0,
			"z0" => z0,
			"τm" => τm,
			"Γm" => Γm,
			"x0" => x0,
			"η" => ηs,
			"φ" => φ,
			"td" => tds,
			"dt" => dt,
			"ϕ" => ϕ,
			"T1" => T1,
			"T2" => T2,
			"ideal" => ideal,
			"dualquad" => dualquad
			)

CSV.write(string(exportpath, "parameters.csv"), pars)
