using Distributed

addprocs(2)


@everywhere begin

	N = 100
	exportpath = "data/run4/"
	mkpath(exportpath)

	cd("/Users/sachagreenfield/Desktop/GitHub/QuantumCircuits.jl")

	using Pkg
	Pkg.activate(".")
	Pkg.instantiate()


	# Pkg.add(QuantumCircuits)

	using QuantumCircuits

	using Random
	using Statistics
	using Distributions
	using DataFrames
	using CSV


	include("../utilities.jl")
	ideal = false


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
	η =  0.6		# efficiency
	φ = 0			# angle
	td = 200 * 1e-3 # time delay for feedback


	"simulation timescales"

	T = (0, 8τm) # simulation duration
	dt = 0.5e-3  # integration time-step


	"feedback target parameters"

	θs = 3π/10 		# target angle on Bloch sphere
	Rs = 1 	# 0.64 			# radius of target
	ϕ = π 					# fixes plane of oscillations
	vec = (sin(θs) * cos(π/2), sin(θs) * sin(π/2), cos(θs))
	σϕ = cos(ϕ)*σx + sin(ϕ)*σy

	"feedback drive parameters"

	Δ0 = ideal ? -sin(2θs)/(4τm) : -sin(2θs)/(4τm*Rs^2)
	Δ1 = ideal ? sin(θs)/τm : sin(θs)/(τm*Rs)


	"decay times and rates"

	T1 = 40 	# energy decay time
	Γ1 = 1/(2T1)# energy decay rate
	T2 = 60 	# environmental dephasing time
	Γ2 = 1/T2 	# environemntal dephasing rate


	"Kraus operators -------------------------------------------------------------"

	H = Hf(t, r) = (Δ0 + Δ1*r[1]) * σϕ/2
	J = ideal ? [(σz, ((1-η)*Γm))] : [(σz, ((1-η)*Γm)), (σm, Γ1), (σz, Γ2)]
	C = [(exp(im * φ) * σz, τm, η)]

	"Bayesian simulation ---------------------------------------------------------"
	#
	# trajs = pmap(m -> begin
	# 		sol = bayesian(T, ρ0, H, J, C; dt=dt, td=td)
	# 		return traj(sol)
	# 	end, 1:N; batch_size=10)

	@time trajs = pmap(1:N; batch_size=10) do m
				println(string("Running trajectory ", m, "..."))
				sol = bayesian(T, ρ0, H, J, C; dt=dt, td=td)
				return traj(sol)
			end

	# loop over N simulations
	# trajs = []
	#
	# for i in 1:N
	# 	sol = bayesian(T, ρ0, H, J, C; dt=dt, td=td)
	# 	tr = traj(sol)
	# 	push!(trajs, tr)
	# end

	t = collect(T[1]:dt:T[2])   # list of times

	dfx = DataFrame([tr.x for tr in trajs], [string("x", i) for i in 1:N])
	dfy = DataFrame([tr.y for tr in trajs], [string("y", i) for i in 1:N])
	dfz = DataFrame([tr.z for tr in trajs], [string("z", i) for i in 1:N])
	dfr = DataFrame([tr.r[1] for tr in trajs], [string("r", i) for i in 1:N])
	dft = DataFrame([t], ["t"])


	CSV.write(string(exportpath, "x.csv"), dfx)
	CSV.write(string(exportpath, "y.csv"), dfy)
	CSV.write(string(exportpath, "z.csv"), dfz)
	CSV.write(string(exportpath, "r.csv"), dfr)
	CSV.write(string(exportpath, "t.csv"), dft)


	"Parameter dictionary"

	pars = Dict("N" => N,
				"x0" => x0,
				"y0" => y0,
				"z0" => z0,
				"τm" => τm,
				"Γm" => Γm,
				"x0" => x0,
				"η" => η,
				"φ" => φ,
				"td" => td,
				"dt" => dt,
				"θs" => θs,
				"Rs" => Rs,
				"ϕ" => ϕ,
				"Δ0" => Δ0,
				"Δ1" => Δ1,
				"Τ1" => T1,
				"Τ2" => T2,
				"ideal" => ideal
				)

	CSV.write(string(exportpath, "parameters.csv"), pars)
end
