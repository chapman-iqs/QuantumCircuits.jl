path_to_repo = Sys.isapple() ?
				"/Users/sachagreenfield/Desktop/GitHub/QC-notebooks/" :
				"C:\\Users\\lfl\\Documents\\Sacha\\Excitation-feedback\\QC-notebooks"
cd(path_to_repo)

using Pkg
Pkg.activate(".")
using Statistics
using LsqFit
using CSV
using DataFrames
using Dates
using QuantumCircuits
using LaTeXStrings
using ProgressMeter

using Distributed
while true
	print(string("There are ", Sys.CPU_THREADS, " logical cores available. \nHow many processes do you want to add? "))
	procs = readline()
	println()
	try
		procs = parse(Int64, procs)
		addprocs(procs; exeflags=string("--project=", path_to_repo))
		break
	catch e
    	@warn string("Please enter an integer between 1 and ", Sys.CPU_THREADS - 1)
	end
end


# include Hilbert space operators
if Sys.isapple()
	include(string(path_to_repo, "/utilities/single-qubit-operators.jl"))
	include(string(path_to_repo, "/utilities/two-qubit-operators.jl"))
else
	include(string(path_to_repo, "\\utilities\\single-qubit-operators.jl"))
	include(string(path_to_repo, "\\utilities\\two-qubit-operators.jl"))
end

# Set export path and activate packages
@everywhere begin

	path_to_data = Sys.isapple() ?
					"/Users/sachagreenfield/Desktop/Physics/Research/2021-Excitation-feedback/data/" :
					"\\Users\\lfl\\Documents\\Sacha\\Excitation-feedback\\data\\"
	path_to_repo = Sys.isapple() ?
					"/Users/sachagreenfield/Desktop/GitHub/QC-notebooks/" :
					"C:\\Users\\lfl\\Documents\\Sacha\\Excitation-feedback\\QC-notebooks"

	function exportpath(foldername)
		if Sys.isapple()
			ep = string(path_to_data, foldername, "/", Dates.now(), "/")
			replace(ep, ":" => "-")
		else
			ep = string(path_to_data, foldername, "\\", Dates.now(), "\\")
			string("C:", replace(ep, ":" => "-"))
		end
	end

	using Pkg
	Pkg.activate(".")
	using QuantumCircuits

end


"""
	singlequbit_scan(; write_exp_values=true, N=1000, batch_size=10, Γ=1.0, η=1.0, ΩR=π/2, dt=1e-3, tf=10.0, comment=""

Keyword Arguments:
N :: Int64, 				-- number of trajectories
write_exp_values :: Bool	-- whether or not to output expectation values as csv
batch_size :: Int64			-- number of trajectories allocated to a worker at a time
Γ :: Rate					-- measurement rate,
η :: Efficiency				-- quantum efficiency,
ΩR :: Rate					-- Rabi rate,
dt :: Timescale 			-- integration time step
tf :: Timescale				-- simulation duration
σs :: StepRangeLen,			-- σ (noise std. dev.) values to scan over
comment :: String; 			-- comment for parameter file for future reference
"""
function singlequbit_scan(; write_exp_values=true, N=1000, batch_size=10, Γ=1.0, η=1.0, ΩR=π/2, dt=1e-3, tf=10.0, comment="")
		"Export parameters -------------------------------------------------------------------------------------"

		ep = exportpath("single-qubit")
		mkpath(ep)
		println(string("Writing to path ", ep))

		pars = Dict("N" => N,
					"batch_size" => batch_size,
					"Γ" => Γ,
					"η" => η,
					"ΩR" => ΩR,
					"dt" => dt,
					"tf" => tf,
					"comment" => comment
					)

		CSV.write(string(ep, "parameters.csv"), pars)

		# initial state
		ρ0 = g


		"Kraus operators -------------------------------------------------------------------------------------------------"
		H = ΩR * σy
		J = η == 1.0 ? [] : [(n, ((1 - η) * Γ))]
		C = [(σz, Γ, η)]

		"Run parallelized code -------------------------------------------------------------------------------------------------"
		ops = [σx, σy, σz]
		op_labels = ["x", "y", "z"]
		exps = ensemble(bayesian, (0.0, tf), ρ0, H, J, C; dt=dt, N=N, batch_size=batch_size, ops=ops)

		"Extract and write expectation values -----------------------------------------------------------------------------------"
		if write_exp_values
			# write data
			for (arr, label) in zip(exps, op_labels)
				df = DataFrame(arr, :auto)
				CSV.write(string(ep, label, "-trajectory.csv"), df)
			end

		end
end # function singlequbit_scan

"""
	singlequbit_T2_scan(; N=1000, batch_size=10, dt=1e-3, tf=20.0, σs=0.5:0.5:10.0, comment="")

Keyword Arguments:
write_exp_values :: Bool	-- whether or not to output expectation values as csv
N :: Int64, 			-- number of trajectories
batch_size :: Int64		-- number of trajectories allocated to a worker at a time
dt :: Timescale 		-- integration time step
tf :: Timescale			-- simulation duration
σs :: StepRangeLen,		-- σ (noise std. dev.) values to scan over
comment :: String; 		-- comment for parameter file for future reference
"""
function singlequbit_T2_scan(; write_exp_values=false, N=1000, batch_size=10, dt=1e-3, tf=20.0, σs=0.5:0.5:10.0, comment="")
		"Export parameters -------------------------------------------------------------------------------------"

		ep = exportpath("single-qubit-T2-scan")
		mkpath(ep)
		println(string("Writing to path ", ep))

		pars = Dict("N" => N,
					"batch_size" => batch_size,
					"dt" => dt,
					"tf" => tf,
					"comment" => comment
					)

		CSV.write(string(ep, "parameters.csv"), pars)

		# initial state
		ρ0 = normalize(e + g)


		"Kraus operators -------------------------------------------------------------------------------------------------"
		H(σ) = t -> σ * randn() * σz

		"Run parallelized code -------------------------------------------------------------------------------------------------"
		ops = [σx, σy]
		op_labels = ["x", "y"]

		# exponential model for fitting T2
		model(t, p) = p[1] * exp.(-p[2] * t)
		p0 = [0.5, 0.5]
		t = collect(0.0:dt:tf)

		T2s = map(σs) do σ
			println(string("Running σ = ", σ, " ..."))
			# get ensemble x-expectations
			(x, y) = ensemble(bayesian, (0.0, tf), ρ0, H(σ), [], []; dt=dt, N=N, batch_size=batch_size, ops=ops)
			xavg = [mean([x[i][j] for i in 1:N]) for j in 1:length(t)]
			yavg = [mean([y[i][j] for i in 1:N]) for j in 1:length(t)]

			# write expectation values
			if write_exp_values
				dfx = DataFrame(x, :auto)
				dfy = DataFrame(y, :auto)
				CSV.write(string(ep, "σ-", σ, "-x-trajectory.csv"), dfx)
				CSV.write(string(ep, "σ-", σ, "-y-trajectory.csv"), dfy)
			end

			# fit model
			fit = curve_fit(model, t, xavg, p0)
			param = fit.param
			T2x = 1/(2 * param[2])

			fit = curve_fit(model, t, yavg, p0)
			param = fit.param
			T2y = 1/(2 * param[2])

			return (T2x + T2y) / 2
		end # T2s = map(σs) do σ


		df = DataFrame([collect(σs), T2s], ["σs", "T2s"])
		CSV.write(string(ep, "T2-scan.csv"), df)
end # function singlequbit_T2_scan

"""
twoqubitfeedback_scan(; write_exp_values=true, N = 1000, batch_size = 10, Γ = 1.0, φ = π/2, σ = 1.0,
									η = 1.0, dt = 1e-3, tf = 10.0, td = 300e-3, Δt = 100e-3, T1 = 40, comment = "")s

Keyword Arguments:
N :: Int64, 				-- number of trajectories
batch_size :: Int64			-- number of trajectories allocated to a worker at a time
Γ :: Rate					-- measurement rate,
φ :: Float64				-- feedback drive angle
σ :: Float64				-- noise std. dev.
η :: Efficiency				-- quantum efficiency
dt :: Timescale 			-- integration time step
tf :: Timescale				-- simulation duration
td :: Timescale				-- time delay for feedback
Δt :: Timescale				-- Rabi pulse duration
T1 :: Timescale,			-- energy decay time
comment :: String; 			-- comment for parameter file for future reference
"""
function twoqubitfeedback_scan(; write_exp_values=true, N = 1000, batch_size = 10, Γ = 1.0, φ = π/2, σ = 1.5,
									η = 1.0, dt = 1e-3, tf = 10.0, td = 300e-3, Δt = 100e-3, T1 = 40, comment = "")
		"Export parameters -------------------------------------------------------------------------------------"

		ep = exportpath("two-qubit")
		mkpath(ep)
		println(string("Writing to path ", ep))

		pars = Dict("N" => N,
					"Γ" => Γ,
					"σ" => σ,
					"φ" => φ,
					"η" => η,
					"dt" => dt,
					"tf" => tf,
					"td" => td,
					"Δt" => Δt,
					"T1" => T1,
					"comment" => comment
					)

		CSV.write(string(ep, "parameters.csv"), pars)

		"System parameters -------------------------------------------------------------------------------------------------"
		# initial state
		ρ0 = g ⊗ g
		ρtarget = ket1p
		targetstr = "ket1p"


		"Kraus operators -------------------------------------------------------------------------------------------------"
		function θ(ρ::State)
			αβ = real(expect(ket1p ⊗ ketm', ρ))
			return 0.5 * asin(2αβ)
		end

		Hf(t::Timescale, ρ::State) = let
				 Ω =
					if t < Δt
						π / (4Δt)
					elseif real(expect(dm(ρtarget), ρ)) < 0.01 && t > td + Δt
						π / (4td) # to get length of pulse right
					elseif td == 0.0
						- θ(ρ) / (2dt) # try θ^3
					else
						- θ(ρ) / (8td)
					end

				return Ω * (cos(φ) * (σx1 + σx2) + sin(φ) * (σy1 + σy2))

			end

		Hz(t, σ) = let
				R1, R2 = σ * randn(), σ * randn()
				R1 * σz1 + R2 * σz2
		end

		Hs(σ) = (t::Timescale, ρ::State) -> Hf(t, ρ) + Hz(t, σ)

		J = η == 1.0 ? [(σm1, 1/T1), (σm2, 1/T1)] : [(n, ((1 - η) * Γ)), (σm1, 1/T1), (σm2, 1/T1)]
		C = [(n, Γ, η)]

		"Run parallelized code -------------------------------------------------------------------------------------------------"
		ops = dm.([ketp, ketm, ket1p, ket1m])
		op_labels = ["ketp", "ketm", "ket1p", "ket1m"]
		(system_exps, filter_exps) = ensemble(bayesian, (0.0, tf), ρ0, (Hs(σ), Hf), J, C; dt=dt, N=N, batch_size=batch_size, ops=ops)

		"Extract and write expectation values -----------------------------------------------------------------------------------"
		if write_exp_values
			# write data
			for (arr, label) in zip(system_exps, op_labels)
				df = DataFrame(arr, :auto)
				CSV.write(string(ep, "system-", label, "-trajectory.csv"), df)
			end

			for (arr, label) in zip(filter_exps, op_labels)
				df = DataFrame(arr, :auto)
				CSV.write(string(ep, "filter-", label, "-trajectory.csv"), df)
			end
		end
end # function twoqubitfeedback_scan

"""
twoqubitfeedback_parameter_sweep(; write_exp_values=true, N = 1000, batch_size = 10, Γs = 0.5:0.5:5.0, Ωmaxs = 0.1:1.0:10.0,
									 φ = π/2, σ = 1.5, ηs = 0.3:0.2:1.0, dt = 1e-3, tf = 10.0, td = 500e-3, Δt = 100e-3,
									 comment = "")

									 DOCUMENTATION NEEDS UPDATE

Keyword Arguments:
N :: Int64, 				-- number of trajectories
batch_size :: Int64			-- number of trajectories allocated to a worker at a time
Γ :: Rate					-- measurement rate,
φ :: Float64				-- feedback drive angle
σ :: Float64				-- noise std. dev.
η :: Efficiency				-- quantum efficiency
dt :: Timescale 			-- integration time step
tf :: Timescale				-- simulation duration
td :: Timescale				-- time delay for feedback
Δt :: Timescale				-- Rabi pulse duration
T1 :: Timescale,			-- energy decay time
comment :: String; 			-- comment for parameter file for future reference
"""
function twoqubitfeedback_parameter_sweep(; write_exp_values=true, N = 1000, batch_size = 10, Γs = 0.5:0.5:5.0, Ωmaxs = 0.1:1.0:10.0,
									 φ = π/2, σ = 1.5, ηs = 0.3:0.2:1.0, dt = 1e-3, tf = 10.0, td = 500e-3, Δt = 100e-3,
									 comment = "")
		"Export parameters -------------------------------------------------------------------------------------"

		ep = exportpath("two-qubit-sweep")
		mkpath(ep)
		println(string("Writing to path ", ep))

		pars = Dict("N" => N,
					"Γs" => Γs,
					"Ωmaxs" => Ωmaxs,
					"σ" => σ,
					"φ" => φ,
					"ηs" => ηs,
					"dt" => dt,
					"tf" => tf,
					"td" => td,
					"Δt" => Δt,
					"comment" => comment
					)

		CSV.write(string(ep, "parameters.csv"), pars)

		"System parameters -------------------------------------------------------------------------------------------------"
		# initial state
		ρ0 = g ⊗ g
		ρtarget = ket1p
		targetstr = "ket1p"


		"Kraus operators -------------------------------------------------------------------------------------------------"
		function θ(ρ::State)
			αβ = real(expect(ket1p ⊗ ketm', ρ))
			return 0.5 * asin(2αβ)
		end

		Ω(t, ρ, Ωmax) =
		   if t < Δt
			   π / (4Δt)
		   elseif real(expect(dm(ρtarget), ρ)) < 0.01 && t > td + Δt
			   π / (4td) # to get length of pulse right
		   elseif td == 0.0
			   - θ(ρ) / (2dt) # try θ^3
		   else
			   - θ(ρ) * Ωmax / (π/4)
		   end

		Hf(Ωmax) = (t::Timescale, ρ::State) -> Ω(t, ρ, Ωmax) * (cos(φ) * (σx1 + σx2) + sin(φ) * (σy1 + σy2))

		Hz(t) = let
				R1, R2 = σ * randn(), σ * randn()
				R1 * σz1 + R2 * σz2
		end

		Hs(Ωmax) = (t::Timescale, ρ::State) -> Hf(Ωmax)(t, ρ) + Hz(t)

		J(Γ, η) = (η == 1.0) ? [] : [(n, ((1 - η) * Γ))]
		C(Γ, η) = [(n, Γ, η)]

		"Loop over parameters -------------------------------------------------------------------------------------------------"
		ops = dm.([ketp, ketm, ket1p, ket1m])
		op_labels = ["ketp", "ketm", "ket1p", "ket1m"]

		sys_data = [Vector{Float64}[] for op in ops]
		fil_data = [Vector{Float64}[] for op in ops]

		for η in ηs
			for Γ in Γs
				for Ωmax in Ωmaxs
					println(string("Running η = ", η, ", Γ = ", Γ, ", Ωmax = ", Ωmax, " ..."))
					(system_exps, filter_exps) = ensemble(bayesian, (0.0, tf), ρ0, (Hs(Ωmax), Hf(Ωmax)), J(Γ, η), C(Γ, η); dt=dt, N=N, batch_size=batch_size, ops=ops)
					sys_mean = mean.(system_exps)
					fil_mean = mean.(filter_exps)

					# save data
					for i in 1:length(ops)# iterates over mean trajectories corresponding to `ops` expectation values
						push!(sys_data[i], vcat([η, Γ, Ωmax], sys_mean[i]))
						push!(fil_data[i], vcat([η, Γ, Ωmax], fil_mean[i]))
					end

				end # for Ωmax
			end # for Γσ
		end # for η


		"Write data  -----------------------------------------------------------------------------------"
		for (arr, label) in zip(sys_data, op_labels)
			df = DataFrame(arr, :auto)
			CSV.write(string(ep, "system-", label, "-mean-trajectory.csv"), df)
		end

		for (arr, label) in zip(fil_data, op_labels)
			df = DataFrame(arr, :auto)
			CSV.write(string(ep, "filter-", label, "-mean-trajectory.csv"), df)
		end
end # function twoqubitfeedback_parameter_sweep


function twoqubit_fluctuators(; write_exp_values=true, N = 1000, batch_size = 10, dt = 1e-3, tf = 10.0, Ω = (2π) * 0.2,
	 							ΩR = (2π) * 25.0, nfluctuators = 10, φd1 = π/2, φd2 = π/2, Γm = 0.0,
								measurementtype = "excitation number", foldername = "two-qubit-fluctuators", comment = "")
		"Export parameters -------------------------------------------------------------------------------------"

		ep = exportpath(foldername)
		mkpath(ep)
		println(string("Writing to path ", ep))

		pars = Dict("N" => N, # number of trajectories
					"dt" => dt, # simulation time step
					"tf" => tf, # simulation duration
					"Ω" => Ω, # fluctuator strength
					"ΩR" => ΩR, # decoupling strength
					"φd1" => φd1, # decoupling axis, qubit 1
					"φd2" => φd2, # decoupling axis, qubit 2
					"Γm" => Γm, # measurement rate
					"nfluctuators" => nfluctuators, # number of fluctuators
					"measurementtype" => measurementtype,
					"comment" => comment
					)

		CSV.write(string(ep, "parameters.csv"), pars)

		"System parameters -------------------------------------------------------------------------------------------------"
		# initial state
		ψ0 = Ψp

		"Setting up the fluctuators ---------------------------------------------------------------------------------------"
		# sample log-uniformly
		(ωmin, ωmax) = 1/tf, 1/dt
		ωlog = range(log(ωmin), log(ωmax), length=nfluctuators)
		ωi = exp.(ωlog)
		τi = map(ω -> 1/ω, ωi) # fluctuators and their signs for all $nfluctuators switching times

		"Hamiltonians -------------------------------------------------------------------------------------------------"
		# fluctuator Hamiltonians for each qubit
		Hf1(τi, tf) = let si = map(τ -> sgn(τ, tf), τi)
						t -> map(s -> (Ω/2) * σz1 * s(t), si) end

		Hf2(τi, tf) = let si = map(τ -> sgn(τ, tf), τi)
						t -> map(s -> (Ω/2) * σz2 * s(t), si) end

		# decoupling Hamiltonian; turn off decoupling by setting ΩR = 0
		H2q = (ΩR/2) * (cos(φd1) * σx1 + sin(φd1) * σy1 + cos(φd2) * σx2 + sin(φd2) * σy2)

		# total Hamiltonian -- a function of a particular realization of fluctuator Hamiltonian Hf
		H(Hf1, Hf2) = t -> H2q + sum(Hf1(t)) + sum(Hf2(t))

		"Measurement -------------------------------------------------------------------------------------------------"
		c = if (measurementtype == "excitation number")
				n
			elseif (measurementtype == "full parity")
				σx1 * σx2
			elseif (measurementtype == "half parity, x1 - x2")
				σx1 - σx2
			elseif (measurementtype == "half parity, x1 + x2")
				σx1 + σx2
			else
				error("Please enter a valid measurement type: 'excitation number', 'full parity', 'half parity, x1 - x2', or 'half parity, x1 + x2'.")
			end

		C = (Γm == 0.0) ? [] : [(c, Γm, 1.0)]

		"Loop over realizations -------------------------------------------------------------------------------------------------"
		solutions = ensemble(bayesian, (0.0, tf), ψ0, H(Hf1(τi, tf), Hf2(τi, tf)), [], C; dt=dt, N=N, batch_size=batch_size)

		# solutions = @showprogress pmap(i -> begin
		#
		# 	Hi1 = Hf1(τi, tf) # new realization of fluctuator Hamiltonian (single qubit)
		# 	Hi2 = Hf1(τi, tf) # new realization of fluctuator Hamiltonian (single qubit)
		# 	return bayesian((0, tf), ψ0, H(Hi1, Hi2), [], C; dt=dt)
		#
	    # end, 1:N; batch_size=batch_size)

		"Process and export data"
		# get expectation values
		ops = bell_basis
		op_labels = bell_basis_strings
		exps = map(op -> [expectations(sol, op) for sol in solutions], ops)
		pur = average_purity(solutions)

		# export data
		if write_exp_values
			for (arr, label) in zip(exps, op_labels)
				df = DataFrame(arr, :auto)
				CSV.write(string(ep, label, "-trajectory.csv"), df)
			end
			df = DataFrame([pur], :auto)
			CSV.write(string(ep, "purity-trajectory.csv"), df)
		end
end # function twoqubit_fluctuators

"""
sgn(τ::Timescale, tf::Timescale)
generates a function returning sign of fluctuator of switching time τ at time t before tf
"""
function sgn(τ::Timescale, tf::Timescale)
	nflips = Int64(floor(tf / τ)) # number of times fluctuator CAN flip in simulation
	series = rand((-1, 1), (nflips + 1)) # string of random signs of fluctuator

	t -> let
		index = Int64(floor(t/tf * nflips)) + 1
		return series[index]
	end
end


"""
twoqubitfeedback_fluctuators(; write_exp_values=true, N = 1000, batch_size = 10, Γ = (2π) * 0.2, dt = 1e-3, tf = 10.0,
									td = 0.0, Δt = 100e-3, tdd = 0.0, Ω = (2π) * 0.2, ΩR = (2π) * 25.0, nfluctuators = 10,
									φd1 = π/2, φd2 = 3π/2, comment = "")

Keyword Arguments:
N :: Int64	 				-- number of trajectories
batch_size :: Int64			-- number of trajectories allocated to a worker at a time
Γ :: Rate					-- measurement rate,
dt :: Timescale 			-- integration time step
tf :: Timescale				-- simulation duration
td :: Timescale				-- time delay for feedback
Δt :: Timescale				-- Rabi pulse duration
tdd :: Timescale			-- time at which dd turns on
Ω :: Rate					-- fluctuator rate
ΩR :: Rate					-- dynamical decoupling rate
nfluctuators :: Int64		-- number of fluctuators
φd1 :: Float64				-- decoupling axis on qubit 1
φd2 :: Float64				-- decoupling axis on qubit 2
comment :: String; 			-- comment for parameter file for future reference
"""
function twoqubitfeedback_fluctuators(; write_exp_values=true, N = 1000, batch_size = 10, Γ = (2π) * 0.2, dt = 1e-3, tf = 10.0,
									td = 0.0, Δt = 100e-3, tdd = 0.0, Ω = (2π) * 0.2, ΩR = (2π) * 25.0, nfluctuators = 10,
									φd1 = π/2, φd2 = 3π/2, comment = "")
		"Export parameters -------------------------------------------------------------------------------------"

		ep = exportpath("two-qubit-feedback-fluctuators")
		mkpath(ep)
		println(string("Writing to path ", ep))

		pars = Dict("N" => N,
					"Γ" => Γ,
					"dt" => dt,
					"tf" => tf,
					"td" => td,
					"Δt" => Δt,
					"Ω" => Ω, # fluctuator strength
					"ΩR" => ΩR, # decoupling strength
					"φd1" => φd1, # decoupling axis, qubit 1
					"φd2" => φd2, # decoupling axis, qubit 2
					"nfluctuators" => nfluctuators, # number of fluctuators
					"comment" => comment
					)

		CSV.write(string(ep, "parameters.csv"), pars)

		"System parameters -------------------------------------------------------------------------------------------------"
		# initial state
		ρ0 = g ⊗ g
		ρtarget = Ψp
		targetstr = "Ψp"

		"Setting up the fluctuators ---------------------------------------------------------------------------------------"
		# sample log-uniformly
		(ωmin, ωmax) = 1/tf, 1/dt
		ωlog = range(log(ωmin), log(ωmax), length=nfluctuators)
		ωi = exp.(ωlog)
		τi = map(ω -> 1/ω, ωi) # fluctuators and their signs for all $nfluctuators switching times

		"Hamiltonians -------------------------------------------------------------------------------------------------"
		"Fluctuator Hamiltonians"
		# fluctuator Hamiltonians for each qubit
		Hf1(τi, tf) = let si = map(τ -> sgn(τ, tf), τi)
						t -> map(s -> (Ω/2) * σz1 * s(t), si) end

		Hf2(τi, tf) = let si = map(τ -> sgn(τ, tf), τi)
						t -> map(s -> (Ω/2) * σz2 * s(t), si) end

		# decoupling Hamiltonian; turn off decoupling by setting ΩR = 0
		H2q(t) = (t < tdd) ? Iq ⊗ Iq : (ΩR/2) * (cos(φd1) * σx1 + sin(φd1) * σy1 + cos(φd2) * σx2 + sin(φd2) * σy2)

		# total Hamiltonian -- a function of a particular realization of fluctuator Hamiltonian Hf
		Hz(Hf1, Hf2) = t -> sum(Hf1(t)) + sum(Hf2(t))

		"Feedback Hamiltonian"
		function θ(ρ::State)
			αβ = real(expect(Φp ⊗ Φm', ρ))
			return 0.5 * asin(2αβ)
		end

		Hf(t::Timescale, ρ::State) = let
				 Ω =
					if t < Δt
						π / (4Δt)
					elseif td == 0.0
						- θ(ρ) / (2dt) # try θ^3
					else
						- θ(ρ) / (8td)
					end

				Ω * (σy1 + σy2) + H2q(t)

			end

		"Total system Hamiltonian"
		Hs(Hf1, Hf2) = (t::Timescale, ρ::State) -> Hf(t, ρ) + Hz(Hf1, Hf2)(t) + H2q(t)
		C = [(n, Γ, 1.0)]

		"Run parallelized code -------------------------------------------------------------------------------------------------"
		ops = bell_basis
		op_labels = bell_basis_strings
		(system_exps, filter_exps) = ensemble(bayesian, (0.0, tf), ρ0, (Hs(Hf1(τi, tf), Hf2(τi, tf)), Hf), [], C; dt=dt, N=N, batch_size=batch_size, ops=ops)

		"Extract and write expectation values -----------------------------------------------------------------------------------"
		if write_exp_values
			# write data
			for (arr, label) in zip(system_exps, op_labels)
				df = DataFrame(arr, :auto)
				CSV.write(string(ep, "system-", label, "-trajectory.csv"), df)
			end

			for (arr, label) in zip(filter_exps, op_labels)
				df = DataFrame(arr, :auto)
				CSV.write(string(ep, "filter-", label, "-trajectory.csv"), df)
			end
		end
end # function twoqubitfeedback_fluctuators
"""
twoqubitfeedback_fluctuators_fixed(; write_exp_values=true, N = 1000, batch_size = 10, Γ = (2π) * 0.2, dt = 1e-3, tf = 10.0,
									td = 0.0, Δt = 100e-3, tdd = 0.0, Ω = (2π) * 0.2, ΩR = (2π) * 25.0, nfluctuators = 10,
									φd1 = π/2, φd2 = 3π/2, comment = "")

Modified compared to twoqubitfeedback_fluctuators, which may be calculating θ incorrectly, this function calculates
according to Leigh's paper.

Keyword Arguments:
N :: Int64	 				-- number of trajectories
batch_size :: Int64			-- number of trajectories allocated to a worker at a time
Γ :: Rate					-- measurement rate,
dt :: Timescale 			-- integration time step
tf :: Timescale				-- simulation duration
td :: Timescale				-- time delay for feedback
Δt :: Timescale				-- Rabi pulse duration
tdd :: Timescale			-- time at which dd turns on
Ω :: Rate					-- fluctuator rate
ΩR :: Rate					-- dynamical decoupling rate
nfluctuators :: Int64		-- number of fluctuators
φd1 :: Float64				-- decoupling axis on qubit 1
φd2 :: Float64				-- decoupling axis on qubit 2
comment :: String; 			-- comment for parameter file for future reference
"""
function twoqubitfeedback_fluctuators_fixed(; write_exp_values=true, N = 1000, batch_size = 10, Γ = (2π) * 0.2, dt = 1e-3, tf = 10.0,
									td = 0.0, Δt = 100e-3, tdd = 0.0, Ω = (2π) * 0.2, ΩR = (2π) * 25.0, nfluctuators = 10,
									φd1 = π/2, φd2 = 3π/2, feedback=true, comment = "")
		"Export parameters -------------------------------------------------------------------------------------"

		ep = exportpath("two-qubit-feedback-fluctuators-fixed")
		mkpath(ep)
		println(string("Writing to path ", ep))

		pars = Dict("N" => N,
					"Γ" => Γ,
					"dt" => dt,
					"tf" => tf,
					"td" => td,
					"Δt" => Δt,
					"Ω" => Ω, # fluctuator strength
					"ΩR" => ΩR, # decoupling strength
					"φd1" => φd1, # decoupling axis, qubit 1
					"φd2" => φd2, # decoupling axis, qubit 2
					"nfluctuators" => nfluctuators, # number of fluctuators
					"comment" => comment
					)

		CSV.write(string(ep, "parameters.csv"), pars)

		"System parameters -------------------------------------------------------------------------------------------------"
		# initial state
		ρ0 = normalize(g + e) ⊗ normalize(g + e) # Ψp
		ρtarget = Ψp
		targetstr = "Ψp"

		"Setting up the fluctuators ---------------------------------------------------------------------------------------"
		# sample log-uniformly
		(ωmin, ωmax) = 1/tf, 1/dt
		ωlog = range(log(ωmin), log(ωmax), length=nfluctuators)
		ωi = exp.(ωlog)
		τi = map(ω -> 1/ω, ωi) # fluctuators and their signs for all $nfluctuators switching times

		"Hamiltonians -------------------------------------------------------------------------------------------------"
		"Fluctuator Hamiltonians"
		# fluctuator Hamiltonians for each qubit
		function sgn(τ::Timescale, tf::Timescale)
			nflips = Int64(floor(tf / τ)) # number of times fluctuator CAN flip in simulation
			series = rand((-1, 1), (nflips + 1)) # string of random signs of fluctuator

			t -> let
				index = Int64(floor(t/tf * nflips)) + 1
				return series[index]
			end
		end

		Hf1(τi, tf) = let si = map(τ -> sgn(τ, tf), τi)
						t -> map(s -> (Ω/2) * σz1 * s(t), si) end

		Hf2(τi, tf) = let si = map(τ -> sgn(τ, tf), τi)
						t -> map(s -> (Ω/2) * σz2 * s(t), si) end

		# decoupling Hamiltonian; turn off decoupling by setting ΩR = 0
		H2q(t) = (t < tdd) ? Iq ⊗ Iq : (ΩR/2) * (cos(φd1) * σx1 + sin(φd1) * σy1 + cos(φd2) * σx2 + sin(φd2) * σy2)

		# total Hamiltonian -- a function of a particular realization of fluctuator Hamiltonian Hf
		Hz(Hf1, Hf2) = t -> sum(Hf1(t)) + sum(Hf2(t))

		"Feedback Hamiltonian"
		function arctan(y, x)
			if x > 0
				atan(y / x)
			elseif y >= 0 && x < 0
				atan(y / x) + π
			elseif y < 0 && x < 0
				atan(y / x) - π
			else
				sign(y) * π/2
			end
		end

		function θ(ρ::State)
			tm = g ⊗ g
			tp = e ⊗ e
			s = Ψm
			t0 = Ψp

			tmt0 = projector(tm, dagger(t0))
			t0tp = projector(t0, dagger(tp))
			ss = projector(s)
			tmtp = projector(tm, dagger(tp))

			F = expect(dm(t0), ρ)
			ρtmt0 = expect(tmt0, ρ)
			ρt0tp = expect(t0tp, ρ)
			ρss = expect(ss, ρ)
			ρtmtp = expect(tmtp, ρ)

			y = real(√8 * (ρtmt0 - ρt0tp))
			x = real(3F + ρss + 2ρtmtp - 1)

			return 0.5 * arctan(y, x)
		end

		Hf(t::Timescale, ρ::State) = feedback ? - θ(ρ) / (2dt) * (σy1 + σy2) + H2q(t) : H2q(t)

		"Total system Hamiltonian"
		Hs(Hf1, Hf2) = (t::Timescale, ρ::State) -> 	Hf(t, ρ) + Hz(Hf1, Hf2)(t) + H2q(t)
		C = feedback ? [(n, Γ, 1.0)] : []

		"Run parallelized code -------------------------------------------------------------------------------------------------"
		ops = bell_basis
		op_labels = bell_basis_strings

		solutions = @showprogress pmap(i -> begin

			Hi1 = Hf1(τi, tf) # new realization of fluctuator Hamiltonian (single qubit)
			Hi2 = Hf1(τi, tf) # new realization of fluctuator Hamiltonian (single qubit)
			return bayesian((0.0, tf), ρ0, (Hs(Hi1, Hi2), Hf), [], C; dt=dt)
		end, 1:N; batch_size=batch_size)

		# (system_exps, filter_exps) = ensemble(bayesian, (0.0, tf), ρ0, (Hs(Hf1(τi, tf), Hf2(τi, tf)), Hf), [], C; dt=dt, N=N, batch_size=batch_size, ops=ops)

		system_exps = map(op -> [expectations(sol[1], op) for sol in solutions], ops)
		filter_exps = map(op -> [expectations(sol[2], op) for sol in solutions], ops)

		"Extract and write expectation values -----------------------------------------------------------------------------------"
		if write_exp_values
			# write data
			for (arr, label) in zip(system_exps, op_labels)
				df = DataFrame(arr, :auto)
				CSV.write(string(ep, "system-", label, "-trajectory.csv"), df)
			end

			for (arr, label) in zip(filter_exps, op_labels)
				df = DataFrame(arr, :auto)
				CSV.write(string(ep, "filter-", label, "-trajectory.csv"), df)
			end
		end
end # function twoqubitfeedback_fluctuators



"""
twoqubitfeedback_fluctuators_CR(; write_exp_values=true, N = 1000, batch_size = 10, Γ = (2π) * 0.2, dt = 1e-3, tf = 10.0,
									td = 0.0, Δt = 100e-3, tdd = 0.0, Ω = (2π) * 0.2, ΩR = (2π) * 25.0, nfluctuators = 10,
									φd1 = π/2, φd2 = 3π/2, comment = "")

Calculates optimal rotation angles according to Leigh's paper draft (θ and Δθ)

Keyword Arguments:
N :: Int64	 				-- number of trajectories
batch_size :: Int64			-- number of trajectories allocated to a worker at a time
ΓH :: Rate					-- measurement rate,
dt :: Timescale 			-- integration time step
tf :: Timescale				-- simulation duration
td :: Timescale				-- time delay for feedback
Ω :: Rate					-- fluctuator rate (per fluctuator)
Ωmax :: Rate				-- maximum counterrotating drive
nfluctuators :: Int64		-- number of fluctuators
comment :: String; 			-- comment for parameter file for future reference
"""
function twoqubitfeedback_fluctuators_CR(; write_trajectories=true, counterrotating=true, N = 1000, batch_size = 10,
									ΓH = 0.5, dt = 1e-3, tf = 10.0, td = 0.0, Ω = (2π) * 0.2, Ωmax = (2π) * 50.0,
									nfluctuators = 10, comment = "")
		"Export parameters -------------------------------------------------------------------------------------"

		ep = exportpath("two-qubit-feedback-fluctuators-CR")
		mkpath(ep)
		println(string("Writing to path ", ep))

		pars = Dict("N" => N,
					"ΓH" => ΓH, # measurement strength
					"dt" => dt, # simulation time step
					"tf" => tf, # simulation duration
					"td" => td, # time delay for feedback
					"Ω" => Ω, # fluctuator strength (per fluctuator)
					"Ωmax" => Ωmax, # maximum Rabi drive
					"nfluctuators" => nfluctuators, # number of fluctuators
					"write_trajectories" => write_trajectories,
					"comment" => comment
					)

		CSV.write(string(ep, "parameters.csv"), pars)

		"System parameters -------------------------------------------------------------------------------------------------"
		# initial state
		ρ0 = normalize((g + e) ⊗ (g + e))

		"Setting up the fluctuators ---------------------------------------------------------------------------------------"
		# sample log-uniformly
		(ωmin, ωmax) = 1/tf, 1/dt
		ωlog = range(log(ωmin), log(ωmax), length=nfluctuators)
		ωi = exp.(ωlog)
		τi = map(ω -> 1/ω, ωi) # fluctuators and their signs for all $nfluctuators switching times

		"Fluctuator Hamiltonians"
		# fluctuator Hamiltonians for each qubit
		function sgn(τ::Timescale, tf::Timescale)
			nflips = Int64(floor(tf / τ)) # number of times fluctuator CAN flip in simulation
			series = rand((-1, 1), (nflips + 1)) # string of random signs of fluctuator

			t -> let
				index = Int64(floor(t/tf * nflips)) + 1
				return series[index]
			end
		end

		Hf1(τi, tf) = let si = map(τ -> sgn(τ, tf), τi)
						t -> map(s -> (Ω/2) * σz1 * s(t), si) end

		Hf2(τi, tf) = let si = map(τ -> sgn(τ, tf), τi)
						t -> map(s -> (Ω/2) * σz2 * s(t), si) end

		# total fluctuator Hamiltonian -- a function of a particular realization of fluctuator Hamiltonian Hf
		Hz(Hf1, Hf2) = t -> sum(Hf1(t)) + sum(Hf2(t))

		"Rotating and counterrotating Hamiltonians ---------------------------------------------------------------------------------------"
		"Helpers"

		"Arctan(y, x) Gives the angle defined by coordinates (x,y) in the correct quadrant."
		function arctan(y, x)
			if x > 0
				atan(y / x)
			elseif y >= 0 && x < 0
				atan(y / x) + π
			elseif y < 0 && x < 0
				atan(y / x) - π
			else
				sign(y) * π/2
			end
		end

		"two-qubit projectors"
		ψpψp = projector(Ψp)
		ϕmϕm = projector(Φm)
		ϕmψp = projector(Φm, dagger(Ψp))

		ϕpϕp = projector(Φp)
		ψmψm = projector(Ψm)
		ψmϕp = projector(Ψm, dagger(Φp))

		"Rotation angle"
		function θ(ρ::State)
			ρψpψp = expect(ψpψp, ρ)
			ρϕmϕm = expect(ϕmϕm, ρ)
			ρψpϕm = expect(ϕmψp, ρ)

			x = real(ρψpψp - ρϕmϕm)
			y = 2 * real(ρψpϕm)

			return 0.5 * arctan(y, x)
		end

		magθ(ρ) = 	let inst = - θ(ρ) / (2dt)
						(abs(inst) < Ωmax) ? inst : sign(inst) * Ωmax
					end

		"Counterrotation angle"
		function Δθ(ρ::State)
			ρϕpϕp = expect(ϕpϕp, ρ)
			ρψmψm = expect(ψmψm, ρ)
			ρϕpψm = expect(ψmϕp, ρ)

			y = real(ρϕpϕp - ρψmψm)
			x = 2 * real(ρϕpψm)

			return 0.5 * (arctan(y, x) + π)
		end

		magΔθ(ρ) = 	let inst = - Δθ(ρ) / (2dt)
						(abs(inst) < Ωmax) ? inst : sign(inst) * Ωmax
					end

		"Feedback (rotating) Hamiltonian"
		HF(t::Timescale, ρ::State) = (magθ(ρ)/2) * (σy1 + σy2)

		"Draining (counterrotating) Hamiltonian"
		HS(ρ::State) = counterrotating ? (magΔθ(ρ)/2) * (σy1 - σy2) : Iq ⊗ Iq


		"Total Hamiltonians -------------------------------------------------------------------------------------------------"
		Hfil(t::Timescale, ρ::State) = HF(t, ρ) + HS(ρ) # filter
		Hsys(Hf1, Hf2) = (t::Timescale, ρ::State) -> Hfil(t, ρ) + Hz(Hf1, Hf2)(t) # system

		C = [((σz1 + σz2)/√2, ΓH, 1.0)]

		"Run parallelized code -------------------------------------------------------------------------------------------------"
		ops = bell_basis
		op_labels = bell_basis_strings

		solutions = @showprogress pmap(i -> begin

			Hi1 = Hf1(τi, tf) # new realization of fluctuator Hamiltonian (single qubit)
			Hi2 = Hf1(τi, tf) # new realization of fluctuator Hamiltonian (single qubit)

			return bayesian((0.0, tf), ρ0, (Hsys(Hi1, Hi2), Hfil), [], C; dt=dt)

		end, 1:N; batch_size=batch_size)

		sys_sols = map(sol -> sol[1], solutions)
		fil_sols = map(sol -> sol[2], solutions)

		sys_avg = ensemble_average(sys_sols)
		fil_avg = ensemble_average(fil_sols)

		"Extract and write expectation values -----------------------------------------------------------------------------------"
		"Ensemble averages"
		for (avg, label) in zip([sys_avg, fil_avg], ["system", "filter"])

			fids = map(ρ -> fidelity(ρ, Ψp), avg)
			exps = [map(ρ -> real(expect(ρ, op)), avg) for op in ops]

			df_fids = DataFrame([fids], :auto)
			df_exps = DataFrame(exps, :auto)

			CSV.write(string(ep, label, "-fidelity-trajectory.csv"), df_fids)
			CSV.write(string(ep, label, "-average-populations.csv"), df_exps)

		end

		"Trajectories"
		if write_trajectories

			system_exps = map(op -> [expectations(sol, op) for sol in sys_sols], ops)
			filter_exps = map(op -> [expectations(sol, op) for sol in fil_sols], ops)

			# write data
			for (arr, label) in zip(system_exps, op_labels)
				df = DataFrame(arr, :auto)
				CSV.write(string(ep, "system-", label, "-trajectory.csv"), df)
			end

			for (arr, label) in zip(filter_exps, op_labels)
				df = DataFrame(arr, :auto)
				CSV.write(string(ep, "filter-", label, "-trajectory.csv"), df)
			end
		end
end # function twoqubitfeedback_fluctuators


function twoqubitfeedback_fluctuators_Zeno(; write_trajectories=true, counterrotating = false, N = 1000, batch_size = 10,
									ΓH = 0.5, ΓZ = 10.0, dt = 1e-3, tf = 10.0, td = 0.0, tZ = 10.0,
									Ω = (2π) * 0.2, Ωmax = (2π) * 50.0,
									ηH = 1.0, ηZ = 1.0,
									nfluctuators = 10, measurementtype ="full parity",
									comment = "")
		"Export parameters -------------------------------------------------------------------------------------"

		ep = exportpath("two-qubit-feedback-fluctuators-Zeno")
		mkpath(ep)
		println(string("Writing to path ", ep))

		pars = Dict("N" => N,
					"ΓH" => ΓH, # feedback measurement strength
					"ΓZ" => ΓZ, # Zeno measurement strength
					"dt" => dt, # simulation time step
					"tf" => tf, # simulation duration
					"td" => td, # time delay for feedback
					"tZ" => tZ, # time delay for turning on Zeno measurement
					"Ω" => Ω, # fluctuator strength (per fluctuator)
					"Ωmax" => Ωmax, # maximum Rabi drive
					"ηH" => ηH, # feedback measurement quantum efficiency
					"ηZ" => ηZ, # Zeno measurement quantum efficiency
					"nfluctuators" => nfluctuators, # number of fluctuators
					"write_trajectories" => write_trajectories,
					"measurementtype" => measurementtype, # operator measured for testing Zeno effect
					"counterrotating" => counterrotating,
					"comment" => comment
					)

		CSV.write(string(ep, "parameters.csv"), pars)

		"System parameters -------------------------------------------------------------------------------------------------"
		# initial state
		ρ0 = normalize((g + e) ⊗ (g + e))

		"Setting up the fluctuators ---------------------------------------------------------------------------------------"
		# sample log-uniformly
		(ωmin, ωmax) = 1/tf, 1/dt
		ωlog = range(log(ωmin), log(ωmax), length=nfluctuators)
		ωi = exp.(ωlog)
		τi = map(ω -> 1/ω, ωi) # fluctuators and their signs for all $nfluctuators switching times

		"Fluctuator Hamiltonians"
		# fluctuator Hamiltonians for each qubit
		function sgn(τ::Timescale, tf::Timescale)
			nflips = Int64(floor(tf / τ)) # number of times fluctuator CAN flip in simulation
			series = rand((-1, 1), (nflips + 1)) # string of random signs of fluctuator

			t -> let
				index = Int64(floor(t/tf * nflips)) + 1
				return series[index]
			end
		end

		Hf1(τi, tf) = let si = map(τ -> sgn(τ, tf), τi)
						t -> map(s -> (Ω/2) * σz1 * s(t), si) end

		Hf2(τi, tf) = let si = map(τ -> sgn(τ, tf), τi)
						t -> map(s -> (Ω/2) * σz2 * s(t), si) end

		# total fluctuator Hamiltonian -- a function of a particular realization of fluctuator Hamiltonian Hf
		Hz(Hf1, Hf2) = t -> sum(Hf1(t)) + sum(Hf2(t))

		"Rotating and counterrotating Hamiltonians ---------------------------------------------------------------------------------------"
		"Helpers"

		"Arctan(y, x) Gives the angle defined by coordinates (x,y) in the correct quadrant."
		function arctan(y, x)
			if x > 0
				atan(y / x)
			elseif y >= 0 && x < 0
				atan(y / x) + π
			elseif y < 0 && x < 0
				atan(y / x) - π
			else
				sign(y) * π/2
			end
		end

		"Rotation angle"
		function θ(ρ::State)
			tm = g ⊗ g
			tp = e ⊗ e
			s = Ψm
			t0 = Ψp

			tmt0 = projector(tm, dagger(t0))
			t0tp = projector(t0, dagger(tp))
			ss = projector(s)
			tmtp = projector(tm, dagger(tp))

			F = expect(dm(t0), ρ)
			ρtmt0 = expect(tmt0, ρ)
			ρt0tp = expect(t0tp, ρ)
			ρss = expect(ss, ρ)
			ρtmtp = expect(tmtp, ρ)

			y = real(√8 * (ρtmt0 - ρt0tp))
			x = real(3F + ρss + 2ρtmtp - 1)

			return 0.5 * arctan(y, x)
		end

		magθ(ρ) = 	let inst = - θ(ρ) / (2dt)
						(abs(inst) < Ωmax) ? inst : sign(inst) * Ωmax
					end


		"Feedback (rotating) Hamiltonian"
		HF(t::Timescale, ρ::State) = (magθ(ρ)/2) * (σy1 + σy2)

		"Total Hamiltonians -------------------------------------------------------------------------------------------------"
		Hfil(t::Timescale, ρ::State) = HF(t, ρ) # filter
		Hsys(Hf1, Hf2) = (t::Timescale, ρ::State) -> Hfil(t, ρ) + Hz(Hf1, Hf2)(t) # system

		"Measurement -----------------------------"

		c =
			if (measurementtype == "full parity")
				σx1 * σx2
			elseif (measurementtype == "half parity, x1 - x2")
				σx1 - σx2
			elseif (measurementtype == "half parity, x1 + x2")
				σx1 + σx2
			else
				error("Please enter a valid measurement type: 'excitation number', 'full parity', 'half parity, x1 - x2', or 'half parity, x1 + x2'.")
			end

		cZ(t) = (t < tZ) ? Iq ⊗ Iq : c


		C = [	((σz1 + σz2)/√2, ΓH, ηH), 	# excitation number measurement for feedback
				(cZ, ΓZ, ηZ)] 				# (half) parity measurement for Zeno pinning

		J = []	# Lindblad dissipation due to inefficient measurement
		for (c, Γ, η) in C
			if η < 1.0
				push!(J, (c, (1 - η) * Γ))
			end
		end

		"Run parallelized code -------------------------------------------------------------------------------------------------"
		ops = bell_basis
		op_labels = bell_basis_strings

		solutions = @showprogress pmap(i -> begin

			Hi1 = Hf1(τi, tf) # new realization of fluctuator Hamiltonian (single qubit)
			Hi2 = Hf1(τi, tf) # new realization of fluctuator Hamiltonian (single qubit)

			return bayesian((0.0, tf), ρ0, (Hsys(Hi1, Hi2), Hfil), J, C; dt=dt)

		end, 1:N; batch_size=batch_size)

		sys_sols = map(sol -> sol[1], solutions)
		fil_sols = map(sol -> sol[2], solutions)

		sys_avg = ensemble_average(sys_sols)
		fil_avg = ensemble_average(fil_sols)

		"Extract and write expectation values -----------------------------------------------------------------------------------"
		"Ensemble averages"
		for (avg, label) in zip([sys_avg, fil_avg], ["system", "filter"])

			fids = map(ρ -> fidelity(ρ, Ψp), avg)
			exps = [map(ρ -> real(expect(ρ, op)), avg) for op in ops]

			df_fids = DataFrame([fids], :auto)
			df_exps = DataFrame(exps, :auto)

			CSV.write(string(ep, label, "-fidelity-trajectory.csv"), df_fids)
			CSV.write(string(ep, label, "-average-populations.csv"), df_exps)

		end

		"Trajectories"
		if write_trajectories

			system_exps = map(op -> [expectations(sol, op) for sol in sys_sols], ops)
			filter_exps = map(op -> [expectations(sol, op) for sol in fil_sols], ops)

			# write data
			for (arr, label) in zip(system_exps, op_labels)
				df = DataFrame(arr, :auto)
				CSV.write(string(ep, "system-", label, "-trajectory.csv"), df)
			end

			for (arr, label) in zip(filter_exps, op_labels)
				df = DataFrame(arr, :auto)
				CSV.write(string(ep, "filter-", label, "-trajectory.csv"), df)
			end
		end
end # function twoqubitfeedback_fluctuators_Zeno




function average_purity(solutions::Vector{Solution})
	return map(ρ -> real(tr(ρ * ρ)), ensemble_average(solutions))
end

function ensemble_average(solutions::Vector{Solution})

	if typeof(solutions[1].ρ[1]) <: Ket
		return [mean(map(sol -> dm(sol.ρ[i]), solutions)) for i in 1:length(solutions[1].t)]
	else
		return [mean(map(sol -> sol.ρ[i], solutions)) for i in 1:length(solutions[1].t)]
	end

end


function twoqubitfeedback_fluctuators_Zeno_feedback(; write_trajectories=true, N = 1000, batch_size = 10,
									ΓH = 0.5, ΓZ = 20.0, dt = 1e-3, tf = 20.0, td = 0.0, tZ = 10.0,
									Ω = (2π) * 0.2, Ωmax = (2π) * 50.0, ΩZ = (2π) * 0.8, ηH = 1.0, ηZ = 1.0,
									readout_integration_bins = 500, nfluctuators = 10, comment = "")
		"Export parameters -------------------------------------------------------------------------------------"

		ep = exportpath("two-qubit-feedback-fluctuators-Zeno-feedback")
		mkpath(ep)
		println(string("Writing to path ", ep))

		pars = Dict("N" => N,
					"ΓH" => ΓH, # feedback measurement strength
					"ΓZ" => ΓZ, # Zeno measurement strength
					"dt" => dt, # simulation time step
					"tf" => tf, # simulation duration
					"td" => td, # time delay for feedback
					"tZ" => tZ, # time delay for turning on Zeno measurement
					"Ω" => Ω, # fluctuator strength (per fluctuator)
					"Ωmax" => Ωmax, # maximum Rabi drive
					"ΩZ" => ΩZ, # Zeno feedback strength
					"ηH" => ηH, # feedback measurement quantum efficiency
					"ηZ" => ηZ, # Zeno measurement quantum efficiency
					"readout_integration_bins" => readout_integration_bins, # number of bins to smooth readout over in feedback
					"nfluctuators" => nfluctuators, # number of fluctuators
					"write_trajectories" => write_trajectories,
					"comment" => comment
					)

		CSV.write(string(ep, "parameters.csv"), pars)

		"System parameters -------------------------------------------------------------------------------------------------"
		# initial state
		ρ0 = normalize((g + e) ⊗ (g + e))

		"Setting up the fluctuators ---------------------------------------------------------------------------------------"
		# sample log-uniformly
		(ωmin, ωmax) = 1/tf, 1/dt
		ωlog = range(log(ωmin), log(ωmax), length=nfluctuators)
		ωi = exp.(ωlog)
		τi = map(ω -> 1/ω, ωi) # fluctuators and their signs for all $nfluctuators switching times

		"Fluctuator Hamiltonians"
		# fluctuator Hamiltonians for each qubit
		function sgn(τ::Timescale, tf::Timescale)
			nflips = Int64(floor(tf / τ)) # number of times fluctuator CAN flip in simulation
			series = rand((-1, 1), (nflips + 1)) # string of random signs of fluctuator

			t -> let
				index = Int64(floor(t/tf * nflips)) + 1
				return series[index]
			end
		end

		Hf1(τi, tf) = let si = map(τ -> sgn(τ, tf), τi)
						t -> map(s -> (Ω/2) * σz1 * s(t), si) end

		Hf2(τi, tf) = let si = map(τ -> sgn(τ, tf), τi)
						t -> map(s -> (Ω/2) * σz2 * s(t), si) end

		# total fluctuator Hamiltonian -- a function of a particular realization of fluctuator Hamiltonian Hf
		Hz(Hf1, Hf2) = t -> sum(Hf1(t)) + sum(Hf2(t))

		"Rotating and counterrotating Hamiltonians ---------------------------------------------------------------------------------------"
		"Helpers"

		"Arctan(y, x) Gives the angle defined by coordinates (x,y) in the correct quadrant."
		function arctan(y, x)
			if x > 0
				atan(y / x)
			elseif y >= 0 && x < 0
				atan(y / x) + π
			elseif y < 0 && x < 0
				atan(y / x) - π
			else
				sign(y) * π/2
			end
		end

		"Rotation angle"
		function θ(ρ::State)
			tm = g ⊗ g
			tp = e ⊗ e
			s = Ψm
			t0 = Ψp

			tmt0 = projector(tm, dagger(t0))
			t0tp = projector(t0, dagger(tp))
			ss = projector(s)
			tmtp = projector(tm, dagger(tp))

			F = expect(dm(t0), ρ)
			ρtmt0 = expect(tmt0, ρ)
			ρt0tp = expect(t0tp, ρ)
			ρss = expect(ss, ρ)
			ρtmtp = expect(tmtp, ρ)

			y = real(√8 * (ρtmt0 - ρt0tp))
			x = real(3F + ρss + 2ρtmtp - 1)

			return 0.5 * arctan(y, x)
		end

		magθ(ρ) = 	let inst = - θ(ρ) / (2dt)
						(abs(inst) < Ωmax) ? inst : sign(inst) * Ωmax
					end

		"Feedback (rotating) Hamiltonian"
		HF(t::Timescale, ρ::State) = (magθ(ρ)/2) * (σy1 + σy2)

		"Zeno Hamiltonian"
		HZ(t::Timescale, r::Record) = t < tZ ? Iq ⊗ Iq : (ΩZ/2) * (r[2] - 1) * σz1


		"Total Hamiltonians -------------------------------------------------------------------------------------------------"
		Hfil(t::Timescale, ρ::State, r::Record) = HF(t, ρ) + HZ(t, r) # filter
		Hsys(Hf1, Hf2) = (t::Timescale, ρ::State, r::Record) -> Hfil(t, ρ, r) + Hz(Hf1, Hf2)(t) # system

		"Measurement -----------------------------"

		c = σx1 * σx2
		cZ(t) = (t < tZ) ? Iq ⊗ Iq : c

		C = [	((σz1 + σz2)/√2, ΓH, ηH), 	# excitation number measurement for feedback
						(cZ, ΓZ, ηZ)] 				# (half) parity measurement for Zeno pinning


		J = []	# Lindblad dissipation due to inefficient measurement
		for (c, Γ, η) in C
			if η < 1.0
				push!(J, (c, (1 - η) * Γ))
			end
		end


		"Run parallelized code -------------------------------------------------------------------------------------------------"
		ops = bell_basis
		op_labels = bell_basis_strings

		solutions = @showprogress pmap(i -> begin

			Hi1 = Hf1(τi, tf) # new realization of fluctuator Hamiltonian (single qubit)
			Hi2 = Hf1(τi, tf) # new realization of fluctuator Hamiltonian (single qubit)

			return bayesian((0.0, tf), ρ0, (Hsys(Hi1, Hi2), Hfil), J, C; dt=dt, readout_integration_bins = readout_integration_bins)

		end, 1:N; batch_size=batch_size)

		sys_sols = map(sol -> sol[1], solutions)
		fil_sols = map(sol -> sol[2], solutions)

		sys_avg = ensemble_average(sys_sols)
		fil_avg = ensemble_average(fil_sols)

		"Extract and write expectation values -----------------------------------------------------------------------------------"
		"Ensemble averages"
		for (avg, label) in zip([sys_avg, fil_avg], ["system", "filter"])

			fids = map(ρ -> fidelity(ρ, Ψp), avg)
			exps = [map(ρ -> real(expect(ρ, op)), avg) for op in ops]

			df_fids = DataFrame([fids], :auto)
			df_exps = DataFrame(exps, :auto)

			CSV.write(string(ep, label, "-fidelity-trajectory.csv"), df_fids)
			CSV.write(string(ep, label, "-average-populations.csv"), df_exps)

		end

		"Trajectories"
		if write_trajectories

			system_exps = map(op -> [expectations(sol, op) for sol in sys_sols], ops)
			filter_exps = map(op -> [expectations(sol, op) for sol in fil_sols], ops)

			# write data
			for (arr, label) in zip(system_exps, op_labels)
				df = DataFrame(arr, :auto)
				CSV.write(string(ep, "system-", label, "-trajectory.csv"), df)
			end

			for (arr, label) in zip(filter_exps, op_labels)
				df = DataFrame(arr, :auto)
				CSV.write(string(ep, "filter-", label, "-trajectory.csv"), df)
			end
		end
end # function twoqubitfeedback_fluctuators_CR_Zeno
