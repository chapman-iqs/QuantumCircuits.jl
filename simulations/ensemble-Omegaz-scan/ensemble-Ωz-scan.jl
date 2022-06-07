N = 1000
write_trajectories = true

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
				ep = string("simulations/ensemble-Ωz-scan/data/", Dates.now(), "/")
				replace(ep, ":" => "-") end
mkpath(exportpath)


include("/Users/sachagreenfield/Desktop/GitHub/QC-notebooks/utilities.jl")

const Time = Float64
const Rate = Float64
const Efficiency = Float64
const Record = Vector{Float64}
const QOp = AbstractOperator
State = Union{Ket, QOp}
const Solution = QuantumCircuits.solution


"Qubit Hilbert space operators -------------------------------------------------------------------------------------"

# Basis
q = SpinBasis(1//2)
Iq = identityoperator(q)
I = Iq ⊗ Iq

# qubit operators, using convention that |-z> is ground state
σx = sigmax(q)
σy = sigmay(q)
σz = sigmaz(q)
σp = sigmap(q)
σm = sigmam(q)
nq = σp * σm

g = spindown(q)
e = spinup(q)

# two-qubit operators
σx1 = σx ⊗ Iq
σx2 = Iq ⊗ σx
σy1 = σy ⊗ Iq
σy2 = Iq ⊗ σy
σz1 = σz ⊗ Iq
σz2 = Iq ⊗ σz
σp1 = σp ⊗ Iq
σm1 = σm ⊗ Iq
σp2 = Iq ⊗ σp
σm2 = Iq ⊗ σm

# number operators
n1 = nq ⊗ Iq
n2 = Iq ⊗ nq
n = n1 + n2

# basis states
ket0 = g ⊗ g
ket1p = normalize(g ⊗ e + e ⊗ g)
ket1m = normalize(g ⊗ e - e ⊗ g)
ket2 = e ⊗ e

ketp = normalize(ket0 + ket2)
ketm = normalize(ket0 - ket2)

"""
System parameters -------------------------------------------------------------------------------------------------
all times given in μs
"""

# initial state
ρ0 = g ⊗ g
ρtarget = ket1p
targetstr = "ket1p"


"measurement / feedback parameters"

Γ = 1.0 # measurement rate
φ = π/2 # feedback drive angle
η = 1.0


"simulation timescales"

dt = 1e-3  # integration time-step
tf = 10.0 # simulation duration
td = 300e-3 # time delay for feedback
Δt = 100e-3 # time duration of initial π/2 pulse

"Kraus operators"
function θ(ρ::State)
	αβ = real(expect(ket1p ⊗ ketm', ρ))
	return 0.5 * asin(2αβ)
end

Hf(t::Time, ρ::State) = let
		 Ω =
			if t < Δt
				π / (4Δt)
			elseif td == 0.0
				- θ(ρ) / (2dt)
			else
				- θ(ρ) / (8td)
			end

		return Ω * (cos(φ) * (σx1 + σx2) + sin(φ) * (σy1 + σy2))

	end

Hz(t, Ωz) = let
		r1, r2 = rand() - 0.5, rand() - 0.5
		Ωz * (r1 * σz1 + r2 * σz2)
end

Hs(Ωz) = (t::Time, ρ::State) -> Hf(t, ρ) + Hz(t, Ωz)


J = η == 1.0 ? [] : [(n, ((1 - η) * Γ))]
C = [(n, Γ, η)]

"Loop over Ωz -------------------------------------------------------------"
Ωzs = [0.0, 5.0, 10.0, 40.0]# range(0.0, 40.0, step=0.5)
basis_states = dm.([ketp, ketm, ket1p, ket1m])
basis_strings = ["ketp", "ketm", "ket1p", "ket1m"]

function trajectory(sol::Solution, op::QOp)
	map(state -> real(expect(op, state)), sol.ρ)
end

### for handling fidelities
system_fid_avg = []
filter_fid_avg = []

# use the built-in QuantumOptics.fidelity function, but include cases to handle Kets
fidelity(ρ::Ket, σ::Ket) = real(ρ' * σ)^2
fidelity(ρ::Ket, σ::QOp) = fidelity(dm(ρ), σ)
fidelity(ρ::QOp, σ::Ket) = fidelity(ρ, dm(σ))

for Ωz in Ωzs

	println(string("Running Ωz = ", Ωz, " ..."))

	system_exps = [Vector{Float64}[] for op in basis_states]
	filter_exps = [Vector{Float64}[] for op in basis_states]
	system_fid = Float64[]
	filter_fid = Float64[]

	@showprogress for i in 1:N
		(sys, fil) = bayesian((0.0, tf), ρ0, (Hs(Ωz), Hf), J, C; dt=dt, td=td)

		# vector of 4 trajectories
		if write_trajectories
			for (arr, op) in zip(system_exps, basis_states)
				push!(arr, trajectory(sys, op)) end

			for (arr, op) in zip(filter_exps, basis_states)
				push!(arr, trajectory(fil, op)) end

		end

		push!(system_fid, fidelity(last(sys.ρ), ρtarget))
		push!(filter_fid, fidelity(last(fil.ρ), ρtarget))

	end

	push!(system_fid_avg, mean(system_fid))
	push!(filter_fid_avg, mean(filter_fid))


	"write trajectory data"
	if write_trajectories

		@debug "type of system_exps is " typeof(system_exps)

		# write
		for (exps, basis_label) in zip(system_exps, basis_strings)
			df = DataFrame(exps, :auto)
			CSV.write(string(exportpath, "Ωz-", Ωz, "_", basis_label, "-system-trajectory.csv"), df)
		end

		for (exps, basis_label) in zip(filter_exps, basis_strings)
			df = DataFrame(exps, :auto)
			CSV.write(string(exportpath, "Ωz-", Ωz, "_", basis_label, "-filter-trajectory.csv"), df)
		end

	end # write trajectories


end # loop over Ωz



"Write data --------------------------------------------------------------------"

"fidelities"
df = DataFrame([Ωzs, system_fid_avg, filter_fid_avg], ["Ωz", "system", "filter"])
CSV.write(string(exportpath, "fidelities-with-target-", targetstr, ".csv"), df)


"Parameter dictionary"

pars = Dict("N" => N,
			"Γ" => Γ,
			"Ωzs" => Ωzs,
			"φ" => φ,
			"η" => η,
			"dt" => dt,
			"tf" => tf,
			"td" => td,
			"Δt" => Δt
			)

CSV.write(string(exportpath, "parameters.csv"), pars)
