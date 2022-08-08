# Solution holds solutions to bayesian or rouchon (rouchon not implemented)#################################
# Abstract types and type aliases
abstract type QObj <: Any end

# Note implementation is light-weight, using mostly type aliases
# for the features in Base that do all the real work
const Timescale = Float64
const Times = Vector{Timescale}
const Rate = Float64
const Efficiency = Float64
const Record = Vector{Float64}	# Measurement record [ra_1, ra_2, ..., ra_t, ... , ra_f] corresponding to measurement
								# outcomes of a single observable a **over time**
const Records = Vector{Record}
const Readout = Vector{Float64}	# Readout [ra_t, rb_t, rc_t, ... ] corresponding to measurement outcomes of measured
								# observables {a, b, c}, etc. at a **single time** t
								# Note that the transpose of a Vector of Records is a Vector of Readouts, and vice versa
const Readouts = Vector{Readout}
const QOp = AbstractOperator
const State = Union{Ket, QOp}
const States = Vector{State}

# types for testing function argument types using `applicable`
const tt = Timescale(0.0)
const rr = Readout([0.0])
const RR = Record([0.0])
const ψψ = Ket(SpinBasis(1//2))
const ρρ = identityoperator(SpinBasis(1//2)) # State


const BasisName = Union{Symbol, Vector{Symbol}}

Base.@kwdef mutable struct Solution

	t::Times
	ρ
	r::Records

	exps = missing
	basis = missing
	basisname::BasisName = missing

	fs::Union{Float64, Missing} = missing
	functions::Union{Function, Missing} = missing
	functionnames::Union{Symbol, Vector{Symbol}, Missing} = missing
end
Solution(t::Times, ρ) = Solution(t, ρ, Records(), [], [], Symbol())
Solution(t::Times, ρ, r::Records) = Solution(t, ρ, r, [], [], Symbol())
Solution(t::Times, ρ, r::Records, basis, bn::BasisName = Symbol()) =
	let exps = [expectations(ρ, op) for op in basis]
	return Solution(t, ρ, r, exps, basis, bn)
end
Solution(t::Times, ρ, basis, bn::BasisName) = Solution(t, ρ, r, basis, bn)
Solution(sol::Solution, basis, bn::BasisName = Symbol()) = Solution(sol.t, sol.ρ, sol.r, basis, bn)


"""
	mutable struct Ensemble(
				t::Times, 	-- times of simulation
				sols::Vector{Solution}, -- vector of `QuantumCircuits.Solution`s containing density matrices and measurement records
				exps,					-- expectation values of basis operators
				basis::Symbol) 			-- symbol representing basis of expectation values

Wrapper for `QuantumCircuits.Solution`s and their expectation values in a given basis.
Dispatches on arguments to provide a simple wrapper for a Vector of `Solution`s, or
if a basis is provided, calculates expectation values in that basis.

###
Ensemble(sols::Vector{Solution}) 				-- returns Ensemble(sols[1].t, sols, (), Symbol())
Ensemble(sols::Vector{Solution}, basis::Symbol) -- returns Ensemble(sols[1].t, sols, exps, basis)
Ensemble(ens::Ensemble, basis::Symbol)			-- returns Ensemble(sols[1].t, sols, exps, basis)

"""
mutable struct Ensemble
	t::Times
	sols::Vector{Solution}
	exps
	basis
	basisname::Union{Symbol, Vector{Symbol}}
end
Ensemble(sols::Vector{Solution}) = Ensemble(sols[1].t, sols, [], [], Symbol())
Ensemble(sols::Vector{Solution}, basis, basisname::Union{Symbol, Vector{Symbol}}) =
	let exps = [map(sol -> expectations(sol, op), sols) for op in basis]
	return Ensemble(sols[1].t, sols, exps, basis, basisname)
end
Ensemble(sols::Vector{Solution}, basis) = Ensemble(sols, basis, Symbol())
Ensemble(ens::Ensemble, basis) = Ensemble(ens.sols, basis)
Ensemble(ens::Ensemble, basis, basisname::Union{Symbol, Vector{Symbol}}) = Ensemble(ens.sols, basis, basisname)
