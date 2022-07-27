"""
	fidelity(ρ::State, σ::State)

Gives the fidelity of two quantum states, ρ and σ, using the standard formula Tr(√(√σ ρ √σ)).
Uses the built-in QuantumOpticsBase.fidelity, but catches cases where ρ or σ may be a Ket,
and makes the appropriate conversion.

### Returns:
  - Tr(√(√σ ρ √σ)) :: Float64
"""
fidelity(ρ::Ket, σ::Ket) = real(ρ' * σ)^2
fidelity(ρ::Ket, σ::QOp) = real(fidelity(dm(ρ), σ))
fidelity(ρ::QOp, σ::Ket) = real(fidelity(ρ, dm(σ)))
fidelity(ρ::QOp, σ::QOp) = real(QuantumOpticsBase.fidelity(DenseOperator(ρ), DenseOperator(σ)))

"""
	expectations(sol::Solution, op::State)
	expectations(ρs::Vector{State}, op::State)

Gives a time series of expectation values of op for the density matrices sol.ρ or ρs.

### Returns :: Vector{Float64}
"""
function expectations(ρs::Vector, op::QOp)
	map(state -> real(expect(op, state)), ρs)
end
function expectations(ρs::Vector, op::Ket)
	map(state -> real(expect(dm(op), state)), ρs)
end
function expectations(sol::Solution, op)
	expectations(sol.ρ, op)
end

"Kronecker delta function"
δ(i,j) = Int(i == j);


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
