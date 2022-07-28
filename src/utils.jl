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


purity(x::Real, y::Real, z::Real) = 0.5(1 + x^2 + y^2 + z^2)
purity(ψ::Ket) = 1.0
purity(ρ::QOp) = real(tr(ρ * ρ))

purities(x::Vector, y::Vector, z::Vector) =  0.5(1 .+ x.^2 .+ y.^2 + z.^2)

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


"""
    coarse_grain(fine; <keyword arguments>)
Argument: fine :: time-series to be coarse-grained
Keyword Argument: n :: (default n=2), number of elements over which to smooth input time-series
Returns: coarse :: smoothed time-series of the same length as fine, but having taken a
					moving-average over n elements
"""
function coarse_grain(fine::Array; n=2)
    coarse = []
	for i in 1:(n-1)
		push!(coarse, mean(fine[1:i])) end

	for i in n:length(fine)
		push!(coarse, mean(fine[i-(n-1):i])) end

    coarse
end

function subselect(a=[]; n=2)
    a[filter(x -> x%n==0, eachindex(a))]
end



"""
trans(mat)

### Returns:
	- Vector{<:Vector{<:Real}}: transpose of mat"
"""
function trans(mat::Vector{<:Vector{<:Real}})
	nrows = length(mat)
	ncols = length(mat[1])
	return [[mat[i][j] for i in 1:nrows] for j in 1:ncols]
end

"Kronecker delta function"
δ(i,j) = Int(i == j);
