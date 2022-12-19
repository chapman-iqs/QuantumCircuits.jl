
"""
	gausskraus(dt::Timescale,
	 	    m::Function/QOp, 	-- measurement operator
			Γ::Function/Rate,	-- measurement rate
			η::Efficiency		-- quantum efficiency
			)

Returns function generating measurement Kraus operator based on readout r,
for measurement of operator m with rate Γ and quantum efficiency η. m and Γ
may be either constant QOps or Functions (of time), and m need not be
Hermitian (observable).

### Returns:
  - (t::Timescale, r) -> SparseOperator(exp(dense((r*v)*m - v*mo2)))

"""
function gausskraus(dt::Timescale, m::QOp, Γ::Rate, η::Efficiency)
	τ = 1/(2Γ * η)
	mo = (m .+ m') / 2
	mo2 = mo^2 / 2
	v = dt/(2τ)
	(t::Timescale, r) -> SparseOperator(exp(dense((r*v)*m - v*mo2)))
end
function gausskraus(dt::Timescale, m::Function, Γ::Rate, η::Efficiency)
	τ = 1/(2Γ * η)
	v = dt/(2τ)
	(t::Timescale, r) -> let mo = (m(t) .+ m(t)') / 2
						mo2 = mo^2 / 2
						SparseOperator(exp(dense((r*v) * m(t) - v*mo2))) end
end
function gausskraus(dt::Timescale, m::QOp, Γ::Function, η::Efficiency)
	mo = (m .+ m') / 2
	mo2 = mo^2 / 2
	(t::Timescale, r) -> let τ = 1/(2Γ(t) * η)
						v = dt/(2τ)
						SparseOperator(exp(dense((r*v)*m - v*mo2))) end
end
function gausskraus(dt::Timescale, m::Function, Γ::Function, η::Efficiency)
	(t::Timescale, r) -> let mo = (m(t) .+ m(t)') / 2
						mo2 = mo^2 / 2
						τ = 1/(2Γ(t) * η)
						v = dt/(2τ)
						SparseOperator(exp(dense((r*v) * m(t) - v*mo2))) end
end
