"""
		readout(dt::Timescale,
	 	    m::Function/QOp, 	-- measurement operator
			Γ::Function/Rate,	-- measurement rate
			η::Efficiency		-- quantum efficiency
			)

Returns function generating random (simulated) readout based on measurement of operator m
with rate Γ and quantum efficiency η. m and Γ may be either constant QOps or Functions (of time).
Note that m need not be Hermitian (observable); only the Hermitian part of the operator --
mo = (m + m')/2 -- will be used to generate the readout.

### Returns:
  - (t::Timescale, ρ) -> σ*randn() + real(expect(mo, ρ))  -- normally distributed random variable with mean
  														real(expect(ρ, mo)) and variance σ^2 = τ/dt
"""
function readout(dt::Timescale, m::Function, Γ::Function, η::Efficiency)
	(t, ρ) -> let mo = (m(t) .+ m(t)')/2;
					τ = 1/(2Γ(t) * η)
					σ = sqrt(τ/dt);
					σ*randn() + real(expect(mo, ρ)) end
end
function readout(dt::Timescale, m::Function, Γ::Rate, η::Efficiency)
	τ = 1/(2Γ * η)
	σ = sqrt(τ/dt)
	(t, ρ) -> let mo = (m(t) .+ m(t)')/2;
					σ*randn() + real(expect(mo, ρ)) end
end
function readout(dt::Timescale, m::QOp, Γ::Function, η::Efficiency)
	mo = (m .+ m')/2
	(t, ρ) -> let τ = 1/(2Γ(t) * η)
					σ = sqrt(τ/dt);
					σ*randn() + real(expect(mo, ρ)) end
end
function readout(dt::Timescale, m::QOp, Γ::Rate, η::Efficiency)
	τ = 1/(2Γ * η)
	mo = (m .+ m')/2
	σ = sqrt(τ/dt)
	(t, ρ) -> σ*randn() + real(expect(mo, ρ))
end
