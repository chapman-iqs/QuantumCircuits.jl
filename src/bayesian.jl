"""
	bayesian(
		(t0, tf)::Tuple{Timescale, Timescale}		-- time interval of simulation
		ρ::State, 									-- initial state
		H::Function/QOp,							-- Hamiltonian, that may be a function of time H(t)
		J::Vector{Tuple{Qop, Rate}},				-- vector of tuples of Lindblad operators and their rates
		C::Vector{Tuple{Qop, Rate, Efficiency}}; 	-- vector of tuples of measurement operators and their rates
		fn::Function = ρ->ρ, 						-- function of ρ to output
		dt::Timescale = 1e-4,						-- simulation time step
		records = Record[],							-- vector of input measurement records
		td::Timescale = 0.0)						-- time delay for feedback

### Returns:
  - Solution object, a tuple (t::Vector{Timescale}, 	-- time series vector
							  ρ::Vector{State},			-- vector of states at each time in t
							  records::Vector{Record})	-- vector of output measurement records
"""
function bayesian((t0, tf), ρ, H0::Union{Function, QOp}, J, C; fn=ρ->ρ, dt=1e-4, records=Record[], td::Timescale = 0.0)

	ts = range(t0, tf, step=dt)
	r0::Readout = [0.0 for i in 1:length(C)] # initial, empty Readout object
	sample = (length(records) == 0)

	H = (typeof(H0) <: Function) ? convertham(H0) : H0

	# define superoperator functions for state update
	m = meas(dt, C; sample=sample) 	# measurement superoperator.
									# If sample=true, meas will generate a random record r at time t and returns a function
									# of (t, ρ), using the randomly generated r in the state update.
									# If sample=false, meas returns a function of (t, ρ, r) and uses the input r in the state update.
	h = ham(dt, H) # unitary superoperator
	l = lind(dt, J) # decay superoperator

	if length(J) > 0 && typeof(ρ) <: Ket
		ρ = dm(ρ)
	end

	# return trajectory
	if applicable(h, tt, ρρ) # no feedback
		if sample
			trajectory(m, h, l, ts, ρ; fn=fn, r0=r0, dt=dt)
		else
			trajectory(m, h, l, ts, ρ, records; fn=fn, dt=dt)
		end
	else # feedback
		if sample
            trajectory(SimpleTimeDelay(), m, h, l, ts, ρ; fn=fn, r0=r0, dt=dt, td=td)
		else
			error(string("Feedback using a predetermined measurement record is not currently ",
							"supported. Please rerun the simulation without an input record."))
		end
	end
end # bayesian


"""

	trajectory(	m::Function,				-- measurement superoperator
				h::Function,				-- unitary superoperator
				l::Function,				-- lindblad superoperator
				ts,							-- simulation times
				ρ[,							-- initial density matrix
				records::Vector{Record};	-- OPTIONAL: measurement record (if not included, will simulate record)
				fn::Function=ρ->ρ, 			-- function of ρ to output
				dt=1e-4 [, 					-- simulation time step
				r0=[0.0]]					-- first Readout, to signal number of measurement
			)

Uses the "jump no-jump" method to efficiently approximate the exact
Lindblad propagator as a composition of Hamiltonian evolution, jumps,
and no-jump informational backaction. Assumes small dt.
[Physical Review A **92**, 052306 (2015)]


### Returns:
	- Solution(ts::Vector{Timescale}, ρs::Vector{State})


NOTE: temporarily took away @inline macro ...

"""
# non-feedback trajectory generating its own measurement records
function trajectory(m::Function, h::Function, l::Function, ts, ρ; fn::Function=ρ->ρ, dt=1e-4, r0=[0.0])

	# initialize
	ρ0 = ρ
	fnρs = [fn(ρ0)]
	readouts = Readouts([r0])


	for t in ts[2:end]
		ρ1 = h(t, ρ) # unitary evolution
		ρ2 = l(t, ρ1) # lindblad evolution
		ρ, ro = m(t, ρ2) # measurement evolution

		push!(fnρs, fn(ρ))
		push!(readouts, ro)

	end

	records = trans(readouts)

	return Solution(collect(ts), fnρs, records)
end

# non-feedback trajectory generated from input measurement records
function trajectory(m::Function, h::Function, l::Function, ts, ρ, records::Vector{Record}; fn::Function=ρ->ρ, dt=1e-4)

	# initialize
	ρ0 = ρ
	fnρs = [fn(ρ0)]
	readouts = trans(records)

	for (t, ro) in zip(ts[2:end], readouts[2:end])
		 ρ1 = h(t, ρ) # unitary evolution
		 ρ2 = l(t, ρ1) 	# lindblad evolution ... is it appropriate to apply before measurement? shouldn't matter right?
		 				# did this for consistency with ftrajectory
		 ρ = m(t, ρ2, ro) # measurement evolution

		 push!(fnρs, fn(ρ))
	end

	return Solution(collect(ts), fnρs, records)
end
