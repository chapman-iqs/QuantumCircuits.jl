module QuantumCircuits

using Reexport
@reexport using QuantumOpticsBase
using Distributions, Distributed, ProgressMeter

#################################
# Abstract types and type aliases
abstract type QObj <: Any end

# Note implementation is light-weight, using mostly type aliases
# for the features in Base that do all the real work
const Timescale = const Rate = const Efficiency = Real
const Record = Vector{Float64}	# Measurement record [ra_1, ra_2, ..., ra_t, ... , ra_f] corresponding to measurement
								# outcomes of a single observable a **over time**
const Readout = Vector{Float64}	# Readout [ra_t, rb_t, rc_t, ... ] corresponding to measurement outcomes of measured
								# observables {a, b, c}, etc. at a **single time** t
								# Note that the transpose of a Vector of Records is a Vector of Readouts, and vice versa
const QOp = AbstractOperator
const State = Union{Ket, QOp}

# types for testing function argument types using `applicable`
const tt = 0.0 # Timescale
const rr = [0.0] # Readout
const ρρ = identityoperator(SpinBasis(1//2)) # State



# holds solutions to bayesian or rouchon (rouchon not implemented)
struct Solution
	t::Vector{Timescale}
	ρ::Vector{State}
	r::Vector{Record}
end

include("utils.jl")



"""
	update(a::QOp, ρ::State)

Update state ρ using operator a.

### Returns:
  - ρ::Ket : a * ρ
  - ρ::QOp : a * ρ * a'
"""

update(a::QOp, ψ::Ket) = a * ψ
update(a::QOp, ρ::QOp) = a * ρ * a'


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


"""
	ensemble(solve, T, ρ0, H, J, C; <keyword arguments>)

solve :: Function, 							-- solver (bayesian or rouchon; may not yet work for rouchon)
(t0,tf) :: Tuple							-- time interval of simulation
ρ0 :: State 								-- initial state
H :: Function/QOp,							-- Hamiltonian, that may be a function of time H(t)
J :: Vector{Tuple{Qop, Rate}},				-- vector of tuples of Lindblad operators and their rates
C :: Vector{Tuple{Qop, Rate, Efficiency}}; 	-- vector of tuples of measurement operators and their rates

Keyword Arguments:
dt :: Timescale = 1e-4,							-- simulation time step
records = Vector{Record}[],					-- vector of vector of input measurement records, s.t. records[i] gives the records
												for the ith trajectory, and records[i][m]::Record gives the record corresponding
												to measurement of the mth measurement operator
N :: Int64 = 10								-- number of trajectories


"""

function ensemble(solve, (t0, tf), ρ, H::Union{Function, QOp}, J, C; dt=1e-4, records=Vector{Record}[], N=10, batch_size=10, ops=[], kwargs...)

    solutions = @showprogress pmap(i -> begin
		rs = length(records) > 0 ? records[i] : Record[]
		return solve((t0, tf), ρ, H, J, C; dt=dt, records=rs, kwargs...)

    end, 1:N; batch_size=batch_size)

	return length(ops) == 0 ?
					solutions :
					map(op -> [expectations(sol, op) for sol in solutions], ops)
end
"Alternatively, ensemble can take a vector of Hamiltonians. It will run bayesian once for each Hamiltonian."
function ensemble(solve, (t0, tf), ρ, Hlist::Vector, J, C; dt=1e-4, records=Vector{Record}[], batch_size=10, ops=[], kwargs...)

    solutions = @showprogress pmap(H -> begin
		rs = length(records) > 0 ? records[i] : Record[]
		return solve((t0, tf), ρ, H, J, C; dt=dt, records=rs, kwargs...)

    end, Hlist; batch_size=batch_size)

	return length(ops) == 0 ?
					solutions :
					map(op -> [expectations(sol, op) for sol in solutions], ops)
end
function ensemble(solve, (t0, tf), ρ, (Hs, Hf), J, C; dt=1e-4, N=10, batch_size=10, ops=[], kwargs...)

    solutions = @showprogress pmap(i -> begin
		return solve((t0, tf), ρ, (Hs, Hf), J, C; dt=dt, kwargs...)
    end, 1:N; batch_size=batch_size)
	# returns a vector of tuples: Vector{Tuple{Solution, Solution}} corresponding to system and filter solutions

	if length(ops) == 0
		return solutions
	else
		sys_exps = map(op -> [expectations(sys, op) for (sys, _) in solutions], ops)
		fil_exps = map(op -> [expectations(fil, op) for (_, fil) in solutions], ops)
		return (sys_exps, fil_exps)
	end
end
function ensemble(solve, (t0, tf), ρ, Htuples::Vector{Tuple}, J, C; dt=1e-4, N=10, batch_size=10, ops=[], kwargs...)

    solutions = @showprogress pmap(Htuple -> begin
		return solve((t0, tf), ρ, Htuple, J, C; dt=dt, kwargs...)
    end, Htuples; batch_size=batch_size)
	# returns a vector of tuples: Vector{Tuple{Solution, Solution}} corresponding to system and filter solutions

	if length(ops) == 0
		return solutions
	else
		sys_exps = map(op -> [expectations(sys, op) for (sys, _) in solutions], ops)
		fil_exps = map(op -> [expectations(fil, op) for (_, fil) in solutions], ops)
		return (sys_exps, fil_exps)
	end
end



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
			ftrajectory(m, h, l, ts, ρ; fn=fn, r0=r0, dt=dt, td=td)
		else
			error(string("Feedback using a predetermined measurement record is not currently ",
							"supported. Please rerun the simulation without an input record."))
		end
	end
end # bayesian
function bayesian((t0, tf), ρ, (Hs0, Hf0), J, C; fn=ρ->ρ, dt=1e-4, td::Timescale = 0.0, readout_integration_bins=1)



	if !(typeof(Hs0) <: Function) || !(typeof(Hf0) <: Function) || applicable(Hs0, tt) || applicable(Hf0, tt)
		str = string("Non-feedback system-filter Hamiltonians are not yet supported. Please input Hamiltonians ",
						"as a tuple of functions (Hs(args1...), Hf(args2...)), where args1, args2 = (t::Timescale, ρ::State),",
						"(t::Timescale, r::Record), or (t::Timescale, ρ::State, r::Record), and Hs and Hf ",
						"are the system and the filter Hamiltonians, respectively. This feature should ",
						"be added in the future...."
						)
		error(str)

	end

	Hs = convertham(Hs0)
	Hf = convertham(Hf0)

	ts = range(t0, tf, step=dt)
	r0::Readout = [0.0 for i in 1:length(C)] # initial, empty Readout object

	# define superoperator functions for state update
	ms = meas(dt, C; sample=true) 	# "system" measurement superoperator (generates a readout based on system state)
	mf = meas(dt, C; sample=false) 	# "filter" measurement superoperator (ALSO generates a readout based on system state
									# -- same superoperator, but needs to be able to take readout as argument)
	hs = ham(dt, Hs) # unitary superoperator for experimentally "known" Hamiltonian (applied to both "system" and "filter")
	hf = ham(dt, Hf) # unitary superoperator for experimentally "unknown" Hamiltonian (applied to only "system")
	l = lind(dt, J) # decay superoperator

	if length(J) > 0 && typeof(ρ) <: Ket
		ρ = dm(ρ)
	end

	ftrajectory((ms, mf), (hs, hf), l, ts, ρ; fn=fn, r0=r0, dt=dt, td=td, readout_integration_bins=readout_integration_bins)
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
function trajectory(m::Function, h::Function, l::Function, ts, ρ; fn::Function=ρ->ρ, dt=1e-4, r0=[0.0])

	# initialize
	ρ0 = ρ
	fnρs = [fn(ρ0)]
	readouts = [r0]


	for t in ts[2:end]
		ρ1 = h(t, ρ) # unitary evolution
		ρ2 = l(t, ρ1) # lindblad evolution
		ρ, ro = m(t, ρ2) # measurement evolution

		push!(fnρs, fn(ρ))
		push!(readouts, ro)

	end

	records = trans(readouts)

	return Solution(ts, fnρs, records)
end
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

	return Solution(ts, fnρs, records)
end
# currently, feeding in a pre-determined measurement record is not supported for feedback,
# since it's not clear when this would be appropriate.
# NOTE: chose to apply unitary, then lindblad, then measurement, so that feedback is applied directly after measurement obtained
# but should check if this is correct.
function ftrajectory(m::Function, h::Function, l::Function, ts, ρ; fn::Function=ρ->ρ, dt=1e-4, r0=[0.0], td=0.0)

	# initialize
	ρ0 = ρ
	ρs = [ρ0]
	fnρs = [fn(ρ0)]
	readouts = [r0]

	# if there is no time delay, delay = 1, and we use the record / state from most recent update in the
	# feedback Hamiltonian
	delay = Int64(round(td / dt))

	# before feedback delay buffer is filled, we must feed back on "fake" record or state
	for t in ts[2:(delay + 1)]

		# measurement evolution
		ρ1, ro = m(t, ρ)
		push!(readouts, ro)
		push!(ρs, ρ1)
		push!(fnρs, fn(ρ1))

		# unitary evolution
		ρ2 = h(t, ρ1, ρ0, r0)
		# ρ2 = apply(h, t, ρ1, ρ0, r0)

		# lindblad evolution
		ρ = l(t, ρ2)




	end

	# after feedback delay buffer is filled, we can feed back on acquired readout / filtered state
	for i in (delay + 2):length(ts)
		t = ts[i]

		# measurement evolution
		ρ1, ro = m(t, ρ)
		push!(readouts, ro)
		push!(ρs, ρ1)
		push!(fnρs, fn(ρ1))

		# unitary evolution
		ρd = ρs[i - delay]
		rd = readouts[i - delay]
		ρ2 = h(t, ρ1, ρd, rd)

		# rd = readouts[i - delay]
		# ρd = ρs[i - delay]
		# ρ2 = apply(h, t, ρ1, ρd, rd)

		# lindblad evolution
		ρ = l(t, ρ2)




    end

	records = trans(readouts)

	return Solution(ts, fnρs, records)
end
function ftrajectory((ms, mf), (hs, hf), l::Function, ts, ρ0; fn::Function=ρ->ρ, dt=1e-4, r0=[0.0], td=0.0, readout_integration_bins=1)

	# initialize system
	ρs = ρ0
	ρss = [ρs]
	fnρs = [fn(ρs)]

	# initialize filter
	ρf = ρ0
	ρfs = [ρf]
	fnρf = [fn(ρ0)]

	# shared readout for system and filter
	readouts = [r0]

	# if there is no time delay, delay = 1, and we use the record / state from most recent update in the
	# feedback Hamiltonian
	delay = Int64(round(td / dt)) + 1


	# feed back on filter state; obtain record using system state
	# before feedback delay buffer is filled, we must feed back on "fake" record or state
	for t in ts[2:delay]

		### system evolution
		ρs1 = hs(t, ρs, ρ0, r0)# unitary
		ρs2 = l(t, ρs1) # lindblad
		ρs, ro = ms(t, ρs2) # measurement

		push!(ρss, ρs)
		push!(readouts, ro)
		push!(fnρs, fn(ρs))

		### filter evolution
		ρf1 = hf(t, ρf, ρ0, r0) # unitary
		ρf2 = l(t, ρf1) # lindblad
		ρf = mf(t, ρf2, ro) # measurement (based on filtered state expectation values, but with system readout)

		push!(ρfs, ρf)
		push!(fnρf, fn(ρf))



	end

	# after feedback delay buffer is filled, we can feed back on acquired readout / filtered state
	for i in (delay + 1):length(ts)
		t = ts[i]
		j = i - delay
		ρd = ρfs[j] # state used in feedback is the filtered state from time td ago
		startsmooth = max(1, j - (readout_integration_bins - 1)) # don't index readouts to smooth from an index less than 1
		rd = mean(readouts[startsmooth:j]) # mean(trans(readouts[startsmooth:j])) 	# readout for feedback is an average of the readout_integration_bins readouts
													# before the delayed time index j

		### system evolution
		ρs1 = hs(t, ρs, ρd, rd) # unitary
		ρs2 = l(t, ρs1) # lindblad
		ρs, ro = ms(t, ρs2) # measurement

		push!(ρss, ρs)
		push!(readouts, ro)
		push!(fnρs, fn(ρs))

		### filter evolution
		ρf1 = hf(t, ρf, ρd, rd) # unitary
		ρf2 = l(t, ρf1) # lindblad
		ρf = mf(t, ρf2, ro) # measurement (based on filtered state expectation values, but with system readout)

		push!(ρfs, ρf)
		push!(fnρf, fn(ρf))

    end

	records = trans(readouts)

	return (Solution(ts, fnρs, records), Solution(ts, fnρf, Record[]))
end



"""
   meas(dt::Timescale								-- simulation time step
		C::Array{Tuple{Qop, Rate, Efficiency}}; -- array of tuples of measurement operators and their rates
		sample::Boolean = true					-- if sample, generate simulated readout using random number generator
		)

### Returns:
	- sample=true : (t::Timescale, ρ::State) --> normalize(ρ1)::State, readout::Readout)
	(updated density matrix after measurement, generated readout)
  	- sample=false : (t::Timescale, ρ::State, r::Readout) --> normalize(ρ1)::State
	(updated density matrix after measurement)


"""
function meas(dt::Timescale, C; sample=true)

	# Assemble readout generating functions ros and Kraus operators gks
	ros = Function[]
	gks = Function[]

	for (m, Γ, η) in C
		if sample
			push!(ros, readout(dt, m, Γ, η)) end
		push!(gks, gausskraus(dt,  m, Γ, η))
	end

	if sample
		# Increment that samples each readout, applies all Kraus operators
		(t, ρ) -> begin
			# evaluate readout generating functions at time t for state ρ
			readout = map(ro -> ro(t, ρ), ros)

			# evaluate each gausskraus operator at time t for the acquired readout
			gs = [g(t, r) for (g, r) in zip(gks, readout)]

			# apply each element of gs via `update` function (sandwich) to ρ
			ρ1 = foldr(update, gs; init=ρ);
			return (normalize(ρ1), readout)
		end

	else
		# same as above, but takes readout as argument
		(t, ρ, readout) -> begin

			gs = [g(t, r) for (g, r) in zip(gks, readout)]
			ρ1 = foldr(update, gs; init=ρ);
			return normalize(ρ1)
		end

	end
end


# Hamiltonian propagation
"""
	ham(dt, H)

Return increment function for Hamiltonian evolution generated
by `H` over a time step `dt`.

Uses an exact (dense) matrix exponential, assuming no time-dependence.

### Returns:
  - typeof(H) <: QOp 		: (t, ρ) -> update(u, ρ)			-- increment function for Hamiltonian-evolved density matrix
  - typeof(H) <: Function 	: (t, ρ) -> ham(dt, H(t))(t, ρ)		-- increment function for Hamiltonian-evolved density matrix,
   																	evaluated at time t

"""
function ham(dt, H::QOp)
	u::QOp = exp( -im * dt * DenseOperator(H))
	(t, ρ) -> update(u, ρ)
end
function ham(dt, H::Function)
	# evaluate H for specific inputs (using test constants tt::Timescale, rr::Record, ρρ::State)
	if applicable(H, tt)
		(t::Timescale, ρ::State) -> ham(dt, H(t))(t, ρ)
	else
		(t::Timescale, ρ::State, ρd::State, r::Readout) -> ham(dt, H(t, ρd, r))(t, ρ)
	end

end


# Jump-nojump Lindblad propagator
"""
	lind(dt[,
	 	  J::Vector{Tuple{QOp, Rate}}];
		  clist=QOp[],
		  flist=Function[]
		  )

Returns increment function over a time step dt for Lindblad dissipative
evolution generated by dissipative operators √γ * j, corresponding to
pairs (j::QOp, γ::Rate) in J.

If J is input as argument to lind, the function sorts the pairs into two
lists: a list `clist` of constant dissipative operators and a list `flist`
of time-dependent functions returning dissipative operators.

### Returns:
  - when J input : (t, ρ) -> lind(dt; clist=clist, flist=flist)
  - else : (t, ρ(t)::QOp) -> L(ρ), where L denotes the Lindblad superoperator
"""
function lind(dt, J)
	# preprocess lindblad arguments
	clist = []
	flist = []

	for (j, γ) in J
		if typeof(γ) <: Function
			push!(flist, (j, γ))
		else
			push!(clist, (j, γ))
		end
	end
	(t, ρ) -> lind(dt; clist=clist, flist=flist)(t, ρ)
end
function lind(dt; clist=QOp[], flist=Function[])
	ns = Function[]
	ds = Function[]
	# Construct operations for constant operators
	if isempty(clist)
		push!(ns, (t, ρ) -> ρ)
	else
		Id = identityoperator(first(clist[1]).basis_l)
		op = DenseOperator(Id - dt * (mapreduce(+, clist) do (j, γ) γ * j' * j end))
		n::Operator = SparseOperator(op.basis_l, op.basis_r, sqrt(op.data))
		push!(ns, (t, ρ) -> n * ρ * n)
		push!(ds, (t, ρ) -> dt * (mapreduce(+, clist) do (j, γ) γ * j * ρ * j' end))
	end
	# Construct operations for time-dependent operators
	if isempty(flist)
		push!(ns, (t, ρ) -> ρ)
	else
		function nf(t)
			Id = identityoperator(first(flist)[1].basis_l)
			op = DenseOperator(Id - dt * (mapreduce(+, flist) do (j, γ) γ(t) * j' * j end))
			return SparseOperator(op.basis_l, op.basis_r, sqrt(op.data))
		end
		push!(ns, (t, ρ) -> nf(t) * ρ * nf(t))
		push!(ds, (t, ρ) -> dt * (mapreduce(+, flist) do (j, γ) γ(t) * j * ρ * j' end))
	end
	push!(ds, (t, ρ) -> last(ns)(t, first(ns)(t, ρ)))
	(t, ρ) -> mapreduce(f -> f(t, ρ), +, ds)
end



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




"""
    rouchon(T, ρ0, H, J, C; <keyword arguments>)

Arguments:

T :: tuple (ti,tf)
ρ0 :: initial density matrix
H :: (time-dependent) Hamiltonian
J :: array of tuples (j, Γ) representing decoherence channel j at rate Γ
C :: array of tuples (c, τ, η) representing measurement of c with timescale τ and collection efficiency η

Keyword Arguments:

dt :: time step; default dt=1e-4
r :: record; default r=[], i.e. simulation generates record by randomly sampling distribution.
            record should be input in the shape [r_1,...,r_Nc] given
            length(C)=Nc collapse operators. Records should have shape
            r_m[i] indexing the mth trajectory at time ts[i], where
			ts = range(first(T), last(T), step=dt)
fn : ρ → Any :: eval function (e.g. to return expectation values instead of density matrices)

Returns: (ts, ρs, r)

ts :: list of simulation times
ρs :: fn(ρ) at each simulation time
r :: input OR simulated record, depending on value of keyword argument r
"""

function rouchon((t0, tf), ρ, H0, J0, Ctups; fn=ρ->ρ, dt=1e-4, r=[])
    ts = range(t0, tf, step=dt)
	ρ = (typeof(ρ) <: Ket) ? dm(ρ) : ρ # eventually, include ket simulation functionality for rouchon
    Id = identityoperator(ρ.basis_l)
    Nt = length(ts)
    Nj = length(J0)
    Nc = length(Ctups)

    ρs = Any[]

    H = length(methods(H0)) > 0 ? H0 : t -> H0
    # J = map(j -> length(methods(j)) > 0 ? j : t -> j, J0)
	J = []
	for (j, Γ) in J0
		push!(J, length(methods(j)) > 0 ? √Γ * j : (t -> √Γ * j) ) end
	C = []
	for (i, (c, τm, η)) in enumerate(Ctups)
		m = c * sqrt(η / 2τm)
		push!(C, length(methods(m)) > 0 ? m : (t -> m))
		# push!(C, length(methods(c)) > 0 ? (c,τm,η) : (t -> c,τm,η))
	end

	# sample / prepare dy array
	sim = (length(r) == 0)
	dy = []
	dW = []
	dist = Normal(0, √dt)

	for (i, (c, τm, η)) in enumerate(Ctups)
		if sim # randomly sample noise time series for EACH stochastic collapse operator C
			push!(dW, rand(dist, Nt))
			push!(dy, zeros(Nt))
		else
			push!(dy, r[i] .* dt ./ sqrt(τm))	end
	end


    for (n, t) in enumerate(ts)
        M = Id - im * H(t) * dt

        # iterate over deterministic collapse operators J
        D = 0Id
        for j in 1:Nj
            M += -0.5J[j](t)' * J[j](t) * dt
            D += J[j](t) * ρ * J[j](t)' * dt
        end

        # initialize dy value, if simulation
        if sim
            for (i, (c, τm, η)) in enumerate(Ctups)
                dy[i][n] = (real(tr(C[i](t) * ρ + ρ * C[i](t)') * dt) + dW[i][n])
            end
        end

        # iterate over stochastic collapse operators C
        for c in 1:Nc
            M += C[c](t) * dy[c][n]
            # D += -C[c](t)*ρ*C[c](t)'*dt

            # nested sum
            for s in 1:Nc
                M += 0.5C[c](t) * C[s](t) * (dy[c][n] * dy[s][n] - δ(c,s) * dt)
            end
        end

        # update ρ according to Rouchon method
        ρ = M*ρ*M' + D
        ρ = ρ / tr(ρ)
        push!(ρs, fn(ρ))
    end

	# revert to r (normalized) version of record
	if sim
		for (i, (c, τm, η)) in enumerate(Ctups)
			push!(r, dy[i] .* sqrt(τm) ./ dt)
		end
	end


    return Solution(ts, ρs, r)
end # rouchon


"""
	Applies Function h to State ρ with correct arguments.
"""
function apply(h::Function, t::Timescale, ρ::State, ρd::State, rd::Record)
	arglist = [(t, ρ, rd), (t, ρ, ρd), (t, ρ, ρd, rd)]
	index = findfirst(map(args -> applicable(h, args...), arglist))
	args = arglist[index]

	return h(args...)
end

"""
	Converts hamiltonian function to have correct argument order and type.
"""
function convertham(H::Function)

	if applicable(H, tt) # non-feedback case
		return H
	else # feedback case
		return (t::Timescale, ρ::State, r::Record) -> let
			arglist = [(t, ρ), (ρ, t), (t, r), (r, t), (t, ρ, r), (t, r, ρ), (r, t, ρ), (r, ρ, t), (ρ, r, t)]
			index = findfirst(map(args -> applicable(H, args...), arglist))
			args = 	try
						arglist[index]
					catch e
						error("Improper argument types input to Hamiltonian.")
					end

		 	H(args...)
		end
	end
end # convertham


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

export bayesian, rouchon, ensemble
export δ, fidelity, coarse_grain, subselect, expectations
export Timescale, Rate, Efficiency, Record, Readout, QOp, State, Solution


end # module
