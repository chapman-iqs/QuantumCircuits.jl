"""
	update(a::QOp, ρ::State)

Update state ρ using operator a.

### Returns:
  - ρ::Ket : a * ρ
  - ρ::QOp : a * ρ * a'
"""

update(a::QOp, ψ::Ket) = a * ψ
update(a::QOp, ρ::QOp) = a * ρ * a'
update(a::Function, ρ::State) = f(ρ) # try this for general state update


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
	elseif applicable(H, tt, ρρ)
		(t::Timescale, ρ::State, ρd::State, r::Readout) -> ham(dt, H(t, ρd))(t, ρ)
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
	fρlist = []

	for (j, γ) in J
		if typeof(γ) <: Function
			if applicable(γ, tt)
				push!(flist, (j, γ))
			elseif applicable(γ, tt, ρρ)
				push!(fρlist, (j, γ))
			else
				error("Please enter γ(t::Timescale) or γ(t::Timescale, ρ::State)")
			end
		else
			push!(clist, (j, γ))
		end
	end
	(t, ρ) -> lind(dt; clist=clist, flist=flist, fρlist=fρlist)(t, ρ)
end
function lind(dt; clist=QOp[], flist=Function[], fρlist=Function[])
	ns = Function[]
	ds = Function[]
	# Construct operations for constant operators
	if isempty(clist)
		push!(ns, (t, ρ) -> ρ)
	else
		n = nullop(dt, clist...)
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
		# nf = nullop(dt, flist...)
		push!(ns, (t, ρ) -> nf(t) * ρ * nf(t))
		push!(ds, (t, ρ) -> dt * (mapreduce(+, flist) do (j, γ) γ(t) * j * ρ * j' end))
	end

	# Construct operations for time- and state-dependent operators
	if isempty(fρlist)
		push!(ns, (t, ρ) -> ρ)
	else
		function nff(t, ρ)
			Id = identityoperator(first(fρlist)[1].basis_l)
			op = DenseOperator(Id - dt * (mapreduce(+, fρlist) do (j, γ) γ(t, ρ) * j' * j end))
			return SparseOperator(op.basis_l, op.basis_r, sqrt(op.data))
		end
		# nff = nullop(dt, fρlist...)
		push!(ns, (t, ρ) -> nff(t, ρ) * ρ * nff(t, ρ))
		push!(ds, (t, ρ) -> dt * (mapreduce(+, fρlist) do (j, γ) γ(t, ρ) * j * ρ * j' end))
	end

	ntotal = (t, ρ) -> let ρ1 = ρ
							for n in ns
								ρ1 = n(t, ρ1)
							end
							return ρ1
						end

	push!(ds, ntotal)

	# push!(ds, (t, ρ) -> last(ns)(t, first(ns)(t, ρ))) # combines into one null map
	(t, ρ) -> mapreduce(f -> f(t, ρ), +, ds)
end



function nullop(dt::Timescale, (j, γ)::Tuple{QOp, Rate}, tuples...) 
	# concatenate all the tuples together; first is distinguished for dispatch
	tuples = [(j, γ), tuples...] 

	# get identity operator from left basis element of j
	Id = identityoperator(j.basis_l)

	# get the operator as √(I - dt * γ1 * j1' * j1 - dt * γ2 * j2' * j2 - ...)
	# square root is taken on op.data at the end
	op = DenseOperator(Id - dt * (mapreduce(+, tuples) do (j, γ) γ * j' * j end))
	return SparseOperator(op.basis_l, op.basis_r, sqrt(op.data))
end

function nullop(dt::Timescale, (j, γ)::Tuple{QOp, Function}, tuples...) 
	# concatenate all the tuples together; first is distinguished for dispatch
	tuples = [(j, γ), tuples...] 

	# get identity operator from left basis element of j
	Id = identityoperator(j.basis_l)

	# get the operator as √(I - dt * γ1 * j1' * j1 - dt * γ2 * j2' * j2 - ...)
	# square root is taken on op.data at the end 
	# args can be (t) or (t, ρ)
	op(args) = DenseOperator(Id - dt * (mapreduce(+, tuples) do (j, γ) γ(args...) * j' * j end))
	return args -> SparseOperator(op(args).basis_l, op(args).basis_r, sqrt(op(args).data))
end
