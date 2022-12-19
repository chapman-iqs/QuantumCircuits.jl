"""
Old utilities that are not currently in use by bayesian or rouchon solvers, but that I'm
not quite ready to delete yet.
"""

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

# function nullop(dt, J...)

# 	# get identity operator from left basis element of the first j
# 	(j1, γ1) = J[1]
# 	Id = identityoperator(j1.basis_l)

# 	# get the operator as √(I - dt * γ1 * j1' * j1 - dt * γ2 * j2' * j2 - ...)
# 	# square root is taken on op.data at the end
# 	op = DenseOperator(Id - dt * (mapreduce(+, J) do (j, γ) γ * j' * j end))
# 	return SparseOperator(op.basis_l, op.basis_r, sqrt(op.data))
# end




# function nullop(dt, J...)

# 	# get identity operator from left basis element of the first j
# 	(j1, γ1) = J[1]
# 	Id = identityoperator(j1.basis_l)

# 	# get the operator as √(I - dt * γ1 * j1' * j1 - dt * γ2 * j2' * j2 - ...)
# 	# square root is taken on op.data at the end
# 	op = DenseOperator(Id - dt * (mapreduce(+, J) do (j, γ) γ * j' * j end))
# 	return SparseOperator(op.basis_l, op.basis_r, sqrt(op.data))
# end
"""
function lind(dt, J)

	# if J is empty, we just return the identity map
	if isempty(J)
		return (t, ρ) -> ρ
	end

	# otherwise, create a map evaluating the J's at time t
	Jconst(t, ρ) = 	map(J) do (j, γ)
						if applicable(γ, tt, ρρ)
							(j, γ(t, ρ))
						elseif applicable(γ, tt)
							(j, γ(t))
						else
							(j, γ)
						end

					end

	# create a time-discrete CP map emulating Lindblad evolution
	return (t, ρ) -> let

		# extract the J applicable at this time
		Jt = Jconst(t, ρ)

		# get the corresponding null operator
		n = nullop(dt, Jt...)
		nullmap = ρ -> n * ρ * n

		# get the event operators
		eventmap = ρ -> dt * (mapreduce(+, Jt) do (j, γ) γ * j * ρ * j' end)
		
		# return the result of applying the operators 
		return nullmap(ρ) + eventmap(ρ)
	end
end"""


# function nullop(dt::Timescale, (j, γ)::Tuple{QOp, Rate}, tuples...)
# 	# concatenate all the tuples together; first is distinguished for dispatch
# 	tuples = [(j, γ), tuples...] 

# 	# get identity operator from left basis element of j
# 	Id = identityoperator(j.basis_l)

# 	# get the operator as √(I - dt * γ1 * j1' * j1 - dt * γ2 * j2' * j2 - ...)
# 	# square root is taken on op.data at the end
# 	op = DenseOperator(Id - dt * (mapreduce(+, tuples) do (j, γ) γ * j' * j end))
# 	return SparseOperator(op.basis_l, op.basis_r, sqrt(op.data))

# end

# function posmap(dt::Timescale, (j, γ)::Tuple{QOp, Rate}, tuples...)
# 	# concatenate all the tuples together; first is distinguished for dispatch
# 	tuples = [(j, γ), tuples...] 

# function nullmap(dt::Timescale, (j, γ)::Tuple{QOp, Function}, tuples...) 
# 	# concatenate all the tuples together; first is distinguished for dispatch
# 	tuples = [(j, γ), tuples...] 

# 	# get identity operator from left basis element of j
# 	Id = identityoperator(j.basis_l)

# 	# get the operator as √(I - dt * γ1 * j1' * j1 - dt * γ2 * j2' * j2 - ...)
# 	# square root is taken on op.data at the end
# 	op(t) = DenseOperator(Id - dt * (mapreduce(+, tuples) do (j, γ) γ(t) * j' * j end))
# 	return t -> SparseOperator(op.basis_l, op.basis_r, sqrt(op(t).data))
# end

# function nullmap(dt::Timescale, (j, γ)::Tuple{QOp, Function}, tuples...) 
# 	# concatenate all the tuples together; first is distinguished for dispatch
# 	tuples = [(j, γ), tuples...] 

# 	# get identity operator from left basis element of j
# 	Id = identityoperator(j.basis_l)

# 	# get the operator as √(I - dt * γ1 * j1' * j1 - dt * γ2 * j2' * j2 - ...)
# 	# square root is taken on op.data at the end
# 	op(t) = DenseOperator(Id - dt * (mapreduce(+, tuples) do (j, γ) γ(t) * j' * j end))
# 	return t -> SparseOperator(op.basis_l, op.basis_r, sqrt(op(t).data))
# end



