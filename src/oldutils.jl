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
