
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
# comment
function ensemble(solve, (t0, tf), ρ, H::Union{Function, QOp}, J, C; dt=1e-4, records=Vector{Record}[], N=10, batch_size=10, ops=[], showprogress=true, kwargs...)

	p = Progress(N, enabled=showprogress)
	solutions = progress_pmap(1:N, progress=p; batch_size=batch_size) do i
		rs = length(records) > 0 ? records[i] : Record[]
		return solve((t0, tf), ρ, H, J, C; dt=dt, records=rs, kwargs...)
	end

	return length(ops) == 0 ?
					solutions :
					map(op -> [expectations(sol, op) for sol in solutions], ops)
end
function ensemble(solve, (t0, tf), ρ, (Hs, Hf)::Tuple, J, C; dt=1e-4, N=10, batch_size=10, ops=[], showprogress=true, kwargs...)

	p = Progress(N, enabled=showprogress)
	solutions = progress_pmap(1:N, progress=p; batch_size=batch_size) do i
		return solve((t0, tf), ρ, (Hs, Hf), J, C; dt=dt, kwargs...)
	end

	if length(ops) == 0
		return solutions
	else
		sys_exps = map(op -> [expectations(sys, op) for (sys, _) in solutions], ops)
		fil_exps = map(op -> [expectations(fil, op) for (_, fil) in solutions], ops)
		return (sys_exps, fil_exps)
	end
end
"Alternatively, ensemble can take a vector of Hamiltonians. It will run bayesian once for each Hamiltonian. Does not yet implement
inputting records."
function ensemble(solve, (t0, tf), ρ, Hlist::Vector, J, C; dt=1e-4, N=10, batch_size=10, ops=[], showprogress=true, kwargs...)

	if length(Hlist) != N
		warn(string(	"Your number of trajectories is ", N,
											", but your vector of Hamiltonians has length ", length(Hlist),
											". Please ensure these lengths match by inputting N as a keyword argument.",
											" (This is a style convention but helps to avoid errors",
											" when comparing feedback and non-feedback Hamiltonians.)"
										  ))
	end

	p = Progress(length(Hlist), enabled=showprogress)
	solutions = progress_pmap(Hlist, progress=p; batch_size=batch_size) do H
	    return solve((t0, tf), ρ, H, J, C; dt=dt, kwargs...)
	end

	return length(ops) == 0 ?
					solutions :
					map(op -> [expectations(sol, op) for sol in solutions], ops)
end
