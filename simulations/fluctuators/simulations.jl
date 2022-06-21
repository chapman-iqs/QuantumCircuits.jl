const fluctuator_test_defaultpars = OrderedDict("N" => 100, "tf" => 10.0, "dt" => 20e-3)

function fluctuator_test(pars = fluctuator_test_defaultpars)

	N, tf, dt = pars.vals
	ts = collect(0.0:dt:tf)
	χtotal = makefluctuators(N, tf, dt)
	χtotal_traj = χtotal.(ts)

	return (ts, χtotal_traj)
end

"""
	T2simulation(v, N, tf, dt)

Keyword Arguments:
Ωz :: Real 		-- total strength of fluctuator Hamiltonian
Nf :: Integer	-- number of fluctuators
Nt :: Integer	-- number of trajectories
tf :: Real 		-- simulation duration
dt :: Real 		-- time step (1/dt is highest frequency fluctuator)
"""

const T2simulation_defaultpars = OrderedDict("Ωz" => 2,
											"Nf" => 20,
											"Nt" => 500,
											"tf" => 20.0,
											"dt" => 20e-3)

function T2simulation(pars::OrderedDict = T2simulation_defaultpars)

	Ωz, Nf, Nt, tf, dt = pars.vals

	# define time-dependent fluctuator Hamiltonian
	Hz(χtotal::Function) = t -> (Ωz/2) * χtotal(t) * σz

	# bayesian simulation
	ops, oplabels = qbasis, qlabels
	ψ0 = normalize(g + e)

	solutions = @showprogress map(i -> begin
		χtotal = makefluctuators(Nf, tf, dt)
		H = Hz(χtotal)
		return bayesian((0.0, tf), ψ0, H, [], []; dt=dt)
	end, 1:Nt)

	return solutions

	# exps = map(op -> [expectations(sol, op) for sol in solutions], ops)
	# return exps

end # function T2simulation
