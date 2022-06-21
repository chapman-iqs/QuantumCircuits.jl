"calculate switching for a fluctuator of frequency γ, up to final time tf"
function switching_times(γ, tf)

	dist = Exponential(1/γ)
	T = 0.0
	Ts = []
	while T < tf
		T += rand(dist)
		push!(Ts, T)
	end
	return Ts
end

"fluctuator sign as a function of time, for fluctuator of frequency γ, up to final time tf"
function fluctuator_sign(γ, tf)

	χ0 = rand([-1,1]) # initial sign
	Ts = switching_times(γ, tf)

	return t -> let j = (findfirst(t .< Ts) - 1)
					return χ0 * (-1)^j
	end
end


"""
	makefluctuators(N, tf, dt)

samples fluctuators between 1/tf and 1/(20dt), where dt is the smallest resolvable timescale
of the simulation.

Keyword Arguments:
N :: Integer	-- number of fluctuators
tf :: Real 		-- simulation duration
dt :: Real 		-- time step (1/dt is highest frequency fluctuator)
"""
function makefluctuators(N, tf, dt)

	"frequency range of fluctuators"
	(fmin, fmax) = 1/tf, 1/(20dt)
	flog = range(log(fmin), log(fmax), length=N)
	fi = exp.(flog)

	"setting up fluctuators"
	χs = [] # array of fluctuator signs
	for f in fi
		γ = 2π * f
		χ = fluctuator_sign(γ, tf)
		push!(χs, χ)
	end

	χtotal(t) = sum(map(χ -> χ(t), χs)) / N

	return χtotal
end # function makefluctuators
