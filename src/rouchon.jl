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


    return Solution(collect(ts), ρs, r)
end # rouchon
