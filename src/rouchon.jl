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
# non-feedback trajectory generating its own measurement record
function rouchon((t0, tf), ρ, H0, J0, C0; fn=ρ->ρ, dt=1e-4, records=Record[])

    ### initialize
    ts = range(t0, tf, step=dt)                         # time series
    # r0::Readout = [0.0 for i in 1:length(C0)]            # initial, empty Readout object
    H = isa(H0, Function) ? H0 : t -> H0

    # convert ρ to a density matrix
    if typeof(ρ) <: Ket
        ρ = dm(ρ)
    end

    # loop elements
    ρ0 = ρ
    fnρs = [fn(ρ0)]
    # readouts = Readouts([r0])

    # initalize rouchon
    # Id = isa(ρ, Ket) ? identityoperator(ρ.basis) : identityoperator(ρ.basis_l)
    Id = identityoperator(ρ.basis_l)
    Nt = length(ts)
    Nj = length(J0)
    Nc = length(C0)

    # modify J
    J = map(J0) do (j, Γ)
            isa(Γ, Function) ? (t -> √Γ(t) * j) : (t -> √Γ * j)
        end
	# J = []
	# for (j, Γ) in J0
	# 	push!(J, isa(J, Function) ? √Γ * j : (t -> √Γ * j) ) 
    # end

    # modify C and sample random numbers for readout
    sim = isempty(records)
    dist = Normal(0, √dt)
    dys = sim ? map(ctup -> rand(dist, Nt), C0) : records_to_dys(records, C0, dt)
	C = []

	for (i, (c, Γ, η)) in enumerate(C0)
		m = √Γ * η * c
		push!(C, isa(m, Function) ? m : (t -> m))
		# push!(C, length(methods(c)) > 0 ? (c,τm,η) : (t -> c,τm,η))
	end

    # loop over times
    for (n, t) in enumerate(ts[2:end])
        M = Id - im * H(t) * dt

        # iterate over deterministic collapse operators J
        D = 0Id
        for j in 1:Nj
            M += -0.5J[j](t)' * J[j](t) * dt
            D += J[j](t) * ρ * J[j](t)' * dt
        end

        # initialize dy value, if simulation
        if sim
            for (i, c) in enumerate(C)
                dys[i][n] += real(tr(c(t) * ρ + ρ * c(t)') * dt)/sqrt(2)
            end
        end

        # iterate over stochastic collapse operators C
        for c in 1:Nc
            # println("length(dy) = ", length(dy))
            M += C[c](t) * dys[c][n]
            # D += -C[c](t)*ρ*C[c](t)'*dt

            # nested sum
            for s in 1:Nc
                M += 0.5C[c](t) * C[s](t) * (dys[c][n] * dys[s][n] - δ(c,s) * dt)
            end
        end

        # update ρ according to Rouchon method
        ρ = M*ρ*M' + D
        ρ = ρ / tr(ρ)
        push!(fnρs, fn(ρ))
    end

	# revert to r (normalized) version of record
    records = sim ? dys_to_records(dys, C0, dt) : records
    # for (i, (c, Γ, η)) in enumerate(C0)
    #     τm = 1/(2Γ*η)
    #     push!(records, dy[i] .* sqrt(τm) ./ dt)
    # end

    return Solution(collect(ts), fnρs, records)
end # rouchon

function record_to_dy(record::Record, (c, Γ, η), dt)
    τm = 1/(2Γ*η)
    return record .* dt ./ sqrt(τm)
end
function dy_to_record(dy, (c, Γ, η), dt)::Record
    τm = 1/(2Γ*η)
    return dy .* sqrt(τm) ./ dt
end
function records_to_dys(records::Vector{Record}, C0, dt)
    map(zip(records, C0)) do (record, ctup)
        record_to_dy(record, ctup, dt)
    end
end
function dys_to_records(dys, C0, dt)::Vector{Record}
    map(zip(dys, C0)) do (dy, ctup)
        dy_to_record(dy, ctup, dt)
    end
end

    

# old rouchon method
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

"""