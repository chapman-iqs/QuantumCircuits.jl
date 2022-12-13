# Utilities
##

"""
    comm(a, b)
Commutator product of `a` and `b`
### Returns:
    - Anti-Hermitian Operator: a * b - b * a
"""
comm(a,b) = a * b' - b * a'

"""
    scomm(a)
Superoperator for commutator with operator `a`.
Assumes Hermitian superket.
### Returns:
    - Superoperator: scomm(a) * superket(b) == superket(a * b - b * a')
"""
scomm(a) = spre(a) - spost(a')

"""
    acomm(a, b) = a ⊕ b
Anti-commutator of operators `a` and `b`.
### Returns:
    - Hermitian operator: a * b' + b * a'
"""
acomm(a, b) = a * b' + b * a'

"""
    sacomm(a)
Superoperator for anticommutator with operator `a`.
Assumes Hermitian superket.
### Returns:
    - Superoperator: scomm(a) * superket(b) == superket(a * b + b * a')
"""
sacomm(a) = spre(a) + spost(a')

"""
    sand(a, b)
Sandwich `b` operator with `a`.
### Returns:
    - Operator: a * b * a'
"""
sand(a, b) = a * b * a'

"""
    ssand(a)
Superoperator for sandwich with operator `a`.
### Returns:
    - Superoperator: ssand(a) * superket(b) == superket(a * b * a')
"""
ssand(a) = spre(a) * spost(a')

"""
    diss(a)
Dissipation function for `a` action.
### Returns:
    - Function: ρ -> sand(a, ρ) - acomm(at*a, ρ)/2
"""
diss(a) = ρ -> sand(a, ρ) - acomm(a'*a, ρ)/2

"""
    sdiss(a)
Dissipation superoperator for `a` action.
### Returns:
    - ssand(a) - sacomm(at*a)/2
"""
sdiss(a) = ssand(a) - sacomm(a'*a)/2


supdate(A::SuperOperator, ρ::State) = A * ρ
supdate(a::QOp, ρ::QOp) = supdate(ssand(a), ρ)
# supdate(a::QOp, ρ::Ket) = supdate(ssand(a), ρ)


# Superoperator Hamiltonian evolution
"""
    sham(dt::Timescale, H::QOp)

Return increment function using a superoperator for Hamiltonian
evolution generated by `H` over a Timescale step `dt`.
Uses an exact (dense) matrix exponential, assuming no Timescale-dependence.

### Returns:
    - U::SuperOperator : Evolution superoperator
    or
    - (t, ρ) -> U(t, ρ)
    or
    - t -> U(t)
"""
function sham(dt::Timescale, H::QOp)::SuperOperator
    u::QOp = SparseOperator(exp( -im * dt * DenseOperator(H)))
    U = ssand(u)
end
function sham(dt::Timescale, H::Function)
    if applicable(H, tt, ρρ)
        (t::Timescale, ρ::State) ->  sham(dt, H(t, ρ))
    elseif applicable(H, tt)
        t::Timescale -> sham(dt, H(t))
    else
        error("Hamiltonian arguments do not match input.")
    end
end


# Jump-nojump Lindblad propagator
"""
	lind(dt[,
	 	  J::Vector{Tuple{QOp, Rate}}];
		  clist=QOp[],
		  flist=Function[]
		  )

Returns increment function over a Timescale step dt for Lindblad dissipative
evolution generated by dissipative operators √γ * j, corresponding to
pairs (j::QOp, γ::Rate) in J.

If J is input as argument to lind, the function sorts the pairs into two
lists: a list `clist` of constant dissipative operators and a list `flist`
of Timescale-dependent functions returning dissipative operators.

### Returns:
  - when J input : (t, ρ) -> lind(dt; clist=clist, flist=flist)
  - else : (t, ρ(t)::QOp) -> L(ρ), where L denotes the Lindblad superoperator
"""
function slind(dt, J)

    if isempty(J)
        return (t, ρ) -> spre(identityoperator(ρ.basis_l))
    end

    # extract identityoperator
    (j1, γ1) = J[1]
    id = identityoperator(j1.basis_l)

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
	(t, ρ) -> slind(dt; id=id, clist=clist, flist=flist, fρlist=fρlist)(t, ρ)
end
function slind(dt; id = identityoperator(SpinBasis(1//2)), clist=QOp[], flist=Function[], fρlist=Function[])
    # first, construct the null map from all operations
    nc, nf, nfρ = map(list -> snullop(dt, list...; id=id), [clist, flist, fρlist])
    ntotal(t, ρ) = nc(t, ρ) * nf(t, ρ) * nfρ(t, ρ)

    # then, construct the null superoperator
    Ns(t, ρ) = spre(ntotal(t, ρ)) * spost(ntotal(t, ρ))

    # second, construct the event maps from all operations
    jc, jf, jfρ = map(list -> seventop(dt, list...; id=id), [clist, flist, fρlist])

    # construct the events superoperator
    Js(t, ρ) = jc(t, ρ)+ jf(t, ρ) + jfρ(t, ρ)

    # the total superoperator is the sum of the two
    return (t, ρ) -> Ns(t, ρ) + Js(t, ρ)
end

snullop(dt::Timescale; id=identityoperator(SpinBasis(1//2))) = (t, ρ) -> id

function snullop(dt::Timescale, J...; id=identityoperator(SpinBasis(1//2))) 

    # if isempty(J)
    #     return (t, ρ) -> id
    # end

    (j, γ) = J[1]

    if applicable(γ, tt, ρρ)
        op(t, ρ) = DenseOperator(id - dt * (mapreduce(+, J) do (j, γ) γ(t, ρ) * j' * j end))
        return (t, ρ) -> SparseOperator(op(t, ρ).basis_l, op(t, ρ).basis_r, sqrt(op(t, ρ).data))

    elseif applicable(γ, tt)
        op(t) = DenseOperator(id - dt * (mapreduce(+, J) do (j, γ) γ(t) * j' * j end))
        return (t, ρ) -> SparseOperator(op(t).basis_l, op(t).basis_r, sqrt(op(t).data))

    elseif typeof(γ) <: Rate
        op = DenseOperator(id - dt * (mapreduce(+, J) do (j, γ) γ * j' * j end))
        return (t, ρ) -> SparseOperator(op.basis_l, op.basis_r, sqrt(op.data))

    else
        error("Lindblad decay rate received improper argument types. Define γ::Rate, γ(t::Timescale), or γ(t::Timescale, ρ::State).")
    end
end

function seventop(dt::Timescale, J...; id=identityoperator(SpinBasis(1//2)))

    if isempty(J)
        return (t, ρ) -> 0 * spre(id)
    end

    (j, γ) = J[1]

    if applicable(γ, tt, ρρ)
        return (t, ρ) -> dt * (mapreduce(+, J) do (j, γ) γ(t, ρ) * ssand(j) end)

    elseif applicable(γ, tt)
        return (t, ρ) -> dt * (mapreduce(+, J) do (j, γ) γ(t) * ssand(j) end)

    elseif typeof(γ) <: Rate
        return (t, ρ) -> dt * (mapreduce(+, J) do (j, γ) γ * ssand(j) end)

    else
        error("Lindblad decay rate received improper argument types. Define γ::Rate, γ(t::Timescale), or γ(t::Timescale, ρ::State).")
    
    end

end

function smeas(dt::Timescale, C; sample=true)

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
			# evaluate readout generating functions at Timescale t for state ρ
			readout = map(ro -> ro(t, ρ), ros)

			# evaluate each gausskraus operator at Timescale t for the acquired readout
			gs = [g(t, r) for (g, r) in zip(gks, readout)]

			# construct a superoperator from all the Kraus operators 
            g = foldr(*, gs)
            return ssand(g), readout
		end

	else
		# same as above, but takes readout as argument
		(t, ρ, readout) -> begin

			gs = [g(t, r) for (g, r) in zip(gks, readout)]
            g = foldr(*, gs)
            return ssand(g)
		end

	end
end


# Superoperator Lindblad propagator
"""
    slind(dt::Timescale, H::QOp; clist=QOp[], flist=Function[])
Return increment function over a Timescale step `dt` for Lindblad dissipative
evolution generated by (an optional) Hamiltonian `H` (operator or function),
a list `clist` of constant dissipative operators, and a list `flist` of
Timescale-dependent functions returning dissipative operators.
Uses direct matrix exponentiation of the total superoperator for
the increment.
### Returns:
    - (t::Timescale, ρvec) -> u * ρvec : Superoperator for total evolution over dt
    """
sslind(dt::Timescale, H) = sham(dt, H)
function sslind(dt::Timescale, H::QOp; clist=QOp[], flist=Function[], fρlist=Function[])
    h = -im*scomm(H)
    init = DenseSuperOperator(h.basis_l, h.basis_r)
    if isempty(flist) && isempty(fρlist)
        let l = DenseSuperOperator(h.basis_l, h.basis_r, (h + mapreduce(sdiss, +, clist; init=init)).data)
            return exp(dt*l);
        end
    elseif isempty(fρlist)
        let l(t::Timescale) = DenseSuperOperator(h.basis_l, h.basis_r, (h + mapreduce(sdiss, +, clist; init=init) +
                                mapreduce(f -> sdiss(f(t)), +, flist; init=init)).data)
            return t::Timescale -> exp(dt*l(t))
        end
    else
        let l(t::Timescale, ρ::State) = DenseSuperOperator(h.basis_l, h.basis_r, (h + mapreduce(sdiss, +, clist; init=init) +
                                mapreduce(f -> sdiss(f(t)), +, flist; init=init) + mapreduce(f -> sdiss(f(t, ρ)), +, fρlist; init=init)).data)
            return (t::Timescale, ρ::State) -> exp(dt*l(t, ρ))
        end
    end
end
function sslind(dt::Timescale, H::Function; clist=QOp[], flist=Function[], fρlist=Function[])
    if applicable(H, tt, ρρ)
        if !isempty(fρlist)
            (t::Timescale, ρh::State, ρl::State) -> sslind(dt, H(t, ρh), clist=clist, flist=flist, fρlist=fρlist)(t, ρl)
        elseif !isempty(flist)
            (t::Timescale, ρh::State, ρl::State) -> sslind(dt, H(t, ρh), clist=clist, flist=flist, fρlist=fρlist)(t)
        else
            (t::Timescale, ρh::State, ρl::State) -> sslind(dt, H(t, ρh), clist=clist, flist=flist, fρlist=fρlist)
        end
    elseif applicable(H, tt)
        if !isempty(fρlist)
            (t::Timescale, ρh::State, ρl::State) -> sslind(dt, H(t), clist=clist, flist=flist, fρlist=fρlist)(t, ρl)
        elseif !isempty(flist)
            (t::Timescale, ρh::State, ρl::State) -> sslind(dt, H(t), clist=clist, flist=flist, fρlist=fρlist)(t)
        else
            (t::Timescale, ρh::State, ρl::State) -> sslind(dt, H(t), clist=clist, flist=flist, fρlist=fρlist)
        end        
    else
        error("Invalid argument types to Hamiltonian.")
    end
end
function sslind(dt, H, J)

        # preprocess lindblad arguments
    clist = []
    flist = []
    fρlist = []

    for (j, γ) in J
        if applicable(γ, tt, ρρ)
            push!(fρlist, (t, ρ) -> √γ(t, ρ) * j)
        elseif applicable(γ, tt)
            push!(flist, t -> √γ(t) * j)
        elseif typeof(γ) <: Rate
            push!(clist, √γ * j)
        else
            error("Improper argument types to γ.")
        end
    end
    return sslind(dt, H; clist=clist, flist=flist, fρlist=fρlist)
end

"""
    strajectory(Meas::Function, H::SuperOperator, L::SuperOperator, ts, ρ; dt=1e-3)

TBW
"""
function strajectory(ts, ρ, Meas::Function, U::SuperOperator, L::SuperOperator; dt=1e-3, r0=[0.0])

	# initialize
    ρs = [ρ]
	readouts = Readouts([r0])

	for t in ts[2:end]
        M, ro = Meas(t, ρ)
        ρ = normalize(foldr(supdate, [U, L, M]; init=ρ))
		push!(ρs, ρ)
        push!(readouts, ro)
	end
    
    records = trans(readouts)

	return Solution(collect(ts), ρs, records)
end

function strajectory(ts, ρ, Meas::Function, L::SuperOperator; dt=1e-3, r0=[0.0])

	# initialize
    ρs = [ρ]
	readouts = Readouts([r0])

	for t in ts[2:end]
        M, ro = Meas(t, ρ)
        ρ = normalize(L * M * ρ)
		push!(ρs, ρ)
        push!(readouts, ro)
	end
    
    records = trans(readouts)

	return Solution(collect(ts), ρs, records)
end

"""
HSys is a (time / state-dep) Hamiltonian map for the system, taking the forward-estimated state as input.
LFil is a (time / state-dep) map containing Lindblad and Hamiltonian terms for the filter, taking the forward-estimated state as input.

"""
function forwardtrajectory(ts, (ρs0, ρf0), Meas::Function, USys::Function, LFil::Function; dt=1e-3, r0=[0.0], td = 0.0)

	# initialize system
    ρs = ρs0
    ρs_list = [ρs]

    # initialize filter
    ρf = ρf0
    ρf_list = [ρf]

    # initialize forward-estimated state
    ρe = ρf0
    ρe_list = [ρe]

    # delay and buffer maps
    delay = Int64(round(td / dt)) + 1
    buffermaps = [spre(identityoperator(ρs0.basis_l)) for i in 1:(delay-1)]

    # initialize readout
	readouts = Readouts([r0])
    

	for t in ts[2:delay]

        # get the forward-estimated state as a function of the initial state and the intervening maps
        ρe = foldr(*, buffermaps; init=ρf0)

        # get the maps
        M, ro = Meas(t, ρs) # measure the SYSTEM to get the bayesian update map
        U = USys(t, ρe)        # get SYSTEM update map based on forward-estimated state
        L = LFil(t, ρe, ρf)        # get FILTER update map based on forward-estimated state

        # update the buffermaps to remember the Hamiltonian operation that was chosen
        pushfirst!(buffermaps, sparse(L))
        pop!(buffermaps)
        

        # state updates
        ρs = normalize(U * M * ρs)    # SYSTEM ρ based on the measurement and Hamiltonian maps
        ρf = normalize(L * M * ρf)   # FILTER ρ based on measurement and Lindblad maps
		push!(ρs_list, ρs)
        push!(ρf_list, ρf)
        push!(ρe_list, ρe)
        push!(readouts, ro)
	end

    for i in (delay + 1):length(ts)
		t = ts[i]

        # get the forward-estimated state as a function of the filter state td ago and the intervening maps
        ρe = foldr(*, buffermaps; init=ρf_list[i - delay])

        # get the maps
        M, ro = Meas(t, ρs) # measure the SYSTEM to get the bayesian update map
        U = USys(t, ρe)        # get SYSTEM update map based on forward-estimated state
        L = LFil(t, ρe, ρf)        # get FILTER update map based on forward-estimated state

        # update the buffermaps to remember the Hamiltonian operation that was chosen
        push!(buffermaps, sparse(L))
        popfirst!(buffermaps)
        

        # state updates
        ρs = normalize(U * M * ρs)    # SYSTEM ρ based on the measurement and Hamiltonian maps
        ρf = normalize(L * M * ρf)   # FILTER ρ based on measurement and Lindblad maps
		push!(ρs_list, ρs)
        push!(ρf_list, ρf)
        push!(ρe_list, ρe)
        push!(readouts, ro)
	end
    
    records = trans(readouts)

	return Solution(collect(ts), ρs_list, records), Solution(collect(ts), ρf_list), Solution(collect(ts), ρe_list)
end

function forwardbayesian((t0, tf), ρ, (Hs, Hf)::Tuple, Jf, C; dt=1e-3, td::Timescale = 0.0)
	# Hs, Hf = convertham.(Hpair)

	ts = range(t0, tf, step=dt)
	r0::Readout = [0.0 for i in 1:length(C)] # initial, empty Readout object

	# define superoperator functions for state update
    USys = sham(dt, Hs)
    LFil = sslind(dt, Hf, Jf)
    Meas = smeas(dt, C; sample=true)

    # type conversion
    ρs0 = ρ
    ρf0 = ρ
	if length(Jf) > 0 && typeof(ρ) <: Ket
		ρf0 = dm(ρf0)
	end

    forwardtrajectory(ts, (ρs0, ρf0), Meas, USys, LFil; dt=1e-3, r0=r0, td=td)
end # bayesian

function sbayesian((t0, tf), ρ, H::QOp, J, C; dt=1e-3)

    r0::Readout = [0.0 for i in 1:length(C)] 
    if length(J) > 0 && typeof(ρ) <: Ket
		ρ = dm(ρ)
	end

    U = sham(dt, H)
    L = slind(dt, J)(tt, ρρ) # assumes const lindblad
    Meas = smeas(dt, C)
ρe
    return strajectory(collect(t0:dt:tf), ρ, Meas, U, L; dt=dt)
end

function ssbayesian((t0, tf), ρ, H::QOp, J, C; dt=1e-3)

    r0::Readout = [0.0 for i in 1:length(C)] 
    if length(J) > 0 && typeof(ρ) <: Ket
		ρ = dm(ρ)
	end

    L = sslind(dt, H, J) # assumes const lindblad
    Meas = smeas(dt, C)

    return strajectory(collect(t0:dt:tf), ρ, Meas, L; dt=dt)
end

