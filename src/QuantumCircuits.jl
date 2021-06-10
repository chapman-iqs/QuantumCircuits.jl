module QuantumCircuits

using Reexport
@reexport using QuantumOpticsBase
using Distributions
using Distributed

δ(i,j) = Int(i == j);
sand(a,b) = a*b*a';

"""
    meas(dt, H; mclist=Tuple{Operator,Time,Float64}[],
                      mflist=Tuple{Function,Time,Float64}[],
                      clist=Operator[], flist=Function[])

Return increment function over a time step `dt` for diffusive monitored
evolution generated by Hamiltonian `H`, a list `mclist` of tuples specifying
constant measurement operators, and a list `mflist` of tuples specifying
time-dependent measurement operators. The tuples have the form `(m, τ, η)`,
where `m` is a (generally non-Hermitian) operator (or time-dependent function)
that specifies the measurement backaction (see below), `τ` is the timescale
of the measurement collapse, and `η` is the measurement efficiency (bounded
between 0 and 1).

These quantities are related as follows:
  - Γ = 1/(2τη) : ensemble measurement dephasing rate
  - Γm = 1/(2τ) : dephasing rate produced by averaging the collected signal
  - γ = Γ - Γm = (1-η)/(2τη) : residual dephasing rate from signal loss
  - ``m = m_o - i m_p``
  - ``m_o = (m + m')/2``  : measured observable (Hermitian part of `m`)
  - ``m_p = i(m - m')/2`` : phase-backaction generator (anti-Hermitian part of `m`)

The backaction takes the form of a three-step process:
  1. Sample `r` from a Gaussian distribution with mean ``⟨m_o⟩`` and variance τ/dt
  2. Conditionally update the state with a purity-preserving Kraus operator:
  - ``M_r = (dt/2πτ)^{1/4} e^{ -i m_p r dt/(2τ) - dt(r - m_o)^2/(4τ) }``
  3. Apply residual Lindblad dephasing evolution, including
  - `m` with rate ``γ = (1-η)/(2τη)``
  - the dephasing operators in `clist` and `flist`
  - natural Hamiltonian evolution `H`

Uses the "jump no-jump" method to efficiently approximate the residual
Lindblad dephasing as a composition of Hamiltonian evolution, jumps,
and no-jump informational backaction. Assumes no time-dependence in `alist`,
and small dt.  [Physical Review A **92**, 052306 (2015)]

### Returns:
  - (t, ρ(t)::Operator) -> (ρ(t+dt)::Operator, rlist::Float64...)

"""
function meas(dt::Float64, H0, J0::Array, C0::Array; rdo=Array[], ts=[])
    # Assemble readout generating functions and Kraus operators
    ros = Function[]
    gks = Function[]

    H = length(methods(H0)) > 0 ? H0 : t -> H0
    J = map(j -> length(methods(j)) > 0 ? j : t -> j, J0)
    C = map(c -> length(methods(c)) > 0 ? c : t -> c, C0)

    for c in C
        if length(rdo) == 0 
            push!(ros, readout(dt, c)) 
        end
        push!(gks, gausskraus(dt, c))
    end

    L = lind(dt, H, clist=[], flist=J)

    # Increment that samples each readout, applies all Kraus operators
    # then applies Lindblad dephasing (including Hamiltonian evolution)
    (t, ρ) -> begin
        rs = length(rdo) > 0 ? map(ro -> ro[argmin(abs.(ts .- t))], rdo) : map(ro -> ro(t, ρ), ros)
        gs = map(z -> z[2](t, z[1]), zip(rs,gks))
        ρ1 = foldr(sand, gs; init=ρ);
        return (L(t, ρ1/tr(ρ1)), rs) 
    end
end

function readout(dt, m::Function)
    σ = sqrt(1/dt)
    (t, ρ) -> let mo = (m(t) .+ m(t)')/2;
                        σ*randn() + real(expect(ρ, mo)) end
end

function gausskraus(dt, m::Function)
    v = dt/2
    (t, r) -> let mo = (m(t) .+ m(t)') / 2
                        mo2 = mo^2 / 2
                        exp(DenseOperator((r*v) * m(t) - v*mo2)) end
end

@inline function trajectory(inc::Function, ts, ρ; fn::Function=ρ->ρ, dt=1e-4)
    
    # probe record size
    _, rs1 = inc(ts[1], ρ)

    # init
    ρ0 = ρ
    dy0 = [0.0 for i in rs1]
    ρs = [fn(ρ0)]
    dy = [dy0]

    for t in ts[2:end]
        ρ, rs = inc(t, ρ)
        push!(ρs, fn(ρ))
        push!(dy, rs)
    end
    
    dy = collect(eachrow(hcat(dy...)))

    return (ts, ρs, dy)
end

function bayesian(tstep::Tuple, ρ, H0, J0::Array, C0::Array; fn=ρ->ρ, dt=1e-4, dy=[])
    ts = range(first(tstep), last(tstep), step=dt)
    return trajectory(meas(dt, H0, J0, C0; rdo=dy, ts=ts), ts, ρ; fn=fn, dt=dt)
end 

# Jump-nojump Lindblad propagator
"""
    lind(dt[, H]; clist=QOp[], flist=Function[])

Return increment function over a time step `dt` for Lindblad dissipative
evolution generated by (an optional) Hamiltonian `H` (operator or function),
a list `clist` of constant dissipative operators, and a list `flist` of
time-dependent functions returning dissipative operators.

Uses the "jump no-jump" method to efficiently approximate the exact
Lindblad propagator as a composition of Hamiltonian evolution, jumps,
and no-jump informational backaction. Assumes small dt.
[Physical Review A **92**, 052306 (2015)]

### Returns:
  - (t, ρ(t)::Operator) -> ρ(t+dt)
"""
function lind(dt; clist=QOp[], flist=Function[])
    ns = Function[]
    ds = Function[]
    # Construct operations for constant operators
    if isempty(clist)
        push!(ns, (t, ρ) -> ρ)
    else
        Id = identityoperator(first(clist).basis_l)
        op = DenseOperator(Id - dt * mapreduce(a -> a' * a, +, clist))
        n::Operator = SparseOperator(op.basis_l, op.basis_r, sqrt(op.data))
        push!(ns, (t, ρ) -> n * ρ * n)
        push!(ds, (t, ρ) -> mapreduce(a -> a * ρ * a', +, clist) * dt)
    end
    # Construct operations for time-dependent operators
    if isempty(flist)
        push!(ns, (t, ρ) -> ρ)
    else
        function nf(t)
            Id = identityoperator(first(flist)(t).basis_l)
            op = DenseOperator(Id - dt * mapreduce(a -> a(t)' * a(t), +, flist))
            return SparseOperator(op.basis_l, op.basis_r, sqrt(op.data))
        end
        push!(ns, (t, ρ) -> nf(t) * ρ * nf(t))
        push!(ds, (t, ρ) -> mapreduce(a -> a(t) * ρ * a(t)', +, flist) * dt)
    end
    push!(ds, (t, ρ) -> last(ns)(t, first(ns)(t, ρ)))
    (t, ρ) -> mapreduce(f -> f(t, ρ), +, ds)
end
function lind(dt, H; clist=QOp[], flist=Function[])
    # Rely on Hamiltonian to specify type of H
    h = ham(dt, H)
    # Apply Hamiltonian first, then the Lindblad increment
    (t, ρ) -> lind(dt, clist=clist, flist=flist)(t, h(t, ρ))
end

# Hamiltonian propagation
"""
    ham(dt, H::Operator; ket=false)

Return increment function for Hamiltonian evolution generated
by `H` over a time step `dt`.

Uses an exact (dense) matrix exponential, assuming no time-dependence.

### Returns:
  - ket=true  : (t, ψ::QKet) -> u * ψ
  - ket=false : (t, ρ::Operator)  -> u * ρ * u'

"""
function ham(dt, H::Operator)
    u::Operator = exp( -im * dt * DenseOperator(H))
    ut = u'
    (t, ρ::Operator) -> u * ρ * ut
end
function ham(dt, H::Function)
    (t, state) -> ham(dt, H(t))(t, state)
end

"""
    rouchon(T, ρ0, H, J, C; <keyword arguments>)

Arguments:

T :: tuple (ti,tf)
ρ0 :: initial density matrix
H :: time-dependent Hamiltonian
J :: deterministic collapse operators
C :: stochastic collapse operators

Keyword Arguments:

dt :: time step; default dt=1e-4
dy :: record; default dy=[], i.e. simulation generates record time series. 
            record should be input in the shape [dy_1,...,dy_Nc] given 
            length(C)=Nc collapse operators. Records should have shape
            dy_c[t] indexing the mth trajectory at time T[n]
fn : ρ → Any :: eval function (e.g. to return expectation values instead of density matrices)
"""

function rouchon(T0, ρ, H0, J0::Array, C0::Array; fn=ρ->ρ, dt=1e-4, dy=[])
    T = range(T0[1], T0[2], step=dt)
    Id = identityoperator(ρ.basis_l)
    Nn = length(T)
    Nj = length(J0)
    Nc = length(C0)
    
    ρs = Any[]

    H = length(methods(H0)) > 0 ? H0 : t -> H0
    J = map(j -> length(methods(j)) > 0 ? j : t -> j, J0)
    C = map(c -> length(methods(c)) > 0 ? c : t -> c, C0)
    
    # randomly sample noise time series for EACH stochastic collapse operator C
    if length(dy) == 0
        sim = true
        dist = Normal(0, sqrt(dt))
        dW=[]
        for c in 1:Nc
            push!(dW,rand(dist,Nn))
            push!(dy,zeros(Nn))
        end
    else
        sim = false
    end

    for n in 1:Nn
        t = T[n]
        M = Id - im*H(t)*dt

        # iterate over deterministic collapse operators J
        D = 0Id
        for j in 1:Nj
            M += -0.5J[j](t)'*J[j](t)*dt
            D += J[j](t)*ρ*J[j](t)'*dt
        end
            
        # initialize dy value, if simulation
        if sim==true
            for c in 1:Nc
                dy[c][n] = real(tr(C[c](t)*ρ + ρ*C[c](t)')*dt) + dW[c][n]
            end
        end

        # iterate over stochastic collapse operators C
        for c in 1:Nc
            M += C[c](t)*conj.(dy[c][n])
            D += -C[c](t)*ρ*C[c](t)'*dt
                
            # nested sum    
            for s in 1:Nc
                M += 0.5C[c](t)*C[s](t)*(conj.(dy[c][n])*conj.(dy[s][n]) - δ(c,s)*dt)
            end
        end

        # update ρ according to Rouchon method
        ρ = M*ρ*M' + D
        ρ = ρ/tr(ρ)
        push!(ρs, fn(ρ))
    end
        
    return (T,ρs,dy)
end

"""
    ensemble(solve, T, ρ0, H, J, C; <keyword arguments>)

solve :: solver function
T :: tuple (ti,tf)
ρ0 :: initial density matrix
H :: time-dependent Hamiltonian
J :: deterministic collapse operators
C :: stochastic collapse operators

Keyword Arguments:

dt :: time step (default dt=1e-4)
record :: (default record=[]), i.e. simulation generates record time series. 
            record should be input in the shape [rec_1,...,rec_Nc] given 
            length(C)=Nc collapse operators. Records should have shape
            rec[m][c,n] indexing the mth trajectory at time T[n]
N :: number of trajectories (default N=10)
"""

function ensemble(solve, T, ρ0, H, J, C; dt=1e-4, record=[], N=10, onstart=x->x, kwargs...)
    data = pmap(m -> begin
        onstart(m)
        dy = length(record) >= m ? record[m] : []
        tt, ρs, dy = solve(T, ρ0, H, J, C; dt=dt, dy=dy, kwargs...)
        return (ρs, dy, tt)
    end, 1:N)

    trajectories = collect(ρs for (ρs, dy) in data)
    record = collect(dy for (ρs, dy) in data)
    ts = data[1][3]

    return (ts, trajectories, record)
end

function coarse_grain(fine=[]; n=2)
    coarse = []
    for i in 1:length(fine)
        if i < n
            push!(coarse, mean(fine[1:i]))
        else
            push!(coarse, mean(fine[i-(n-1):i]))
        end
    end
    coarse
end

function subselect(a=[]; n=2)
    a[filter(x -> x%n==0, eachindex(a))]
end

export δ, rouchon, ensemble, meas, trajectory, bayesian, coarse_grain, subselect

end # module
