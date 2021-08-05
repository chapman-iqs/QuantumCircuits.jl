### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 2c9f5f4e-9236-4625-b7f3-e27a0747895d
begin
	using Random
	using Statistics
	using PyPlot
	using Distributions
	using QuantumOptics 	# note: for some reason, QuantumCircuits has to be 
							# "used" last
end

# ╔═╡ ffd8cd8a-d08e-11eb-1b1e-6b0db734faf5


# ╔═╡ d0c661c3-8423-43aa-9997-fe29fb49becd
md" # Linear feedback stabilization"

# ╔═╡ b871f40b-ee15-486d-aa83-ccf6c07ebb72
md"""
In this interactive notebook, we implement the linear feedback stabilization described in [1].
"""

# ╔═╡ 8a5d25cb-3038-43bc-9ae2-83b68ff11340
md"""
## Problem setup

We begin by defining a two-level Hilbert space for the system. `QuantumCircuits.jl` uses `QuantumOptics.jl` as its backend: `SpinBasis`, `sigmax`, `identityoperator` and so forth are `QuantumOptics.jl` functions.
"""

# ╔═╡ 52152b81-68f9-4b37-9e9d-4df1a0fef4b1
begin
	# Basis
	q = SpinBasis(1//2)

	# Operators, using convention that |-z> is ground state
	σx = sigmax(q)
	σy = sigmay(q)
	σz = sigmaz(q)
	σp = sigmap(q)
	σm = sigmam(q)
	id = identityoperator(q)
end

# ╔═╡ 31edbe57-d620-4a53-a7f3-8256906b0411
md"""
### Description of the evolution

The Hamiltonian is a Rabi drive modulated by the feedback from the measurement readout:

$\hat H_c = \hbar \big[\Delta_0 + \Delta_1 \tilde{r}(t - T_d) \big] \frac{\hat \sigma_\phi}2$

where

$\hat \sigma_\phi = \cos \phi \hat \sigma_x + \sin \phi \hat \sigma_y,$

and $\hat \sigma_x$ and $\hat \sigma_y$ are the corresponding Pauli operators. Thus there is a base Rabi frequency $\Delta_0$ modulated by a linear feedback term $\Delta_1 \tilde{r} (t-T_d)$.

"""

# ╔═╡ 8394e44b-c1ac-438d-847e-55c722bd3cd8
md" ## Functions"

# ╔═╡ 0602e66f-701c-471d-9c1d-8770b61944c1
md" Look at how code already handles readout-dependent measurement operators; then do something similar for Hamiltonian."

# ╔═╡ 87faf926-9862-4f4e-8b38-f5a6aac11026
function gausskraus(dt, m::Function)
    v = dt/2
    (t, r) -> let mo = (m(t) .+ m(t)') / 2
                        mo2 = mo^2 / 2
                        exp(DenseOperator(conj(r) * v * m(t) - v*mo2)) end
end

# ╔═╡ b566b767-8aa3-44ed-900e-7762f7130c6d
function readout(dt, m::Function)
    dist = Normal(0, sqrt(dt))
    (t, ρ) -> expect(ρ, m(t)) + (rand(dist) + im*rand(dist))/(dt * √2)
end

# ╔═╡ 03c74991-d45f-455f-b0a8-81957d55f97a
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

# ╔═╡ 979f6729-253a-4769-8971-5ed7d96e9eea
function lindf(dt; clist=QOp[], flist=Function[])
    ns = Function[]
    ds = Function[]
    # Construct operations for constant operators
    if isempty(clist)
        push!(ns, (t, ρ, r) -> ρ)
    else
        Id = identityoperator(first(clist).basis_l)
        op = DenseOperator(Id - dt * mapreduce(a -> a' * a, +, clist))
        n::Operator = SparseOperator(op.basis_l, op.basis_r, sqrt(op.data))
        push!(ns, (t, ρ, r) -> n * ρ * n)
        push!(ds, (t, ρ, r) -> mapreduce(a -> a * ρ * a', +, clist) * dt)
    end
    # Construct operations for time-dependent operators
    if isempty(flist)
        push!(ns, (t, ρ, r) -> ρ)
    else
        function nf(t)
            Id = identityoperator(first(flist)(t).basis_l)
            op = DenseOperator(Id - dt * mapreduce(a -> a(t)' * a(t), +, flist))
            return SparseOperator(op.basis_l, op.basis_r, sqrt(op.data))
        end
        push!(ns, (t, ρ, r) -> nf(t) * ρ * nf(t))
        push!(ds, (t, ρ, r) -> mapreduce(a -> a(t) * ρ * a(t)', +, flist) * dt)
    end
    push!(ds, (t, ρ, r) -> last(ns)(t, first(ns)(t, ρ, r)))
    (t, ρ, r) -> mapreduce(f -> f(t, ρ, r), +, ds)
end


# ╔═╡ 09e06118-9102-4ff4-9854-79ccccce95ee
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

# ╔═╡ f3cc4146-dd3d-4993-9804-7b1d48240d9a
function ham(dt, H::Function)
    (t, state) -> ham(dt, H(t))(t, state)
end

# ╔═╡ 888950c5-2f5e-441d-a067-4addc4bd1026
function lind(dt, H; clist=QOp[], flist=Function[])
    # Rely on Hamiltonian to specify type of H
    h = ham(dt, H)
    # Apply Hamiltonian first, then the Lindblad increment
    (t, ρ) -> lind(dt, clist=clist, flist=flist)(t, h(t, ρ))
end

# ╔═╡ c28bb0e9-1951-408b-932f-c1f1be42e174
function hamf(dt, H::Function)
    (t, state, r) -> ham(dt, H(t, r))(t, state, r)
end

# ╔═╡ 90bf5d08-61f8-46b4-9af1-80f1792afbcd
function lindf(dt, H; clist=QOp[], flist=Function[])
    # Rely on Hamiltonian to specify type of H
    h = hamf(dt, H)
    # Apply Hamiltonian first, then the Lindblad increment
    (t, ρ, r) -> lindf(dt, clist=clist, flist=flist)(t, h(t, ρ, r))
end

# ╔═╡ 202d6a71-4c81-4fa8-bc12-21fad3266e57
function meas(dt::Float64, H0, J0::Array, C0::Array; rdo=Array[], ts=[], td=0, sample=true)

	feedback = applicable(H0, ts[1], [1])
	
	# Assemble readout generating functions and Kraus operators
	ros = Function[]
	gks = Function[]
	
	H = length(methods(H0)) > 0 ? H0 : t -> H0 # if H0 is already a function of t, set H = H0; if H0 = const, define h(t) a constant function of t
	J = map(j -> length(methods(j)) > 0 ? j : t -> j, J0) # same, for each element of J0
	C = map(c -> length(methods(c)) > 0 ? c : t -> c, C0) # same, for each element of C0
	
	for c in C
		if length(rdo) == 0 && sample
			push!(ros, readout(dt, c))
		end
		push!(gks, gausskraus(dt, c))
	end
	
	L = feedback ? lindf(dt, H, clist=[], flist=J) : lind(dt, H, clist=[], flist=J)
	
# Increment that samples each readout, applies all Kraus operators
	# then applies Lindblad dephasing (including Hamiltonian evolution)
	(t, ρ) -> begin
		# for each experimental record (corresponding to each observable), get the value closest to t
		# else, for each simulated record, get the value corresponding to (t,ρ)
		# gives an array (?) of readout values for each observable, at time t.
		rs = length(rdo) > 0 ? map(ro -> ro[argmin(abs.(ts .- t))], rdo) : map(ro -> ro(t, ρ), ros)
		
		if feedback
			rd = length(rdo) > 0 ? map(ro -> ro[argmin(abs.(ts .- (t-td)))], rdo) : 					map(ro -> ro(t-td, ρ), ros) end
	
		# for each pair in `zip(rs, gks)`, get the gausskraus operator `z[2]` as a function of `(t,z[1])`, 
		# where `z[1]` is the readout value at `t`
		gs = map(z -> z[2](t, z[1]), zip(rs, gks))
	
		# apply each element of gs via `sand` function (sandwich) to ρ
		ρ1 = foldr(sand, gs; init=ρ);
		return feedback ? (L(t, ρ1/tr(ρ1), rd), rs) : (L(t, ρ1/tr(ρ1)), rs)
	end
end

# ╔═╡ e818500b-58d7-47a1-bae1-413960bfb1ca
function trajectory(inc::Function, ts, ρ; fn::Function=ρ->ρ, dt=1e-4, td=td)

    # probe record size
    _, rs1 = inc(ts[1], ρ)

    # init
    ρ0 = ρ
    dy0 = [0.0im for i in rs1]
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

# ╔═╡ f3ec5c67-50d2-4744-b7a3-bf1ba7f7bb2e
function bayesian(tstep::Tuple, ρ, H0, J0::Array, C0::Array; fn=ρ->ρ, dt=1e-4, dy=[], td=0, sample=true)
    ts = range(first(tstep), last(tstep), step=dt)
    return trajectory(meas(dt, H0, J0, C0; rdo=dy, ts=ts, sample=sample), ts, ρ; fn=fn, dt=dt, td=td)
end

# ╔═╡ c52046a0-dd05-4905-b51b-4134232875b9
md" ## Parameters and operators"

# ╔═╡ af7ddcfe-45aa-4baf-b919-317aece1297a
begin
	ρ0 = dm(spindown(q)) # initial state
	dt = 1e-2  # integration time-step
	td = dt*15
	
	Δ0 = 2π
	Δ1 = Δ0/30
	τ = 3.0
	η = 0.3
	Γ = 1/(2τ)
end

# ╔═╡ 86ae37a7-51c3-44d0-b401-ea0b20310e09
begin
	dist = Normal(0, sqrt(dt))
	rand(dist)
end

# ╔═╡ 1e66d35d-fa19-4cf0-a70a-1e94d0d93799
begin
	# Hamiltonian defined at a certain time t, 
	# takes in readout r in some observable with time delay td
	# H0(t, r::Array) = (Δ0 + Δ1*r[1])*σx 
	H0 = Δ0*σx
	J0 = [√Γ*σz]
	C0 = [√(Γ*η)*σz];
end

# ╔═╡ a7dfc96a-6f51-4762-9472-46f0ab142a53
applicable(H0, 1, [1])

# ╔═╡ 995d0307-2be2-4671-87cb-7f2bef8e952e
md" ## Testing"

# ╔═╡ cf030661-bbb7-4959-bc7d-9238bc1669d8
begin
	Random.seed!(1)
	sol1 = bayesian((0, 4τ), ρ0, H0, [J0], [C0]; dt=dt, td=td)
end

# ╔═╡ a0f7774a-47ca-4121-9b18-be491e8c6314
md"""
### References

[1] T. L. Patti, A. Chantasri, L. P. García-Pintos, A. N. Jordan, and J. Dressel, Linear Feedback Stabilization of a Dispersively Monitored Qubit, Phys. Rev. A 96, 022311 (2017).

"""

# ╔═╡ Cell order:
# ╠═ffd8cd8a-d08e-11eb-1b1e-6b0db734faf5
# ╠═2c9f5f4e-9236-4625-b7f3-e27a0747895d
# ╠═d0c661c3-8423-43aa-9997-fe29fb49becd
# ╠═b871f40b-ee15-486d-aa83-ccf6c07ebb72
# ╠═8a5d25cb-3038-43bc-9ae2-83b68ff11340
# ╠═52152b81-68f9-4b37-9e9d-4df1a0fef4b1
# ╠═31edbe57-d620-4a53-a7f3-8256906b0411
# ╟─8394e44b-c1ac-438d-847e-55c722bd3cd8
# ╟─0602e66f-701c-471d-9c1d-8770b61944c1
# ╠═86ae37a7-51c3-44d0-b401-ea0b20310e09
# ╠═87faf926-9862-4f4e-8b38-f5a6aac11026
# ╠═b566b767-8aa3-44ed-900e-7762f7130c6d
# ╠═03c74991-d45f-455f-b0a8-81957d55f97a
# ╠═888950c5-2f5e-441d-a067-4addc4bd1026
# ╠═979f6729-253a-4769-8971-5ed7d96e9eea
# ╠═90bf5d08-61f8-46b4-9af1-80f1792afbcd
# ╠═09e06118-9102-4ff4-9854-79ccccce95ee
# ╠═f3cc4146-dd3d-4993-9804-7b1d48240d9a
# ╠═c28bb0e9-1951-408b-932f-c1f1be42e174
# ╠═202d6a71-4c81-4fa8-bc12-21fad3266e57
# ╠═f3ec5c67-50d2-4744-b7a3-bf1ba7f7bb2e
# ╠═e818500b-58d7-47a1-bae1-413960bfb1ca
# ╟─c52046a0-dd05-4905-b51b-4134232875b9
# ╠═af7ddcfe-45aa-4baf-b919-317aece1297a
# ╠═1e66d35d-fa19-4cf0-a70a-1e94d0d93799
# ╠═a7dfc96a-6f51-4762-9472-46f0ab142a53
# ╟─995d0307-2be2-4671-87cb-7f2bef8e952e
# ╠═cf030661-bbb7-4959-bc7d-9238bc1669d8
# ╟─a0f7774a-47ca-4121-9b18-be491e8c6314
