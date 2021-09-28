### A Pluto.jl notebook ###
# v0.16.0

using Markdown
using InteractiveUtils

# ╔═╡ 7d43b95c-d3a4-11eb-16d3-47148906d472
begin
	cd("/Users/sachagreenfield/Desktop/GitHub/QuantumCircuits.jl")
	import Pkg
	Pkg.activate(".")
	using Random
	using Statistics
	using Distributions
	using PyPlot
	using QuantumCircuits
end

# ╔═╡ be7d340a-4ef1-4727-a531-7d6e5c568ad9
md" # Feedback: matrix exponentiation"

# ╔═╡ b9bd68f2-df48-4072-8269-1898b7cf1b15
md"""
In this interactive notebook, I test the feedback matrix exponentiation code to compare to bloch equations and verify that `QuantumCircuits.jl` is working as expected.
"""

# ╔═╡ e66e6723-16fe-4d8d-ae1f-f8237f8e8897
md"""
## Problem setup

We begin by defining a two-level Hilbert space for the system. `QuantumCircuits.jl` uses `QuantumOptics.jl` as its backend: `SpinBasis`, `sigmax`, `identityoperator` and so forth are `QuantumOptics.jl` functions.
"""

# ╔═╡ 601a8b56-2f41-4b2c-9ea2-6b26a52b4e26
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

# ╔═╡ 18c29abc-1f35-4b5a-bb27-a491c02cc98f
md"""
### Description of the evolution

The Hamiltonian is a Rabi drive modulated by the feedback from the measurement readout:

$\hat H_c = \hbar \big[\Delta_0 + \Delta_1 \tilde{r}(t - T_d) \big] \frac{\hat \sigma_\phi}2$

where

$\hat \sigma_\phi = \cos \phi \hat \sigma_x + \sin \phi \hat \sigma_y,$

and $\hat \sigma_x$ and $\hat \sigma_y$ are the corresponding Pauli operators. Thus there is a base Rabi frequency $\Delta_0$ modulated by a linear feedback term $\Delta_1 \tilde{r} (t-T_d)$.

"""

# ╔═╡ 01f57775-b647-4fea-8e96-0b8c8ceeff05
md" ### Parameters and operators"

# ╔═╡ 071641ac-a5a8-4d37-b531-88b1e083416c
md" #### Bloch simulation (generate record)"

# ╔═╡ c5de979e-de3b-4a20-9fc4-649851a311fa
ideal = true

# ╔═╡ eae605ed-f411-4f33-8066-bd8f01fc8a2d
begin
	# all times given in μs
	# initial state
	(x0,y0,z0) = (0., 0.3, 0.91)
	ρ0 = DenseOperator(0.5*(id + x0*σx + y0*σy + z0*σz))
	
	dt = 0.5e-3  # integration time-step
	td = 0 # 200e-3 # delay time
	η = 1 # 0.41 # measurement efficiency
	
	θs = 3π/10 # target angle on Bloch sphere
	ϕ = π # fixes plane of oscillations
	σϕ = cos(ϕ)*σx + sin(ϕ)*σy
	
	τm =  0.2 # measurement time
	Γm = 1/(2τm) # measurement rate
	T = (0, 4) #(0, 10τm) # simulation duration
	td = 0 # time delay for feedback
	
	Rs = 1 # 0.64 # radius of target
	
	# feedback drive parameters
	Δ0 = ideal ? -sin(2θs)/(4τm) : -sin(2θs)/(4τm*Rs^2) 
	Δ1 = ideal ? sin(θs)/τm : sin(θs)/(τm*Rs)
		
	T1 = 40 # energy decay time
	T2 = 60 # environmental dephasing time

	# corresponding rates
	Γ1 = 1/(2T1)
	Γ2 = 1/T2

	
end

# ╔═╡ 6294f81c-ea9d-4500-9751-b45f8a348639
begin
	ts = range(first(T), last(T), step=dt)
	ztar = [Rs*cos(θs) for i in 1:length(ts)]
	ytar = [Rs*sin(θs) for i in 1:length(ts)]
	xtar = zeros(length(ts))
end

# ╔═╡ ee541c05-c187-4b43-a803-2255e254efe5
begin
	# Hamiltonian defined at a certain time t, 
	# takes in readout r in some observable with time delay td
	H0(t, r::Array) = (Δ0 + Δ1*r[1])*σϕ/2
	# H0 = Δ0*σϕ/2 
	C0 = [(σz, τm, η)]
	J0 = [√((1-η)*Γm)*σz]
	#J = ideal ? [√((1-η)*Γm)*σz] : [√((1-η) * Γm) * σz, √Γ1*σm, √Γ2*σz] 
end

# ╔═╡ 8aa08bfb-ac91-4e5f-9fb2-dbce02a38b8a
begin
	
	title_str = ideal ? "ideal Rabi oscillation with feedback, general" : "non-ideal Rabi oscillation with feedback, general";
	title_str_hc = ideal ? "ideal Rabi oscillation with feedback, hard-coded" : "non-ideal Rabi oscillation with feedback, hard-coded";
	
end

# ╔═╡ ea558387-68dc-4dff-ae04-f7a877379510
md" ### Bloch simulation from Patti et al."

# ╔═╡ f2892559-e623-46d4-9ff2-f04fd9253734
function blochevolve(; rs=[])
	# readout distribution
	sim = (length(rs) == 0)
	ts = range(first(T), last(T), step=dt)
	
	# first time step
	(xn, yn, zn) = (x0, y0, z0)	
	(xs, ys, zs) = ([x0], [y0], [z0])

	
	if sim
		dist = Normal(0, sqrt(τm/dt))
		push!(rs, 0) end
		

	for i in 2:length(ts)
		
		# sample from distribution		
		if sim
			r = zn + rand(dist)
			push!(rs, r)
		else
			r = rs[i] end
		
 		# update equations
		pn = cosh(r*dt/τm) + zn*sinh(r*dt/τm)
		xn1 = xn/pn
		yn1 = yn/pn
		zn1 = (zn*cosh(r*dt/τm) + sinh(r*dt/τm))/pn
		Δ = (Δ0 + Δ1*r) # the Rabi-drive depends on the readout (instantaneous here)
		
		yn2 = yn1*cos(dt*Δ) + zn1*sin(dt*Δ)
		zn2 = zn1*cos(dt*Δ) - yn1*sin(dt*Δ)
		
		# if ideal 
		xn = xn1*exp(-dt*(1-η)/(2τm*η))
		yn = yn2*exp(-dt*(1-η)/(2τm*η))
		zn = zn2

		# else 
		# 	xn = xn1*exp(-dt/(2T1) - dt/T2 - dt*(1-η)/(2τm*η))
		# 	yn = yn2*exp(-dt/(2T1) - dt/T2 - dt*(1-η)/(2τm*η))
		# 	zn = zn2*exp(-dt/T1) - (1 - exp(-dt/T1)) 
		# end

		# store values
		push!(xs, xn)
		push!(ys, yn)
		push!(zs, zn)
		

	end
	
	ρs = 0.5*(1 .+ xs.^2 + ys.^2 + zs.^2)
	
	return ts, rs, (xs, ys, zs, ρs)
	
	
end

# ╔═╡ 4487690f-349c-461c-a5b1-f9a6be9897bd
begin
	Random.seed!(2)
	(tt, r, (xx, yy, zz, rr)) = blochevolve()
	tt = collect(tt)
end

# ╔═╡ 15613717-87c1-45d8-b194-1c7b1085d4f7
md" ### Matrix exponentiation"

# ╔═╡ 840beeea-a731-4a38-b36b-e2ab848147d9
begin
	H = length(methods(H0)) > 0 ? H0 : t -> H0
	J = map(j -> length(methods(j)) > 0 ? j : t -> j, J0) # same, for each element of J0
	C = []
	for (c, τm, η) in C0
		push!(C, length(methods(c)) > 0 ? (c,τm,η) : (t -> c,τm,η) ) end
end

# ╔═╡ bc8e8a5f-6476-402f-9cc6-0af722a6c033
ρ0.data

# ╔═╡ e7353f43-899d-4899-86b5-1303676c07ec
ρ0.data[2,2]

# ╔═╡ 2c1f165f-cabd-4599-afc7-348799754f5e
function expevolve(; rs=[])
	# readout distribution
	sim = (length(rs) == 0)
	ts = range(first(T), last(T), step=dt)
	
	# first time step
	(xn, yn, zn) = (x0, y0, z0)	
	(xs, ys, zs, ρs) = ([x0], [y0], [z0], [tr(ρ0*ρ0)])
	ρ = ρ0

	
	if sim
		dist = Normal(0, sqrt(τm/dt))
		push!(rs, 0) end
		

	for i in 2:length(ts)
		
		# sample from distribution		
		if sim
			r = zn + rand(dist)
			push!(rs, r)
		else
			r = rs[i] end
		
 		# update equations
		
		# measurement backaction
		M = exp(r' * dt * DenseOperator(σz) / 2τm)
		ρ = M * ρ * M'
		ρ = ρ / tr(ρ) 
		
		# Hamiltonian evolution
		H = (Δ0 + Δ1*r) * σϕ/2
		U = exp(-im * dt * DenseOperator(H))
		ρ = U * ρ * U'
		
		# Lindblad evolution
		if ideal
			ρ.data[1,2] = ρ.data[1,2] * exp(-dt * (1 - η)/(2τm * η))
			ρ.data[2,1] = ρ.data[2,1] * exp(-dt * (1 - η)/(2τm * η))
			
		else
			ρ.data[1,2] = ρ.data[1,2] * exp(-dt/2T1 - dt/T2 - dt*(1 - η)/(2τm * η))
			ρ.data[2,1] = ρ.data[2,1] * exp(-dt/2T1 - dt/T2 - dt*(1 - η)/(2τm * η))
		end
		
		# store values
		push!(xs, real(expect(ρ, σx)))
		push!(ys, real(expect(ρ, σy)))
		push!(zs, real(expect(ρ, σz)))
		push!(ρs, tr(ρ*ρ))
		

	end
	
	return ts, rs, (xs, ys, zs, ρs)
	
	
end

# ╔═╡ 90511e52-c343-490b-81b2-1502c546e16b
(_, _, (xxe, yye, zze, rre)) = expevolve(rs=r)

# ╔═╡ 9e731ac9-e11c-427b-9569-4a925799042a
md" #### Testing QuantumCircuits.jl functions"

# ╔═╡ ab459f70-e40a-499f-bd05-248270d7d4dd
md" ##### Functions"

# ╔═╡ 2487c71c-2d00-4ac0-8cec-7d3f2e80ab3e
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

# ╔═╡ 8c422daa-fbf6-4a9e-83ce-3523948aa069
function ham(dt, H::Operator)

    u::Operator = exp( -im * dt * DenseOperator(H))
    ut = u'
    (t, ρ::Operator) -> u * ρ * ut

end

# ╔═╡ 3ad727fc-6955-4bb9-b9bf-cb34b63c8c70
function ham(dt, H::Function)
	feedback = applicable(H, 1, [1])
	return feedback ?
			(t, state, r) -> ham(dt, H(t, r))(t, state) :
			(t, state) -> ham(dt, H(t))(t, state)
end

# ╔═╡ 5d2e29b6-3f8e-4e66-9442-c325af182a6a
function lind(dt, H; clist=QOp[], flist=Function[])

	feedback = applicable(H, 1, 1)
    # Rely on Hamiltonian to specify type of H
    h = ham(dt, H)
    # Apply Hamiltonian first, then the Lindblad increment
	return feedback ?
			(t, ρ, r) -> lind(dt, clist=clist, flist=flist)(t, h(t, ρ, r)) :
			(t, ρ) -> lind(dt, clist=clist, flist=flist)(t, h(t, ρ))
end

# ╔═╡ f9e49ee3-d8d2-4d92-8ce0-f9d07401d61c
md" ##### Testing: no time-dependence"

# ╔═╡ 4ee435fb-e786-4a9f-8f46-1820e6ce9936
function expevolve2(; rs=[])
	# readout distribution
	sim = (length(rs) == 0)
	ts = range(first(T), last(T), step=dt)
	
	# first time step
	(xn, yn, zn) = (x0, y0, z0)	
	(xs, ys, zs, ρs) = ([x0], [y0], [z0], [tr(ρ0*ρ0)])
	ρ = ρ0

	
	if sim
		dist = Normal(0, sqrt(τm/dt))
		push!(rs, 0) end
		

	for i in 2:length(ts)
		t = ts[i]
		
		# sample from distribution		
		if sim
			r = zn + rand(dist)
			push!(rs, r)
		else
			r = rs[i] end
		
 		# update equations
		
		# measurement backaction
		M = exp(r' * dt * DenseOperator(σz) / 2τm)
		ρ = M * ρ * M'
		ρ = ρ / tr(ρ) 
		
		# Hamiltonian evolution
		H = (Δ0 + Δ1*r) * σϕ/2
		
		# Lindblad evolution
		L = lind(dt, H, clist=[], flist=J)
		
		ρ = L(t, ρ/tr(ρ))
		
		
		# store values
		push!(xs, real(expect(ρ, σx)))
		push!(ys, real(expect(ρ, σy)))
		push!(zs, real(expect(ρ, σz)))
		push!(ρs, tr(ρ*ρ))
		

	end
	
	return ts, rs, (xs, ys, zs, ρs)
	
	
end

# ╔═╡ d7d0ec12-242f-4c1c-875e-f960a2d9c88f
(_, _, (xxe2, yye2, zze2, rre2)) = expevolve2(rs=r)

# ╔═╡ 406317be-0a5a-4895-ad79-70fbc683f79b
md" ##### Testing: functionally defined Hamiltonian"

# ╔═╡ 8f11676f-ab29-4a30-a23c-4e65039f3a44
Hf(t, r) = (Δ0 + Δ1*r[1]) * σϕ/2

# ╔═╡ 9a21aa2e-8034-41b0-a87b-d261efe70ee2
rs = map(el -> [el], r)

# ╔═╡ 295bef3d-8f73-453e-86dc-b8e61cdded83
function expevolve3(; rs=[])
	# readout distribution
	sim = (length(rs) == 0)
	ts = range(first(T), last(T), step=dt)
	
	# first time step
	(xn, yn, zn) = (x0, y0, z0)	
	(xs, ys, zs, ρs) = ([x0], [y0], [z0], [tr(ρ0*ρ0)])
	ρ = ρ0

	
	if sim
		dist = Normal(0, sqrt(τm/dt))
		push!(rs, [0]) end
		

	for i in 2:length(ts)
		t = ts[i]
		
		# sample from distribution		
		if sim
			r = zn + rand(dist)
			push!(rs, [r])
		else
			r = rs[i] end
		
 		# update equations
		
		# measurement backaction
		M = exp(r[1]' * dt * DenseOperator(σz) / 2τm)
		ρ = M * ρ * M'
		ρ = ρ / tr(ρ) 
		
		# Hamiltonian evolution
		
		# Lindblad evolution
		L = lind(dt, Hf, clist=[], flist=J)
		
		ρ = L(t, ρ/tr(ρ), r)
		
		
		# store values
		push!(xs, real(expect(ρ, σx)))
		push!(ys, real(expect(ρ, σy)))
		push!(zs, real(expect(ρ, σz)))
		push!(ρs, tr(ρ*ρ))
		

	end
	
	return ts, rs, (xs, ys, zs, ρs)
	
	
end

# ╔═╡ d429e614-a7d4-43e9-ac9b-d7dcfdfb4efa
tt

# ╔═╡ fcefa57c-d828-47b0-bc98-c107c9078946
argmin(abs.(tt .- td))

# ╔═╡ 5957c035-687a-48d0-8c9a-ad0489c3f240
function expevolve4(; rs=[], td=0)
	# readout distribution
	sim = (length(rs) == 0)
	ts = range(first(T), last(T), step=dt)
	
	# first time step
	(xn, yn, zn) = (x0, y0, z0)	
	(xs, ys, zs, ρs) = ([x0], [y0], [z0], [tr(ρ0*ρ0)])
	ρ = ρ0

	
	if sim
		dist = Normal(0, sqrt(τm/dt))
		push!(rs, [0]) end
		
	td_index = argmin(abs.(ts .- td))
	
		for i in 2:td_index
		t = ts[i]
		rf = r[i - td_index]
		
		# sample from distribution		
		if sim
			r = zn + rand(dist)
			push!(rs, [r])
		else
			r = rs[i] end
		
 		# update equations
		
		# measurement backaction
		M = exp(r[1]' * dt * DenseOperator(σz) / 2τm)
		ρ = M * ρ * M'
		ρ = ρ / tr(ρ) 
		
		# Hamiltonian evolution
		
		# Lindblad evolution
		L = lind(dt, Hf, clist=[], flist=J)
		
		ρ = L(t, ρ/tr(ρ), 0)
		
		
		# store values
		push!(xs, real(expect(ρ, σx)))
		push!(ys, real(expect(ρ, σy)))
		push!(zs, real(expect(ρ, σz)))
		push!(ρs, tr(ρ*ρ))
		

	end

	for i in (td_index+1):length(ts)
		t = ts[i]
		rf = rs[i - (td_index - 1)]
		
		# sample from distribution		
		if sim
			r = zn + rand(dist)
			push!(rs, [r])
		else
			r = rs[i] end
		
 		# update equations
		
		# measurement backaction
		M = exp(r[1]' * dt * DenseOperator(σz) / 2τm)
		ρ = M * ρ * M'
		ρ = ρ / tr(ρ) 
		
		# Hamiltonian evolution
		
		# Lindblad evolution
		L = lind(dt, Hf, clist=[], flist=J)
		
		ρ = L(t, ρ/tr(ρ), rf)
		
		
		# store values
		push!(xs, real(expect(ρ, σx)))
		push!(ys, real(expect(ρ, σy)))
		push!(zs, real(expect(ρ, σz)))
		push!(ρs, tr(ρ*ρ))
		

	end
	
	return ts, rs, (xs, ys, zs, ρs)
	
	
end

# ╔═╡ d3ec6758-b27d-4c1e-a7de-7383f0f9027e
(_, _, (xxe3, yye3, zze3, rre3)) = expevolve4(rs=rs)

# ╔═╡ 75bc7334-aac0-4d5e-895d-894ba49ce492
yye3

# ╔═╡ 3258df22-de38-4ab5-92b7-0a828bc32155
md" ## Utilities "

# ╔═╡ e59530ae-7f27-4836-a352-9d0a08d62451
function rms(ser1::Array, ser2::Array)
	l = min(length(ser1), length(ser2))
	sqrt(sum((ser1[1:l] - ser2[1:l]).^2))/l
end

# ╔═╡ 0fa84d7f-59ff-4a27-821f-2114464129a6
expects = ρ -> collect(real(expect(ρ, s)) for s in [σx,σy,σz,ρ]) # ρ -> [<x>,<y>,<z>,<ρ>]

# ╔═╡ a589fb58-fc99-493d-bca4-fd21dc6a78d8
function blochs(ρs::Array)
	# Get Bloch components
	evs0 = expects.(ρs);
	xx,yy,zz,rr = [map(x -> x[i], evs0) for i in 1:4];

	(xx, yy, zz, rr)
	
end

# ╔═╡ c9a3a847-d910-4433-926a-d1bfd012f248
function blochs(sol)
	(tt, ρt, _) = sol

	# Get Bloch components
	evs0 = expects.(ρt);
	xx,yy,zz,ρρ = [map(x -> x[i], evs0) for i in 1:4];

	(collect(tt), xx, yy, zz, ρρ)
	
end

# ╔═╡ 739984d0-91a6-45d1-bd2a-0561bfb3e4d5
purity(x,y,z) = 0.5*(1 + x^2 + y^2 + z^2)

# ╔═╡ 8f40ea32-6f75-4119-b90c-c72b7749360f
R(x,y,z) = sqrt(x^2 + y^2 + z^2)

# ╔═╡ 8bb2e4ea-e480-4e03-be07-24a5c4d6c780
function subseries(rec, T, dt; scale=2)
	ts = collect(range(first(T), last(T), step=dt))
	tts = subselect(real(coarse_grain(ts; n=scale)); n=scale)
	(ti, tf) = (tts[1], tts[end])
	dtt = dt * scale
	
	subrec = subselect(real(coarse_grain(rec; n=scale)); n=scale)
	
	(tts, subrec)
	
end

# ╔═╡ 16aed026-dcfd-432f-8028-32c7ef3afa5b
# colorscheme
begin
	colorscheme = "Paired"
	cmap = plt.matplotlib.cm.get_cmap(colorscheme)
	colors=collect(map(x -> cmap(x), 0:11))
	md" `colorscheme`"
end

# ╔═╡ 69855c0b-c7e3-4e6c-bffa-7d5e47a55b1b
# Plotting
function plot_solution(sol; plot_title="Rabi Oscillation")
	
	close("all")
	
	(tt, ρt, dys) = sol
    
    # Get Bloch components
    evs0 = expects.(ρt);
    xx,yy,zz,ρρ = [map(x -> x[i], evs0) for i in 1:4];
    
    # Plot Bloch components vs. time
    
    p = plot(tt, xx, color=colors[2], label=L"$x$")
    plot(tt, yy, color=colors[4],label=L"$y$")
    ax = gca()
    ax.set_ylim([-1.1,1.1]) 
    plot(tt, zz, color=colors[6], label=L"$z$")
	plot(tt, ρρ, color=colors[8], label=L"Tr $\rho^2$")
    xlabel(L"$t$")
    ylabel("Bloch coordinates")
    title(plot_title)
	legend()
	gcf()

end

# ╔═╡ b41b51dc-0763-4ad0-90e7-8e54fcc7d916
# Plotting
function plot_blochs((tt, blochs); plot_title="Rabi Oscillation")
	close("all")
	
	(xx,yy,zz) = blochs
	pp = map(i -> purity(xx[i], yy[i], zz[i]), 1:length(tt))
    
    # Plot Bloch components vs. time
    
    p = plot(tt, xx, color=colors[2], label=L"$x$")
    plot(tt, yy, color=colors[4],label=L"$y$")
    ax = gca()
    ax.set_ylim([-1.1,1.1]) 
    plot(tt, zz, color=colors[6], label=L"$z$")
	plot(tt, pp, color=colors[8], label=L"Tr $\rho^2$")
    xlabel(L"$t$")
    ylabel("Bloch coordinates")
    title(plot_title)
    legend()
    gcf()
end

# ╔═╡ 6afc5215-2b51-425f-aecc-b20e15d5128a
function record_histograms(records...; plot_title="record histogram", labels=[]::Array, density=false)
	close("all")
	
	μσs = []
	hist_colors = []
	hist_labels = []
	
	for i in 1:length(records)
		label = i > length(labels) ? i : labels[i]
		dys = records[i]
		
		# get mean and std dev for real part
		(μ, σ) = params(fit(Normal, dys))
		push!(μσs, map(p -> round(p, digits=4), (μ, σ)))

		# make histogram		
		n, bins, patches = hist(dys, 50, density=density, 					 									facecolor=colors[2i], alpha=1, label=label)
		push!(hist_colors, colors[2i])
		

	end
	
	# write down (μ, σ) pairs as text boxes
	μσ_strings = map(μσ -> string("(μ, σ) = (", μσ[1], ", ", μσ[2], ")\n"), μσs)
	ax = gca()
	for i in 1:length(μσ_strings)
		str = μσ_strings[i]

		ax.text(0.05, 1 - 0.05i, str, transform=ax.transAxes, fontsize=10,
			verticalalignment="top", color=hist_colors[i])
		
	end
	
	
	legend()
	xlabel("value (arbitrary units)")
	ylabel(density ? "relative frequency" : "frequency")
	title("record histograms")
	gcf()
		
end

# ╔═╡ a1c5942e-4c27-4cbb-8fa8-b47d9449c523
# Plotting
function plot_timeseries(tt::Array, series...; plot_title="time series", xlab=L"$t$", ylab="arbitrary units", labels=[]::Array, ylims=[], colorpairs=false, kwargs...)
	close("all")
	ser_colors(i) = colorpairs ? colors[i] : colors[2i]
	
	label(i) = i > length(labels) ? i : labels[i]
    
    # Plot records vs. time
	ser = series[1]
	p = plot(tt, ser, color=ser_colors(1), label=label(1), kwargs...)
	ax = gca()
	if length(ylims) > 0
		ax.set_ylim(ylims)  end
	
	if length(series) > 1
		
		for i in 2:length(series)
			ser = series[i]
			plot(tt, ser, color=ser_colors(i),label=label(i), kwargs...) 
		end
		
	end
	
    xlabel(xlab)
    ylabel(ylab)
    title(plot_title)
    legend()
    gcf()
	
end

# ╔═╡ aec50388-2682-44cd-aab8-6474cc863f31
# Plotting
function plot_timeseries(ttseries...; plot_title="time series", xlab=L"$t$", ylab="arbitrary units", labels=[]::Array, colorpairs=false)
	close("all")
	ser_colors(i) = colorpairs ? colors[i] : colors[2i]
	
	label(i) = i > length(labels) ? i : labels[i]
    
    # Plot records vs. time
	(tt, ser) = ttseries[1]
	p = plot(tt, ser, color=ser_colors(1), label=label(1))
	ax = gca()
	
	if length(ttseries) > 1
		
		for i in 2:length(ttseries)
			(tt, ser) = ttseries[i]
			plot(tt, ser, color=ser_colors(i),label=label(i)) 
		end
		
	end
	
    xlabel(xlab)
    ylabel(ylab)
    title(plot_title)
    legend()
    gcf()
	
end

# ╔═╡ a348c59c-acb5-428f-b544-d5e1038c56d6
plot_timeseries(tt, xtar, xx, ytar, yy, ztar, zz, rr; plot_title="hard-coded evolution", xlab=L"$t$", ylab="bloch coordinates", labels=["x target", "x", "y target", "y", "z target", "z", string("Tr ", L"$\rho^2$")], colorpairs=true, ylims=[0.5,0.9])

# ╔═╡ 7bc8459e-1652-497f-87da-e2e45b145df7
plot_timeseries((tt, r); plot_title="blochevolve simulated record")

# ╔═╡ 2a842a56-4d47-4a4e-9c3c-7f7ce5c1ac75
plot_timeseries(tt, xtar, xxe, ytar, yye, ztar, zze, rre; plot_title="hard-coded matrix exp evolution", xlab=L"$t$", ylab="bloch coordinates", labels=["x target", "x", "y target", "y", "z target", "z", string("Tr ", L"$\rho^2$")], colorpairs=true, ylims=[-1.1,1.1])

# ╔═╡ 477756ac-f4ff-471e-b73a-13183585592a
plot_timeseries(tt, xtar, xx, ytar, yy, ztar, zz, rr; plot_title="hard-coded evolution", xlab=L"$t$", ylab="bloch coordinates", labels=["x target", "x", "y target", "y", "z target", "z", string("Tr ", L"$\rho^2$")], colorpairs=true, ylims=[-1.1,1.1])

# ╔═╡ fe17c0d6-7195-4491-84b6-582d907eeb92
plot_timeseries(tt, xtar, xxe2, ytar, yye2, ztar, zze2, rre2; plot_title="matrix exp evolution using functions", xlab=L"$t$", ylab="bloch coordinates", labels=["x target", "x", "y target", "y", "z target", "z", string("Tr ", L"$\rho^2$")], colorpairs=true, ylims=[-1.1,1.1])

# ╔═╡ 43de5e01-e2e7-4302-8914-0fc444645db4
plot_timeseries(tt, xtar, xxe3, ytar, yye3, ztar, zze3, rre3; plot_title="matrix exp evolution using functions", xlab=L"$t$", ylab="bloch coordinates", labels=["x target", "x", "y target", "y", "z target", "z", string("Tr ", L"$\rho^2$")], colorpairs=true, ylims=[-1.1,1.1])

# ╔═╡ c541a88e-0a5a-48d2-b490-91ef7500dbe5
# Plotting
function plot_solutions((sol1,sol2); plot_title="Rabi Oscillation")
    close("all")
    
    tt1 = sol1[1]
    ρt1 = sol1[2]
    tt2 = sol2[1]
    ρt2 = sol2[2]
    
    # Get Bloch components
    evs1 = expects.(ρt1);
    x1,y1,z1,ρ1 = [map(x -> x[i], evs1) for i in 1:4];
    evs2 = expects.(ρt2);
    x2,y2,z2,ρ2 = [map(x -> x[i], evs2) for i in 1:4];
    
    # Plot Bloch components vs. time
	p = plot(tt1, x1, color=colors[2], linewidth=2, label=L"$x_{\eta = 0}$")
    plot(tt1, y1, color=colors[4], linewidth=2, label=L"$y_{\eta = 0}$")
    plot(tt1, z1, color=colors[6], linewidth=2, label=L"$z_{\eta = 0}$")
    plot(tt1, p1, color=colors[8], linewidth=2, label=L"(Tr $\rho^2)_{\eta = 0}$")
    plot(tt2, x2,  color=colors[1], linestyle="dashed", label=L"$x_{avg}$")
    plot(tt2, y2, color=colors[3], linestyle="dashed", label=L"$y_{avg}$")
    plot(tt2, p2, color=colors[7], linewidth=2, linestyle="dashed", label=L"(Tr $\rho^2)_{avg}$")
  
	ax = gca()
	ax.set_ylim([-1.1,1.1]) 
	xlabel(L"$t$")
	ylabel("Bloch coordinates")
	title(plot_title)
	legend()
	gcf()
end

# ╔═╡ c0c7bcf3-3a8e-489d-9bda-4113764190ec
function plot_evals((tt1, evals); α=0.1, linewidth=1, labels=false)
    xxs,yys,zzs,ρρs = [map(x -> x[i], evals) for i in 1:4];
    if labels
        plot(tt1, xxs, color=colors[2], alpha=α, linewidth=linewidth, label=L"$x$")
        plot(tt1, yys, color=colors[4], alpha=α, linewidth=linewidth, label=L"$y$")
        plot(tt1, zzs, color=colors[6], alpha=α, linewidth=linewidth, label=L"$z$")
        plot(tt1, ρρs, color=colors[8], alpha=α, linewidth=linewidth, label=L"Tr $ρ^2$")
    else
        plot(tt1, xxs, color=colors[2], alpha=α, linewidth=linewidth)
        plot(tt1, yys, color=colors[4], alpha=α, linewidth=linewidth)
        plot(tt1, zzs, color=colors[6], alpha=α, linewidth=linewidth)
        plot(tt1, ρρs, color=colors[8], alpha=α, linewidth=linewidth)
        
    end

end

# ╔═╡ abcd28b0-8d41-4953-9bda-800c56a9fa27
function plot_ensemble(sol_ens; α=0.1, linewidth=1, labels=false, average=false)
    close("all")
	tt1 = sol_ens[1]
    evs = collect(map(ρs -> expects.(ρs), sol_ens[2]));

    for i in 1:50
        plot_evals((tt1, evs[i]); α=α, labels=labels, linewidth=linewidth)
    end

    if average
        plot_evals((tt1, mean(evs)), α=1, linewidth=1.5, labels=true)
        title_string = "Trajectories w/ ensemble average"
    else
        plot_evals((tt1, evs[1]), α=1, linewidth=1.5, labels=true)
        title_string = "Trajectories"
    end
    
    ax = gca()
    ax.set_ylim([-1.1,1.1]) 
	legend()

    xlabel(L"$t$")
    ylabel("Bloch coordinates")
    title(title_string)
	
	gcf()
end

# ╔═╡ 36b2bbc2-55ac-45ad-b458-a792c0acddc1
green(text; title="Note") = Markdown.MD(Markdown.Admonition("correct", title, [text]))

# ╔═╡ b5c85a10-9f32-43d4-a4ea-5ec2b5dd37da
red(text; title="Note") = Markdown.MD(Markdown.Admonition("danger", title, [text]))

# ╔═╡ 7bfcd716-0b89-4ff0-85e5-0b7fd65a305d
tan(text; title="Note") = Markdown.MD(Markdown.Admonition("warning", title, [text]))

# ╔═╡ c44d4025-7ce8-4623-bc97-7bbb555914d1
blue(text; title="Note") = Markdown.MD(Markdown.Admonition("note", title, [text]))

# ╔═╡ 3cf94694-4002-4c1c-9593-e61ea5c99b7b
hint(text; title="Hint") = Markdown.MD(Markdown.Admonition("hint", title, [text]))

# ╔═╡ Cell order:
# ╠═7d43b95c-d3a4-11eb-16d3-47148906d472
# ╟─be7d340a-4ef1-4727-a531-7d6e5c568ad9
# ╟─b9bd68f2-df48-4072-8269-1898b7cf1b15
# ╟─e66e6723-16fe-4d8d-ae1f-f8237f8e8897
# ╠═601a8b56-2f41-4b2c-9ea2-6b26a52b4e26
# ╟─18c29abc-1f35-4b5a-bb27-a491c02cc98f
# ╟─01f57775-b647-4fea-8e96-0b8c8ceeff05
# ╠═eae605ed-f411-4f33-8066-bd8f01fc8a2d
# ╟─071641ac-a5a8-4d37-b531-88b1e083416c
# ╠═6294f81c-ea9d-4500-9751-b45f8a348639
# ╠═ee541c05-c187-4b43-a803-2255e254efe5
# ╠═c5de979e-de3b-4a20-9fc4-649851a311fa
# ╟─8aa08bfb-ac91-4e5f-9fb2-dbce02a38b8a
# ╠═ea558387-68dc-4dff-ae04-f7a877379510
# ╠═f2892559-e623-46d4-9ff2-f04fd9253734
# ╠═4487690f-349c-461c-a5b1-f9a6be9897bd
# ╠═a348c59c-acb5-428f-b544-d5e1038c56d6
# ╠═7bc8459e-1652-497f-87da-e2e45b145df7
# ╟─15613717-87c1-45d8-b194-1c7b1085d4f7
# ╠═840beeea-a731-4a38-b36b-e2ab848147d9
# ╠═bc8e8a5f-6476-402f-9cc6-0af722a6c033
# ╠═e7353f43-899d-4899-86b5-1303676c07ec
# ╠═2c1f165f-cabd-4599-afc7-348799754f5e
# ╠═90511e52-c343-490b-81b2-1502c546e16b
# ╠═2a842a56-4d47-4a4e-9c3c-7f7ce5c1ac75
# ╠═477756ac-f4ff-471e-b73a-13183585592a
# ╠═9e731ac9-e11c-427b-9569-4a925799042a
# ╟─ab459f70-e40a-499f-bd05-248270d7d4dd
# ╟─2487c71c-2d00-4ac0-8cec-7d3f2e80ab3e
# ╠═5d2e29b6-3f8e-4e66-9442-c325af182a6a
# ╟─8c422daa-fbf6-4a9e-83ce-3523948aa069
# ╟─3ad727fc-6955-4bb9-b9bf-cb34b63c8c70
# ╟─f9e49ee3-d8d2-4d92-8ce0-f9d07401d61c
# ╟─4ee435fb-e786-4a9f-8f46-1820e6ce9936
# ╠═d7d0ec12-242f-4c1c-875e-f960a2d9c88f
# ╠═fe17c0d6-7195-4491-84b6-582d907eeb92
# ╟─406317be-0a5a-4895-ad79-70fbc683f79b
# ╠═8f11676f-ab29-4a30-a23c-4e65039f3a44
# ╠═9a21aa2e-8034-41b0-a87b-d261efe70ee2
# ╟─295bef3d-8f73-453e-86dc-b8e61cdded83
# ╠═d429e614-a7d4-43e9-ac9b-d7dcfdfb4efa
# ╠═fcefa57c-d828-47b0-bc98-c107c9078946
# ╠═5957c035-687a-48d0-8c9a-ad0489c3f240
# ╠═d3ec6758-b27d-4c1e-a7de-7383f0f9027e
# ╠═75bc7334-aac0-4d5e-895d-894ba49ce492
# ╠═43de5e01-e2e7-4302-8914-0fc444645db4
# ╠═3258df22-de38-4ab5-92b7-0a828bc32155
# ╟─e59530ae-7f27-4836-a352-9d0a08d62451
# ╠═a589fb58-fc99-493d-bca4-fd21dc6a78d8
# ╠═c9a3a847-d910-4433-926a-d1bfd012f248
# ╠═0fa84d7f-59ff-4a27-821f-2114464129a6
# ╟─739984d0-91a6-45d1-bd2a-0561bfb3e4d5
# ╟─8f40ea32-6f75-4119-b90c-c72b7749360f
# ╠═69855c0b-c7e3-4e6c-bffa-7d5e47a55b1b
# ╟─b41b51dc-0763-4ad0-90e7-8e54fcc7d916
# ╟─6afc5215-2b51-425f-aecc-b20e15d5128a
# ╟─8bb2e4ea-e480-4e03-be07-24a5c4d6c780
# ╠═a1c5942e-4c27-4cbb-8fa8-b47d9449c523
# ╟─aec50388-2682-44cd-aab8-6474cc863f31
# ╟─c541a88e-0a5a-48d2-b490-91ef7500dbe5
# ╟─c0c7bcf3-3a8e-489d-9bda-4113764190ec
# ╟─abcd28b0-8d41-4953-9bda-800c56a9fa27
# ╟─16aed026-dcfd-432f-8028-32c7ef3afa5b
# ╟─36b2bbc2-55ac-45ad-b458-a792c0acddc1
# ╟─b5c85a10-9f32-43d4-a4ea-5ec2b5dd37da
# ╟─7bfcd716-0b89-4ff0-85e5-0b7fd65a305d
# ╟─c44d4025-7ce8-4623-bc97-7bbb555914d1
# ╟─3cf94694-4002-4c1c-9593-e61ea5c99b7b
