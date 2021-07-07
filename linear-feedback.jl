### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 7d43b95c-d3a4-11eb-16d3-47148906d472
begin
	using Random
	using Statistics
	using PyPlot
	using Distributions
	using QuantumCircuits 	# note: for some reason, QuantumCircuits has to be 
							# "used" last
end

# ╔═╡ be7d340a-4ef1-4727-a531-7d6e5c568ad9
md" # Linear feedback stabilization"

# ╔═╡ b9bd68f2-df48-4072-8269-1898b7cf1b15
md"""
In this interactive notebook, we implement the linear feedback stabilization described in [1].
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

# ╔═╡ c5de979e-de3b-4a20-9fc4-649851a311fa
ideal = true

# ╔═╡ eae605ed-f411-4f33-8066-bd8f01fc8a2d
begin
	# all times given in μs
	# initial state
	(x0,y0,z0) = (0., 0.3, 0.91)
	ρ0 = DenseOperator(0.5*(id + x0*σx + y0*σy + z0*σz))
	
	dt = 0.5e-4  # integration time-step
	td = 0 # 200e-3 # delay time
	η = 0.41 # measurement efficiency
	
	θs = 3π/10 # target angle on Bloch sphere
	
	τm = 0.2 # measurement time
	Γm = 1/(2τm*η) # measurement rate
	T = (0, 4τm) # simulation duration
	
	Rs = 0.64 # radius of target
	
	# feedback drive parameters
	Δ0 = ideal ? -sin(2θs)/(4τm) : -sin(2θs)/(4τm*Rs^2) 
	Δ1 = 0 # ideal ? sin(θs)/τm : sin(θs)/(τm*Rs)
		
	T1 = 40 # energy decay time
	T2 = 60 # environmental dephasing time

	# corresponding rates
	Γ1 = 1/(2T1)
	Γ2 = 1/T2

	
end

# ╔═╡ ee541c05-c187-4b43-a803-2255e254efe5
begin
	# Hamiltonian defined at a certain time t, 
	# takes in readout r in some observable with time delay td
	H0(t, r::Array) = (Δ0 + Δ1*r[1])*σx/2 
	C = [√Γm*σz]
	J = ideal ? [√Γm*σz] : [√Γ1*σm, √(Γ2+Γm)*σz]
end

# ╔═╡ 8aa08bfb-ac91-4e5f-9fb2-dbce02a38b8a
begin
	
	title_str = ideal ? "ideal Rabi oscillation with feedback, general" : "non-ideal Rabi oscillation with feedback, general";
	title_str_hc = ideal ? "ideal Rabi oscillation with feedback, hard-coded" : "non-ideal Rabi oscillation with feedback, hard-coded";
	
end

# ╔═╡ 436c2d6e-eed8-4313-9c0c-295a9a706344
md" ## Simulation"

# ╔═╡ 40a10a59-80a1-4b21-b6e1-d1a774330dd0
md" ### `QuantumCircuits.jl` general implementation"

# ╔═╡ df33032d-29df-41dc-aedb-6dde844f5f01
begin
	Random.seed!(1)
	sol1 = bayesian(T, ρ0, H0, J, C; dt=dt, td=td, heterodyne=false)
end

# ╔═╡ 0ef74ae9-2b4a-4fd4-abc1-e7bfc2eea931
md" ### Check against hard-coded numerical simulation"

# ╔═╡ 5022461e-6138-4ff6-9292-2820a8ec18d1
# using a `let` block (instead of begin) allows for the necessary
# scope for variables like xn, yn, zn

function simulate()
	# readout distribution
	Random.seed!(1)
	dist = Normal(0, sqrt(τm/dt))
	
	ts = range(first(T), last(T), step=dt)
	(dy, xs, ys, zs) = map(i -> [], 1:4)
	
	# first time step
	xn = x0
	yn = y0
	zn = z0
	push!(dy, 0)
	push!(xs, xn)
	push!(ys, yn)
	push!(zs, zn)

	
	for t in ts[2:end]
		
		# sample from distribution
		r = zn + rand(dist)
		
 		# update equations
		pn = cosh(r*dt/τm) + zn*sinh(r*dt/τm)
		xn1 = xn/pn
		yn1 = yn/pn
		zn1 = (zn*cosh(r*dt/τm) + sinh(r*dt/τm))/pn
		Δ = Δ0 + Δ1*r # the Rabi-drive depends on the readout (instantaneous here)
		
		yn2 = yn1*cos(dt*Δ) + zn1*sin(dt*Δ)
		zn2 = zn1*cos(dt*Δ) - yn1*sin(dt*Δ)
		
		if ideal
			
			xn = xn1*exp(-dt*(1-η)/(2τm*η))
			yn = yn2*exp(-dt*(1-η)/(2τm*η))
			zn = zn2
			
		else
			
			xn = xn1*exp(-dt/(2T1) - dt/T2 - dt*(1-η)/(2τm*η))
			yn = yn2*exp(-dt/(2T1) - dt/T2 - dt*(1-η)/(2τm*η))
			zn = zn2*exp(-dt/T1) - (1 - exp(-dt/T1))
			
		end
		
		# store values
		push!(dy, r)
		push!(xs, xn)
		push!(ys, yn)
		push!(zs, zn)
		

	end
	
	return ts, [dy], (xs, ys, zs)
	
	
end

# ╔═╡ d614159d-c5d4-49cb-9320-115e2b9d59d3
(tt, dys, blochs) = simulate()

# ╔═╡ a74e31c5-a582-4240-a54f-cd09a3b0b720
fit(Normal, convert(Array{Float64}, dys[1]))

# ╔═╡ 0ed59062-6696-4a94-b85b-baa7c93df3a1
(xx,yy,zz) = blochs

# ╔═╡ 2b750028-1c81-4180-a59a-02bcbd331c7a
μ = mean(zz)

# ╔═╡ e2acd522-2add-4aff-ad49-30d7ab0c8292
σ = sqrt(τm/dt)

# ╔═╡ 65d5b466-263b-48d4-a0ce-7ecf32bcacbb
fit(Normal, convert(Array{Float64}, sol1[3][1]))

# ╔═╡ 7d41021a-5ee7-4f06-97c0-4efd490b88f9
dif = abs.(sol1[3] - dys)

# ╔═╡ 79604a4e-372e-4a89-b377-46ac6b284078
collect(sol1[3][1])

# ╔═╡ d37aaf9e-1d74-47ed-a0a4-0876316f6cda
dys[1]

# ╔═╡ 3258df22-de38-4ab5-92b7-0a828bc32155
md" ## Utilities "

# ╔═╡ 1e027b22-90e8-4b48-9eb9-2722e1aa612e
expects = ρ -> collect(real(expect(ρ, s)) for s in [σx,σy,σz,ρ]) # ρ -> [<x>,<y>,<z>,<ρ>]

# ╔═╡ 1bd31808-2978-4168-a40b-831abd09b69a
purity(x,y,z) = 0.5*(1 + x^2 + y^2 + z^2)

# ╔═╡ fbcefdcc-0db5-4676-8ca8-386501f6a790
R(x,y,z) = sqrt(x^2 + y^2 + z^2)

# ╔═╡ 47e449dd-57b5-4bf1-936b-f56a132f541a
# colorscheme
begin
	colorscheme = "Paired"
	cmap = plt.matplotlib.cm.get_cmap(colorscheme)
	colors=collect(map(x -> cmap(x), 0:11))
	md" `colorscheme`"
end

# ╔═╡ 4f018bfb-c03d-410b-bec6-1b7070bb309f
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

# ╔═╡ 6473130f-484e-47a4-8fac-8ebd5733e4a1
plot_solution(sol1; plot_title=title_str)

# ╔═╡ 604d9542-155f-4c29-a5ba-5997ac586ef5
plot_solution(sol1; plot_title=title_str)

# ╔═╡ 5779b365-b466-4fad-90b7-e47df73ea707
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

# ╔═╡ 0d026ed1-10ed-4d5a-b897-6a3e91a78c1e
plot_blochs((tt,blochs); plot_title=title_str_hc)

# ╔═╡ 92affe90-b18d-4f7a-8da0-f708dc0f7bb8
# Plotting
function plot_records(sol; plot_title="record")
	close("all")
	
	(tt, _, dys) = sol
    
    # Plot records vs. time
	p = plot(tt, dys[1], color=colors[2], label=L"record $1$")
	ax = gca()
	
	for i in 2:length(dys)
		dy = dys[i]
    	plot(tt, dy, color=colors[2i],label=L"$record i$")
	end
	
    xlabel(L"$t$")
    ylabel("arbitrary units")
    title(plot_title)
    legend()
    gcf()
end

# ╔═╡ d6734dcc-e5ee-4523-9970-b46542b1372a
plot_records((tt,blochs,dys); plot_title="record, hard-coded")

# ╔═╡ 9fa0d08b-d902-4850-ba87-546ab58a108f
plot_records(sol1; plot_title="record, general")

# ╔═╡ eef64e1d-bf71-4ef9-9e1a-d683cd7be679
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

# ╔═╡ 32aa09f0-3493-4360-814c-0c3928029c94
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

# ╔═╡ daef5add-d4b1-4c68-b273-f5c386f511d0
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

# ╔═╡ 2faba1da-cc7c-408b-9253-10c0185d978d
green(text; title="Note") = Markdown.MD(Markdown.Admonition("correct", title, [text]))

# ╔═╡ 9098967b-9e56-4062-9455-65bbaab3c5e4
red(text; title="Note") = Markdown.MD(Markdown.Admonition("danger", title, [text]))

# ╔═╡ 5cac0386-f12d-47a9-8e3e-8b0e1b131216
tan(text; title="Note") = Markdown.MD(Markdown.Admonition("warning", title, [text]))

# ╔═╡ e8b31676-5bcf-4ee0-87c9-5b88cc095155

tan(md"I don't think the feedback code I wrote in `QuantumCircuits.jl` is working properly. I think when I am calling `dy` I am not using the form of dy that is already filled with the measured values...")


# ╔═╡ 9391f19c-8c1d-4366-b76a-f304fe65ac0d
blue(text; title="Note") = Markdown.MD(Markdown.Admonition("note", title, [text]))

# ╔═╡ 8387e17a-5672-4aed-9bee-130c8743375b
hint(text; title="Hint") = Markdown.MD(Markdown.Admonition("hint", title, [text]))

# ╔═╡ 7c9bcaed-d611-4994-96b4-71ea86867808
md"""
!!! warning "Edit"

	We need to use simple master equation API rather than `bayesian` here. I'm not sure this currently exists in `QuantumCircuits.jl`, but should probably be the `jump-no-jump` method. Or does `bayesian` automatically do that in the absence of measurement?
"""

# ╔═╡ Cell order:
# ╠═7d43b95c-d3a4-11eb-16d3-47148906d472
# ╟─be7d340a-4ef1-4727-a531-7d6e5c568ad9
# ╟─b9bd68f2-df48-4072-8269-1898b7cf1b15
# ╟─e66e6723-16fe-4d8d-ae1f-f8237f8e8897
# ╠═601a8b56-2f41-4b2c-9ea2-6b26a52b4e26
# ╟─18c29abc-1f35-4b5a-bb27-a491c02cc98f
# ╟─01f57775-b647-4fea-8e96-0b8c8ceeff05
# ╠═eae605ed-f411-4f33-8066-bd8f01fc8a2d
# ╠═ee541c05-c187-4b43-a803-2255e254efe5
# ╠═c5de979e-de3b-4a20-9fc4-649851a311fa
# ╟─8aa08bfb-ac91-4e5f-9fb2-dbce02a38b8a
# ╟─436c2d6e-eed8-4313-9c0c-295a9a706344
# ╟─40a10a59-80a1-4b21-b6e1-d1a774330dd0
# ╠═df33032d-29df-41dc-aedb-6dde844f5f01
# ╠═6473130f-484e-47a4-8fac-8ebd5733e4a1
# ╟─0ef74ae9-2b4a-4fd4-abc1-e7bfc2eea931
# ╠═5022461e-6138-4ff6-9292-2820a8ec18d1
# ╠═d614159d-c5d4-49cb-9320-115e2b9d59d3
# ╠═0d026ed1-10ed-4d5a-b897-6a3e91a78c1e
# ╠═604d9542-155f-4c29-a5ba-5997ac586ef5
# ╠═d6734dcc-e5ee-4523-9970-b46542b1372a
# ╠═a74e31c5-a582-4240-a54f-cd09a3b0b720
# ╠═0ed59062-6696-4a94-b85b-baa7c93df3a1
# ╠═2b750028-1c81-4180-a59a-02bcbd331c7a
# ╠═e2acd522-2add-4aff-ad49-30d7ab0c8292
# ╠═65d5b466-263b-48d4-a0ce-7ecf32bcacbb
# ╠═9fa0d08b-d902-4850-ba87-546ab58a108f
# ╠═7d41021a-5ee7-4f06-97c0-4efd490b88f9
# ╠═79604a4e-372e-4a89-b377-46ac6b284078
# ╠═d37aaf9e-1d74-47ed-a0a4-0876316f6cda
# ╟─e8b31676-5bcf-4ee0-87c9-5b88cc095155
# ╟─3258df22-de38-4ab5-92b7-0a828bc32155
# ╟─1e027b22-90e8-4b48-9eb9-2722e1aa612e
# ╟─1bd31808-2978-4168-a40b-831abd09b69a
# ╟─fbcefdcc-0db5-4676-8ca8-386501f6a790
# ╠═4f018bfb-c03d-410b-bec6-1b7070bb309f
# ╠═5779b365-b466-4fad-90b7-e47df73ea707
# ╠═92affe90-b18d-4f7a-8da0-f708dc0f7bb8
# ╟─eef64e1d-bf71-4ef9-9e1a-d683cd7be679
# ╟─32aa09f0-3493-4360-814c-0c3928029c94
# ╟─daef5add-d4b1-4c68-b273-f5c386f511d0
# ╟─47e449dd-57b5-4bf1-936b-f56a132f541a
# ╟─2faba1da-cc7c-408b-9253-10c0185d978d
# ╟─9098967b-9e56-4062-9455-65bbaab3c5e4
# ╠═5cac0386-f12d-47a9-8e3e-8b0e1b131216
# ╟─9391f19c-8c1d-4366-b76a-f304fe65ac0d
# ╟─8387e17a-5672-4aed-9bee-130c8743375b
# ╟─7c9bcaed-d611-4994-96b4-71ea86867808
