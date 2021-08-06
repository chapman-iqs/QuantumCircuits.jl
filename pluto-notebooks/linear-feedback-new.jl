### A Pluto.jl notebook ###
# v0.15.1

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
	C = [(σz, τm, η)]
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

# ╔═╡ 95ebf3b0-afbf-4dfe-869e-f0428fe49d63
md"""
Note `bayesian` outputs a record of the form
```math
R \equiv \frac{r}{\sqrt{\tau_m}} = \frac{z}{\sqrt{\tau_m}} + \frac{dW}{dt},
```
where $dW \sim \mathcal N(0,\sqrt{dt})$ is a Wiener increment. The hard-coded simulation that follows will take the normalized record $r$ as input, so this is calculated below.
"""

# ╔═╡ 304d5f8b-8870-46e9-a293-2df0d2c2895b
begin
	(tb, ρb, rs) = sol1
	rs = collect(rs[1]) # record output from bayesian
end

# ╔═╡ 0ef74ae9-2b4a-4fd4-abc1-e7bfc2eea931
md" ### Check against hard-coded numerical simulation"

# ╔═╡ 5022461e-6138-4ff6-9292-2820a8ec18d1
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
		Δ = -(Δ0 + Δ1*r) # the Rabi-drive depends on the readout (instantaneous here)
		
		yn2 = yn1*cos(dt*Δ) + zn1*sin(dt*Δ)
		zn2 = zn1*cos(dt*Δ) - yn1*sin(dt*Δ)
		
		if ideal 
			xn = xn1*exp(-dt*(1-η)/(2τm*η))
			yn = yn2*exp(-dt*(1-η)/(2τm*η))
			zn = zn2
		
		else 
			xn = xn1*exp(-dt/(2T1) - dt/T2 - dt*(1-η)/(2τm*η))
			yn = yn2*exp(-dt/(2T1) - dt/T2 - dt*(1-η)/(2τm*η))
			zn = zn2*exp(-dt/T1) - (1 - exp(-dt/T1)) end
		
		# store values
		push!(xs, xn)
		push!(ys, yn)
		push!(zs, zn)
		

	end
	
	return ts, [rs], (xs, ys, zs)
	
	
end

# ╔═╡ d614159d-c5d4-49cb-9320-115e2b9d59d3
(tt, r, bloch) = blochevolve(rs=rs)

# ╔═╡ edc0f4eb-755d-4ff3-bd70-634515d50e06
r[1]==rs

# ╔═╡ a74e31c5-a582-4240-a54f-cd09a3b0b720
fit(Normal, convert(Array{Float64}, dys[1]))

# ╔═╡ e2acd522-2add-4aff-ad49-30d7ab0c8292
σ = sqrt(τm/dt)

# ╔═╡ 65d5b466-263b-48d4-a0ce-7ecf32bcacbb
fit(Normal, convert(Array{Float64}, sol1[3][1]))

# ╔═╡ 9fa0d08b-d902-4850-ba87-546ab58a108f
plot_records(sol1; plot_title="record, general")

# ╔═╡ 7d41021a-5ee7-4f06-97c0-4efd490b88f9
dif = abs.(sol1[3] - dys)

# ╔═╡ 79604a4e-372e-4a89-b377-46ac6b284078
collect(sol1[3][1])

# ╔═╡ d37aaf9e-1d74-47ed-a0a4-0876316f6cda
dys[1]

# ╔═╡ 3258df22-de38-4ab5-92b7-0a828bc32155
md" ## Utilities "

# ╔═╡ e59530ae-7f27-4836-a352-9d0a08d62451
function rms(ser1::Array, ser2::Array)
	l = min(length(ser1), length(ser2))
	sqrt(sum((ser1[1:l] - ser2[1:l]).^2))/l
end

# ╔═╡ 0fa84d7f-59ff-4a27-821f-2114464129a6
expects = ρ -> collect(real(expect(ρ, s)) for s in [σx,σy,σz,ρ]) # ρ -> [<x>,<y>,<z>,<ρ>]

# ╔═╡ c9a3a847-d910-4433-926a-d1bfd012f248
function blochs(sol)
	(tt, ρt, _) = sol

	# Get Bloch components
	evs0 = expects.(ρt);
	xx,yy,zz,ρρ = [map(x -> x[i], evs0) for i in 1:4];

	(collect(tt), xx, yy, zz, ρρ)
	
end

# ╔═╡ d6734dcc-e5ee-4523-9970-b46542b1372a
plot_records((tt,blochs,dys); plot_title="record, hard-coded")

# ╔═╡ 0ed59062-6696-4a94-b85b-baa7c93df3a1
(xx,yy,zz) = blochs

# ╔═╡ 2b750028-1c81-4180-a59a-02bcbd331c7a
μ = mean(zz)

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

# ╔═╡ 6473130f-484e-47a4-8fac-8ebd5733e4a1
plot_solution(sol1; plot_title=title_str)

# ╔═╡ 604d9542-155f-4c29-a5ba-5997ac586ef5
plot_solution(sol1; plot_title=title_str)

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

# ╔═╡ 0d026ed1-10ed-4d5a-b897-6a3e91a78c1e
plot_blochs((tt,bloch); plot_title=title_str_hc)

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
function plot_timeseries(tt::Array, series...; plot_title="time series", xlab=L"$t$", ylab="arbitrary units", labels=[]::Array, colorpairs=false, kwargs...)
	close("all")
	ser_colors(i) = colorpairs ? colors[i] : colors[2i]
	
	label(i) = i > length(labels) ? i : labels[i]
    
    # Plot records vs. time
	ser = series[1]
	p = plot(tt, ser, color=ser_colors(1), label=label(1), kwargs...)
	ax = gca()
	
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

# ╔═╡ 369aa5bf-b76c-444f-aa48-4c7403b1e92c
plot_timeseries((tb, rs); plot_title="bayesian simulated record")

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

# ╔═╡ e8b31676-5bcf-4ee0-87c9-5b88cc095155

tan(md"I don't think the feedback code I wrote in `QuantumCircuits.jl` is working properly. I think when I am calling `dy` I am not using the form of dy that is already filled with the measured values...")


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
# ╠═ee541c05-c187-4b43-a803-2255e254efe5
# ╠═c5de979e-de3b-4a20-9fc4-649851a311fa
# ╟─8aa08bfb-ac91-4e5f-9fb2-dbce02a38b8a
# ╟─436c2d6e-eed8-4313-9c0c-295a9a706344
# ╟─40a10a59-80a1-4b21-b6e1-d1a774330dd0
# ╠═df33032d-29df-41dc-aedb-6dde844f5f01
# ╟─95ebf3b0-afbf-4dfe-869e-f0428fe49d63
# ╠═304d5f8b-8870-46e9-a293-2df0d2c2895b
# ╠═6473130f-484e-47a4-8fac-8ebd5733e4a1
# ╠═369aa5bf-b76c-444f-aa48-4c7403b1e92c
# ╟─0ef74ae9-2b4a-4fd4-abc1-e7bfc2eea931
# ╠═5022461e-6138-4ff6-9292-2820a8ec18d1
# ╠═edc0f4eb-755d-4ff3-bd70-634515d50e06
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
# ╟─e59530ae-7f27-4836-a352-9d0a08d62451
# ╟─c9a3a847-d910-4433-926a-d1bfd012f248
# ╟─0fa84d7f-59ff-4a27-821f-2114464129a6
# ╟─739984d0-91a6-45d1-bd2a-0561bfb3e4d5
# ╟─8f40ea32-6f75-4119-b90c-c72b7749360f
# ╟─69855c0b-c7e3-4e6c-bffa-7d5e47a55b1b
# ╟─b41b51dc-0763-4ad0-90e7-8e54fcc7d916
# ╟─6afc5215-2b51-425f-aecc-b20e15d5128a
# ╟─8bb2e4ea-e480-4e03-be07-24a5c4d6c780
# ╟─a1c5942e-4c27-4cbb-8fa8-b47d9449c523
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
