### A Pluto.jl notebook ###
# v0.17.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 7d43b95c-d3a4-11eb-16d3-47148906d472
begin
	cd("/Users/sachagreenfield/Desktop/GitHub/QuantumCircuits.jl")
	import Pkg
	Pkg.activate(".")
	using Random
	using Statistics
	using PyPlot
	using QuantumCircuits
end

# ╔═╡ be7d340a-4ef1-4727-a531-7d6e5c568ad9
md"""
# Record normalization

In this interactive notebook, we do basic checks on records output by our simulations.

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

Consider basic two-level system evolution with Hamiltonian

$\hat H = \frac{\Omega_R}2 \hat \sigma_x,$
with measurement decay rate $\Gamma_m = 1/2 \tau_m$. If non-idealities are included, there is also energy decay with rate $\Gamma_1 = 1/2 T_1$ and environmental dephasing with rate $\Gamma_2 = 1/T_2$.
	
"""

# ╔═╡ 01f57775-b647-4fea-8e96-0b8c8ceeff05
md" ### Parameters and operators"

# ╔═╡ c43c0486-d3fb-4457-b7ff-1927ec105d53
begin
	ΩR  = 2π # Rabi frequency
	τm = 3.0 # Measurement collapse timescale
	Γm = 1/(2τm) # Measurement dephasing rate
	η = 1
	
	T1 = 10 # energy decay time
	T2 = 10 # environmental dephasing time

	# corresponding rates
	Γ1 = 1/(2T1)
	Γ2 = 1/T2
end

# ╔═╡ 09875de9-a56c-436f-a51e-fb5639d4f267
begin
	T = (0,4τm) # simulation duration
	(x0,y0,z0) = (0., 0.3, 0.91)
	ρ0 = DenseOperator(0.5*(id + x0*σx + y0*σy + z0*σz))
	dt = 1e-3 # integration time-step
end

# ╔═╡ f877ba39-0665-4267-b2dd-a84d90d41e35
expect(ρ0, ρ0)

# ╔═╡ 436c2d6e-eed8-4313-9c0c-295a9a706344
md" ## Simulation"

# ╔═╡ dd8df49a-e3b7-48c4-87c5-5e31b6141a3e
md"""
measurement type:  $(@bind meas_type html"<select><option value='heterodyne'>heterodyne</option><option value='homodyne'>homodyne</option></select>")
"""

# ╔═╡ d1184046-7dc1-45b7-898d-00e7a5837dd7
md"""
simulation type:  $(@bind sim_type html"<select><option value='bayesian'>bayesian</option><option value='rouchon'>rouchon</option></select>")
"""

# ╔═╡ de66a6c3-3180-4c54-a8da-5a9976f678e5
md" Simulate using the `QuantumCircuits.jl`'s $sim_type method."

# ╔═╡ c5de979e-de3b-4a20-9fc4-649851a311fa
md"""
**include non-idealities?** $(@bind nonideal html"<input type=checkbox >")
"""

# ╔═╡ eae605ed-f411-4f33-8066-bd8f01fc8a2d
begin
	H = ΩR*σx/2
	J = nonideal ? [√Γ1*σm, √(Γ2+Γm*(1-η))*σz] : [√((1-η)*Γm)*σz] 
	C = [√(Γm*η)*σz]
end

# ╔═╡ f28ef695-b27f-4c7c-8874-290d68e4e4ff
md"Now, let's look at the records output from $sim_type. They are stored in the third element of `sol1`:"

# ╔═╡ efcdb0b0-f054-4d04-a319-cd7fd6336183
md"Plotted as a time series,"

# ╔═╡ c9661f73-7e76-4b60-a5cb-e0233519c492
md"""
In the continuum limit, the record takes the form 
$r = z + \zeta$,
where $z = \langle \hat \sigma_z \rangle$ is the Bloch coordinate being measured, and
$\zeta \sim \mathcal{N}(0,\sqrt{\tau_m/dt})$ is a zero-mean Gaussian noise of standard deviation $\sqrt{\tau_m/dt}$. Note that `bayesian` and `rouchon` output records of different forms.
"""

# ╔═╡ 0ef74ae9-2b4a-4fd4-abc1-e7bfc2eea931
md" ### Check against hard-coded numerical simulation"

# ╔═╡ cfb94e57-06a7-4a52-b42e-7db995582670
md"""
The following code implements a hard-coded simulation of the system, to compare against the `QuantumCircuits.jl` implementation.
"""

# ╔═╡ 5022461e-6138-4ff6-9292-2820a8ec18d1
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
		
		yn2 = yn1*cos(dt*ΩR) - zn1*sin(dt*ΩR)
		zn2 = zn1*cos(dt*ΩR) + yn1*sin(dt*ΩR)
		
		if nonideal
			
			xn = xn1*exp(-dt/(2T1) - dt/T2 - dt*(1-η)/(2τm*η))
			yn = yn2*exp(-dt/(2T1) - dt/T2 - dt*(1-η)/(2τm*η))
			zn = zn2*exp(-dt/T1) - (1 - exp(-dt/T1))
			
		else
			
			xn = xn1*exp(-dt*(1-η)/(2τm*η))
			yn = yn2*exp(-dt*(1-η)/(2τm*η))
			zn = zn2
			
			
		end
		
		# store values
		push!(dy, r)
		push!(xs, xn)
		push!(ys, yn)
		push!(zs, zn)
		

	end
	
	return ts, [dy], (xs, ys, zs)
	
	
end

# ╔═╡ c06f5b1c-20e0-47a8-b80f-3d6f5353c880
(tt, dys, blochs) = simulate()

# ╔═╡ 865e5929-db0d-4eb9-8f40-9dcdb8bd6c30
σ = round(sqrt(τm/dt), digits=4)

# ╔═╡ dd5439f2-5572-4584-b39c-9f18e07a2165
sqrt(1/dt)

# ╔═╡ 3258df22-de38-4ab5-92b7-0a828bc32155
md" ## Utilities "

# ╔═╡ 8bdb2b60-b8b5-4793-a8d0-0a5d86879409
Statistics

# ╔═╡ 8aa08bfb-ac91-4e5f-9fb2-dbce02a38b8a
title_str = nonideal ? "non-ideal Rabi oscillation" : "ideal Rabi oscillation"

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

# ╔═╡ 9f08d06c-437e-47a9-9969-0297c12dbac8
plot_blochs((tt,blochs), plot_title=string(title_str, ", hard-coded"))

# ╔═╡ 8d8b700b-bb00-4279-b64e-eb06bd3cb986
function record_histograms(sols...; plot_title="record histogram", labels=[]::Array)
	close("all")
	
	μσs = []
	hist_colors = []
	hist_labels = []
	
	for i in 1:length(sols)
		label = i > length(labels) ? i : labels[i]
		sol = sols[i]

		(tt, _, dys) = sol
		ty = typeof(dys[1][end])
		comp = (ty == ComplexF64)
		dys = convert(Array{ty}, dys[1])
		
		# get mean and std dev for real part
		(μ, σ) = params(fit(Normal, real.(dys)))
		push!(μσs, map(p -> round(p, digits=4), (μ, σ)))

		# make histogram
		sublabel = comp ? string(label, " (real)") : label
		
		n, bins, patches = hist(real.(dys), 50, density=false, 					 									facecolor=colors[2i], alpha=1, label=sublabel)
		push!(hist_colors, colors[2i])
		
		
		if comp
			
			# get mean and std dev for imaginary part
			(μ, σ) = params(fit(Normal, imag.(dys)))
			push!(μσs, map(p -> round(p, digits=4), (μ, σ)))

			# make histogram
			n, bins, patches = hist(imag.(dys), 50, density=false, 					 									facecolor=colors[2i-1], alpha=0.75, 												label=string(label, " (imag)"))
			push!(hist_colors, colors[2i-1])

		else 
			
		end

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
	ylabel("frequency")
	title("record histograms")
	gcf()
		
end

# ╔═╡ 92affe90-b18d-4f7a-8da0-f708dc0f7bb8
# Plotting
function plot_records(sol; plot_title="record", labels=[]::Array)
	close("all")
	
	label(i) = i > length(labels) ? i : labels[i]
	labr(i) = string(label(i), " (real)")
	labi(i) = string(label(i), " (imag)")
	
	(tt, _, dys) = sol
	comp = (typeof(dys[1][end]) == ComplexF64)
    
    # Plot records vs. time
	if comp
		p = plot(tt, real.(dys[1]), color=colors[2], label=labr(1))
		plot(tt, imag.(dys[1]), color=colors[1], label=labi(1))
		ax = gca()
		
		for i in 2:length(dys)
			dy = dys[i]
			plot(tt, real.(dy), color=colors[2i],label=labr(i))
			plot(tt, imag.(dy), color=colors[2i-1],label=labi(i)) end
		
	else
		p = plot(tt, dys[1], color=colors[2], label=label(1))
		ax = gca()
		
		for i in 2:length(dys)
			dy = dys[i]
			plot(tt, dy, color=colors[2i],label=label(i)) end
		
	end
	
    xlabel(L"$t$")
    ylabel("arbitrary units")
    title(plot_title)
    legend()
    gcf()
end

# ╔═╡ 5e5a01f3-adc0-46fd-b54f-d6542c3766c7
plot_records((tt,blochs,dys); plot_title="record, hard-coded", labels=("bayesian"))

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

# ╔═╡ 863696ee-6787-455b-9891-86865ca0c681
begin 
	bayes = (sim_type == "bayesian")
	heterodyne = (meas_type == "heterodyne")
	
	if !bayes && heterodyne
		tan(md" `rouchon` currently does not implement heterodyne. A homodyne simulation will be run.") end
	
end

# ╔═╡ df33032d-29df-41dc-aedb-6dde844f5f01
begin
	Random.seed!(1)
	sol1 = bayes ? bayesian(T, ρ0, H, J, C; dt=dt, heterodyne=heterodyne) :
					rouchon(T, ρ0, H, J, C; dt=dt)
end

# ╔═╡ 6473130f-484e-47a4-8fac-8ebd5733e4a1
plot_solution(sol1; plot_title=string(title_str, ", $sim_type"))

# ╔═╡ 5fc6f4f2-202e-404d-be38-3d2626df1ed0
sol1[3]

# ╔═╡ 2f1d495c-a6b1-440c-acdc-ad89d36bc8f5
plot_records(sol1; plot_title="record, general", labels=[sim_type])

# ╔═╡ c49b0f3f-fe63-4086-bafb-2f81141e08a6
record_histograms(sol1; labels=[sim_type])

# ╔═╡ 1501fcfb-a944-4ef0-91c7-88f4fa63764e
plot_solution(sol1; plot_title=string(title_str, ", $sim_type"))

# ╔═╡ 1ee7f64b-7be6-4bd6-a187-c35ec0426096
record_histograms((tt,blochs,dys), sol1; labels=["hard-coded", "QuantumCircuits.jl"])

# ╔═╡ ec829366-a768-47e6-aeb9-753b2f05c5ed
if bayes && heterodyne
	md"""
	`bayesian` outputs a record of the form 
	
	$R \equiv \frac{r}{\sqrt{\tau_m}} = \frac{z}{\sqrt{\tau_m}} + \frac{dW}{dt},$
	
	where $dW \sim \mathcal N(0,\sqrt{dt})$ is a Wiener increment .
	
	Since we are doing a heterodyne measurement, $sim_type outputs one complex record. This can be understood as two independent real records, each with 
	
	$\sigma = 1/\sqrt{2 dt}$  
	
	= $(round(1/sqrt(2*dt), digits=4)) (**check this**)
	
	"""
	
elseif bayes && !heterodyne
	
	md"""
	`bayesian` outputs a record of the form 
	
	$R \equiv \frac{r}{\sqrt{\tau_m}} = \frac{z}{\sqrt{\tau_m}} + \frac{dW}{dt},$
	
	where $dW \sim \mathcal N(0,\sqrt{dt})$ is a Wiener increment .
	
	Since we are doing a homodyne measurement, `bayesian`
	outputs one real record of standard deviation
	
	$\sigma = 1/\sqrt{dt}$ 
	
	= $(round(1/sqrt(dt), digits=4)).
	
	"""
		
		
	
elseif !bayes && heterodyne
	md"""
	`rouchon` outputs a record of the form 
	
	$dy \equiv r \frac{dt}{\sqrt{\tau_m}} = z \frac{dt}{\sqrt{\tau_m}}+ dW,$
	
	where $dW \sim \mathcal N(0,\sqrt{dt})$ is a Wiener increment .
	
	Since we are doing a heterodyne measurement, `rouchon` outputs one complex record. This can be understood as two independent real records, each with $\sigma = \sqrt{dt/2}$ = $(round(sqrt(dt/2),digits=4)). (**check this**) **Currently rouchon does not implement heterodyne measurement.**
	
	"""
	
elseif !bayes && !heterodyne
	md"""
	`rouchon` outputs a record of the form 
	
	$dy \equiv r \frac{dt}{\sqrt{\tau_m}} = z \frac{dt}{\sqrt{\tau_m}}+ dW,$
	
	where $dW \sim \mathcal N(0,\sqrt{dt})$ is a Wiener increment .
	
	Since we are doing a homodyne measurement, `rouchon` outputs one real record of standard deviation 
	
	$\sigma = \sqrt{dt}$ 
	
	= $(round(sqrt(dt),digits=4)).
	
	"""
	
end

# ╔═╡ 9391f19c-8c1d-4366-b76a-f304fe65ac0d
blue(text; title="Note") = Markdown.MD(Markdown.Admonition("note", title, [text]))

# ╔═╡ 8387e17a-5672-4aed-9bee-130c8743375b
hint(text; title="Hint") = Markdown.MD(Markdown.Admonition("hint", title, [text]))

# ╔═╡ Cell order:
# ╠═7d43b95c-d3a4-11eb-16d3-47148906d472
# ╟─be7d340a-4ef1-4727-a531-7d6e5c568ad9
# ╟─e66e6723-16fe-4d8d-ae1f-f8237f8e8897
# ╠═601a8b56-2f41-4b2c-9ea2-6b26a52b4e26
# ╟─18c29abc-1f35-4b5a-bb27-a491c02cc98f
# ╟─01f57775-b647-4fea-8e96-0b8c8ceeff05
# ╠═c43c0486-d3fb-4457-b7ff-1927ec105d53
# ╠═eae605ed-f411-4f33-8066-bd8f01fc8a2d
# ╠═09875de9-a56c-436f-a51e-fb5639d4f267
# ╠═f877ba39-0665-4267-b2dd-a84d90d41e35
# ╟─436c2d6e-eed8-4313-9c0c-295a9a706344
# ╟─de66a6c3-3180-4c54-a8da-5a9976f678e5
# ╟─dd8df49a-e3b7-48c4-87c5-5e31b6141a3e
# ╠═d1184046-7dc1-45b7-898d-00e7a5837dd7
# ╟─c5de979e-de3b-4a20-9fc4-649851a311fa
# ╠═863696ee-6787-455b-9891-86865ca0c681
# ╠═df33032d-29df-41dc-aedb-6dde844f5f01
# ╠═6473130f-484e-47a4-8fac-8ebd5733e4a1
# ╟─f28ef695-b27f-4c7c-8874-290d68e4e4ff
# ╠═5fc6f4f2-202e-404d-be38-3d2626df1ed0
# ╟─efcdb0b0-f054-4d04-a319-cd7fd6336183
# ╠═2f1d495c-a6b1-440c-acdc-ad89d36bc8f5
# ╟─c9661f73-7e76-4b60-a5cb-e0233519c492
# ╠═ec829366-a768-47e6-aeb9-753b2f05c5ed
# ╠═c49b0f3f-fe63-4086-bafb-2f81141e08a6
# ╟─0ef74ae9-2b4a-4fd4-abc1-e7bfc2eea931
# ╟─cfb94e57-06a7-4a52-b42e-7db995582670
# ╠═5022461e-6138-4ff6-9292-2820a8ec18d1
# ╠═c06f5b1c-20e0-47a8-b80f-3d6f5353c880
# ╠═9f08d06c-437e-47a9-9969-0297c12dbac8
# ╠═1501fcfb-a944-4ef0-91c7-88f4fa63764e
# ╠═5e5a01f3-adc0-46fd-b54f-d6542c3766c7
# ╠═865e5929-db0d-4eb9-8f40-9dcdb8bd6c30
# ╠═dd5439f2-5572-4584-b39c-9f18e07a2165
# ╠═1ee7f64b-7be6-4bd6-a187-c35ec0426096
# ╟─3258df22-de38-4ab5-92b7-0a828bc32155
# ╠═8bdb2b60-b8b5-4793-a8d0-0a5d86879409
# ╟─8aa08bfb-ac91-4e5f-9fb2-dbce02a38b8a
# ╟─1e027b22-90e8-4b48-9eb9-2722e1aa612e
# ╟─1bd31808-2978-4168-a40b-831abd09b69a
# ╟─fbcefdcc-0db5-4676-8ca8-386501f6a790
# ╟─4f018bfb-c03d-410b-bec6-1b7070bb309f
# ╟─5779b365-b466-4fad-90b7-e47df73ea707
# ╠═8d8b700b-bb00-4279-b64e-eb06bd3cb986
# ╟─92affe90-b18d-4f7a-8da0-f708dc0f7bb8
# ╟─eef64e1d-bf71-4ef9-9e1a-d683cd7be679
# ╟─32aa09f0-3493-4360-814c-0c3928029c94
# ╟─daef5add-d4b1-4c68-b273-f5c386f511d0
# ╟─47e449dd-57b5-4bf1-936b-f56a132f541a
# ╟─2faba1da-cc7c-408b-9253-10c0185d978d
# ╟─9098967b-9e56-4062-9455-65bbaab3c5e4
# ╠═5cac0386-f12d-47a9-8e3e-8b0e1b131216
# ╟─9391f19c-8c1d-4366-b76a-f304fe65ac0d
# ╟─8387e17a-5672-4aed-9bee-130c8743375b
