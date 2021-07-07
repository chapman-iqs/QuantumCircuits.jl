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

$\hat H = \frac{\Omega_R}2 \sigma_y$.

"""

# ╔═╡ 01f57775-b647-4fea-8e96-0b8c8ceeff05
md" ### Parameters and operators"

# ╔═╡ c43c0486-d3fb-4457-b7ff-1927ec105d53
begin
	ΩR  = 2π # Rabi frequency
	τm = 3.0 # Measurement collapse timescale
	Γm = 1/(2τm) # Measurement dephasing rate
	η = 0.3
	
	T1 = 10 # energy decay time
	T2 = 10 # environmental dephasing time

	# corresponding rates
	Γ1 = 1/(2T1)
	Γ2 = 1/T2
end

# ╔═╡ c5de979e-de3b-4a20-9fc4-649851a311fa
ideal = true

# ╔═╡ eae605ed-f411-4f33-8066-bd8f01fc8a2d
begin
	H = ΩR*σx/2
	J = ideal ? [√Γm*σz] : [√Γ1*σm, √(Γ2+Γm)*σz]
	C = [√(Γm*η)*σz]
end

# ╔═╡ 8aa08bfb-ac91-4e5f-9fb2-dbce02a38b8a
title_str = ideal ? "ideal Rabi oscillation" : "non-ideal Rabi oscillation";

# ╔═╡ 09875de9-a56c-436f-a51e-fb5639d4f267
begin
	T = (0,4τm) # simulation duration
	(x0,y0,z0) = (0., 0.3, 0.91)
	ρ0 = DenseOperator(0.5*(id + x0*σx + y0*σy + z0*σz))
	dt = 1e-3  # integration time-step
end

# ╔═╡ 436c2d6e-eed8-4313-9c0c-295a9a706344
md" ## Simulation"

# ╔═╡ 40a10a59-80a1-4b21-b6e1-d1a774330dd0
md" ### `QuantumCircuits.jl` general implementation"

# ╔═╡ df33032d-29df-41dc-aedb-6dde844f5f01
begin
	Random.seed!(1)
	sol1 = bayesian(T, ρ0, H, J, C; dt=dt, heterodyne=false)
end

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
		
		yn2 = yn1*cos(dt*ΩR) + zn1*sin(dt*ΩR)
		zn2 = zn1*cos(dt*ΩR) - yn1*sin(dt*ΩR)
		
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

# ╔═╡ c06f5b1c-20e0-47a8-b80f-3d6f5353c880
(tt, dys, blochs) = simulate()

# ╔═╡ 92f9effa-86a9-49b5-9be7-f81a891b23a8
fit(Normal, convert(Array{Float64}, sol1[3][1]))

# ╔═╡ 865e5929-db0d-4eb9-8f40-9dcdb8bd6c30
σ = round(sqrt(τm/dt), digits=4)

# ╔═╡ 03049154-e003-4e32-9e3f-90e4d11dbcc8
md"""
We expect that the mean of the readout should be σ = $σ, but this does not match the results below. We should check the `master` branch of `QuantumCircuits.jl` to see if this is just an issue on `sacha-dev` branch (some mistake I made while trying to update the code to work with feedback), or if it is an older issue with our code.

Remember this depends on whether the measurement is heterodyne or not, which changes the variance of the distribution by a factor of $\sqrt{2}$; see method `readout`. The below simulations are homodyne measurements.
"""

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

# ╔═╡ 69f19624-075e-42b9-893a-19bc33e5cfbe
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

# ╔═╡ 9f08d06c-437e-47a9-9969-0297c12dbac8
plot_blochs((tt,blochs))

# ╔═╡ 8d8b700b-bb00-4279-b64e-eb06bd3cb986
function record_histograms(sols...; plot_title="record histogram", labels=[]::Array)
	close("all")
	
	μσs = []
	
	for i in 1:length(sols)
		label = i > length(labels) ? i : labels[i]
		sol = sols[i]

		(tt, _, dys) = sol
		ty = typeof(dys[1][end])
		dys = convert(Array{ty}, dys[1])

		# get mean and std dev
		(μ, σ) = params(fit(Normal, dys))
		push!(μσs, map(p -> round(p, digits=4), (μ, σ)))

		# make histogram
		n, bins, patches = hist(dys, 50, density=false, facecolor=colors[2i], 									alpha=0.75, label=label)

	end
	
	# write down (μ, σ) pairs as text boxes
	μσ_strings = map(μσ -> string("(μ, σ) = (", μσ[1], ", ", μσ[2], ")\n"), μσs)
	ax = gca()
	for i in 1:length(μσ_strings)
		str = μσ_strings[i]

		ax.text(0.05, 1 - 0.05i, str, transform=ax.transAxes, fontsize=10,
			verticalalignment="top", color=colors[2i])
		
	end
	
	
	legend()
	xlabel(L"$t$")
	ylabel("frequency")
	title("histograms")
	gcf()
end

# ╔═╡ 1ee7f64b-7be6-4bd6-a187-c35ec0426096
record_histograms((tt,blochs,dys), sol1; labels=["hard-coded", "QuantumCircuits.jl"])

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

# ╔═╡ 5e5a01f3-adc0-46fd-b54f-d6542c3766c7
plot_records((tt,blochs,dys); plot_title="record, hard-coded")

# ╔═╡ bb21dcf4-5c4c-4683-b60d-1c577888ff09
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
# ╟─e66e6723-16fe-4d8d-ae1f-f8237f8e8897
# ╠═601a8b56-2f41-4b2c-9ea2-6b26a52b4e26
# ╟─18c29abc-1f35-4b5a-bb27-a491c02cc98f
# ╟─01f57775-b647-4fea-8e96-0b8c8ceeff05
# ╠═c43c0486-d3fb-4457-b7ff-1927ec105d53
# ╠═eae605ed-f411-4f33-8066-bd8f01fc8a2d
# ╠═c5de979e-de3b-4a20-9fc4-649851a311fa
# ╠═8aa08bfb-ac91-4e5f-9fb2-dbce02a38b8a
# ╠═09875de9-a56c-436f-a51e-fb5639d4f267
# ╟─436c2d6e-eed8-4313-9c0c-295a9a706344
# ╟─40a10a59-80a1-4b21-b6e1-d1a774330dd0
# ╠═df33032d-29df-41dc-aedb-6dde844f5f01
# ╠═6473130f-484e-47a4-8fac-8ebd5733e4a1
# ╟─0ef74ae9-2b4a-4fd4-abc1-e7bfc2eea931
# ╟─cfb94e57-06a7-4a52-b42e-7db995582670
# ╟─5022461e-6138-4ff6-9292-2820a8ec18d1
# ╠═c06f5b1c-20e0-47a8-b80f-3d6f5353c880
# ╠═9f08d06c-437e-47a9-9969-0297c12dbac8
# ╠═69f19624-075e-42b9-893a-19bc33e5cfbe
# ╠═5e5a01f3-adc0-46fd-b54f-d6542c3766c7
# ╠═92f9effa-86a9-49b5-9be7-f81a891b23a8
# ╠═bb21dcf4-5c4c-4683-b60d-1c577888ff09
# ╠═865e5929-db0d-4eb9-8f40-9dcdb8bd6c30
# ╠═03049154-e003-4e32-9e3f-90e4d11dbcc8
# ╠═1ee7f64b-7be6-4bd6-a187-c35ec0426096
# ╟─e8b31676-5bcf-4ee0-87c9-5b88cc095155
# ╟─3258df22-de38-4ab5-92b7-0a828bc32155
# ╟─1e027b22-90e8-4b48-9eb9-2722e1aa612e
# ╟─1bd31808-2978-4168-a40b-831abd09b69a
# ╟─fbcefdcc-0db5-4676-8ca8-386501f6a790
# ╟─4f018bfb-c03d-410b-bec6-1b7070bb309f
# ╟─5779b365-b466-4fad-90b7-e47df73ea707
# ╠═8d8b700b-bb00-4279-b64e-eb06bd3cb986
# ╟─92affe90-b18d-4f7a-8da0-f708dc0f7bb8
# ╠═eef64e1d-bf71-4ef9-9e1a-d683cd7be679
# ╟─32aa09f0-3493-4360-814c-0c3928029c94
# ╟─daef5add-d4b1-4c68-b273-f5c386f511d0
# ╟─47e449dd-57b5-4bf1-936b-f56a132f541a
# ╟─2faba1da-cc7c-408b-9253-10c0185d978d
# ╟─9098967b-9e56-4062-9455-65bbaab3c5e4
# ╠═5cac0386-f12d-47a9-8e3e-8b0e1b131216
# ╟─9391f19c-8c1d-4366-b76a-f304fe65ac0d
# ╟─8387e17a-5672-4aed-9bee-130c8743375b
# ╟─7c9bcaed-d611-4994-96b4-71ea86867808
