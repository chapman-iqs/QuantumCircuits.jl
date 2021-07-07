### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 75236fca-cfac-11eb-2ffe-c394a1c504cf
begin
	using Random
	using Statistics
	using PyPlot
	using QuantumCircuits 	# note: for some reason, QuantumCircuits has to be 
							# "used" last
end

# ╔═╡ 337d1109-0183-4411-b269-34aeb03961b8
md"""
# Ensemble-average convergence

The ensemble-average should be independent of $\eta$, and should converge to the master equation solution ($\eta = 0$). Below we consider the example of the reduced qubit as described in `reduced-qubit-evolution.jl`.

"""

# ╔═╡ eb760d9e-3e4c-4f57-b906-ff9147bd0e74
md" ##### Qubit Hilbert space"

# ╔═╡ 7f920903-a4ff-4e12-8db2-f04db4175ea2
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

# ╔═╡ 8b4d5601-1310-4773-b0d5-699d86789827
md" ##### Operators and definitions"

# ╔═╡ 4e9af06b-8d7e-40fa-8d3e-fb7b12d86591
begin
	ΩR  = 2π # Rabi frequency
	τ = 3.0 # Measurement collapse timescale
	Γ = 1/(2τ) # Measurement dephasing rate
end

# ╔═╡ c6666b04-6972-4aa2-99b0-8fc9f7f25e02
md" ##### Simulation"

# ╔═╡ 15360fec-3fba-48d2-a52c-844d3ae2f725
η = 1.

# ╔═╡ 90d544e3-a262-4209-840a-cb796e113d92
begin
	H = ΩR*σy/2
	J = √Γ*σz
	C = √(Γ*η)*σz
end

# ╔═╡ 1ec1382d-3355-4159-8ef6-a5101d6f6939
begin
	local ρ0 = dm(spindown(q)) # initial state
	local dt = 1e-2  # integration time-step
	
	Random.seed!(1)
	sol_ens = @time ensemble(bayesian, (0,4τ), ρ0, H, [J], [C]; dt=dt, N=150)
	sol_avg = (sol_ens[1],mean(sol_ens[2]))
end

# ╔═╡ 549c46e5-f1d2-463b-b355-ccde14fe4249
md"""
Try changing the value of `η` in the code above to see how it affects the ensemble average.
"""

# ╔═╡ 7c281ad8-963e-43e1-a885-68cb3ac71c62
begin
	local ρ0 = dm(spindown(q)) # initial state
	local dt = 1e-2  # integration time-step
	
	sol_η0 = bayesian((0, 4τ), ρ0, H, [J], []; dt=dt)
end

# ╔═╡ cb42e11d-85c6-425b-a48a-99087b695c5c
md" ## Time-step convergence"

# ╔═╡ a1b4f410-acee-402c-8e2c-1ccf5468e88a
md"""
Choosing the right time-step is critical to getting simulation results that make sense. In weakly measured systems, there is a natural finite timestep of the experiment due to the time-resolution of the detector. Ideally, you should use a simulation timestep that matches the experimental time-resolution. However, you must first ensure that the integration is robust at the desired time-step. Below are several checks in order of increasing complexity to ensure this is the case.

"""

# ╔═╡ 87ba3059-72fe-48a2-869b-3d5bac1f6328
md" ##### ensemble average convergence for different initial conditions "

# ╔═╡ c9db4dd8-ad8d-467f-9f14-a3b4f5c1a238
md"""
The ensemble average should be independent of the initial condition chosen. If it is not, a smaller time-step may be necessary.

"""

# ╔═╡ d4f0260d-3d43-4616-9257-9357caf3ac76
begin
	# initial conditions
	local ρ0 = dm(spindown(q)) 
	local ρ1 = dm(spinup(q))
	local u = (ρ0 + ρ1)/2
	
	dt = 1e-3

	bayes_sol1 = ensemble(bayesian, (0,4τ), ρ0, H, [J], [C]; dt=dt, N=500)
	bayes_sol2 = ensemble(bayesian, (0,4τ), u, H, [J], [C]; dt=dt, N=500);
end

# ╔═╡ 5b49af65-c3f5-408c-9f8b-dd9a409a7ba7
md" ##### single-trajectory convergence across $\Delta t$"

# ╔═╡ bb7c273b-aabe-439e-a858-82fd2d5471f5
md"""
A single trajectory (correpsonding to a given random seed / noise realization) should be robust to changes in $\Delta t$.
"""

# ╔═╡ 685f4b80-a9e5-4c9d-823a-819707f04601
md" ### Utilities"

# ╔═╡ b1a25a71-ba30-4d92-b7fc-cfb54836234e
# colorscheme
begin
	colorscheme = "Paired"
	cmap = plt.matplotlib.cm.get_cmap(colorscheme)
	colors=collect(map(x -> cmap(x), 0:11))
	md" `colorscheme`"
end

# ╔═╡ d0204e86-97e9-4fc5-b124-ca84d6605f27
expects = ρ -> collect(real(expect(ρ, s)) for s in [σx,σy,σz,ρ]) # ρ -> [<x>,<y>,<z>,<ρ>]

# ╔═╡ c93a00a2-4c24-4bd0-9c82-ce8a506ed00f
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

# ╔═╡ d6e770d5-462d-4c16-bd1e-0d0a140c9a3b
begin
	green(text; title="Note") = Markdown.MD(Markdown.Admonition("correct", title, [text]))
	
	red(text; title="Note") = Markdown.MD(Markdown.Admonition("danger", title, [text]))
	
	tan(text; title="Note") = Markdown.MD(Markdown.Admonition("warning", title, [text]))
	
	blue(text; title="Note") = Markdown.MD(Markdown.Admonition("note", title, [text]))
	
	hint(text; title="Hint") = Markdown.MD(Markdown.Admonition("hint", title, [text]))
	
end

# ╔═╡ d09a282e-4b9b-4718-89bf-62c42603f937
hint(md"If there seems to be a mismatch, try increasing $N$ (the number of trajectories) simulated in the ensemble. The larger $\eta$ is, the greater $N$ must be to see good convergence.")

# ╔═╡ 34a3b82f-7722-407a-b703-555ad2efb59c
green(md"If you're not getting ensemble-average convergence in your own simulations, this could be for a few reasons: (i) you need to simulate more trajectories to get convergence (especially if $\eta$ is close to 1)
	; (ii) you are having time-step convergence issues (see next section); or (iii) there is a mismatch between your measurement operators `C` and your dissipation operators `J`.
	"; title="Tip")

# ╔═╡ eb9b4d46-97e3-40d7-b270-b9abf002f838
function plot_ensemble(sol_ens; α=0.1, linewidth=1, labels=false, average=false, plot_title="Trajectories")
    close("all")
	tt1 = sol_ens[1]
    evs = collect(map(ρs -> expects.(ρs), sol_ens[2]));

    for i in 1:50
        plot_evals((tt1, evs[i]); α=α, labels=labels, linewidth=linewidth)
    end

    if average
        plot_evals((tt1, mean(evs)), α=1, linewidth=1.5, labels=true)
    else
        plot_evals((tt1, evs[1]), α=1, linewidth=2, labels=true)
    end
    
    ax = gca()
    ax.set_ylim([-1.1,1.1]) 
	legend()

    xlabel(L"$t$")
    ylabel("Bloch coordinates")
    title(plot_title)
	
	gcf()
end

# ╔═╡ dd7aba68-1d4d-4e38-a155-b0ff1dc5ea6f
 plot_ensemble(sol_ens; average=true, plot_title="Trajectories w/ ensemble average; η = $(η) ")

# ╔═╡ 93b538c3-dcdc-462d-b381-b7fc4075e703
# Plotting
function plot_solutions((sol1,sol2); plot_title="Rabi Oscillation", labels=(L"$\eta = 0$", "avg"))
    close("all")
	
	(l1,l2) = labels
    
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
	p = plot(tt1, x1, color=colors[2], linewidth=2, label=L"$x_l1$")
    plot(tt1, y1, color=colors[4], linewidth=2, label=L"$y_l1$")
    plot(tt1, z1, color=colors[6], linewidth=2, label=L"$z_l1$")
    plot(tt1, ρ1, color=colors[8], linewidth=2, label=L"(Tr $\rho^2)_l1$")
    plot(tt2, x2,  color=colors[1], linestyle="dashed", label=L"$x_l2$")
    plot(tt2, y2, color=colors[3], linestyle="dashed", label=L"$y_l2")
	plot(tt2, z2, color=colors[5], linewidth=2, label=L"$z_l2$")
    plot(tt2, ρ2, color=colors[7], linewidth=2, linestyle="dashed", label=L"(Tr $\rho^2)_l2$")
  
	ax = gca()
	ax.set_ylim([-1.1,1.1]) 
	xlabel(L"$t$")
	ylabel("Bloch coordinates")
	title(plot_title)
	legend()
	gcf()
end

# ╔═╡ 97e22fd8-c3ce-4f08-93b8-7c2f4ac5f38d
plot_solutions((sol_η0, sol_avg), plot_title="η = 0 and ensemble average for η = $(η) ")

# ╔═╡ 57e31f36-bf5a-45ea-b0bd-9a723ec879d8
begin
	bayes_sol1avg = (bayes_sol1[1], mean(bayes_sol1[2]))
	bayes_sol2avg = (bayes_sol2[1], mean(bayes_sol2[2]))
	
	plot_solutions((bayes_sol1avg, bayes_sol2avg), plot_title="title")
end

# ╔═╡ Cell order:
# ╠═75236fca-cfac-11eb-2ffe-c394a1c504cf
# ╟─337d1109-0183-4411-b269-34aeb03961b8
# ╟─eb760d9e-3e4c-4f57-b906-ff9147bd0e74
# ╟─7f920903-a4ff-4e12-8db2-f04db4175ea2
# ╟─8b4d5601-1310-4773-b0d5-699d86789827
# ╠═90d544e3-a262-4209-840a-cb796e113d92
# ╟─4e9af06b-8d7e-40fa-8d3e-fb7b12d86591
# ╟─c6666b04-6972-4aa2-99b0-8fc9f7f25e02
# ╠═1ec1382d-3355-4159-8ef6-a5101d6f6939
# ╠═dd7aba68-1d4d-4e38-a155-b0ff1dc5ea6f
# ╠═97e22fd8-c3ce-4f08-93b8-7c2f4ac5f38d
# ╠═15360fec-3fba-48d2-a52c-844d3ae2f725
# ╟─549c46e5-f1d2-463b-b355-ccde14fe4249
# ╟─d09a282e-4b9b-4718-89bf-62c42603f937
# ╟─7c281ad8-963e-43e1-a885-68cb3ac71c62
# ╟─34a3b82f-7722-407a-b703-555ad2efb59c
# ╟─cb42e11d-85c6-425b-a48a-99087b695c5c
# ╟─a1b4f410-acee-402c-8e2c-1ccf5468e88a
# ╟─87ba3059-72fe-48a2-869b-3d5bac1f6328
# ╟─c9db4dd8-ad8d-467f-9f14-a3b4f5c1a238
# ╠═d4f0260d-3d43-4616-9257-9357caf3ac76
# ╠═57e31f36-bf5a-45ea-b0bd-9a723ec879d8
# ╟─5b49af65-c3f5-408c-9f8b-dd9a409a7ba7
# ╟─bb7c273b-aabe-439e-a858-82fd2d5471f5
# ╟─685f4b80-a9e5-4c9d-823a-819707f04601
# ╟─b1a25a71-ba30-4d92-b7fc-cfb54836234e
# ╟─d0204e86-97e9-4fc5-b124-ca84d6605f27
# ╟─c93a00a2-4c24-4bd0-9c82-ce8a506ed00f
# ╟─d6e770d5-462d-4c16-bd1e-0d0a140c9a3b
# ╟─eb9b4d46-97e3-40d7-b270-b9abf002f838
# ╠═93b538c3-dcdc-462d-b381-b7fc4075e703
