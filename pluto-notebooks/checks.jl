### A Pluto.jl notebook ###
# v0.15.1

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
	local dt = 1e-3  # integration time-step
	
	Random.seed!(1)
	sol_ens = ensemble(bayesian, (0,4τ), ρ0, H, [J], [C]; dt=dt, N=300)
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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
PyPlot = "~2.9.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Conda]]
deps = ["JSON", "VersionParsing"]
git-tree-sha1 = "299304989a5e6473d985212c28928899c74e9421"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.5.2"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "0fb723cd8c45858c22169b2e42269e53271a6df7"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.7"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "438d35d2d95ae2c5e8780b330592b6de8494e779"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.3"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "169bb8ea6b1b143c5cf57df6d34d022a7b60c6db"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.92.3"

[[PyPlot]]
deps = ["Colors", "LaTeXStrings", "PyCall", "Sockets", "Test", "VersionParsing"]
git-tree-sha1 = "67dde2482fe1a72ef62ed93f8c239f947638e5a2"
uuid = "d330b81b-6aea-500a-939a-2ce795aea3ee"
version = "2.9.0"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[VersionParsing]]
git-tree-sha1 = "80229be1f670524750d905f8fc8148e5a8c4537f"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.2.0"
"""

# ╔═╡ Cell order:
# ╠═75236fca-cfac-11eb-2ffe-c394a1c504cf
# ╟─337d1109-0183-4411-b269-34aeb03961b8
# ╟─eb760d9e-3e4c-4f57-b906-ff9147bd0e74
# ╟─7f920903-a4ff-4e12-8db2-f04db4175ea2
# ╟─8b4d5601-1310-4773-b0d5-699d86789827
# ╠═90d544e3-a262-4209-840a-cb796e113d92
# ╠═4e9af06b-8d7e-40fa-8d3e-fb7b12d86591
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
# ╟─5b49af65-c3f5-408c-9f8b-dd9a409a7ba7
# ╟─bb7c273b-aabe-439e-a858-82fd2d5471f5
# ╟─685f4b80-a9e5-4c9d-823a-819707f04601
# ╟─b1a25a71-ba30-4d92-b7fc-cfb54836234e
# ╟─d0204e86-97e9-4fc5-b124-ca84d6605f27
# ╟─c93a00a2-4c24-4bd0-9c82-ce8a506ed00f
# ╟─d6e770d5-462d-4c16-bd1e-0d0a140c9a3b
# ╟─eb9b4d46-97e3-40d7-b270-b9abf002f838
# ╟─93b538c3-dcdc-462d-b381-b7fc4075e703
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
