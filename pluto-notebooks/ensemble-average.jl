### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 4c45fe3e-cd69-11eb-20cd-7bfb98c040cf
begin
	using Random
	using Statistics
	using PyPlot
	using QuantumCircuits 	# note: for some reason, QuantumCircuits has to be 
							# "used" last
end

# ╔═╡ 377a3336-20bd-4baa-a033-af8bbc8668a8
md"""
# Reduced qubit evolution

In this interactive notebook, we'll look at reduced qubit evolution as a basic example to gain intuition for Cavity Quantum Electrodynamics (CQED).
"""

# ╔═╡ 7f1176b8-f6d3-4fd1-a23e-31710dcfff10
md"""
## Problem setup

We begin by defining a two-level Hilbert space for the system. `QuantumCircuits.jl` uses `QuantumOptics.jl` as its backend: `SpinBasis`, `sigmax`, `identityoperator` and so forth are `QuantumOptics.jl` functions.
"""

# ╔═╡ 982565bc-b137-4c10-8e13-5fe37be86823
md" ##### Qubit Hilbert space"

# ╔═╡ 4284173a-be05-4b58-a8d9-7189301344fd
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

# ╔═╡ 85206c80-af15-4e1b-81e9-04bcb1083063
md" ##### Description of the evolution"

# ╔═╡ 629669ac-006a-4774-ae3e-67787bd0fa0a
md"""
We work in the rotating frame, so that $\omega_q = 0$ effectively. Thus the Hamiltonian is

$H_0 = \frac{\Omega_R}{2} \sigma_y.$
"""

# ╔═╡ bbbc92cd-8daa-4ed8-8bf0-d10461af4edd
md"""
The qubit is weakly and continuously monitored in its informational basis, $|+z \rangle$ (excited) and $|-z \rangle$ (ground). This ''weak measurement'' manifests in two primary ways:

* The information *collected* from the measurement creates **measurement backaction**, characterized by measurement collapse timescale $\tau$. This is represented via the dissipation operator $J_0 = \sqrt{\Gamma} \sigma_z$.

* The information *not collected* leads to **dephasing** at rate $\Gamma = 1/(2\tau)$. This is represented via the measurement collapse operator $C_0 = \sqrt{\Gamma \eta}\sigma_z$.

The balance of these two effects is determined by the signal collection efficiency $\eta$, where $\eta = 1$ means all information is collected, and $\eta = 0$ means no information is collected.
"""

# ╔═╡ 5b48435b-da1b-4547-ac7c-72c95e1d07ff
md"""
##### System parameters

The default parameters of our qubit system are chosen as
"""

# ╔═╡ 99f2b351-877a-43c8-a508-d5688f92ae0c
begin
	ΩR0  = 2π # Rabi frequency
	τ0 = 3.0 # Measurement collapse timescale
	Γ0 = 1/(2τ0) # Measurement dephasing rate
	η0 = 0.3
end

# ╔═╡ 259ed22f-b0c0-4fbe-9c6b-97a8730b6b14
H0 = ΩR0*σy/2

# ╔═╡ b290f188-47f9-4af1-a5f0-11d99c29c1e3
begin
	J0 = √Γ0*σz
	C0 = √(Γ0*η0)*σz;
end

# ╔═╡ 677b683b-ba7c-4b01-be5b-40ef97a075cd
md" ##Measurement-free evolution"

# ╔═╡ 1aa85cc1-32d2-4fa9-a29d-b042d4d4b94a
md"""
## Quantum trajectories

A *quantum trajectory* is a particular realization of the qubit evolution, based on the measurement record. The Hamiltonian evolution of the qubit is interleaved with Bayesian measurement updates based on information acquired from measurement. The resulting trajectory is a noisy perturbation of the ensemble-average behavior.

In the following code, we use `QuantumCircuits.jl`'s `bayesian` function to simulate trajectories for the reduced qubit system.

!!! warning "Edit"
	Include description of what a quantum trajectory is and the Bayesian filtering method for simulating one.

"""

# ╔═╡ 139565f9-0219-4939-a470-4252397dd6a0
begin
	ρ0 = dm(spindown(q)) # initial state
	dt = 1e-3  # integration time-step
	
	Random.seed!(1)
	sol1 = bayesian((0, 4τ0), ρ0, H0, [J0], [C0]; dt=dt)
end

# ╔═╡ 62170241-f050-4634-ad6a-8842ac96096c
md" This is one of infinitely many possible trajectories for the system. Let's generate an ensemble of 500 trajectories using `ensemble`:"

# ╔═╡ d1544d60-7de6-4f54-bffb-9ea0c3ddddc1
begin
	Random.seed!(1)
	sol_ens = @time ensemble(rouchon, (0,4τ0), ρ0, H0, [J0], [C0]; dt=dt, N=500)
end

# ╔═╡ 4d42ca32-5348-48b5-8a24-eafa5b303690
md" Let's look at 50 of these trajectories, with our original highlighted: "

# ╔═╡ ee6b82a5-52b0-4467-b73a-d6534de57202
md" While there is some randomness to the trajectories, there is clearly an average behavior that is consistent with all of them. By averaging the Bloch coordinates of all 500 trajectories, for every time step, we can obtain the average qubit evolution: "

# ╔═╡ 0b2ba4be-0213-4b3b-9943-291381c2b4fc
md" η0 = $η0"

# ╔═╡ 0cb6b898-347c-4baf-988f-d1dfac1d6557
md" This is equivalent to setting $\eta = 0$, and reproduces the master equation solution. Intuitively, when no signal is collected, our ''best estimate'' of the quantum state is the Hamiltonian evolution plus additional decoherence due to the uncollected measurement. You can gain intuition for this effect in the next section."

# ╔═╡ 7961752d-faee-4799-99b1-3bd3d17d1f31
md"""
## Measurement efficiency $\eta$ and measurement collapse timescale $\tau$

Let's gain intuition for how $\eta$ and $\tau$ affect our simulation. Play with the values of $\eta$ and $\tau$ below to see what happens.
"""

# ╔═╡ dff8aec2-7216-49e9-a19d-81b7ff79eb4d
md"""
##### Parameters:

"""

# ╔═╡ 6391d48d-4332-4209-bfe1-bec2e60e30c1
md"""

ΩR = 0
$(@bind ΩR_coeff html"<input type=range min=0.1 max=6 step=0.1 value=2>")
ΩR = 10π

η = 0
$(@bind η html"<input type=range min=0. max=1. step=0.02 value=0.3>")
η = 1 | - - - - - - | τ = 0.5 
$(@bind τ html"<input type=range min=0.5 max=10 step=0.5 value=3.0>")
τ=10.0

"""


# ╔═╡ 8638af64-25f9-4793-bc8b-5d22c3ef70fa
ΩR = ΩR_coeff*π

# ╔═╡ 02ce38f3-f5e1-4fea-9e6e-b887baa548e7
md" ##### Simulate"

# ╔═╡ 5e8b047e-b70f-4478-a8b2-a2e48b4a56af
begin
	# other parameters
	Γ = 1/(2τ)
	
	# operators
	H = ΩR*σy/2
	J = √Γ*σz
	C = √(Γ*η)*σz
	md" `operators defined`"
end

# ╔═╡ 2db9b4c4-0256-4e29-990c-2ca3ee9c7ddf
md"""
ΩR = $(round(ΩR/π, digits=3)) π | - - - - - - | η = $η | - - - - - - | τ = $τ | - - - - - - | Γ = 1/(2τ) = $(round(Γ,digits=3))
"""

# ╔═╡ 09b609f6-a441-4ae2-aaff-66c1b8b6ffae
begin
	Random.seed!(1)
	sol2 = @time bayesian((0, 4τ0), ρ0, H, [J], [C]; dt=dt, heterodyne=true)
	print()
end

# ╔═╡ 7e62e6c6-aa3c-4350-901c-15d017b8db42
md" ## Utilities "

# ╔═╡ b2235afc-1bc9-4cf4-ae21-6f535eb37e96
expects = ρ -> collect(real(expect(ρ, s)) for s in [σx,σy,σz,ρ]) # ρ -> [<x>,<y>,<z>,<ρ>]

# ╔═╡ 0a2219d8-2fbd-4b03-867f-966c4a595e78
# colorscheme
begin
	colorscheme = "Paired"
	cmap = plt.matplotlib.cm.get_cmap(colorscheme)
	colors=collect(map(x -> cmap(x), 0:11))
	md" `colorscheme`"
end

# ╔═╡ 31371407-5798-49a4-a83a-2c068d622c4b
# Plotting
function plot_solution(sol; plot_title="Rabi Oscillation")
	
	close("all")
    
    tt = sol[1]
    ρt = sol[2]
    
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

# ╔═╡ 7218d350-0b14-4fed-b87e-65a9eb573780
plot_solution(sol1; plot_title="Monitored Rabi Oscillation")

# ╔═╡ cdbd897e-d689-42fa-9977-2b3169b29d08
begin
	title_string = "bayesian heterodyne; τ = $(τ), η = $(η)"
	plot_solution(sol2; plot_title=title_string)
end

# ╔═╡ 99e90182-1ee7-45ef-8656-cc5a1e01d564
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

# ╔═╡ 53f32791-a721-4d80-96ce-fae9f533dbaa
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

# ╔═╡ 7d4f4a02-f3b8-4369-aad3-03ff4bc7ec0b
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

# ╔═╡ 2b5ae77e-07c0-4750-85dd-138dcf84da1f
 plot_ensemble(sol_ens)

# ╔═╡ 8696fe89-1c19-4852-afe0-2b642f5a9382
 plot_ensemble(sol_ens; average=true)

# ╔═╡ a255f7a5-d617-4444-80a4-4e4e7762b236
green(text; title="Note") = Markdown.MD(Markdown.Admonition("correct", title, [text]))

# ╔═╡ 6e86b5ee-5b38-4245-b555-c208151c0c53
red(text; title="Note") = Markdown.MD(Markdown.Admonition("danger", title, [text]))

# ╔═╡ 4585fc96-cc13-476d-99d6-de839c4ec691
if τ > τ0
	green(md"""As the measurement timescale τ increases, the dephasing rate $\Gamma = 1/(2\tau)$ decreases. Thus the system takes longer to dephase (remains high purity).""", title="τ increased from 3.0")
	
elseif τ < τ0
	red(md"""As the measurement timescale τ decreases, the dephasing rate $\Gamma = 1/(2\tau)$ increases. Thus the system dephases quickly and the purity decreases.""", title="τ decreased from 3.0")
		
end

# ╔═╡ fb32f241-34bb-413b-99d2-c0155b363670
tan(text; title="Note") = Markdown.MD(Markdown.Admonition("warning", title, [text]))

# ╔═╡ e730e11b-6fa5-4121-80eb-69a59d74b852
blue(text; title="Note") = Markdown.MD(Markdown.Admonition("note", title, [text]))

# ╔═╡ 5590bdc6-a943-4614-a24d-335512b0f73f
if η == 0.
	tan(md"""When η = 0, the estimated trajectory is exactly the master equation solution. There is no backaction.""",title="η = 0 (no signal collected)")
	
elseif η == 1
	blue(md"""When η = 1, the whole leaked signal is collected, and the state remains pure.""",title="η = 1 (all signal collected)")
	
elseif η < η0
	red(md"""As the signal collection efficiency η decreases, the backaction strength decreases. The trajectory becomes less noisy.""", title="η decreased from 0.3")

elseif η > η0
	green(md"""As the signal collection efficiency η increases, the backaction strength increases. The trajectory becomes more noisy.""",title="η increased from 0.3")
		
end

# ╔═╡ c8aa6f9a-a2b7-45d7-9e39-344380a3909b
hint(text; title="Hint") = Markdown.MD(Markdown.Admonition("hint", title, [text]))

# ╔═╡ dfcd22dc-5517-4ecb-99d1-2151d8c73005
hint(md"To see the most dramatic effects, try setting one parameter to an extreme value and then change the other parameter.")

# ╔═╡ 6bb8877f-d861-492e-b816-2396fcd6d59c
md"""
!!! warning "Edit"

	We need to use simple master equation API rather than `bayesian` here. I'm not sure this currently exists in `QuantumCircuits.jl`, but should probably be the `jump-no-jump` method. Or does `bayesian` automatically do that in the absence of measurement?
"""

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
# ╠═4c45fe3e-cd69-11eb-20cd-7bfb98c040cf
# ╟─377a3336-20bd-4baa-a033-af8bbc8668a8
# ╟─7f1176b8-f6d3-4fd1-a23e-31710dcfff10
# ╟─982565bc-b137-4c10-8e13-5fe37be86823
# ╠═4284173a-be05-4b58-a8d9-7189301344fd
# ╟─85206c80-af15-4e1b-81e9-04bcb1083063
# ╟─629669ac-006a-4774-ae3e-67787bd0fa0a
# ╟─259ed22f-b0c0-4fbe-9c6b-97a8730b6b14
# ╟─bbbc92cd-8daa-4ed8-8bf0-d10461af4edd
# ╠═b290f188-47f9-4af1-a5f0-11d99c29c1e3
# ╟─5b48435b-da1b-4547-ac7c-72c95e1d07ff
# ╠═99f2b351-877a-43c8-a508-d5688f92ae0c
# ╠═677b683b-ba7c-4b01-be5b-40ef97a075cd
# ╟─1aa85cc1-32d2-4fa9-a29d-b042d4d4b94a
# ╠═139565f9-0219-4939-a470-4252397dd6a0
# ╠═7218d350-0b14-4fed-b87e-65a9eb573780
# ╟─62170241-f050-4634-ad6a-8842ac96096c
# ╠═d1544d60-7de6-4f54-bffb-9ea0c3ddddc1
# ╟─4d42ca32-5348-48b5-8a24-eafa5b303690
# ╠═2b5ae77e-07c0-4750-85dd-138dcf84da1f
# ╟─ee6b82a5-52b0-4467-b73a-d6534de57202
# ╠═8696fe89-1c19-4852-afe0-2b642f5a9382
# ╟─0b2ba4be-0213-4b3b-9943-291381c2b4fc
# ╟─0cb6b898-347c-4baf-988f-d1dfac1d6557
# ╟─7961752d-faee-4799-99b1-3bd3d17d1f31
# ╟─dff8aec2-7216-49e9-a19d-81b7ff79eb4d
# ╟─6391d48d-4332-4209-bfe1-bec2e60e30c1
# ╟─8638af64-25f9-4793-bc8b-5d22c3ef70fa
# ╟─2db9b4c4-0256-4e29-990c-2ca3ee9c7ddf
# ╟─5590bdc6-a943-4614-a24d-335512b0f73f
# ╟─4585fc96-cc13-476d-99d6-de839c4ec691
# ╟─02ce38f3-f5e1-4fea-9e6e-b887baa548e7
# ╠═09b609f6-a441-4ae2-aaff-66c1b8b6ffae
# ╠═cdbd897e-d689-42fa-9977-2b3169b29d08
# ╟─dfcd22dc-5517-4ecb-99d1-2151d8c73005
# ╟─5e8b047e-b70f-4478-a8b2-a2e48b4a56af
# ╟─7e62e6c6-aa3c-4350-901c-15d017b8db42
# ╠═b2235afc-1bc9-4cf4-ae21-6f535eb37e96
# ╠═31371407-5798-49a4-a83a-2c068d622c4b
# ╠═99e90182-1ee7-45ef-8656-cc5a1e01d564
# ╠═53f32791-a721-4d80-96ce-fae9f533dbaa
# ╟─7d4f4a02-f3b8-4369-aad3-03ff4bc7ec0b
# ╠═0a2219d8-2fbd-4b03-867f-966c4a595e78
# ╠═a255f7a5-d617-4444-80a4-4e4e7762b236
# ╠═6e86b5ee-5b38-4245-b555-c208151c0c53
# ╠═fb32f241-34bb-413b-99d2-c0155b363670
# ╠═e730e11b-6fa5-4121-80eb-69a59d74b852
# ╠═c8aa6f9a-a2b7-45d7-9e39-344380a3909b
# ╟─6bb8877f-d861-492e-b816-2396fcd6d59c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
