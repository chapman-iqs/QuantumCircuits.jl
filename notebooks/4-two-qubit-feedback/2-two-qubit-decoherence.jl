### A Pluto.jl notebook ###
# v0.19.5

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

# ╔═╡ df539f73-7c8a-4d65-b909-0d94720f724f
begin
	directory_name = "QuantumCircuits.jl"
	path = let 
			arr = split(pwd(), "/")
			index = findfirst(s -> s == directory_name, arr)
			join(map(a -> string(a, "/"), arr[1:index]))
	end
	cd(path)
	
	import Pkg
	Pkg.activate(".")
	
	using PlutoUI
	using Random
	using LaTeXStrings
	using Statistics
	using Distributions
	using QuantumCircuits
	using Plots
	using Plots.Measures
	using StatsPlots
	using ProgressMeter
	using LsqFit
	
	include("utilities/two-qubit-operators.jl")
	include("utilities/plotting.jl")
	include("utilities/utilities.jl")

	include("notebooks/table-of-contents.jl")
	include("resources.jl")
	
	md" # Packages and julia files"
end

# ╔═╡ e625f413-e93b-48c9-b14d-181c37b0fe0a
md"""
In this notebook, we check what happens to entangled qubit states under different decoherence mechanisms. We also check how the feedback works after a decoherence or decay event has taken place (including knowledge and ignorance of the event).
"""

# ╔═╡ 249f33ea-92aa-4b9c-bdfa-d35fe5254225
mdp(table_of_contents📔)

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents(title="Two-qubit decoherence mechanisms")

# ╔═╡ 7cb4592d-4b47-41c6-9cb3-cca1033b6c9c
md"""
# Decoherence mechanisms

We introduce some common decoherence mechanisms and their action on two-qubit states.
"""

# ╔═╡ 81394981-ea8f-4718-b4b3-138546ca8a71
md"""
#### System parameters
"""

# ╔═╡ 41954632-a8c7-49d1-9c7e-559c0b0c0556
begin
	dt = 1e-3  # integration time-step
	tf = 10.0
end

# ╔═╡ 38cb1efe-da72-4e1e-a03e-36f7926f986f
md"""
## Random phase rotations
"""

# ╔═╡ 74e7c611-b47d-4c76-a2e0-d748ff90d88b
md"""
We assume random $\hat \sigma_z$ rotations of the form

$\hat H_z(t) = R_1(t) \hat \sigma_1^z + R_2(t) \hat \sigma_2^z$

where $R_1, R_2 \sim \mathcal{N}(0, \sigma)$ are i.i.d. Gaussian distributed random variables with mean $0$ and standard deviation $\sigma$.


"""

# ╔═╡ 17ba9107-98f4-4075-8b97-95f7803a2b8d
md"""
Choose the noise level for the system:
$(@bind noise_str Select(["σ = 1.0 (Τ2 > 100 μs)", "σ = 2.5 (T2 ~ 40 μs)", "σ = 5.0 (T2 ~ 10 μs)", "σ = 10.0 (T2 ~ 2 μs)"], default="σ = 2.0 (T2 ~ 40 μs)"))
"""

# ╔═╡ eb4f3f7a-e2d7-4a99-864d-d52e63f94cb2
md"""
Choose the initial subspace:
$(@bind subspace_str Select(["product", "entangled"]))
"""

# ╔═╡ c471d19f-1c59-4244-a38f-b8b6e9bea368
if subspace_str == "entangled"
	md"""
	Choose the initial system state: $\psi_0$ = 
	$(@bind init_str Select(["|Ψ+>", "|Ψ->", "|Φ+>", "|Φ->"]))
	"""

 else
	md"""
	Choose the initial system state: $\psi_0$ = 
	$(@bind init_str1 Select(["|0>", "|1>", "|+>", "|->"])) ⊗ $(@bind init_str2 Select(["|0>", "|1>", "|+>", "|->"]))
	"""

 end


# ╔═╡ 9bbc57d1-ea42-4ff7-ae7c-28093ff8ddd2
begin
	noise_dict = Dict(["σ = 1.0 (Τ2 > 100 μs)" => 1.0, 
					"σ = 2.5 (T2 ~ 40 μs)" => 2.5,
					"σ = 5.0 (T2 ~ 10 μs)" => 5.0, 
					"σ = 10.0 (T2 ~ 2 μs)" => 10.0])
	σ = noise_dict[noise_str]    

	init_dict = Dict([	"|Ψ+>" => Ψp, 
						"|Ψ->" => Ψm, 
						"|Φ+>" => Φp, 
						"|Φ->" => Φm,
						"|0>" => g, 
						"|1>" => e, 
						"|+>" => normalize(g + e), 
						"|->" => normalize(g - e)])
	if subspace_str == "entangled"
		ψ0 = init_dict[init_str]
	else
		ψ0 = init_dict[init_str1] ⊗ init_dict[init_str2] 
	end
end

# ╔═╡ 63d6d9e6-c7c6-408b-b6be-f719f0af9c4e
Hz(t) = let
			R1, R2 = σ * randn(), σ * randn()
			R1 * σz1 + R2 * σz2
		end

# ╔═╡ 4ebd2e2b-c0a1-4fb4-819b-2baab3e3afce
sol1 = bayesian((0, tf), ψ0, Hz, [], []; dt=dt)

# ╔═╡ 66b87a94-5f06-4bb7-b244-3c166c7a133e
md" ### Correspondence to T2 of a single qubit"

# ╔═╡ 04e7b659-0501-4ce2-818e-5bc1d15f7f7d
md"""
We can extract the effective $T_2$ corresponding to the chosen $\sigma$ by fitting the exponential decay for a single qubit:
"""

# ╔═╡ cdea1cd7-9271-4aa0-b972-cd52c8a1caae
H1(σ) = t -> σ * randn() * σz

# ╔═╡ 13b8e451-4af1-4d60-8702-754096555c56
begin

	N = 500 # number of noise realizations
	ψ0q = normalize(e + g) # initialize in in the x-y plane

	sols = map(1:N) do m
		bayesian((0, tf), ψ0q, H1(σ), [], []; dt=dt)
	end

end

# ╔═╡ 4a81e2d4-6d19-41f4-8267-c81d59e309aa
function get_T2(ensemble)

	t = sols[1].t
	xs, ys, zs = [], [], []

	for sol in ensemble
		x, y, z = [map(ρ -> real(expect(op, ρ)), sol.ρ) for op in [σx, σy, σz]]
		for (list, traj) in zip([xs, ys, zs], [x, y, z])
			push!(list, traj)
		end
		
	end

	N = length(xs)
	xavg = [mean([xs[i][j] for i in 1:N]) for j in 1:length(t)]

	# fit model
	model(t, p) = p[1] * exp.(-p[2] * t)
	p0 = [0.5, 0.5]
	fit = curve_fit(model, t, xavg, p0)
	param = fit.param
	T2 = 1/(2 * param[2])

	return T2

end

# ╔═╡ f9a003ed-9f2a-48b0-b035-b61de8c79849
md"""
## Qubit swaps
"""

# ╔═╡ 389fd79b-f5df-4f57-8241-5e4382ea2ed9
md"""
## Spontaneous emission
"""

# ╔═╡ bc68de97-beef-47a5-bd37-5eaca4a740c6
md"""
# Feedback under decoherence
"""

# ╔═╡ 32e30f5f-ec56-4586-aa7d-cfc6d410e899
md"""
## Assuming detection of decohering event
"""

# ╔═╡ 79b0d8f8-0601-41fb-8927-3437ccc1e011
md"""
### Random phase rotations
"""

# ╔═╡ 614eb04d-3eb0-4541-8e22-aaa1024a19e7
md"""
### Qubit swaps
"""

# ╔═╡ 99be4964-9236-4b08-b61d-e797a8fb8bbc
md"""
### Spontaneous emission
"""

# ╔═╡ d4818d0e-fb6d-4005-8c04-df7a73f3c6d5
md"""
## Assuming ignorance of decohering event
"""

# ╔═╡ 9a1b96b4-7f3b-4ddf-8e39-a7b51b7030a3
md"""
### Random phase rotations
"""

# ╔═╡ e9a15c74-2f80-47f4-9b9c-304ec972eb81
md"""
### Qubit swaps
"""

# ╔═╡ 479a6ccb-1d71-49f4-9ef1-03efb30513e7
md"""
### Spontaneous emission
"""

# ╔═╡ 200afa00-0f0f-465e-bb26-009f9da78e68
md"""
### Single noise realization, no measurement
"""

# ╔═╡ 133b6939-10b2-4c8e-acf8-5658ca96a0f9
md" # Utilities"

# ╔═╡ ea04744b-4296-4dc2-9a3c-1f477c96f1ac
md"""
## Plotting
"""

# ╔═╡ bf00661d-719d-45f8-b649-fe93b1174a10
begin
	colors1q = palette(:tab10)
	colors_system = palette(:rainbow)
	colors_filter = palette(:lightrainbow)	
end

# ╔═╡ 944e3623-39b5-448c-a4fd-d879d3838d4f
begin
	p = plot()
	t = sols[1].t
	xs, ys, zs = [], [], []
	
	for sol in sols
		x, y, z = [map(ρ -> real(expect(op, ρ)), sol.ρ) for op in [σx, σy, σz]]
		for (list, traj) in zip([xs, ys, zs], [x, y, z])
			push!(list, traj)
		end
		
	end

	for x in xs[1:50]
		plot!(t, x, alpha=0.1, label=:none, color=colors1q[1])
	end
	plot!(t, xs[51], label="x trajectory", color=colors1q[1])

	xavg = [mean([xs[i][j] for i in 1:N]) for j in 1:length(t)]

	# fit model
	model(t, p) = p[1] * exp.(-p[2] * t)
	p0 = [0.5, 0.5]
	fit = curve_fit(model, t, xavg, p0)
	param = fit.param
	T2 = 1/(2 * param[2])
	

	# plot
	plot!(t, xavg, alpha=1, color=:black, label="x average", linewidth=3)
	plot!(t, map(tt -> model(tt, param), t), alpha=1, color=:red, linestyle=:dash, linewidth=3, label="exponential fit", legend=:right)

end

# ╔═╡ de114db9-c561-4b5e-bbcf-08b8108e69e3
md" The model fit found T2 = $(round(T2, digits=2)) for σ = $σ."

# ╔═╡ bb21a169-8a65-4660-84d2-365563b6dd4b
function single_qubit_plots(sol::Solution)
	colors = colors1q 
	
	# calculate expectation values --------------------------------------------
	q1_basis = [σx1, σy1, σz1]
	q2_basis = [σx2, σy2, σz2]
	
	exps1 = map(op -> expectations(sol, op), q1_basis)
	exps2 = map(op -> expectations(sol, op), q2_basis)

	qlabels = ["x", "y", "z"]

	p1s = 0.5 * (1 .+ exps1[1].^2 .+ exps1[2].^2 .+ exps1[3].^2)
	p2s = 0.5 * (1 .+ exps2[1].^2 .+ exps2[2].^2 .+ exps2[3].^2)
	
	# plot ----------------------------------------------------------------------
	l = @layout [bloch1{0.5h}; bloch2{0.5h}]

	p1 = plot(sol.t, p1s, color=colors[4], label=L"Tr(\rho_1^2)", linestyle=:dash, title="single-qubit states")

	for l in 1:3
		label = qlabels[l]
		color = colors[l]
		exp = exps1[l]
		plot!(sol.t, exp, color=color, label=label, legend=:outerright, ylims=[-1,1])
	end

	p2 = plot(sol.t, p2s, color=colors[4], label=L"Tr(\rho_2^2)", linestyle=:dash, title="")

	for l in 1:3
		label = qlabels[l]
		color = colors[l]
		exp = exps2[l]
		plot!(sol.t, exp, color=color, label=label, legend=:outerright, xlabel="t (μs)", ylims=[-1,1])
	end
	
	plot(p1, p2, layout = l, link=:y, size=(600,300), legendfontsize=8, titlefontsize=12, legend=:outerright)
	
end

# ╔═╡ d2c8faa0-4d6e-49bd-98cb-4447640036d9
single_qubit_plots(sol1)

# ╔═╡ 70975e08-6513-4543-b20c-9a594416a4f0
function bell_plot(sol::Solution)

	basis = bell_basis
	colors = colors_system
	labels = bell_basis_labels
	title = "Bell states"

	exps = map(op -> expectations(sol, dm(op)), basis)
	
	pl = plot(size=(600,300), legendfontsize=12, titlefontsize=12, legend=:outerright, title=title)
	
	for l in 1:length(basis)
		label = labels[l]
		color = colors[l]
		exp = exps[l]
		plot!(sol.t, exps[l], color=color, label=label, legend=:outerright, xlabel="t (μs)", ylims=[0,1])
	end

	pl

end

# ╔═╡ 3488df41-3805-445c-8196-ff4a60dd2efa
function bell_plot(sys::Solution, fil::Solution)

	basis = bell_basis
	labels = bell_basis_labels
	title = "Bell states"

	exps_sys = map(op -> expectations(sys, dm(op)), basis)
	exps_fil = map(op -> expectations(fil, dm(op)), basis)
	
	pl = plot(size=(600,300), legendfontsize=12, titlefontsize=12, legend=:outerright, title=title, xlabel="t (μs)",)
	
	for l in 1:length(basis)
		plot!(sys.t, exps_sys[l], color=colors_system[l], label=labels[l], legend=:outerright, ylims=[0,1])
		plot!(fil.t, exps_fil[l], color=colors_filter[l], label=:none, linestyle=:dash)
	end

	pl

end

# ╔═╡ 2b745af1-3787-4888-b26a-3b278c6693ed
bell_plot(sol1)

# ╔═╡ 575fb11f-3eac-4f9a-bfd0-c74ccc726ab2
md"""
## Old utilities (delete)
"""

# ╔═╡ c617a096-895d-499a-8d1a-a35c21dfbb7d
function plotresults3(sol::QuantumCircuits.solution)
	# calculate expectation values --------------------------------------------
	ns, n1s, n2s, x1s, x2s, y1s, y2s, z1s, z2s = [map(ρ -> real(expect(ρ, op)), sol.ρ) for op in [n, n1, n2, σx1, σx2, σy1, σy2, σz1, σz2]]

	exps = [map(ρ -> real(expect(ρ, op)), sol.ρ) for op in [ket0, ket1p, ket1m, ket2]]

	p1s = 0.5 * (1 .+ x1s.^2 .+ y1s.^2 .+ z1s.^2)
	p2s = 0.5 * (1 .+ x2s.^2 .+ y2s.^2 .+ z2s.^2)

	ps = map(ρ -> real(expect(ρ, ρ)), sol.ρ)

	
	# plot ----------------------------------------------------------------------
	l = @layout [nexp{0.25h}; basisstates{0.25h}; bloch1{0.25h}; bloch2{0.25h}]


	p2 = plot(sol.t, ns, ylabel=L"\langle n \rangle", color=colors2[8], title="excitation number", legend=:false, ylims=[0,2.2])

	p3 = plot()
	for (i, eval) in enumerate(exps)
		plot!(sol.t, eval, color=colors2[i], label=fs_labels[i], legend=:outerright)
	end
	plot!(sol.t, ps, color=colors[4], label=L"Tr(\rho^2)", title="2-qubit states", linestyle=:dash)

	p4 = plot(sol.t, x1s, color=colors[1], label=L"x", legend=:outerright)
	plot!(sol.t, y1s, color=colors[2], label=L"y")
	plot!(sol.t, z1s, color=colors[3], label=L"z")
	plot!(sol.t, p1s, color=colors[4], label=L"Tr(\rho_1^2)", title="qubit 1")

	p5 = plot(sol.t, x2s, color=colors[1], label=L"x", legend=:outerright)
	plot!(sol.t, y2s, color=colors[2], label=L"y")
	plot!(sol.t, z2s, color=colors[3], label=L"z")
	plot!(sol.t, p2s, color=colors[4], label=L"Tr(\rho_2^2)", title="qubit 2")

	plot(p2, p3, p4, p5,  layout = l, link=:y, size=(600,800), legendfontsize=8, titlefontsize=12, legend=:outerright)
	
end

# ╔═╡ ece2cf9e-45c5-4f05-b61a-cb1eb1acfc72
function plotresults4(sol::QuantumCircuits.solution)
	# calculate expectation values --------------------------------------------
	ns, n1s, n2s, x1s, x2s, y1s, y2s, z1s, z2s = [map(ρ -> real(expect(ρ, op)), sol.ρ) for op in [n, n1, n2, σx1, σx2, σy1, σy2, σz1, σz2]]

	exps = [map(ρ -> real(expect(ρ, dm(op))), sol.ρ) for op in [ketp, ket1p, ket1m, ketm]]

	p1s = 0.5 * (1 .+ x1s.^2 .+ y1s.^2 .+ z1s.^2)
	p2s = 0.5 * (1 .+ x2s.^2 .+ y2s.^2 .+ z2s.^2)

	ps = map(ρ -> real(expect(ρ, ρ)), sol.ρ)

	
	# plot ----------------------------------------------------------------------
	l = @layout [nexp{0.25h}; basisstates{0.25h}; bloch1{0.25h}; bloch2{0.25h}]


	p2 = plot(sol.t, ns, ylabel=L"\langle n \rangle", color=colors2[8], title="excitation number", legend=:false, ylims=[0,2.2])

	p3 = plot()
	for (i, eval) in enumerate(exps)
		plot!(sol.t, eval, color=colors2[i], label=rabi_labels[i], legend=:outerright, xlabel="t (μs)")
	end
	plot!(sol.t, ps, color=colors[4], label=L"Tr(\rho^2)", title="2-qubit states", linestyle=:dash)

	p4 = plot(sol.t, x1s, color=colors[1], label=L"x", legend=:outerright)
	plot!(sol.t, y1s, color=colors[2], label=L"y")
	plot!(sol.t, z1s, color=colors[3], label=L"z")
	plot!(sol.t, p1s, color=colors[4], label=L"Tr(\rho_1^2)", title="qubit 1", xlabel="t (μs)")

	p5 = plot(sol.t, x2s, color=colors[1], label=L"x", legend=:outerright)
	plot!(sol.t, y2s, color=colors[2], label=L"y")
	plot!(sol.t, z2s, color=colors[3], label=L"z")
	plot!(sol.t, p2s, color=colors[4], label=L"Tr(\rho_2^2)", title="qubit 2")

	plot(p2, p3, p4, p5,  layout = l, link=:y, size=(600,800), legendfontsize=8, titlefontsize=12, legend=:outerright)
	
end

# ╔═╡ 4efec7f2-8a4e-4c4d-9fe9-53754fe35065
function plotresults_3q_entangled(sol::QuantumCircuits.solution)
	# calculate expectation values --------------------------------------------
	ns, x1s, x2s, x3s, y1s, y2s, y3s, z1s, z2s, z3s = [map(ρ -> real(expect(ρ, op)), sol.ρ) for op in [n_3, σx1_3, σx2_3, σx3_3, σy1_3, σy2_3, σy3_3, σz1_3, σz2_3, σz3_3]]

	exps = [map(ρ -> real(expect(ρ, dm(op))), sol.ρ) for op in [ket0_3, ket1_3(0,0), ket2_3(0,0), ket3_3]]
	labs_3 = [L"| \overline{0}\rangle", L"|\overline{1}_{0,0}\rangle", L"|\overline{2}_{0,0}\rangle", L"|\overline{3}\rangle"]

	p1s = 0.5 * (1 .+ x1s.^2 .+ y1s.^2 .+ z1s.^2)
	p2s = 0.5 * (1 .+ x2s.^2 .+ y2s.^2 .+ z2s.^2)
	p3s = 0.5 * (1 .+ x3s.^2 .+ y3s.^2 .+ z3s.^2)

	ps = map(ρ -> real(expect(ρ, ρ)), sol.ρ)

	
	# plot ----------------------------------------------------------------------
	l = @layout [nexp{0.25h}; basisstates{0.25h}]


	p1 = plot(sol.t, ns, color=colors2[8], label=L"\langle \hat n \rangle", ylims=[0,3.2], legend=:outerright)

	p2 = plot()
	for (i, eval) in enumerate(exps)
		plot!(sol.t, eval, color=colors2[i], label=labs_3[i], legend=:outerright,  legendfontsize=12)
	end
	plot!(sol.t, ps, color=colors[4], label=L"Tr(\rho^2)", linestyle=:dash, xlabel="t (μs)")

	plot(p1, p2,  layout = l, link=:y, size=(600,800), legendfontsize=10, titlefontsize=12)
	
end

# ╔═╡ 3df47627-fc2f-4712-804d-59eeb42c7db3
function plotresults_3q(sol::QuantumCircuits.solution)
	# calculate expectation values --------------------------------------------
	ns, x1s, x2s, x3s, y1s, y2s, y3s, z1s, z2s, z3s = [map(ρ -> real(expect(ρ, op)), sol.ρ) for op in [n_3, σx1_3, σx2_3, σx3_3, σy1_3, σy2_3, σy3_3, σz1_3, σz2_3, σz3_3]]

	exps = [map(ρ -> real(expect(ρ, dm(op))), sol.ρ) for op in basis_states_3]

	p1s = 0.5 * (1 .+ x1s.^2 .+ y1s.^2 .+ z1s.^2)
	p2s = 0.5 * (1 .+ x2s.^2 .+ y2s.^2 .+ z2s.^2)
	p3s = 0.5 * (1 .+ x3s.^2 .+ y3s.^2 .+ z3s.^2)

	ps = map(ρ -> real(expect(ρ, ρ)), sol.ρ)

	
	# plot ----------------------------------------------------------------------
	l = @layout [nexp{0.2h}; basisstates{0.5h}; bloch1{0.1h}; bloch2{0.1h}; bloch3{0.1h}]


	p1 = plot(sol.t, ns, ylabel=L"\langle n \rangle", color=colors2[8], title="excitation number", label="n", ylims=[0,3.2], legend=:outerright)

	# p2 = plot()
	# for (i, eval) in enumerate(exps)
	# 	plot!(sol.t, eval, color=colors2[i], label=fs_labels_3[i], legend=:outerright)
	# end
	# plot!(sol.t, ps, color=colors[4], label=L"Tr(\rho^2)", title="2-qubit states", linestyle=:dash)

	ls = [:solid, :solid, :dash, :dot, :solid, :dash, :dot, :solid]
		p2 = plot()
	for (i, eval) in enumerate(exps)
		plot!(sol.t, eval, color=colors2[i], label=basis_labels_3[i], legend=:outerright,  legendfontsize=12, linestyle=ls[i])
	end
	plot!(sol.t, ps, color=colors[4], label=L"Tr(\rho^2)", title="3-qubit states", linestyle=:dash)

	p3 = plot(sol.t, x1s, color=colors[1], label=L"x", legend=:outerright)
	plot!(sol.t, y1s, color=colors[2], label=L"y")
	plot!(sol.t, z1s, color=colors[3], label=L"z")
	plot!(sol.t, p1s, color=colors[4], label=L"Tr(\rho_1^2)", title="qubit 1")

	p4 = plot(sol.t, x2s, color=colors[1], label=L"x", legend=:outerright)
	plot!(sol.t, y2s, color=colors[2], label=L"y")
	plot!(sol.t, z2s, color=colors[3], label=L"z")
	plot!(sol.t, p2s, color=colors[4], label=L"Tr(\rho_2^2)", title="qubit 2")

	p5 = plot(sol.t, x3s, color=colors[1], label=L"x", legend=:outerright)
	plot!(sol.t, y3s, color=colors[2], label=L"y")
	plot!(sol.t, z3s, color=colors[3], label=L"z")
	plot!(sol.t, p3s, color=colors[4], label=L"Tr(\rho_3^2)", title="qubit 3")

	plot(p1, p2, p3, p4, p5,  layout = l, link=:y, size=(600,800), legendfontsize=10, titlefontsize=12)
	
end

# ╔═╡ 42025655-1b85-4087-9ecc-11c5a412c6ca
ketp ⊗ ket1p'

# ╔═╡ 5d0a5796-c91b-4cb7-ba52-f55ddd34d0a5
function plotresults5(sol::QuantumCircuits.solution)
	# calculate expectation values --------------------------------------------
	ns, n1s, n2s, x1s, x2s, y1s, y2s, z1s, z2s = [map(ρ -> real(expect(ρ, op)), sol.ρ) for op in [n, n1, n2, σx1, σx2, σy1, σy2, σz1, σz2]]

	exps = [map(ρ -> real(expect(ρ, op)), sol.ρ) for op in [ketp, ket1p, ket1m, ketm]]

	cors = [map(ρ -> real(expect(ρ, op)), sol.ρ) for op in [ketp ⊗ ket1p', ketm ⊗ ket1p']]

	p1s = 0.5 * (1 .+ x1s.^2 .+ y1s.^2 .+ z1s.^2)
	p2s = 0.5 * (1 .+ x2s.^2 .+ y2s.^2 .+ z2s.^2)

	ps = map(ρ -> real(expect(ρ, ρ)), sol.ρ)

	
	# plot ----------------------------------------------------------------------
	l = @layout [nexp{0.2h}; basisstates{0.2h}; cors{0.2h}; bloch1{0.2h}; bloch2{0.2h}]


	p1 = plot(sol.t, ns, label=L"\langle \hat n \rangle", color=colors2[8], ylims=[0,2.2])

	p2 = plot()
	for (i, eval) in enumerate(exps)
		plot!(sol.t, eval, color=colors2[i], label=rabi_labels[i])
	end
	plot!(sol.t, ps, color=colors[4], label=L"Tr(\rho^2)", linestyle=:dash, xlabel="t (μs)")

	p3 = plot()
	plot!(sol.t, - θ.(sol.ρ) / (2td), color=:blue, label=L"\Omega = \theta(t) / 2 t_d")

	p4 = plot(sol.t, x1s, color=colors[1], label=L"x", legend=:outerright)
	plot!(sol.t, y1s, color=colors[2], label=L"y")
	plot!(sol.t, z1s, color=colors[3], label=L"z")
	plot!(sol.t, p1s, color=colors[4], label=L"Tr(\rho_1^2)", title="qubit 1")

	p5 = plot(sol.t, x2s, color=colors[1], label=L"x", legend=:outerright)
	plot!(sol.t, y2s, color=colors[2], label=L"y")
	plot!(sol.t, z2s, color=colors[3], label=L"z")
	plot!(sol.t, p2s, color=colors[4], label=L"Tr(\rho_2^2)", title="qubit 2")

	plot(p1, p2, p3, p4, p5, layout = l, link=:y, size=(600,800), legendfontsize=8, titlefontsize=12, legend=:outerright)
	
end

# ╔═╡ b4fd50f9-7073-4233-838a-7702dfa74c27


# ╔═╡ fc3a9a28-42c5-4c9c-bdc9-8745f02cbe83


# ╔═╡ f3ac9b41-b123-4236-8c62-8e2b988463c2
begin
	@userplot MyPlot
	
	@recipe function f(mp::MyPlot; add_marker=false)
		
		x, y = mp.args
		
		linecolor   --> :blue
		seriestype  :=  :path
		markershape --> (add_marker ? :circle : :none)
		legend := :none
		
		@series begin
			x, y
		end
	end
	
end

# ╔═╡ 8b4bc294-989b-4f1c-82b9-e418fdd3c40b
md"""
### Other
"""

# ╔═╡ 87f29d4c-63d5-4d2f-9126-90bbd67fe2bd
f(b::Vector) = "this is a vector"

# ╔═╡ 2b310f40-0132-41b5-b3c2-9c334bf5b63e
f(k::Operator) = "this is an operator"

# ╔═╡ 235cee23-c8af-4de2-a1c3-a2b173156703
expects(ops) = ρ -> collect(real(expect(ρ, s)) for s in vcat(ops, ρ)) # ρ -> [<x>,<y>,<z>,<ρ>]

# ╔═╡ e27bd39c-58b7-4c5c-a677-3fe70f500ee8
xyz(θ, ϕ) = (sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ))

# ╔═╡ 28d8df47-07cf-4a0d-a447-f894d021b2bc
begin
	mutable struct traj
		t::Vector{Float64}
		x::Vector{Float64}
		y::Vector{Float64}
		z::Vector{Float64}
		p::Vector{Float64}
		r
	end
	
	function traj(t, ρ, r)
		x, y, z = [real(expect(σi, ρ)) for σi in (σx, σy, σz)]
		p = real(expect.(ρ, ρ))
		traj(t, x, y, z, p, r)
	end
	
	function traj(sol::QuantumCircuits.solution; resonator=false)
		t, ρ, r = (sol.t, sol.ρ, sol.r)
		x, y, z = resonator ? 
					[real(expect(σi ⊗ id, ρ)) for σi in (σx, σy, σz)] :
					[real(expect(σi, ρ)) for σi in (σx, σy, σz)]
		
		p = 0.5 .* (1 .+ x.^2 .+ y.^2 .+ z.^2)
		traj(t, x, y, z, p, r)
	end
	
end

# ╔═╡ dc3e2dcb-5bf1-492d-8337-f366dcf0170b
φdict = Dict("0" => 0, 
				"π/8" => π/8, 
				"π/4" => π/4,
				"3π/8" => 3π/8,
				"π/2" => π/2)

# ╔═╡ 460379bc-0d49-434a-b0bc-3efcfdc47b5c


# ╔═╡ 01523c93-5737-4d94-87fc-2e5fb731002c
green(text; title="Note") = Markdown.MD(Markdown.Admonition("correct", title, [text]))

# ╔═╡ de310c78-ae02-488f-a939-e4d29faa3651
red(text; title="Note") = Markdown.MD(Markdown.Admonition("danger", title, [text]))

# ╔═╡ f8cfd829-06ef-4971-b216-5a6abe1b072d
tan(text; title="Note") = Markdown.MD(Markdown.Admonition("warning", title, [text]))

# ╔═╡ 87e2a8c9-75b6-486e-95b3-2506dd255992
blue(text; title="Note") = Markdown.MD(Markdown.Admonition("note", title, [text]))

# ╔═╡ Cell order:
# ╟─e625f413-e93b-48c9-b14d-181c37b0fe0a
# ╟─249f33ea-92aa-4b9c-bdfa-d35fe5254225
# ╟─3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─7cb4592d-4b47-41c6-9cb3-cca1033b6c9c
# ╟─81394981-ea8f-4718-b4b3-138546ca8a71
# ╠═41954632-a8c7-49d1-9c7e-559c0b0c0556
# ╟─38cb1efe-da72-4e1e-a03e-36f7926f986f
# ╟─74e7c611-b47d-4c76-a2e0-d748ff90d88b
# ╠═63d6d9e6-c7c6-408b-b6be-f719f0af9c4e
# ╠═17ba9107-98f4-4075-8b97-95f7803a2b8d
# ╟─eb4f3f7a-e2d7-4a99-864d-d52e63f94cb2
# ╟─c471d19f-1c59-4244-a38f-b8b6e9bea368
# ╟─9bbc57d1-ea42-4ff7-ae7c-28093ff8ddd2
# ╟─4ebd2e2b-c0a1-4fb4-819b-2baab3e3afce
# ╟─d2c8faa0-4d6e-49bd-98cb-4447640036d9
# ╟─2b745af1-3787-4888-b26a-3b278c6693ed
# ╟─66b87a94-5f06-4bb7-b244-3c166c7a133e
# ╟─04e7b659-0501-4ce2-818e-5bc1d15f7f7d
# ╠═cdea1cd7-9271-4aa0-b972-cd52c8a1caae
# ╠═13b8e451-4af1-4d60-8702-754096555c56
# ╠═944e3623-39b5-448c-a4fd-d879d3838d4f
# ╟─de114db9-c561-4b5e-bbcf-08b8108e69e3
# ╟─4a81e2d4-6d19-41f4-8267-c81d59e309aa
# ╟─f9a003ed-9f2a-48b0-b035-b61de8c79849
# ╟─389fd79b-f5df-4f57-8241-5e4382ea2ed9
# ╟─bc68de97-beef-47a5-bd37-5eaca4a740c6
# ╟─32e30f5f-ec56-4586-aa7d-cfc6d410e899
# ╟─79b0d8f8-0601-41fb-8927-3437ccc1e011
# ╟─614eb04d-3eb0-4541-8e22-aaa1024a19e7
# ╟─99be4964-9236-4b08-b61d-e797a8fb8bbc
# ╟─d4818d0e-fb6d-4005-8c04-df7a73f3c6d5
# ╟─9a1b96b4-7f3b-4ddf-8e39-a7b51b7030a3
# ╟─e9a15c74-2f80-47f4-9b9c-304ec972eb81
# ╟─479a6ccb-1d71-49f4-9ef1-03efb30513e7
# ╟─200afa00-0f0f-465e-bb26-009f9da78e68
# ╟─133b6939-10b2-4c8e-acf8-5658ca96a0f9
# ╠═ea04744b-4296-4dc2-9a3c-1f477c96f1ac
# ╠═bf00661d-719d-45f8-b649-fe93b1174a10
# ╠═bb21a169-8a65-4660-84d2-365563b6dd4b
# ╠═70975e08-6513-4543-b20c-9a594416a4f0
# ╠═3488df41-3805-445c-8196-ff4a60dd2efa
# ╟─575fb11f-3eac-4f9a-bfd0-c74ccc726ab2
# ╟─c617a096-895d-499a-8d1a-a35c21dfbb7d
# ╠═ece2cf9e-45c5-4f05-b61a-cb1eb1acfc72
# ╟─4efec7f2-8a4e-4c4d-9fe9-53754fe35065
# ╟─3df47627-fc2f-4712-804d-59eeb42c7db3
# ╠═42025655-1b85-4087-9ecc-11c5a412c6ca
# ╠═5d0a5796-c91b-4cb7-ba52-f55ddd34d0a5
# ╠═b4fd50f9-7073-4233-838a-7702dfa74c27
# ╠═fc3a9a28-42c5-4c9c-bdc9-8745f02cbe83
# ╠═f3ac9b41-b123-4236-8c62-8e2b988463c2
# ╟─8b4bc294-989b-4f1c-82b9-e418fdd3c40b
# ╠═87f29d4c-63d5-4d2f-9126-90bbd67fe2bd
# ╠═2b310f40-0132-41b5-b3c2-9c334bf5b63e
# ╟─235cee23-c8af-4de2-a1c3-a2b173156703
# ╟─e27bd39c-58b7-4c5c-a677-3fe70f500ee8
# ╠═28d8df47-07cf-4a0d-a447-f894d021b2bc
# ╟─dc3e2dcb-5bf1-492d-8337-f366dcf0170b
# ╠═460379bc-0d49-434a-b0bc-3efcfdc47b5c
# ╠═01523c93-5737-4d94-87fc-2e5fb731002c
# ╠═de310c78-ae02-488f-a939-e4d29faa3651
# ╠═f8cfd829-06ef-4971-b216-5a6abe1b072d
# ╠═87e2a8c9-75b6-486e-95b3-2506dd255992
# ╠═df539f73-7c8a-4d65-b909-0d94720f724f
