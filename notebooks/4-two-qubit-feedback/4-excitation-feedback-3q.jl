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

# ╔═╡ 4a66d07c-4828-43d2-9c1a-de5475d966c9
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

	include("utilities/three-qubit-operators.jl")
	include("utilities/utilities.jl")
	include("utilities/plotting.jl")

	include("notebooks/table-of-contents.jl")
	include("resources.jl")
	
	md" # Packages and julia files"
end

# ╔═╡ c2f8e743-48bf-4c1b-a199-b03cabba8121
mdp("In this interactive notebook, we show similar examples to ", excitation_feedback📔, ", but for a three-qubit Hilbert space.)")

# ╔═╡ b7976c8f-903b-4700-8cf1-5d9a00bed17f
mdp(table_of_contents📔)

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents(title="Three-qubit excitation feedback")

# ╔═╡ 0b2c3e3c-831a-45df-8adb-33d2af64e210
md"""
# Three qubit example
"""

# ╔═╡ 0dabf623-f0ea-4a77-92ed-4b30ddd3a800
md"""
## System description
"""

# ╔═╡ 6f2d9396-110d-41f8-870e-bd70c948f587
md"""
For the three-qubit example, we can define a standard product eigenbasis:

$\big\{ \ket{000}, \ket{001}, \ket{010}, \ket{100}, \ket{011}, \ket{101}, \ket{110}, \ket{111} \big\}.$

Under measurement of the total excitation (number) operator $\hat n$, we can define number eigenstates

$\ket{\overline{0}} \equiv \ket{000}$

$\ket{\overline{1}}_{\alpha, \beta, \gamma} = \alpha_1 \ket{001} + \beta_1 \ket{010} + \gamma_1 \ket{100}$

$\ket{\overline{2}}_{\alpha, \beta, \gamma} = \alpha_2 \ket{110} + \beta_2 \ket{101} + \gamma_2 \ket{011}$

$\ket{\overline{3}} \equiv \ket{111}$

such that $|\alpha_i|^2 + |\beta_i|^2 + |\gamma_i|^2 = 1$, and $\hat n \ket{\overline i} = i \ket{\overline i}$. The flexibility in choice of $(\alpha_i, \beta_i, \gamma_i)$ reflects the degeneracy of eigenstates of $\hat n$ with eigenvalues $1$ or $2$.

To facilitate plotting, we'll define the POVM $\{\Lambda_i \}_{i=0}^3$ such that

$\Lambda_0 \equiv \ket{000}\bra{000},$

$\Lambda_1 \equiv \ket{001}\bra{001} + \ket{010}\bra{010} + \ket{100}\bra{100},$

$\Lambda_2 \equiv \ket{110}\bra{110} + \ket{101}\bra{101} + \ket{011}\bra{011},$

$\Lambda_3 \equiv \ket{111}\bra{111}.$


"""

# ╔═╡ 830e6b02-0de2-4802-b212-2f8305008dfa
md"""
## Measurement-only evolution
"""

# ╔═╡ ff04e9d3-2890-437f-992f-b55761ed12dd
md"""
## Feedback
"""

# ╔═╡ c67b7c75-5526-4e54-8e96-a0870b3b4d1b
md"""
Choose the target subspace:
$(@bind target_str Select(["Λ0", "Λ1", "Λ2", "Λ3"], default="Λ1"))
"""

# ╔═╡ f54e0863-2dd3-4719-8f69-fee15883a924
begin
	 target_dict = Dict("Λ0" => Λ0, "Λ1" => Λ1, "Λ2" => Λ2, "Λ3" => Λ3)
	 target = target_dict[target_str]
end

# ╔═╡ b6598bb0-c52b-46fe-acce-604c92e288a0
md"""
### On optimal Rabi axis
"""

# ╔═╡ cfc89850-d76d-4228-9326-78ac85a0e80d
md"""
### On expectation of n
"""

# ╔═╡ 501439e9-0343-4ba5-8df4-6a87b456ba5f
md"""
Further consideration is needed to understand what is the optimal protocol -- is there an analogue to the "effective qubit" $\theta$ feedback used in the two-qubit example, `excitation-feedback.jl`?
"""

# ╔═╡ 133b6939-10b2-4c8e-acf8-5658ca96a0f9
md" # Utilities"

# ╔═╡ ea04744b-4296-4dc2-9a3c-1f477c96f1ac
md"""
## Plotting
"""

# ╔═╡ 3d84b923-42fa-422b-a9e2-4ed26880e1e9
begin
	colors1q = palette(:tab10)
	colors3q_number = palette(:lightrainbow)
	colors3q_prod = palette(:okabe_ito)
	colors2q_bell = palette(:rainbow)
end

# ╔═╡ 4efec7f2-8a4e-4c4d-9fe9-53754fe35065
function number_plot_3q(sol::Solution)
	# calculate expectation values --------------------------------------------
	basis = number_POVM
	labels = number_labels
	colors = colors3q_number
	title = "3-qubit number states"

	exps = map(op -> expectations(sol, op), basis)

	pl = plot(size=(600,300), legendfontsize=12, titlefontsize=12, legend=:outerright, title=title)
	
	for l in 1:length(basis)
		label = labels[l]
		color = colors[l]
		exp = exps[l]
		plot!(sol.t, exps[l], color=color, label=label, legend=:outerright, xlabel="t (μs)", ylims=[0,1])
	end

	pl
	
end

# ╔═╡ 4525b42b-5421-4416-a5e7-57369d2e90cb
let
	ψ0 = normalize(ket100 + ket011)
	dt = 1e-3  # integration time-step

	Γ = 1.0 # Ensemble measurement dephasing rate (MHz)	
	η = 1.0 # collection efficiency
	tf = 10.0
	
	# Kraus operators --------------------------------------------------------------
	H = 0.0 * I
	
	J = [(n, ((1 - η) * Γ))]
	C = [(n, Γ, η)]
	
	# Random.seed!(s)
	sol = bayesian((0.0, tf), ψ0, H, J, C; dt=dt)

	number_plot_3q(sol)
end

# ╔═╡ 3d38025e-8595-42cd-92b1-958a6a52ddb4
let
	ρ0 = ket000
	dt = 1e-3  # integration time-step

	Γ = 1.0 # Ensemble measurement dephasing rate (MHz)	
	η = 1.0 # collection efficiency
	ΩR = π/2
	tf = 6.0
	
	# Kraus operators --------------------------------------------------------------
	H(t::Timescale, ρ::State) = 
			let
				h(ϕ) = ΩR * (cos(ϕ) * (σy1 + σy2 + σy3) + sin(ϕ) * (σx1 + σx2 + σx3))
				u(ϕ) = exp( -im * dt * DenseOperator(h(ϕ)))
				ρnext(ϕ) = u(ϕ) * ρ * u(ϕ)'

				ϕs = range(0, 2π, step=π/10)
				exps = [real(expect(target, ρn)) for ρn in ρnext.(ϕs)]

				h(ϕs[findmax(exps)[2]])

			end	
	
	J = [(n, ((1 - η) * Γ))]
	C = [(n, Γ, η)]
	
	global sol1 = bayesian((0.0, tf), ρ0, H, J, C; dt=dt)
	
	number_plot_3q(sol1)

	
end

# ╔═╡ f2273cd2-6808-4d4d-ac33-7825deffb72f
let
	ρ0 = ket000
	dt = 1e-3  # integration time-step

	Γ = 1.0 # Ensemble measurement dephasing rate (MHz)	
	η = 1.0 # collection efficiency
	ΩR = π/2

	tf = 6.0
	ntarget = expect(n, normalize(target))
	
	# Kraus operators --------------------------------------------------------------
	H(t::Timescale, ρ::State) = ΩR * (ntarget - real(expect(n, ρ))) * (σx1 + σx2 + σx3)
	
	J = [(n, ((1 - η) * Γ))]
	C = [(n, Γ, η)]
	
	global sol2 = bayesian((0.0, tf), ρ0, H, J, C; dt=dt)
	
	number_plot_3q(sol2)

	
end

# ╔═╡ 11f2dc87-601b-4ca0-9357-7aa1458d01a6
function product_plot_3q(sol::Solution)
	# calculate expectation values --------------------------------------------
	basis = basis_3q
	labels = basis_3q_labels
	colors = colors3q_prod
	title = "3-qubit product states"

	exps = map(op -> expectations(sol, op), basis)

	pl = plot(size=(600,300), legendfontsize=12, titlefontsize=12, legend=:outerright, title=title)
	
	for l in 1:length(basis)
		label = labels[l]
		color = colors[l]
		exp = exps[l]
		plot!(sol.t, exps[l], color=color, label=label, legend=:outerright, xlabel="t (μs)", ylims=[0,1])
	end

	pl
	
end

# ╔═╡ 268ba231-3d96-4c42-b89f-4406db276c27
product_plot_3q(sol1)

# ╔═╡ 9f1ee82c-11d6-446e-a3b7-09a0e13b8377
product_plot_3q(sol2)

# ╔═╡ Cell order:
# ╟─c2f8e743-48bf-4c1b-a199-b03cabba8121
# ╟─b7976c8f-903b-4700-8cf1-5d9a00bed17f
# ╟─3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─0b2c3e3c-831a-45df-8adb-33d2af64e210
# ╟─0dabf623-f0ea-4a77-92ed-4b30ddd3a800
# ╟─6f2d9396-110d-41f8-870e-bd70c948f587
# ╟─830e6b02-0de2-4802-b212-2f8305008dfa
# ╟─4525b42b-5421-4416-a5e7-57369d2e90cb
# ╟─ff04e9d3-2890-437f-992f-b55761ed12dd
# ╟─c67b7c75-5526-4e54-8e96-a0870b3b4d1b
# ╟─f54e0863-2dd3-4719-8f69-fee15883a924
# ╟─b6598bb0-c52b-46fe-acce-604c92e288a0
# ╠═3d38025e-8595-42cd-92b1-958a6a52ddb4
# ╠═268ba231-3d96-4c42-b89f-4406db276c27
# ╟─cfc89850-d76d-4228-9326-78ac85a0e80d
# ╟─f2273cd2-6808-4d4d-ac33-7825deffb72f
# ╟─9f1ee82c-11d6-446e-a3b7-09a0e13b8377
# ╟─501439e9-0343-4ba5-8df4-6a87b456ba5f
# ╟─133b6939-10b2-4c8e-acf8-5658ca96a0f9
# ╠═ea04744b-4296-4dc2-9a3c-1f477c96f1ac
# ╟─3d84b923-42fa-422b-a9e2-4ed26880e1e9
# ╟─4efec7f2-8a4e-4c4d-9fe9-53754fe35065
# ╟─11f2dc87-601b-4ca0-9357-7aa1458d01a6
# ╠═4a66d07c-4828-43d2-9c1a-de5475d966c9
