### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

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
	
	# include("utilities/two-qubit-operators.jl")
	include("utilities/plotting.jl")
	include("utilities/utilities.jl")

	include("notebooks/table-of-contents.jl")
	include("resources.jl")
	include("notebook-utilities.jl")
	
	md" # Packages and julia files"
end

# ╔═╡ e625f413-e93b-48c9-b14d-181c37b0fe0a
md"""
In this notebook, we model PMME using a bipartite representation, based on $Budini_2014.
"""

# ╔═╡ 249f33ea-92aa-4b9c-bdfa-d35fe5254225
mdp(table_of_contents📔)

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents(title="Fluctuator bath")

# ╔═╡ 97e21ba9-f148-46b5-81a2-53fb00c0cf6e
md"""
# Hilbert space
"""

# ╔═╡ d340648b-06ce-432e-b6e5-932ee9193ee5
begin
	"Single-qubit operators ------------------------------------------------"
	
	# Basis
	q = SpinBasis(1//2)
	Is = identityoperator(q)
	Ia = identityoperator(q)
	
	# qubit operators, using convention that |-z> is ground state (aligned with QuantumOpticsBase convention)
	σx = sigmax(q)
	σy = sigmay(q)
	σz = sigmaz(q)
	σ = sigmam(q)
	nq = σ' * σ
	
	g = spindown(q) # "minus" |-> or |a0> state in Budini
	m = g
	a0 = g
	
	e = spinup(q) # "plus" |+> or |a1> state in Budini
	p = e
	a1 = e
	
end

# ╔═╡ 5b9c647d-da6c-4f1f-8872-c70143d1195e
md"""
# Definitions
"""

# ╔═╡ d75ea611-5c86-4ff0-8b31-777b66a7914a
md"""
## Functions
"""

# ╔═╡ 2fefce63-2f69-46e2-87af-01f2b761a5ef
begin
	com(A, B) = A * B - B * A
	acom(A, B) = A * B + B * A
end

# ╔═╡ af4416fc-9c30-49a0-a8ec-9b3798cdcc75
md"""
## Operators and superoperators
"""

# ╔═╡ 0f5a2f11-838e-4f23-a286-b54e58dfcb1e
Rate

# ╔═╡ c38dce45-7573-4da5-b492-b22972cb0c40
md"""
$V_\alpha = \ket{r_\alpha} \bra{u}$ is the "system operator" correpsonding to transition to state $\ket{r_\alpha}$; in this case there is only one operator, where we choose $\ket u = \ket +$ and $\{\ket{r_\alpha}\} = \ket{-}$, leading to

$V = \sigma = \ket{-} \bra{+}.$

Then the diagonalized Lindblad superoperator is given by

$\mathcal{C}_s = \mathcal{D}_s + \mathcal{J}_s$


"""

# ╔═╡ bf9f4f40-f5e8-4175-86a7-b5f19a465b8e
md"""

where

$\mathcal{D}_s[\rho] = - \frac12 \sum_\alpha \gamma_\alpha \{V_\alpha^\dagger V_\alpha, \rho \} = - \frac12 \gamma \{\ket{u}\bra{u}, \rho \} = -\frac12 \gamma \{\ket + \bra +, \rho \}$
"""

# ╔═╡ 4e117236-b419-4f47-a2ba-b7546e262c38
md"""

is the dissipation superoperator and

$\mathcal{J}_s[\rho] = \sum_\alpha \gamma_\alpha V_\alpha \rho V_\alpha^\dagger = \gamma \overline{\rho}_s \bra u \rho \ket{u} = \gamma \bra + \rho \ket + \ket - \bra -$
"""

# ╔═╡ f9a65621-d076-48f5-b438-f8281c3d3798
md"""
is the jump superoperator, and  $\overline{\rho}_s$ is the system "resetting state":

$\overline{\rho}_s = \sum_\alpha p_\alpha \ket{r_\alpha} \bra{r_\alpha} = \ket - \bra -.$
"""

# ╔═╡ ec0d8380-3c0f-4012-934a-56fdcb5ac6b6
ρ̃ₛ = dm(m)

# ╔═╡ 8b0a2916-47c9-4289-8605-9e290f9e290d
md"""
We also define an initial state $\rho_0^s = \ket - \bra -$.
"""

# ╔═╡ bb839574-bc32-447b-9741-ab50e1ba7a15
ρ₀ˢ = dm(m)

# ╔═╡ 5e236a20-2b4c-44ed-b50f-f82a76fcfc87
md"""
We write the kernel in the Laplace domain as 

$\hat k_I(u) = \frac{\phi}{u + \phi + \varphi}.$
"""

# ╔═╡ 692a8268-bfb8-4f20-8bc9-a7790207c54b
md"""
Then the dissipation and jump superoperators in the Laplace domain become

$\hat {\mathcal{D}}_s(u) = \mathcal{D}_s \hat{k}_I(u - \mathcal L_s - \mathcal D_s)$
"""

# ╔═╡ fd604412-7d11-4177-8eb7-e6ebef6c5725
md"""
and 

$\hat {\mathcal{J}}_s(u) = \mathcal{J}_s \hat{k}_I(u - \mathcal L_s - \mathcal D_s)$
"""

# ╔═╡ 8b0aabe9-1538-48ad-adbd-c1d9479ea61d
md"""
We define the "unnormalized conditional propagator" as 

$\hat{\mathcal{T}}(u) = \frac1{u - \mathcal{L}_s - \hat{\mathcal{D}}_s(u)}$

"""

# ╔═╡ faeb3219-d9bf-4a58-a969-56d1aa9cc167
md"""
We also define the propagator

$\hat{\mathcal{G}}'(u) = \frac1{u - \mathcal{L}_s - \mathcal{C}_s \hat{k}_{\text{II}}'(u - \mathcal{L}_s)}$
"""

# ╔═╡ d6343434-3f7f-4f19-8e74-0cfc3ee5285c
md"""
with kernel 

$\hat{k}_\text{II}'(u) = 1 - \frac{\varphi}{u + \phi + \varphi}$
"""

# ╔═╡ 9d15048f-f41e-4eec-a681-130fc8f420b9
md"""
We define additional propagators and superoperators

$\hat{\mathcal{T}}'(u) = \frac1{u - \mathcal{L}_s - \hat{\mathcal{D}}_s'(u)}$
"""

# ╔═╡ 21d8de14-b36e-48ad-bbfa-68a50214cd7b
md"""
$\hat{\mathcal{D}}_s'(u) = \mathcal{D}_s \hat{k}_\text{II}'(u - \mathcal{L}_s)$
"""

# ╔═╡ 2aa1c26b-87cb-466f-8820-71c28ed95ede
md"""
$\hat{\mathcal{J}}_s'(u) = \mathcal{J}_s \hat{k}_\text{II}'(u - \mathcal{L}_s)$
"""

# ╔═╡ f1385a4f-5aeb-4364-8466-4de1390bb44f
md"""
These allow us to calculate the waiting time probability density 

$\hat{w}(u) = \text{Tr}_s [\hat{\mathcal{J'}}_s(u) \hat{\mathcal{T'}}(u)\overline{\rho}_s]$
"""

# ╔═╡ 2d5008a3-aa5b-432b-a071-7590c9421a49
md"""
and the initial waiting time probability density

$\hat{w}_\text{in}(u) = \text{Tr}_s [\hat{\mathcal{J}}_s(u) \hat{\mathcal{T}}(u)\rho_0^s]$
"""

# ╔═╡ 7e098d1e-89d4-489f-9e98-7e362890c858
md"""
**Will need to multiplex functions to take QOp or frequency inputs.**
"""

# ╔═╡ 045fa364-8be3-4df7-9ff3-22483c15b60e
md"""
## Parameters
"""

# ╔═╡ 0cc3cda9-3fd6-4549-9c62-148137f239a0
begin
	γ = 1
	Ω = 0.15
	φ = 0.01
	ϕ = 0.01	
	ωs = 1 		# system energy splitting
end

# ╔═╡ 0587d065-64aa-4a24-a8e5-45710247f67b
begin
	ρ0s = dm(m) 		# initial system state
	ρ0sa = ρ0s ⊗ dm(a0) # initial system-ancilla state

	# operators
	A = Is ⊗ σ' 		# ancilla operator
	T = σ ⊗ dm(a1) 		# system decay operator (coupling term)
	Hₛ = (ωs/2) * σz 	# system Hamiltonian

	# superoperators
	Lₛ(ρ) = -im * (Ω/2) * com(σx, ρ)
	I(ρ) = ρ # identity
	# Cs(ρ) = (γ/2) * (com(σ, ρ * σ') + com(σ * ρ, σ'))

	# density matrix evolution
	ddt_ρtsa(ρ) = -im * (Ω/2) * com(σx ⊗ Ia, ρ) + (γ/2) * (com(T, ρ * T') + com(T * ρ, T')) + (ϕ/2) * (com(A, ρ * A') + com(A * ρ, A')) + (φ/2) * (com(A', ρ * A) + com(A' * ρ, A))

end

# ╔═╡ aa34e59a-3eea-4493-a348-fb5f5459aff9
Dₛ(ρ) = -0.5 * γ * acom(dm(p), ρ)

# ╔═╡ e8c731b9-2e7e-477b-93d5-6dbf11458330
Jₛ(ρ) = γ * (p' * ρ * p) * dm(m)

# ╔═╡ b6f8035c-b023-4da4-96a9-85a67d2bffe0
Cₛ(ρ) = Dₛ(ρ) + Jₛ(ρ)

# ╔═╡ 3b1a73b2-7508-4ae3-bf3b-861cf5e46374
begin
	k̂₁(u::Rate) = ρ -> I(ρ) *  (ϕ / (u + ϕ + φ))
	k̂₁(u::QOp) = ρ -> I(ρ) * ϕ / (u + I(ρ) * (ϕ + φ))
end

# ╔═╡ 454b43fb-7045-4f50-bd3b-377c58e237a0
begin
	D̂ₛ(u::Rate) = ρ -> Dₛ(ρ) * k̂₁(u * I(ρ) - Lₛ(ρ) - Dₛ(ρ))
	D̂ₛ(u::QOp) = ρ -> Dₛ(ρ) * k̂₁(u - Lₛ(ρ) - Dₛ(ρ))
end

# ╔═╡ 41a073ea-5de2-4f6d-9c55-7c0e3abaa537
begin
	T̂(u::Rate) = ρ -> 1/(u * I(ρ) - Lₛ(ρ) - D̂ₛ(u)(ρ))
	T̂(u::QOp) = ρ -> 1/(u - Lₛ(ρ) - D̂ₛ(u)(ρ))
end

# ╔═╡ ea360261-aa88-43b4-9037-84e83e067b57
begin
	Ĵₛ(u::Rate) = ρ -> Jₛ(ρ) * k̂₁(u * I(ρ) - Lₛ(ρ) - Dₛ(ρ))
	Ĵₛ(u::QOp) = ρ -> Jₛ(ρ) * k̂₁(u - Lₛ(ρ) - Dₛ(ρ))
end

# ╔═╡ e6a1a375-b267-4370-9db7-1ab5770bb126
ŵᵢ(u) = tr(Ĵₛ(u)(T̂(u)(ρ₀ˢ)))

# ╔═╡ 5de56302-75ca-42d7-9ec1-20209c550571
begin
	k̂₂(u::Rate) = ρ -> I(ρ) * (1 - φ/(u + ϕ + φ))
	k̂₂(u::QOp) = ρ -> I(ρ) - I(ρ) * φ / (u + I(ρ) * (ϕ + φ))
end

# ╔═╡ bb7b2dc4-4893-4f8b-a49c-2a8623edeceb
begin
	Ĝp(u::Rate) = ρ -> 1/(u * I(ρ) - Lₛ(ρ) - Cₛ(ρ) * k̂₂(u - Lₛ(ρ)))
	Ĝp(u::QOp) = ρ -> 1/(u - Lₛ(ρ) - Cₛ(ρ) * k̂₂(u - Lₛ(ρ)))
end

# ╔═╡ f46c322c-f592-40b2-a8dc-51a8e5c76b1a
begin
	D̂ₛp(u::Rate) = ρ -> Dₛ(ρ) * k̂₂(u * I(ρ) - Lₛ(ρ))
	D̂ₛp(u::QOp) = ρ -> Dₛ(ρ) * k̂₂(u - Lₛ(ρ))
end

# ╔═╡ 19f748ae-f947-4073-b787-359fbe25fcbe
begin
	T̂p(u::Rate) = ρ -> 1 / (u * I(ρ) - Lₛ(ρ) - D̂ₛp(u)(ρ))
	T̂p(u::QOp) = ρ -> 1 / (u - Lₛ(ρ) - D̂ₛp(u)(ρ))
end

# ╔═╡ 44a1a11b-40bd-4ef8-a737-8ffbfc5a6dfb
begin
	Ĵₛp(u::Rate) = ρ -> Jₛ(ρ) * k̂₂(u * I(ρ) - Lₛ(ρ))
	Ĵₛp(u::QOp) = ρ -> Jₛ(ρ) * k̂₂(u - Lₛ(ρ))
end

# ╔═╡ a20cc119-e99b-4bfc-8db3-4b5a93d47d00
ŵ(u::Rate) = tr(Ĵₛp(u)(T̂p(u)(ρ̃ₛ)))

# ╔═╡ 9c48ad25-ba45-4b52-a5c5-acbf34526e89
ŵ(0.5)

# ╔═╡ 0f5189fb-f05c-4aa5-98ae-218b13acb4b4
md"""
# Utilities
"""

# ╔═╡ 1b7e3334-9825-4832-a0b8-d58816fd5236
begin
	@userplot PlotFluctuators
	@recipe function f(pf::PlotFluctuators)
		ωi, times, Hi, Ω = pf.args
	
		labels = permutedims(map(ω -> string(ω, " MHz"), ωi))
	
		# Plot time series --------------------------------------------------------
	
		xlabel --> "t (μs)"
		ylabel --> "fluctuator value"
		legend --> :outerright
		label --> labels
	
		palette := :tab10
		linealpha --> 1
	
		legendfontsize --> 6
		titlefontsize --> 12
		xtickfontsize --> 8
		ytickfontsize --> 8
		xguidefontsize --> 10
		yguidefontsize --> 10
		size --> (600,300)
		linewidth --> 1.5
	
		for i in reverse(1:length(ωi))
	
			fluctuators = map(t -> real(expect(Hi(t)[i], σz)) / Ω, times)
	
			@series begin
				times, fluctuators
			end
	
		end
	end
end

# ╔═╡ 2e072a88-7166-4d5d-96e6-52fcbe668505
function bloch_plots(sols::Vector{Solution}; alpha=0.1, N=50, kwargs...)
	colors = palette(:tab10) 
	
	# calculate expectation values --------------------------------------------
	t = sols[1].t
	xs, ys, zs = [], [], []
	
	for sol in sols
		x, y, z = map(op -> expectations(sol, op), qbasis)
		for (list, traj) in zip([xs, ys, zs], [x, y, z])
			push!(list, traj)
		end
		
	end
	
	
	# plot ----------------------------------------------------------------------
	function bloch(os; color=colors1q[1], xlabel="", ylabel="")
		
		po = plot(ylims = [-1,1], xlabel=xlabel, ylabel=ylabel)
		
		for o in os[1:min(N, 50)]
			plot!(t, o, alpha=alpha, label=:none, color=color)
		end
		
		oavg = [mean([os[i][j] for i in 1:N]) for j in 1:length(t)]
		plot!(t, oavg, alpha=1, color=color, linewidth=3)

		po
		
	end

	
	l = @layout [xplot{0.33h}; yplot{0.33h}; zplot{0.33h}]
	
	px = bloch(xs, color=colors[1], ylabel="x")
	py = bloch(ys, color=colors[2], ylabel="y")
	pz = bloch(zs, color=colors[3], ylabel="z", xlabel="t (μs)")
	
	plot(px, py, pz, layout = l, link=:y, size=(800,500), legendfontsize=8, titlefontsize=12, legend=:none, kwargs...)
	
end

# ╔═╡ Cell order:
# ╟─e625f413-e93b-48c9-b14d-181c37b0fe0a
# ╟─249f33ea-92aa-4b9c-bdfa-d35fe5254225
# ╟─3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─97e21ba9-f148-46b5-81a2-53fb00c0cf6e
# ╠═d340648b-06ce-432e-b6e5-932ee9193ee5
# ╟─5b9c647d-da6c-4f1f-8872-c70143d1195e
# ╟─d75ea611-5c86-4ff0-8b31-777b66a7914a
# ╠═2fefce63-2f69-46e2-87af-01f2b761a5ef
# ╟─af4416fc-9c30-49a0-a8ec-9b3798cdcc75
# ╠═0f5a2f11-838e-4f23-a286-b54e58dfcb1e
# ╠═0587d065-64aa-4a24-a8e5-45710247f67b
# ╟─c38dce45-7573-4da5-b492-b22972cb0c40
# ╟─b6f8035c-b023-4da4-96a9-85a67d2bffe0
# ╟─bf9f4f40-f5e8-4175-86a7-b5f19a465b8e
# ╠═aa34e59a-3eea-4493-a348-fb5f5459aff9
# ╟─4e117236-b419-4f47-a2ba-b7546e262c38
# ╠═e8c731b9-2e7e-477b-93d5-6dbf11458330
# ╟─f9a65621-d076-48f5-b438-f8281c3d3798
# ╠═ec0d8380-3c0f-4012-934a-56fdcb5ac6b6
# ╟─8b0a2916-47c9-4289-8605-9e290f9e290d
# ╠═bb839574-bc32-447b-9741-ab50e1ba7a15
# ╟─5e236a20-2b4c-44ed-b50f-f82a76fcfc87
# ╠═3b1a73b2-7508-4ae3-bf3b-861cf5e46374
# ╟─692a8268-bfb8-4f20-8bc9-a7790207c54b
# ╠═454b43fb-7045-4f50-bd3b-377c58e237a0
# ╟─fd604412-7d11-4177-8eb7-e6ebef6c5725
# ╠═ea360261-aa88-43b4-9037-84e83e067b57
# ╟─8b0aabe9-1538-48ad-adbd-c1d9479ea61d
# ╠═41a073ea-5de2-4f6d-9c55-7c0e3abaa537
# ╟─faeb3219-d9bf-4a58-a969-56d1aa9cc167
# ╠═bb7b2dc4-4893-4f8b-a49c-2a8623edeceb
# ╟─d6343434-3f7f-4f19-8e74-0cfc3ee5285c
# ╠═5de56302-75ca-42d7-9ec1-20209c550571
# ╟─9d15048f-f41e-4eec-a681-130fc8f420b9
# ╠═19f748ae-f947-4073-b787-359fbe25fcbe
# ╟─21d8de14-b36e-48ad-bbfa-68a50214cd7b
# ╠═f46c322c-f592-40b2-a8dc-51a8e5c76b1a
# ╟─2aa1c26b-87cb-466f-8820-71c28ed95ede
# ╠═44a1a11b-40bd-4ef8-a737-8ffbfc5a6dfb
# ╟─f1385a4f-5aeb-4364-8466-4de1390bb44f
# ╠═a20cc119-e99b-4bfc-8db3-4b5a93d47d00
# ╟─2d5008a3-aa5b-432b-a071-7590c9421a49
# ╟─e6a1a375-b267-4370-9db7-1ab5770bb126
# ╠═9c48ad25-ba45-4b52-a5c5-acbf34526e89
# ╟─7e098d1e-89d4-489f-9e98-7e362890c858
# ╟─045fa364-8be3-4df7-9ff3-22483c15b60e
# ╠═0cc3cda9-3fd6-4549-9c62-148137f239a0
# ╟─0f5189fb-f05c-4aa5-98ae-218b13acb4b4
# ╠═1b7e3334-9825-4832-a0b8-d58816fd5236
# ╠═2e072a88-7166-4d5d-96e6-52fcbe668505
# ╠═df539f73-7c8a-4d65-b909-0d94720f724f
