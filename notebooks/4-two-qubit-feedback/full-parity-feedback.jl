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

# ╔═╡ 897b9f97-3dab-4561-afed-be6e67efe047
begin
	directory_name = "QC-notebooks"
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

	using DataFrames
	using CSV
	
	include("utilities/two-qubit-operators.jl")
	include("utilities/plotting.jl")
	include("utilities/utilities.jl")

	include("notebooks/table-of-contents.jl")
	include("resources.jl")

	include("notebook-utilities.jl")
	
	md" # Packages and julia files"
end

# ╔═╡ ed642b58-e611-4ca2-a491-4d43f9f7ce29
TableOfContents(title="Full parity feedback")

# ╔═╡ 48faa1c1-c7b1-47e9-8e07-11cc57b1c0e9
md"""
# Setup
"""

# ╔═╡ 2931a779-d095-49a5-96bb-52e46573e266
md"""
## System parameters
"""

# ╔═╡ 7cb94127-7392-4953-a488-5d9fa573f372
begin
	ΓH = 0.5
	ΓZ = 40.0
	dt = 1e-3
	tf = 20.0
	td = 0.0
	tZ = 10.0
	Ω = (2π) * 0.2
	Ωmax = (2π) * 50.0
	ηH = 1.0
	ηZ = 1.0
	nfluctuators = 10
	ρ0 = normalize((g + e) ⊗ (g + e))
end

# ╔═╡ 1be3fe91-2f5d-4a93-bc16-0e8acf2eaa02
md"""
## Setting up the fluctuators
"""

# ╔═╡ 805fc3ba-55bd-467c-abda-21556c306c92
begin
	# sample log-uniformly
	(ωmin, ωmax) = 1/tf, 1/dt
	ωlog = range(log(ωmin), log(ωmax), length=nfluctuators)
	ωi = exp.(ωlog)
	τi = map(ω -> 1/ω, ωi) # fluctuators and their signs for all $nfluctuators switching times
end

# ╔═╡ 305ce5e8-b73f-40d5-9965-1316c6b59635
md"""
### Fluctuator Hamiltonians
"""

# ╔═╡ 384d4c66-63d4-4a7d-9092-83b4b399b8d9
function sgn(τ::Timescale, tf::Timescale)
	nflips = Int64(floor(tf / τ)) # number of times fluctuator CAN flip in simulation
	series = rand((-1, 1), (nflips + 1)) # string of random signs of fluctuator

	t -> let
		index = Int64(floor(t/tf * nflips)) + 1
		return series[index]
	end
end

# ╔═╡ 7b684ee3-e7c1-4b1a-b412-1a3e87a081f2
begin
	Hf1(τi, tf) = let si = map(τ -> sgn(τ, tf), τi)
					t -> map(s -> (Ω/2) * σz1 * s(t), si) end

	Hf2(τi, tf) = let si = map(τ -> sgn(τ, tf), τi)
					t -> map(s -> (Ω/2) * σz2 * s(t), si) end

	# total fluctuator Hamiltonian -- a function of a particular realization of fluctuator Hamiltonian Hf
	Hz(Hf1, Hf2) = t -> sum(Hf1(t)) + sum(Hf2(t))
end

# ╔═╡ 8c4c6702-c37f-46ab-b24f-60415bd512e2
md"""
## Feedback Hamiltonians
"""

# ╔═╡ b68e9606-fab7-4689-ad7b-daccc8ed6135
md"""
### Helpers
"""

# ╔═╡ e2a21ce6-e082-44bd-9d50-d7af2144af75
function arctan(y, x)
	if x > 0
		atan(y / x)
	elseif y >= 0 && x < 0
		atan(y / x) + π
	elseif y < 0 && x < 0
		atan(y / x) - π
	else
		sign(y) * π/2
	end
end

# ╔═╡ eb52e9d1-6b55-403c-90c5-778eaed71817
function θ(ρ::State)
	tm = g ⊗ g
	tp = e ⊗ e
	s = Ψm
	t0 = Ψp

	tmt0 = projector(tm, dagger(t0))
	t0tp = projector(t0, dagger(tp))
	ss = projector(s)
	tmtp = projector(tm, dagger(tp))

	F = expect(dm(t0), ρ)
	ρtmt0 = expect(tmt0, ρ)
	ρt0tp = expect(t0tp, ρ)
	ρss = expect(ss, ρ)
	ρtmtp = expect(tmtp, ρ)

	y = real(√8 * (ρtmt0 - ρt0tp))
	x = real(3F + ρss + 2ρtmtp - 1)

	return 0.5 * arctan(y, x)
end

# ╔═╡ 133fc9dc-9b77-4c44-af21-6cd083f2f467
magθ(ρ) = 	let inst = - θ(ρ) / (2dt)
				(abs(inst) < Ωmax) ? inst : sign(inst) * Ωmax
			end

# ╔═╡ 65d796cc-80c8-47f5-a594-3cf44c83079a
md"""
### Excitation feedback Hamiltonian
"""

# ╔═╡ dc824332-854e-4035-af37-e24e87c85db5
HF(t::Timescale, ρ::State) = (magθ(ρ)/2) * (σy1 + σy2)

# ╔═╡ d7b48569-b67f-4390-948d-45c23a9538bf
md"""
### Parity feedback Hamiltonian
"""

# ╔═╡ 31eff151-0ce3-4aa6-8f93-6bd9e969d610
threshold = 0.05

# ╔═╡ 4c455a91-412d-4e7a-9308-db3c90248537
md"""
## Total Hamiltonians (system and filter)
"""

# ╔═╡ 91d24100-8268-461e-a585-f9bb36f89d7b
md"""
## Measurement
"""

# ╔═╡ 61bf8483-a10f-4140-b50f-2351ecce69e1
begin
	c = σx1 * σx2
	cZ(t) = (t < tZ) ? Iq ⊗ Iq : c
	
	C = [	((σz1 + σz2)/√2, ΓH, ηH), 	# excitation number measurement for feedback
					(cZ, ΓZ, ηZ)] 				# (half) parity measurement for Zeno pinning
	
	
	J = []	# Lindblad dissipation due to inefficient measurement
	for (c, Γ, η) in C
		if η < 1.0
			push!(J, (c, (1 - η) * Γ))
		end
	end
end

# ╔═╡ 41122eaa-dbbe-11ec-0cd6-fb171d50e120
md"""
# Simulation
"""

# ╔═╡ f2a9c32d-17a9-4dbf-9e21-6b7c3ba1e3c0
s=7

# ╔═╡ ef6f8efb-76e6-4624-b332-156dcbfd23cb
md"""
ΩZ = 0.0 MHz $(@bind ΩZ Slider((2π).*(0.0:0.1:2.0))) 2.0 MHz
"""

# ╔═╡ 9c487fcf-b1e2-4815-94aa-14f0fb0b5220
# HZ(t::Timescale, r::Record) = (ΩZ/2) * r[2] * σz1
HZ(t::Timescale, r::Record) = t < tZ ? Iq ⊗ Iq : (-ΩZ/2) * (r[2] - 1) * σz1
# HZ(t::Timescale, ρ::State) = 	if t > tZ && (real(expect(dm(Ψp), ρ)) < threshold)
# 									(π/2dt) * σz1
# 								else
# 									Iq ⊗ Iq
# 								end

# ╔═╡ b07c9329-9e78-4ed7-b7ef-78705cfcbccb
begin
	Hfil(t::Timescale, ρ::State, r::Record) = HF(t, ρ) + HZ(t, r) # filter
	Hsys(Hf1, Hf2) = (t::Timescale, ρ::State, r::Record) -> Hfil(t, ρ, r) + Hz(Hf1, Hf2)(t) # system
end

# ╔═╡ 0384ca43-9ea1-4aff-9514-1fc79c1fdfaa
ΩZ/2π

# ╔═╡ 87c9e502-d0ed-4039-a706-0e737ca19132
readout_integration_bins = 500

# ╔═╡ 93dc4b96-0dd1-4177-81d5-179789a135b5
begin
	Random.seed!(s)
	Hi1 = Hf1(τi, tf) # new realization of fluctuator Hamiltonian (single qubit)
	Hi2 = Hf1(τi, tf) # new realization of fluctuator Hamiltonian (single qubit)

	sol = bayesian((0.0, tf), ρ0, (Hsys(Hi1, Hi2), Hfil), J, C; dt=dt, readout_integration_bins = readout_integration_bins)
end

# ╔═╡ 8e37e639-069d-4a28-b71d-37a2340b4b70
function trans(mat::Vector{<:Vector{<:Real}})
	nrows = length(mat)
	ncols = length(mat[1])
	return [[mat[i][j] for i in 1:nrows] for j in 1:ncols]
end

# ╔═╡ 437075b9-d429-48e7-a18c-f8cf3c6e95ad
md"""
## Ensemble
"""

# ╔═╡ 67d3e834-5352-47f1-a30c-16fc7a3671b1
N = 100

# ╔═╡ 97936682-9e5e-4dac-9e44-4d5dcba9a65a
# begin
# 	sols = map(1:N) do n
		
# 		Hi1 = Hf1(τi, tf) # new realization of fluctuator Hamiltonian (single qubit)
# 		Hi2 = Hf1(τi, tf) # new realization of fluctuator Hamiltonian (single qubit)

# 		bayesian((0.0, tf), ρ0, (Hsys(Hi1, Hi2), Hfil), J, C; dt=dt, readout_integration_bins = readout_integration_bins)

# 	end
# end

# ╔═╡ 3a728387-32d9-4c53-9669-a546e231dd3e
times = collect(0.0:dt:tf)

# ╔═╡ 51c409d4-37bc-4e2c-8b55-f17d1add60a7
typeof(times)

# ╔═╡ 135cadb3-5f46-4120-bf5f-50b1685dc808
ΩZ/2π

# ╔═╡ 50fe781c-40fe-47c8-8235-adf28309ec13
function bell_plot(exps::Vector{Vector{Float64}}, t::Vector{Float64})

	basis = bell_basis
	colors = palette(:rainbow)
	labels = bell_basis_labels
	title = "Bell states"
	
	pl = plot(size=(600,300), legendfontsize=12, titlefontsize=12, legend=:outerright, title=title)
	
	for l in 1:length(basis)
		label = labels[l]
		color = colors[l]
		exp = exps[l]
		plot!(t, exps[l], color=color, label=label, legend=:outerright, xlabel="t (μs)", ylims=[0,1])
	end

	pl

end

# ╔═╡ ef6c09c8-eef4-47be-b114-9217f3124bd5
bell_plot(sol...)

# ╔═╡ 0ed67bd2-1eba-4128-b611-87a85bd2cba2
bp = bell_plot(sol[1])

# ╔═╡ 9d46b037-a439-4b17-8d7b-2d9b83a16a81
begin
	rp = plot(sol[1].t,coarse_grain(sol[1].r[2], n = readout_integration_bins) .-1, legend=:outerright, xlabel="t (μs)", ylabel="smoothed record", label="record  ", ylims=[-4,1])
	plot(bp, rp, layout=grid(2,1))
end

# ╔═╡ cb1d93f0-4792-4720-be61-bad0ba4c0d23
bell_plot(sol[2])

# ╔═╡ b22495a7-1733-48eb-ad3e-8b0d0560408c
function ensemble_average(solutions::Vector{Solution})

	if typeof(solutions[1].ρ[1]) <: Ket
		return [mean(map(sol -> dm(sol.ρ[i]), solutions)) for i in 1:length(solutions[1].t)]
	else
		return [mean(map(sol -> sol.ρ[i], solutions)) for i in 1:length(solutions[1].t)]
	end

end

# ╔═╡ 57ed4275-eba1-4924-8a62-853d6670eee2
begin
	ops = bell_basis
	op_labels = bell_basis_strings
	
	sys_sols = map(sol -> sol[1], sols)
	fil_sols = map(sol -> sol[2], sols)
	
	sys_avg = ensemble_average(sys_sols)
	fil_avg = ensemble_average(fil_sols)

	fids_sys = map(ρ -> fidelity(ρ, Ψp), sys_avg)
	exps_sys = [map(ρ -> real(expect(ρ, op)), sys_avg) for op in ops]

	fids_fil = map(ρ -> fidelity(ρ, Ψp), fil_avg)
	exps_fil = [map(ρ -> real(expect(ρ, op)), fil_avg) for op in ops]

end

# ╔═╡ db928ce6-6286-4ddf-b242-599312a91a3c
bell_plot(exps_sys, times)

# ╔═╡ d163d238-4653-4cb1-be7e-0dd5a96e1b54
exps_sys

# ╔═╡ 4e7ebe7d-4924-4c3e-953b-c56d1b99023f
md"""
# Testing
"""

# ╔═╡ a90ca242-d51f-49db-a459-92ea10d84496
begin
	const tt = 0.0 # Timescale
	const rr = [0.0] # Readout
	const ρρ = identityoperator(SpinBasis(1//2)) # State
end

# ╔═╡ a1b640f2-35d5-4527-98cc-ae0f27ffbcf1


# ╔═╡ bd469f04-85a8-4b1a-be45-715fce80f522
arglist = [(tt, ρρ, rr), (tt, ρρ, ρρ), (tt, ρρ, ρρ, rr)]

# ╔═╡ 21a0079b-8776-4864-a4a1-242fd8d01f66
function ham(dt, H::QOp)
	u::QOp = exp( -im * dt * DenseOperator(H))
	(t, ρ) -> update(u, ρ)
end

# ╔═╡ 39d888ec-213f-41db-988b-000da9a8c895
function ham(dt, H::Function)
	# evaluate H for specific inputs (using test constants tt::Timescale, rr::Record, ρρ::State)
	if applicable(H, tt, rr)
		(t::Timescale, ρ::State, r::Readout) -> ham(dt, H(t, r))(t, ρ)

	elseif applicable(H, tt, ρρ)
		(t::Timescale, ρ::State, ρd::State) -> ham(dt, H(t, ρd))(t, ρ)

	elseif applicable(H, tt, ρρ, rr)
		(t::Timescale, ρ::State, ρd::State, r::Readout) -> ham(dt, H(t, ρd, r))(t, ρ)

	elseif applicable(H, tt)
		(t::Timescale, ρ::State) -> ham(dt, H(t))(t, ρ)

	else
		message = string("Hamiltonian does not have correct argument types. Define a function ",
							"H(t::Float64, r::Vector{Float64}) for readout feedback, H(t::Float64, ",
							"ρ::AbstractOperator) for state feedback, or H(t::Float64) for simple time dependence.")
		error(message)
	end
end

# ╔═╡ 6944b675-fe92-4450-abc2-d001f78bc7d2
hf = ham(dt, Hfil)

# ╔═╡ 8efad48d-7170-4945-a3f8-5fe2a3da616e
function apply(h::Function, t::Timescale, ρ::State, ρd::State, rd::Record)
	arglist = [(t, ρ, rd), (tt, ρ, ρd), (tt, ρ, ρd, rd)]
	index = findfirst(map(args -> applicable(hf, args...), arglist))
	args = arglist[index]

	return h(args...)
end

# ╔═╡ 901bac3f-6bc4-4416-a307-f5061abcf947
md"""
# Utilities
"""

# ╔═╡ Cell order:
# ╟─ed642b58-e611-4ca2-a491-4d43f9f7ce29
# ╟─48faa1c1-c7b1-47e9-8e07-11cc57b1c0e9
# ╟─2931a779-d095-49a5-96bb-52e46573e266
# ╠═7cb94127-7392-4953-a488-5d9fa573f372
# ╟─1be3fe91-2f5d-4a93-bc16-0e8acf2eaa02
# ╠═805fc3ba-55bd-467c-abda-21556c306c92
# ╟─305ce5e8-b73f-40d5-9965-1316c6b59635
# ╠═384d4c66-63d4-4a7d-9092-83b4b399b8d9
# ╠═7b684ee3-e7c1-4b1a-b412-1a3e87a081f2
# ╟─8c4c6702-c37f-46ab-b24f-60415bd512e2
# ╟─b68e9606-fab7-4689-ad7b-daccc8ed6135
# ╠═e2a21ce6-e082-44bd-9d50-d7af2144af75
# ╠═eb52e9d1-6b55-403c-90c5-778eaed71817
# ╠═133fc9dc-9b77-4c44-af21-6cd083f2f467
# ╟─65d796cc-80c8-47f5-a594-3cf44c83079a
# ╠═dc824332-854e-4035-af37-e24e87c85db5
# ╟─d7b48569-b67f-4390-948d-45c23a9538bf
# ╠═31eff151-0ce3-4aa6-8f93-6bd9e969d610
# ╠═9c487fcf-b1e2-4815-94aa-14f0fb0b5220
# ╟─4c455a91-412d-4e7a-9308-db3c90248537
# ╠═b07c9329-9e78-4ed7-b7ef-78705cfcbccb
# ╟─91d24100-8268-461e-a585-f9bb36f89d7b
# ╠═61bf8483-a10f-4140-b50f-2351ecce69e1
# ╟─41122eaa-dbbe-11ec-0cd6-fb171d50e120
# ╠═93dc4b96-0dd1-4177-81d5-179789a135b5
# ╠═ef6c09c8-eef4-47be-b114-9217f3124bd5
# ╠═0ed67bd2-1eba-4128-b611-87a85bd2cba2
# ╠═f2a9c32d-17a9-4dbf-9e21-6b7c3ba1e3c0
# ╠═ef6f8efb-76e6-4624-b332-156dcbfd23cb
# ╠═0384ca43-9ea1-4aff-9514-1fc79c1fdfaa
# ╠═87c9e502-d0ed-4039-a706-0e737ca19132
# ╠═9d46b037-a439-4b17-8d7b-2d9b83a16a81
# ╠═8e37e639-069d-4a28-b71d-37a2340b4b70
# ╠═cb1d93f0-4792-4720-be61-bad0ba4c0d23
# ╟─437075b9-d429-48e7-a18c-f8cf3c6e95ad
# ╠═67d3e834-5352-47f1-a30c-16fc7a3671b1
# ╠═97936682-9e5e-4dac-9e44-4d5dcba9a65a
# ╠═57ed4275-eba1-4924-8a62-853d6670eee2
# ╠═3a728387-32d9-4c53-9669-a546e231dd3e
# ╠═51c409d4-37bc-4e2c-8b55-f17d1add60a7
# ╠═db928ce6-6286-4ddf-b242-599312a91a3c
# ╠═135cadb3-5f46-4120-bf5f-50b1685dc808
# ╠═50fe781c-40fe-47c8-8235-adf28309ec13
# ╠═d163d238-4653-4cb1-be7e-0dd5a96e1b54
# ╠═b22495a7-1733-48eb-ad3e-8b0d0560408c
# ╟─4e7ebe7d-4924-4c3e-953b-c56d1b99023f
# ╠═a90ca242-d51f-49db-a459-92ea10d84496
# ╠═6944b675-fe92-4450-abc2-d001f78bc7d2
# ╠═a1b640f2-35d5-4527-98cc-ae0f27ffbcf1
# ╠═bd469f04-85a8-4b1a-be45-715fce80f522
# ╠═8efad48d-7170-4945-a3f8-5fe2a3da616e
# ╠═21a0079b-8776-4864-a4a1-242fd8d01f66
# ╠═39d888ec-213f-41db-988b-000da9a8c895
# ╠═901bac3f-6bc4-4416-a307-f5061abcf947
# ╠═897b9f97-3dab-4561-afed-be6e67efe047
