### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

# ╔═╡ df539f73-7c8a-4d65-b909-0d94720f724f
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

# ╔═╡ e625f413-e93b-48c9-b14d-181c37b0fe0a
md"""
In this notebook, we model dephasing using a fluctuator bath and look at how the counterrotating drive interacts with this.
"""

# ╔═╡ 249f33ea-92aa-4b9c-bdfa-d35fe5254225
mdp(table_of_contents📔)

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents(title="Fluctuator bath")

# ╔═╡ be587838-24f3-4b4b-8a51-789168089aa7
md"""
# Setup
"""

# ╔═╡ bb91d819-871b-48e3-92fa-7a494a67c836
md"""
## Parameters
"""

# ╔═╡ d172ea75-59de-475f-9c0f-365cc1f09b85
begin	
	N = 1000
	nfluctuators = 10
	
	ΓH = 0.5
	dt = 1e-3
	tf = 20.0
	td = 0.0
	Ω = (2π) * 0.2

	ρ0 = normalize((g + e) ⊗ (g + e))#g ⊗ g
end

# ╔═╡ cf4feaf2-0fa6-449c-86d1-ac988b5df569
md"""
## Setting up the fluctuators
"""

# ╔═╡ 3cddac9b-5ab2-4120-9355-00411c76dce3
begin
	# sample frequencies log-uniformly
	(ωmin, ωmax) = 1/tf, 1/dt
	ωlog = range(log(ωmin), log(ωmax), length=nfluctuators)
	ωi = exp.(ωlog)
	
	# fluctuators and their signs for all $nfluctuators switching times
	τi = map(ω -> 1/ω, ωi) 
end

# ╔═╡ 41dc7c32-dc78-41ea-b48a-3a99415da603
function sgn(τ::Timescale, tf::Timescale)
	nflips = Int64(floor(tf / τ)) # number of times fluctuator CAN flip in simulation
	series = rand((-1, 1), (nflips + 1)) # string of random signs of fluctuator

	t -> let
		index = Int64(floor(t/tf * nflips)) + 1
		return series[index]
	end
end

# ╔═╡ 5c1c9b7c-a16d-453c-966c-55d31ef64066
md"""
## Hamiltonians
"""

# ╔═╡ dabf9a58-907f-43ac-a45f-9fa902536f6a
md"""
### Fluctuator Hamiltonians
"""

# ╔═╡ a7d7312f-3dab-4f10-8509-2b8f90d8074f
begin
	Hf1(τi, tf) = let si = map(τ -> sgn(τ, tf), τi)
					t -> map(s -> (Ω/2) * σz1 * s(t), si) end
	
	Hf2(τi, tf) = let si = map(τ -> sgn(τ, tf), τi)
					t -> map(s -> (Ω/2) * σz2 * s(t), si) end

	# total Hamiltonian 
	# a function of a particular realization of fluctuator Hamiltonian Hf
	Hz(Hf1, Hf2) = t -> sum(Hf1(t)) + sum(Hf2(t))
end

# ╔═╡ c81d1356-120b-49dd-947e-e66ee6952a21
md"""
### Draining (counterrotating) Hamiltonian
"""

# ╔═╡ 54a89faa-4c16-46f4-8852-aa363b3fa101
md"""
The optimal angles for the feedback (rotating) and draining (counterrotating) drives are, respectively,

$\theta_{\text{opt}} = \text{atan}(2\text{Re} [\rho_{\psi+,\phi-}], \rho_{\psi+,\psi+} - \rho_{\phi-,\phi-})/2$

$\Delta \theta_{\text{opt}} = \text{atan}(\rho_{\phi+,\phi+} - \rho_{\psi-,\psi-}, 2\text{Re} [\rho_{\phi+,\psi-}])/2$
"""

# ╔═╡ e07021db-e4d4-400a-9b83-61e0ffd19159
begin
	ψpψp = projector(Ψp)
	ϕmϕm = projector(Φm)
	# ψpϕm = projector(Ψp, dagger(Φm))
	ϕmψp = projector(Φm, dagger(Ψp))
	
	ϕpϕp = projector(Φp)
	ψmψm = projector(Ψm)
	# ϕpψm = projector(Φp, dagger(Ψm))
	ψmϕp = projector(Ψm, dagger(Φp))

	md" 🚩 projectors"
end

# ╔═╡ a33e0d94-d99b-4643-8e6a-d4136237d48c
md"""

with the arctangent defined in the appropriate quadrant as:

$\arctan[y,x] = \begin{cases}
		\arctan[y/x], & x > 0 \\
		\arctan[y/x] + \pi, & y \geq 0, x < 0 \\
		\arctan[y/x] - \pi, & y < 0, x < 0 \\
		\text{sign}(y) \pi/2, & x = 0.
	\end{cases}$
"""

# ╔═╡ 5bc80650-4fb3-41ed-bc4b-15f51f937ed5
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

# ╔═╡ 5ae2203c-4a4a-4265-932e-be0deb01c97e
function θ(ρ::State)
	ρψpψp = expect(ψpψp, ρ)
	ρϕmϕm = expect(ϕmϕm, ρ)
	ρψpϕm = expect(ϕmψp, ρ) #expect(ψpϕm, ρ)

	x = real(ρψpψp - ρϕmϕm)
	y = 2 * real(ρψpϕm)

	return 0.5 * arctan(y, x)
end

# function θ(ρ::State)
# 	tm = g ⊗ g
# 	tp = e ⊗ e
# 	s = Ψm
# 	t0 = Ψp

# 	tmt0 = projector(tm, dagger(t0))
# 	t0tp = projector(t0, dagger(tp))
# 	ss = projector(s)
# 	tmtp = projector(tm, dagger(tp))

# 	F = expect(dm(t0), ρ)
# 	ρtmt0 = expect(tmt0, ρ)
# 	ρt0tp = expect(t0tp, ρ)
# 	ρss = expect(ss, ρ)
# 	ρtmtp = expect(tmtp, ρ)

# 	y = real(√8 * (ρtmt0 - ρt0tp))
# 	x = real(3F + ρss + 2ρtmtp - 1)

# 	return 0.5 * arctan(y, x)
# end

# ╔═╡ 8aab4b6c-2ce4-4667-a122-1fa4d8ae37e7
function Δθ(ρ::State)
	ρϕpϕp = expect(ϕpϕp, ρ)
	ρψmψm = expect(ψmψm, ρ)
	ρϕpψm = expect(ψmϕp, ρ) #expect(ϕpψm, ρ)

	y = real(ρϕpϕp - ρψmψm)
	x = 2 * real(ρϕpψm)

	return 0.5 * (arctan(y, x) + π)
end

# ╔═╡ 84dc0aef-6230-423f-bb6d-6a5bf66652fc
md"""
Then the Hamiltonian is 

$H_S(t) = \Gamma_S c(t) (\sigma_{y1} - \sigma_{y2})/2$

where $\Gamma_S$ is a constant rate and $c(t)$ is its modulation based on the angle. To achieve the appropriate rotation we choose $\Delta \theta(t) dt = \Gamma_S c(t)/2$ so that 

$H_S(t) = \frac{\Delta \theta(t)}{dt} (\sigma_{y1} - \sigma_{y2}).$
"""

# ╔═╡ a3554988-fa10-425e-b535-690fdf668615
md"""
### Feedback (rotating) Hamiltonian
"""

# ╔═╡ f8a1debc-8dff-4638-969d-926c756dc2d0
md"""
The feedback Hamiltonian is given by

$H_F(t) = \sqrt{\Gamma_H} P(t) (\sigma_{y,1} + \sigma_{y,2}) / \sqrt2$

where $\Gamma_H = 0.5$ MHz is the measurement rate corresponding to continuous measurement of the observable $X_H = (\sigma_{z,1} + \sigma_{z,2} )/ \sqrt{2}$. To get the appropriate rotation angle of $\theta_\text{opt}$ we thus want $\sqrt{\Gamma_H / 2} P(t) = \theta_\text{opt}/dt$, i.e.

$H_F(t) = \frac{\theta(t)}{dt} (\sigma_{y,1} + \sigma_{y,2})$
"""

# ╔═╡ 3cca9018-c461-45df-8c5c-dc7cef7dfa69
md"""
### Total Hamiltonian
"""

# ╔═╡ 56426fae-e189-4ddc-bcc2-2a598c08f844
md"""
The total known Hamiltonian (labeled $H_\text{fil}$ for "filter") is thus

$H_\text{fil}(t) = H_F(t) + H_S(t),$
"""

# ╔═╡ c88ff48f-ddd0-4b30-8149-01f9303144ab
md"""
while the total true Hamiltonian (labeled $H_\text{sys}$ for "system") is 

$H_\text{sys}(t) = H_\text{fil}(t) + H_z(t)$
"""

# ╔═╡ a3e273ff-6ce8-430d-a4f9-7432fab9cce2
md"""
## Measurement
"""

# ╔═╡ 5ea7a206-94e7-4122-8cf1-cb9e0b48e0fb
md"""
We measure  $X_H = (\sigma_{z,1} + \sigma_{z,2} )/ \sqrt{2}$ at rate $\Gamma_H$ with unit quantum efficiency $\eta = 1.0$:
"""

# ╔═╡ c1aa9c52-2251-4d35-8c31-a7509e8f1da8
C = [((σz1 + σz2)/√2, ΓH, 1.0)]

# ╔═╡ 49ba376b-bd44-4544-b641-d81b0f1a6a46
md"""
# Simulation (QuantumCircuits.jl)
"""

# ╔═╡ a26a8f10-3eb4-4e1c-98b8-e15652df6ce3
md"""
## Single trajectory (with fluctuators)
"""

# ╔═╡ 1730ce6e-af0c-4125-89c7-e2bb5f66e7c2
# plot(sol.t, map(ρ -> real(expect(dm(Ψm), ρ)), sol.ρ))

# ╔═╡ 425eb1c0-5e87-41bf-971a-ad77f0085de7
md"""
## Single trajectory (system-filter)
"""

# ╔═╡ 87f49f41-7d2c-4410-b1ed-027758032dc7
dtt = 50*dt

# ╔═╡ b06c4db4-29fa-4b5f-aa1c-09568020b0c8
begin
	thresholded = true
	threshold = (2π) * 50.0
end

# ╔═╡ e88781df-61ae-4553-ae7f-a2549f3d9b8b
magΔθ(ρ) = 	if thresholded
				inst = - Δθ(ρ) / (2dt) 
				(abs(inst) < threshold) ? inst : sign(inst) * threshold
			else
				- Δθ(ρ) / (2dtt)
			end

# ╔═╡ 64913d45-adca-436b-b0e8-153ec5a1a6f5
HS(ρ::State) = (magΔθ(ρ)/2) * (σy1 - σy2)

# ╔═╡ df90c78f-52fa-48c6-a487-1091acd3129f
magθ(ρ) = 	if thresholded
				inst = - θ(ρ) / (2dt) 
				(abs(inst) < threshold) ? inst : sign(inst) * threshold
			else
				- θ(ρ) / (2dtt)
			end

# ╔═╡ 8232efab-e666-4459-a8b1-797b075d6d9d
HF(t::Timescale, ρ::State) = (magθ(ρ)/2) * (σy1 + σy2)

# ╔═╡ f6a02c48-62b3-4d69-9ff1-185088aefe4f
Hfil(t::Timescale, ρ::State) = HF(t, ρ) + HS(ρ)

# ╔═╡ db474913-ae34-48ed-b100-f4feb5087ec3
Hsys(Hf1, Hf2) = (t::Timescale, ρ::State) -> Hfil(t, ρ) + Hz(Hf1, Hf2)(t)

# ╔═╡ c62ffa10-3eae-43b2-9925-f44e3378e560
begin
	Hi1 = Hf1(τi, tf) # new realization of fluctuator Hamiltonian (single qubit)
	Hi2 = Hf2(τi, tf) # new realization of fluctuator Hamiltonian (single qubit)
	
	(sys, fil) = bayesian((0.0, tf), ρ0, (Hsys(Hi1, Hi2), Hfil), [], C; dt=dt)
end

# ╔═╡ c0a8e3e6-c136-4287-9967-73fd9b7ec88e
sol = bayesian((0.0, tf), ρ0, Hsys(Hi1, Hi2), [], C; dt=dt)

# ╔═╡ cb094b31-ae08-427d-8ad5-4f333a81c405
bell_plot(sol)

# ╔═╡ 73a9d4dd-1447-4763-9d1d-a190be90b4b6
begin
	plot(sol.t, Δθ.(sol.ρ) /(2π * dtt), xlabel="t (μs)", label=L"|H_S| (counter)")
	plot!(sol.t, θ.(sol.ρ) / (2π * dtt), label=L"|H_F| (rotating)",legendfontsize=10)
end

# ╔═╡ 983e1c7e-63a6-4a23-8cd4-2076ab1bec9c
begin
	ops = bell_basis
	op_labels = bell_basis_strings
	
	system_exps = map(op -> expectations(sys, op), ops)
	filter_exps = map(op -> expectations(fil, op), ops)	
end

# ╔═╡ ca36b0a2-3718-459d-a637-2065b1a68a2c
bell_plot(sys)

# ╔═╡ b41a5c71-a16b-4f86-adda-e1807681edde
bell_plot(sys, fil)

# ╔═╡ 52cbc301-9ca3-43c2-9b74-c89e3691ec10
begin
	plot(sol.t, magθ.(sol.ρ) / (2π), xlabel="t (μs)", label=string(L"|H_F|/2", " (rotating)"))
	plot!(sol.t, magΔθ.(sol.ρ) / (2π), label=string(L"|H_S|/2", " (counterrotating)"),legendfontsize=10, ylabel="Hamiltonian in MHz")
end

# ╔═╡ 0f5189fb-f05c-4aa5-98ae-218b13acb4b4
md"""
# Utilities
"""

# ╔═╡ f8eb195f-e0d9-48ef-bf22-e1af09cd47a1
xyz(θ, ϕ) = (sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ))

# ╔═╡ 9ac6937c-245e-4ee5-a92b-84ae6432cc82
function purities(sol::Solution)
	ρs = (typeof(sol.ρ[1]) <: Ket) ? dm.(sol.ρ) : sol.ρ
	map(ρ -> real(tr(ρ * ρ)), ρs)
end

# ╔═╡ 0a6f05f6-6bf3-4ca9-8543-3633b6819e87
purities(exps) = 0.5 .* (1 .+ exps[1].^2 .+ exps[2].^2 .+ exps[3].^2)

# ╔═╡ 1fc7a91b-8303-4a24-8cd9-3430033cf316
function ensembleavg(sols::Vector{Solution}; ops = qbasis)
	exparrs = map(op -> [], ops)

	for sol in sols
		exps = map(op -> expectations(sol, op), ops)
		for (list, traj) in zip(exparrs, exps)
			push!(list, traj)
		end
	end

	return map(arr -> mean(arr), exparrs)

end

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

# ╔═╡ 06ace002-be3d-4478-860d-7dc23838c9de
function bloch_vectors(vecs...; mesh=30, ax=false, viewϕ=0, connecting=false, labels=[], size=(400,400))
	bcolors = palette(:seaborn_bright)

	# Wire frame coordinates ---------------------------------------------------

	x(θ, ϕ) = sin(θ) * cos(ϕ + viewϕ)
	y(θ, ϕ) = sin(θ) * sin(ϕ + viewϕ)
	z(θ, ϕ) = cos(θ)

	θs = range(0, 2π, length=mesh)
	ϕs = range(0, π, length=mesh) .+ viewϕ


	# Plot wireframe -----------------------------------------------------------
	wf = plot()
	
	# Longitudes
	for ϕ in ϕs
		plot!([x(θ, ϕ) for θ in θs], [y(θ, ϕ) for θ in θs], [z(θ, ϕ) for θ in θs], label=:none, linecolor="steelblue", linewidth=0.5, linealpha=1, seriestype=path3d, aspect_ratio=1.0, size=size)
			
	end

	# Latitudes
	for θ in θs
		plot!([x(θ, ϕ) for ϕ in ϕs], [y(θ, ϕ) for ϕ in ϕs], [z(θ, ϕ) for ϕ in ϕs], label=:none, linecolor="steelblue", linewidth=0.5, linealpha=1, seriestype=path3d, aspect_ratio=1.0, size=size)
	end


	# Plot reference axes ------------------------------------------------------
	
	colors = [palette(:tab10)[i] for i in 1:3]

	if ax
		plot!([0, cos(viewϕ)], [0, sin(viewϕ)], [0, 0], linecolor=colors[1], linewidth=3.0, label=:none)
		plot!([0, -sin(viewϕ)], [0, cos(viewϕ)], [0, 0], linecolor=colors[2], linewidth=3.0, label=:none)
		plot!([0, 0], [0, 0], [0, 1], linecolor=colors[3], linewidth=3.0, label=:none)
	end


	# Plot Bloch vectors input by user --------------------------------

	for (i, vec) in enumerate(vecs)

		(xvv, yvv, zvv) = vec

		xv = xvv * cos(viewϕ) - yvv * sin(viewϕ)
		yv = xvv .* sin(viewϕ) .+ yvv * cos(viewϕ)
		zv = zvv

		if connecting
			plot!([0, xv], [0, yv], [0, zv], label=:none, linewidth=2, linecolor=bcolors[i])
		end
	
		plot!([xv], [yv], [zv], legend=:outerright, legendfontsize=10, marker=(:circle, 5), markercolor=bcolors[mod(i, length(bcolors)) + 1], label=try labels[i] catch e "" end)
	end

	return wf


end

# ╔═╡ bb58f61d-334f-424e-80c1-4730abb59050
mod(10, 10)

# ╔═╡ Cell order:
# ╟─e625f413-e93b-48c9-b14d-181c37b0fe0a
# ╟─249f33ea-92aa-4b9c-bdfa-d35fe5254225
# ╟─3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─be587838-24f3-4b4b-8a51-789168089aa7
# ╟─bb91d819-871b-48e3-92fa-7a494a67c836
# ╠═d172ea75-59de-475f-9c0f-365cc1f09b85
# ╟─cf4feaf2-0fa6-449c-86d1-ac988b5df569
# ╠═3cddac9b-5ab2-4120-9355-00411c76dce3
# ╠═41dc7c32-dc78-41ea-b48a-3a99415da603
# ╟─5c1c9b7c-a16d-453c-966c-55d31ef64066
# ╟─dabf9a58-907f-43ac-a45f-9fa902536f6a
# ╠═a7d7312f-3dab-4f10-8509-2b8f90d8074f
# ╟─c81d1356-120b-49dd-947e-e66ee6952a21
# ╟─54a89faa-4c16-46f4-8852-aa363b3fa101
# ╠═e07021db-e4d4-400a-9b83-61e0ffd19159
# ╠═5ae2203c-4a4a-4265-932e-be0deb01c97e
# ╠═8aab4b6c-2ce4-4667-a122-1fa4d8ae37e7
# ╟─a33e0d94-d99b-4643-8e6a-d4136237d48c
# ╠═5bc80650-4fb3-41ed-bc4b-15f51f937ed5
# ╟─84dc0aef-6230-423f-bb6d-6a5bf66652fc
# ╠═64913d45-adca-436b-b0e8-153ec5a1a6f5
# ╠═e88781df-61ae-4553-ae7f-a2549f3d9b8b
# ╟─a3554988-fa10-425e-b535-690fdf668615
# ╟─f8a1debc-8dff-4638-969d-926c756dc2d0
# ╠═8232efab-e666-4459-a8b1-797b075d6d9d
# ╠═df90c78f-52fa-48c6-a487-1091acd3129f
# ╟─3cca9018-c461-45df-8c5c-dc7cef7dfa69
# ╟─56426fae-e189-4ddc-bcc2-2a598c08f844
# ╠═f6a02c48-62b3-4d69-9ff1-185088aefe4f
# ╟─c88ff48f-ddd0-4b30-8149-01f9303144ab
# ╠═db474913-ae34-48ed-b100-f4feb5087ec3
# ╟─a3e273ff-6ce8-430d-a4f9-7432fab9cce2
# ╟─5ea7a206-94e7-4122-8cf1-cb9e0b48e0fb
# ╠═c1aa9c52-2251-4d35-8c31-a7509e8f1da8
# ╟─49ba376b-bd44-4544-b641-d81b0f1a6a46
# ╠═a26a8f10-3eb4-4e1c-98b8-e15652df6ce3
# ╠═c0a8e3e6-c136-4287-9967-73fd9b7ec88e
# ╠═cb094b31-ae08-427d-8ad5-4f333a81c405
# ╠═1730ce6e-af0c-4125-89c7-e2bb5f66e7c2
# ╠═73a9d4dd-1447-4763-9d1d-a190be90b4b6
# ╠═425eb1c0-5e87-41bf-971a-ad77f0085de7
# ╠═c62ffa10-3eae-43b2-9925-f44e3378e560
# ╠═983e1c7e-63a6-4a23-8cd4-2076ab1bec9c
# ╠═ca36b0a2-3718-459d-a637-2065b1a68a2c
# ╠═87f49f41-7d2c-4410-b1ed-027758032dc7
# ╠═b06c4db4-29fa-4b5f-aa1c-09568020b0c8
# ╠═b41a5c71-a16b-4f86-adda-e1807681edde
# ╠═52cbc301-9ca3-43c2-9b74-c89e3691ec10
# ╟─0f5189fb-f05c-4aa5-98ae-218b13acb4b4
# ╠═f8eb195f-e0d9-48ef-bf22-e1af09cd47a1
# ╠═9ac6937c-245e-4ee5-a92b-84ae6432cc82
# ╠═0a6f05f6-6bf3-4ca9-8543-3633b6819e87
# ╠═1fc7a91b-8303-4a24-8cd9-3430033cf316
# ╠═1b7e3334-9825-4832-a0b8-d58816fd5236
# ╠═2e072a88-7166-4d5d-96e6-52fcbe668505
# ╠═06ace002-be3d-4478-860d-7dc23838c9de
# ╠═bb58f61d-334f-424e-80c1-4730abb59050
# ╠═df539f73-7c8a-4d65-b909-0d94720f724f
