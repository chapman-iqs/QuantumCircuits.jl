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

# ╔═╡ 9b379ec4-cdfe-4110-ac35-39a6653cfa2b
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
	using LaTeXStrings
	using Random
	using Statistics
	using Distributions
	using QuantumCircuits
	using Plots

	include("utilities/qubit-resonator-operators.jl")
	include("utilities/plotting.jl")

	include("notebooks/table-of-contents.jl")
	include("resources.jl")
	
	md" # Packages and julia files"
end

# ╔═╡ 377a3336-20bd-4baa-a033-af8bbc8668a8
mdp("In this interactive notebook, we explore qubit-resonator simulations where the resonator is simulated as a harmonic oscillator. This gives intuition for the subtleties of readout that may not be captured by the coherent state approximation described in ", coherent_states📔, ".")

# ╔═╡ 3b27a77c-4efc-46f7-a1fd-f22522ed2d46
md"""
**Note that this notebook has an unresolved bug where the dual-quadrature measurement leads to unusually high photon numbers...**
"""

# ╔═╡ 00ae069c-3d04-4937-9e12-51637a237aaf
mdp(table_of_contents📔)

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents(title="Qubit-resonator dynamics")

# ╔═╡ ede93004-ef9e-461f-b4ea-062f9d7879f3
md"""
# Qubit-resonator description
"""

# ╔═╡ f1ed08db-97a9-4439-8f06-619ff9bffcc1
md" $(@bind nonideal CheckBox()) Include T1, T2 effects "

# ╔═╡ fa6e2db2-0f58-41cd-9193-5ad05e501e81
md"""
Measurement type: $(@bind quad Select(["single quadrature", "dual quadrature"]))
"""

# ╔═╡ dcd54469-41a8-4505-b0ed-15d0545f8742
begin
	single_quadrature = (quad == "single quadrature")
	title = string(quad, " measurement")
	ideal = !nonideal
end

# ╔═╡ 49492b7b-10b4-4a93-bd3e-29835ceec60f
md" $(@bind anim_plot1 CheckBox()) Animate cavity centers "

# ╔═╡ bab3e2b6-2796-4167-8393-e3f04568abf8
md"""
# Shifting the resonator frequency
"""

# ╔═╡ afc18faf-9a07-4974-8b9e-bffc44d5429e
1/1.56

# ╔═╡ 8af5af9f-a6cd-41a8-8c34-db981bee7096


# ╔═╡ 7e62e6c6-aa3c-4350-901c-15d017b8db42
md" # Utilities "

# ╔═╡ b4ca38d2-e56b-4977-893d-3b3e856e238d
md"""
## Plotting
"""

# ╔═╡ 1aec6816-7205-4651-a733-a5b7f01b1bcf
function cavity_plot(ts, εt, (αp_list, αm_list); xlims=[-0.5,0.5], ylims=[-0.5,0.5], steadystate=false)
	return cavity_plot(length(ts), ts, εt, (αp_list, αm_list); xlims=xlims, ylims=ylims, steadystate=steadystate)
end

# ╔═╡ 159bb702-fb0e-4aaf-91f7-86c2e2b8d9a3
function cavity_plot(ts, εt, (αp_list, αm_list), (αp_list2, αm_list2); xlims=[-0.5,0.5], ylims=[-0.5,0.5])
	return cavity_plot(length(ts), ts, εt, (αp_list, αm_list), (αp_list2, αm_list2); xlims=xlims, ylims=ylims)
end

# ╔═╡ 2a1eb3ba-64ff-4bb2-b9dd-fdc2a9f09acb
function qubit_plot(sol::Solution; record=false, title="", legendpos=:bottomleft, rlims=:default, rlabels=["r1", "r2"])

	basis = qbasis

	t = sol.t
	exps = map(op -> expectations(sol, op), basis)
	rs = record ? sol.r : []

	return qubit_plot((t, exps, rs); title=title, legendpos=legendpos, rlims=rlims, rlabels=rlabels)
end

# ╔═╡ c5b251a5-7a70-466f-a9d6-91d9704ce039
function bloch_plots(sols::Vector{Solution}, sol_η0::Solution; alpha=0.1, N=50)
	colors = colors1q 
	
	# calculate expectation values --------------------------------------------
	t = sols[1].t
	xs, ys, zs = [], [], []
	
	for sol in sols
		x, y, z = map(op -> expectations(sol, op), qbasis)
		for (list, traj) in zip([xs, ys, zs], [x, y, z])
			push!(list, traj)
		end
		
	end

	# η = 0 solution
	xη0, yη0, zη0 = map(op -> expectations(sol_η0, op), qbasis)

	
	
	# plot ----------------------------------------------------------------------
	function bloch(os, oη0; color=colors1q[1], xlabel="", ylabel="")
		
		po = plot(ylims = [-1,1], xlabel=xlabel, ylabel=ylabel)
		
		for o in os[1:min(N, 50)]
			plot!(t, o, alpha=alpha, label=:none, color=color)
		end
		
		oavg = [mean([os[i][j] for i in 1:N]) for j in 1:length(t)]
		plot!(t, oavg, alpha=1, color=color, label="average", linewidth=3)
		plot!(t, oη0, alpha=1, color=:black, label="η = 0", linewidth=2)

		po
		
	end

	
	l = @layout [xplot{0.33h}; yplot{0.33h}; zplot{0.33h}]
	
	px = bloch(xs, xη0, color=colors[1], ylabel="x")
	py = bloch(ys, yη0, color=colors[2], ylabel="y")
	pz = bloch(zs, zη0, color=colors[3], ylabel="z", xlabel="t (μs)")
	
	plot(px, py, pz, layout = l, link=:y, size=(800,500), legendfontsize=8, titlefontsize=12, legend=:outerright)
	
end

# ╔═╡ be064ca9-ac02-4b74-ba51-ce7ec361939b
function animate_plot(tt, plot, args...; step=100, fps=15)
	anim = @animate for i ∈ range(1, length(tt), step=step)
			plot(i, args...) end
	gif(anim, fps = fps)
end

# ╔═╡ 46605307-01f2-4f35-b0d0-fc4229e24038
begin
	colors1 = palette(:rainbow)
	colors2 = palette(:lightrainbow)
	colorsq1 = palette(:tab20)[1:2:20]
	colorsq2 = palette(:tab20)[2:2:20]

	mixed_colors = let
		cols = []
		for i in 1:length(colors1)
			push!(cols, colors1[i])
			push!(cols, colors2[i])
		end
		cols
	end
end

# ╔═╡ 5934b164-fab7-4a9c-ab1f-a06296fb8b58
function cavity_plot(i::Int64, ts, εt, (αp_list, αm_list); xlims=[-0.5,0.5], ylims=[-0.5,0.5], steadystate=false)

	colors=colors1
	
	αps = αp_list[1:i]
	αms = αm_list[1:i]

	labels = steadystate ? 
				[L"\alpha_+^{s.s.}", L"\alpha_-^{s.s.}"] :
				[L"\alpha_+", L"\alpha_-"]

	ls = steadystate ? :dash : :solid

	pα = plot(xlims=xlims, ylims=ylims, xlabel=string("Re ", L"\alpha_\pm"), ylabel=string("Im ", L"\alpha_\pm"), legend=:right, title="coherent state wavepacket centers", titlefontsize=12)
	
	plot!(real.(αps), imag.(αps), color=colors[2], linestyle=ls, label=:none)
	plot!([real(last(αps))], [imag(last(αps))], color=colors[2], marker="o", label=labels[1])
	plot!(real.(αms), imag.(αms), color=colors[4], linestyle=ls, label=:none)
	plot!([real(last(αms))], [imag(last(αms))], color=colors[4], marker="o", label=labels[2])

	plot!(ts, εt, color=colors[5], label="", xlabel = "t (μs)", ylabel = "ε (MHz)", legend=:none, inset = (1, bbox(0.05, 0.05, 0.4, 0.25, :top, :right)), subplot=2, title="cavity drive", titlefontsize=10, axisfontsize=10)
	plot!([ts[i]], [εt[i]], color=colors[5], marker="o", label=:none, subplot=2)
	
end

# ╔═╡ 224bae70-b220-4f8e-9532-fddd8ca486e2
function cavity_plot(i::Int64, ts, εt, (αpss, αmss), (αp_list2, αm_list2); xlims=[-1.5, 1.5], ylims=[-1.5,1.5])

	colors=colors1
	
	αps1 = αpss[1:i]
	αms1 = αmss[1:i]
	αps2 = αp_list2[1:i]
	αms2 = αm_list2[1:i]

	pα = plot(xlims=xlims, ylims=ylims, xlabel=string("Re ", L"\alpha_\pm"), ylabel=string("Im ", L"\alpha_\pm"), legend=:outerbottomright, title="coherent state wavepacket centers", titlefontsize=12)
	
	plot!(real.(αps1), imag.(αps1), color=colors[2], linestyle=:dash, label=:none)
	plot!([real(last(αps1))], [imag(last(αps1))], color=colors[2], marker=:cross, label=string(L"\alpha_+", " (s.s.)"))
	plot!(real.(αms1), imag.(αms1), color=colors[4], linestyle=:dash, label=:none)
	plot!([real(last(αms1))], [imag(last(αms1))], color=colors[4], marker=:cross, label=string(L"\alpha_-", " (s.s.)"))

	plot!(real.(αps2), imag.(αps2), color=colors[2], label=:none)
	plot!([real(last(αps2))], [imag(last(αps2))], color=colors[2], marker="o", label=L"\alpha_+")
	plot!(real.(αms2), imag.(αms2), color=colors[4], label=:none)
	plot!([real(last(αms2))], [imag(last(αms2))], color=colors[4], marker="o", label=string(L"\alpha_-"))

	plot!(ts, εt, color=colors[5], label=:none, xlabel = "t (μs)", ylabel = "ε (MHz)", legend=:none, inset = (1, bbox(0.05, 0.05, 0.4, 0.25, :top, :right)), subplot=2, title="cavity drive", titlefontsize=10, axisfontsize=10)
	plot!([ts[i]], [εt[i]], color=colors[5], marker="o", label=:none, subplot=2)
	
end

# ╔═╡ b2a4380b-b49c-4600-8fe8-7341faf6a649
function qubit_plot((t, exps, rs); title="", legendpos=:bottomleft, rlims=:default, rlabels=["r1", "r2"])

	record = (rs != [])

	basis = qbasis
	colors = colorsq1
	labels = qlabels

	p = 0.5 .* (1 .+ exps[1].^2 .+ exps[2].^2 .+ exps[3].^2)

	pl = plot(size=(600,300), legendfontsize=10, titlefontsize=12, legend=:outerright, ylabel="bloch coordinates", xlabel = record ? "" : "t (μs)", linewidth=1.5, title=title)

	for l in 1:length(basis)
		label = labels[l]
		color = colors[l]
		exp = exps[l]
		plot!(t, exps[l], color=color, label=label, legend=legendpos, ylims=[-1,1], linewidth=1.5)
	end

	plot!(t, p, color=colors[4], label=L"Tr(\rho_q^2)")

	if !record
		return pl
	else
		l = @layout [blochs{0.6h}; record{0.4h}]
		pr = plot(xlabel="t (μs)", ylabel="record", label=:none, legend=legendpos, title="", linewidth=0.8, ylims=rlims)
		for (i, r) in enumerate(rs)
			plot!(t, r, color=mixed_colors[i], label=rlabels[i])
		end
		return plot(pl, pr, layout = l, link=:y)
	end
end

# ╔═╡ 0184ea80-f1b8-41f5-9d88-c70f4bf7fbd0
function qubit_plot(sol1::Solution, sol2::Solution; record=false, title="", color1=colorsq1, color2=colorsq2, l1="", l2="", rlims=:default)

	basis = qbasis

	t1 = sol1.t
	exps1 = map(op -> expectations(sol1, op), basis)
	r1 = record ? sol1.r[1] : []
	
	t2 = sol2.t
	exps2 = map(op -> expectations(sol2, op), basis)
	r2 = record ? sol2.r[1] : []
	
	return qubit_plot((t1, exps1, r1), (t2, exps2, r2); title=title, color1=color1, color2=color2, l1=l1, l2=l2, rlims=rlims)
	
end

# ╔═╡ a16bbccc-6cba-4832-99ce-680e50296959
function qubit_plot((t1, exps1, r1), (t2, exps2, r2); title="", color1=colorsq1, color2=colorsq2, l1="", l2="", rlims=:default)

	basis = qbasis
	labels = qlabels

	record = ((r1 != []) || (r2 != []))

	p1 = 0.5 .* (1 .+ exps1[1].^2 .+ exps1[2].^2 .+ exps1[3].^2)
	p2 = 0.5 .* (1 .+ exps2[1].^2 .+ exps2[2].^2 .+ exps2[3].^2)

	pl = plot(size=(600,300), legendfontsize=10, titlefontsize=12, legend=:outerright, ylabel="bloch coordinates", xlabel = record ? "" : "t (μs)", linewidth=1.2, title=title)

	for l in 1:length(basis)
		plot!(t1, exps1[l], color=color1[l], label=labels[l], legend=:outerright, ylims=[-1,1], linewidth=1.2)
	end
	plot!(t1, p1, color=color1[4], label=L"Tr(\rho_q^2)", linewidth=1.2)
	
	for l in 1:length(basis)
		plot!(t2, exps2[l], color=color2[l], label=:none, legend=:outerright, ylims=[-1,1], linestyle=:dash, linewidth=2)
	end
	plot!(t2, p2, color=color2[4], label=:none, linestyle=:dash, title=string(l1, " ---, ", l2, " - - -"), linewidth=2)


	if !record
		return pl
	else
		l = @layout [blochs{0.6h}; record{0.4h}]
		pr = plot()
		if r1 != []	
			plot!(t1, r1, color=mixed_colors[1], xlabel="t (μs)", ylabel="record", label=:none, legend=:outerright, title="", ylims=rlims)
		end
		if r2 != []
			plot!(t2, r2, color=mixed_colors[2], xlabel="t (μs)", ylabel="record", label=:none, legend=:outerright, title="", ylims=rlims)
		end
		return plot(pl, pr, layout = l, link=:both)
	end
end

# ╔═╡ 278694c3-6035-4603-8498-6c5126442f8a
md"""
## Other
"""

# ╔═╡ 69602c2b-84f6-45cd-a7c4-6a13ccaf2d33
function findmax(times, series; inrange=())
	if inrange == ()
		range = 1:length(times)
	else
		i1 = findfirst(x -> x == first(inrange), times)
		i2 = findfirst(x -> x == last(inrange), times)
		range = i1:i2
	end
	t = times[range]
	ser = series[range]
	i = findfirst(x -> x == maximum(ser), ser)
	return (t[i], ser[i])
end

# ╔═╡ 3d1472b2-30ac-4237-ace3-00ad5f5768b6
xyz(θ, ϕ) = (sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ))

# ╔═╡ 76d5307b-f58b-4ab1-8ba5-b43491c54552
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
	
	function traj(sol::Solution)
		t, ρ, r = (sol.t, sol.ρ, sol.r)
		σs = (size(ρ[1]) == (2, 2)) ? (σxq, σyq, σzq) : (σx, σy, σz)
		x, y, z =  [real(expect(σi, ρ)) for σi in σs] 
		p = real(expect.(ρ, ρ))
		traj(t, x, y, z, p, r)
	end
	
	
end

# ╔═╡ a823cf8d-8637-44b9-af2a-8feefd22d986
begin
	
	mutable struct qr_traj
		t::Vector{Float64}
		
		x::Vector{Float64}
		y::Vector{Float64}
		z::Vector{Float64}
		
		n::Vector{Float64}
		ap::Vector{ComplexF64}
		am::Vector{ComplexF64}
		
		p::Vector{Float64}
		pq::Vector{Float64}
		r
	end
	
	
	function qr_traj(sol::Solution)
		t, ρ, r = (sol.t, sol.ρ, sol.r)
		
		# get expectation values
		x, y, z =  [real(expect(σi, ρ)) for σi in (σx, σy, σz)] 
		n, zpp, zmm = [real(expect(op, ρ)) for op in ((a' * a), zp, zm)] 
		αpp, αmm = [expect(op, ρ) for op in (αp, αm)] 
		p = typeof(ρ[1]) <: Ket ? [1.0 for tt in t] : real(expect.(ρ, ρ))

		# calculate functions of exp. values
		pq = 0.5 .* (1 .+ x.^2 .+ y.^2 .+ z.^2)
		ap = αpp ./ zpp
		am = αmm ./ zmm
		
		qr_traj(t, x, y, z, n, ap, am, p, pq, r)
	end
	
end

# ╔═╡ 55b946e9-20e2-4ea8-85e2-455dd9b6b8a1
let
	# Parameters -------------------------------------------------------------------
	
	dt = 1e-3 # integration time step (μs)
	tf = 10.0 # total time duration (μs)
	
	# carrier frequencies - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	# 	ωr = 2π * (6.67e3) # rad MHz; bare resonator frequency
	# 	ωε = 2π * (6.67e3)# rad MHz; resonator readout pulse carrier frequency
	# 	ωq = 2π * (5.56e3) # rad MHz; bare qubit frequency
	#   ωR = 2π * (5.56e3) # 2π*(5.559871e3) # rad MHz; qubit drive carrier frequency
		
	
	# detunings - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	Δc = 0.0 	# ωr - ωε # cavity readout pulse from bare cavity frequency (MHz)
	Δq = 0.0 # 0.129 	# ωq - ωR: Rabi drive from bare qubit frequency (MHz)
	
	
	# envelopes - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	ΩR = 0.0 #2π * (0.4) 	# rad MHz; Rabi drive
	ε = 2π * (0.7) 		# rad MHz; cavity readout drive envelope
	χ = 2π * (-0.47) 	# rad MHz; dispersive shift / 2

	global εt = [ε for t in 0.0:dt:tf]


	# measurement parameters - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	θR = π/2 	# Rabi-drive axis angle in x-y plane, w.r.t. x-axis
	φ = 0 		# measurement quadrature
	η = 1.0 		# signal collection efficiency
	
	
	# rates - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	T1 = 160 		# qubit energy decay time (μs)
	γ1 = 1/T1 		# qubit energy decay rate
	
	Tϕ = 16 		# cavity-induced dephasing time (μs)
	γϕ = 1/Tϕ 		# cavity-induced dephasing rate
	
	κ = 2π * (1.56) # cavity linewidth / decay rate (rad Hz)
	τ = 1/(2η * κ) 	# measurement collapse time
	
	# initial state  - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	θ0 = π/2 	# initial qubit polar angle
	ϕ0 = 0 		# initial qubit azimuthal angle
	
	ψq = normalize(cos(θ0/2) * g + exp(im * ϕ0) * sin(θ0/2) * e) # initial qubit state
	ψc = fockstate(f, 0) 	# initial cavity state (vacuum)
	ψ0 = ψq ⊗ ψc 			# initial qubit-resonator state
	
	
	# Kraus operators -----------------------------------------------------------
	
	# Hamiltonian - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	Hc = Δc * (a' * a) 								# cavity energy
	Hq =  Δq * σz / 2 								# qubit energy
	Hqc = χ * (a' * a) * σz 						# cavity-qubit dispersive coupling
	Hε = (ε / 2) * (a + a') 						# cavity readout tone
	HR = (ΩR / 2) * (cos(θR) * σx + sin(θR) * σy) 	# qubit Rabi drive
	
	H = Hc + Hq + Hqc + Hε + HR
	# H = Hqc + Hε
	
	
	# Lindblad operators - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	J =	[] #ideal ? [(a, (1 - η) * κ/2)] : [(a, (1 - η) * κ/2), (σm, γ1), (σz, γϕ/2)]
	
	# Measurement operators - - - - - - - - - - - - - - - - - - - - - - - - - - 
	C = single_quadrature ? 
			[(exp(im * φ) * a, κ, η)] : 
			[(exp(im * φ) * a, κ, η/2), (exp(im * (φ + π/2)) * a, κ, η/2)]

	
	# Simulation ----------------------------------------------------------------
	
	global sol1 = bayesian((0.0, tf), ψ0, H, J, C; dt=dt)
	global traj1 = qr_traj(sol1)
	
end

# ╔═╡ b658c1fa-ac7c-4cf3-8fb7-13abff46af26
qubit_plot(sol1, record=true, legendpos=:bottomright, rlabels=["I", "Q"], title=title)

# ╔═╡ 36fcc1a1-b8e0-4ac3-b79e-eb3d288bee0c
begin
	p1 = plot(sol1.t, traj1.n, ylabel="photon number", label=:none, title=title, titlefontsize=12)
	p2 = plot(sol1.t, traj1.p, color=:red, label="qubit-resonator", xlabel="t (μs)", ylims=[0.5,1], ylabel="purities")
	plot!(sol1.t, traj1.pq, color=:red, linestyle=:dash, label="reduced qubit", legend=:right)
	plot(p1, p2, layout=grid(2,1))
end

# ╔═╡ 433b3ddd-dfe9-44c9-a673-bbcf72cadbff
if anim_plot1
		animate_plot(sol1.t, cavity_plot, sol1.t, (εt ./ 2π), (traj1.ap, traj1.am); step=150, fps=10)
	else
		cavity_plot(sol1.t, (εt ./ 2π), (traj1.ap, traj1.am))
	end

# ╔═╡ 3aca6d8d-1c65-4585-9578-9ddbc8dab5cb
let
	# Parameters -------------------------------------------------------------------
	
	dt = 1e-3 # integration time step (μs)
	tf = 5.0 # total time duration (μs)

	# detunings - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	Δc = 0.0 	# ωr - ωε # cavity readout pulse from bare cavity frequency (MHz)	

	# rates - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	φ = 0 		# measurement quadrature
	η = 1.0 		# signal collection efficiency
	
	κ = 2π * (1.56) # cavity linewidth / decay rate (rad Hz)
	τ = 1/(2η * κ) 	# measurement collapse time
	
	# envelopes - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	ε = 2π * (0.7) 		# rad MHz; cavity readout drive envelope

	ωε = κ/2
	global εfunc(t) = ε * (cos(ωε * t))^2


	# measurement parameters - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	

	
	# initial state  - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	
	ψc = fockstate(f, 0) 	# initial cavity state (vacuum)
	ψ0 = g ⊗ ψc 			# initial qubit-resonator state
	
	
	# Kraus operators -----------------------------------------------------------
	
	# Hamiltonian - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	Hc = Δc * (a' * a) 								# cavity energy
	Hε(t) = (εfunc(t) / 2) * (a + a') 						# cavity readout tone
	
	H(t) = Hc + Hε(t)
	
	# Lindblad operators - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	J =	[] #ideal ? [(a, (1 - η) * κ/2)] : [(a, (1 - η) * κ/2), (σm, γ1), (σz, γϕ/2)]
	
	# Measurement operators - - - - - - - - - - - - - - - - - - - - - - - - - - 
	C = [(exp(im * φ) * a, κ, η/2), (exp(im * (φ + π/2)) * a, κ, η/2)]

	
	# Simulation ----------------------------------------------------------------
	
	global sol2 = bayesian((0.0, tf), ψ0, H, J, C; dt=dt)
	global traj2 = qr_traj(sol2)
	
end

# ╔═╡ 64bcd706-8096-4cb6-83e3-9602638a43ba
npeak = findmax(sol2.t, traj2.n; inrange=(3.2, 3.6))

# ╔═╡ c7bd4723-2f27-4233-8fad-6505adf9a55a
εpeak = findmax(sol2.t, εfunc.(sol2.t); inrange=(3.0, 3.4))

# ╔═╡ ebf20139-da0f-47de-834e-9987b0c53e9e
δt = npeak[1] - εpeak[1]

# ╔═╡ 678a2b09-e84e-4e33-928e-239b4f335049
δt * 1.56

# ╔═╡ 4bf1f0c9-3d5f-4976-b298-18bdba0c77ff
@bind index Slider(1:length(sol2.t))

# ╔═╡ 16801edb-ca17-49c5-926b-e61c9780908f
@bind index2 Slider(1:length(sol2.t))

# ╔═╡ 31ab44c5-586d-41ba-88ec-955a88ec5009
let
	p1 = plot(sol2.t, traj2.n, ylabel="photon number", label="photon number",title="dual quadrature measurement", titlefontsize=12)
	plot!(sol2.t, εfunc.(sol2.t)/15, label="ε(t)", legend=:bottomright, color=:red, linestyle=:dash, xlabel="t (μs)")
	plot!([sol2.t[index], sol2.t[index]], [0, 0.35], color=:black, linestyle=:solid, label=string("t = ", round(sol2.t[index], digits=3), " μs"))
	plot!([sol2.t[index2], sol2.t[index2]], [0, 0.35], color=:green, linestyle=:solid, label=string("t = ", round(sol2.t[index2], digits=3), " μs"))
end

# ╔═╡ b136e6e0-6e88-4008-94bb-7756716e913e
φdict = Dict("0" => 0, 
				"π/8" => π/8, 
				"π/4" => π/4,
				"3π/8" => 3π/8,
				"π/2" => π/2)

# ╔═╡ dfb8774e-bb14-4f8b-9492-ecd2220e393d
md"""
## Old utilities
"""

# ╔═╡ a255f7a5-d617-4444-80a4-4e4e7762b236
green(text; title="Note") = Markdown.MD(Markdown.Admonition("correct", title, [text]))

# ╔═╡ 6e86b5ee-5b38-4245-b555-c208151c0c53
red(text; title="Note") = Markdown.MD(Markdown.Admonition("danger", title, [text]))

# ╔═╡ fb32f241-34bb-413b-99d2-c0155b363670
tan(text; title="Note") = Markdown.MD(Markdown.Admonition("warning", title, [text]))

# ╔═╡ e730e11b-6fa5-4121-80eb-69a59d74b852
blue(text; title="Note") = Markdown.MD(Markdown.Admonition("note", title, [text]))

# ╔═╡ c8aa6f9a-a2b7-45d7-9e39-344380a3909b
hint(text; title="Hint") = Markdown.MD(Markdown.Admonition("hint", title, [text]))

# ╔═╡ Cell order:
# ╟─377a3336-20bd-4baa-a033-af8bbc8668a8
# ╟─3b27a77c-4efc-46f7-a1fd-f22522ed2d46
# ╠═00ae069c-3d04-4937-9e12-51637a237aaf
# ╟─3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─ede93004-ef9e-461f-b4ea-062f9d7879f3
# ╟─f1ed08db-97a9-4439-8f06-619ff9bffcc1
# ╟─fa6e2db2-0f58-41cd-9193-5ad05e501e81
# ╟─dcd54469-41a8-4505-b0ed-15d0545f8742
# ╠═55b946e9-20e2-4ea8-85e2-455dd9b6b8a1
# ╠═b658c1fa-ac7c-4cf3-8fb7-13abff46af26
# ╠═36fcc1a1-b8e0-4ac3-b79e-eb3d288bee0c
# ╟─49492b7b-10b4-4a93-bd3e-29835ceec60f
# ╠═433b3ddd-dfe9-44c9-a673-bbcf72cadbff
# ╠═bab3e2b6-2796-4167-8393-e3f04568abf8
# ╠═3aca6d8d-1c65-4585-9578-9ddbc8dab5cb
# ╠═31ab44c5-586d-41ba-88ec-955a88ec5009
# ╠═64bcd706-8096-4cb6-83e3-9602638a43ba
# ╠═c7bd4723-2f27-4233-8fad-6505adf9a55a
# ╠═ebf20139-da0f-47de-834e-9987b0c53e9e
# ╠═afc18faf-9a07-4974-8b9e-bffc44d5429e
# ╠═678a2b09-e84e-4e33-928e-239b4f335049
# ╠═4bf1f0c9-3d5f-4976-b298-18bdba0c77ff
# ╠═16801edb-ca17-49c5-926b-e61c9780908f
# ╠═8af5af9f-a6cd-41a8-8c34-db981bee7096
# ╟─7e62e6c6-aa3c-4350-901c-15d017b8db42
# ╟─b4ca38d2-e56b-4977-893d-3b3e856e238d
# ╟─5934b164-fab7-4a9c-ab1f-a06296fb8b58
# ╟─224bae70-b220-4f8e-9532-fddd8ca486e2
# ╟─1aec6816-7205-4651-a733-a5b7f01b1bcf
# ╟─159bb702-fb0e-4aaf-91f7-86c2e2b8d9a3
# ╟─2a1eb3ba-64ff-4bb2-b9dd-fdc2a9f09acb
# ╟─b2a4380b-b49c-4600-8fe8-7341faf6a649
# ╟─0184ea80-f1b8-41f5-9d88-c70f4bf7fbd0
# ╟─a16bbccc-6cba-4832-99ce-680e50296959
# ╟─c5b251a5-7a70-466f-a9d6-91d9704ce039
# ╟─be064ca9-ac02-4b74-ba51-ce7ec361939b
# ╟─46605307-01f2-4f35-b0d0-fc4229e24038
# ╟─278694c3-6035-4603-8498-6c5126442f8a
# ╟─69602c2b-84f6-45cd-a7c4-6a13ccaf2d33
# ╠═3d1472b2-30ac-4237-ace3-00ad5f5768b6
# ╟─76d5307b-f58b-4ab1-8ba5-b43491c54552
# ╠═a823cf8d-8637-44b9-af2a-8feefd22d986
# ╠═b136e6e0-6e88-4008-94bb-7756716e913e
# ╟─dfb8774e-bb14-4f8b-9492-ecd2220e393d
# ╟─a255f7a5-d617-4444-80a4-4e4e7762b236
# ╟─6e86b5ee-5b38-4245-b555-c208151c0c53
# ╟─fb32f241-34bb-413b-99d2-c0155b363670
# ╟─e730e11b-6fa5-4121-80eb-69a59d74b852
# ╟─c8aa6f9a-a2b7-45d7-9e39-344380a3909b
# ╠═9b379ec4-cdfe-4110-ac35-39a6653cfa2b
