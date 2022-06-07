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

# ╔═╡ a5beb57b-9a1a-421e-80c1-1a9fca19b6e0
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
	using LaTeXStrings
	using Random
	using Statistics
	using Distributions
	using Plots
	using DataFrames
	using Dates
	using QuantumCircuits

	include("utilities/single-qubit-operators.jl")
	include("utilities/utilities.jl")
	include("utilities/plotting.jl")

	include("notebooks/table-of-contents.jl")
	include("resources.jl")

	# where to save exported gifs (optional)
	datapath = "/Users/sachagreenfield/Desktop/Physics/Research/2021-Linear-feedback/data/";

	md" # Packages and julia files"
	
	
end

# ╔═╡ 377a3336-20bd-4baa-a033-af8bbc8668a8
md"""
In this interactive notebook, we'll study some basic principles of measurement backaction in a qubit.
"""

# ╔═╡ 696c4757-cebb-4d96-a837-0e9e5155e5a7
mdp(table_of_contents📔)

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents(title="Measurement backaction principles")

# ╔═╡ 8f7c6440-eac8-49d3-bd23-c2545ec16830
md"""
# Measurement backaction
"""

# ╔═╡ 727ced31-7a51-4892-b1d9-eff4fcdc9767
md"""
Weak measurement of a qubit observable $\hat A$ generates a stochastic record

$\tilde r_t= \braket{\hat A}_t + \zeta_t$

where $\braket{\hat A}_t  \equiv \text{Tr} \big(\rho_q(t) \hat A \big)$ is the expectation value of the observable at time $t$, and 

$\zeta_t \sim \mathcal{N} \Big(0, \sqrt{\tau_m/ dt} \Big)$

is a zero-mean, Gaussian distributed random variable with variance $\sigma = \tau_m / dt$. Here, $\tau_m$ is the measurement timescale determined by the strength of the measurement.
"""

# ╔═╡ 1aaa9cbe-2f98-41e1-af41-598ba6578818
md"""
In this notebook, we will consider measurement of the qubit state $\ket0$ or $\ket1$, such that $\hat A = e^{i \varphi} \hat \sigma_z$. The phase $\varphi$ is called the "measurement angle" and correpsonds to the quadrature of amplification of the signal. It determines the type of backaction experienced by the qubit.
"""

# ╔═╡ 9c821306-43ac-4577-9687-d8b0a4e4cecb
md"""
## Informational vs. phasal backaction
"""

# ╔═╡ 875adb1b-22f7-4375-83c5-2a17c9a8e46c
md"""

Rotate view: `ϕv = 0`
$(@bind ϕvIB html"<input type=range min=0 max=32 step=1 value=0>") 
`ϕv = 2π`

"""

# ╔═╡ 3531558a-f9b7-4f6b-9800-884fe0b04712
md" $(@bind animIB_bs CheckBox()) Animate Bloch sphere "

# ╔═╡ 71604112-0704-4832-8a80-e886b6155e34
md"""
`φ =` $(@bind φIB Select(["0", "π/8", "π/4", "3π/8", "π/2"]))
"""

# ╔═╡ d3e48ac4-973f-4397-bf06-89150a0c13ed
md"""
Informational ("spooky") backaction occurs when $\varphi = 0$. Phase ("realistic") backaction occurs when $\varphi = \pi/2$. Intermediate angles result in a combination of both types of backaction.

In the following simulation, the state is initialized in the state $\ket{+x}$. By default, the measurement angle is $\varphi = 0$, leading to informational backaction. The measurement angle is currently $\varphi$ = $φIB.

Interact with the plots below and try different parameter settings to gain intuition for informational and phase backaction.
"""

# ╔═╡ faf0b11c-4339-4a82-a23c-9d35eb9d10b4
md" $(@bind animIB CheckBox(default=false)) Animate informational backaction"

# ╔═╡ 67e9c47a-d688-4724-8b21-90472311951b
begin
	exportpath = let
					ep = string(datapath, "animations/", Dates.now())
					replace(ep, ":" => "-") end
	md" $(@bind exportgifs CheckBox(default=false)) Export gifs"
end

# ╔═╡ 8f8606a1-bf7e-469d-bf63-7d93cee714c4
string(exportpath, "-projections.gif")

# ╔═╡ 3ff551a9-3a07-4a2f-928e-880b7e3ba1fc
if 	φIB == "0"
	
	md""" **Informational backaction**  The qubit undergoes a biased random walk towards one of the poles. For $\eta = 1$, the walk is confined to the circle $x^2 + z^2 = 1$. Once it reaches a pole, it is pinned there. As an aside, this demonstrates how continuous weak monitoring can be understand as a projective measurement in the limit of measuring for a long time. The statistics of final states will follow that of projective measurement, in the absence of Hamiltonian evolution.

While the initial state creates a bias towards the pole it is closer to, it is possible to stabilize to either pole so long as the initial state is a superposition.
	"""

elseif φIB == "π/2"
	
	md""" **Phase backaction**  The qubit undergoes a unbiased random walk in the $x$ -- $y$ plane.  For $\eta = 1$, the walk is confined to the circle $x^2 + y^2 = 1$. Note that this effect is preserved for a state initialized outside of the $x$ -- $y$ plane; the random walk will then be confined to its original plane: $x^2 + y^2 = 1 - z^2$.

	"""

else
	
	md""" **Informational and phase backaction**  The qubit undergoes a biased random walk towards one of the poles.  For $\eta = 1$, the walk is confined to the sphere $x^2 + y^2 + z^2 = 1$. 

	"""

end

# ╔═╡ 911d4aa1-41cd-4eea-8e0f-1af6f8024290
@bind go Button("Rerun simulation")

# ╔═╡ 3ae402b0-476a-4c11-823c-d0016793adab
md"""
Change the measurement angle to see how backaction depends on the measurement quadrature.
"""

# ╔═╡ 8d7bed45-4431-41df-9226-f689d1ca1977
md"""
## Measurement collapse in the ensemble
"""

# ╔═╡ a70b0d26-fa32-4b94-9e06-73bd10742c5e
md" $(@bind animate_collapse CheckBox()) Animate collapse "

# ╔═╡ ab188d7c-9eb4-48f6-86e6-53c057ce67ae
md"""
## Linear feedback
"""

# ╔═╡ 5b83ae62-9486-485f-821b-eb127f01e487
md"""
Linear feedback takes advantage of the fact that informational (phase) backaction looks like random rotations about the $x$ - or $y$ -axes ($z$ -axis). Such rotations can be mirrored by Hamiltonian evolution. Thus, feeding back the measurement record into the Hamiltonian can be used to "erase" certain effects of measurement backaction.
"""

# ╔═╡ 71addbb9-4362-470d-80d1-63f3eed4b7c3
md"""
Informational linear feedback, implemented in $Patti_et_al, erases the effects of measurement backaction while changing the effective measurement pole. 
"""

# ╔═╡ 133b6939-10b2-4c8e-acf8-5658ca96a0f9
md" # Utilities"

# ╔═╡ a7347fdf-723b-4ba9-9431-7c2ab61a1106
md"""
## Plotting
"""

# ╔═╡ 5f00c276-33ff-4d37-a51a-831b496c7e1c
function qubit_plot(sol::Solution; record=false, title="", legendpos=:bottomleft)

	basis = qbasis

	t = sol.t
	exps = map(op -> expectations(sol, op), basis)
	r = record ? sol.r[1] : []

	return qubit_plot((t, exps, r); title=title, legendpos=legendpos)
end

# ╔═╡ d46fe705-2d2f-4ecb-93f5-39d5de028417
# Plotting
function plot_timeseries(series...; plot_title="time series", xlabel="t (μs)", ylabel="arbitrary units", labels=[], colors = palette(:lightrainbow), kwargs...)

	label(i) = i > length(labels) ? i : labels[i]

	p = plot(size=(600,300))

	# Plot records vs. time
	for (i, (tt, ser)) in enumerate(series)
		plot!(tt, ser, color=colors[i], label=label(i), xlabel=xlabel, ylabel=ylabel, title=plot_title, legend=:outerright, kwargs...)
	end
	
    p
	
end

# ╔═╡ b8295b27-47d3-4d93-b2fc-4f0a8319bf02
# Plotting
function plot_timeseries(tt::Vector{Timescale}, series...; plot_title="time series", xlabel=L"$t$", ylabel="arbitrary units", labels=[], colors = palette(:lightrainbow), kwargs...)

	ser_colors(i) = colorpairs ? colors[i] : colors[2i]
	label(i) = i > length(labels) ? i : labels[i]

	p = plot(size=(600,300))

	# Plot records vs. time
	for (i, ser) in enumerate(series)
		plot!(tt, ser, color=ser_colors(i), label=label(i), xlabel=xlabel, ylabel=ylabel, title=plot_title, legend=:outerright, kwargs)
	end
	
    p
	
end

# ╔═╡ 7d8d7147-10b5-42c9-a68b-fc98cc0ddb7e
# Plotting
function series_histogram(tt, series, index; title="time series", xlabel="t (μs)", ylabel="arbitrary units", labels=[], color=:blue)

	t = tt[index]
	ser = series[1]
	p = plot([t, t], [-1, 1], linestyle=:dash, color=:black, size=(600,300), linewidth=1.5)

	# Plot records vs. time
	histdata = []
	for (i, ser) in enumerate(series)
		plot!(tt, ser, color=color, xlabel=xlabel, ylabel=ylabel, title=title, legend=:none, alpha=0.5, linewidth=0.7)
		push!(histdata, ser[index])
	end

	h =	histogram(histdata, bins=(-1):0.1:1, normalize=false, color=color, alpha=0.9, orientation=:h, legend=:none, xlabel="occurences", xlims=[0, round(length(series)/2)])
		
	l = @layout [timeseries{0.7w} histogram{0.3w}]
	return plot(p, h, layout=l, link=:y)
	
end

# ╔═╡ f837190b-d6e9-4f20-a641-f10f8e2c678c
# Plotting
function plot_records(series; plot_title="time series", xlabel="t (μs)", ylabel="arbitrary units", labels=[], colors = palette(:lightrainbow), histograms=false)

	label(i) = i > length(labels) ? i : labels[i]

	p = plot(size=(600,300))

	# Plot records vs. time
	for (i, (tt, ser)) in enumerate(series)
		plot!(tt, ser, color=colors[i], xlabel=xlabel, ylabel=ylabel, title=plot_title, legend=:none)
	end

	if !histograms
	    return p
		
	else
		
		h =	histogram(title="record histograms", xlabel="probability", orientation=:h)
		for (i, (tt, ser)) in enumerate(series)
			(μ, σ) = params(fit(Normal, ser))
			histogram!(ser, bins = :scott, normalize=true, color=colors[i], alpha=0.65, label=string(labels[i], ", σ = ", round(σ, digits=3)), orientation=:h, legend=:bottomright)
		end

		l = @layout [timeseries{0.5w} histogram{0.5w}]
		return plot(p, h, layout=l, link=:y)
		
	end
	
end

# ╔═╡ 8a6cceb1-31fa-44b4-a95d-5c095a57555c
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

# ╔═╡ 60513c9e-719e-4bbe-ac93-2d6a557d4982
md"""
#### Colors
"""

# ╔═╡ 15424721-26c1-4f40-898a-41edebc886e2
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

# ╔═╡ d0e342fc-d3b3-4ef8-9ebe-da8bd0621643
function qubit_plot((t, exps, r); title="", legendpos=:bottomleft)

	record = (r != [])

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

	plot!(t, p, color=colors[4], label=L"Tr(\rho^2)")

	if !record
		return pl
	else
		l = @layout [blochs{0.6h}; record{0.4h}]
		pr = plot(t, r, color=colors1[1], xlabel="t (μs)", ylabel="record", label=:none, legend=legendpos, title="", linewidth=0.8)
		return plot(pl, pr, layout = l, link=:y)
	end
end

# ╔═╡ 226b3d09-d81e-4fcc-ac75-ec39f7e0b037
function qubit_plot(sol1::Solution, sol2::Solution; record=false, title="", color1=colorsq1, color2=colorsq2, l1="", l2="")

	basis = qbasis

	t1 = sol1.t
	exps1 = map(op -> expectations(sol1, op), basis)
	r1 = record ? sol1.r[1] : []
	
	t2 = sol2.t
	exps2 = map(op -> expectations(sol2, op), basis)
	r2 = record ? sol2.r[1] : []
	
	return qubit_plot((t1, exps1, r1), (t2, exps2, r2); title=title, color1=color1, color2=color2, l1=l1, l2=l2)
	
end

# ╔═╡ 7c618a96-6770-4853-abae-185e426bbc15
function qubit_plot((t1, exps1, r1), (t2, exps2, r2); title="", color1=colorsq1, color2=colorsq2, l1="", l2="")

	basis = qbasis
	labels = qlabels

	record = ((r1 != []) || (r2 != []))

	p1 = 0.5 .* (1 .+ exps1[1].^2 .+ exps1[2].^2 .+ exps1[3].^2)
	p2 = 0.5 .* (1 .+ exps2[1].^2 .+ exps2[2].^2 .+ exps2[3].^2)

	pl = plot(size=(600,300), legendfontsize=10, titlefontsize=12, legend=:outerright, ylabel="bloch coordinates", xlabel = record ? "" : "t (μs)", linewidth=1.2, title=title)

	for l in 1:length(basis)
		plot!(t1, exps1[l], color=color1[l], label=labels[l], legend=:outerright, ylims=[-1,1], linewidth=1.2)
	end
	plot!(t1, p1, color=color1[4], label=L"Tr(\rho^2)", linewidth=1.2)
	
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
			plot!(t1, r1, color=mixed_colors[1], xlabel="t (μs)", ylabel="record", label=:none, legend=:outerright, title="")
		end
		if r2 != []
			plot!(t2, r2, color=mixed_colors[2], xlabel="t (μs)", ylabel="record", label=:none, legend=:outerright, title="")
		end
		return plot(pl, pr, layout = l, link=:both)
	end
end

# ╔═╡ 99fec4d3-d0e2-4a53-b460-d5f46d34e97e
md"""
## Other
"""

# ╔═╡ 2150a7fe-2403-4b2f-9e03-c9f374557457
begin
	mutable struct traj
	  t::Vector{Float64}
	  x::Vector{Float64}
	  y::Vector{Float64}
	  z::Vector{Float64}
	  p::Vector{Float64}
	  r
	end
	
	function traj(sol::Solution)
	  t, ρ, r = (sol.t, sol.ρ, sol.r)
	  x, y, z = map(op -> expectations(sol, op), [σx, σy, σz])
	  p = (typeof(ρ[1]) <: Ket) ? [1.0 for el in ρ] : real(expect.(ρ, ρ))
	  traj(t, x, y, z, p, r)
	end
end

# ╔═╡ ccbcf668-d948-4ec6-a5f7-39a178d54c29
xyz(θ, ϕ) = (sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ))

# ╔═╡ 618168dc-53dd-4562-8061-67a0b56587aa
φdict = Dict("0" => 0, 
				"π/8" => π/8, 
				"π/4" => π/4,
				"3π/8" => 3π/8,
				"π/2" => π/2)

# ╔═╡ 86314c80-9a75-4268-8f39-8fe621079580
begin
	# System parameters ------------------------------------------------------------
	# all times given in μs
	
	# initial state 
	ψ0 = normalize(g + e)
	
	# measurement parameters
	Γ = 0.5				# rate
	η = 1.0 				# efficiency
	φ = get(φdict, φIB, 0) 	# angle
	
	# simulation timescales
	tf = 8.0
	T = (0.0, tf) # simulation duration
	dt = 1e-3  # integration time-step
	ΩR = 0.0
end

# ╔═╡ 49065a29-93ac-458d-ac0c-84f2f2151f9d
let
	# Kraus operators -------------------------------------------------------------
	H = (ΩR / 2) * σy
	J = [(σz, ((1 - η) * Γ/2))]
	C = [(exp(im * φ) * σz, Γ, η)]
	
	# Bayesian simulation ---------------------------------------------------------
	go
	global solB = bayesian(T, ψ0, H, J, C; dt=dt)
	
	global IB = traj(solB)
	
	md" ###### 🔻 Simulation "
	
end

# ╔═╡ f9eaaf70-4e0f-4503-8a64-2380682354ce
let 
	ϕv = ϕvIB * (π/16)
	tt, xx, yy, zz = (IB.t, IB.x, IB.y, IB.z)
	
	if animIB_bs
		
		anim = @animate for i ∈ range(1, length(tt), step=100)
			blochsphere(xx[1:i], yy[1:i], zz[1:i], linewidth=1., linealpha=0.85, ax=true, viewϕ = ϕv, blochmark=true) end

		if exportgifs
			gif(anim, string(exportpath, "-blochsphere.gif"), fps = 15)
		else
			gif(anim, fps = 15)
		end
		
	else
		blochsphere(xx, yy, zz, linewidth=1., linealpha=0.85, ax=true, viewϕ = ϕv, blochmark=true)
		
	end
	
end

# ╔═╡ 6bc3e01f-b3cb-4a32-ab6b-e5dcc967b07f
let 
	tt, xx, yy, zz = (IB.t, IB.x, IB.y, IB.z)
	
	if animIB
		anim = @animate for i ∈ range(1, length(tt), step=100)
			blochprojections(xx[1:i], yy[1:i], zz[1:i], blochmark=true) end
		if exportgifs
			gif(anim, string(exportpath, "-projections.gif"), fps = 15)
		else
			gif(anim, fps = 15)
		end
	else
		blochprojections(xx, yy, zz, blochmark=true)
	end
end

# ╔═╡ 34ac25b0-daff-4e71-9ef7-43c13692a4bd
qubit_plot(solB, record=true, legendpos=:bottomright)

# ╔═╡ b9258bba-17cf-4b93-ac0d-6b2d517cdadb
md"""
Looking at an ensemble of trajectories, we see that the statistics of continuous measurement collapse resemble those of projective measurement: initializing in $\ket{+}$, we achieve 50/50 statistics for measurement $\ket{\pm z}$. Most collapses occur within a few measurement times, $\tau_m = 1/2\Gamma$ = $(1/(2Γ)) μs.
"""

# ╔═╡ 5af97094-fcdd-4f5d-a176-5b3e7a1fc2a8
begin
	N1 = 200 # number of realizations
	
	# Kraus operators -----------------------------------------
	H = (ΩR / 2) * σy
	J(η) = (η == 1.0) ? [] : [(σz, (1 - η) * Γ/2)]
	C(η) = (η == 0.0) ? [] : [(exp(im * φ) * σz, Γ, η)]


	sols = map(1:N1) do m
		bayesian((0, 1.5tf), ψ0, H, J(η), C(η); dt=dt)
	end

	η0_sol = bayesian((0, 1.5tf), ψ0, H, J(0.0), C(0.0); dt=dt)

	collapse_zs = map(sol -> expectations(sol, σz), sols[1:50])
end

# ╔═╡ ab361db2-4ba5-47cf-83fe-cd8075a79de1
bloch_plots(sols, η0_sol, alpha=0.15, N=N1)

# ╔═╡ 1c6f257d-0810-4d2f-9e15-09f01bb213e5
if animate_collapse
	
	anim = @animate for i ∈ range(1, length(sols[1].t), step=300)
	series_histogram(sols[1].t, collapse_zs, i, color=colorsq1[3], ylabel="z", title="measurement collapse") end
	gif(anim, fps = 5)
	
else
	series_histogram(sols[1].t, collapse_zs, length(sols[1].t), color=colorsq1[3], ylabel="z", title="measurement collapse")
end

# ╔═╡ 8abc2b3f-e3ee-484f-a73b-344e7f4c3988
md"""
## Paths
"""

# ╔═╡ 9c6e0b0b-513d-4bad-954c-36302418c562
refpath = string(path, "notebooks/5-linear-feedback/2-linear-feedback.jl")

# ╔═╡ cc115e97-479c-4123-92eb-1b9dae3f8a80
Markdown.parse("""
You can study these effects in more detail in [Linear feedback📔](./open?path=$refpath)."""
)

# ╔═╡ Cell order:
# ╟─377a3336-20bd-4baa-a033-af8bbc8668a8
# ╟─696c4757-cebb-4d96-a837-0e9e5155e5a7
# ╟─3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─8f7c6440-eac8-49d3-bd23-c2545ec16830
# ╟─727ced31-7a51-4892-b1d9-eff4fcdc9767
# ╟─1aaa9cbe-2f98-41e1-af41-598ba6578818
# ╟─9c821306-43ac-4577-9687-d8b0a4e4cecb
# ╟─d3e48ac4-973f-4397-bf06-89150a0c13ed
# ╠═86314c80-9a75-4268-8f39-8fe621079580
# ╠═49065a29-93ac-458d-ac0c-84f2f2151f9d
# ╟─875adb1b-22f7-4375-83c5-2a17c9a8e46c
# ╟─3531558a-f9b7-4f6b-9800-884fe0b04712
# ╟─71604112-0704-4832-8a80-e886b6155e34
# ╟─f9eaaf70-4e0f-4503-8a64-2380682354ce
# ╟─faf0b11c-4339-4a82-a23c-9d35eb9d10b4
# ╠═67e9c47a-d688-4724-8b21-90472311951b
# ╟─6bc3e01f-b3cb-4a32-ab6b-e5dcc967b07f
# ╟─8f8606a1-bf7e-469d-bf63-7d93cee714c4
# ╠═34ac25b0-daff-4e71-9ef7-43c13692a4bd
# ╟─3ff551a9-3a07-4a2f-928e-880b7e3ba1fc
# ╟─911d4aa1-41cd-4eea-8e0f-1af6f8024290
# ╟─3ae402b0-476a-4c11-823c-d0016793adab
# ╟─8d7bed45-4431-41df-9226-f689d1ca1977
# ╟─b9258bba-17cf-4b93-ac0d-6b2d517cdadb
# ╠═5af97094-fcdd-4f5d-a176-5b3e7a1fc2a8
# ╠═ab361db2-4ba5-47cf-83fe-cd8075a79de1
# ╟─a70b0d26-fa32-4b94-9e06-73bd10742c5e
# ╠═1c6f257d-0810-4d2f-9e15-09f01bb213e5
# ╟─ab188d7c-9eb4-48f6-86e6-53c057ce67ae
# ╟─5b83ae62-9486-485f-821b-eb127f01e487
# ╟─71addbb9-4362-470d-80d1-63f3eed4b7c3
# ╟─cc115e97-479c-4123-92eb-1b9dae3f8a80
# ╟─133b6939-10b2-4c8e-acf8-5658ca96a0f9
# ╟─a7347fdf-723b-4ba9-9431-7c2ab61a1106
# ╟─5f00c276-33ff-4d37-a51a-831b496c7e1c
# ╟─d0e342fc-d3b3-4ef8-9ebe-da8bd0621643
# ╟─226b3d09-d81e-4fcc-ac75-ec39f7e0b037
# ╟─7c618a96-6770-4853-abae-185e426bbc15
# ╟─d46fe705-2d2f-4ecb-93f5-39d5de028417
# ╟─b8295b27-47d3-4d93-b2fc-4f0a8319bf02
# ╟─7d8d7147-10b5-42c9-a68b-fc98cc0ddb7e
# ╟─f837190b-d6e9-4f20-a641-f10f8e2c678c
# ╟─8a6cceb1-31fa-44b4-a95d-5c095a57555c
# ╟─60513c9e-719e-4bbe-ac93-2d6a557d4982
# ╠═15424721-26c1-4f40-898a-41edebc886e2
# ╟─99fec4d3-d0e2-4a53-b460-d5f46d34e97e
# ╠═2150a7fe-2403-4b2f-9e03-c9f374557457
# ╠═ccbcf668-d948-4ec6-a5f7-39a178d54c29
# ╠═618168dc-53dd-4562-8061-67a0b56587aa
# ╟─8abc2b3f-e3ee-484f-a73b-344e7f4c3988
# ╟─9c6e0b0b-513d-4bad-954c-36302418c562
# ╠═a5beb57b-9a1a-421e-80c1-1a9fca19b6e0
