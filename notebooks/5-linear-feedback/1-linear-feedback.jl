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

# ╔═╡ 77d623fb-7dbc-406c-bf46-05385a049530
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
	datapath = "/Users/sachagreenfield/Desktop/Physics/Research/2021-Linear-feedback/data/"

	md" # Packages and julia files"
	
	
end

# ╔═╡ 377a3336-20bd-4baa-a033-af8bbc8668a8
md"""
In this interactive notebook, we'll look at dynamical linear feedback stabilization of quantum trajectories, based on $Patti_et_al.
"""

# ╔═╡ e77d8c4b-73fc-4ba1-8054-f5da97a36748
mdp("This notebook relies on the concepts of informational and phasal backaction, introduced in ",  measurement_backaction📔, ".")

# ╔═╡ 9014ea2e-c244-4aeb-bc56-811dcfb6b70c
mdp(table_of_contents📔)

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents(title="Linear feedback stabilization")

# ╔═╡ ab188d7c-9eb4-48f6-86e6-53c057ce67ae
md"""
# Linear feedback
"""

# ╔═╡ 5b83ae62-9486-485f-821b-eb127f01e487
md"""
Linear feedback takes advantage of the fact that informational (phase) backaction looks like random rotations about the $x$ - or $y$ -axes ($z$ -axis). Such rotations can be mirrored by Hamiltonian evolution. Thus, feeding back the measurement record into the Hamiltonian can be used to "erase" certain effects of measurement backaction.

Informational linear feedback, introduced by $(Patti_et_al), erases the effects of measurement backaction while changing the effective measurement pole. First, we will gain intuition for how this works in the next section. Then, we will attempt to apply similar principles to stabilize phase backaction.
"""

# ╔═╡ 1eb3c79e-2efb-4a4e-9fb2-07d78e9cf21a
md"""
# Informational linear feedback
"""

# ╔═╡ d724fb1e-bd93-4a69-8082-eadaa50b31a0
md"""
## Description of the problem
"""

# ╔═╡ 34e1740c-20d6-4099-ab7a-0c64160957c4
md"""
A transmon qubit is dispersively coupled to a cavity. The cavity bare resonant frequency undergoes a dispersive shift $\pm \chi$ depending on whether the qubit is in the ground or excited state. Thus, measuring the frequency of the cavity yields qubit information, which is encoded in a measurement record $\tilde r$.

Linear feedback stabilization feeds the measurement record back into qubit Hamiltonian evolution to dynamically stabilize a particular state $\ket{\psi_\text{target}}$ on the Bloch sphere. In the ideal scenario, the state is pinned to $\ket{\psi_\text{target}}$ with fidelity $\rightarrow 1$ after a transient evolution. 

As non-idealities are introduced (such as $\eta < 1$, or $T_1, T_2 < \infty$, or finite feedback delay $T_d > 0$), the fidelity of stabilization is reduced.
"""

# ╔═╡ 7ff6b7dc-ba8b-4ce3-8a1c-0308e871986d
md" ##### Weak monitoring"

# ╔═╡ 7fb819fc-f0ae-49ac-a702-8b5329adccd3
md"""
The qubit is weakly monitored due to the small cavity population ($\langle a^\dagger a \rangle \approx 0.2$ photons). Continuous monitoring of the cavity yields a time-dependent stochastic record of the form 

$\tilde r (t) = z(t) + \zeta(t)$

where $z(t) \equiv \text{Tr} \big( \rho_q(t) \sigma_z \big)$ is the $z$ Bloch coordinate of the qubit, and 

$\zeta(t) \sim \mathcal{N} \Big(0, \sqrt{\tau_m/ dt} \Big)$

is a zero-mean, Gaussian distributed random variable with variance $\sigma = \tau_m / dt$.
"""

# ╔═╡ 18c9725c-ebad-49b8-841c-57893558251f
md" ##### Measurement backaction"

# ╔═╡ 1edc57b3-de6a-4546-abab-78b780922081
md"""
The weak measurement of the qubit state creates informational ($\sigma_z$) measurement backaction. This means that the qubit state evolves according to readout-dependent measurement Kraus operators

$M_{\tilde r} \propto e^{- \tilde{r} dt / 2\tau_m} \ket{0} \bra{0} + e^{\tilde{r} dt / 2\tau_m} \ket{1} \bra{1} = \exp(\tilde r dt \sigma_z / 2\tau_m).$

In other words, measurement of $\tilde r$ leads to an outcome-dependent state update. Critically, in the Bloch picture, note that $M_{\tilde r}$ will only modify the $z$ coordinate (see [1], Eq. 7). This has two important consequences:

* In the absence of a $\sigma_z$ (qubit energy) term in the Hamiltonian, dynamics will be constrained to the plane of Rabi oscillations.

* The effect of $M_\tilde{r}$ can be counteracted by an appropriate feedback into Rabi oscillations, which rotate $z$ into $x$ and/or $y$.


"""

# ╔═╡ e338a516-c4c8-4571-8e84-d85dcc68c985
md"""
Note that this scheme relies on the backaction being **informational** (towards $\pm z$). This is achieved via phase-sensitive amplification of the readout pulse along the *informational* axis. This discriminates between cavity states $a \ket{0} \bra{0}$ and $a \ket{1}\bra{1}$ by phase. However, if the readout pulse is amplified in quadrature, this effectively measures the cavity photon number, creating backaction in the qubit Stark shift and thus creating random rotations about $\sigma_z$. This is called **phase** backaction.

Later in this notebook, we will look into whether a similar scheme might be used to counteract phase backaction. This would allow for stabilization of qubits measured via phase-preserving amplification, which results in both informational and phase backaction.

"""

# ╔═╡ a7a38426-33b0-4057-b125-0488d4bcb279
md" ##### Feedback Hamiltonian"

# ╔═╡ a87a6bce-c5b6-4b04-9d14-d9656969d87b
md"""

The Hamiltonian is a Rabi drive modulated by the feedback from the measurement readout:

$\hat H = \hbar \big[\Delta_0 + \Delta_1 \tilde{r}(t - T_d) \big] \frac{\hat \sigma_\phi}2$

where

$\hat \sigma_\phi = \cos \phi \hat \sigma_x + \sin \phi \hat \sigma_y,$

and $\hat \sigma_x$ and $\hat \sigma_y$ are the corresponding Pauli operators. Thus there is a base Rabi frequency $\Delta_0$ modulated by a linear feedback term $\Delta_1 \tilde{r} (t-T_d)$. The feedback is delayed by $T_d$ due to finite speed of transmission of the signal through transmission line.

The parameters $\Delta_0$, $\Delta_1$ are chosen to balance the predicted (informational) measurement backaction with a Rabi drive in the opposite direction. 

Following the analysis in [1], this notebook focuses on clockwise rotations in the $y$ -- $z$ plane, which corresponds to $\phi = \pi$.

"""

# ╔═╡ fb578156-094b-43b1-8c64-740f64b193bc
md" ## Simulation "

# ╔═╡ 5e76a75a-3041-49f6-b5b6-96047b1d5bd4
ideal = true

# ╔═╡ b5ddea46-9aee-4ce9-8230-902611228ca7
begin
	x0 = 0
	y0 = 0.3
	z0 = sqrt(1 - y0^2)
	ρ0 = DenseOperator(0.5*(Iq + x0*σx + y0*σy + z0*σz))
end

# ╔═╡ e486fbed-157c-4139-914c-ada6bb88d7b4
md" ##### Options: "

# ╔═╡ 4898ab97-4058-4e70-a959-2962641d9611
md" $(@bind nonideal CheckBox(default=false)) Include T1, T2 effects "

# ╔═╡ 6b93bb84-0ba1-44e8-910e-612793b3df3b
md"""
Change signal collection efficiency `η =` $(@bind ηstr Select(["1.0", "0.9", "0.8", "0.7", "0.6", "0.5", "0.4", "0.3", "0.2", "0.1", "0"]))
"""

# ╔═╡ ce79057b-eb31-49b1-9db7-bfd793a447c5
md"""
Add time delay `Td =` $(@bind Tdstr Select(["0", "1", "2", "4", "8", "20", "40", "80", "200", "400"])) ns
"""

# ╔═╡ ca7a2351-cff5-4d77-ba16-00de384b8b7c
md"""
##### Plots and animations
"""

# ╔═╡ f3e7794e-a9a1-4103-8469-32a8e37d2d82
md"""

Rotate view: `ϕv = 0`
$(@bind ϕvc html"<input type=range min=0 max=32 step=1 value=0>") 
`ϕv = 2π`

"""

# ╔═╡ a12cdb8c-e9a1-4c2d-9811-cff266e152d8
md" $(@bind show_gif CheckBox()) Animate Bloch sphere "

# ╔═╡ 0a7f28c9-1d84-43e6-b62c-711a231a3972
md" $(@bind show_gif2 CheckBox()) Animate Bloch series "

# ╔═╡ d9f2f00f-4ee2-45b5-91d6-6552d6d5b6c1
md" $(@bind show_gif3 CheckBox()) Animate cross-sections "

# ╔═╡ e28cc7a9-ccc7-418a-94fa-068befe4f151
md"""
##### Discussion
"""

# ╔═╡ 3679d9db-97aa-40cd-ab1b-2329051c7156
md"""
###### Optimizing the fidelity
"""

# ╔═╡ 1159cb26-2a33-46af-a2cb-3cc667b79c0d
md""" $T_1$ and $T_2$ effects are mild enough to not be visible on the simulation timescale for small $\tau_m$. However, the stabilization fidelity rapidly decreases with $\eta<1$. There is no straightforward way to counteract this effect.

Even for $\eta = 1$, a realistic measurement delay time ($T_d \approx 200$ ns) rapidly reduces stabilization fidelity. This can be alleviated by increasing the ensemble-average measurement collapse timescale $\tau_m$. In the bad-cavity limit, this is approximately

$\tau_m = 1/2 \Gamma \approx \frac{\kappa}{16 \chi^2 \bar n}.$

Thus, the finite feedback delay can be alleviated by taking a very weak measurement, e.g. by making the average cavity population $\bar n$ very small.
"""

# ╔═╡ 08f808ab-13f7-4799-951c-72042d71e1be
md"""
Increase $\tau_m$ to see how it affects simulations: `τm =` $(@bind τmstr Select(["0.2", "0.4", "0.8", "2", "4", "8"])) μs
"""

# ╔═╡ 237967b2-cd1b-4a77-95b1-d3972f593e2f
md"""
In practice, possibilities for increasing $\tau_m$ are probably bounded by $T_1, T_2$.
"""

# ╔═╡ d842a68e-e5ef-4d93-8c96-da9b5445bdfa
md" ###### Measurement angle"

# ╔═╡ ac1b9f43-39bf-4538-a588-aa2172ba6de6
md"""
We can change the measurement angle to see how it affects the simulation.
"""

# ╔═╡ 7da8d44e-841e-42a3-b42a-0991569fa424
md"""
`φ =` $(@bind φstr Select(["0", "π/8", "π/4", "3π/8", "π/2"]))
"""

# ╔═╡ f8e6a7c8-9d4a-434f-b96a-5f27d5a69b23
md"""
Evidently, stabilization effectiveness is highest for informational backaction. Stabilization fidelity appears to go to zero for $\varphi = \pi/2$, but this should be checked rigorously. 

This raises the question: Can we adapt this method to similarly stabilize phase backaction? This is the topic of the next section.

"""

# ╔═╡ e15a6a43-31fe-488d-aaab-aff9f50d1b8a
md"""
# Phasal linear feedback
"""

# ╔═╡ 6a5deff1-d8ff-47ea-9efa-ecc3b8c6dec0
md"""
Because phase backaction consists of random rotations about $\sigma_z$, the counteracting Hamiltonian term will be a readout-dependent ac Stark shift term. Physically, this involves modulating the cavity drive amplitude, thus limiting backaction on the photon number (which would create the qubit $\sigma_z$ rotations).

This is more complicated, because the Stark shift *and* the measurement timescale will both depend on $\bar n$. A full treatment would look at this in terms of qubit-cavity simulations. However, for simplicity, I will assume feedback variations in $\tau_m$ are small enough that we can model this as a fluctuating Stark shift.

Note that the goal here is not to stabilize to an arbitrary state, but simply to counteract phase backaction. In the absence of Rabi oscillations, this will look like reducing the RMS distance of the final state from the initial state. Including Rabi oscillations, the feedback should confine the state to the Rabi plane.

Then, in a dual-quadrature measurement, the two stabilization techniques could presumably be used in conjunction to increase fidelity which would otherwise be lost due to phase backaction (check it is lost).
"""

# ╔═╡ fe7af10f-ab48-4a49-9bb1-e368ec94d928
md"""
## Pure phase backaction
"""

# ╔═╡ 39679392-2d6a-42f1-89d7-4bd20a953090
md" `ΔS = ` 0 $(@bind ΔS Slider(0:0.1:2)) 2 MHz" 

# ╔═╡ afb6751a-31af-44cf-8922-0f1b816755ca
md" `ΔS = ` $ΔS MHz"

# ╔═╡ 85a3d2ae-f740-43e6-80e3-9ba76c91fb68
@bind rs Button("New random seed")

# ╔═╡ 16ccad7b-1f63-4373-ab81-3f77484a3545
begin
	rs
	seed = abs(rand(Int))
end

# ╔═╡ 3c0e6012-0a8b-4975-baf7-bfe00571933f
md"""
Evidently, $\Delta_S = 1$ MHz leads to perfect erasure of phasal measurement backaction. Now, let's experiment with adding a Rabi drive:
"""

# ╔═╡ 175af755-b2e2-4a3e-9789-855bdec353ae
md" $(@bind rabi CheckBox()) Add Rabi oscillations "

# ╔═╡ ebb1ef94-f30a-4203-8862-5d0294588d30
md"""
The qubit is monitored in quadrature, leading to stochastic record

$\tilde r_Q(t) = \zeta_Q(t)$

where $\zeta(t) \sim \mathcal{N}\big(0, \sqrt{\tau_m/dt} \big)$.

The qubit undergoes Hamiltonian evolution

$\hat H = \frac{\Delta_S}2 \tilde r_Q(t - T_d) \hat \sigma_z + \frac{\Omega_R}2 \hat \sigma_y.$

In the below simulations, the Rabi drive is set to $\Omega_R = 0$ by default. Currently `ΩR = ` $(rabi ? "2π" : 0).
"""

# ╔═╡ 44cfbce3-0f4e-4641-90b8-f4853ccd68ea
md"""
## Dual quadrature measurement
"""

# ╔═╡ c414534a-01d8-422b-9286-ce218a39bee8
idealDQM = true

# ╔═╡ 1dc2ee28-6a9e-4ca3-9c7e-390388882bb5
md" $(@bind info_stabilize CheckBox()) Stabilize target state "

# ╔═╡ 4ada2450-72ce-4960-9278-a6ec125d86a2
md" $(@bind phase_stabilize CheckBox()) Stabilize phase backaction "

# ╔═╡ 4396ff4d-2768-464e-8e90-ddef35cbef5e
md" `ΔS = ` 0.2 $(@bind ΔS2 Slider(0.25:0.005:0.35)) 0.3 MHz"

# ╔═╡ da94dffd-a398-4d72-9742-295f4837f120
md" `ΔS = ` $ΔS2 MHz"

# ╔═╡ dbf13ff0-482b-4f28-9860-e0f1d556450c
@bind rs2 Button("New random seed")

# ╔═╡ f58ca757-4de6-4561-bb4e-4522c51e2b3e
begin
	rs2
	seed2 = abs(rand(Int))
end

# ╔═╡ 1d4b8d69-3c7c-4b9c-9610-489641a619d5
md"""

Rotate view: `ϕv = 0`
$(@bind ϕvDQM html"<input type=range min=0 max=32 step=1 value=0>") 
`ϕv = 2π`

"""

# ╔═╡ 523e81ca-81ea-4d16-9098-7067e6ffa75d
md" $(@bind BlochDQM CheckBox()) Animate Bloch sphere "

# ╔═╡ 9e98e449-859e-4cf7-9e28-b76390d961c9
md"""
# References

 $Patti_et_al T. L. Patti, A. Chantasri, L. P. García-Pintos, A. N. Jordan, and J. Dressel, Linear Feedback Stabilization of a Dispersively Monitored Qubit, Phys. Rev. A 96, 022311 (2017).

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

# ╔═╡ 09bcbe71-e239-4a7c-ac14-99a487c1f9a4
let
	
	# System parameters ------------------------------------------------------------
	# all times given in μs
	
	# initial state 
	x0, y0, z0 = xyz(π/2, 0)
	ρ0 = DenseOperator(0.5*(Iq + x0*σx + y0*σy + z0*σz))
	
	# measurement parameters
	τm = 1.					# time
	Γm = 1/(2τm) 			# rate
	η =  1.					# efficiency
	
	# simulation timescales
	T = (0, 8τm) # simulation duration
	dt = 0.5e-3  # integration time-step
	
	# Hamiltonian parameters
	ΩR = rabi ? 2π : 0
	
	
	# Kraus operators -------------------------------------------------------------
	
	H(t::Timescale, r::Record) = (ΔS/2) * r[1] * σz + (ΩR / 2) * σy
	J = [(σz, ((1-η)*Γm))]
	C = [(im * σz, Γm, η)]
	

	# Bayesian simulation ---------------------------------------------------------
	Random.seed!(seed)
	sol = bayesian(T, ρ0, H, J, C; dt=dt)
	
	global PLF = traj(sol)
	
	md" ###### 🔻 Simulation "
	
end

# ╔═╡ c337e3e0-95cf-4c18-987e-59216bf54419
let
	
	sim = PLF
	
	blochtimeseries(sim.t, sim.x, sim.y, sim.z, title = "Monitored qubit", tf=last(sim.t), ylims=[-1.1,1.1])
	
end

# ╔═╡ 88ec8109-87ed-457d-b9ab-b374176150b1
let
	
	# System parameters ------------------------------------------------------------
	# all times given in μs
	
	# initial state 
	x0, y0, z0 = xyz(π/2, π/2)
	ρ0 = DenseOperator(0.5*(Iq + x0*σx + y0*σy + z0*σz))
	
	# measurement parameters
	τm = 2 					# time
	Γm = 1/(2τm) 			# rate
	η =  1.0				# efficiency
	td = 0.0
	
	T1 = 40 	# energy decay time
	Γ1 = 1/(2T1)# energy decay rate
	T2 = 60 	# environmental dephasing time
	Γ2 = 1/T2 	# environemntal dephasing rate

	# simulation timescales
	T = (0, 25τm) # simulation duration
	dt = 0.5e-3  # integration time-step
	
	# Hamiltonian parameters
	θs = 3π/10 # target angle on Bloch sphere
	ϕ = π # fixes plane of oscillations
	global vecDQM = (sin(θs) * cos(π/2), sin(θs) * sin(π/2), cos(θs))
	σϕ = cos(ϕ)*σx + sin(ϕ)*σy
	
	Δ0 = info_stabilize ? -sin(2θs)/(4τm) : 0
	Δ1 = info_stabilize ? sin(θs)/τm : 0
	ΔS = phase_stabilize ? ΔS2 : 0
	
	
	# Kraus operators -------------------------------------------------------------

	H = 0. * Iq
	# H(t::Time, r::Record) = (Δ0 + Δ1*r[1]) * σϕ/2 + ΔS * r[2] * σz/2 
	J = idealDQM ? [(σz, ((1-η)*Γm))] : [(σz, ((1-η)*Γm + Γ2)), (σm, Γ1)]
	C = [(σz, Γm, η/2), (im * σz, Γm, η/2)]
	

	# Bayesian simulation ---------------------------------------------------------
	Random.seed!(seed2)
	sol = bayesian(T, ρ0, H, J, C; dt=dt)

	global DQM = traj(sol)
	
	md" ###### 🔻 Simulation "
	
end

# ╔═╡ 5615a4a2-768f-47f3-98ff-ae78a4e14413
let
	
	sim = DQM
	vec = vecDQM
	
	blochtimeseries(sim.t, sim.x, sim.y, sim.z, title = "Monitored qubit", tf=last(sim.t), ylims=[-1.1,1.1], vec=vec)
	
end

# ╔═╡ 5c8ee263-061e-41b5-9138-3977e2c6dd09
let 
	ϕv = ϕvDQM * (π/16)
	sim = DQM 
	vec = vecDQM
	
	if BlochDQM
		anim = @animate for i ∈ range(1, length(sim.t), step=100)
			blochsphere(sim.x[1:i], sim.y[1:i], sim.z[1:i], linewidth=1., linealpha=0.85, ax=true, vec=vec, viewϕ = ϕv, blochmark=true) end
		gif(anim, fps = 15)
	else
		blochsphere(sim.x, sim.y, sim.z, linewidth=1., linealpha=0.85, ax=true, vec = vec, viewϕ = ϕv, blochmark=true)
	end
end

# ╔═╡ 660d4a39-1818-4f27-8371-acdbae557b97
let
	
	# System parameters ------------------------------------------------------------
	# all times given in μs
	
	# initial state 
	x0, y0, z0 = xyz(π/2, π/2)
	ρ0 = DenseOperator(0.5*(Iq + x0*σx + y0*σy + z0*σz))
	
	# measurement parameters
	τm = 1 					# time
	Γm = 1/(2τm) 			# rate
	η =  0.5			# efficiency
	
	# simulation timescales
	T = (0, 25τm) # simulation duration
	dt = 0.5e-3  # integration time-step
	
	# Hamiltonian parameters
	θs = 3π/10 # target angle on Bloch sphere
	ϕ = π # fixes plane of oscillations
	global vecSQM = (sin(θs) * cos(π/2), sin(θs) * sin(π/2), cos(θs))
	σϕ = cos(ϕ)*σx + sin(ϕ)*σy
	
	
	Δ0 = info_stabilize ? -sin(2θs)/(4τm) : 0
	Δ1 = info_stabilize ? sin(θs)/τm : 0
	
	
	# Kraus operators -------------------------------------------------------------
	
	H(t, r) = (Δ0 + Δ1*r[1]) * σϕ/2
	J = [(σz, ((1-η)*Γm))]
	C = [(σz, τm, η)]
	

	# Bayesian simulation ---------------------------------------------------------
	Random.seed!(seed2)
	sol = bayesian(T, ρ0, H, J, C; dt=dt)

	global SQM = traj(sol)
	
	md" ###### 🔻 Single-quadrature comparison "
	
end

# ╔═╡ c1ac49f3-c6ba-4f89-b351-49cfee9bb8f8
let
	
	sim = SQM
	vec = vecSQM
	
	blochtimeseries(sim.t, sim.x, sim.y, sim.z, title = "Monitored qubit", tf=last(sim.t), ylims=[-1.1,1.1], vec=vec)
	
end

# ╔═╡ 618168dc-53dd-4562-8061-67a0b56587aa
φdict = Dict("0" => 0, 
				"π/8" => π/8, 
				"π/4" => π/4,
				"3π/8" => 3π/8,
				"π/2" => π/2)

# ╔═╡ 13e19d59-5132-45f3-bac1-3a3b3a7a12b2
let
	
	# System parameters ------------------------------------------------------------
	# all times given in μs
	
	# initial state 
	x0 = 0
	y0 = 0.5 # 0.3
	z0 = sqrt(1 - y0^2)
	ρ0 = DenseOperator(0.5*(Iq + x0*σx + y0*σy + z0*σz))
	
	
	# measurement parameters
	
	η =  parse(Float64, ηstr)			# efficiency
	Γm = 2. # parse(Float64, Γmstr)			# measurement dephasing rate
	τm = 1/(2Γm) 						# measurement collapse time
	φ = get(φdict, φstr, 0) 			# angle
	td = parse(Float64, Tdstr) * 1e-3 	# time delay for feedback

	
	# simulation timescales
	
	T = (0, 4) # simulation duration
	dt = 0.5e-3  # integration time-step
	
	
	# feedback target parameters
	
	global θs = 3π/10 		# target angle on Bloch sphere
	Rs = 1 # 0.64 			# radius of target
	ϕ = π 					# fixes plane of oscillations
	global vec = (sin(θs) * cos(π/2), sin(θs) * sin(π/2), cos(θs))
	σϕ = cos(ϕ)*σx + sin(ϕ)*σy
	
	
	# feedback drive parameters
	
	Δ0 = ideal ? -sin(2θs)/(4τm) : -sin(2θs)/(4τm*Rs^2) 
	Δ1 = ideal ? sin(θs)/τm : sin(θs)/(τm*Rs)
	
	
	# decay times and rates
	
	T1 = 40 	# energy decay time
	Γ1 = 1/(2T1)# energy decay rate
	T2 = 60 	# environmental dephasing time
	Γ2 = 1/T2 	# environemntal dephasing rate
	
	
	# Kraus operators -------------------------------------------------------------

	H(t::Timescale, r::Readout) = (Δ0 + Δ1*r[1]) * σϕ/2
	J = ideal ? [(σz, ((1-η)*Γm))] : [(σz, ((1-η)*Γm)), (σm, Γ1), (σz, Γ2)]
	C = [(exp(im * φ) * σz, Γm, η)]
	
	
	# Bayesian simulation ---------------------------------------------------------
	sol = bayesian(T, ρ0, H, J, C; dt=dt, td=td)
	global ILF = traj(sol)
		
	md" ###### 🔻 Simulation "
	
end

# ╔═╡ c9d3c54e-acff-4a84-8f95-937ee1602350
md"""
In this simulation, the stabilization target is chosen to be $\theta_s =$ $(round(θs/π, digits=3)) π in the plane of Rabi oscillations. In other words, the state will be stabilized to the point (x, y, z) = $(round.(vec, digits=3)). 

In plots, the target is indicated by a red vector or dashed lines.
"""

# ╔═╡ 725dc4c3-cc74-4400-819c-2cffd06fbbf9
let 
	ϕv = ϕvc * (π/16)
	sim = ILF
	
	if show_gif
		anim = @animate for i ∈ range(1, length(sim.t), step=100)
			blochsphere(sim.x[1:i], sim.y[1:i], sim.z[1:i], linewidth=1., linealpha=0.85, ax=true, vec=vec, viewϕ = ϕv, blochmark=true) end
		gif(anim, fps = 15)
	else
		blochsphere(sim.x, sim.y, sim.z, linewidth=1., linealpha=0.85, ax=true, vec = vec, viewϕ = ϕv, blochmark=true)
	end
end

# ╔═╡ 34a700bb-5809-4755-a7fa-def102c5fd4c
let 
	sim = ILF
	
	if show_gif2
		anim = @animate for i ∈ range(1, length(tt), step=100)
			blochtimeseries(sim.t[1:i], sim.x[1:i], sim.y[1:i], sim.z[1:i], vec=vec, title = "Monitored Rabi oscillations", tf=last(sim.t)) end
		gif(anim, fps = 15)
	else
		blochtimeseries(sim.t, sim.x, sim.y, sim.z, vec=vec, title = "Monitored Rabi oscillations", tf=last(sim.t), ylims=[-1.1,1.5])
	end
end

# ╔═╡ bb5f3187-2773-4647-807a-63141e16c2b4
let 
	sim = ILF
	
	if show_gif3
		anim = @animate for i ∈ range(1, length(sim.t), step=100)
			blochprojections(sim.x[1:i], sim.y[1:i], sim.z[1:i], blochmark=true, blochmarkcolor="white", vec=vec) end
		gif(anim, fps = 15)
	else
		blochprojections(sim.x, sim.y, sim.z, blochmark=true, vec=vec, blochmarkcolor="white")
	end
end



# ╔═╡ 9d98ff41-28d5-4f3d-b672-49131dafd5db
md"""
## Paths
"""

# ╔═╡ fef481d9-5151-4d1e-81eb-41231d0febe6
refpath = string(path, "notebooks/1-quantum-mechanics-intro/4-measurement-backaction.jl")

# ╔═╡ Cell order:
# ╟─377a3336-20bd-4baa-a033-af8bbc8668a8
# ╟─e77d8c4b-73fc-4ba1-8054-f5da97a36748
# ╟─9014ea2e-c244-4aeb-bc56-811dcfb6b70c
# ╟─3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─ab188d7c-9eb4-48f6-86e6-53c057ce67ae
# ╟─5b83ae62-9486-485f-821b-eb127f01e487
# ╟─1eb3c79e-2efb-4a4e-9fb2-07d78e9cf21a
# ╟─d724fb1e-bd93-4a69-8082-eadaa50b31a0
# ╟─34e1740c-20d6-4099-ab7a-0c64160957c4
# ╟─7ff6b7dc-ba8b-4ce3-8a1c-0308e871986d
# ╟─7fb819fc-f0ae-49ac-a702-8b5329adccd3
# ╟─18c9725c-ebad-49b8-841c-57893558251f
# ╟─1edc57b3-de6a-4546-abab-78b780922081
# ╟─e338a516-c4c8-4571-8e84-d85dcc68c985
# ╟─a7a38426-33b0-4057-b125-0488d4bcb279
# ╟─a87a6bce-c5b6-4b04-9d14-d9656969d87b
# ╟─fb578156-094b-43b1-8c64-740f64b193bc
# ╟─c9d3c54e-acff-4a84-8f95-937ee1602350
# ╟─5e76a75a-3041-49f6-b5b6-96047b1d5bd4
# ╠═b5ddea46-9aee-4ce9-8230-902611228ca7
# ╠═13e19d59-5132-45f3-bac1-3a3b3a7a12b2
# ╟─e486fbed-157c-4139-914c-ada6bb88d7b4
# ╟─4898ab97-4058-4e70-a959-2962641d9611
# ╟─6b93bb84-0ba1-44e8-910e-612793b3df3b
# ╟─ce79057b-eb31-49b1-9db7-bfd793a447c5
# ╟─ca7a2351-cff5-4d77-ba16-00de384b8b7c
# ╟─f3e7794e-a9a1-4103-8469-32a8e37d2d82
# ╟─a12cdb8c-e9a1-4c2d-9811-cff266e152d8
# ╟─725dc4c3-cc74-4400-819c-2cffd06fbbf9
# ╟─0a7f28c9-1d84-43e6-b62c-711a231a3972
# ╠═34a700bb-5809-4755-a7fa-def102c5fd4c
# ╟─d9f2f00f-4ee2-45b5-91d6-6552d6d5b6c1
# ╟─bb5f3187-2773-4647-807a-63141e16c2b4
# ╟─e28cc7a9-ccc7-418a-94fa-068befe4f151
# ╟─3679d9db-97aa-40cd-ab1b-2329051c7156
# ╟─1159cb26-2a33-46af-a2cb-3cc667b79c0d
# ╟─08f808ab-13f7-4799-951c-72042d71e1be
# ╟─237967b2-cd1b-4a77-95b1-d3972f593e2f
# ╟─d842a68e-e5ef-4d93-8c96-da9b5445bdfa
# ╟─ac1b9f43-39bf-4538-a588-aa2172ba6de6
# ╟─7da8d44e-841e-42a3-b42a-0991569fa424
# ╟─f8e6a7c8-9d4a-434f-b96a-5f27d5a69b23
# ╟─e15a6a43-31fe-488d-aaab-aff9f50d1b8a
# ╟─6a5deff1-d8ff-47ea-9efa-ecc3b8c6dec0
# ╟─fe7af10f-ab48-4a49-9bb1-e368ec94d928
# ╟─ebb1ef94-f30a-4203-8862-5d0294588d30
# ╠═09bcbe71-e239-4a7c-ac14-99a487c1f9a4
# ╠═39679392-2d6a-42f1-89d7-4bd20a953090
# ╟─afb6751a-31af-44cf-8922-0f1b816755ca
# ╟─c337e3e0-95cf-4c18-987e-59216bf54419
# ╟─85a3d2ae-f740-43e6-80e3-9ba76c91fb68
# ╟─16ccad7b-1f63-4373-ab81-3f77484a3545
# ╟─3c0e6012-0a8b-4975-baf7-bfe00571933f
# ╟─175af755-b2e2-4a3e-9789-855bdec353ae
# ╟─44cfbce3-0f4e-4641-90b8-f4853ccd68ea
# ╠═c414534a-01d8-422b-9286-ce218a39bee8
# ╠═88ec8109-87ed-457d-b9ab-b374176150b1
# ╟─1dc2ee28-6a9e-4ca3-9c7e-390388882bb5
# ╟─4ada2450-72ce-4960-9278-a6ec125d86a2
# ╟─4396ff4d-2768-464e-8e90-ddef35cbef5e
# ╟─da94dffd-a398-4d72-9742-295f4837f120
# ╠═5615a4a2-768f-47f3-98ff-ae78a4e14413
# ╟─dbf13ff0-482b-4f28-9860-e0f1d556450c
# ╟─f58ca757-4de6-4561-bb4e-4522c51e2b3e
# ╟─1d4b8d69-3c7c-4b9c-9610-489641a619d5
# ╠═523e81ca-81ea-4d16-9098-7067e6ffa75d
# ╠═5c8ee263-061e-41b5-9138-3977e2c6dd09
# ╟─660d4a39-1818-4f27-8371-acdbae557b97
# ╠═c1ac49f3-c6ba-4f89-b351-49cfee9bb8f8
# ╟─9e98e449-859e-4cf7-9e28-b76390d961c9
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
# ╟─9d98ff41-28d5-4f3d-b672-49131dafd5db
# ╠═fef481d9-5151-4d1e-81eb-41231d0febe6
# ╠═77d623fb-7dbc-406c-bf46-05385a049530
