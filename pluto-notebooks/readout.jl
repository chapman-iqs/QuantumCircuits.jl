### A Pluto.jl notebook ###
# v0.17.1

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

# ╔═╡ 4c45fe3e-cd69-11eb-20cd-7bfb98c040cf
begin
	cd("/Users/sachagreenfield/Desktop/GitHub/QuantumCircuits.jl")
	import Pkg
	Pkg.activate(".")
	
	using PlutoUI
	using LaTeXStrings
	using Random
	using Statistics
	using Distributions
	using QuantumCircuits
	using Plots
	using Plots.Measures
	
	include("plotting.jl")
	
	md" ### 🔶 Packages and julia files"
end

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents()

# ╔═╡ 377a3336-20bd-4baa-a033-af8bbc8668a8
md"""
# Readout

In this interactive notebook, we'll understand how readout works and how it affects quantum state evolution.
"""

# ╔═╡ 3be820e3-8541-424c-b042-facf686f644b
md"""
# Rabi flopping
"""

# ╔═╡ 66d52a3a-edfd-4bdd-8c18-9f1a27041876
md"""
To better observe the effects of continuous measurement, we'll look at a so-called "Rabi flopping" regime, where the qubit is being driven at a significant Rabi rate ($\Omega_R = 3$ MHz), while also being measured relatively strongly at rate $\Gamma = 40$ MHz. 

The corresponding measurement collapse timescale -- the time the state takes to collapse to a measurement eigenstate -- is $\tau = 1/2\Gamma = 0.0125$ μs. Thus, the state will generally be pinned to an eigenstate; however, the Rabi drive allows for quick transitions or "jumps" to the opposite eigenstate.
"""

# ╔═╡ 282851f0-fa1d-42b1-a243-5620339129cc
md"""
### Simulation
"""

# ╔═╡ 259978cd-7888-4f38-bf83-e40719d9777b
md"""
You can select from single-quadrature or dual quadrature measurement below, and observe the effects of changing the measurement angle $\varphi$.
"""

# ╔═╡ 8d749aee-829a-479c-ba65-b66b2d18bee3
md"""
When measuring two quadratures, the efficiency on each is halved, and the measurement time is doubled. **check this -- should input parameters to Bayesian update really be 2τ, or is this already captured in η/2?**
"""

# ╔═╡ ad636be5-a60b-4f36-9470-e33d5d3c891a
md"""
Change the measurement angle below:
"""

# ╔═╡ 2662b242-1d39-4223-a600-8b5a5d6a148b
md"""

ϕ = 0 
$(@bind ϕcd html"<input type=range min=0 max=16 step=1 value=0>") 
ϕ = π/2

"""

# ╔═╡ db19a712-aa39-4691-93dc-bd9d18f5db7e
begin
	ϕ = ϕcd * (π/32)
	md" ϕ = $(ϕcd/32) π "	
end

# ╔═╡ cf5721a2-6f73-47e9-9b11-1c9562c29a03
md"""

`tf = 0`
$(@bind tfactor html"<input type=range min=0.05 max=1 step=0.05 value=1>") 
`tf = 2 μs`

"""

# ╔═╡ d05919e8-48c5-40dd-95fe-151d793049d7
md" $(@bind dualquad CheckBox()) Dual quadrature "

# ╔═╡ dbf9e12b-06b0-42cb-968b-d49f7d1c4e10
md"""
In the above plot, the upper time series represents the readout with arbitrary y-axis units, and time as given in the lower plot. The circle to the right indicates the measurement angle. Purple (gray) indicates the measurement is informational (non-informational). The lower time series is the bloch coordinate evolution. 
"""

# ╔═╡ 0396d908-8123-49db-a1d6-a637651f81b5
md"""
## Measurement record sampling
"""

# ╔═╡ 1105dac1-61ff-44e9-a82a-fc9823ea7757
md"""
When we do readout, we are sampling from a distribution centered around the expected value of $z$. We then use this information to update our estimate of the quantum state.
"""

# ╔═╡ 6145101c-9ce3-4938-bf10-83fa14ff84db
md"""
The records plotted above are sampled in this way. Thus, as $z$ evolves, the mean of the distribution shifts slightly along with the quantum state. Depending on the strength of measurement, this shift can be almost imperceptible compared to the variance. However, this small shift accounts for state-dependent backaction.
"""

# ╔═╡ 49debdd0-b3db-41c3-b6e0-a42750787fed
md"""
You can play with the plot below to understand the relationship between the semi-random readout and the quantum state.
"""

# ╔═╡ f21f2ca8-83a9-4b78-8053-f6cfc9536338
md"""

`t = 0`
$(@bind ti html"<input type=range min=2 max=1001 step=1 value=1>") 
`t = 1 μs`

"""

# ╔═╡ 3c55abb4-5114-4183-ba06-d12e2f4d8a01
md"""
The first time series is the readout -- what is observed in the lab, up to shifting and rescaling. The Gaussian curve is the distribution the data are sampled from, which moves with the expectation value of $z$. The green time series represents the readout's tracking of the $z$ coordinate, compared to the black series, which is the true $z$ value.
"""

# ╔═╡ a968952f-a64a-4882-a151-8a2d9a1723c7
md"""
**In the above, should we expect tracking efficacy to be halved for dual-quadrature measurement? Or is this already captured in $\tau \mapsto 2\tau$?**
"""

# ╔═╡ d0957eb9-0df7-4015-9720-5fdb4d1df3b3
md"""
This demonstrates how the quantum state evolution affects the sampling of the measurement record. However, the connection goes both ways: the measurement itself affects quantum state evolution. To see this in action, we will next consider Bayes' rule and its application in Bayesian update.
"""

# ╔═╡ 37dac3f2-f394-467c-ba4f-cf3f9d6d9577
md"""
# Bayesian update
"""

# ╔═╡ 553dec63-d7c0-4105-b5dc-f2c36b11eafe
md"""
Bayesian update describes how the quantum state updates based on a measurement outcome. Essentially, it is responsible for measurement backaction. It is founded upon Bayes rule.
"""

# ╔═╡ 32ee829c-6a9c-4ae8-a8ed-8d3516d1afec
md"""
## Bayes' rule
"""

# ╔═╡ 5b40fc7a-4280-49f2-b2f3-000e9cacc947
md"""
In its canonical formulation, Bayes rule is written

$P(A | B) = \frac{P(B | A) \cdot P(A)}{P(B)},$

where $A$ and $B$ are two events. It reads: the probability of $A$, *given* $B$ is true -- $P(A | B)$ -- is proportional to the probability of $B$, *given* $A$ is true -- $P(B | A)$. The scaling factor is the ratio of the independent probabilities of each event: $P(A) / P(B)$.

Suppose we are keeping track of the probability of $A$ over time. Then as we gain information about $A$ *indirectly* -- through gaining information about $B$, which is *correlated* with $A$ -- we can update our estimate of $A$. More precisely, this temporal formulation can be written

$P(A_{t + \Delta t}) = P(A_{t} | B_{t}) = \frac{P(B_t | A_t) \cdot P(A_t)}{P(B_t)},$

where the first equality assumes that gaining information via $B$ is the only thing affecting $A$.

In the case of qubit readout, the event $A$ is measuring $\ket{0}$. The event $B$ is obtaining a measurement outcome via our weak measurement. Thus,

* ` $P(A_{t + \Delta t}) = P(A_t | B_t)$ is our updated estimate of $P(\ket{0}) = 0.5(1 + z)$, based on our measurement of $B_t$

* ` $B_t$ is our measurement outcome, $\tilde{r} (t)$, and $P(B_t)$ the unconditioned probability of measuring it

* ` $A_t$ is our prior estimate of $P(\ket{0})$

* ` $P(B_t | A_t)$ is the conditional probability of measuring $\tilde{r}(t)$ given $\ket{0}$ at time $t$

"""

# ╔═╡ 34a2807c-bd8c-4563-a616-10e0ca3d2e52
md"""
## Quantum mechanical formalism
"""

# ╔═╡ 9a9c5acb-6db5-4a71-8174-75a64d7cb2d0
md"""
### Simulation: pure measurement collapse
"""

# ╔═╡ dbf58f2c-9926-4a57-a0f4-0c8a1332df80
md"""
As above, we will assume that the measurement is the only thing affecting the quantum state. Then we can consider the evolution due purely to measurement backaction. This will make the effect clearer.

Since Bayesian update makes no assumption of small $dt$, we can take $dt$ to be whatever we want -- in this case, we use the common experimental value of $dt = 40$ ns.
"""

# ╔═╡ b5c18b03-2543-4d47-a488-93b2c43663fd
md"""
We can see that measurement collapses the state to $z = 1$, i.e. to the state $\ket{0}$, in a finite amount of time: approximately $1 $ μs $= 4\tau$ (recall $2\tau$ is approximately the amount of time we expect for collapse). Now, let's pull apart what is going on here.
"""

# ╔═╡ 7e444185-4e30-412a-9a38-cadb547e8372
md"""
##### t = 0

Our initial state is $\ket{\psi}_0 = \frac{1}{\sqrt{2}}(\ket{0} + \ket{1}).$ So
* ` $P_0(\ket{0}) = 0.5$
* ` $P_0(\tilde{r}) \sim \mathcal{N}(0, \sqrt{\tau/dt}) \Big|_{\tilde r}$
* ` $P_0(\tilde{r} \big| \ket{\psi}_0) \sim \mathcal{N}(z_0, \sqrt{\tau/dt}) = \mathcal{N}(0, \sqrt{\tau/dt})  \Big|_{\tilde r}$

"""

# ╔═╡ b4068526-178a-42d9-994d-e4e3c879bb6f
length(range(0,1,step=0.05)) .* 1e-3

# ╔═╡ 133b6939-10b2-4c8e-acf8-5658ca96a0f9
md" # Utilities"

# ╔═╡ c1d743fb-356c-4d7a-a290-110c658e20dd
getclosest(array, val) = argmin(abs.(val .- array))

# ╔═╡ a2be6ad8-1ed5-4869-9741-106a348a9d82
colors = palette(:tab10)

# ╔═╡ c63be186-2510-448b-a1e6-fb4c9b96b28f
begin
	phase(φ) = exp(im * φ)
	angles = range(0, 2.2π, step=2π/100)
	phases = phase.(angles)
end

# ╔═╡ 5cef3447-784d-4be1-8e0e-44737997b4ea
begin
	@userplot BlochReadout3
	
	@recipe function f(bts::BlochReadout3; vec=nothing, tf=nothing)
		ϕ, sim = bts.args
		ts, xs, ys, zs, ps, r = sim.t, sim.x, sim.y, sim.z, sim.p, sim.r
	
		# Plot time series --------------------------------------------------------
		
		legend := [:none :none :bottomright]
		label --> [L"$x$" L"$y$" L"$z$" L"$Tr ( \rho^2 )$"]
		xlabel --> "t (μs)"
	    link := :both
	    framestyle := [:none :none :axes]
	    grid := false
	    layout := @layout [readout1{1.0w, 1.0h}  _
							readout2{1.0w, 1.0h}  _
	                       blochseries ]
	
		palette := :tab10
		linealpha --> 1
	
		legendfontsize --> 10
		titlefontsize --> 12
		xtickfontsize --> 10
		ytickfontsize --> 10
		xguidefontsize --> 10
		yguidefontsize --> 10
		size --> (500,400)
		linewidth --> 1.5
	
		tf = (tf == nothing) ? last(ts) : tf
		xlims --> [first(ts), tf]
		ylims --> [-1, 1]
		
	
	
		for bs in (xs, ys, zs, ps)
	
			@series begin
				# top_margin := -5mm
				right_margin := 10mm
				subplot := 3
				ts, bs
			end
	
		end
		
		
		ylims := [minimum(r[1]), maximum(r[1])]
		linewidth := 1
		
		
		color := [RGB(200/255,abs(sin(ϕ))*200/255,200/255) RGB(200/255,abs(cos(ϕ))*200/255,200/255)]
		
		
		@series begin
			top_margins := -20mm
			bottom_margins := -20mm
			label := :none
			subplot := 1
			ts, r[1]
		end
		
		if length(r) == 1
			
			@series begin
				top_margins := -20mm
				bottom_margins := -20mm
				label := :none
				linealpha := 0
				subplot := 2
				ts, r[1]
			end
			
		else
			
			@series begin
				top_margins := -20mm
				bottom_margins := -20mm
				label := :none
				subplot := 2
				ts, r[2]
			end
			
		end
		
		
		
			
		
	

	
	end
	
	md" 💧 Bloch readout macro"
end

# ╔═╡ 2546f296-ef35-4702-9dd1-13c1b308529a
begin
	@userplot BlochReadout2
	
	@recipe function f(bts::BlochReadout2; vec=nothing, tf=nothing)
		ϕ, sim = bts.args
		ts, xs, ys, zs, ps, r = sim.t, sim.x, sim.y, sim.z, sim.p, sim.r
	
		# Plot time series --------------------------------------------------------
		
		legend := [:none :none :bottomright]
		label --> [L"$x$" L"$y$" L"$z$" L"$Tr ( \rho^2 )$"]
		xlabel --> "t (μs)"
	    link := :both
	    framestyle := [:none :none :axes]
	    grid := false
	    layout := @layout [readout1{1.0w, 1.0h}  _
							readout2{1.0w, 1.0h}  _
	                       blochseries ]
	
		palette := :tab10
		linealpha --> 1
	
		legendfontsize --> 10
		titlefontsize --> 12
		xtickfontsize --> 10
		ytickfontsize --> 10
		xguidefontsize --> 10
		yguidefontsize --> 10
		size --> (500,400)
		linewidth --> 1.5
	
		tf = (tf == nothing) ? last(ts) : tf
		xlims --> [first(ts), tf]
		ylims --> [-1, 1]
		
	
	
		for bs in (xs, ys, zs, ps)
	
			@series begin
				# top_margin := -5mm
				right_margin := 10mm
				subplot := length(r) + 1
				ts, bs
			end
	
		end
		
		
		ylims := [minimum(r[1]), maximum(r[1])]
		linewidth := 0.25
		
		
		color := [RGB(200/255,abs(sin(ϕ))*200/255,200/255) RGB(200/255,abs(cos(ϕ))*200/255,200/255)]
		
	
		for (i, rs) in enumerate(r)
		
			@series begin
				top_margins := -20mm
				bottom_margins := -20mm
				label := :none
				subplot := i
				ts, r[i]
			end
			
		end
	

	
	end
	
	md" 💧 Bloch readout macro"
end

# ╔═╡ 3f0bf48d-a27d-492c-b2d7-2291222c4607
colorangle(ϕ) = RGB(200/255,abs(sin(ϕ))*200/255,200/255)

# ╔═╡ a427e58a-fa75-4d50-acab-da11e49dd1bf
function plot_phases3(sim, ϕ, tfactor)
	
	index = Int64(floor(length(sim.t) * tfactor))
	
	
	br = blochreadout3(ϕ, sim, tf=sim.t[index], ylims=[-1.1,1.5])
	
	plot!([(real(phases), imag(phases)), ([0, cos(ϕ)], [0, sin(ϕ)])],  inset = bbox(0, 0.1, 0.1, 0.1, :right), subplot=4, legend = :none, frame = :none, aspect_ratio = :equal, color = [:black colorangle(ϕ)], linewidth=[1 2])
	
	
	plot!([(real(phases), imag(phases)), ([0, cos(ϕ + π/2)], [0, sin(ϕ + π/2)])],  inset = bbox(0, 0.38, 0.1, 0.1, :right), subplot=5, legend = :none, frame = :none, aspect_ratio = :equal, color = [:black colorangle(ϕ + π/2)], linewidth=[1 2])
		

	
end
	

# ╔═╡ 2d2e530d-5e8a-4045-93bd-97ff99aecaa1
begin
	@userplot BlochReadout
	
	@recipe function f(bts::BlochReadout; vec=nothing, tf=nothing)
		ϕ, sim = bts.args
		ts, xs, ys, zs, ps, r = sim.t, sim.x, sim.y, sim.z, sim.p, sim.r
	
		# Plot time series --------------------------------------------------------
		
		legend := [:right :bottomright]
		label --> [L"$x$" L"$y$" L"$z$" L"$Tr ( \rho^2 )$"]
		xlabel --> "t (μs)"
	    link := :both
	    framestyle := [:none :axes]
	    grid := false
	    layout := @layout [readout{1.0w, 1.0h}  _
	                       blochseries ]
	
		palette := :tab10
		linealpha --> 1
	
		legendfontsize --> 10
		titlefontsize --> 12
		xtickfontsize --> 10
		ytickfontsize --> 10
		xguidefontsize --> 10
		yguidefontsize --> 10
		size --> (400,300)
		linewidth --> 1.5
	
		tf = (tf == nothing) ? last(ts) : tf
		xlims --> [first(ts), tf]
		ylims --> [-1, 1]
		
		top_margin := -20mm
		
	
	
		for bs in (xs, ys, zs, ps)
	
			@series begin
				right_margin := 10mm
				subplot := 2
				ts, bs
			end
	
		end
		
		
		ylims := [minimum(r[1]), maximum(r[1])]
		linewidth := 0.5
		
	
		
		@series begin
			label := :none #string("ϕ = ", ϕ)
			color := colorangle(ϕ)
			subplot := 1
			ts, r[1]
		end
	

	
	end
	
	md" 💧 Bloch readout macro"
end

# ╔═╡ 57e71783-8e44-453e-a1d9-42377ffe41e9
function plot_phases(sim, ϕ, tfactor)
	
	index = Int64(floor(length(sim.t) * tfactor))
	
	
	br = (length(sim.r) == 1) ? 
			blochreadout(ϕ, sim, tf=sim.t[index], ylims=[-1.1,1.5]) :
			blochreadout2(ϕ, sim, tf=sim.t[index], ylims=[-1.1,1.5])
	
	if length(sim.r) == 1
		
		plot!([(real(phases), imag(phases)), ([0, cos(ϕ)], [0, sin(ϕ)])],  inset = bbox(0, 0.15, 0.12, 0.12, :right), subplot=3, legend = :none, frame = :none, aspect_ratio = :equal, color = [:black colorangle(ϕ)], linewidth=[1 2])
		
	else

		plot!([(real(phases), imag(phases)), ([0, cos(ϕ)], [0, sin(ϕ)])],  inset = bbox(0, 0.1, 0.1, 0.1, :right), subplot=4, legend = :none, frame = :none, aspect_ratio = :equal, color = [:black colorangle(ϕ)], linewidth=[1 2])

		plot!([(real(phases), imag(phases)), ([0, cos(ϕ + π/2)], [0, sin(ϕ + π/2)])],  inset = bbox(0, 0.38, 0.1, 0.1, :right), subplot=5, legend = :none, frame = :none, aspect_ratio = :equal, color = [:black colorangle(ϕ + π/2)], linewidth=[1 2])
		
	end
	
end
	

# ╔═╡ e43e5329-bd96-41ce-a183-1bd206204f65
begin
	Nfock = 15
	
	# Basis
	q = SpinBasis(1//2)
	f = FockBasis(Nfock)
	Iq = identityoperator(q)
	If = identityoperator(f)

	# qubit operators, using convention that |-z> is ground state
	σx = sigmax(q)
	σy = sigmay(q)
	σz = sigmaz(q)
	σp = sigmap(q)
	σm = sigmam(q)
	
	
	ground = spindown(q)
	excited = spinup(q)


	# Qubit-resonator operators
	Id = Iq ⊗ If
	a = Iq ⊗ destroy(f)
	
	X = σx ⊗ If
	Y = σy ⊗ If
	Z = σz ⊗ If
	P = σp ⊗ If
	M = σm ⊗ If
	
	# projectors
	ZP = dm(excited) ⊗ If
	ZM = dm(ground) ⊗ If
	
	αP = a * ZP
	αM = a * ZM
	
	md"### 🔶 Hilbert space operators "
end

# ╔═╡ e7b46d6f-2f2e-4b60-8684-ee65adcf1431
begin
	ρ0 = dm((ground + excited)/√2) # initial state
	dt = 1e-3  # integration time-step

	ΩR  = 2π * 3 # Rabi frequency (rad * MHz)
	Γ = 40 #0.15 # Measurement dephasing rate (MHz)
	τ = 1/(2Γ) #  Measurement collapse timescale (μs)
	η = 1 # collection efficiency
end

# ╔═╡ 6843e772-74b5-48e1-95c8-5648a93f0787
begin
	fnorm(μ, σ) = x ->  1/(σ * √(2π)) * exp(-0.5 * (x - μ)^2  / σ^2)
	fun(z, τ) = fnorm(z, √(τ/dt))
	md" 🌀 Define normal distribution functions"
end

# ╔═╡ ea04744b-4296-4dc2-9a3c-1f477c96f1ac
md"""
### Plotting
"""

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

# ╔═╡ 4e50b40b-82fc-4eac-a96c-ac10546eecfe
let
	
	# Parameters -------------------------------------------------------------------
	
# 	ρ0 = dm((ground + excited)/2) # initial state
# 	dt = 1e-3  # integration time-step
	
# 	ΩR  = 0 # Rabi frequency (rad * MHz)
# 	Γ = 0.15 # Measurement dephasing rate (MHz)
# 	τ = 1/(2Γ) #  Measurement collapse timescale (μs)
# 	η = 1 # collection efficiency
	
	# Kraus operators --------------------------------------------------------------
	
	H = ΩR * σy/2
	J = [(σz, ((1-η) * Γ))]
	C = [(exp(im * ϕ) * σz, τ, η)]
	
	Random.seed!(1)
	sol = bayesian((0, 2), ρ0, H, J, C; dt=dt)
	
	global Rsq = traj(sol)
	

	
	md" ###### 🔻 Bayesian simulation (single-quadrature)"
end

# ╔═╡ 9b91683a-e5a4-41e9-81fb-b84cb781cd36
let
	
	# Kraus operators --------------------------------------------------------------
	
	H = ΩR * σy/2
	J = [(σz, ((1-η) * Γ))]
	C = [(exp(im * ϕ) * σz, 2τ, η/2), (exp(im * (ϕ + π/2)) * σz, 2τ, η/2)]
	
	Random.seed!(1)
	sol = bayesian((0, 2), ρ0, H, J, C; dt=dt)
	
	global Rdq = traj(sol)
	

	
	md" ###### 🔻 Bayesian simulation (dual-quadrature)"
end

# ╔═╡ fac6aa1f-28ec-4e67-b513-6ba485e88eee
begin
	sim = dualquad ? Rdq : Rsq
	plot_phases3(sim, ϕ, tfactor)
end

# ╔═╡ e2443859-dda5-48b4-bf3f-312deeff7a8e
tf = getclosest(sim.t, 1)

# ╔═╡ 082f9721-a5eb-4b78-93cc-a1dc9cc20ad4
md" `t = ` $(sim.t[ti]) μs"

# ╔═╡ 64fef35f-b894-48db-9ad3-9939c57fcc66
if dualquad

	let

		l = @layout [readout{0.65w} righthist{0.25w} zseries{0.1w}]
		ts = sim.t[1:ti] * 1000
		rs1 = sim.r[1][1:ti]
		rs2 = sim.r[2][1:ti]
		zs = sim.z[1:ti]
		
		values1 = range(floor(minimum(sim.r[1])) - 1, floor(maximum(sim.r[1])), step=0.1)
		values2 = range(floor(minimum(sim.r[2])) - 1, floor(maximum(sim.r[2])), step=0.1)
		
		density1 = fun(last(zs) * cos(ϕ), 2τ).(values1)
		density2 = fun(last(zs) * cos(ϕ + π/2), 2τ).(values2)


		p1 = plot(ts, rs1, xlims=[first(sim.t), sim.t[tf] * 1000], color = colorangle(ϕ), xlabel="t (ns)", ylabel="arbitrary units", left_margin=1mm)
		plot!(ts, rs2, color = colorangle(ϕ + π/2))
		plot!([last(ts)], [last(rs1)], marker=true, color = colorangle(ϕ))
		plot!([last(ts)], [last(rs2)], marker=true, color = colorangle(ϕ + π/2))

		p2 = plot(density1, values1, frame = :none, color = colorangle(ϕ), linewidth=2)
		ind1 = getclosest(values1, last(rs1))
		plot!(density2, values2, frame = :none, color = colorangle(ϕ + π/2), linewidth=2)
		ind2 = getclosest(values2, last(rs2))
		indz1 = getclosest(values1, last(zs))
		indz2 = getclosest(values2, last(zs))
		plot!([density1[ind1]], [last(rs1)], marker=true, color = colorangle(ϕ))
		plot!([density2[ind2]], [last(rs2)], marker=true, color = colorangle(ϕ + π/2))
		plot!([density1[indz1]], [last(zs) * cos(ϕ)], marker=true, color = colors[3])
		plot!([density2[indz2]], [last(zs) * cos(ϕ + π/2)], marker=true, color = colors[3])


		scalefactor = last(ts)/(sim.t[tf] * 1000)
		p3 = plot(-ts, zs, color=:black, linewidth=0.5)
		plot!(-ts, zs .* cos(ϕ), color=colors[3], frame=:true, ymirror=:true, framestyle = :origin, yticks=[-1.0, 0, 1.0], xticks=:none, showaxis=:hide, right_margin=5mm, annotations=[(-450 * scalefactor, -2.5, text("t", 10)), (900 * scalefactor, 0, text("z", 10))], linewidth=1.5)
		plot!(-ts, zs .* cos(ϕ + π/2), color=colors[3])

		# plot!([last(-ts)], [last(zs) * cos(ϕ)], marker=true, color = colors[3])


		plot(p1, p2, p3, layout = l, link=:y, size=(600,500), margins=-2.2mm,legends=:none)
	end
end

# ╔═╡ 6da939ad-8625-43c3-9537-169eab3e97f6
if !dualquad
	
	let


		l = @layout [readout{0.65w} righthist{0.25w} zseries{0.1w}]
		ts = sim.t[1:ti] * 1000
		rs = sim.r[1][1:ti]
		zs = sim.z[1:ti]
		values = range(floor(minimum( sim.r[1])) - 1, floor(maximum( sim.r[1])), step=0.1)
		density = fun(last(zs) * cos(ϕ), 2τ).(values)


		p1 = plot(ts, rs, xlims=[first(sim.t), sim.t[tf] * 1000], color = colorangle(ϕ), xlabel="t (ns)", ylabel="arbitrary units", left_margin=1mm)
		plot!([last(ts)], [last(rs)], marker=true, color = colorangle(ϕ))

		p2 = plot(density, values, frame = :none, color = colorangle(ϕ), linewidth=2)
		ind = getclosest(values, last(rs))
		indz = getclosest(values, last(zs))
		plot!([density[ind]], [last(rs)], marker=true, color = colorangle(ϕ))
		plot!([density[indz]], [last(zs) * cos(ϕ)], marker=true, color = colors[3])


		scalefactor = last(ts)/(sim.t[tf] * 1000)
		p4 = plot(-ts, zs, color=:black, linewidth=0.5)
		plot!(-ts, zs .* cos(ϕ), color=colors[3], frame=:true, ymirror=:true, framestyle = :origin, yticks=[-1.0, 0, 1.0], xticks=:none, showaxis=:hide, right_margin=5mm, annotations=[(-450 * scalefactor, -2.5, text("t", 10)), (900 * scalefactor, 0, text("z", 10))], linewidth=1.5)



		plot(p1, p2, p4, layout = l, link=:y, size=(600,500), margins=-2.0mm,legends=:none)
		
	end
end

# ╔═╡ d67dc66f-93cc-4c36-a305-01ee7a6e897f
let
	
	# Parameters -------------------------------------------------------------------
	
	ρ0 = dm((ground + excited)/√2) # initial state
	dt = 40e-3  # integration time-step
	
	ΩR  = 0 # Rabi frequency (rad * MHz)
	Γ = 2 # Measurement dephasing rate (MHz)
	τ = 1/(2Γ) #  Measurement collapse timescale (μs)
	η = 1 # collection efficiency
	global ϕpm = 0
	
	
	# Kraus operators --------------------------------------------------------------
	
	H = ΩR * σy/2
	J = [(σz, ((1-η) * Γ))]
	C = [(exp(im * ϕpm) * σz, τ, η)]
	
	Random.seed!(1)
	sol = bayesian((0, 3), ρ0, H, J, C; dt=dt)
	
	global Rpm = traj(sol)
	

	
	md" ###### 🔻 Bayesian simulation (single-quadrature)"
end

# ╔═╡ 700cde86-804b-409a-9c2f-03eff2163be8
plot_phases3(Rpm, ϕpm, 1)

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

# ╔═╡ 7a785533-08b8-48ce-94f1-ec219ac64196
if dualquad
	blue(md"""When measuring two quadratures, notice that the Rabi flopping rate is decreased compared to single-quadrature informational measurement.""",title="Dual quadrature measurement")
	
else
	blue(md"""When measuring one quadrature, notice that Rabi flopping exists for informational readout ($\phi = 0$) and goes away for non-informational readout ($\phi = \pi/2$).""",title="Single quadrature measurement")

end

# ╔═╡ Cell order:
# ╠═3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─377a3336-20bd-4baa-a033-af8bbc8668a8
# ╟─3be820e3-8541-424c-b042-facf686f644b
# ╟─66d52a3a-edfd-4bdd-8c18-9f1a27041876
# ╠═e7b46d6f-2f2e-4b60-8684-ee65adcf1431
# ╟─282851f0-fa1d-42b1-a243-5620339129cc
# ╟─259978cd-7888-4f38-bf83-e40719d9777b
# ╟─8d749aee-829a-479c-ba65-b66b2d18bee3
# ╟─4e50b40b-82fc-4eac-a96c-ac10546eecfe
# ╟─9b91683a-e5a4-41e9-81fb-b84cb781cd36
# ╟─ad636be5-a60b-4f36-9470-e33d5d3c891a
# ╟─2662b242-1d39-4223-a600-8b5a5d6a148b
# ╟─db19a712-aa39-4691-93dc-bd9d18f5db7e
# ╟─cf5721a2-6f73-47e9-9b11-1c9562c29a03
# ╟─d05919e8-48c5-40dd-95fe-151d793049d7
# ╟─fac6aa1f-28ec-4e67-b513-6ba485e88eee
# ╟─dbf9e12b-06b0-42cb-968b-d49f7d1c4e10
# ╟─7a785533-08b8-48ce-94f1-ec219ac64196
# ╠═0396d908-8123-49db-a1d6-a637651f81b5
# ╟─1105dac1-61ff-44e9-a82a-fc9823ea7757
# ╟─6145101c-9ce3-4938-bf10-83fa14ff84db
# ╟─49debdd0-b3db-41c3-b6e0-a42750787fed
# ╟─6843e772-74b5-48e1-95c8-5648a93f0787
# ╟─e2443859-dda5-48b4-bf3f-312deeff7a8e
# ╟─082f9721-a5eb-4b78-93cc-a1dc9cc20ad4
# ╟─f21f2ca8-83a9-4b78-8053-f6cfc9536338
# ╟─64fef35f-b894-48db-9ad3-9939c57fcc66
# ╠═6da939ad-8625-43c3-9537-169eab3e97f6
# ╟─3c55abb4-5114-4183-ba06-d12e2f4d8a01
# ╟─a968952f-a64a-4882-a151-8a2d9a1723c7
# ╟─d0957eb9-0df7-4015-9720-5fdb4d1df3b3
# ╟─37dac3f2-f394-467c-ba4f-cf3f9d6d9577
# ╟─553dec63-d7c0-4105-b5dc-f2c36b11eafe
# ╟─32ee829c-6a9c-4ae8-a8ed-8d3516d1afec
# ╟─5b40fc7a-4280-49f2-b2f3-000e9cacc947
# ╟─34a2807c-bd8c-4563-a616-10e0ca3d2e52
# ╟─9a9c5acb-6db5-4a71-8174-75a64d7cb2d0
# ╟─dbf58f2c-9926-4a57-a0f4-0c8a1332df80
# ╟─d67dc66f-93cc-4c36-a305-01ee7a6e897f
# ╟─700cde86-804b-409a-9c2f-03eff2163be8
# ╟─b5c18b03-2543-4d47-a488-93b2c43663fd
# ╟─7e444185-4e30-412a-9a38-cadb547e8372
# ╠═b4068526-178a-42d9-994d-e4e3c879bb6f
# ╟─133b6939-10b2-4c8e-acf8-5658ca96a0f9
# ╠═c1d743fb-356c-4d7a-a290-110c658e20dd
# ╠═a2be6ad8-1ed5-4869-9741-106a348a9d82
# ╟─c63be186-2510-448b-a1e6-fb4c9b96b28f
# ╠═a427e58a-fa75-4d50-acab-da11e49dd1bf
# ╟─57e71783-8e44-453e-a1d9-42377ffe41e9
# ╠═5cef3447-784d-4be1-8e0e-44737997b4ea
# ╟─2546f296-ef35-4702-9dd1-13c1b308529a
# ╟─2d2e530d-5e8a-4045-93bd-97ff99aecaa1
# ╟─3f0bf48d-a27d-492c-b2d7-2291222c4607
# ╟─4c45fe3e-cd69-11eb-20cd-7bfb98c040cf
# ╠═e43e5329-bd96-41ce-a183-1bd206204f65
# ╟─ea04744b-4296-4dc2-9a3c-1f477c96f1ac
# ╠═f3ac9b41-b123-4236-8c62-8e2b988463c2
# ╟─8b4bc294-989b-4f1c-82b9-e418fdd3c40b
# ╟─235cee23-c8af-4de2-a1c3-a2b173156703
# ╟─e27bd39c-58b7-4c5c-a677-3fe70f500ee8
# ╠═28d8df47-07cf-4a0d-a447-f894d021b2bc
# ╟─dc3e2dcb-5bf1-492d-8337-f366dcf0170b
# ╠═460379bc-0d49-434a-b0bc-3efcfdc47b5c
# ╠═01523c93-5737-4d94-87fc-2e5fb731002c
# ╠═de310c78-ae02-488f-a939-e4d29faa3651
# ╠═f8cfd829-06ef-4971-b216-5a6abe1b072d
# ╠═87e2a8c9-75b6-486e-95b3-2506dd255992
