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
# Excitation feedback

In this interactive notebook, we'll explore controlling onto an excitation state using measurement feedback.
"""

# ╔═╡ 3be820e3-8541-424c-b042-facf686f644b
md"""
# Two qubit example
"""

# ╔═╡ 009c334a-74cc-466c-9eb8-0f8aebfc9164
md"""
## System description
"""

# ╔═╡ 66d52a3a-edfd-4bdd-8c18-9f1a27041876
md"""
We'll start from the two qubit example, where we'll define 

$\ket{\overline{0}} \equiv \ket{00}$

$\ket{\overline{1}} \equiv \frac{\ket{01} + \ket{10}}{\sqrt{2}}$

$\ket{\overline{2}} \equiv \ket{11}.$

We can count the number of excitations using the total excitation operator 

$\hat n = \hat n_1 + \hat n_2 = (\hat \sigma_+ \hat \sigma_-)_1 \otimes I_2 + I_1 \otimes (\hat \sigma_+ \hat \sigma_-)_2.$
"""

# ╔═╡ deb2a0f3-32d1-46f9-bfbe-2e197ead243d
md"""
We weakly measure the qubits collectively by measuring the number of dispersive shifts of the resonator. This corresponds to a measurement of $\hat n$. The qubits will be effectively identical in all respects, including resonant frequency and coupling to the cavity. We assume the qubits do not interact with one another.

The measurement of $\hat n$ results in a noisy measurement record

$\tilde r(t) = \braket{\hat n} + \zeta(t), \hspace{5mm} \zeta(t) \sim \mathcal N(0, \tau/dt) \label{blash}$

where $\tau$ is the measurement collapse timescale. **I'm guessing this is the same collectively as it is for an individual qubit, but not sure.**
"""

# ╔═╡ 58004c67-ca73-4095-b38f-e6658f1f998f
md"""
Each qubit receives an identical time-dependent Rabi drive modulated by the stochastic readout, offset by $r_{target}$, and amplified with some gain $\Omega_R$. Thus the full Hamiltonian in the rotating frame is

$H = \Omega_R (\tilde{r}(t) - r_{target}) (\hat \sigma_1^x \otimes I_2 + I_1 \otimes \hat \sigma_2^x).$
"""

# ╔═╡ 5a01f315-4db7-4dac-8b98-891b5b8eef22
md"""
## Feedback protocol
"""

# ╔═╡ f3cc4e10-afad-4f28-8915-8132de900b4f
md"""
Assuming homodyne measurement, there will be three DC voltages corresponding to the number eigenstates $\ket{\overline 0}$, $\ket{\overline 1}$, and $\ket {\overline 2}$ which we will call $V_0$, $V_1$, and $V_2$. In simulation, we model the measurement record as a shifted and stretched version of the experimental voltage trace, such that it has the form given above as $\tilde{r}(t)$. Then for each $V_i$ we can associate an $r_i = \langle \overline{i} | \hat n | \overline i \rangle$, i.e. $r_0 = 0$, $r_1 = 1$, $r_2 = 2$. 

To control onto a specific target state, we set $r_{target} = r_i$. The Rabi drive effectively cancels the backaction when the system is in the target state, creating a point of stability there that the system will be drawn towards.
"""

# ╔═╡ e43e5329-bd96-41ce-a183-1bd206204f65
begin
	# Basis
	q = SpinBasis(1//2)
	Iq = identityoperator(q)
	I = Iq ⊗ Iq

	# qubit operators, using convention that |-z> is ground state
	σx = sigmax(q)
	σy = sigmay(q)
	σz = sigmaz(q)
	σp = sigmap(q)
	σm = sigmam(q)
	nq = σp * σm

	g = spindown(q)
	e = spinup(q)

	# two-qubit operators
	σx1 = σx ⊗ Iq
	σx2 = Iq ⊗ σx
	σy1 = σy ⊗ Iq
	σy2 = Iq ⊗ σy
	σz1 = σz ⊗ Iq
	σz2 = Iq ⊗ σz
	σp1 = σp ⊗ Iq
	σm1 = σm ⊗ Iq
	σp2 = Iq ⊗ σp
	σm2 = Iq ⊗ σm

	# number operators
	n1 = nq ⊗ Iq
	n2 = Iq ⊗ nq
	n = n1 + n2


	# projectors
	gg = dm(g ⊗ g)
	ge = dm(g ⊗ e)
	eg = dm(e ⊗ g)
	ee = dm(e ⊗ e)
	
	md"### 🔶 Hilbert space operators "
end

# ╔═╡ 282851f0-fa1d-42b1-a243-5620339129cc
md"""
## Simulation
"""

# ╔═╡ cff9fded-2ecb-48fa-899d-202de36f74a9
begin
	function smooth(fine::Array; n=2)
	    coarse = []
		m = Int64(floor(n/2))
		for i in 1:(m-1)
			push!(coarse, mean(fine[1:i+m])) end
	
		for i in m:(length(fine) - m)
			push!(coarse, mean(fine[i-(m-1):i+m])) end
	
		for i in (length(fine) - m + 1):length(fine)
			push!(coarse, mean(fine[i-(m-1):length(fine)])) end
		
	    coarse
	end
	function smooth(fineSA::SubArray; n=2)
	fine = parent(fineSA)
    coarse = []
	m = Int64(floor(n/2))
	for i in 1:(m-1)
		push!(coarse, mean(fine[1:i+m])) end

	for i in m:(length(fine) - m)
		push!(coarse, mean(fine[i-(m-1):i+m])) end

	for i in (length(fine) - m + 1):length(fine)
		push!(coarse, mean(fine[i-(m-1):length(fine)])) end
	
    coarse
	end
end

# ╔═╡ c15a1717-4c3a-4e49-8523-d6653b6e76b4
md" ` seed = 1` $(@bind s Slider(1:1:100)) `100`"

# ╔═╡ afa6e410-5c7a-48fa-b724-6419d04f39af
md" s = $s "

# ╔═╡ 467cb6d1-cbc2-482e-96fd-5432c7d18c6f
lims = [-0.5,5]

# ╔═╡ 9f076306-27b0-4573-9385-a75964d4a670
N = 300

# ╔═╡ 016972b6-895f-4295-8a85-33377cf158fb
md"""
$ρ_0 \propto \epsilon \ket{11} + \gamma \ket{00} + (\ket{01} + \ket{10})$
"""

# ╔═╡ e64e47a2-410e-43de-bd8f-d9d2ff011bea
)

# ╔═╡ 12d22513-478c-4bb0-8077-513da10004e0
tf = 5

# ╔═╡ f2670e6e-f37e-4eff-ad7e-05181772aec7
2*(n - I)

# ╔═╡ 6bf29e9a-3f01-4e0f-b9c6-1007c8007ffb
md"""
### Single simulation
"""

# ╔═╡ e7b46d6f-2f2e-4b60-8684-ee65adcf1431
begin
	β = 2
	ε, γ = (0, β)
	ρ0 = normalize(dm(ε*(e ⊗ e) + γ*(g ⊗ g) +(e ⊗ g + g ⊗ e)))
	dt = 1e-3  # integration time-step

	Γ = 1.0 # Ensemble measurement dephasing rate (MHz)	
	η = 1.0 # collection efficiency
	ΩR = 0 # Rabi frequency (rad * MHz)

	md" ### 🌀 Parameters"
end

# ╔═╡ 020ce994-b653-4024-9445-851da8af7aca
let
	Random.seed!()
	# Kraus operators --------------------------------------------------------------
	H = 0 * σx1
	J = [(n, ((1 - η) * Γ))]
	C = [(n, Γ, η)]

	ρ00(ε, γ) = normalize(dm(ε*(e ⊗ e) + γ*(g ⊗ g) +(e ⊗ g + g ⊗ e)))

	global nfε = map(1:N) do m
					ε, γ = (β, 0)
					sol = bayesian((0, tf), ρ00(ε, γ), H, J, C; dt=dt)
					ρf = last(sol.ρ)
					return round(real(expect(ρf, n)), digits=3)
				end

	global nfγ = map(1:N) do m
					ε, γ = (0, β)
					sol = bayesian((0, tf), ρ00(ε, γ), H, J, C; dt=dt)
					ρf = last(sol.ρ)
					return round(real(expect(ρf, n)), digits=3)
				end
	

	md" ###### 🔻 Bayesian simulation"
end

# ╔═╡ 9da84b37-932c-4d17-94a5-6288d876f396
nfε[nfε .> 1.9]

# ╔═╡ 1ad5981d-e66f-4b71-b83d-21898723d9c6
let
	binrange = -0.2:0.1:2.2
	histogram(round.(nfε, digits=1), bins=binrange, alpha=0.4, label="(ε, γ) = (β, 0)", title=string("final n for N = ", N, " runs, tf = ", tf, " μs, β = ", β))
	histogram!(round.(nfγ, digits=1), bins=binrange, alpha=0.4, label="(ε, γ) = (0, β)", legend=:left)
	# histogram!(n1f, bins=0:0.05:2.2, alpha=0.4, label="n1", barwidth=0.1)
	# histogram!(n2f, bins=0:0.05:2.2, alpha=0.4, label="n2", barwidth=0.1)
end

# ╔═╡ 4e50b40b-82fc-4eac-a96c-ac10546eecfe
let
	# Kraus operators --------------------------------------------------------------
	δ = 5
	m = n - δ * I
	
	H = ΩR * (σx1 + σx2)
	J = [(m, ((1 - η) * Γ))]
	C = [(m, Γ, η)]
	
	Random.seed!(s)
	global sol = bayesian((0, tf), ρ0, H, J, C; dt=dt)
	
	md" ###### 🔻 Bayesian simulation"
end

# ╔═╡ f6634cd0-990e-4e9f-bd21-f001719c52c4
smoothed = smooth(sol.r[1], n=2000)

# ╔═╡ 2d3627cf-1a6e-4461-a747-b8d5b8c1ef2a
begin
	ns, n1s, n2s, x1s, x2s, y1s, y2s, z1s, z2s = [map(ρ -> real(expect(ρ, op)), sol.ρ) for op in [n, n1, n2, σx1, σx2, σy1, σy2, σz1, σz2]]

	ggs, ges, egs, ees = [map(ρ -> real(expect(ρ, op)), sol.ρ) for op in [gg, ge, eg, ee]]

	p1s = 0.5 * (1 .+ x1s.^2 .+ y1s.^2 .+ z1s.^2)
	p2s = 0.5 * (1 .+ x2s.^2 .+ y2s.^2 .+ z2s.^2)

	ps = map(ρ -> real(expect(ρ, ρ)), sol.ρ)

end

# ╔═╡ 133b6939-10b2-4c8e-acf8-5658ca96a0f9
md" # Utilities"

# ╔═╡ c1d743fb-356c-4d7a-a290-110c658e20dd
getclosest(array, val) = argmin(abs.(val .- array))

# ╔═╡ a2be6ad8-1ed5-4869-9741-106a348a9d82
colors = palette(:tab10)

# ╔═╡ d71e8206-8cb0-4399-8379-aba003e32ea2
colors2 = palette(:tol_muted)

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

# ╔═╡ 7cb61b30-803b-4412-ad5a-65e3e10bda0a
let
	# l = @layout [readout{0.2h}; nexp{0.2h}; basisstates{0.2h}; bloch1{0.2h}; bloch2{0.2h}]
	l = @layout [readout{0.3h}; nexp{0.15h}; basisstates{0.25h}; bloch1{0.15h}; bloch2{0.15h}]

	p1 = plot(sol.t, smoothed, color = :black, ylims=lims, legend=:false, title=string("readout; seed = ", s))
	plot!(twinx(), sol.t, sol.r, color = colorangle(0), left_margin=1mm, foreground_color_axis=colorangle(0), foreground_color_border=colorangle(0), legend=:false)
	plot!(twinx(), sol.t, smoothed, color = :black, ylims=lims, frame=:none, legend=:false)
	plot!(twinx(), [0, last(sol.t)], [1,1], color = :black, linestyle=:dash, ylims=lims, frame=:none, legend=:false, linewidth=2)
	plot!(twinx(), [0, last(sol.t)], [0,0], color = :black, linestyle=:dash, ylims=lims, frame=:none, legend=:false, linewidth=2)
	plot!(twinx(), [0, last(sol.t)], [2,2], color = :black, linestyle=:dash, ylims=lims, frame=:none, legend=:false, linewidth=2)

	p2 = plot(sol.t, ns, ylabel=L"\langle n \rangle", color=colors2[8], title="excitation number", legend=:false, ylims=[0,2.2])

	p3 = plot(sol.t, ggs, color=colors2[1], label=L"|00 \rangle")
	plot!(sol.t, ges, color=colors2[2], label=L"|01 \rangle")
	plot!(sol.t, egs, color=colors2[3], label=L"|10 \rangle")
	plot!(sol.t, ees, color=colors2[4], label=L"|11 \rangle")
	plot!(sol.t, ps, color=colors[4], label=L"Tr(\rho^2)", title="2-qubit states")

	p4 = plot(sol.t, x1s, color=colors[1], label=L"x")
	plot!(sol.t, y1s, color=colors[2], label=L"y")
	plot!(sol.t, z1s, color=colors[3], label=L"z")
	plot!(sol.t, p1s, color=colors[4], label=L"Tr(\rho_1^2)", title="qubit 1")

	p5 = plot(sol.t, x2s, color=colors[1], label=L"x")
	plot!(sol.t, y2s, color=colors[2], label=L"y")
	plot!(sol.t, z2s, color=colors[3], label=L"z")
	plot!(sol.t, p2s, color=colors[4], label=L"Tr(\rho_2^2)", title="qubit 2")

	plot(p1, p2, p3, p4, p5,  layout = l, link=:y, size=(600,800), legendfontsize=10, titlefontsize=12)
	
end

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

# ╔═╡ 87f29d4c-63d5-4d2f-9126-90bbd67fe2bd
f(b::Vector) = "this is a vector"

# ╔═╡ 2b310f40-0132-41b5-b3c2-9c334bf5b63e
f(k::Operator) = "this is an operator"

# ╔═╡ bcbb2d05-5278-49a0-a7f5-2624f04b782f
f(ρ0)

# ╔═╡ 23dfa7e1-a9d5-49e4-9674-21300fb34b24
f(sol.r[1])

# ╔═╡ 6df040fd-3585-4839-a63e-5e057f6818df
typeof(sol.r[1]) <: Vector

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

# ╔═╡ 614813fb-069b-404d-9ec5-7e7e4a8fe6a3
md"""
### Equation numbering
"""

# ╔═╡ Cell order:
# ╠═3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─377a3336-20bd-4baa-a033-af8bbc8668a8
# ╟─3be820e3-8541-424c-b042-facf686f644b
# ╟─009c334a-74cc-466c-9eb8-0f8aebfc9164
# ╟─66d52a3a-edfd-4bdd-8c18-9f1a27041876
# ╟─deb2a0f3-32d1-46f9-bfbe-2e197ead243d
# ╟─58004c67-ca73-4095-b38f-e6658f1f998f
# ╟─5a01f315-4db7-4dac-8b98-891b5b8eef22
# ╟─f3cc4e10-afad-4f28-8915-8132de900b4f
# ╟─e43e5329-bd96-41ce-a183-1bd206204f65
# ╟─282851f0-fa1d-42b1-a243-5620339129cc
# ╟─cff9fded-2ecb-48fa-899d-202de36f74a9
# ╠═f6634cd0-990e-4e9f-bd21-f001719c52c4
# ╟─c15a1717-4c3a-4e49-8523-d6653b6e76b4
# ╟─afa6e410-5c7a-48fa-b724-6419d04f39af
# ╠═467cb6d1-cbc2-482e-96fd-5432c7d18c6f
# ╟─7cb61b30-803b-4412-ad5a-65e3e10bda0a
# ╠═9f076306-27b0-4573-9385-a75964d4a670
# ╠═020ce994-b653-4024-9445-851da8af7aca
# ╟─016972b6-895f-4295-8a85-33377cf158fb
# ╠═9da84b37-932c-4d17-94a5-6288d876f396
# ╠═e64e47a2-410e-43de-bd8f-d9d2ff011bea
# ╠═1ad5981d-e66f-4b71-b83d-21898723d9c6
# ╠═12d22513-478c-4bb0-8077-513da10004e0
# ╠═f2670e6e-f37e-4eff-ad7e-05181772aec7
# ╟─6bf29e9a-3f01-4e0f-b9c6-1007c8007ffb
# ╠═4e50b40b-82fc-4eac-a96c-ac10546eecfe
# ╠═e7b46d6f-2f2e-4b60-8684-ee65adcf1431
# ╠═2d3627cf-1a6e-4461-a747-b8d5b8c1ef2a
# ╟─133b6939-10b2-4c8e-acf8-5658ca96a0f9
# ╠═c1d743fb-356c-4d7a-a290-110c658e20dd
# ╠═a2be6ad8-1ed5-4869-9741-106a348a9d82
# ╠═d71e8206-8cb0-4399-8379-aba003e32ea2
# ╟─c63be186-2510-448b-a1e6-fb4c9b96b28f
# ╟─a427e58a-fa75-4d50-acab-da11e49dd1bf
# ╟─57e71783-8e44-453e-a1d9-42377ffe41e9
# ╟─5cef3447-784d-4be1-8e0e-44737997b4ea
# ╟─2546f296-ef35-4702-9dd1-13c1b308529a
# ╟─2d2e530d-5e8a-4045-93bd-97ff99aecaa1
# ╟─3f0bf48d-a27d-492c-b2d7-2291222c4607
# ╟─4c45fe3e-cd69-11eb-20cd-7bfb98c040cf
# ╟─ea04744b-4296-4dc2-9a3c-1f477c96f1ac
# ╠═f3ac9b41-b123-4236-8c62-8e2b988463c2
# ╟─8b4bc294-989b-4f1c-82b9-e418fdd3c40b
# ╠═87f29d4c-63d5-4d2f-9126-90bbd67fe2bd
# ╠═2b310f40-0132-41b5-b3c2-9c334bf5b63e
# ╠═bcbb2d05-5278-49a0-a7f5-2624f04b782f
# ╠═23dfa7e1-a9d5-49e4-9674-21300fb34b24
# ╠═6df040fd-3585-4839-a63e-5e057f6818df
# ╟─235cee23-c8af-4de2-a1c3-a2b173156703
# ╟─e27bd39c-58b7-4c5c-a677-3fe70f500ee8
# ╠═28d8df47-07cf-4a0d-a447-f894d021b2bc
# ╟─dc3e2dcb-5bf1-492d-8337-f366dcf0170b
# ╠═460379bc-0d49-434a-b0bc-3efcfdc47b5c
# ╠═01523c93-5737-4d94-87fc-2e5fb731002c
# ╠═de310c78-ae02-488f-a939-e4d29faa3651
# ╠═f8cfd829-06ef-4971-b216-5a6abe1b072d
# ╠═87e2a8c9-75b6-486e-95b3-2506dd255992
# ╟─614813fb-069b-404d-9ec5-7e7e4a8fe6a3
