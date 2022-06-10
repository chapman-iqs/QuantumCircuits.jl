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

# ╔═╡ e49291f6-9ab8-43f7-af88-95b6c9c2181c
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
	using Plots.Measures
	
	include("utilities/single-qubit-operators.jl")
	include("utilities/utilities.jl")
	include("utilities/plotting.jl")

	include("notebooks/table-of-contents.jl")
	include("resources.jl")

	md" # Packages and julia files"
	
	
end

# ╔═╡ 377a3336-20bd-4baa-a033-af8bbc8668a8
md"""
In this interactive notebook, we'll look at the Quantum Zeno effect to help gain intuition for measurement and readout.
"""

# ╔═╡ f0b0cb6e-4112-481c-89bb-35d0d6ad1d96
mdp(table_of_contents📔)

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents(title="Quantum Zeno effect")

# ╔═╡ 3be820e3-8541-424c-b042-facf686f644b
md"""
# Quantum Zeno effect
"""

# ╔═╡ 66d52a3a-edfd-4bdd-8c18-9f1a27041876
md"""
To better observe the effects of continuous measurement, we'll look at the so-called Quantum Zeno effect ($Itano_et_al), described in terms of a superconducting qubit experiment in $(Slichter_et_al).
"""

# ╔═╡ e7b46d6f-2f2e-4b60-8684-ee65adcf1431
begin
	ψ0 = normalize(g + e) # initial state
	dt = 1e-3  # integration time-step
	tf = 5.0

	ΩR  = 2π * 3.6 # Rabi frequency (rad * MHz)
	Γ = 134. #0.15 # Measurement dephasing rate (MHz)
	η = 1 # collection efficiency
	τ = 1/(2Γ*η) #  Measurement collapse timescale (μs)
end

# ╔═╡ bc4f647d-037f-41c4-b31e-4b21bf883236
md"""
The qubit is driven at a significant Rabi rate ($\Omega_R$ = $(round(ΩR/(2π), digits=3)) MHz), 

while also being measured strongly at rate $\Gamma$ = $Γ MHz. This imitates parameters used in the experiment in $Slichter_et_al. 

Note that the strong measurement rate chosen here violates assumptions of the `bayesian` and `rouchon` solvers, namely that the readout can be modeled as sampling from a single Gaussian with shifted mean. However, the solvers still reproduce the important qualitative features of the effect (particularly since the state primarily is constrained to the two eigenstates), so we show it here as a demonstration to give intuition about measurement.
"""

# ╔═╡ ff0ac02b-a19e-434a-bb2e-bcc30c415efb
md"""


The corresponding measurement collapse timescale -- the time the state takes to collapse to a measurement eigenstate -- is $\tau = 1/2\Gamma \eta$ = $(round(τ*1000, digits=2)) ns. Thus, the state will generally be pinned to an eigenstate; however, the Rabi drive allows for quick transitions or "jumps" to the opposite eigenstate.
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
When measuring two quadratures, the efficiency $\eta$ on each is halved. The effective measurement time $\tau = 1/(2\Gamma\eta)$ is doubled; however, the measurement rate $\Gamma$ input to the solvers remains the same, since this corresponds to an ensemble-averaged measurement dephasing rate.
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

`tf = 0.05 μs`
$(@bind tfactor html"<input type=range min=0.05 max=5 step=0.05 value=1>") 
`tf = 5 μs`

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
When we do readout, we are sampling from a distribution centered around the expected value of $z$. (Remember that this is an approximation in the limit of weak measurement, as discussed on p. 4-5 of $Jacobs_Steck_2006.) We then use this information to update our estimate of the quantum state.
"""

# ╔═╡ 6145101c-9ce3-4938-bf10-83fa14ff84db
md"""
The records plotted above are sampled in this way. Thus, as $z$ evolves, the mean of the distribution shifts along with the quantum state. Depending on the strength of measurement, this shift can be almost imperceptible compared to the variance. However, this small shift accounts for state-dependent backaction.

In this toy example, we have chosen large $\Gamma$ to exagerrate this shift so that it is perceptible.
"""

# ╔═╡ 49debdd0-b3db-41c3-b6e0-a42750787fed
md"""
You can play with the plot below to understand the relationship between the semi-random readout and the quantum state.
"""

# ╔═╡ 6843e772-74b5-48e1-95c8-5648a93f0787
begin
	fnorm(μ, σ) = x ->  1/(σ * √(2π)) * exp(-0.5 * (x - μ)^2  / σ^2)
	fun(z, τ) = fnorm(z, √(τ/dt))
	md" 🌀 Define normal distribution functions"
end

# ╔═╡ f21f2ca8-83a9-4b78-8053-f6cfc9536338
md"""

`t = 0`
$(@bind ti html"<input type=range min=5 max=2000 step=5 value=5>") 
`t = 2 μs`

"""

# ╔═╡ 3c55abb4-5114-4183-ba06-d12e2f4d8a01
md"""
The first time series is the readout -- what is observed in the lab, up to shifting and rescaling. The Gaussian curve is the distribution the data are sampled from, which moves with the expectation value of $z$. The green time series represents the readout's tracking of the $z$ coordinate, compared to the black series, which is the true $z$ value.
"""

# ╔═╡ f0f48d8f-997d-4021-8923-1e2bc954a539
mdp("This demonstrates how the quantum state evolution affects the sampling of the measurement record. However, the connection goes both ways: the measurement itself affects quantum state evolution. To see this in action, take a look at ",  measurement_backaction📔, ".")

# ╔═╡ 133b6939-10b2-4c8e-acf8-5658ca96a0f9
md" # Utilities"

# ╔═╡ ea04744b-4296-4dc2-9a3c-1f477c96f1ac
md"""
## Plotting
"""

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
	@userplot BlochReadout
	
	@recipe function f(bts::BlochReadout; vec=nothing, tf=nothing)
		ϕ, sim = bts.args
		ts, xs, ys, zs, ps, r = sim.t, sim.x, sim.y, sim.z, sim.p, sim.r
	
		# Plot time series --------------------------------------------------------
		
		legend := [:none :none :outerright]
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

# ╔═╡ 3f0bf48d-a27d-492c-b2d7-2291222c4607
colorangle(ϕ) = RGB(200/255,abs(sin(ϕ))*200/255,200/255)

# ╔═╡ a427e58a-fa75-4d50-acab-da11e49dd1bf
function plot_phases(sim, ϕ, tfactor)
	
	index = Int64(floor(length(sim.t) * tfactor / tf))
	
	
	br = blochreadout(ϕ, sim, tf=sim.t[index], ylims=[-1.1,1.5])
	
	plot!([(real(phases), imag(phases)), ([0, cos(ϕ)], [0, sin(ϕ)])],  inset = bbox(0, 0.1, 0.1, 0.1, :right), subplot=4, legend = :none, frame = :none, aspect_ratio = :equal, color = [:black colorangle(ϕ)], linewidth=[1 2])
	
	
	plot!([(real(phases), imag(phases)), ([0, cos(ϕ + π/2)], [0, sin(ϕ + π/2)])],  inset = bbox(0, 0.38, 0.1, 0.1, :right), subplot=5, legend = :none, frame = :none, aspect_ratio = :equal, color = [:black colorangle(ϕ + π/2)], linewidth=[1 2])
		

	
end
	

# ╔═╡ 8b4bc294-989b-4f1c-82b9-e418fdd3c40b
md"""
## Other
"""

# ╔═╡ c1d743fb-356c-4d7a-a290-110c658e20dd
getclosest(array, val) = argmin(abs.(val .- array))

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
	
	function traj(sol::Solution)
	  t, ρ, r = (sol.t, sol.ρ, sol.r)
	  x, y, z = map(op -> expectations(sol, op), [σx, σy, σz])
	  p = (typeof(ρ[1]) <: Ket) ? [1.0 for el in ρ] : real(expect.(ρ, ρ))
	  traj(t, x, y, z, p, r)
	end
end

# ╔═╡ 4e50b40b-82fc-4eac-a96c-ac10546eecfe
let	
	# Kraus operators --------------------------------------------------------------
	
	H = ΩR * σy/2
	J = [(σz, ((1-η) * Γ))]
	C = [(exp(im * ϕ) * σz, Γ, η)]
	
	sol = bayesian((0, tf), ψ0, H, J, C; dt=dt)
	
	global Rsq = traj(sol)
	
	md" ###### 🔻 Bayesian simulation (single-quadrature)"
end

# ╔═╡ 9b91683a-e5a4-41e9-81fb-b84cb781cd36
let
	# Kraus operators --------------------------------------------------------------
	
	H = ΩR * σy/2
	J = [(σz, ((1-η) * Γ))]
	C = [(exp(im * ϕ) * σz, Γ, η/2), (exp(im * (ϕ + π/2)) * σz, Γ, η/2)]
	
	Random.seed!(1)
	sol = bayesian((0, tf), ψ0, H, J, C; dt=dt)
	
	global Rdq = traj(sol)
	

	
	md" ###### 🔻 Bayesian simulation (dual-quadrature)"
end

# ╔═╡ fac6aa1f-28ec-4e67-b513-6ba485e88eee
begin
	sim = dualquad ? Rdq : Rsq
	plot_phases(sim, ϕ, tfactor)
end

# ╔═╡ 082f9721-a5eb-4b78-93cc-a1dc9cc20ad4
md" `t = ` $(sim.t[ti]) μs"

# ╔═╡ 3670a0e6-7a61-4c89-87a3-a98792c40950
tff = getclosest(sim.t, 1)

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

		# xlims=[first(sim.t), sim.t[tf] * 1000]


		p1 = plot(ts, rs1, xlims=[first(sim.t), sim.t[tff] * 2000], color = colorangle(ϕ), xlabel="t (ns)", ylabel="arbitrary units", left_margin=1mm)
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


		scalefactor = last(ts)/(sim.t[tff] * 2000)
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


		p1 = plot(ts, rs, xlims=[first(sim.t), sim.t[tff] * 2000], color = colorangle(ϕ), xlabel="t (ns)", ylabel="arbitrary units", left_margin=1mm)
		plot!([last(ts)], [last(rs)], marker=true, color = colorangle(ϕ))

		p2 = plot(density, values, frame = :none, color = colorangle(ϕ), linewidth=2)
		ind = getclosest(values, last(rs))
		indz = getclosest(values, last(zs))
		plot!([density[ind]], [last(rs)], marker=true, color = colorangle(ϕ))
		plot!([density[indz]], [last(zs) * cos(ϕ)], marker=true, color = colors[3])


		scalefactor = last(ts)/(sim.t[tff] * 2000)
		p4 = plot(-ts, zs, color=:black, linewidth=0.5)
		plot!(-ts, zs .* cos(ϕ), color=colors[3], frame=:true, ymirror=:true, framestyle = :origin, yticks=[-1.0, 0, 1.0], xticks=:none, showaxis=:hide, right_margin=5mm, annotations=[(-450 * scalefactor, -2.5, text("t", 10)), (900 * scalefactor, 0, text("z", 10))], linewidth=1.5)



		plot(p1, p2, p4, layout = l, link=:y, size=(600,500), margins=-2.0mm,legends=:none)
		
	end
end

# ╔═╡ dc3e2dcb-5bf1-492d-8337-f366dcf0170b
φdict = Dict("0" => 0, 
				"π/8" => π/8, 
				"π/4" => π/4,
				"3π/8" => 3π/8,
				"π/2" => π/2)

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
	blue(md"""When measuring two quadratures, notice that the transition rate $\ket{+z} \leftrightarrow \ket{-z}$ is decreased compared to single-quadrature informational measurement.""",title="Dual quadrature measurement")
	
else
	blue(md"""When measuring one quadrature, notice that transitions $\ket{+z} \leftrightarrow \ket{-z}$ exist for informational readout ($\phi = 0$) and go away for non-informational readout ($\phi = \pi/2$).""",title="Single quadrature measurement")

end

# ╔═╡ Cell order:
# ╟─377a3336-20bd-4baa-a033-af8bbc8668a8
# ╟─f0b0cb6e-4112-481c-89bb-35d0d6ad1d96
# ╟─3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─3be820e3-8541-424c-b042-facf686f644b
# ╟─66d52a3a-edfd-4bdd-8c18-9f1a27041876
# ╟─bc4f647d-037f-41c4-b31e-4b21bf883236
# ╟─ff0ac02b-a19e-434a-bb2e-bcc30c415efb
# ╠═e7b46d6f-2f2e-4b60-8684-ee65adcf1431
# ╟─282851f0-fa1d-42b1-a243-5620339129cc
# ╟─259978cd-7888-4f38-bf83-e40719d9777b
# ╟─8d749aee-829a-479c-ba65-b66b2d18bee3
# ╠═4e50b40b-82fc-4eac-a96c-ac10546eecfe
# ╟─9b91683a-e5a4-41e9-81fb-b84cb781cd36
# ╟─ad636be5-a60b-4f36-9470-e33d5d3c891a
# ╟─2662b242-1d39-4223-a600-8b5a5d6a148b
# ╟─db19a712-aa39-4691-93dc-bd9d18f5db7e
# ╠═cf5721a2-6f73-47e9-9b11-1c9562c29a03
# ╟─d05919e8-48c5-40dd-95fe-151d793049d7
# ╟─fac6aa1f-28ec-4e67-b513-6ba485e88eee
# ╟─dbf9e12b-06b0-42cb-968b-d49f7d1c4e10
# ╟─7a785533-08b8-48ce-94f1-ec219ac64196
# ╟─0396d908-8123-49db-a1d6-a637651f81b5
# ╟─1105dac1-61ff-44e9-a82a-fc9823ea7757
# ╟─6145101c-9ce3-4938-bf10-83fa14ff84db
# ╟─49debdd0-b3db-41c3-b6e0-a42750787fed
# ╟─6843e772-74b5-48e1-95c8-5648a93f0787
# ╟─082f9721-a5eb-4b78-93cc-a1dc9cc20ad4
# ╟─f21f2ca8-83a9-4b78-8053-f6cfc9536338
# ╟─3670a0e6-7a61-4c89-87a3-a98792c40950
# ╟─64fef35f-b894-48db-9ad3-9939c57fcc66
# ╟─6da939ad-8625-43c3-9537-169eab3e97f6
# ╟─3c55abb4-5114-4183-ba06-d12e2f4d8a01
# ╟─f0f48d8f-997d-4021-8923-1e2bc954a539
# ╟─133b6939-10b2-4c8e-acf8-5658ca96a0f9
# ╟─ea04744b-4296-4dc2-9a3c-1f477c96f1ac
# ╠═a2be6ad8-1ed5-4869-9741-106a348a9d82
# ╟─c63be186-2510-448b-a1e6-fb4c9b96b28f
# ╟─a427e58a-fa75-4d50-acab-da11e49dd1bf
# ╟─5cef3447-784d-4be1-8e0e-44737997b4ea
# ╟─3f0bf48d-a27d-492c-b2d7-2291222c4607
# ╟─8b4bc294-989b-4f1c-82b9-e418fdd3c40b
# ╟─c1d743fb-356c-4d7a-a290-110c658e20dd
# ╟─28d8df47-07cf-4a0d-a447-f894d021b2bc
# ╟─dc3e2dcb-5bf1-492d-8337-f366dcf0170b
# ╟─01523c93-5737-4d94-87fc-2e5fb731002c
# ╠═de310c78-ae02-488f-a939-e4d29faa3651
# ╟─f8cfd829-06ef-4971-b216-5a6abe1b072d
# ╟─87e2a8c9-75b6-486e-95b3-2506dd255992
# ╠═e49291f6-9ab8-43f7-af88-95b6c9c2181c
