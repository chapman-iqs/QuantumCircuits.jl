### A Pluto.jl notebook ###
# v0.16.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
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
	
	include("plotting.jl")
	
	md" ###### 🔶 Packages and julia files"
end

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents()

# ╔═╡ 377a3336-20bd-4baa-a033-af8bbc8668a8
md"""
# Bloch sphere

In this interactive notebook, we'll look at quantum trajectories on the Bloch sphere.
"""

# ╔═╡ e43e5329-bd96-41ce-a183-1bd206204f65
begin
	# Basis
	q = SpinBasis(1//2)

	# Operators, using convention that |-z> is ground state
	σx = sigmax(q)
	σy = sigmay(q)
	σz = sigmaz(q)
	σp = sigmap(q)
	σm = sigmam(q)
	I = identityoperator(q)
	
	ground = spindown(q)
	excited = spinup(q)

	md"###### 🔶 Qubit Hilbert space operators "
end

# ╔═╡ 02d2ca0b-20bf-4dc2-a8f8-e72b5ae40ed2
md"""

Change measurement angle: `ϕ = 0 `
$(@bind ϕc html"<input type=range min=0 max=16 step=1 value=0>") 
`ϕ = 2π`

"""

# ╔═╡ f3e7794e-a9a1-4103-8469-32a8e37d2d82
md"""

Rotate view: `ϕv = 0`
$(@bind ϕvc html"<input type=range min=0 max=32 step=1 value=0>") 
`ϕv = 2π`

"""

# ╔═╡ 01b79d3d-a6ba-4523-a668-b78a455279cb
begin
	ϕ = ϕc * (π/8)
	ϕv = ϕvc * (π/16)
	md" `ϕ =` $(ϕc/8) π,  `ϕv =` $(ϕvc/16) π"
end

# ╔═╡ a12cdb8c-e9a1-4c2d-9811-cff266e152d8
md" $(@bind show_gif CheckBox()) Animate "

# ╔═╡ 0a7f28c9-1d84-43e6-b62c-711a231a3972
md" $(@bind show_gif2 CheckBox()) Animate Bloch series "

# ╔═╡ d9f2f00f-4ee2-45b5-91d6-6552d6d5b6c1
md" $(@bind show_gif3 CheckBox()) Animate cross-sections "

# ╔═╡ 133b6939-10b2-4c8e-acf8-5658ca96a0f9
md" # Utilities"

# ╔═╡ 235cee23-c8af-4de2-a1c3-a2b173156703
expects(ops) = ρ -> collect(real(expect(ρ, s)) for s in vcat(ops, ρ)) # ρ -> [<x>,<y>,<z>,<ρ>]

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
	
	function traj(sol::QuantumCircuits.solution)
		t, ρ, r = (sol.t, sol.ρ, sol.r)
		x, y, z = [real(expect(σi, ρ)) for σi in (σx, σy, σz)]
		p = real(expect.(ρ, ρ))
		traj(t, x, y, z, p, r)
	end
	
end

# ╔═╡ ecc19d12-e53f-4904-b13f-7a2ccc4912d7
let
	# simulation parameters
	ρ0 = dm(ground) # initial state
	dt = 1e-3  # integration time-step
	
	# system parameters
	ΩR  = 2π # Rabi frequency (rad * MHz)
	Γ = 0.15 # Measurement dephasing rate (MHz)
	τ = 1/(2Γ) #  Measurement collapse timescale (μs)
	η = 0.4 # collection efficiency
	
	# Kraus operators
	H = ΩR * σy/2
	J = [(σz, ((1 - η) * Γ))]
	C = [(exp(im * ϕ) * σz, τ, η)]
	
	Random.seed!(1)
	sol = bayesian((0, 4τ), ρ0, H, J, C; dt=dt)
	global B1 = traj(sol)
	
	md" ###### 🔻 Bayesian simulation (single-quadrature)"
end

# ╔═╡ 725dc4c3-cc74-4400-819c-2cffd06fbbf9
let 
	sim = B1
	tt, xx, yy, zz = (sim.t, sim.x, sim.y, sim.z)
	
	if show_gif
		anim = @animate for i ∈ range(1, length(tt), step=100)
			blochsphere(xx[1:i], yy[1:i], zz[1:i], linewidth=1., linealpha=0.85, ax=true, viewϕ = ϕv) end
		gif(anim, fps = 15)
	else
		blochsphere(xx, yy, zz, linewidth=1., linealpha=0.85, ax=true, viewϕ = ϕv)
	end
end

# ╔═╡ 34a700bb-5809-4755-a7fa-def102c5fd4c
let 
	sim = B1
	tt, xx, yy, zz = (sim.t, sim.x, sim.y, sim.z)
	
	if show_gif2
		anim = @animate for i ∈ range(1, length(tt), step=100)
			blochtimeseries(tt[1:i], xx[1:i], yy[1:i], zz[1:i], title = "Monitored Rabi oscillations", xlims = [0, last(tt)]) end
		gif(anim, fps = 15)
	else
		blochtimeseries(tt, xx, yy, zz, title = "Monitored Rabi oscillations")
	end
end

# ╔═╡ bb5f3187-2773-4647-807a-63141e16c2b4
let 
	sim = B1
	tt, xx, yy, zz = (sim.t, sim.x, sim.y, sim.z)
	
	if show_gif3
		anim = @animate for i ∈ range(1, length(tt), step=100)
			blochprojections(xx[1:i], yy[1:i], zz[1:i]) end
		gif(anim, fps = 15)
	else
		blochprojections(xx, yy, zz)
	end
end



# ╔═╡ dc3e2dcb-5bf1-492d-8337-f366dcf0170b
φdict = Dict("0" => 0, 
				"π/8" => π/8, 
				"π/4" => π/4,
				"3π/8" => 3π/8,
				"π/2" => π/2)

# ╔═╡ Cell order:
# ╠═4c45fe3e-cd69-11eb-20cd-7bfb98c040cf
# ╠═3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─377a3336-20bd-4baa-a033-af8bbc8668a8
# ╟─e43e5329-bd96-41ce-a183-1bd206204f65
# ╟─ecc19d12-e53f-4904-b13f-7a2ccc4912d7
# ╟─02d2ca0b-20bf-4dc2-a8f8-e72b5ae40ed2
# ╟─f3e7794e-a9a1-4103-8469-32a8e37d2d82
# ╟─01b79d3d-a6ba-4523-a668-b78a455279cb
# ╟─a12cdb8c-e9a1-4c2d-9811-cff266e152d8
# ╟─725dc4c3-cc74-4400-819c-2cffd06fbbf9
# ╟─0a7f28c9-1d84-43e6-b62c-711a231a3972
# ╟─34a700bb-5809-4755-a7fa-def102c5fd4c
# ╟─d9f2f00f-4ee2-45b5-91d6-6552d6d5b6c1
# ╟─bb5f3187-2773-4647-807a-63141e16c2b4
# ╟─133b6939-10b2-4c8e-acf8-5658ca96a0f9
# ╠═235cee23-c8af-4de2-a1c3-a2b173156703
# ╠═f3ac9b41-b123-4236-8c62-8e2b988463c2
# ╠═e27bd39c-58b7-4c5c-a677-3fe70f500ee8
# ╠═28d8df47-07cf-4a0d-a447-f894d021b2bc
# ╠═dc3e2dcb-5bf1-492d-8337-f366dcf0170b
