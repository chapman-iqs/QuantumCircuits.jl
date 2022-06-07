### A Pluto.jl notebook ###
# v0.17.4

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

# ╔═╡ 0a9177a9-5d9c-47b5-999f-a5559d31a7d7
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
	using QuantumCircuits
	using Plots

	include("notebooks/table-of-contents.jl")
	include("notebooks/resources.jl")

	include("utilities/single-qubit-operators.jl")
	include("utilities/plotting.jl")
	
	md" # Packages and julia files"
end

# ╔═╡ 377a3336-20bd-4baa-a033-af8bbc8668a8
md"""
In this interactive notebook, we'll look at quantum trajectories on the Bloch sphere.
"""

# ╔═╡ f5244b50-85a1-4222-955e-5b455758ee25
mdp(table_of_contents📔)

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents(title="Bloch sphere")

# ╔═╡ 1cb8036a-3d95-439b-afdd-71de52f4286a
md"""
# Simulations
"""

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

# ╔═╡ 67741543-18d0-44df-9bc3-dd2df6bc4ff2
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

# ╔═╡ ecc19d12-e53f-4904-b13f-7a2ccc4912d7
let
	# simulation parameters
	ψ0 = g # initial state
	dt = 1e-3  # integration time-step
	
	# system parameters
	ΩR  = 2π # Rabi frequency (rad * MHz)
	Γ = 0.15 # Measurement dephasing rate (MHz)
	η = 0.4 # collection efficiency
	τ = 1/(2Γ * η) #  Measurement collapse timescale (μs)
	
	# Kraus operators
	H = (ΩR/2) * σy
	J = [(σz, (1 - η) * Γ)]
	C = [(exp(im * ϕ) * σz, Γ, η)]
	
	Random.seed!(1)
	sol = bayesian((0, 4τ), ψ0, H, J, C; dt=dt)
	global B1 = traj(sol)
	
	md" ###### 🔻 Bayesian simulation (single-quadrature)"
end

# ╔═╡ 725dc4c3-cc74-4400-819c-2cffd06fbbf9
let 
	sim = B1
	t = sim.t
	
	tt, xx, yy, zz = (sim.t, sim.x, sim.y, sim.z)
	
	if show_gif
		anim = @animate for i ∈ range(1, length(tt), step=100)
			blochsphere(x[1:i], y[1:i], z[1:i], linewidth=1., linealpha=0.85, ax=true, viewϕ = ϕv) end
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
		blochtimeseries(tt, xx, yy, zz, title = "Monitored Rabi oscillations", size=(600,300))
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



# ╔═╡ Cell order:
# ╟─377a3336-20bd-4baa-a033-af8bbc8668a8
# ╟─f5244b50-85a1-4222-955e-5b455758ee25
# ╟─3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─1cb8036a-3d95-439b-afdd-71de52f4286a
# ╠═ecc19d12-e53f-4904-b13f-7a2ccc4912d7
# ╟─02d2ca0b-20bf-4dc2-a8f8-e72b5ae40ed2
# ╟─f3e7794e-a9a1-4103-8469-32a8e37d2d82
# ╟─01b79d3d-a6ba-4523-a668-b78a455279cb
# ╟─a12cdb8c-e9a1-4c2d-9811-cff266e152d8
# ╠═725dc4c3-cc74-4400-819c-2cffd06fbbf9
# ╟─0a7f28c9-1d84-43e6-b62c-711a231a3972
# ╟─34a700bb-5809-4755-a7fa-def102c5fd4c
# ╟─d9f2f00f-4ee2-45b5-91d6-6552d6d5b6c1
# ╟─bb5f3187-2773-4647-807a-63141e16c2b4
# ╟─133b6939-10b2-4c8e-acf8-5658ca96a0f9
# ╠═67741543-18d0-44df-9bc3-dd2df6bc4ff2
# ╠═0a9177a9-5d9c-47b5-999f-a5559d31a7d7
