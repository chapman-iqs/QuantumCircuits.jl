### A Pluto.jl notebook ###
# v0.16.0

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
	σxq = sigmax(q)
	σyq = sigmay(q)
	σzq = sigmaz(q)
	σpq = sigmap(q)
	σmq = sigmam(q)
	Iq = identityoperator(q)
	
	ground = spindown(q)
	excited = spinup(q)

	md"###### 🔶 Qubit Hilbert space operators "
end

# ╔═╡ d9c286d7-fa58-4a0c-ba63-3b7f2b04b598
begin
	ρ0 = dm(spindown(q)) # initial state
	dt = 1e-3  # integration time-step
	md" ###### 🌀 Simulation parameters"
end

# ╔═╡ 9ed5fc83-eacd-4571-8806-42899fefa4bf
begin
	ΩR0  = 2π # Rabi frequency (rad * MHz)
	Γ0 = 0.15 # Measurement dephasing rate (MHz)
	τ0 = 1/(2Γ0) #  Measurement collapse timescale (μs)
	η0 = 0.4 # collection efficiency
	md" ###### 🌀 System parameters"
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

# ╔═╡ 23d6bc90-332d-47ca-a826-17d2fdeeaf51
begin
	H0 = ΩR0*σyq/2
	J0 = [(σzq, ((1-η0)*Γ0))]
	C0s = [(exp(im * ϕ) * σzq, τ0, η0)]
	md" ###### 💢 Kraus operators"
end

# ╔═╡ a12cdb8c-e9a1-4c2d-9811-cff266e152d8
md" $(@bind show_gif CheckBox()) Animate "

# ╔═╡ 0a7f28c9-1d84-43e6-b62c-711a231a3972
md" $(@bind show_gif2 CheckBox()) Animate Bloch series "

# ╔═╡ d9f2f00f-4ee2-45b5-91d6-6552d6d5b6c1
md" $(@bind show_gif3 CheckBox()) Animate cross-sections "

# ╔═╡ 133b6939-10b2-4c8e-acf8-5658ca96a0f9
md" # Utilities"

# ╔═╡ 8e85754f-d66b-477b-8153-b162519edb7c
expects(ops) = ρ -> collect(real(expect(ρ, s)) for s in vcat(ops, ρ)) # ρ -> [<x>,<y>,<z>,<ρ>]

# ╔═╡ ecc19d12-e53f-4904-b13f-7a2ccc4912d7
begin
	Random.seed!(1)
	sol1 = bayesian((0, 4τ0), ρ0, H0, J0, C0s; dt=dt)
	
	# collect outputs
	tt = sol1[1]
    ρt = sol1[2]
	r = collect(sol1[3][1])
	
	# get expectation values
	evs0 = expects([σxq, σyq, σzq]).(ρt);
    xx,yy,zz,ρρ = [map(x -> x[i], evs0) for i in 1:4];
	
	md" ###### 🔻 Bayesian simulation (single-quadrature)"
end

# ╔═╡ 725dc4c3-cc74-4400-819c-2cffd06fbbf9
let 
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
	if show_gif3
		anim = @animate for i ∈ range(1, length(tt), step=100)
			blochprojections(xx[1:i], yy[1:i], zz[1:i]) end
		gif(anim, fps = 15)
	else
		blochprojections(xx, yy, zz)
	end
end



# ╔═╡ Cell order:
# ╠═4c45fe3e-cd69-11eb-20cd-7bfb98c040cf
# ╠═3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─377a3336-20bd-4baa-a033-af8bbc8668a8
# ╟─e43e5329-bd96-41ce-a183-1bd206204f65
# ╟─d9c286d7-fa58-4a0c-ba63-3b7f2b04b598
# ╟─9ed5fc83-eacd-4571-8806-42899fefa4bf
# ╟─23d6bc90-332d-47ca-a826-17d2fdeeaf51
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
# ╟─8e85754f-d66b-477b-8153-b162519edb7c
