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

# ╔═╡ d3c8c931-efb2-4bb3-9ff1-b32ea59753cb
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
	
	include("utilities/plotting.jl")
	include("utilities/utilities.jl")
	include("utilities/single-qubit-operators.jl")

	
	include("notebooks/table-of-contents.jl")
	include("resources.jl")

	
	md" # Packages and julia files"
end

# ╔═╡ 64dd4d12-4473-4199-92f0-56ba1afd38dc
md"""
In this interactive notebook, we'll explore adiabatically changing the Rabi axis during evolution.
"""

# ╔═╡ 8166c912-cbd1-4f80-bacf-cb91820ae077
mdp(table_of_contents📔)

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents(title="Rabi dragging")

# ╔═╡ 2f02e4be-09b6-43ce-af0c-556626fa074d
md"""
# System description
"""

# ╔═╡ d529b5c0-18e8-41df-a345-b29c967fc49a
md"""
The state evolves according to the Hamiltonian

$H(t) = \Omega_R \hat \sigma_\phi(t)$

where $\hat \sigma_\phi(t) = cos(\omega t) \hat \sigma_x + sin(\omega t) \hat \sigma_y.$ This constitutes a drive of constant Rabi rate with axis rotating at rate $\omega$.

To satisfy adiabaticity, we require $\omega \ll \Omega_R$.

We'll initialize the state in $\ket{+x}$, so that it is an energy eigenstate under the Rabi drive. As we change the Rabi axis, we expect the state to get "dragged" along.
"""

# ╔═╡ fcc4e9a1-a063-409c-83ed-6152d9cf3f03
md"""
Note that the time-dependent Hamiltonian doesn't self-commute at different times: $[H(t_i), H(t_j)] \neq 0$ for $i \neq j$. So we must calculate unitary evolution in small steps, even in the absence of measurement or decoherence. `bayesian` implements this using unitaries $U_{\Delta t}(t) = e^{-i H(t) \Delta t}$, and evolves the state according to

$\ket{\psi(t_f)} = \prod_{t = t_0}^{t_f} U_{\Delta t}(t) \ket{\psi(0)}$

for pure states, or

$\rho(t_f) = \prod_{t = t_0}^{t_f} \mathcal{U_{\Delta t, t}} \big(\rho(0)\big)$

for mixed states, where $\mathcal{U_{\Delta t, t}} = U_{\Delta t}(t) \hspace{1mm} \cdot \hspace{1mm} U_{\Delta t}(t)^\dagger$.
"""

# ╔═╡ b06f3a5d-4691-42d9-a794-883d8655f593
md" $(@bind ideal CheckBox(default=true)) ideal "

# ╔═╡ dacb20f4-6960-43b0-a4cc-62b5047fcd9b
md"""
ω = 2π / d, d = 1 $(@bind denom Slider(1:1:100, default=70)) 100
"""

# ╔═╡ d08be14b-8389-4416-bc4f-6c884d9c449e
let
	# parameters -----------------------------------------------------------------
	ψ0 = normalize(g + e)
	dt = 1e-3  # integration time-step
	ΩR = 1 * 2π # Rabi frequency (rad * MHz)
	global ω = 2π / denom # rotation rate of Rabi axis -- `denom` on slider below

	global tf = 0.85 * (2π) / ω

	global T1 = 40
	global T2 = 60
	
	# Kraus operators --------------------------------------------------------------

	H(t) = ΩR * (cos(ω * t) * σx + sin(ω * t) * σy)
	J = ideal ? [] : [(σm, 1/T1), (σz, 1/(2T2))]
	C = []
	
	global solb = bayesian((0, tf), ψ0, H, J, C; dt=dt)
	
	
	md" ###### 🔻 Bayesian simulation"
end

# ╔═╡ 31d0c9ef-0ca0-4360-8e61-0a3eee2c6049
begin
	ops = (σx, σy, σz)
	(x, y, z) = map(op -> expectations(solb, op), ops)
	p = 0.5 .* (1 .+ x.^2 + y.^2 + z.^2)
	t = solb.t
end	

# ╔═╡ 987fc381-900f-4b17-9fc3-8cc66662b331
md"""
ω = 2π / $denom
"""

# ╔═╡ f3e7794e-a9a1-4103-8469-32a8e37d2d82
md"""

Rotate view: `ϕv = 0`
$(@bind ϕvc html"<input type=range min=0 max=32 step=1 value=0>") 
`ϕv = 2π`

"""

# ╔═╡ 01b79d3d-a6ba-4523-a668-b78a455279cb
begin
	ϕv = ϕvc * (π/16)
	md" `ϕv =` $(ϕvc/16) π"
end

# ╔═╡ a12cdb8c-e9a1-4c2d-9811-cff266e152d8
md" $(@bind show_gif CheckBox()) Animate "

# ╔═╡ 725dc4c3-cc74-4400-819c-2cffd06fbbf9
let 	
	if show_gif
		anim = @animate for i ∈ range(1, length(t), step=100)
			blochsphere(x[1:i], y[1:i], z[1:i], linewidth=1., linealpha=0.85, ax=true, viewϕ = ϕv, vec=(cos(ω * t[i]), sin(ω * t[i]), 0), blochmark = true) end
		gif(anim, fps = 40)
	else
		blochsphere(x, y, z, linewidth=1., linealpha=0.85, ax=true, viewϕ = ϕv, vec=(cos(ω * tf), sin(ω * tf), 0), blochmark = true)
	end
end

# ╔═╡ 0a7f28c9-1d84-43e6-b62c-711a231a3972
md" $(@bind show_gif2 CheckBox()) Animate Bloch series "

# ╔═╡ d9f2f00f-4ee2-45b5-91d6-6552d6d5b6c1
md" $(@bind show_gif3 CheckBox()) Animate cross-sections "

# ╔═╡ bb5f3187-2773-4647-807a-63141e16c2b4
if show_gif3
	anim = @animate for i ∈ range(1, length(t), step=100)
		blochprojections(x[1:i], y[1:i], z[1:i], vec=(cos(ω * t[i]), sin(ω * t[i]), 0)) end
	gif(anim, fps = 15)
else
	blochprojections(x, y, z, vec=(cos(ω * tf), sin(ω * tf), 0))
end

# ╔═╡ 133b6939-10b2-4c8e-acf8-5658ca96a0f9
md" # Utilities"

# ╔═╡ 3eab3404-fbf0-45c6-a72b-b5e43af4380a
colors = palette(:tab10)

# ╔═╡ 4da1b0ce-c067-412f-8228-0037ac1d02b3
if ideal

	plot(t, x, color=colors[1], label=L"x", legendfontsize=12, xlabel="t (μs)", title=string("ω = ΩR / ", denom, ", T1 = ", (ideal ? "∞" : T1 ), " μs, T2 = ", (ideal ? "∞" : T2), " μs"))
	plot!(t, (t -> cos(ω * t)).(t), color=:black, linestyle=:dot, linewidth=2, label=L"cos(\omega t)")
	plot!(t, y, color=colors[2], label=L"y")
	plot!(t, (t -> sin(ω * t)).(t), color=:black, linestyle=:dot, linewidth=2, label=L"sin(\omega t)")
	plot!(t, z, color=colors[3], label=L"z")
	plot!(t, p, color=colors[4], label=L"Tr(\rho^2)")

else
	plot(t, x, color=colors[1], label=L"x", legendfontsize=12, xlabel="t (μs)", title=string("ω = ΩR / ", denom, ", T1 = ", (ideal ? "∞" : T1 ), " μs, T2 = ", (ideal ? "∞" : T2), " μs"))
	plot!(t, (t -> exp(-t/(2T2) - t/T1) * cos(ω * t)).(t), color=:black, linestyle=:dot, linewidth=2, label=L"cos(\omega t)")
	plot!(t, y, color=colors[2], label=L"y")
	plot!(t, (t -> exp(-t/(2T2) - t/T1) * sin(ω * t)).(t), color=:black, linestyle=:dot, linewidth=2, label=L"sin(\omega t)")
	plot!(t, z, color=colors[3], label=L"z")
	plot!(t, p, color=colors[4], label=L"Tr(\rho^2)")
	
end

# ╔═╡ e27bd39c-58b7-4c5c-a677-3fe70f500ee8
xyz(θ, ϕ) = (sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ))

# ╔═╡ Cell order:
# ╟─64dd4d12-4473-4199-92f0-56ba1afd38dc
# ╟─8166c912-cbd1-4f80-bacf-cb91820ae077
# ╟─3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─2f02e4be-09b6-43ce-af0c-556626fa074d
# ╟─d529b5c0-18e8-41df-a345-b29c967fc49a
# ╟─fcc4e9a1-a063-409c-83ed-6152d9cf3f03
# ╟─b06f3a5d-4691-42d9-a794-883d8655f593
# ╠═d08be14b-8389-4416-bc4f-6c884d9c449e
# ╠═31d0c9ef-0ca0-4360-8e61-0a3eee2c6049
# ╟─dacb20f4-6960-43b0-a4cc-62b5047fcd9b
# ╟─987fc381-900f-4b17-9fc3-8cc66662b331
# ╟─f3e7794e-a9a1-4103-8469-32a8e37d2d82
# ╟─01b79d3d-a6ba-4523-a668-b78a455279cb
# ╟─a12cdb8c-e9a1-4c2d-9811-cff266e152d8
# ╟─725dc4c3-cc74-4400-819c-2cffd06fbbf9
# ╟─0a7f28c9-1d84-43e6-b62c-711a231a3972
# ╟─4da1b0ce-c067-412f-8228-0037ac1d02b3
# ╟─d9f2f00f-4ee2-45b5-91d6-6552d6d5b6c1
# ╟─bb5f3187-2773-4647-807a-63141e16c2b4
# ╟─133b6939-10b2-4c8e-acf8-5658ca96a0f9
# ╠═3eab3404-fbf0-45c6-a72b-b5e43af4380a
# ╠═e27bd39c-58b7-4c5c-a677-3fe70f500ee8
# ╠═d3c8c931-efb2-4bb3-9ff1-b32ea59753cb
