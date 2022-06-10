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

# ╔═╡ ce94b9f1-4612-4bbd-af04-0c9a7e433939
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
	using Random
	using LaTeXStrings
	using Statistics
	using Distributions
	using QuantumCircuits
	using Plots
	using Plots.Measures
	using StatsPlots
	
	include("utilities/two-qubit-operators.jl")
	include("utilities/utilities.jl")
	include("utilities/plotting.jl")

	include("notebooks/table-of-contents.jl")
	include("resources.jl")
	
	md" # Packages and julia files"
end

# ╔═╡ acb801e9-5ed1-40cf-8b14-56def620df00
mdp("This notebooks tests the system-filter functionality in `QuantumCircuits.jl`, using the excitation feedback as an example. (See ", excitation_feedback📔, " for a description of this system.)")

# ╔═╡ f53eaaa9-e5f2-4964-93f6-32f3968b8f5f
md"""
The terminology of "system / filter" comes from $Szigeti_et_al. We keep that terminology here, though "system / estimate" may be clearer, since both the "system" state and the "filter" state arise through Bayesian filtering of a measurement record.
"""

# ╔═╡ 920e897f-c383-4cb9-9f18-ee42230f36c7
mdp(table_of_contents📔)

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents(title="System / filter")

# ╔═╡ ac42fe33-b241-4a68-b167-f58b78a634b4
md"""
#### Simulation parameters
"""

# ╔═╡ ac0008af-d797-411e-b329-272d3f7a47e0
md" `Γ = 0.1 ` MHz $(@bind Γ Slider(0.1:0.1:3.0, default=1.0)) 3.0 MHz"

# ╔═╡ c40b9f1c-86a8-45e6-bf3f-3d344267d20d
md" `Ωmax = 0.1` rad MHz $(@bind Ωmax Slider(0.1:0.1:3.0, default=1.57)) 3.0 rad MHz"

# ╔═╡ 91937177-d607-418a-9753-1a3f5f18d9e8
md" `td = 0` ns $(@bind td Slider(0:1e-3:300e-3, default=0)) 300 ns"

# ╔═╡ 54a42a51-4be4-4f77-b7a5-a45e1358f54d
md" `η = 0.01` $(@bind η Slider(0.01:0.01:1.0, default=1.0)) 1.0"

# ╔═╡ df7e9778-d3e2-4fea-854b-314382b595c9
md""" $\Omega_{max}$ = $Ωmax rad MHz, 

 $\Gamma$ = $Γ MHz,
 
 $t_d$ = $(td*1000) ns,
 
 $\eta$ = $η
"""

# ╔═╡ e47c00dc-09b6-4b3b-a132-c201434014da
begin
	tf = 20.0   
	ψ0 = ket00
	dt = 1e-3
	target = Ψp
	ntarget = round(real(expect(n, target)))
	ψ0 = ket00
	Δt = 100e-3 # time duration of initial π/2 pulse
	margin = 0.01
end

# ╔═╡ 582ffe51-cb4e-4b28-a633-acc9a05f74e1
md"""
# "System" and "filter" Hamiltonians identical

In this case, we expect that the trajectories should be exactly identical, since they share a measurement record.
"""

# ╔═╡ 0e48d9ef-e8d3-45cc-96d2-29d4930c331b
md"""
### Angle feedback
"""

# ╔═╡ a772504d-23c4-4a7c-93e8-571fa8438120
θ(ρ::State) = let 
			αβ = real(expect(Ψp ⊗ Φm', ρ))
			return 0.5 * asin(2αβ)
end

# ╔═╡ c08cda77-0ea0-4c4e-8a2b-87acf0c14cda
let
	
	# Kraus operators --------------------------------------------------------------
	H(t::Timescale, ρ::State) = let
		 Ω = 
			if t < Δt
				π / (4Δt)
			elseif fidelity(target, ρ) < margin && t > td + Δt
				π / (4Δt) 
				# "instantaneous" π/2 pulse to avoid Zeno pinning to unwanted state
			elseif td == 0.0
				- θ(ρ) / (2dt)
			else
				- θ(ρ) * Ωmax / π
			end

		return Ω * (σy1 + σy2)

	end

	
	J = η == 1.0 ? [] : [(n, ((1 - η) * Γ))]
	C = [(n, Γ, η)]

	global (sys1, fil1) = bayesian((0.0, tf), ψ0, (H, H), J, C; dt=dt, td=td)
end

# ╔═╡ 822fa001-6e65-4e04-b128-c4006944f601
sys1.ρ == fil1.ρ

# ╔═╡ 40683f6e-96f3-4832-be4a-83ead9e86532
md"""
### Number feedback
"""

# ╔═╡ d02cd7c4-ebf5-4592-baa2-01232c63627c
let
	# Kraus operators --------------------------------------------------------------
	H(t::Timescale, ρ::QOp) = Ωmax * (real(expect(n, ρ)) - ntarget) * (σy1 + σy2)
	J = [(n, ((1 - η) * Γ))]
	C = [(n, Γ, η)]
	
	global (sys2, fil2) = bayesian((0.0, tf), ψ0, (H, H), J, C; dt=dt, td=td)
	
end

# ╔═╡ 7616c011-35ba-44c0-9ab5-019b3bc3beb3
sys2.ρ == fil2.ρ

# ╔═╡ e8aa97bb-00d2-4dc1-b946-9c9c390d43c3
md"""
# Different Hamiltonians

Now, let's look at what happens when we give them different Hamiltonians, which is the point of the system-filter funcitonality. In the simulations below, we define a filter (feedback) Hamiltonian $H_f(t, \rho)$ identically to the two previous examples. However, we also define a system Hamiltonian $H_s(t, \rho) = H_f(t, \rho) + H_z(t)$ such that

$H_z(t) = R_1 \hspace{1mm} \hat \sigma_1^z + R_2 \hspace{1mm}  \hat \sigma_2^z,$ 

where $R_1, R_2 \sim \mathcal{N}(0, \sigma)$ are i.i.d. Gaussian-distributed random variables with standard deviation $\sigma$.

This models a scenario in which an experimental controller outputs a feedback drive based on limited knowledge of the actual evolution. The controller assumes the two qubits evolve only according to the controller-output Hamiltonian $H_f(t, \rho)$. The state known by the controller is called the "filtered" state $\rho_f$. However, the "actual" system state $\rho_s$ is subject to additional phase noise $H_z(t)$, in addition to the controller drive $H_f(t, \rho)$. In both cases, the readout (and thus the measurement backaction) is determined by the system state, i.e.

$r_t = \text{Tr}(\hat n \rho_s) + \zeta_t, \hspace{5mm} \zeta_t \sim \mathcal{N}(0, \sqrt{\tau/dt}).$

This method of modeling system and filter allows us to evaluate how certain types of experimentally un-measurable noise may disrupt the feedback process. While a Hamiltonian of the type $H_z(t)$ models $T_2$ (dephasing) effects in the ensemble-average, the standard way of modeling such effects using $\sigma_z$ Lindblad operators fails to capture the effect of such noise on feedback efficacy for individual trajectories, since Lindblads only capture ensemble-averaged effects.

"""

# ╔═╡ 8c2af690-5161-4e31-a9bf-566773dfd38f
md"""
Choose the noise level for the system:
$(@bind noise_str Select(["σ = 1.0 (Τ2 ~ 100 μs)", "σ = 1.5 (T2 ~ 40 μs)", "σ = 5.0 (T2 ~ 10 μs)", "σ = 10.0 (T2 ~ 2 μs)"]))
"""

# ╔═╡ 7b109ad4-9be0-4266-b0c1-09a28c0772e3
begin
	noise_dict = Dict(["σ = 1.0 (Τ2 ~ 100 μs)" => 1.0, 
					"σ = 1.5 (T2 ~ 40 μs)" => 1.5,
					"σ = 5.0 (T2 ~ 10 μs)" => 5.0, 
					"σ = 10.0 (T2 ~ 2 μs)" => 10.0])
	σ = noise_dict[noise_str]    
end

# ╔═╡ 3c763b47-0206-403c-93cd-78671226dcfd
md"""
### Angle feedback
"""

# ╔═╡ d7f43d14-97dd-47d9-8e7a-173403fbfdee
begin
	
	Hf(t::Timescale, ρ::State) = let
			 Ω = 
				if t < Δt
					π / (4Δt)
				elseif fidelity(target, ρ) < margin && t > td + Δt
					π / (4Δt) 
					# "instantaneous" π/2 pulse to avoid Zeno pinning to unwanted state
				elseif td == 0.0
					- θ(ρ) / (2dt)
				else
					- θ(ρ) * Ωmax / π
				end
	
			return Ω * (σy1 + σy2)
	
		end
	
	Hz(t, σ) = let
			R1, R2 = σ * randn(), σ * randn()
			R1 * σz1 + R2 * σz2
	end
	
	Hs(σ) = (t::Timescale, ρ::State) -> Hf(t, ρ) + Hz(t, σ)
	
	J = η == 1.0 ? [] : [(n, ((1 - η) * Γ))]
	C = [(n, Γ, η)]

	global (sys3, fil3) = bayesian((0.0, tf), ψ0, (Hs(σ), Hf), J, C; dt=dt, td=td)
end

# ╔═╡ 573858ac-31f5-41f4-9dad-b3673223e4f5
md"""
### Number feedback
"""

# ╔═╡ 6f30b61b-b12a-4e37-ab75-021b7303a2f2
let
	# Kraus operators --------------------------------------------------------------
	Hf(t::Timescale, ρ::State) = Ωmax * (real(expect(n, ρ)) - ntarget) * (σy1 + σy2)

	Hz(t, σ) = let
			R1, R2 = σ * randn(), σ * randn()
			R1 * σz1 + R2 * σz2
	end
	
	Hs(σ) = (t::Timescale, ρ::State) -> Hf(t, ρ) + Hz(t, σ)
	
	J = η == 1.0 ? [] : [(n, ((1 - η) * Γ))]
	C = [(n, Γ, η)]

	global (sys4, fil4) = bayesian((0.0, tf), ψ0, (Hs(σ), Hf), J, C; dt=dt, td=td)
	
end

# ╔═╡ 133b6939-10b2-4c8e-acf8-5658ca96a0f9
md" # Utilities"

# ╔═╡ ea04744b-4296-4dc2-9a3c-1f477c96f1ac
md"""
## Plotting
"""

# ╔═╡ 8cb229d1-9303-4d53-b5e6-6c5255b66a68
begin
	colors1q = palette(:tab10)
	colors_system = palette(:rainbow)
	colors_filter = palette(:lightrainbow)	
end

# ╔═╡ 71208ed0-f538-4c32-9c5e-aa2e65008b23
function bell_plot(sol::Solution)

	basis = bell_basis
	colors = colors2q_bell
	labels = bell_basis_labels
	title = "Bell states"

	exps = map(op -> expectations(sol, dm(op)), basis)
	
	pl = plot(size=(600,300), legendfontsize=12, titlefontsize=12, legend=:outerright, title=title)
	
	for l in 1:length(basis)
		label = labels[l]
		color = colors[l]
		exp = exps[l]
		plot!(sol.t, exps[l], color=color, label=label, legend=:outerright, xlabel="t (μs)", ylims=[0,1])
	end

	pl

end

# ╔═╡ e688adcc-43c4-4c9c-a541-c867bf9cbf96
function bell_plot(sys::Solution, fil::Solution)

	basis = bell_basis
	labels = bell_basis_labels
	title = "Bell states"

	exps_sys = map(op -> expectations(sys, dm(op)), basis)
	exps_fil = map(op -> expectations(fil, dm(op)), basis)
	
	pl = plot(size=(600,300), legendfontsize=12, titlefontsize=12, legend=:outerright, title=title, xlabel="t (μs)",)
	
	for l in 1:length(basis)
		plot!(sys.t, exps_sys[l], color=colors_system[l], label=labels[l], legend=:outerright, ylims=[0,1])
		plot!(fil.t, exps_fil[l], color=colors_filter[l], label=:none, linestyle=:dash)
	end

	pl

end

# ╔═╡ f1a39f79-5b3c-46a8-bb09-c28a7fb21f6e
bell_plot(sys1, fil1)

# ╔═╡ 460e5547-a89d-4213-948f-892e68da8c9c
bell_plot(sys2, fil2)

# ╔═╡ 76627d20-7db7-4ee4-baab-c3d93c774a55
bell_plot(sys3, fil3)

# ╔═╡ bd38356c-d166-468a-8c6a-b5f53c926bf0
bell_plot(sys4, fil4)

# ╔═╡ Cell order:
# ╟─acb801e9-5ed1-40cf-8b14-56def620df00
# ╟─f53eaaa9-e5f2-4964-93f6-32f3968b8f5f
# ╟─920e897f-c383-4cb9-9f18-ee42230f36c7
# ╟─3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─ac42fe33-b241-4a68-b167-f58b78a634b4
# ╟─ac0008af-d797-411e-b329-272d3f7a47e0
# ╟─c40b9f1c-86a8-45e6-bf3f-3d344267d20d
# ╟─91937177-d607-418a-9753-1a3f5f18d9e8
# ╟─54a42a51-4be4-4f77-b7a5-a45e1358f54d
# ╟─df7e9778-d3e2-4fea-854b-314382b595c9
# ╠═e47c00dc-09b6-4b3b-a132-c201434014da
# ╟─582ffe51-cb4e-4b28-a633-acc9a05f74e1
# ╟─0e48d9ef-e8d3-45cc-96d2-29d4930c331b
# ╟─a772504d-23c4-4a7c-93e8-571fa8438120
# ╠═c08cda77-0ea0-4c4e-8a2b-87acf0c14cda
# ╟─f1a39f79-5b3c-46a8-bb09-c28a7fb21f6e
# ╠═822fa001-6e65-4e04-b128-c4006944f601
# ╟─40683f6e-96f3-4832-be4a-83ead9e86532
# ╠═d02cd7c4-ebf5-4592-baa2-01232c63627c
# ╠═460e5547-a89d-4213-948f-892e68da8c9c
# ╠═7616c011-35ba-44c0-9ab5-019b3bc3beb3
# ╟─e8aa97bb-00d2-4dc1-b946-9c9c390d43c3
# ╟─8c2af690-5161-4e31-a9bf-566773dfd38f
# ╟─7b109ad4-9be0-4266-b0c1-09a28c0772e3
# ╟─3c763b47-0206-403c-93cd-78671226dcfd
# ╠═d7f43d14-97dd-47d9-8e7a-173403fbfdee
# ╠═76627d20-7db7-4ee4-baab-c3d93c774a55
# ╟─573858ac-31f5-41f4-9dad-b3673223e4f5
# ╠═6f30b61b-b12a-4e37-ab75-021b7303a2f2
# ╠═bd38356c-d166-468a-8c6a-b5f53c926bf0
# ╟─133b6939-10b2-4c8e-acf8-5658ca96a0f9
# ╠═ea04744b-4296-4dc2-9a3c-1f477c96f1ac
# ╟─8cb229d1-9303-4d53-b5e6-6c5255b66a68
# ╟─71208ed0-f538-4c32-9c5e-aa2e65008b23
# ╟─e688adcc-43c4-4c9c-a541-c867bf9cbf96
# ╠═ce94b9f1-4612-4bbd-af04-0c9a7e433939
