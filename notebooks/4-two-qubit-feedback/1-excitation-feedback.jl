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

# ╔═╡ 5093db05-a7f6-4708-b8e7-b939a6fa9d88
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

# ╔═╡ 377a3336-20bd-4baa-a033-af8bbc8668a8
md"""
In this interactive notebook, we'll explore generating two-qubit entangled states (Bell states) using weak measurement feedback on the total excitation number.
"""

# ╔═╡ 3030ad81-2a28-4f7f-b18b-80c757bf0575
mdp(table_of_contents📔)

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents(title="Two-qubit excitation feedback")

# ╔═╡ 009c334a-74cc-466c-9eb8-0f8aebfc9164
md"""
# System description
"""

# ╔═╡ be257103-3716-4ac0-aa94-845eb4892b55
md"""
##### Basis states
"""

# ╔═╡ 66d52a3a-edfd-4bdd-8c18-9f1a27041876
md"""
Recall that the two-qubit Hilbert space is spanned by four entangled states, called the Bell states:

$\ket{\Phi_\pm} \equiv \frac1{\sqrt2}\big(\ket{00} \pm \ket{11} \big)$

$\ket{\Psi_\pm} \equiv \frac1{\sqrt2}\big(\ket{01} \pm \ket{10} \big).$

These states have the special property that there are Rabi drives that swap between two states and leave the other two invariant, particularly:

 $\hat H_{+x+x} \equiv \Omega_R (\hat \sigma_1^x + \hat \sigma_2^x)$ swaps $\ket{\Phi_+} \leftrightarrow \ket{\Psi_+}$

 $\hat H_{+y+y} \equiv \Omega_R (\hat \sigma_1^y + \hat \sigma_2^y)$ swaps $\ket{\Phi_-} \leftrightarrow \ket{\Psi_+},$

which will be important for our feedback protocol later.
"""

# ╔═╡ 8bfa2445-5d97-45be-beee-a4361d2db190
md"""
##### Measurement
"""

# ╔═╡ 7190793d-6104-45cb-abb3-f1a174f882a6
md"""


We can count the total number of qubit excitations using the operator 

$\hat n = \hat n_1 + \hat n_2,$

where $n_i = (\hat \sigma_i^+ \hat \sigma_i^-)$ is a single-qubit number operator such that $n_i \ket{0} = 0$ and $n_i \ket{1} = \ket{1}$.
"""

# ╔═╡ deb2a0f3-32d1-46f9-bfbe-2e197ead243d
md"""
By weakly measuring the total number of dispersive shifts of the resonator coupled to both qubits, we can make a weak measurement of $\hat n$. We assume that the qubits do not interact with each other.

The measurement of $\hat n$ results in a noisy measurement record

$\tilde r(t) = \braket{\hat n} + \zeta(t), \hspace{5mm} \zeta(t) \sim \mathcal N(0, \tau/dt)$

where $\tau$ is the measurement collapse timescale determined by the resonator linewidth $\kappa$, the dispersive shift $\chi$, and the resonator population $\overline{n}$.
"""

# ╔═╡ 57b54539-03b6-49e3-84b2-291533cdac12
md"""
This measurement has eigenstates $\big\{\ket{00}, \alpha \ket{01} + \beta \ket{10}, \ket{11}) \big\}$, where $\alpha, \beta \in \mathbb{C}, |\alpha|^2 + |\beta|^2 = 1$ are complex coefficients reflecting the degeneracy of measuring $\hat n$: all these states have $1$ excitation, and are indistinguishable by the measurement.
"""

# ╔═╡ e8ee0728-e3c9-4307-9f81-1ffb1c3258b5
md"""
#### Expected behaviors
"""

# ╔═╡ 9cab60f8-d444-4786-951f-534666951f0e
md"""
First, we'll do a few checks to make sure the simulation is working as we expect:
* Under **measurement-only evolution**, the state should collapse to one of the measurement eigenstates: $\big\{\ket{00}, \alpha \ket{01} + \beta \ket{10}, \ket{11}) \big\}$
* Under **drive-only evolution**, product states should remain product states, with each qubit should oscillating between $\ket{0}$ and $\ket{1}$ independently; entangled states should remain entangled, and swap depending on the relative Rabi axes of both qubits.
"""

# ╔═╡ d708b8da-1d0f-4107-bfb6-35e132679cdd
md"""
# Checks
"""

# ╔═╡ 6bf29e9a-3f01-4e0f-b9c6-1007c8007ffb
md"""
## Measurement-only evolution
"""

# ╔═╡ f8b6a349-6c9b-45c6-9f72-61e7a4da67e8
md"""
### Measurement collapse
"""

# ╔═╡ 98328f7a-8b90-48fb-823e-e691c7edb1d5
md"""
If we initialize the state slightly away from a measurement eigenstate $| \Psi_+ \rangle$, it collapses to that state due to the measurement with high probability.
"""

# ╔═╡ ec88adf4-7e87-4326-b4f4-197b29cb6c2a
md"""
Arbitrary superpositions of $|01\rangle$ and $|10\rangle$ are also measurement eigenstates that are unchanged by the measurement:

$\ket{\psi} = \cos(\theta / 2) \ket{01} + e^{i \phi} \sin(\theta/2) \ket{10}$
"""

# ╔═╡ b9f11243-36c2-4570-bad3-0482ea923056
md"""
θ0 = 0 $(@bind θ0 Slider(0:0.1:3.14)) π

ϕ0 = 0 $(@bind ϕ0 Slider(0:0.1:6.28)) 2π
"""

# ╔═╡ fb8d78b4-ae13-4dd2-9df9-2c32dd134ee8
md"""
θ0 = $(round(θ0/π, digits=3)) π, ϕ0 = $(round(ϕ0/π, digits=3)) π
"""

# ╔═╡ 0e3196af-eaa0-4f60-bf5a-88340d91eeb7
md"""
### Collapse statistics
"""

# ╔═╡ 80e4d0d3-d8d1-4601-b49b-4657c4c7584e
md"""
Measurement collapses the state in a way that preserves the relative phase of $\ket{01}$ and $\ket{10}$ of the initial superposition, and also matches probabilities of the measurement eigenstates $\big\{\ket{00}, \alpha \ket{01} + \beta \ket{10}, \ket{11}) \big\}$.
"""

# ╔═╡ dce38821-9948-41b7-94f7-d78cd77d7617
function bulk_measure(ψ0::State; trajectory=false, N::Int64=50, tf=10.0)
	dt = 1e-3  # integration time-step

	Γ = 1.0 # Ensemble measurement dephasing rate (MHz)	
	η = 1.0 # collection efficiency
	ΩR = 0 # Rabi frequency (rad * MHz)
	
	Random.seed!()
	# Kraus operators --------------------------------------------------------------
	H = 0 * I
	J = [(n, ((1 - η) * Γ))]
	C = [(n, Γ, η)]

	map(1:N) do m
		sol = bayesian((0, tf), ψ0, H, J, C; dt=dt)
		return trajectory ? sol : last(sol.ρ) 
	end

	
end

# ╔═╡ afc6f9e7-e450-4936-9fce-5ec8a40da27a
begin
	labels = [L"\rho_0 = |0 0 \rangle",
			L"\rho_0 = (|0\rangle + |1 \rangle) \otimes (|0\rangle + |1 \rangle)",
			L"\rho_0 = (|0\rangle - |1 \rangle) \otimes (|0\rangle - |1 \rangle)", 
			L"\rho_0 = (|0\rangle + |1 \rangle) \otimes (|0\rangle - |1 \rangle)",
			L"\rho_0 = (|0\rangle - |1 \rangle) \otimes (|0\rangle + |1 \rangle)",
			L"\rho_0 = (|0 0\rangle - |1 0 \rangle - |0 1 \rangle - |1 1 \rangle",
			L"\rho_0 = (|1 1 \rangle"]
	
	init_states = [dm(g ⊗ g), 
					normalize(dm((g + e) ⊗ (g + e))), 
					normalize(dm((g - e) ⊗ (g - e))), 
					normalize(dm((g + e) ⊗ (g - e))), 
					normalize(dm((g - e) ⊗ (g + e))),
					normalize(dm(g ⊗ g - e ⊗ g - g ⊗ e - e ⊗ e)),
					dm(e ⊗ e)]

	fs_labels = [L"|00\rangle", L"\Psi_+\rangle", L"\Psi_-\rangle", L"|11\rangle"]

	rabi_labels = [L"| \Phi_+ \rangle",
				L"| \Psi_+ \rangle",
				L"| \Psi_- \rangle",
				L"| \Phi_{-} \rangle"]

	finalstates = [ket00, Ψp, Ψm, ket11]
end

# ╔═╡ b2703b3e-e4f4-4cc2-8f86-e38481fd4367
length(init_states)

# ╔═╡ b1bb0c2c-9eeb-4a1d-a5be-60a7bbe481c3
function statefreq(ρs, ket; threshold=0.99)
	evals = [real(expect(ρ, ket)) for ρ in ρs]
	length(evals[evals .> threshold])
end

# ╔═╡ 5dbe2d6a-5112-42d4-b5b2-8196e8f146dc
counted(sim) = sum([statefreq(sim, finalstate) for finalstate in finalstates])

# ╔═╡ 17f3a778-0d1c-40ca-9289-c1a7dbc51cb0
function convergencelist(ρfs, ket; threshold=0.99)
	evals = [real(expect(ρ, ket)) for ρ in ρfs]
	(evals .> threshold)
end

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

# ╔═╡ 33fc313f-5085-4ac8-bde6-0f8e7523e1a4
md"""
## Drive-only evolution
"""

# ╔═╡ 572ff0ac-ae60-44f0-9ef7-bd5b801522af
md"""
The Rabi drive creates transitions based on the axis ϕ. The drive below is $\Omega_R (\sigma_1^x + \sigma_2^x)$.
"""

# ╔═╡ 83f27327-336f-4988-9dd8-1f2601f7c773
md"""
Choose the inital subspace:
$(@bind initial Select(["entangled subspace", "product subspace"]))
"""

# ╔═╡ 40714981-bbfe-4e4b-b9b0-b4a07f63d8cd
if initial == "entangled subspace"
	
	md"When the state is initialized in $\ket{\Psi_+}$ (maximally entangled state), the Rabi drive rotates the state in the entangled subspace, while the reduced qubit states remain fully mixed. "

else

	md"When the state is initialized in $\ket{00}$ (product state), the Rabi drive rotates the state in the product subspace. Reduced qubit states remain pure and oscillate independently."
	
end

# ╔═╡ 8e418c51-3af3-4fe8-8f86-1987b2f6ecee
md"""
# Feedback
"""

# ╔═╡ 5c5d6110-c5d3-4a7c-adfb-bc266a868cbc
md"""
Below, we compare four different feedback protocols, corresponding to three different Hamiltonians:

0.  $\hat H = \pm \Omega_{max} (\hat \sigma_1^x + \hat \sigma_2^x)$, with the sign chosen to reduce the fidelity with $\psi_{target}$

1.  $\hat H = \Omega_{max} (\braket{n} - n_{target}) (\hat \sigma_1^x + \hat \sigma_2^x)$

2.  $\displaystyle \hat H = \Omega_{max} \Big( \frac{\theta(\rho)}{\pi} \Big)  (\hat \sigma_1^x + \hat \sigma_2^x)$, where $\theta(\rho) = \frac12 \text{arcsin} \Big\{ 2 \text{Tr}\big(\rho \ket{\Psi_+} \bra{\Phi_-}\big) \Big\}$ is the effective "angle" between $\ket{\Psi_+}$ and $\ket{\Phi_-}$

3.  $\hat H = \Omega_{max} \big(1 - F(\rho, \rho_{target})\big) (\hat \sigma_1^x + \hat \sigma_2^x)$, where $F(\rho, \rho_{target}) = \big(\text{Tr}\sqrt{\sqrt \rho \hspace{1mm} \rho_{target} \sqrt \rho}\big)^2$ is the fidelity with the target state $\rho_{target} = \ket{\psi_{target}} \bra{\psi_{target}}$
"""

# ╔═╡ 8baddd19-a12d-439c-8b34-d8e2a4a89ba6
md"""
using the following
##### simulation parameters:
"""

# ╔═╡ 832d3598-ce0c-46fb-af0d-8ab274789312
md" `Γ = 0.1 ` MHz $(@bind Γ Slider(0.1:0.1:3.0, default=1.0)) 3.0 MHz"

# ╔═╡ fb157192-c029-4766-b7d0-848e2351a68b
md" `Ωmax = 0.1` rad MHz $(@bind Ωmax Slider(0.1:0.1:3.0, default=1.57)) 3.0 rad MHz"

# ╔═╡ 5162c87e-2b32-4d7f-8d1a-3a8d7a5f36a3
md" `td = 0` ns $(@bind td Slider(0:1e-3:300e-3, default=0)) 300 ns"

# ╔═╡ 34863769-409d-408c-8d9f-e860a6adaf7e
md" `η = 0.01` $(@bind η Slider(0.01:0.01:1.0, default=1.0)) 1.0"

# ╔═╡ 50a8afd8-2a7d-4a71-883d-89b7139a119d
begin
	tf = 20.0   
	ψ0 = ket00
	dt = 1e-3
	target = Ψp
	ntarget = round(real(expect(n, target)))
end

# ╔═╡ 3578c71e-d386-4280-995a-d27889475edc
begin
	tf0 = 10.0
	N = 10
	sims = [bulk_measure(ψ0, N=N, tf=tf) for ψ0 in init_states]
end

# ╔═╡ ba421651-9ef4-4948-888b-a2efda39d49a
let
	mn = vcat([[statefreq(ρs, op) for op in finalstates] for ρs in sims]...)
	
	sx = repeat(labels, inner=length(fs_labels))
	
	nam = repeat(fs_labels, outer=length(sims))
	
	groupedbar(nam, mn, group = sx, ylabel = "counts", 
	        title = string("final state given initial state, tf = ", tf0, " μs, N = ", N), bar_width = 0.67, legendfontsize=11,tickfontsize=12, size=(800,400), legend=:outertopright, legendlabel="Initial state")
end

# ╔═╡ a94d1749-3420-4eb7-9369-ff4754fdb025
md" Ωmax = $Ωmax rad MHz, Γ = $Γ MHz, td = $(td*1000) ns, η = $η"

# ╔═╡ e36c59c6-4835-4875-92cc-627acb577310
md"""
## Protocols
"""

# ╔═╡ 18627c90-f5a8-4a25-af33-01fe51621922
md"""
### Optimizing sign only
"""

# ╔═╡ 384577f4-a4f5-46ea-b74e-466d687dbf8c
md"""
The sign of the Rabi drive is chosen to improve the fidelity with the target state on the next time step. The Rabi amplitude is not modulated. This works well for $t_d = 0$, but quickly breaks down for finite time delay.
"""

# ╔═╡ 77ac96ca-e2aa-4c6c-96b9-23810e113b0c
md"""
##### Simulation
"""

# ╔═╡ 314fc98c-ba13-40dc-9437-0ed735bf3698
let
	
	# Kraus operators --------------------------------------------------------------
	H(t::Float64, ρ::State) = 
		let
			h(sign) = sign * Ωmax * (σy1 + σy2)
			u(sign) = sign * exp( -im * dt * DenseOperator(h(sign)))
			ρnext(sign) = u(sign) * ρ * u(sign)'

			signs = [-1, 0, 1]
			sign = signs[findmax([fidelity(ρnext(sign), target) for sign in signs])[2]]

			h(sign)
		end

	
	J = [(n, ((1 - η) * Γ))]
	C = [(n, Γ, η)]

	global sol5 = bayesian((0.0, tf), ψ0, H, J, C; dt=dt, td=td)
end

# ╔═╡ ba2181b1-7a7f-4c07-9615-bb30e682556f
md"""
### On angle in Rabi plane
"""

# ╔═╡ 0c10665f-d1fc-4302-bcea-d1af9ad16e5e
md"""
$\displaystyle \hat H = \Omega_{max} \Big( \frac{\theta(\rho)}{\pi} \Big)  (\hat \sigma_1^x + \hat \sigma_2^x),$ where $\theta(\rho) = \frac12 \text{arcsin} \Big\{ 2 \text{Tr}\big(\rho \ket{\Psi_+} \bra{\Phi_-}\big) \Big\}$ is the effective "angle" between $\ket{\Psi_+}$ and $\ket{\Phi_-}$.

This protocol is optimal in the limit of zero time delay $t_d \rightarrow 0$, but seems about as effective as the $\braket n$ protocol for finite $t_d$.
"""

# ╔═╡ 8eb60364-c3f2-4515-9f07-44ec344ea41e
md"""
###### Simulation
"""

# ╔═╡ 5006da3f-1296-497c-8a2c-919f80afeefd
θ(ρ::State) = let 
			αβ = real(expect(Ψp ⊗ Φm', ρ))
			return 0.5 * asin(2αβ)
end

# ╔═╡ 8b2a6b25-5b5d-43ac-9d3d-ba2aaafed366
let
	Δt = 100e-3 # time duration of initial π/2 pulse
	margin = 0.01
	
	# Kraus operators --------------------------------------------------------------
	H(t::Float64, ρ::State) = let
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

	global sol3 = bayesian((0.0, tf), ψ0, H, J, C; dt=dt, td=td)
end

# ╔═╡ 5a01f315-4db7-4dac-8b98-891b5b8eef22
md"""
### On n average
"""

# ╔═╡ b5ac11c6-5a11-4806-b5d9-28fd70e5bf5c
md"""
$\hat H = \Omega_{max} (\braket{n} - n_{target}) (\hat \sigma_1^x + \hat \sigma_2^x).$

This protocol is intuitive to understand and works well even for finite time delay $t_d$. It is probably the easiest to compute and implement.
"""

# ╔═╡ e6c7e4e6-c523-4125-8dbc-82b23830bcfd
md"""
###### Simulation
"""

# ╔═╡ c931b23a-2925-4d91-b3a5-189188263649
let
	# Kraus operators --------------------------------------------------------------
	H(t::Timescale, ρ::QOp) = Ωmax * (real(expect(ρ, n)) - ntarget) * (σy1 + σy2)
	J = [(n, ((1 - η) * Γ))]
	C = [(n, Γ, η)]
	
	global sol2 = bayesian((0.0, tf), ψ0, H, J, C; dt=dt, td=td)
	
end

# ╔═╡ fa8a7c1b-8f25-4a96-8d82-bcfa9ce3b574
md"""
### On fidelity with target state
"""

# ╔═╡ 0ebfd0a0-bf2e-4a01-9919-b0ceb5a69c17
md"""
$\hat H = \Omega_{max} \big(1 - F(\rho, \rho_{target})\big) (\hat \sigma_1^x + \hat \sigma_2^x),$ where $F(\rho, \rho_{target}) = \big(\text{Tr}\sqrt{\sqrt \rho \hspace{1mm} \rho_{target} \sqrt \rho}\big)^2$ is the fidelity with the target state $\rho_{target} = \ket{\psi_{target}} \bra{\psi_{target}}.$

This is less effective on average, since it cannot take into account the optimal sign of the Rabi drive.
"""

# ╔═╡ 84a86253-e2f5-42cf-a986-c357a35053c4
md"""
###### Simulation
"""

# ╔═╡ cead7077-6c7e-4d7a-bcd2-90ca5f80a9a3
let	
	# Kraus operators --------------------------------------------------------------
	H(t::Timescale, ρ::QOp) = Ωmax * (1 - fidelity(ρ, target)) * (σy1 + σy2)
	J = [(n, ((1 - η) * Γ))]
	C = [(n, Γ, η)]
	
	# Random.seed!(s)
	global sol4 = bayesian((0.0, tf), ψ0, H, J, C; dt=dt, td=td)
end

# ╔═╡ 133b6939-10b2-4c8e-acf8-5658ca96a0f9
md" # Utilities"

# ╔═╡ ea04744b-4296-4dc2-9a3c-1f477c96f1ac
md"""
### Plotting
"""

# ╔═╡ 5a0be813-c725-4674-98c2-169701d14a40
begin
	colors1q = palette(:tab10)
	colors3q_number = palette(:lightrainbow)
	colors3q_prod = palette(:okabe_ito)
	colors2q_bell = palette(:rainbow)
end

# ╔═╡ 9050cb31-715a-45ee-8a0b-3136b4775742
function single_qubit_plots(sol::Solution)
	colors = colors1q 
	
	# calculate expectation values --------------------------------------------
	q1_basis = [σx1, σy1, σz1]
	q2_basis = [σx2, σy2, σz2]
	
	exps1 = map(op -> expectations(sol, op), q1_basis)
	exps2 = map(op -> expectations(sol, op), q2_basis)

	qlabels = ["x", "y", "z"]

	p1s = 0.5 * (1 .+ exps1[1].^2 .+ exps1[2].^2 .+ exps1[3].^2)
	p2s = 0.5 * (1 .+ exps2[1].^2 .+ exps2[2].^2 .+ exps2[3].^2)
	
	# plot ----------------------------------------------------------------------
	l = @layout [bloch1{0.5h}; bloch2{0.5h}]

	p1 = plot(sol.t, p1s, color=colors[4], label=L"Tr(\rho_1^2)", linestyle=:dash, title="single-qubit states")

	for l in 1:3
		label = qlabels[l]
		color = colors[l]
		exp = exps1[l]
		plot!(sol.t, exp, color=color, label=label, legend=:outerright, ylims=[-1,1])
	end

	p2 = plot(sol.t, p2s, color=colors[4], label=L"Tr(\rho_2^2)", linestyle=:dash, title="")

	for l in 1:3
		label = qlabels[l]
		color = colors[l]
		exp = exps2[l]
		plot!(sol.t, exp, color=color, label=label, legend=:outerright, xlabel="t (μs)", ylims=[-1,1])
	end
	
	plot(p1, p2, layout = l, link=:y, size=(600,300), legendfontsize=8, titlefontsize=12, legend=:outerright)
	
end

# ╔═╡ 3cbd88c3-f3e8-4692-8dc1-57e1197c1551
single_qubit_plots(sol5)

# ╔═╡ 815fab6b-8aa5-4fc8-9fe3-943f23239f88
single_qubit_plots(sol3)

# ╔═╡ 5bc8fe4c-85bb-47a6-95f2-8c99177956da
single_qubit_plots(sol2)

# ╔═╡ f3391975-49e7-4ff9-9f3f-2bdc64f329f0
single_qubit_plots(sol4)

# ╔═╡ c24c1980-714d-478d-8b1c-63af392b1534
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

# ╔═╡ 74a6a381-6a07-4682-90e4-40ea4f3ef6c6
let
	# parameters -----------------------------------------------------------------
	ψ0 = normalize(Ψp + 0.3 * ket00)
	dt = 1e-3  # integration time-step

	Γ = 1.0 # Ensemble measurement dephasing rate (MHz)	
	η = 1.0 # collection efficiency
	ΩR = 0.0
	tf = 10.0
	
	# Kraus operators --------------------------------------------------------------
	H = 0 * I
	J = [(n, ((1 - η) * Γ))]
	C = [(n, Γ, η)]
	
	sol = bayesian((0, tf), ψ0, H, J, C; dt=dt)

	bell_plot(sol)
end

# ╔═╡ 1de32aba-cc55-4026-b5ef-3a7896fff445
let
	# parameters -----------------------------------------------------------------
	ψ0 = cos(θ0/2) * (g ⊗ e) + sin(θ0/2)*exp(im * ϕ0) * (e ⊗ g)
	dt = 1e-3  # integration time-step

	Γ = 1.0 # Ensemble measurement dephasing rate (MHz)	
	η = 1.0 # collection efficiency
	ΩR = 0.0

	tf = 10.0
	
	# Kraus operators --------------------------------------------------------------

	H = 0 * I
	J = [(n, ((1 - η) * Γ))]
	C = [(n, Γ, η)]
	
	sol = bayesian((0, tf), ψ0, H, J, C; dt=dt)

	bell_plot(sol)
end

# ╔═╡ 9062da9f-167a-42d0-9337-ab130e0871d3
let
	# parameters -----------------------------------------------------------------
	dt = 1e-3  # integration time-step
	ψ0 = (initial == "entangled subspace") ? Ψp : ket00
	ΩR = π/2

	# Hamiltonian
	H = ΩR * (σx1 + σx2)

	global sol1 = bayesian((0, 5), ψ0, H, [], []; dt=dt)
	bell_plot(sol1)
end

# ╔═╡ ef5ad6d8-8de5-41b6-8847-a3ecde6031ad
single_qubit_plots(sol1)

# ╔═╡ b84ed4f6-5a27-486b-a842-e65726fcab14
bell_plot(sol5)

# ╔═╡ 4ce49e96-2bcc-44fc-a926-3006f726c92a
bell_plot(sol3)

# ╔═╡ 65b0dc75-a132-45e4-8820-d71a3b75cc93
bell_plot(sol2)

# ╔═╡ 1227c183-2724-47ee-9af4-1e6c5217c3e1
bell_plot(sol4)

# ╔═╡ Cell order:
# ╟─377a3336-20bd-4baa-a033-af8bbc8668a8
# ╟─3030ad81-2a28-4f7f-b18b-80c757bf0575
# ╟─3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─009c334a-74cc-466c-9eb8-0f8aebfc9164
# ╟─be257103-3716-4ac0-aa94-845eb4892b55
# ╟─66d52a3a-edfd-4bdd-8c18-9f1a27041876
# ╟─8bfa2445-5d97-45be-beee-a4361d2db190
# ╟─7190793d-6104-45cb-abb3-f1a174f882a6
# ╟─deb2a0f3-32d1-46f9-bfbe-2e197ead243d
# ╟─57b54539-03b6-49e3-84b2-291533cdac12
# ╟─e8ee0728-e3c9-4307-9f81-1ffb1c3258b5
# ╟─9cab60f8-d444-4786-951f-534666951f0e
# ╟─d708b8da-1d0f-4107-bfb6-35e132679cdd
# ╟─6bf29e9a-3f01-4e0f-b9c6-1007c8007ffb
# ╟─f8b6a349-6c9b-45c6-9f72-61e7a4da67e8
# ╟─98328f7a-8b90-48fb-823e-e691c7edb1d5
# ╟─74a6a381-6a07-4682-90e4-40ea4f3ef6c6
# ╟─ec88adf4-7e87-4326-b4f4-197b29cb6c2a
# ╟─b9f11243-36c2-4570-bad3-0482ea923056
# ╟─fb8d78b4-ae13-4dd2-9df9-2c32dd134ee8
# ╟─1de32aba-cc55-4026-b5ef-3a7896fff445
# ╟─0e3196af-eaa0-4f60-bf5a-88340d91eeb7
# ╟─80e4d0d3-d8d1-4601-b49b-4657c4c7584e
# ╟─dce38821-9948-41b7-94f7-d78cd77d7617
# ╠═b2703b3e-e4f4-4cc2-8f86-e38481fd4367
# ╠═3578c71e-d386-4280-995a-d27889475edc
# ╠═ba421651-9ef4-4948-888b-a2efda39d49a
# ╠═afc6f9e7-e450-4936-9fce-5ec8a40da27a
# ╟─b1bb0c2c-9eeb-4a1d-a5be-60a7bbe481c3
# ╟─5dbe2d6a-5112-42d4-b5b2-8196e8f146dc
# ╟─17f3a778-0d1c-40ca-9289-c1a7dbc51cb0
# ╟─cff9fded-2ecb-48fa-899d-202de36f74a9
# ╟─33fc313f-5085-4ac8-bde6-0f8e7523e1a4
# ╟─572ff0ac-ae60-44f0-9ef7-bd5b801522af
# ╟─83f27327-336f-4988-9dd8-1f2601f7c773
# ╟─40714981-bbfe-4e4b-b9b0-b4a07f63d8cd
# ╟─9062da9f-167a-42d0-9337-ab130e0871d3
# ╟─ef5ad6d8-8de5-41b6-8847-a3ecde6031ad
# ╟─8e418c51-3af3-4fe8-8f86-1987b2f6ecee
# ╟─5c5d6110-c5d3-4a7c-adfb-bc266a868cbc
# ╟─8baddd19-a12d-439c-8b34-d8e2a4a89ba6
# ╟─832d3598-ce0c-46fb-af0d-8ab274789312
# ╟─fb157192-c029-4766-b7d0-848e2351a68b
# ╟─5162c87e-2b32-4d7f-8d1a-3a8d7a5f36a3
# ╟─34863769-409d-408c-8d9f-e860a6adaf7e
# ╟─50a8afd8-2a7d-4a71-883d-89b7139a119d
# ╟─a94d1749-3420-4eb7-9369-ff4754fdb025
# ╟─e36c59c6-4835-4875-92cc-627acb577310
# ╟─18627c90-f5a8-4a25-af33-01fe51621922
# ╟─384577f4-a4f5-46ea-b74e-466d687dbf8c
# ╟─77ac96ca-e2aa-4c6c-96b9-23810e113b0c
# ╠═314fc98c-ba13-40dc-9437-0ed735bf3698
# ╟─b84ed4f6-5a27-486b-a842-e65726fcab14
# ╟─3cbd88c3-f3e8-4692-8dc1-57e1197c1551
# ╟─ba2181b1-7a7f-4c07-9615-bb30e682556f
# ╟─0c10665f-d1fc-4302-bcea-d1af9ad16e5e
# ╟─8eb60364-c3f2-4515-9f07-44ec344ea41e
# ╟─5006da3f-1296-497c-8a2c-919f80afeefd
# ╠═8b2a6b25-5b5d-43ac-9d3d-ba2aaafed366
# ╟─4ce49e96-2bcc-44fc-a926-3006f726c92a
# ╠═815fab6b-8aa5-4fc8-9fe3-943f23239f88
# ╟─5a01f315-4db7-4dac-8b98-891b5b8eef22
# ╟─b5ac11c6-5a11-4806-b5d9-28fd70e5bf5c
# ╟─e6c7e4e6-c523-4125-8dbc-82b23830bcfd
# ╟─c931b23a-2925-4d91-b3a5-189188263649
# ╟─65b0dc75-a132-45e4-8820-d71a3b75cc93
# ╟─5bc8fe4c-85bb-47a6-95f2-8c99177956da
# ╟─fa8a7c1b-8f25-4a96-8d82-bcfa9ce3b574
# ╟─0ebfd0a0-bf2e-4a01-9919-b0ceb5a69c17
# ╟─84a86253-e2f5-42cf-a986-c357a35053c4
# ╟─cead7077-6c7e-4d7a-bcd2-90ca5f80a9a3
# ╟─1227c183-2724-47ee-9af4-1e6c5217c3e1
# ╟─f3391975-49e7-4ff9-9f3f-2bdc64f329f0
# ╟─133b6939-10b2-4c8e-acf8-5658ca96a0f9
# ╟─ea04744b-4296-4dc2-9a3c-1f477c96f1ac
# ╠═5a0be813-c725-4674-98c2-169701d14a40
# ╠═9050cb31-715a-45ee-8a0b-3136b4775742
# ╠═c24c1980-714d-478d-8b1c-63af392b1534
# ╠═5093db05-a7f6-4708-b8e7-b939a6fa9d88
