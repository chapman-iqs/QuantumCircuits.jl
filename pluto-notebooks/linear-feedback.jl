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
	using DataFrames
	
	include("plotting.jl")
	
	md" ###### 🔶 Packages and julia files"
end

# ╔═╡ 88b77211-fe49-4afb-9c2a-c64791c20395
using CSV

# ╔═╡ 377a3336-20bd-4baa-a033-af8bbc8668a8
md"""
# Linear feedback stabilization

In this interactive notebook, we'll look at dynamical linear feedback stabilization of quantum trajectories, based on [1].
"""

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents()

# ╔═╡ 8f7c6440-eac8-49d3-bd23-c2545ec16830
md"""
## Measurement backaction principles
"""

# ╔═╡ ff721d9a-c37c-43a1-adee-015221ddd133
md"""
To understand how linear feedback stabilization works, we must first look at how measurement creates backaction in a qubit.
"""

# ╔═╡ 727ced31-7a51-4892-b1d9-eff4fcdc9767
md"""
Weak measurement of a qubit observable $\hat A$ generates a stochastic record

$\tilde r (t) = \braket{\hat A}(t) + \zeta(t)$

where $\braket{\hat A}(t)  \equiv \text{Tr} \big(\rho_q(t) \hat A \big)$ is the expectation value of the observable at time $t$, and 

$\zeta(t) \sim \mathcal{N} \Big(0, \sqrt{\tau_m/ dt} \Big)$

is a zero-mean, Gaussian distributed random variable with variance $\sigma = \tau_m / dt$. Here, $\tau_m$ is the measurement timescale determined by the strength of the measurement.
"""

# ╔═╡ 1aaa9cbe-2f98-41e1-af41-598ba6578818
md"""
In this notebook, we will consider measurement of the qubit state $\ket0$ or $\ket1$, such that $\hat A = e^{i \varphi} \hat \sigma_z$. The phase $\varphi$ is called the "measurement angle" and correpsonds to the quadrature of amplification of the signal. It determines the type of backaction experienced by the qubit.
"""

# ╔═╡ 1ee143ba-2f0d-4eee-a76c-4bf58469e0be
md"""
### Simulation
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
	id = identityoperator(q)
	
	ground = spindown(q)
	excited = spinup(q)

	md"###### 🔶 Qubit Hilbert space operators "
end

# ╔═╡ 875adb1b-22f7-4375-83c5-2a17c9a8e46c
md"""

Rotate view: `ϕv = 0`
$(@bind ϕvIB html"<input type=range min=0 max=32 step=1 value=0>") 
`ϕv = 2π`

"""

# ╔═╡ 3531558a-f9b7-4f6b-9800-884fe0b04712
md" $(@bind animIB_bs CheckBox()) Animate Bloch sphere "

# ╔═╡ faf0b11c-4339-4a82-a23c-9d35eb9d10b4
md" $(@bind animIB CheckBox(default=false)) Animate informational backaction"

# ╔═╡ 911d4aa1-41cd-4eea-8e0f-1af6f8024290
@bind go Button("Rerun simulation")

# ╔═╡ 3ae402b0-476a-4c11-823c-d0016793adab
md"""
Change the measurement angle to see how backaction depends on the measurement quadrature.
"""

# ╔═╡ 71604112-0704-4832-8a80-e886b6155e34
md"""
`φ =` $(@bind φIB Select(["0", "π/8", "π/4", "3π/8", "π/2"]))
"""

# ╔═╡ d3e48ac4-973f-4397-bf06-89150a0c13ed
md"""
Informational backaction occurs when $\varphi = 0$. Phase backaction occurs when $\varphi = \pi/2$. Intermediate angles result in a combination of both types of backaction.

In the following simulation, the state is initialized in the state $\ket{+x}$. By default, the measurement angle is $\varphi = 0$, leading to informational backaction. The measurement angle is currently `φ =` $φIB.

Interact with the following plots and try different parameter settings (below) to gain intuition for informational and phase backaction.
"""

# ╔═╡ 3ff551a9-3a07-4a2f-928e-880b7e3ba1fc
if 	φIB == "0"
	
	md""" **Informational backaction**  The qubit undergoes a biased random walk towards one of the poles. For $\eta = 1$, the walk is confined to the circle $x^2 + z^2 = 1$. Once it reaches a pole, it is pinned there. As an aside, this demonstrates how continuous weak monitoring can be understand as a projective measurement in the limit of measuring for a long time. The statistics of final states will follow that of projective measurement, in the absence of Hamiltonian evolution.

While the initial state creates a bias towards the pole it is closer to, it is possible to stabilize to either pole so long as the initial state is a superposition. Click the button below to see this:
	"""

elseif φIB == "π/2"
	
	md""" **Phase backaction**  The qubit undergoes a unbiased random walk in the $x$ -- $y$ plane.  For $\eta = 1$, the walk is confined to the circle $x^2 + y^2 = 1$. Note that this effect is preserved for a state initialized outside of the $x$ -- $y$ plane; the random walk will then be confined to its original plane: $x^2 + y^2 = 1 - z^2$.

	"""

else
	
	md""" **Informational and phase backaction**  The qubit undergoes a biased random walk towards one of the poles.  For $\eta = 1$, the walk is confined to the sphere $x^2 + y^2 + z^2 = 1$. 

	"""

end

# ╔═╡ ab188d7c-9eb4-48f6-86e6-53c057ce67ae
md"""
## Linear feedback
"""

# ╔═╡ 5b83ae62-9486-485f-821b-eb127f01e487
md"""
Linear feedback takes advantage of the fact that informational (phase) backaction looks like random rotations about the $x$ - or $y$ -axes ($z$ -axis). Such rotations can be straightforwardly implemented via Hamiltonian evolution. Thus, feeding back the measurement record into the Hamiltonian can be used to "erase" certain effects of measurement backaction.

Informational linear feedback, implemented in [1], erases the effects of measurement backaction while changing the effective measurement pole. First, we will gain intuition for how this works in the next section. Then, we will attempt to apply similar principles to stabilize phase backaction.
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

Later in this notebook, we will look into whether a similar scheme might be used to stabilize phase backaction. This would allow for stabilization of qubits measured via phase-preserving amplification, which results inevitably in both informational and phase backaction.

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
ideal = false

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

# ╔═╡ 219487b2-c85c-4b64-8387-5d8a90801dc6
md"""
### Ensemble
"""

# ╔═╡ 98baaa59-213b-4f6e-929b-db10ad4fd718
N = 10

# ╔═╡ d7e4182d-46ba-40f5-97bd-ced14022b3f0
md"""
###### Export DataFrame
"""

# ╔═╡ 9638226c-7c76-4230-9cc9-cc6700232e20
testd = Dict([("A", 1), ("B", 2)])

# ╔═╡ 2ef48cde-d72e-4d91-8de3-5fcb52d764eb
testd."A"

# ╔═╡ a8b375ae-5203-4c4f-96c8-e594c5eb098d
md"""
###### Import dataframe
"""

# ╔═╡ e8a25893-452f-461f-9050-ee16b6f80f54
newx = DataFrame(CSV.File("data/xdata.csv"))

# ╔═╡ 9279175b-0bb8-4f4d-b1ed-6289aacdaacf
newz = DataFrame(CSV.File("data/zdata.csv"))

# ╔═╡ 849589c4-d260-4f24-a1c5-cf4887ba9701
md"""
### Scan fidelities over η
"""

# ╔═╡ 84275f2c-a215-4d80-ba81-326ed0eab8c3
N_ηs = 100

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
idealDQM = false

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

# ╔═╡ 723b3eab-d717-4823-a7eb-da39202f6f7b
md"""
### Scan fidelities over ΔS
"""

# ╔═╡ 1f017209-c017-4d77-9a94-8edde2c02fdb
NΔS = 100

# ╔═╡ 9e98e449-859e-4cf7-9e28-b76390d961c9
md"""
# References

[1] T. L. Patti, A. Chantasri, L. P. García-Pintos, A. N. Jordan, and J. Dressel, Linear Feedback Stabilization of a Dispersively Monitored Qubit, Phys. Rev. A 96, 022311 (2017).

"""

# ╔═╡ 133b6939-10b2-4c8e-acf8-5658ca96a0f9
md" # Utilities"

# ╔═╡ 8e85754f-d66b-477b-8153-b162519edb7c
expects(ops) = ρ -> collect(real(expect(ρ, s)) for s in vcat(ops, ρ)) # ρ -> [<x>,<y>,<z>,<ρ>]

# ╔═╡ fdcfcea3-e2b6-4939-ac6c-eada7421f3dd
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

# ╔═╡ ccbcf668-d948-4ec6-a5f7-39a178d54c29
xyz(θ, ϕ) = (sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ))

# ╔═╡ 2b35485b-155b-4fb1-868d-b431cc867d61
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

# ╔═╡ 09bcbe71-e239-4a7c-ac14-99a487c1f9a4
let
	
	# System parameters ------------------------------------------------------------
	# all times given in μs
	
	# initial state 
	x0, y0, z0 = xyz(π/2, 0)
	ρ0 = DenseOperator(0.5*(id + x0*σx + y0*σy + z0*σz))
	
	# measurement parameters
	τm = 1 					# time
	Γm = 1/(2τm) 			# rate
	η =  1 					# efficiency
	
	# simulation timescales
	T = (0, 8τm) # simulation duration
	dt = 0.5e-3  # integration time-step
	
	# Hamiltonian parameters
	ΩR = rabi ? 2π : 0
	
	
	# Kraus operators -------------------------------------------------------------
	
	H(t, r) = (ΔS/2) * r[1] * σz + (ΩR / 2) * σy
	J = [(σz, ((1-η)*Γm))]
	C = [(im * σz, τm, η)]
	

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
	ρ0 = DenseOperator(0.5*(id + x0*σx + y0*σy + z0*σz))
	
	# measurement parameters
	τm = 2 					# time
	Γm = 1/(2τm) 			# rate
	η =  0.6					# efficiency
	td = 0.2
	
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
	
	H(t, r) = (Δ0 + Δ1*r[1]) * σϕ/2 + ΔS * r[2] * σz/2 
	J = idealDQM ? [(σz, ((1-η)*Γm))] : [(σz, ((1-η)*Γm + Γ2)), (σm, Γ1)]
	C = [(σz, τm, η/2), (im * σz, τm, η/2)]
	

	# Bayesian simulation ---------------------------------------------------------
	Random.seed!(seed2)
	sol = bayesian(T, ρ0, H, J, C; dt=dt, td=td)

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
	ρ0 = DenseOperator(0.5*(id + x0*σx + y0*σy + z0*σz))
	
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

# ╔═╡ f48235c9-2261-4f37-9d78-62e181fa4d43
let
	
	# System parameters ------------------------------------------------------------
	# all times given in μs
	
	# initial state 
	x0, y0, z0 = xyz(π/2, π/2)
	ρ0 = DenseOperator(0.5*(id + x0*σx + y0*σy + z0*σz))
	
	# measurement parameters
	τm = 1 					# time
	Γm = 1/(2τm) 			# rate
	η = 1 					# efficiency
	
	# simulation timescales
	T = (0, 15τm) # simulation duration
	dt = 0.5e-3  # integration time-step
	
	# Hamiltonian parameters
	θs = 3π/10 # target angle on Bloch sphere
	ϕ = π # fixes plane of oscillations
	vecDQM = (sin(θs) * cos(π/2), sin(θs) * sin(π/2), cos(θs))
	σϕ = cos(ϕ)*σx + sin(ϕ)*σy
	
	Δ0 = -sin(2θs)/(4τm)
	Δ1 = sin(θs)/τm
	global ΔSs = 0.5:0.025:1
	
	
	# Scan over ΔS ---------------------------------------------------------------
	
	global ps_ΔS = []
	
	for ΔS in ΔSs

		# Kraus operators - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		H(t, r) = (Δ0 + Δ1*r[1]) * σϕ/2 + ΔS * r[2] * σz/2 
		J = [(σz, ((1-η)*Γm))]
		C = [(σz, τm, η/2), (im * σz, τm, η/2)]

		
		# loop over N simulations - - - - - - - - - - - - - - - - - - - - - - - - -
		
		xs, ys, zs = ([], [], [])
		
		for i in 1:NΔS
			sol = bayesian(T, ρ0, H, J, C; dt=dt)
			tr = traj(sol)
			push!(xs, last(tr.x))
			push!(ys, last(tr.y))
			push!(zs, last(tr.z))
		end

		x, y, z = mean.((xs, ys, zs))
		p = 0.5 * (1 + x^2 + y^2 + z^2) # state purity
		push!(ps_ΔS, p)
	end


	
	md" ###### 🔻 Simulation "
	
end

# ╔═╡ e1fbdf59-c77f-4975-8a11-1bddc0a5b14d
begin
	myplot(ΔSs, ps_ΔS; add_marker=true, xguide="ΔS", yguide="purity", title="purity of ensemble average (N = 100) vs. ΔS", xticks=0:0.1:1)
end

# ╔═╡ 618168dc-53dd-4562-8061-67a0b56587aa
φdict = Dict("0" => 0, 
				"π/8" => π/8, 
				"π/4" => π/4,
				"3π/8" => 3π/8,
				"π/2" => π/2)

# ╔═╡ 49065a29-93ac-458d-ac0c-84f2f2151f9d
let
	
	# System parameters ------------------------------------------------------------
	# all times given in μs
	
	# initial state 
	x0, y0, z0 = xyz(π/2, 0)
	ρ0 = DenseOperator(0.5*(id + x0*σx + y0*σy + z0*σz))
	
	# measurement parameters
	τm = 1 					# time
	Γm = 1/(2τm) 			# rate
	η =  1 					# efficiency
	φ = get(φdict, φIB, 0) 	# angle
	
	# simulation timescales
	T = (0, 8τm) # simulation duration
	dt = 0.5e-3  # integration time-step
	
	
	# Kraus operators -------------------------------------------------------------
	
	H = 0 * id
	J = [(σz, ((1-η)*Γm))]
	C = [(exp(im * φ) * σz, τm, η)]
	
	
	
	# Bayesian simulation ---------------------------------------------------------
	go
	sol = bayesian(T, ρ0, H, J, C; dt=dt)
	
	global IB = traj(sol)
	
	md" ###### 🔻 Simulation "
	
end

# ╔═╡ f9eaaf70-4e0f-4503-8a64-2380682354ce
let 
	ϕv = ϕvIB * (π/16)
	tt, xx, yy, zz = (IB.t, IB.x, IB.y, IB.z)
	
	if animIB_bs
		
		anim = @animate for i ∈ range(1, length(tt), step=100)
			blochsphere(xx[1:i], yy[1:i], zz[1:i], linewidth=1., linealpha=0.85, ax=true, viewϕ = ϕv, blochmark=true) end
		
		gif(anim, fps = 15)
		
	else
		blochsphere(xx, yy, zz, linewidth=1., linealpha=0.85, ax=true, viewϕ = ϕv, blochmark=true)
		
	end
	
end

# ╔═╡ 6bc3e01f-b3cb-4a32-ab6b-e5dcc967b07f
let 
	tt, xx, yy, zz = (IB.t, IB.x, IB.y, IB.z)
	
	if animIB
		anim = @animate for i ∈ range(1, length(tt), step=100)
			blochprojections(xx[1:i], yy[1:i], zz[1:i], blochmark=true) end
		gif(anim, fps = 15)
	else
		blochprojections(xx, yy, zz, blochmark=true)
	end
end

# ╔═╡ df97d34b-16a7-49a0-a143-939f18248f48
blochtimeseries(IB.t, IB.x, IB.y, IB.z, title = "Monitored qubit", tf=last(IB.t), ylims=[-1.1,1.5])

# ╔═╡ 62472483-7fae-4adc-976c-9275e9d5ebfc
myplot(IB.t, coarsegrain(IB.r, n=50))

# ╔═╡ 13e19d59-5132-45f3-bac1-3a3b3a7a12b2
let
	
	# System parameters ------------------------------------------------------------
	# all times given in μs
	
	# initial state 
	(x0,y0,z0) = (0., 0.3, 0.91)
	ρ0 = DenseOperator(0.5*(id + x0*σx + y0*σy + z0*σz))
	
	
	# measurement parameters
	
	τm = parse(Float64, τmstr)			# time
	Γm = 1/(2τm) 						# rate
	η =  parse(Float64, ηstr)			# efficiency
	φ = get(φdict, φstr, 0) 			# angle
	td = parse(Float64, Tdstr) * 1e-3 	# time delay for feedback

	
	# simulation timescales
	
	T = (0, 8τm) # simulation duration
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
	
	H = Hf(t, r) = (Δ0 + Δ1*r[1]) * σϕ/2
	J = ideal ? [(σz, ((1-η)*Γm))] : [(σz, ((1-η)*Γm)), (σm, Γ1), (σz, Γ2)]
	C = [(exp(im * φ) * σz, τm, η)]
	
	
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



# ╔═╡ 6affc985-2426-42ed-a7be-4a18d21ccf27
let
	
	# System parameters ------------------------------------------------------------
	# all times given in μs
	
	# initial state 
	(x0,y0,z0) = (0., 0.3, 0.91)
	ρ0 = DenseOperator(0.5*(id + x0*σx + y0*σy + z0*σz))
	
	
	# measurement parameters
	
	τm = parse(Float64, τmstr)			# time
	Γm = 1/(2τm) 						# rate
	η =  parse(Float64, ηstr)			# efficiency
	φ = get(φdict, φIB, 0) 				# angle
	td = parse(Float64, Tdstr) * 1e-3 	# time delay for feedback

	
	# simulation timescales
	
	T = (0, 8τm) # simulation duration
	dt = 0.5e-3  # integration time-step
	
	
	# feedback target parameters
	
	θs = 3π/10 		# target angle on Bloch sphere
	Rs = 1 # 0.64 			# radius of target
	ϕ = π 					# fixes plane of oscillations
	vec = (sin(θs) * cos(π/2), sin(θs) * sin(π/2), cos(θs))
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
	
	H = Hf(t, r) = (Δ0 + Δ1*r[1]) * σϕ/2
	J = ideal ? [(σz, ((1-η)*Γm))] : [(σz, ((1-η)*Γm)), (σm, Γ1), (σz, Γ2)]
	C = [(exp(im * φ) * σz, τm, η)]
	
	
	# Bayesian simulation ---------------------------------------------------------
	trajs = []
	
	# loop over N simulations
	for i in 1:N
		sol = bayesian(T, ρ0, H, J, C; dt=dt, td=td)
		tr = traj(sol)
		push!(trajs, tr)
	end
	
# 	xm, ym, zm = (mean(xs), mean(ys), mean(zs))
	
# 	p = 0.5 .* (1 .+ xm.^2 + ym.^2 + zm.^2) # state purity
	t = collect(T[1]:dt:T[2])   # list of times
	
	global dfx = DataFrame([tr.x for tr in trajs], [string("x", i) for i in 1:N])
	global dfy = DataFrame([tr.y for tr in trajs], [string("y", i) for i in 1:N])
	global dfz = DataFrame([tr.z for tr in trajs], [string("z", i) for i in 1:N])
	global dfr = DataFrame([tr.r[1] for tr in trajs], [string("r", i) for i in 1:N])
	global dft = DataFrame([t], ["t"])
	
	# global ILFens = traj(t, xm, ym, zm, p, ())
	md" ###### 🔻 Simulation "
	
end

# ╔═╡ 682c16f9-80e9-40a3-96c9-085876d62ed6
CSV.write("xdata.csv", dfx)

# ╔═╡ 42cfa26c-eb57-4014-897b-9fab1d991db6
begin
	CSV.write("data/ydata.csv", dfy)
	CSV.write("data/zdata.csv", dfz)
	CSV.write("data/rdata.csv", dfr)
	CSV.write("data/tdata.csv", dfr)
end

# ╔═╡ 11a4a023-0cc6-4531-96d7-0764f226334b
let
	
	# System parameters ------------------------------------------------------------
	# all times given in μs
	

	
	# initial state 
	(x0,y0,z0) = (0., 0.3, 0.91)
	ρ0 = DenseOperator(0.5*(id + x0*σx + y0*σy + z0*σz))
	
	
	# measurement parameters
	
	τm = parse(Float64, τmstr)			# time
	Γm = 1/(2τm) 						# rate
	global ηs =  0.1:0.1:1.0 			# efficiency
	φ = get(φdict, φIB, 0) 				# angle
	td = parse(Float64, Tdstr) * 1e-3 	# time delay for feedback

	
	# simulation timescales
	
	T = (0, 8τm) # simulation duration
	dt = 0.5e-3  # integration time-step
	
	
	# feedback target parameters
	
	θs = 3π/10 		# target angle on Bloch sphere
	Rs = 1 # 0.64 			# radius of target
	ϕ = π 					# fixes plane of oscillations
	vec = (sin(θs) * cos(π/2), sin(θs) * sin(π/2), cos(θs))
	σϕ = cos(ϕ)*σx + sin(ϕ)*σy
	
	
	# feedback drive parameters
	
	Δ0 = ideal ? -sin(2θs)/(4τm) : -sin(2θs)/(4τm*Rs^2) 
	Δ1 = ideal ? sin(θs)/τm : sin(θs)/(τm*Rs)
	
	
	# decay times and rates
	
	T1 = 40 	# energy decay time
	Γ1 = 1/(2T1)# energy decay rate
	T2 = 60 	# environmental dephasing time
	Γ2 = 1/T2 	# environemntal dephasing rate
	
	
	# loop over η values -----------------------------------------------------------
	global ps = []
	
	for η in ηs

		# Kraus operators - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

		H(t, r) = (Δ0 + Δ1*r[1]) * σϕ/2
		J = ideal ? [(σz, ((1-η)*Γm))] : [(σz, ((1-η)*Γm)), (σm, Γ1), (σz, Γ2)]
		C = [(exp(im * φ) * σz, τm, η)]


		# Bayesian simulation - - - - - - - - - - - - - - - - - - - - - - - - - - - 
		
		xs, ys, zs = ([], [], [])

		# loop over N simulations
		for i in 1:N_ηs
			sol = bayesian(T, ρ0, H, J, C; dt=dt, td=td)
			tr = traj(sol)
			push!(xs, last(tr.x))
			push!(ys, last(tr.y))
			push!(zs, last(tr.z))
		end

		x, y, z = mean.((xs, ys, zs))
		p = 0.5 * (1 + x^2 + y^2 + z^2) # state purity
		push!(ps, p)

	end
	
	md" ###### 🔻 Simulation "
	
end

# ╔═╡ 5de35dba-fb82-49a4-9eec-c92750e7a4e9
begin
	myplot(ηs, ps; add_marker=true, xguide="η", yguide="purity", title="purity of ensemble average (N = 100) vs. η")
end

# ╔═╡ Cell order:
# ╠═4c45fe3e-cd69-11eb-20cd-7bfb98c040cf
# ╟─377a3336-20bd-4baa-a033-af8bbc8668a8
# ╟─3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─8f7c6440-eac8-49d3-bd23-c2545ec16830
# ╟─ff721d9a-c37c-43a1-adee-015221ddd133
# ╟─727ced31-7a51-4892-b1d9-eff4fcdc9767
# ╟─1aaa9cbe-2f98-41e1-af41-598ba6578818
# ╟─1ee143ba-2f0d-4eee-a76c-4bf58469e0be
# ╟─d3e48ac4-973f-4397-bf06-89150a0c13ed
# ╟─e43e5329-bd96-41ce-a183-1bd206204f65
# ╠═49065a29-93ac-458d-ac0c-84f2f2151f9d
# ╟─875adb1b-22f7-4375-83c5-2a17c9a8e46c
# ╟─3531558a-f9b7-4f6b-9800-884fe0b04712
# ╠═f9eaaf70-4e0f-4503-8a64-2380682354ce
# ╟─faf0b11c-4339-4a82-a23c-9d35eb9d10b4
# ╟─6bc3e01f-b3cb-4a32-ab6b-e5dcc967b07f
# ╠═df97d34b-16a7-49a0-a143-939f18248f48
# ╠═62472483-7fae-4adc-976c-9275e9d5ebfc
# ╟─3ff551a9-3a07-4a2f-928e-880b7e3ba1fc
# ╟─911d4aa1-41cd-4eea-8e0f-1af6f8024290
# ╟─3ae402b0-476a-4c11-823c-d0016793adab
# ╠═71604112-0704-4832-8a80-e886b6155e34
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
# ╠═c9d3c54e-acff-4a84-8f95-937ee1602350
# ╠═5e76a75a-3041-49f6-b5b6-96047b1d5bd4
# ╠═13e19d59-5132-45f3-bac1-3a3b3a7a12b2
# ╟─e486fbed-157c-4139-914c-ada6bb88d7b4
# ╟─4898ab97-4058-4e70-a959-2962641d9611
# ╟─6b93bb84-0ba1-44e8-910e-612793b3df3b
# ╟─ce79057b-eb31-49b1-9db7-bfd793a447c5
# ╟─ca7a2351-cff5-4d77-ba16-00de384b8b7c
# ╟─f3e7794e-a9a1-4103-8469-32a8e37d2d82
# ╟─a12cdb8c-e9a1-4c2d-9811-cff266e152d8
# ╠═725dc4c3-cc74-4400-819c-2cffd06fbbf9
# ╟─0a7f28c9-1d84-43e6-b62c-711a231a3972
# ╟─34a700bb-5809-4755-a7fa-def102c5fd4c
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
# ╟─219487b2-c85c-4b64-8387-5d8a90801dc6
# ╠═98baaa59-213b-4f6e-929b-db10ad4fd718
# ╠═6affc985-2426-42ed-a7be-4a18d21ccf27
# ╟─d7e4182d-46ba-40f5-97bd-ced14022b3f0
# ╠═88b77211-fe49-4afb-9c2a-c64791c20395
# ╠═682c16f9-80e9-40a3-96c9-085876d62ed6
# ╠═42cfa26c-eb57-4014-897b-9fab1d991db6
# ╠═9638226c-7c76-4230-9cc9-cc6700232e20
# ╠═2ef48cde-d72e-4d91-8de3-5fcb52d764eb
# ╟─a8b375ae-5203-4c4f-96c8-e594c5eb098d
# ╠═e8a25893-452f-461f-9050-ee16b6f80f54
# ╠═9279175b-0bb8-4f4d-b1ed-6289aacdaacf
# ╟─849589c4-d260-4f24-a1c5-cf4887ba9701
# ╠═84275f2c-a215-4d80-ba81-326ed0eab8c3
# ╠═11a4a023-0cc6-4531-96d7-0764f226334b
# ╟─5de35dba-fb82-49a4-9eec-c92750e7a4e9
# ╟─e15a6a43-31fe-488d-aaab-aff9f50d1b8a
# ╟─6a5deff1-d8ff-47ea-9efa-ecc3b8c6dec0
# ╟─fe7af10f-ab48-4a49-9bb1-e368ec94d928
# ╟─ebb1ef94-f30a-4203-8862-5d0294588d30
# ╟─09bcbe71-e239-4a7c-ac14-99a487c1f9a4
# ╠═39679392-2d6a-42f1-89d7-4bd20a953090
# ╟─afb6751a-31af-44cf-8922-0f1b816755ca
# ╟─c337e3e0-95cf-4c18-987e-59216bf54419
# ╟─85a3d2ae-f740-43e6-80e3-9ba76c91fb68
# ╟─16ccad7b-1f63-4373-ab81-3f77484a3545
# ╟─3c0e6012-0a8b-4975-baf7-bfe00571933f
# ╟─175af755-b2e2-4a3e-9789-855bdec353ae
# ╟─44cfbce3-0f4e-4641-90b8-f4853ccd68ea
# ╟─c414534a-01d8-422b-9286-ce218a39bee8
# ╟─88ec8109-87ed-457d-b9ab-b374176150b1
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
# ╟─723b3eab-d717-4823-a7eb-da39202f6f7b
# ╠═1f017209-c017-4d77-9a94-8edde2c02fdb
# ╠═f48235c9-2261-4f37-9d78-62e181fa4d43
# ╠═e1fbdf59-c77f-4975-8a11-1bddc0a5b14d
# ╟─9e98e449-859e-4cf7-9e28-b76390d961c9
# ╟─133b6939-10b2-4c8e-acf8-5658ca96a0f9
# ╟─8e85754f-d66b-477b-8153-b162519edb7c
# ╠═fdcfcea3-e2b6-4939-ac6c-eada7421f3dd
# ╠═ccbcf668-d948-4ec6-a5f7-39a178d54c29
# ╠═2b35485b-155b-4fb1-868d-b431cc867d61
# ╠═618168dc-53dd-4562-8061-67a0b56587aa
