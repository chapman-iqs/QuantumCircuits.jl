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

# ╔═╡ 1aa3dd28-34e8-4ebe-bc07-8ec41527bb1d
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
In this interactive notebook, we'll look at properties and dynamics of coherent states, which are key to understanding superconducting qubit readout. This notebook is based on unpublished notes ($Dressel_notes_2019).
"""

# ╔═╡ 13f42433-1bc7-4f5e-b97c-a4a2664e1a71
mdp(table_of_contents📔)

# ╔═╡ 3dfd0a4f-ed67-4c19-9fdb-98b3b96c62d1
TableOfContents(title="Coherent state dynamics")

# ╔═╡ 1c6ac0e4-690d-4eff-9a88-68883abc5f67
md"""
# Problem setup
"""

# ╔═╡ e4ffe1eb-c206-4ea9-9536-8942d2ae7e36
md"""
### Qubit-resonator interaction

A qubit with bare frequency $\omega_q \sim 4-9$ GHz is dispersively coupled to a microwave resonator with bare frequency $\omega_r$ detuned by $\Delta \equiv \omega_q - \omega_r \sim 1$ GHz from the qubit. The resonator frequency shifts by $2\chi$ depending on the qubit state, with $\chi \sim 1$ MHz, i.e. the resonator has dressed frequency $(\omega_r)^{\pm} = \omega_r \pm \chi$, depending on whether the qubit is excited ($+$) or ground ($-$). The qubit is also coupled to a (time-dependent) Rabi drive with frequency $\Omega_R(t) \ll \omega_q$. The Hamiltonian describing this interaction is

```math
\hat H_{qr} / \hbar = \omega_q \hat \sigma_z/2 + (\omega_r + \chi \sigma_z) \hat a^\dagger \hat a+ \Omega_R(t) \hat \sigma_y /2.
```
"""

# ╔═╡ db847000-8785-4bcb-9f93-db82f6ead06d
md"""
### Resonator-transmission line interaction

The resonator is also coupled to a transmission line with traveling input and output fields $\hat c_\text{in}$ and $\hat c_\text{out}$. The resonator has energy leakage at rate $\kappa \approx \omega_r (Z_r / Z_t)$, where $Z_r / Z_t$ is the ratio of impedances of the resonator and transmission line. This interaction is described by

```math
\hat H_{rt} / \hbar = i (\sqrt{\kappa}/2) (\hat a + \hat a^\dagger) (\hat c_\text{in} - \hat  c_\text{in}^\dagger - \hat c_\text{out} + \hat c_\text{out}^\dagger)
```

with boundary condition

```math
\sqrt{\kappa} \hat a = \hat c_\text{in} + \hat c_\text{out}.
```

"""

# ╔═╡ b7be6342-01f9-4afd-9250-27f0479ed5c5
md"""
The resonator field therefore evolves according to

```math
\dot{\hat{a}} = \frac{i}{\hbar} [\hat H_{qr} + \hat H_{rt}, \hat a]

=-i(\omega_r + \chi \hat \sigma_z) \hat a + (\sqrt{\kappa}/2)(\hat c_\text{in} - \hat  c_\text{in}^\dagger - \hat c_\text{out} + \hat c_\text{out}^\dagger).
```

With the boundary condition this becomes

```math
\dot{\hat{a}} = -i(\omega_r + \chi \hat \sigma_z) \hat a + \sqrt{\kappa} (\hat c_\text{in} - \hat  c_\text{in}^\dagger) - \frac{\kappa}2 (\hat a - \hat a^\dagger)
```
"""

# ╔═╡ 306f55d1-fd48-4830-ba0a-b54ca5dccf82
md"""
### Rotating-wave and coherent-field approximations
"""

# ╔═╡ 56f69b9a-af7e-4e50-bf6a-0d4e436a7f95
md"""
We assume:
* coherent input fields, such that $\langle \hat c_\text{in} \rangle = -i \epsilon(t) e^{-i \omega_d t}/\sqrt{\kappa}$ with $\omega_d \sim \omega_r$.

* slowly varying drive envelope $\epsilon(t)/\sqrt{\kappa}$ compared to drive frequency $\omega_d$

* resonator field adiabatically follows qubit state (slow enough $\Omega_R$)

* resonator state corresponding to definite qubit state ($\ket{0}$ / $\ket{1}$) is the coherent state $e^{-i \phi_\pm(t)} \ket{\alpha_\pm(t) e^{-i \omega_d t}}$, so that the cavity expectation values are $\braket{\hat a}_\pm = \alpha_\pm(t) e^{-i \omega_d t}$ when $\braket{\hat \sigma_z} = \pm 1$.

Then the resonator-transmission line interaction in the rotating frame simplifies to

$$\hat H_{rt}/\hbar = \hat a^\dagger \epsilon(t) e^{-i \omega_d t} + \hat a \epsilon^*(t) e^{i \omega_d t}$$

and the equation for $\dot{\hat{a}}$ can be expressed simply in terms of the expected coherent state amplitudes

$$\dot{\alpha}_\pm(t) = -i(\omega_r - \omega_d \pm \chi) \alpha_\pm(t) - i \epsilon(t) - \frac{\kappa}2 \alpha_\pm(t),$$

again in the rotating frame.
"""

# ╔═╡ c681b1d2-a32b-4c10-bc80-6447368f13eb
md"""
# Characterizing the resonator
"""

# ╔═╡ d5b3318e-a1a8-48c2-baf7-92b398902aff
md"""
## Steady-state characteristics
"""

# ╔═╡ ce3f9a90-70ed-4d73-b7ac-e61791eb1aa9
md"""
In the case where $\epsilon(t) = \epsilon$ is constant, setting the differential equation for $\dot{\alpha}_\pm(t) = 0$ immediately yields steady state amplitudes of the resonator:

$$\alpha_\pm^{(s.s.)} \equiv \frac{2 \epsilon}{\kappa} \frac{-i}{1 + i[2(\Delta \pm \chi)/\kappa]},$$

where $\Delta \equiv \omega_r - \omega_d$ is the drive detuning from the bare resonator frequency. The steady-state amplitudes allow us to characterize the cavity (and the measurement) in several important ways.

"""

# ╔═╡ 543b95ab-589b-4523-9bf5-b955d249699e
begin
	ε0 = 2π * 0.4 # MHz
	κ = 7.5 # MHz
	χ = 0.25 # MHz
	Δs = range(-10, 10, step=0.2) # MHz
	md" 🌀 parameters"
end

# ╔═╡ 9e816960-3731-4b98-96c8-9aa689bdf1ac
begin
	αp(Δ) = (2ε0/κ) * -im / (1 + im * (2(Δ + χ)/κ))
	αm(Δ) = (2ε0/κ) * -im / (1 + im * (2(Δ - χ)/κ))
	md" 🌀 function definitions"
end

# ╔═╡ 83f4d601-33f2-4c3d-92c4-37231d55cae3
md"""
### State distinguishability
"""

# ╔═╡ 54e9278d-7757-4af9-968b-c78dfbe04758
md"""
Superconducting qubit readout depends on the distinguishability of cavity coherent states conditioned on the state of the qubit. Since resonator amplitudes are complex numbers, 

$\alpha_\pm = |\alpha_\pm| e^{i \phi_\pm},$


they can be distinguished in amplitude or in phase. To better visualize this, they can be plotted in the complex plane via their real and imaginary parts (called "quadratures"). The quadratures are the experimentally observable quantities used to infer the amplitude and phase.
"""

# ╔═╡ 9c860ffd-7a2d-45df-84a8-1ddcf0f46298
md"""

 $\Delta$ = -4 MHz
$(@bind Δ Slider(-4:0.1:4, default=0.1))
4 MHz

"""



# ╔═╡ 897d6364-90ee-41cc-ab70-d5ef0c537bd2
md"""
Currently, $\Delta$ = $Δ MHz
"""

# ╔═╡ e73546a2-2b9e-4062-89f2-5a0913209988
md"""
### Adiabatic resonator evolution
"""

# ╔═╡ 05def7cc-6c0e-496e-9429-1fd457575b31
md"""
So long as we assume $\epsilon(t)$ changes very slowly (we'll discuss what is "slow enough" in the following section), we can assume the resonator remains in steady-state as the drive amplitude changes. Let's look at how the resonator states evolve with the drive amplitude. As an example, we consider a drive amplitude with sinusoidal ramp-up / ramp-down.
"""

# ╔═╡ 38702e03-a6ff-4c7a-bbb6-4dd0bcf89ba8
begin
	ΩR0 = 2π * 0.4 # rad MHz
	Δ0 = 0
	dt = 0.01
	md" 🌀 parameters"
end

# ╔═╡ 00257bca-f1bd-4c4b-af0b-db0bafd5e406
md" $(@bind anim_plot1 CheckBox()) Animate cavity population "

# ╔═╡ f5719cfb-1acf-4f87-b12e-a0ed65f635d5
md"""

Δ : -4 MHz
$(@bind Δ1 Slider(-4:0.2:4, default=0.0))
4 MHz

"""

# ╔═╡ ca2e2708-789b-4a95-ba3f-0402f93cb84e
if Δ1 == 0 
	md"""
	The cavity drive detuning is currently $\Delta$ = $Δ1 MHz. This means the cavity is being driven on resonance. Adjust the detuning below to see what changes.
	"""
else
	md" The cavity drive detuning is currently $\Delta$ = $Δ1 MHz. The cavity is being driven off-resonance."
end

# ╔═╡ cf3fd908-d968-4581-a8b4-5aac46c4afb0
md"""
## Non-adiabatic evolution
"""

# ╔═╡ c5141459-60d0-47ed-90f3-ac4e741e5707
md"""
Since we have assumed that $\alpha_+, \alpha_-$ couple independently to the qubit states $\ket{e}, \ket{g}$, we can consider their evolution given a time-dependent drive without reference to qubit evolution. In the following, we solve the differential equations in $\alpha_\pm$,

$$\dot{\alpha}_\pm(t) = -i(\Delta \pm \chi) \alpha_\pm(t) - i \epsilon(t) - \frac{\kappa}2 \alpha_\pm(t),$$

using the $4^{th}$ order Runge-Kutta integration method:


"""

# ╔═╡ 71ea6e65-46b3-4c6c-a919-87dd31a6674f
begin
	α̇p(t, αp, ε) = -im * (Δ1 + χ) * αp - im * ε(t) - (κ/2) * αp
	α̇m(t, αm, ε) = -im * (Δ1 - χ) * αm - im * ε(t) - (κ/2) * αm
	md" 🌀 function definitions"
end

# ╔═╡ 83e516e3-20dd-4cbc-9f6a-84e1662d0dcb
md" $(@bind anim_plot2 CheckBox()) Animate cavity population "

# ╔═╡ a9dc195d-fe8f-4c6b-8078-597585b2ee3f
md"""

ramp time = $((1/0.01)/2) μs
$(@bind fε html"<input type=range min=0.01 max=0.12 step=0.01 value=0.01>")
$(round((1/0.12)/2, digits=3)) μs

"""

# ╔═╡ 3ca15002-9b64-416c-94e3-a7904766a87b
md"""
Recall our second assumption, namely that the drive envelope $\epsilon(t)/\sqrt{\kappa}$ varies slowly compared to the drive frequency $\omega_d$. We can see already that there is a visible discrepancy between the steady state values $\alpha_+^{s.s.}, \alpha_-^{s.s.}$ and the time-dependent $\alpha_+, \alpha_-$. Currently, the ramp is $(round((1/fε)/2, digits=3)) μs long.

We'll can see an even greater discrepancy if we ramp-up the drive faster. Try changing the ramp time below:
"""

# ╔═╡ 14cccb18-0187-4c91-ba53-b920ce0c3704
md"""
In real qubit readout, the adiabatic approximation is clearly inadequate since pulses can last just 20 ns, with $\epsilon(t)$ at much larger values than we have considered.
"""

# ╔═╡ e8299f27-e74b-47d1-8349-7e28c3ea4495
md"""
# Reduced-qubit dynamics
"""

# ╔═╡ adb901f4-b296-4549-a92f-b075b40e8f8f
md"""
The time-dependency of cavity states has several effects on the qubit evolution. 
"""

# ╔═╡ 69055a05-e2f1-4c8a-ad2b-1ca8ac0dadc8
md"""
##### Qubit rates
"""

# ╔═╡ e4269abf-9b63-4fbe-8001-397a449dd0e6
md"""
The time-evolving cavity photon number weak value and state separation

$n_w(t) \equiv \alpha_+^*(t) \alpha_-(t)$

$\Delta \alpha(t) \equiv \alpha_+(t) - \alpha_-(t)$

lead to a time-dependent AC-Stark shift and time-dependent dephasing rate

$\Delta \tilde \omega_q(t) \equiv 2 \chi \text{Re}(n_w(t)) - \kappa \eta \text{Im}(n_w(t)) - \eta \text{Re}(\epsilon^*\Delta \alpha(t))$

$\tilde{\Gamma}(t) \equiv 2 \chi \text{Im} (n_w(t)) - \frac{\kappa \eta}2 |\Delta \alpha(t)|^2.$


"""

# ╔═╡ 58e7e646-eaf8-4ba0-8204-9819a054fe65
md"""
##### Unitary evolution, decoherence, and measurement
"""

# ╔═╡ 02d2a1f6-12af-4345-b42a-db2b11e2a883
md"""
The AC-Stark shift leads to a time-dependent qubit Hamiltonian in the rotating frame

$\hat H_q(t) = \frac{\Delta \tilde \omega_q(t)}2 \hat \sigma_z$

which gives the discrete time-evolution operator

$\hat U_{\Delta t}(t) = \exp(-i \Delta t \hat H_q(t)).$

"""

# ╔═╡ 498bfa13-4278-4808-b1d6-279eabcd9ed8
md"""
The time-dependent dephasing rate leads to decoherence

$\rho_{01}^{(q)} \mapsto \rho_{01}^{(q)} e^{-i \Delta t \tilde \Gamma}.$

In the `QuantumCircuits.jl` `bayesian` or `rouchon` implementations, this is expressed using a Lindblad operator  $\hat \sigma_z$ with rate $\tilde \Gamma(t)$. (Note that since `bayesian` and `rouchon` are actually CP maps and do not involve true Lindbladians, they can handle negativity of $\tilde \Gamma(t)$. This can be interpretted as a non-Markovian "re-coherence", and is a special feature of the non-adiabatic qubit-resonator evolution.)

Finally, measurement of the random complex signals $r_\varphi$ and $r_{\varphi + \pi/2}$ leads to measurement backaction, expressed through the measurement operator
"""

# ╔═╡ 9423dbb2-c325-415a-aeca-31cc8b205173
html"""<a class="anchor" id="measurement_operator"></a>"""

# ╔═╡ 50ac7363-0591-49dc-9bdf-dc37dcf43f60
md"""
$\hat M_{r_\theta}(t) \equiv \exp \Bigg[\frac{r_\theta \Delta t}{2 \tau_m} e^{i(\theta + \text{angle}(\Delta \alpha))}\hat \sigma_z \Bigg] = \Bigg[\frac{r_\theta \Delta t}{2 \tau_m} e^{i \theta} \hat \sigma_z \Bigg]$
"""

# ╔═╡ c7eaf3d6-73ac-4291-ab9c-5857cb77cdaf
md"""
(The second equality follows from our choosing $\text{angle} (\Delta \alpha) = 0$, which is an experimental freedom to choose the definition of quadratures.) This is equivalent to measuring $e^{i \theta} \hat \sigma_z$ with measurement collapse time $\tau_m(t) = 1/2\Gamma_m(t)\eta$ and efficiency $\eta$, where $\Gamma_m(t)$ is the ensemble measurement dephasing rate

$\Gamma_m(t) = \frac{\kappa}2 |\Delta \alpha(t)|^2.$
"""

# ╔═╡ d0528468-7e59-41a1-a81a-5d8933864f49
md"""
##### Measurement record
"""

# ╔═╡ f568802e-e754-4bad-a221-5d6102e2068f
md"""
Here, we emphasize the dual-quadrature measurement record $\tilde r = r_\varphi + i r_{\varphi + \pi/2}$, where

$r_\theta(t) = \Big\langle{\frac{\hat m_\theta(t) + \hat m_\theta(t)^\dagger}{2}} \Big\rangle + d\zeta(t)$

where $d\zeta(t) \sim \mathcal{N}(0, \tau_m/dt)$ is a normally distributed zero-mean random variable with variance $\tau_m / dt$, and 

$\hat m_\theta(t) = e^{i \theta} \hat \sigma_z + \frac{\Delta \tilde{\omega}_q(t)}{\Gamma_m(t)} \hat I$

is the measured operator. The phase on $\hat \sigma_z$ expresses the quadrature of measurement, with $\theta = 0$ being informational and $\theta = π/2$ being non-informational (often called "quadrature", meaning "in quadrature to the informational axis"). The $(\Delta \tilde{\omega}_q(t) /\Gamma_m(t)) \hat I$ term creates an overall bias in the record, and is due to the dispersive shift of the *resonator* away from its undressed frequency.

It is interesting to note that it is possible (and, in this physical scenario, inevitable) for $\hat m_\theta(t)$ to be non-Hermitian. This will sound like a red flag for those used to thinking of measured operators as corresponding to observables, which must be Hermitian. However, the non-Hermitian part of $\hat m_\theta(t)$ is not observed -- as shown in the expression for $r_\theta$, which includes on the Hermitian part $(\hat m_\theta(t) + \hat m_\theta(t)^\dagger)/2$. It only comes into play in the measurement operator.
"""

# ╔═╡ fca59554-0975-415a-8143-56b86d3de5bc
mdp("A non-Hermitian measurement operator is used to model backaction due to a non-informational measurement. See ", measurement_backaction📔, " for a description of this.")

# ╔═╡ 0c326351-6d79-4e16-9bf0-8170fd92a671
md"""
## Cavity steady state dynamics
"""

# ╔═╡ 4e368417-1107-45af-9af8-10fc16eaabe5
md"""
In the following discussion, we assume the cavity is driven on resonance ($\Delta = 0$) and look at steady state values of qubit evolution parameters. When the cavity is at steady state in the bad-cavity (Markovian) limit, this leads to several simplifications relative to the previous discussion.
"""

# ╔═╡ ebfe6948-64be-4190-8470-d5e196ee43ec
md"""
$n_w^{s.s.} = (\alpha_+^{s.s.})^* \alpha_-^{s.s.},$

produces a weak dispersive shift 

$2\chi n_w = \omega_S + i \Gamma.$

The real part of the weak photon number is itself the effective photon number in the cavity:

$\bar n = \text{Re} (n_w).$ 

"""

# ╔═╡ 23b8def8-5019-41ab-8e55-b70ee857b7a6
md"""


This determines the ensemble-averaged AC-Stark shift

$\omega_S = 2\chi \bar n \approx 2 \chi \Bigg| \frac{2 \epsilon}{\kappa} \Bigg|^2$ 

and the ensemble-averaged measurement dephasing rate 

$\Gamma_m = 2 \chi \text{Im} (n_w) \approx \frac{8 \chi^2 \bar n}{\kappa}.$

The approximations are valid when $|\kappa| \gg |\Delta|, |\chi|$ (bad-cavity regime) and when the resonator is at steady-state, and lead to an effective qubit frequency

$\tilde{\omega}_q \approx \omega_q + (1 - \eta) 2 \chi \bar{n}$

and effective dephasing rate

$\tilde \Gamma \approx (1 - \eta)\Gamma.$

"""

# ╔═╡ 67c05b7c-4376-495a-8a21-b5d4f749c0d6
md"""
$$\alpha_\pm^{(s.s.)} \equiv \frac{2 \epsilon}{\kappa} \frac{-i}{1 + i[2(\Delta \pm \chi)/\kappa]},$$


"""

# ╔═╡ 3c50b081-f1d0-40a4-994d-d1d4813b98f4
md"""
### Single trajectory simulation
"""

# ╔═╡ fc628bbf-47cf-4603-bba6-52bcb4eed93f
begin
	local ε = ΩR0 # MHz
	local Δ = 0.0
	
	ϕ = 0
	η = 0.6
	ΩR = 0.0

	ωs = 2χ * abs(2ε / κ)^2
	Γm = 8 * (χ^2) * abs(2ε / κ)^2/ κ
	
	Δω̃q(η) = (1 - η) * ωs
	Γ̃(η) = (1 - η) * Γm

	md" 🌀 Reduced qubit steady-state parameters"
end

# ╔═╡ a0b3eef2-9d43-441a-b541-e7ef430d86be
let
	ψ0 = normalize(g + e)
	dt = 1e-3  # integration time-step
	tf = 10.0

	# evolution operators -------------------
	J = [(σz, Γ̃(η))]
	C = [(exp(im * ϕ) * (σz + (ωs / Γm) * Iq), Γm, η/2), (exp(im * (ϕ + π/2)) * (σz + (ωs / Γm) * Iq), Γm, η/2)]
	# C = [(exp(im * ϕ) * σz, Γ, η/2), (exp(im * (ϕ + π/2)) * σz, Γ, η2/)]
	H = (ΩR / 2) * σy + (Δω̃q(η) / 2) * σz

	global solss = bayesian((0.0, tf), ψ0, H, J, C; dt=dt)
	
end

# ╔═╡ 6af54bde-991c-4810-b5a4-7be74f743ee3
md"""
### Check: ensemble average
"""

# ╔═╡ c76e3314-e3a4-44d3-844e-cf141b1f550a
N = 50

# ╔═╡ ab601ea6-69cb-4079-9cb7-a57ecfd15b12
let	
	ψ0 = normalize(g + e)
	dt = 1e-3  # integration time-step
	tf = 100.0

	# evolution operators -------------------
	J(η) = [(σz, Γ̃(η)/2)]
	C(η) = [(exp(im * ϕ) * σz - (ωs / Γm) * Iq, Γm, η/2), (exp(im * (ϕ + π/2)) * σz - (ωs / Γm) * Iq, Γm, η/2)]
	# C = [(exp(im * ϕ) * σz, Γ, η/2), (exp(im * (ϕ + π/2)) * σz, Γ, η2/)]
	H(η) = (ΩR / 2) * σy + (Δω̃q(η) / 2) * σz

	global ss_ensemble = map(1:N) do m
		bayesian((0.0, tf), ψ0, H(η), J(η), C(η); dt=dt)
	end

	global ss_η0 = bayesian((0.0, tf), ψ0, H(0), J(0), []; dt=dt)
	
end

# ╔═╡ 2b8f2e1c-0f25-4be0-ab40-2bb0d96ba5e8
md"""
## Non-steady state dynamics
"""

# ╔═╡ 6cb5d42f-b074-45c4-8f19-a70f04f4fad2
md"""
### In bad-cavity regime
"""

# ╔═╡ 44b8a3f8-b220-4567-a486-576908bd4255
md"""
Below, we check that the exact dynamics for setting $\epsilon$ = $(round(ΩR0, digits=3)) MHz (as used in the steady state solution above) match the steady-state solutions, given the same record. This is another check that the bad cavity approximation is valid in this regime.
"""

# ╔═╡ e13483e5-41b3-438e-9504-2b836d544c00
md"""
### In general
"""

# ╔═╡ 1fc48bb8-d254-470d-b3bc-76acc342c6b6
md"""
Notice that the informational and quadrature records have a transient behavior where they have very large variance. This occurs during the transient of $\Gamma_m(t)$ before it stabilizes to its steady state value. Because of the $1/2\tau_m$ term in the [measurement operator](#measurement_operator), the backaction during this time will actually be reduced.
"""

# ╔═╡ 7e62e6c6-aa3c-4350-901c-15d017b8db42
md" # Utilities "

# ╔═╡ 04f6e40d-2dfa-46fd-8e43-58c3f647b5b5
md"""
## Plotting
"""

# ╔═╡ 62b5300d-32c5-477f-8150-89d90b5ae31c
function cavity_plot(ts, εt, (αp_list, αm_list); xlims=[-1.5,1.5], ylims=[-1.5,1.5])
	return cavity_plot(length(ts), ts, εt, (αp_list, αm_list); xlims=xlims, ylims=ylims)
end

# ╔═╡ 6038bb81-8e47-4be6-a5be-715101349309
function cavity_plot(ts, εt, (αp_list, αm_list), (αp_list2, αm_list2); xlims=[-1.5,1.5], ylims=[-1.5,1.5])
	return cavity_plot(length(ts), ts, εt, (αp_list, αm_list), (αp_list2, αm_list2); xlims=xlims, ylims=ylims)
end

# ╔═╡ 29b15810-41e0-4e2d-b676-7e7378cd261a
function qubit_plot(sol::Solution; record=false, title="", legendpos=:bottomleft, rlims=:default, rlabels=["r1", "r2"])

	basis = qbasis

	t = sol.t
	exps = map(op -> expectations(sol, op), basis)
	rs = record ? sol.r : []

	return qubit_plot((t, exps, rs); title=title, legendpos=legendpos, rlims=rlims, rlabels=rlabels)
end

# ╔═╡ 7af4d4ff-30ad-4ec8-9aea-ca5c2b431f57
function bloch_plots(sols::Vector{Solution}, sol_η0::Solution; alpha=0.1, N=50)
	colors = colors1q 
	
	# calculate expectation values --------------------------------------------
	t = sols[1].t
	xs, ys, zs = [], [], []
	
	for sol in sols
		x, y, z = map(op -> expectations(sol, op), qbasis)
		for (list, traj) in zip([xs, ys, zs], [x, y, z])
			push!(list, traj)
		end
		
	end

	# η = 0 solution
	xη0, yη0, zη0 = map(op -> expectations(sol_η0, op), qbasis)

	
	
	# plot ----------------------------------------------------------------------
	function bloch(os, oη0; color=colors1q[1], xlabel="", ylabel="")
		
		po = plot(ylims = [-1,1], xlabel=xlabel, ylabel=ylabel)
		
		for o in os[1:min(N, 50)]
			plot!(t, o, alpha=alpha, label=:none, color=color)
		end
		
		oavg = [mean([os[i][j] for i in 1:N]) for j in 1:length(t)]
		plot!(t, oavg, alpha=1, color=color, label="average", linewidth=3)
		plot!(t, oη0, alpha=1, color=:black, label="η = 0", linewidth=2)

		po
		
	end

	
	l = @layout [xplot{0.33h}; yplot{0.33h}; zplot{0.33h}]
	
	px = bloch(xs, xη0, color=colors[1], ylabel="x")
	py = bloch(ys, yη0, color=colors[2], ylabel="y")
	pz = bloch(zs, zη0, color=colors[3], ylabel="z", xlabel="t (μs)")
	
	plot(px, py, pz, layout = l, link=:y, size=(800,500), legendfontsize=8, titlefontsize=12, legend=:outerright)
	
end

# ╔═╡ c9efabea-2802-4bcd-adbf-3608155fba58
bloch_plots(ss_ensemble, ss_η0, alpha=0.15, N=N)

# ╔═╡ d88511e7-af21-4b0b-a654-348060892498
function animate_plot(tt, plot, args...; step=100, fps=15)
	anim = @animate for i ∈ range(1, length(tt), step=step)
			plot(i, args...) end
	gif(anim, fps = fps)
end	

# ╔═╡ b07e6b7b-b96e-4d0d-ac42-e3b717d5163a
begin
	colors1 = palette(:rainbow)
	colors2 = palette(:lightrainbow)
	colorsq1 = palette(:tab20)[1:2:20]
	colorsq2 = palette(:tab20)[2:2:20]

	mixed_colors = let
		cols = []
		for i in 1:length(colors1)
			push!(cols, colors1[i])
			push!(cols, colors2[i])
		end
		cols
	end
end

# ╔═╡ 67d7fdb7-145f-42f3-8dec-cb6fc8bb3018
let
	colors=colors1
	
	αp_list = αp.(Δs)
	αm_list = αm.(Δs)
	
	# Plot phase vs. detuning Δ ---------------------------------------------------
	
	p1 = plot(ylabel="amplitude", legend=:bottomright)
	plot!(Δs, abs.(αp_list ), color=colors[2], label=L"\alpha_+")
	plot!(Δs, abs.(αm_list ), color=colors[4], label=L"\alpha_-")
	plot!([Δ,Δ], [0.2,1.0], linestyle=:dash, color="red", label=:none)

	# Plot amplitude vs. detuning Δ ------------------------------------------------
	
	p2 = plot(ylabel="phase", xlabel=string(L"$\Delta$", " (MHz)"), legend=:bottomright)
	plot!(Δs, angle.(αp_list), color=colors[2], label=L"\alpha_+")
	plot!(Δs, angle.(αm_list), color=colors[4], label=L"\alpha_-")
	plot!([Δ,Δ], [-3.0,-0.2], linestyle=:dash, color="red", label=:none)

	
	# Plot α in complex plane ---------------------------------------------------
	
	p3 = plot(xlims = [-1.1, 1.1], ylims = [-1.1, 1.1], xlabel=string("Re ", L"\alpha_\pm"), ylabel=string("Im ", L"\alpha_\pm"), legend=:topright, title="coherent state wavepacket centers", titlefontsize=10, legendfontsize=10)
	plot!(real.(αp_list), imag.(αp_list), color=:gray, label=:none)
	plot!([real(αp(Δ))], [imag(αp(Δ))], color=colors[2], marker="o", label=L"\alpha_+")
	plot!([real(αm(Δ))], [imag(αm(Δ))], color=colors[4], marker="o", label=L"\alpha_-")

	l = @layout [a{0.5w, 1.0h} grid(2,1)]

	plot(p3, p1, p2, layout=l, size=(600,300))
	

	
end

# ╔═╡ 4029ecf7-7983-465f-bc8c-b8b6fb25e47e
function cavity_plot(i::Int64, ts, εt, (αp_list, αm_list); xlims=[-1.5,1.5], ylims=[-1.5,1.5])

	colors=colors1
	
	αps = αp_list[1:i]
	αms = αm_list[1:i]

	pα = plot(xlims=xlims, ylims=ylims, xlabel=string("Re ", L"\alpha_\pm"), ylabel=string("Im ", L"\alpha_\pm"), legend=:right, title="coherent state wavepacket centers", titlefontsize=12)
	
	plot!(real.(αps), imag.(αps), color=colors[2], linestyle=:dash, label=:none)
	plot!([real(last(αps))], [imag(last(αps))], color=colors[2], marker="o", label=L"\alpha_+^{s.s.}")
	plot!(real.(αms), imag.(αms), color=colors[4], linestyle=:dash, label=:none)
	plot!([real(last(αms))], [imag(last(αms))], color=colors[4], marker="o", label=L"\alpha_-^{s.s.}")

	plot!(ts, εt, color=colors[5], label="", xlabel = "t (μs)", ylabel = "ε (MHz)", legend=:none, inset = (1, bbox(0.05, 0.05, 0.4, 0.25, :top, :right)), subplot=2, title="cavity drive", titlefontsize=10, axisfontsize=10)
	plot!([ts[i]], [εt[i]], color=colors[5], marker="o", label=:none, subplot=2)
	
end

# ╔═╡ b5e8f30e-1573-4ec0-a20a-8673ab582d46

function cavity_plot(i::Int64, ts, εt, (αpss, αmss), (αp_list2, αm_list2); xlims=[-1.5, 1.5], ylims=[-1.5,1.5])

	colors=colors1
	
	αps1 = αpss[1:i]
	αms1 = αmss[1:i]
	αps2 = αp_list2[1:i]
	αms2 = αm_list2[1:i]

	pα = plot(xlims=xlims, ylims=ylims, xlabel=string("Re ", L"\alpha_\pm"), ylabel=string("Im ", L"\alpha_\pm"), legend=:outerbottomright, title="coherent state wavepacket centers", titlefontsize=12)
	
	plot!(real.(αps1), imag.(αps1), color=colors[2], linestyle=:dash, label=:none)
	plot!([real(last(αps1))], [imag(last(αps1))], color=colors[2], marker=:cross, label=string(L"\alpha_+", " (s.s.)"))
	plot!(real.(αms1), imag.(αms1), color=colors[4], linestyle=:dash, label=:none)
	plot!([real(last(αms1))], [imag(last(αms1))], color=colors[4], marker=:cross, label=string(L"\alpha_-", " (s.s.)"))

	plot!(real.(αps2), imag.(αps2), color=colors[2], label=:none)
	plot!([real(last(αps2))], [imag(last(αps2))], color=colors[2], marker="o", label=L"\alpha_+")
	plot!(real.(αms2), imag.(αms2), color=colors[4], label=:none)
	plot!([real(last(αms2))], [imag(last(αms2))], color=colors[4], marker="o", label=string(L"\alpha_-"))

	plot!(ts, εt, color=colors[5], label=:none, xlabel = "t (μs)", ylabel = "ε (MHz)", legend=:none, inset = (1, bbox(0.05, 0.05, 0.4, 0.25, :top, :right)), subplot=2, title="cavity drive", titlefontsize=10, axisfontsize=10)
	plot!([ts[i]], [εt[i]], color=colors[5], marker="o", label=:none, subplot=2)
	
end

# ╔═╡ 45ad1a01-dc16-464c-a44d-2259542d4004
function qubit_plot((t, exps, rs); title="", legendpos=:bottomleft, rlims=:default, rlabels=["r1", "r2"])

	record = (rs != [])

	basis = qbasis
	colors = colorsq1
	labels = qlabels

	p = 0.5 .* (1 .+ exps[1].^2 .+ exps[2].^2 .+ exps[3].^2)

	pl = plot(size=(600,300), legendfontsize=10, titlefontsize=12, legend=:outerright, ylabel="bloch coordinates", xlabel = record ? "" : "t (μs)", linewidth=1.5, title=title)

	for l in 1:length(basis)
		label = labels[l]
		color = colors[l]
		exp = exps[l]
		plot!(t, exps[l], color=color, label=label, legend=legendpos, ylims=[-1,1], linewidth=1.5)
	end

	plot!(t, p, color=colors[4], label=L"Tr(\rho^2)")

	if !record
		return pl
	else
		l = @layout [blochs{0.6h}; record{0.4h}]
		pr = plot(xlabel="t (μs)", ylabel="record", label=:none, legend=legendpos, title="", linewidth=0.8, ylims=rlims)
		for (i, r) in enumerate(rs)
			plot!(t, r, color=mixed_colors[i], label=rlabels[i])
		end
		return plot(pl, pr, layout = l, link=:y)
	end
end

# ╔═╡ 9b7b91e0-921d-4178-b7d4-57a5d0659db5
function qubit_plot(sol1::Solution, sol2::Solution; record=false, title="", color1=colorsq1, color2=colorsq2, l1="", l2="", rlims=:default)

	basis = qbasis

	t1 = sol1.t
	exps1 = map(op -> expectations(sol1, op), basis)
	r1 = record ? sol1.r[1] : []
	
	t2 = sol2.t
	exps2 = map(op -> expectations(sol2, op), basis)
	r2 = record ? sol2.r[1] : []
	
	return qubit_plot((t1, exps1, r1), (t2, exps2, r2); title=title, color1=color1, color2=color2, l1=l1, l2=l2, rlims=rlims)
	
end

# ╔═╡ 748a5310-3c5e-44fc-b3f9-0d994fea53eb
function qubit_plot((t1, exps1, r1), (t2, exps2, r2); title="", color1=colorsq1, color2=colorsq2, l1="", l2="", rlims=:default)

	basis = qbasis
	labels = qlabels

	record = ((r1 != []) || (r2 != []))

	p1 = 0.5 .* (1 .+ exps1[1].^2 .+ exps1[2].^2 .+ exps1[3].^2)
	p2 = 0.5 .* (1 .+ exps2[1].^2 .+ exps2[2].^2 .+ exps2[3].^2)

	pl = plot(size=(600,300), legendfontsize=10, titlefontsize=12, legend=:outerright, ylabel="bloch coordinates", xlabel = record ? "" : "t (μs)", linewidth=1.2, title=title)

	for l in 1:length(basis)
		plot!(t1, exps1[l], color=color1[l], label=labels[l], legend=:outerright, ylims=[-1,1], linewidth=1.2)
	end
	plot!(t1, p1, color=color1[4], label=L"Tr(\rho^2)", linewidth=1.2)
	
	for l in 1:length(basis)
		plot!(t2, exps2[l], color=color2[l], label=:none, legend=:outerright, ylims=[-1,1], linestyle=:dash, linewidth=2)
	end
	plot!(t2, p2, color=color2[4], label=:none, linestyle=:dash, title=string(l1, " ---, ", l2, " - - -"), linewidth=2)


	if !record
		return pl
	else
		l = @layout [blochs{0.6h}; record{0.4h}]
		pr = plot()
		if r1 != []	
			plot!(t1, r1, color=mixed_colors[1], xlabel="t (μs)", ylabel="record", label=:none, legend=:outerright, title="", ylims=rlims)
		end
		if r2 != []
			plot!(t2, r2, color=mixed_colors[2], xlabel="t (μs)", ylabel="record", label=:none, legend=:outerright, title="", ylims=rlims)
		end
		return plot(pl, pr, layout = l, link=:both)
	end
end

# ╔═╡ 6089fbf5-9764-4b20-bd5d-4c134482da64
qubit_plot(solss, record=true, legendpos=:bottomright, rlabels=["I", "Q"])

# ╔═╡ 30f58c5b-6869-4949-ae7c-01ca263b6ee3
md"""
## Functions
"""

# ╔═╡ d7ff20d4-674d-4248-a493-f63bff94aefe
function envelope(t, funcs...)
	for (f, T) in funcs
		if t >= T[1] && t < T[2]
			return f(t)
		end
	end
	
	return 0
	
end

# ╔═╡ 18c5a83d-cf4f-4303-842d-c800193e52b1
let
	fε = 0.02
	ωε = 2π * fε # MHz
	
	global t0 = 0 # μs
	t1 = (1/fε)/2 # μs
	global t3 = 100 # μs
	t2 = t3 - (1/fε)/2 # μs
	
	ε1(t) = (ΩR0/2) * (1 - cos(ωε * t))
	ε2(t) = ΩR0
	ε3(t) = (ΩR0/2) * (1 - cos(ωε * (t - t2) + π))
	
	ε(t) = envelope(t, (ε1, (t0, t1)), (ε2, (t1, t2)), (ε3, (t2, t3)))
	
	global ts = range(t0, t3, step=dt)
	global εt = ε.(ts)
	
	md" 🌀 drive definition"
end

# ╔═╡ c728af0d-72e3-4596-9f27-01f82fe86437
begin
	αpss(ε, Δ, χ, κ) = (2ε/κ) * -im / (1 + im * (2(Δ + χ)/κ))
	αmss(ε, Δ, χ, κ) = (2ε/κ) * -im / (1 + im * (2(Δ - χ)/κ))
	
	αp_list = αpss.(εt, Δ1, χ, κ)
	αm_list = αmss.(εt, Δ1, χ, κ)
	md" 🌀 function definitions"
end

# ╔═╡ 8e19c7cf-bdf8-494b-9e2e-eeff6c655f09
# RUNGE KUTTA INTEGRATION OF αp, αm

function cavity_solution((t0, tf), ε::Function, Δ, χ, κ; dt=1e-3)

	ts = range(t0, tf, step=dt)

	α̇p(t, αp, ε) = -im * (Δ + χ) * αp - im * ε(t) - (κ/2) * αp
	α̇m(t, αm, ε) = -im * (Δ - χ) * αm - im * ε(t) - (κ/2) * αm

	
	αp = αpss(ε(t0), Δ, χ, κ)
	αm = αmss(ε(t0), Δ, χ, κ)

	αplist = [αp]
	αmlist = [αm]	

	for i in 2:length(ts)
		
		t = ts[i]

		# update αp
		K1p = dt * α̇p(t, αp, ε)
		K2p = dt * α̇p(t + dt/2, αp + K1p/2, ε)
		K3p = dt * α̇p(t + dt/2, αp + K2p/2, ε)
		K4p = dt * α̇p(t + dt, αp + K3p, ε)

		αp = αp + (K1p + 2*K2p + 2*K3p + K4p)/6
		push!(αplist, αp)

		# update αm
		K1m = dt * α̇m(t, αm, ε)
		K2m = dt * α̇m(t + dt/2, αm + K1m/2, ε)
		K3m = dt * α̇m(t + dt/2, αm + K2m/2, ε)
		K4m = dt * α̇m(t+ dt, αm + K3m, ε)

		αm = αm + (K1m + 2*K2m + 2*K3m + K4m)/6
		push!(αmlist, αm)

	end

	return (αplist, αmlist)

md" 🔻 solving Runge-Kutta equations"
end

# ╔═╡ 3f41b58d-8ce0-4140-8061-a5413771e73e
md"""
Note that here, we are working in the bad cavity limit since 

 $\kappa$ = $(κ) MHz, 

 $\chi$ = $χ MHz, and

 $\Delta$ =  0 MHz. 
 
 In this regime, the cavity photon number approximation $\overline n \approx |2\epsilon / \kappa|^2$ = $(round(abs(2ΩR0 / κ)^2, digits=3)) is approximately equal to the exact 

 $\overline n = \text{Re} (\alpha_+^* \alpha_-)$ = $(round(real(αpss(ΩR0, 0, χ, κ)' * αmss(ΩR0, 0, χ, κ)), digits=3)). As a result, the measurement rates are approximately equal, and the approximation is roughly valid.
"""

# ╔═╡ abd37631-e506-49e2-ab38-5c06d380ddcf
if anim_plot1
	animate_plot(ts, cavity_plot, ts, εt, (αp_list, αm_list); step=150, fps=10)
else
	cavity_plot(ts, εt, (αp_list, αm_list); xlims=[-0.1,0.1], ylims=[-0.7,0.7])
end

# ╔═╡ 2b5dc7be-d3c8-40b8-b13a-d36b51f114d8
md"""

t : $(first(ts)) μs
$(@bind j html"<input type=range min=1 max=10001 step=10 value=1>")
$(last(ts)) μs

"""



# ╔═╡ e0db1833-9b68-4254-8556-c76d0162a3d1
md"""
t = $(ts[j]) μs
"""

# ╔═╡ 0d85d373-e094-4864-8a1a-83936c3362d5
let
	ωε = 2π * fε # MHz
	
	t0 = 0 # μs
	t1 = (1/fε)/2 # μs
	t3 = 100 # μs
	t2 = t3 - (1/fε)/2 # μs
	
	ε1(t) = (ΩR0/2) * (1 - cos(ωε * t))
	ε2(t) = ΩR0
	ε3(t) = (ΩR0/2) * (1 - cos(ωε * (t - t2) + π))

	global ε(t) = envelope(t, (ε1, (t0, t1)), (ε2, (t1, t2)), (ε3, (t2, t3)))
	
	md" 🌀 drive and function definitions"
end

# ╔═╡ f6456101-0876-4445-b0af-a00370077460
(αplist, αmlist) = cavity_solution((t0, t3), ε, Δ1, χ, κ; dt=dt)

# ╔═╡ 81a7ea85-d3aa-4c67-8ec3-26644ae04d5b
let
	εt = ε.(ts)
	αp_ss = αpss.(εt, Δ1, χ, κ)
	αm_ss = αmss.(εt, Δ1, χ, κ)

	if anim_plot2
		animate_plot(ts, cavity_plot, ts, εt, (αp_ss, αm_ss), (αplist, αmlist); step=150, fps=10)
	else
		cavity_plot(ts, εt, (αp_ss, αm_ss), (αplist, αmlist); xlims=[-0.1,0.1], ylims=[-0.7,0.7])
	end

end

# ╔═╡ fe014ea3-d8d5-4437-81b6-2a688ac9a41b
getclosest(array, val) = argmin(abs.(val .- array))

# ╔═╡ 44e14f46-1940-4e84-b6ad-c42c398e7da9
function attime(array, time, ts)
	i = getclosest(ts, time)
	return array[i]
end

# ╔═╡ fd18414e-6737-4f08-8104-2f21a3c1eebf
let
	ψ0 = normalize(g + e)
	dt = 1e-3  # integration time-step
	tf = 10.0
	ts = range(0.0, tf, step=dt)

	# ε0 = 1 # MHz
	κ = 7.5 # MHz
	χ = 0.25 # MHz
	Δ = 0 # MHz

	ε = t -> ΩR0
	
	(αplist, αmlist) = cavity_solution((0.0, tf), ε, Δ, χ, κ; dt=dt)
	
	αp(t) = attime(αplist, t, ts)
	αm(t) = attime(αmlist, t, ts)
	nw(t) = αp(t)' * αm(t)
	Δα(t) = αp(t) - αm(t)

	

	Γ̃(t) = 2χ * imag(nw(t)) - (κ * η / 2) * abs(Δα(t))^2
	Δω̃q(t) = 2χ * real(nw(t)) - (κ * η) * imag(nw(t)) - η * real(ε(t)' * Δα(t))
	Γm(t) = (κ / 2) * abs(Δα(t))^2


	# evolution operators -------------------
	J = [(σz, Γ̃)]
	C = [(exp(im * ϕ) * (σz + ωs * Iq), Γm, η/2), (exp(im * (ϕ + π/2)) * (σz + ωs * Iq), Γm, η/2)]
	# C = [(exp(im * ϕ) * σz, Γ, η/2), (exp(im * (ϕ + π/2)) * σz, Γ, η2/)]
	H(t) = (ΩR / 2) * σy + (Δω̃q(t) / 2) * σz

	global solns = bayesian((0.0, tf), ψ0, H, J, C; dt=dt, records=solss.r)
	
	
	
end

# ╔═╡ c2f2d30f-1290-4ce0-8b90-a12648bfc76d
qubit_plot(solns, record=true, legendpos=:bottomright, rlabels=["I", "Q"])

# ╔═╡ 699bd693-cc49-4ab2-a07f-7553acd613e7
sine_envelope(t, rt, tf, εmax) = let ω = π/rt
	if t > 0.0 && t < rt
		(εmax/2) * (1 - cos(ω * t))	
	elseif t >= rt && t < (tf - rt)
		εmax
	elseif t >= (tf - rt) && t < tf
		(εmax/2) * (1 - cos(ω * (t - (tf - rt)) + π))
	else
		0.0
	end
end

# ╔═╡ 9af66e2b-012f-44ec-8b8f-aedb2a773c7c
let
	ψ0 = normalize(g + e)
	dt = 1e-3  # integration time-step
	tf = 20.0
	ts = range(0.0, tf, step=dt)

	κ = 7.5 # MHz
	χ = 0.25 # MHz
	Δ = 0 # MHz

	# cavity solutions ------------------------------------------------------------
	ε(t) = sine_envelope(t, 4.0, 20.0, ΩR0)		
	
	(αplist, αmlist) = cavity_solution((0.0, tf), ε, Δ, χ, κ; dt=dt)
	αpsslist = map(εt -> αpss(εt, Δ, χ, κ), ε.(ts))
	αmsslist = map(εt -> αmss(εt, Δ, χ, κ), ε.(ts))
	nsslist = map(zip(αpsslist, αmsslist)) do (αp, αm)
				real(αp' * αm)
	end
	
	αp(t) = attime(αplist, t, ts)
	αm(t) = attime(αmlist, t, ts)
	nw(t) = αp(t)' * αm(t)
	Δα(t) = αp(t) - αm(t)

	Γ̃(t) = 2χ * imag(nw(t)) - (κ * η / 2) * abs(Δα(t))^2
	Δω̃q(t) = 2χ * real(nw(t)) - (κ * η) * imag(nw(t)) - η * real(ε(t)' * Δα(t))
	Γm(t) = (κ / 2) * abs(Δα(t))^2

	# cavity plots ------------------------------------------------------------
	p0 = plot(ts, Γ̃.(ts), color=colors1[1], ylabel=L"\tilde{\Gamma}", legend=:none)
	p1 = plot(ts, Γm.(ts), color=colors1[2], ylabel=L"\Gamma_m", legend=:none)
	p2 = plot(ts, real.(nw.(ts)), color=colors1[3], ylabel=string("Re", L"[n_w(t)]"), label="actual")
	plot!(ts, nsslist, color=colors1[3], linestyle=:dash, ylabel=string("Re ", L"[n_w(t)]"), label="steady state", legend=:bottomright, xlabel="t (μs)")
	p3 = cavity_plot(ts, ε.(ts), (αpsslist, αmsslist), (αplist, αmlist); xlims=[-0.1,0.1], ylims=[-0.7,0.7])

	global cavity_plots = (p0, p1, p2, p3)


	# evolution operators ---------------------------------------------------------
	J = [(σz, Γ̃)]
	C = [(exp(im * ϕ) * (σz + ωs * Iq), Γm, η/2), (exp(im * (ϕ + π/2)) * (σz + ωs * Iq), Γm, η/2)]
	# C = [(exp(im * ϕ) * σz, Γ, η/2), (exp(im * (ϕ + π/2)) * σz, Γ, η2/)]
	H(t) = (ΩR / 2) * σy + (Δω̃q(t) / 2) * σz

	# solution --------------------------------------------------------------------
	global sol3 = bayesian((0.0, tf), ψ0, H, J, C; dt=dt)
	
	
	
end

# ╔═╡ 81ebe638-7c9f-493f-8fa2-953cbec48551
plot(cavity_plots[1], cavity_plots[2], cavity_plots[3], layout=grid(3,1))

# ╔═╡ 418628f1-27df-4f1c-858e-b445ea865cd6
qubit_plot(sol3, record=true, legendpos=:bottomright, rlims=[-1000,1000], rlabels=["I","Q"])

# ╔═╡ df300433-9581-4b88-83e2-de1742b497fd
md"""
## Other
"""

# ╔═╡ d304e60c-f47f-4b75-b34d-4ef5dc7b5750
green(text; title="Note") = Markdown.MD(Markdown.Admonition("correct", title, [text]))

# ╔═╡ a7332fed-6089-4ffd-a93e-8893024d87e8
red(text; title="Note") = Markdown.MD(Markdown.Admonition("danger", title, [text]))

# ╔═╡ dfb8cf89-86f5-4ab9-aac3-5867130fd6bd
tan(text; title="Note") = Markdown.MD(Markdown.Admonition("warning", title, [text]))

# ╔═╡ b01d4d05-7b32-4bd4-b554-e85abca61595
blue(text; title="Note") = Markdown.MD(Markdown.Admonition("note", title, [text]))

# ╔═╡ 91e5577b-fe11-4d8e-b37e-813a81621385
if Δ == 0.
	green(md"""When the cavity is driven on its bare resonance, coherent states are maximally distinguishable in phase.""",title="Driving bare resonance")
	
elseif Δ == χ
	blue(md"""A detuning of Δ = 1 MHz = χ corresponds to driving on resonance with the α+ coherent state.""",title="Driving α+ state.")
	
elseif Δ == -χ
	blue(md"""A detuning of Δ = 1 MHz = -χ  corresponds to driving on resonance with the α- coherent state.""", title="Driving α- state.")
	
elseif abs(Δ) > 3χ 
	red(md"""When the detuning is large, the coherent states are not easily distinguished in either amplitude or phase.""", title="Driving off resonance.")

else
	tan(md"""Change the value of Δ (resonator-drive detuning) using the slider to see how it affects the conditional coherent states.""")

		
end

# ╔═╡ 7a17cf5d-db67-499c-886f-2c56683719f9
hint(text; title="Hint") = Markdown.MD(Markdown.Admonition("hint", title, [text]))

# ╔═╡ Cell order:
# ╟─377a3336-20bd-4baa-a033-af8bbc8668a8
# ╟─13f42433-1bc7-4f5e-b97c-a4a2664e1a71
# ╟─3dfd0a4f-ed67-4c19-9fdb-98b3b96c62d1
# ╟─1c6ac0e4-690d-4eff-9a88-68883abc5f67
# ╟─e4ffe1eb-c206-4ea9-9536-8942d2ae7e36
# ╟─db847000-8785-4bcb-9f93-db82f6ead06d
# ╟─b7be6342-01f9-4afd-9250-27f0479ed5c5
# ╟─306f55d1-fd48-4830-ba0a-b54ca5dccf82
# ╟─56f69b9a-af7e-4e50-bf6a-0d4e436a7f95
# ╟─c681b1d2-a32b-4c10-bc80-6447368f13eb
# ╟─d5b3318e-a1a8-48c2-baf7-92b398902aff
# ╟─ce3f9a90-70ed-4d73-b7ac-e61791eb1aa9
# ╟─543b95ab-589b-4523-9bf5-b955d249699e
# ╟─9e816960-3731-4b98-96c8-9aa689bdf1ac
# ╟─83f4d601-33f2-4c3d-92c4-37231d55cae3
# ╟─54e9278d-7757-4af9-968b-c78dfbe04758
# ╟─91e5577b-fe11-4d8e-b37e-813a81621385
# ╟─9c860ffd-7a2d-45df-84a8-1ddcf0f46298
# ╟─897d6364-90ee-41cc-ab70-d5ef0c537bd2
# ╟─67d7fdb7-145f-42f3-8dec-cb6fc8bb3018
# ╟─e73546a2-2b9e-4062-89f2-5a0913209988
# ╟─05def7cc-6c0e-496e-9429-1fd457575b31
# ╟─38702e03-a6ff-4c7a-bbb6-4dd0bcf89ba8
# ╟─18c5a83d-cf4f-4303-842d-c800193e52b1
# ╟─c728af0d-72e3-4596-9f27-01f82fe86437
# ╟─00257bca-f1bd-4c4b-af0b-db0bafd5e406
# ╠═abd37631-e506-49e2-ab38-5c06d380ddcf
# ╟─ca2e2708-789b-4a95-ba3f-0402f93cb84e
# ╠═f5719cfb-1acf-4f87-b12e-a0ed65f635d5
# ╟─cf3fd908-d968-4581-a8b4-5aac46c4afb0
# ╟─c5141459-60d0-47ed-90f3-ac4e741e5707
# ╟─8e19c7cf-bdf8-494b-9e2e-eeff6c655f09
# ╟─0d85d373-e094-4864-8a1a-83936c3362d5
# ╟─71ea6e65-46b3-4c6c-a919-87dd31a6674f
# ╠═f6456101-0876-4445-b0af-a00370077460
# ╟─2b5dc7be-d3c8-40b8-b13a-d36b51f114d8
# ╟─e0db1833-9b68-4254-8556-c76d0162a3d1
# ╟─83e516e3-20dd-4cbc-9f6a-84e1662d0dcb
# ╠═81a7ea85-d3aa-4c67-8ec3-26644ae04d5b
# ╟─3ca15002-9b64-416c-94e3-a7904766a87b
# ╟─a9dc195d-fe8f-4c6b-8078-597585b2ee3f
# ╟─14cccb18-0187-4c91-ba53-b920ce0c3704
# ╟─e8299f27-e74b-47d1-8349-7e28c3ea4495
# ╟─adb901f4-b296-4549-a92f-b075b40e8f8f
# ╟─69055a05-e2f1-4c8a-ad2b-1ca8ac0dadc8
# ╟─e4269abf-9b63-4fbe-8001-397a449dd0e6
# ╟─58e7e646-eaf8-4ba0-8204-9819a054fe65
# ╟─02d2a1f6-12af-4345-b42a-db2b11e2a883
# ╟─498bfa13-4278-4808-b1d6-279eabcd9ed8
# ╟─9423dbb2-c325-415a-aeca-31cc8b205173
# ╟─50ac7363-0591-49dc-9bdf-dc37dcf43f60
# ╟─c7eaf3d6-73ac-4291-ab9c-5857cb77cdaf
# ╟─d0528468-7e59-41a1-a81a-5d8933864f49
# ╟─f568802e-e754-4bad-a221-5d6102e2068f
# ╟─fca59554-0975-415a-8143-56b86d3de5bc
# ╟─0c326351-6d79-4e16-9bf0-8170fd92a671
# ╟─4e368417-1107-45af-9af8-10fc16eaabe5
# ╟─ebfe6948-64be-4190-8470-d5e196ee43ec
# ╟─23b8def8-5019-41ab-8e55-b70ee857b7a6
# ╟─67c05b7c-4376-495a-8a21-b5d4f749c0d6
# ╟─3c50b081-f1d0-40a4-994d-d1d4813b98f4
# ╟─fc628bbf-47cf-4603-bba6-52bcb4eed93f
# ╟─3f41b58d-8ce0-4140-8061-a5413771e73e
# ╟─a0b3eef2-9d43-441a-b541-e7ef430d86be
# ╟─6089fbf5-9764-4b20-bd5d-4c134482da64
# ╟─6af54bde-991c-4810-b5a4-7be74f743ee3
# ╠═c76e3314-e3a4-44d3-844e-cf141b1f550a
# ╠═ab601ea6-69cb-4079-9cb7-a57ecfd15b12
# ╠═c9efabea-2802-4bcd-adbf-3608155fba58
# ╟─2b8f2e1c-0f25-4be0-ab40-2bb0d96ba5e8
# ╟─6cb5d42f-b074-45c4-8f19-a70f04f4fad2
# ╟─44b8a3f8-b220-4567-a486-576908bd4255
# ╟─fd18414e-6737-4f08-8104-2f21a3c1eebf
# ╟─c2f2d30f-1290-4ce0-8b90-a12648bfc76d
# ╟─e13483e5-41b3-438e-9504-2b836d544c00
# ╟─9af66e2b-012f-44ec-8b8f-aedb2a773c7c
# ╠═81ebe638-7c9f-493f-8fa2-953cbec48551
# ╟─418628f1-27df-4f1c-858e-b445ea865cd6
# ╟─1fc48bb8-d254-470d-b3bc-76acc342c6b6
# ╟─7e62e6c6-aa3c-4350-901c-15d017b8db42
# ╟─04f6e40d-2dfa-46fd-8e43-58c3f647b5b5
# ╟─4029ecf7-7983-465f-bc8c-b8b6fb25e47e
# ╠═b5e8f30e-1573-4ec0-a20a-8673ab582d46
# ╟─62b5300d-32c5-477f-8150-89d90b5ae31c
# ╟─6038bb81-8e47-4be6-a5be-715101349309
# ╟─29b15810-41e0-4e2d-b676-7e7378cd261a
# ╟─45ad1a01-dc16-464c-a44d-2259542d4004
# ╟─9b7b91e0-921d-4178-b7d4-57a5d0659db5
# ╟─748a5310-3c5e-44fc-b3f9-0d994fea53eb
# ╟─7af4d4ff-30ad-4ec8-9aea-ca5c2b431f57
# ╟─d88511e7-af21-4b0b-a654-348060892498
# ╠═b07e6b7b-b96e-4d0d-ac42-e3b717d5163a
# ╟─30f58c5b-6869-4949-ae7c-01ca263b6ee3
# ╟─44e14f46-1940-4e84-b6ad-c42c398e7da9
# ╟─d7ff20d4-674d-4248-a493-f63bff94aefe
# ╟─fe014ea3-d8d5-4437-81b6-2a688ac9a41b
# ╟─699bd693-cc49-4ab2-a07f-7553acd613e7
# ╟─df300433-9581-4b88-83e2-de1742b497fd
# ╟─d304e60c-f47f-4b75-b34d-4ef5dc7b5750
# ╟─a7332fed-6089-4ffd-a93e-8893024d87e8
# ╟─dfb8cf89-86f5-4ab9-aac3-5867130fd6bd
# ╟─b01d4d05-7b32-4bd4-b554-e85abca61595
# ╟─7a17cf5d-db67-499c-886f-2c56683719f9
# ╠═1aa3dd28-34e8-4ebe-bc07-8ec41527bb1d
