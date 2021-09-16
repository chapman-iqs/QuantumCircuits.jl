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
	using Random
	using Statistics
	using Distributions
	using PyPlot
	using QuantumCircuits
end

# ╔═╡ 377a3336-20bd-4baa-a033-af8bbc8668a8
md"""
# Coherent states of the resonator

In this interactive notebook, we'll look at properties and dynamics of coherent states, which are key to understanding superconducting qubit readout. This notebook is based on unpublished notes [1].
"""

# ╔═╡ db847000-8785-4bcb-9f93-db82f6ead06d
md"""
## Problem setup

##### Qubit-resonator interaction

A qubit with bare frequency $\omega_q \sim 4-9$ GHz is dispersively coupled to a microwave resonator with bare frequency $\omega_r$ detuned by $\Delta \equiv \omega_q - \omega_r \sim 1$ GHz from the qubit. The resonator frequency shifts by $2\chi$ depending on the qubit state, with $\chi \sim 1$ MHz, i.e. the resonator has dressed frequency $(\omega_r)^{\pm} = \omega_r \pm \chi$, depending on whether the qubit is excited ($+$) or ground ($-$). The qubit is also coupled to a (time-dependent) Rabi drive with frequency $\Omega_R(t) \ll \omega_q$. The Hamiltonian describing this interaction is

```math
\hat H_{qr} / \hbar = \omega_q \hat \sigma_z/2 + (\omega_r + \chi \sigma_z) \hat a^\dagger \hat a+ \Omega_R(t) \hat \sigma_y /2.
```

##### Resonator-transmission line interaction

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

=-i(\omega_r + \chi \hat \sigma_z) \hat a + (\sqrt{\kappa}/2)(\hat c_\text{in} - \hat  c_\text{in}^\dagger - \hat c_\text{out} + \hat c_\text{out}^\dagger),
```

which with the boundary condition becomes

```math
\dot{\hat{a}} = -i(\omega_r + \chi \hat \sigma_z) \hat a + \sqrt{\kappa} (\hat c_\text{in} - \hat  c_\text{in}^\dagger) - \frac{\kappa}2 (\hat a - \hat a^\dagger)
```
"""

# ╔═╡ 306f55d1-fd48-4830-ba0a-b54ca5dccf82
md"""
##### Rotating-wave and coherent-field approximations
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
## Characterizing the resonator
"""

# ╔═╡ d5b3318e-a1a8-48c2-baf7-92b398902aff
md"""
#### Steady-state characteristics
"""

# ╔═╡ ce3f9a90-70ed-4d73-b7ac-e61791eb1aa9
md"""
In the case where $\epsilon(t) = \epsilon$ is constant, setting the differential equation for $\dot{\alpha}_\pm(t) = 0$ immediately yields steady state amplitudes of the resonator:

$$\alpha_\pm^{(s.s.)} \equiv \frac{2 \epsilon}{\kappa} \frac{-i}{1 + i[2(\Delta \pm \chi)/\kappa]},$$

where $\Delta \equiv \omega_r - \omega_d$ is the drive detuning from the bare resonator frequency. The steady-state amplitudes allow us to characterize the cavity (and the measurement) in several important ways.

"""

# ╔═╡ 543b95ab-589b-4523-9bf5-b955d249699e
begin
	ε0 = 1 # MHz
	κ = 2 # MHz
	χ = 1 # MHz
	Δs = range(-5, 5, step=0.1) # MHz
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
##### State distinguishability
"""

# ╔═╡ 54e9278d-7757-4af9-968b-c78dfbe04758
md"""
Superconducting qubit readout depends on the distinguishability of cavity coherent states conditioned on the state of the qubit. Since resonator amplitudes are complex numbers, 

$\alpha_\pm = |\alpha_\pm| e^{i \phi_\pm},$


they can be distinguished in amplitude or in phase. Alternatively, they can be plotted in the complex plane via their real and imaginary parts (called "quadratures"). 
"""

# ╔═╡ 9c860ffd-7a2d-45df-84a8-1ddcf0f46298
md"""

Δ : -4
$(@bind Δ html"<input type=range min=-4 max=4 step=0.1 value=0.1>")
4

"""



# ╔═╡ 897d6364-90ee-41cc-ab70-d5ef0c537bd2
md"""
Δ = $Δ
"""

# ╔═╡ 1c0539b3-e6ca-4869-aa65-e268c6842d22
md"""
##### Cavity photon number and measurement rate
"""

# ╔═╡ e73546a2-2b9e-4062-89f2-5a0913209988
md"""
##### Adiabatic evolution of coherent states with time-dependent drive
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

# ╔═╡ 11838257-1e11-484d-831e-10eeb28abb6d
md" Move the slider to evolve the system in time."

# ╔═╡ ca2e2708-789b-4a95-ba3f-0402f93cb84e
md"""
The code above assumes that the cavity is driven on resonance. Adjust the detuning below to see what changes.
"""

# ╔═╡ f5719cfb-1acf-4f87-b12e-a0ed65f635d5
md"""

Δ : -4 MHz
$(@bind Δ1 html"<input type=range min=-4 max=4 step=0.2 value=0>")
4 MHz

"""



# ╔═╡ e16b9185-418e-4846-96dc-c076ad781e70
md"""
Δ = $Δ1 MHz
"""

# ╔═╡ cf3fd908-d968-4581-a8b4-5aac46c4afb0
md"""
#### Non-adiabatic evolution of coherent states with time-dependent drive
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

# ╔═╡ 3ca15002-9b64-416c-94e3-a7904766a87b
md"""
Recall our second assumption, namely that the drive envelope $\epsilon(t)/\sqrt{\kappa}$ varies slowly compared to the drive frequency $\omega_d$. We can see already that there is a visible discrepancy between the steady state values $\alpha_+^{s.s.}, \alpha_-^{s.s.}$ and the time-dependent $\alpha_+, \alpha_-$. We'll can see an even greater discrepancy if we ramp-up the drive faster. Try changing fε:
"""

# ╔═╡ a9dc195d-fe8f-4c6b-8078-597585b2ee3f
md"""

fε : 0.01 MHz
$(@bind fε html"<input type=range min=0.01 max=0.12 step=0.01 value=0.01>")
0.12 MHz

"""

# ╔═╡ 17d6c0a1-9191-4da8-a298-e63fd4008453
md"""
fε = $fε MHz
"""

# ╔═╡ 14cccb18-0187-4c91-ba53-b920ce0c3704
md"""
In real qubit readout, the adiabatic approximation is clearly inadequate since pulses can last just 20 ns, with $\epsilon(t)$ at much larger values than we have considered.
"""

# ╔═╡ e8299f27-e74b-47d1-8349-7e28c3ea4495
md"""
## Reduced-qubit dynamics
"""

# ╔═╡ ec64b4d3-7631-4a1e-af3e-cf895f2fd922
md"""
##### (Cavity) steady-state dynamics
"""

# ╔═╡ 4e368417-1107-45af-9af8-10fc16eaabe5
md"""
In the following discussion, we assume the cavity is driven on resonance ($\Delta = 0$) and look at steady state values of qubit evolution parameters.
"""

# ╔═╡ 61e73b8f-1ee6-4a5c-aa11-beb85bdd80d1


# ╔═╡ 5a86ea71-0684-476f-8d5d-ee065b43c92f
let
	ε0 = 1 # MHz
	κ = 5 # MHz
	χ = 0.2 # MHz
	Δ = 0
	global η = 1
	global ΩR = 0
	
	global ωs = 2χ * abs(2ε0 / κ)^2
	global Γ = 8 * (χ^2) * abs(2ε0 / κ)^2/ κ
	
	global ω̃q = (1 - η) * ωs
	global Γ̃ = (1 - η) * Γ
	global τ = 1/(2Γ)
	
	global αpss = (2ε0/κ) * -im / (1 + im * (2(Δ + χ)/κ))
	global αmss = (2ε0/κ) * -im / (1 + im * (2(Δ - χ)/κ))

	md" 🌀 qubit / cavity parameters"
end

# ╔═╡ a3ec8af8-9f28-468a-8d9f-286c66e69164
let
	ε0 = 1 # MHz
	κ = 5 # MHz
	χ = 0.2 # MHz
	Δ = 0
	
	κ * imag(αpss' * αmss) + real(ε0 * (αpss - αmss))
end

# ╔═╡ 0ec0004b-92dc-451b-bc74-54a8d44d85a5
ωs

# ╔═╡ ebfe6948-64be-4190-8470-d5e196ee43ec
md"""
The time-dependent cavity states has several effects on the qubit evolution. The cavity photon number weak value,

$n_w^{s.s.} = (\alpha_+^{s.s.})^* \alpha_-^{s.s.},$

produces a weak dispersive shift 

$2\chi n_w = \omega_S + i \Gamma.$

The real part of the weak photon number is itself the effective photon number in the cavity:

$\bar n = \text{Re} (n_w).$ 

This determines the ensemble-averaged AC-Stark shift

$\omega_S = 2\chi \bar n \approx 2 \chi \Bigg| \frac{2 \epsilon}{\kappa} \Bigg|^2$ 

and the ensemble-averaged dephasing rate 

$\Gamma = 2 \chi \text{Im} (n_w) \approx \frac{8 \chi^2 \bar n}{\kappa}.$

The approximations are valid when $|\kappa| \gg |\Delta|, |\chi|$ (bad-cavity regime) and when the resonator is at steady-state, and lead to an effective qubit frequency

$\tilde{\omega}_q \approx \omega_q + (1 - \eta) 2 \chi \bar{n}$

and effective dephasing rate

$\tilde \Gamma \approx (1 - \eta)\Gamma$


"""

# ╔═╡ 67c05b7c-4376-495a-8a21-b5d4f749c0d6
md"""
$$\alpha_\pm^{(s.s.)} \equiv \frac{2 \epsilon}{\kappa} \frac{-i}{1 + i[2(\Delta \pm \chi)/\kappa]},$$


"""

# ╔═╡ 2b8f2e1c-0f25-4be0-ab40-2bb0d96ba5e8
md"""
##### Non-steady state dynamics
"""

# ╔═╡ e4269abf-9b63-4fbe-8001-397a449dd0e6
md"""
The time-dependency of cavity states has several effects on the qubit evolution. The time-evolving cavity photon number weak value and state separation

$n_w \equiv \alpha_+^* \alpha_-$

$\Delta \alpha \equiv \alpha_+ - \alpha_-$

lead to a time-dependent AC-Stark shift and time-dependent dephasing rate

$\Delta \tilde \omega_q \equiv 2 \chi \text{Re}(n_w) - \kappa \eta \text{Im}(n_w) - \eta \text{Re}(\epsilon^*\Delta \alpha)$

$\tilde{\Gamma} \equiv 2 \chi \text{Im} (n_w) - \frac{\kappa \eta}2 |\Delta \alpha|^2.$

(Here, $\alpha_+$, $\alpha_-$, $n_w$, $\Delta \alpha$,  and $\epsilon$ are all implicitly time-dependent.)

"""

# ╔═╡ e1345e44-5104-42bb-81a6-1e85db1ccf91
begin
	nwss = αpss' * αmss
	nss = real(nwss)
	Δαss = αpss - αmss
	
	md" 🌀 cavity characteristics"
end

# ╔═╡ 7780e6ad-5241-4bde-bf43-a7410d4c1650
κ

# ╔═╡ 202fde97-a8ed-44d9-a02e-6a7f564b348d
Γ

# ╔═╡ 07629b53-a846-4341-9333-0c442ad2fc96
(5 / 2) * abs(Δαss)^2

# ╔═╡ 8d353b07-2dc9-4611-bb94-fd413dabf120
2 * 0.2 * imag(nwss)

# ╔═╡ 64da2aa4-96ac-49b8-8028-4dad4698353e
ωs

# ╔═╡ b464fa05-4b99-49a7-838b-93826da4f0e6
2 * 0.2 * real(nwss)

# ╔═╡ dbaaafbb-18b0-43ff-9b20-8d337ca745a9
ω̃q

# ╔═╡ c2085fde-2ec3-4d27-8aca-9b6ecc8218b1
2 * 0.2 * real(nwss) - 5 * η * imag(nwss) - η * real(Δαss)

# ╔═╡ c0e189dc-26cc-4bbe-9b72-aa681dd57e46
η

# ╔═╡ dc9478a2-0532-4b18-bd6e-2a0942eef585
Γ̃

# ╔═╡ ec327412-51bd-47f8-b40e-b89b4679dc90
let
	ε = 1
	κ = 5
	χ = 0.2
	
	2χ * (2ε/κ)^2 * (4χ/κ) / ((1 - 4 * (χ/κ)^2 )^2 + 16 * (χ/κ)^2)
	
	8 * χ^2 * real(nwss) / κ
end

# ╔═╡ 02d2a1f6-12af-4345-b42a-db2b11e2a883
md"""
The AC-Stark shift leads to a time-dependent Hamiltonian in the rotating frame

$\hat H_q = \frac{\Delta \tilde \omega_q}2 \hat \sigma_z + \frac{\Omega_R(t)}2 \hat \sigma_y,$

which gives the discrete time-evolution operator

$\hat U_{\Delta t} = \exp(-i \Delta t \hat H_q).$

"""

# ╔═╡ 498bfa13-4278-4808-b1d6-279eabcd9ed8
md"""
The time-dependent dephasing rate leads to coherence decay

$\rho_{01}^{(q)} \mapsto \rho_{01}^{(q)} e^{-i \Delta t \tilde \Gamma}.$
This can also be expressed using a Lindblad operator $J = \sqrt{\tilde \Gamma} \hat \sigma_z$.

Finally, the random complex signal $\tilde r = r_\varphi + i r_{\varphi + \pi/2}$ leads to measurement backaction

$\hat M_{\tilde{r}} \equiv \exp \Bigg[\frac{\tilde r \Delta t}{2 \tau_m} e^{i(\varphi + \text{angle}(\Delta \alpha))}\hat \sigma_z \Bigg] = \Bigg[\frac{\tilde r \Delta t}{2 \tau_m} \hat \sigma_z \Bigg]$
"""

# ╔═╡ c7eaf3d6-73ac-4291-ab9c-5857cb77cdaf
md"""
This is equivalent to measuring $\hat \sigma_z$ with measurement collapse time $\tau_m$ and efficiency $\eta$. The second equality follows from our choosing $\varphi = - \text{angle} (\Delta \alpha)$.

The measurement collapse time is directly related to the ensemble measurement-dephasing rate

$\Gamma_m = \frac{\kappa}2 |\Delta \alpha|^2, \hspace{5mm} \tau_m = \frac1{2\Gamma_m}.$
"""

# ╔═╡ 7f1176b8-f6d3-4fd1-a23e-31710dcfff10
md"""
#### Operators

We begin by defining a two-level Hilbert space for the system. `QuantumCircuits.jl` uses `QuantumOptics.jl` as its backend: `SpinBasis`, `sigmax`, `identityoperator` and so forth are `QuantumOptics.jl` functions.
"""

# ╔═╡ 982565bc-b137-4c10-8e13-5fe37be86823
md" ##### Qubit Hilbert space"

# ╔═╡ 4284173a-be05-4b58-a8d9-7189301344fd
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
end

# ╔═╡ 682b7883-fb96-4ebb-a196-5d85c5f91418
begin
	# η = 0 evolution
	J0ss = [(σz, Γ)]
	C0ss = []
	H0ss = (ΩR / 2) * σy + (ωs / 2) * σz
	
	# stochastic evolution
	Jss = [(σz, Γ̃)]
	Css = [(σz, τ, η), (im*σz, τ, η)];
	Hss = (ΩR / 2) * σy + (ω̃q / 2) * σz
	md" 🌀 evolution Kraus operators (`H`, `C`, `J`)"
end

# ╔═╡ cd1f1146-0b21-4b01-89de-e8bf67ed7867
let
	ε0 = 1 # MHz
	κ = 5 # MHz
	χ = 0.2 # MHz
	Δ = 0
	
	ωs = 2χ * real(nwss)
	Δω̃q = ωs - κ * η * imag(nwss) - η * real(ε0' * Δαss)
	
	Γm = (κ / 2) * abs(Δαss)^2
	τ = 1/(2Γm)
	
	Γ = 2χ * imag(nwss)
	Γ̃ = Γ - η * Γm
	
	# η = 0 evolution
	global J0 = [(σz, Γ)]
	global C0 = []
	global H0 = (ΩR / 2) * σy + (ωs / 2) * σz
	
	# stochastic evolution
	global J = [(σz, Γ̃)]
	global C = [(σz, τ, η)];
	global H = (ΩR / 2) * σy + (Δω̃q / 2) * σz
	
	md" 🌀 parameters and Kraus operators"

end

# ╔═╡ e51ca599-2452-4bfb-9f64-f27751d80aa0
md"""
### Citations

[1] Dressel, J. Phase-sensitive Qubit Monitoring with Resonator Transients. Unpublished notes, 2019.
"""

# ╔═╡ 7e62e6c6-aa3c-4350-901c-15d017b8db42
md" ## Utilities "

# ╔═╡ d7ff20d4-674d-4248-a493-f63bff94aefe
function envelope(t, funcs...)
	for (f, T) in funcs
		if t >= T[1] && t < T[2]
			return f(t)
		end
	end
	
	return 0
	
end

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

# ╔═╡ 18c5a83d-cf4f-4303-842d-c800193e52b1
let
	fε = 0.02
	ωε = 2π * fε # MHz
	
	t0 = 0 # μs
	t1 = (1/fε)/2 # μs
	t3 = 100 # μs
	t2 = t3 - (1/fε)/2 # μs
	
	ε1(t) = (ΩR0/2) * (1 - cos(ωε * t))
	ε2(t) = ΩR0
	ε3(t) = (ΩR0/2) * (1 - cos(ωε * (t - t2) + π))

	
	# ε(t) = 	if 		t <= t1 	ε1(t)
	# 		elseif 	t <= t2 	ε2(t)
	# 		elseif  t <= t3		ε3(t)
	# 		else  				0
	# end
	
	ε(t) = envelope(t, (ε1, (t0, t1)), (ε2, (t1, t2)), (ε3, (t2, t3)))
	
	global ts = range(t0, t3, step=dt)
	global εt = ε.(ts)
	
	md" 🌀 drive definition"
end

# ╔═╡ c728af0d-72e3-4596-9f27-01f82fe86437
begin
	αp(Δ1, ε) = (2ε/κ) * -im / (1 + im * (2(Δ1 + χ)/κ))
	αm(Δ1, ε) = (2ε/κ) * -im / (1 + im * (2(Δ1 - χ)/κ))
	
	αp_list = αp.(Δ1, εt)
	αm_list = αm.(Δ1, εt)
	md" 🌀 function definitions"
end

# ╔═╡ a86b6204-1a3a-45e8-96e6-21dd9369f059
md"""

t : $(first(ts)) μs
$(@bind i html"<input type=range min=1 max=10001 step=10 value=1>")
$(last(ts)) μs

"""


# ╔═╡ 0f1146c2-2675-492f-883d-d6a353aab727
md" t = $(ts[i]) μs "

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

# ╔═╡ 5fe7e35f-c2d1-425d-96b4-bc3e2c377f61
# RUNGE KUTTA INTEGRATION OF αp, αm

begin

	αplist = [αp(Δ1, ε(ts[1]))]
	αmlist = [αm(Δ1, ε(ts[2]))]

	for i in 2:length(ts)
		
		t = ts[i]

		# update αp
		K1p = dt * α̇p(t, αplist[i-1], ε)
		K2p = dt * α̇p(t + dt/2, αplist[i-1] + K1p/2, ε)
		K3p = dt * α̇p(t + dt/2, αplist[i-1] + K2p/2, ε)
		K4p = dt * α̇p(t + dt, αplist[i-1] + K3p, ε)
		push!(αplist, αplist[i-1] + (K1p + 2*K2p + 2*K3p + K4p)/6)

		# update αm
		K1m = dt * α̇m(t, αmlist[i-1], ε)
		K2m = dt * α̇m(t + dt/2, αmlist[i-1] + K1m/2, ε)
		K3m = dt * α̇m(t + dt/2, αmlist[i-1] + K2m/2, ε)
		K4m = dt * α̇m(t+ dt, αmlist[i-1] + K3m, ε)
		push!(αmlist, αmlist[i-1] + (K1m + 2*K2m + 2*K3m + K4m)/6)

	end

md" 🔻 solving Runge-Kutta equations"
end

# ╔═╡ f3c97871-bcf3-45f2-aae1-2c5d72461d4d
function rms(ser1::Array, ser2::Array)
	l = min(length(ser1), length(ser2))
	sqrt(sum((ser1[1:l] - ser2[1:l]).^2))/l
end

# ╔═╡ 3274ad19-1b27-4090-8aee-ef8c148b5e8e
expects = ρ -> collect(real(expect(ρ, s)) for s in [σx,σy,σz,ρ]) # ρ -> [<x>,<y>,<z>,<ρ>]

# ╔═╡ 3e789d78-38f8-41c0-b2ae-b230849be534
function blochs(sol)
	(tt, ρt, _) = sol

	# Get Bloch components
	evs0 = expects.(ρt);
	xx,yy,zz,ρρ = [map(x -> x[i], evs0) for i in 1:4];

	(collect(tt), xx, yy, zz, ρρ)
	
end

# ╔═╡ 3cc5e5ea-7ad4-4dc3-b478-06f5ada3cb1a
begin
	ρ0 = dm((spindown(q) + spinup(q))/√2) # initial state
	Δt = 1e-3  # integration time-step
	
	Random.seed!(1)
	solss = bayesian((0, 4τ), ρ0, H0ss, J0ss, C0ss; dt=Δt)
	
	(tη0,xη0,yη0,zη0,ρη0) = blochs(solss)
	
	md"🔻 simulate η = 0"
end

# ╔═╡ f45d3d31-decc-4517-9bfd-494604299d75
begin
	Random.seed!(1)
	sol_ens = ensemble(bayesian, (0,4τ), ρ0, Hss, Jss, Css; dt=Δt, N=10)
	
	tens = sol_ens[1]
	evs = mean(collect(map(ρs -> expects.(ρs), sol_ens[2])))
	(xens,yens,zens,ρens) = [map(x -> x[i], evs) for i in 1:4]
	md" 🔻 Simulate ensemble"
end

# ╔═╡ 8a9c8844-0554-4992-bf92-016f46c0e2df
begin
	Random.seed!(1)
	sol = bayesian((0, 4τ), ρ0, H0, J0, C0; dt=Δt, heterodyne=true)
	
	(ttη0,xxη0,yyη0,zzη0,ρρη0) = blochs(solss)
	
	md"🔻 simulate η = 0"
end

# ╔═╡ 457bd92a-d80c-47c1-878f-3ea9c7e5511b
purity(x,y,z) = 0.5*(1 + x^2 + y^2 + z^2)

# ╔═╡ 4646ff94-42f3-40e0-9119-23dde83b499f
R(x,y,z) = sqrt(x^2 + y^2 + z^2)

# ╔═╡ 9a4d431a-ead6-48eb-9a18-5d56fd7b56d0
function subseries(rec, T, dt; scale=2)
	ts = collect(range(first(T), last(T), step=dt))
	tts = subselect(real(coarse_grain(ts; n=scale)); n=scale)
	(ti, tf) = (tts[1], tts[end])
	dtt = dt * scale
	
	subrec = subselect(real(coarse_grain(rec; n=scale)); n=scale)
	
	(tts, subrec)
	
end

# ╔═╡ e4306041-7bfe-4363-ba25-60556a03d4d2
# colorscheme
begin
	colorscheme = "Paired"
	cmap = plt.matplotlib.cm.get_cmap(colorscheme)
	colors=collect(map(x -> cmap(x), 0:11))
	md" `colorscheme`"
end

# ╔═╡ 67d7fdb7-145f-42f3-8dec-cb6fc8bb3018
let
	close("all")
	subplot(2, 2, 1)
	
	αp_list = αp.(Δs)
	αm_list = αm.(Δs)
	
	p = plot(Δs, angle.(αp_list), color=colors[2], label=L"\alpha_+")
	plot(Δs, angle.(αm_list), color=colors[4], label=L"\alpha_-")
	plot([Δ,Δ], [-3.0,-0.2], linestyle="dashed", color="red")
	


	ax = gca()
	ax.set_xticks(range(-4,4,step=1))
	ax.grid()
	
    ax.set_xlabel(string(L"$\Delta$", " (MHz)"))
    title(string("Phase: ArcTan( Re ", L"\alpha_\pm", "/ Im ", L"\alpha_\pm", ")"))
    legend()
    gcf()
	
	subplot(2, 2, 2)
	plot(Δs, abs.(αp_list ), color=colors[2], label=L"\alpha_+")
	plot(Δs, abs.(αm_list ), color=colors[4], label=L"\alpha_-")
	plot([Δ,Δ], [0.2,1.0], linestyle="dashed", color="red")
	
	ax = gca()
	ax.set_xticks(range(-4,4,step=1))
	ax.grid()
	
	ax.set_xlabel(string(L"$\Delta$", " (MHz)"))
	
    # xlabel(string(L"$\Delta$", " (MHz)"))
    # ylabel(string("Amplitude ", L"|\alpha_\pm|"), loc="top")
    title(string("Amplitude ", L"|\alpha_\pm|"))
    legend()
    gcf()
	
	subplot(2, 2, 3) 
	plot(real.(αp_list), imag.(αp_list), color="gray")
	plot([real(αp(Δ))], [imag(αp(Δ))], color=colors[2], marker="o", label=L"\alpha_+")
	plot([real(αm(Δ))], [imag(αm(Δ))], color=colors[4], marker="o", label=L"\alpha_-")
	
			tight_layout()
	
	ax = gca()
	# ax.set_xticks(range(-2,2,step=1))
	ax.grid()
	
    xlabel(string("Re", L"\alpha_\pm"))
    ylabel(string("Im", L"\alpha_\pm"))
    title("")
    legend()
    gcf()
	

	
end

# ╔═╡ 57a11094-5c9b-44ad-9636-e729bc407ebe
let
	close("all")
	
	subplot(2, 2, 1)
	
	p = plot(ts, εt, color=colors[8], label="")
	plot([ts[i]], [εt[i]], color=colors[8], marker="o")

	ax = gca()
	# ax.set_xticks(range(-4,4,step=1))
	ax.grid()
	
    xlabel(string(L"$t$", " (μs)"))
	ylabel(string(L"$\epsilon(t)$", " (MHz ?)"))
    gcf()
	
	# # # # # 
	
	subplot(2, 2, 2)
	
	αps = αp_list[1:i]
	αms = αm_list[1:i]
	
	plot(real.(αps), imag.(αps), color=colors[2], linestyle="dashed")
	plot([real(last(αps))], [imag(last(αps))], color=colors[2], marker="o", label=L"\alpha_+^{s.s.}")
	plot(real.(αms), imag.(αms), color=colors[4], linestyle="dashed")
	plot([real(last(αms))], [imag(last(αms))], color=colors[4], marker="o", label=L"\alpha_-^{s.s.}")
	
	tight_layout()
	
	ax = gca()
	ax.set_xticks(range(-1.5,1.5,step=0.5))
	ax.set_yticks(range(-2.4, 0.2, step=0.4))
	ax.grid()
	
	ylim([-2.8,0.2])
	
    xlabel(string("Re", L"\alpha_\pm"))
    ylabel(string("Im", L"\alpha_\pm"))
    title("")
    legend()
    gcf()
	
end

# ╔═╡ 941b53b6-02e3-4cbc-ba8a-6e4b6b764f51
let
	close("all")
	
	subplot(2, 2, 1)
	
	p = plot(ts, ε.(ts), color=colors[8], label="")
	plot([ts[j]], [ε(ts[j])], color=colors[8], marker="o")

	ax = gca()
	# ax.set_xticks(range(-4,4,step=1))
	ax.grid()
	
    xlabel(string(L"$t$", " (μs)"))
	ylabel(string(L"$\epsilon(t)$", " (MHz ?)"))
    gcf()
	
	# # # # # 
	
	subplot(2, 2, 2)
	
	αps1 = αp_list[1:j]
	αms1 = αm_list[1:j]
	αps2 = αplist[1:j]
	αms2 = αmlist[1:j]
	
	plot(real.(αps1), imag.(αps1), color=colors[2], linestyle="dashed")
	plot([real(last(αps1))], [imag(last(αps1))], color=colors[2], marker="o", label=L"\alpha_+^{s.s.}")
	plot(real.(αms1), imag.(αms1), color=colors[4], linestyle="dashed")
	plot([real(last(αms1))], [imag(last(αms1))], color=colors[4], marker="o", label=L"\alpha_-^{s.s.}")
	
	plot(real.(αps2), imag.(αps2), color=colors[6], linestyle="dashed")
	plot([real(last(αps2))], [imag(last(αps2))], color=colors[6], marker="o", label=L"\alpha_+")
	plot(real.(αms2), imag.(αms2), color=colors[8], linestyle="dashed")
	plot([real(last(αms2))], [imag(last(αms2))], color=colors[8], marker="o", label=L"\alpha_-")
	
	tight_layout()
	
	ax = gca()
	ax.set_xticks(range(-1.5,1.5,step=0.5))
	ax.set_yticks(range(-2.4, 0.2, step=0.4))
	ax.grid()
	
	ylim([-2.8,0.2])
	
    xlabel(string("Re", L"\alpha_\pm"))
    ylabel(string("Im", L"\alpha_\pm"))
    title("")
    legend()
    gcf()
	
end

# ╔═╡ edd523bb-7ac3-406a-b9a0-41c0967ccb69
# Plotting
function plot_solution(sol; plot_title="Rabi Oscillation")
	
	close("all")
	
	(tt, ρt, dys) = sol
    
    # Get Bloch components
    evs0 = expects.(ρt);
    xx,yy,zz,ρρ = [map(x -> x[i], evs0) for i in 1:4];
    
    # Plot Bloch components vs. time
    
    p = plot(tt, xx, color=colors[2], label=L"$x$")
    plot(tt, yy, color=colors[4],label=L"$y$")
    ax = gca()
    ax.set_ylim([-1.1,1.1]) 
    plot(tt, zz, color=colors[6], label=L"$z$")
	plot(tt, ρρ, color=colors[8], label=L"Tr $\rho^2$")
    xlabel(L"$t$")
    ylabel("Bloch coordinates")
    title(plot_title)
	legend()
	gcf()

end

# ╔═╡ 5eb98909-4d69-4e4e-a7b0-70a1bc7e8874
plot_solution(solss; plot_title="Monitored Rabi Oscillation")

# ╔═╡ 79803729-7f3c-4473-a84f-1f02b9246e1e
plot_solution(sol; plot_title="Monitored Rabi Oscillation")

# ╔═╡ 9e44da9d-5efb-4d9c-9b06-f648b0f37973
# Plotting
function plot_blochs((tt, blochs); plot_title="Rabi Oscillation")
	close("all")
	
	(xx,yy,zz) = blochs
	pp = map(i -> purity(xx[i], yy[i], zz[i]), 1:length(tt))
    
    # Plot Bloch components vs. time
    
    p = plot(tt, xx, color=colors[2], label=L"$x$")
    plot(tt, yy, color=colors[4],label=L"$y$")
    ax = gca()
    ax.set_ylim([-1.1,1.1]) 
    plot(tt, zz, color=colors[6], label=L"$z$")
	plot(tt, pp, color=colors[8], label=L"Tr $\rho^2$")
    xlabel(L"$t$")
    ylabel("Bloch coordinates")
    title(plot_title)
    legend()
    gcf()
end

# ╔═╡ 1d5608cb-e976-477d-a79d-beb5f9b65722
function record_histograms(records...; plot_title="record histogram", labels=[]::Array, density=false)
	close("all")
	
	μσs = []
	hist_colors = []
	hist_labels = []
	
	for i in 1:length(records)
		label = i > length(labels) ? i : labels[i]
		dys = records[i]
		
		# get mean and std dev for real part
		(μ, σ) = params(fit(Normal, dys))
		push!(μσs, map(p -> round(p, digits=4), (μ, σ)))

		# make histogram		
		n, bins, patches = hist(dys, 50, density=density, 					 									facecolor=colors[2i], alpha=1, label=label)
		push!(hist_colors, colors[2i])
		

	end
	
	# write down (μ, σ) pairs as text boxes
	μσ_strings = map(μσ -> string("(μ, σ) = (", μσ[1], ", ", μσ[2], ")\n"), μσs)
	ax = gca()
	for i in 1:length(μσ_strings)
		str = μσ_strings[i]

		ax.text(0.05, 1 - 0.05i, str, transform=ax.transAxes, fontsize=10,
			verticalalignment="top", color=hist_colors[i])
		
	end
	
	
	legend()
	xlabel("value (arbitrary units)")
	ylabel(density ? "relative frequency" : "frequency")
	title("record histograms")
	gcf()
		
end

# ╔═╡ c2d5c6c4-f321-4399-a01d-dd559867ff9e
# Plotting
function plot_timeseries(tt::Array, series...; plot_title="time series", xlab=L"$t$", ylab="arbitrary units", labels=[]::Array, colorpairs=false, kwargs...)
	close("all")
	ser_colors(i) = colorpairs ? colors[i] : colors[2i]
	
	label(i) = i > length(labels) ? i : labels[i]
    
    # Plot records vs. time
	ser = series[1]
	p = plot(tt, ser, color=ser_colors(1), label=label(1), kwargs...)
	ax = gca()
	
	if length(series) > 1
		
		for i in 2:length(series)
			ser = series[i]
			plot(tt, ser, color=ser_colors(i),label=label(i), kwargs...) 
		end
		
	end
	
    xlabel(xlab)
    ylabel(ylab)
    title(plot_title)
    legend()
    gcf()
	
end

# ╔═╡ 98e5f842-148c-4682-9e36-917b43b6e5d5
# Plotting
function plot_timeseries(ttseries...; plot_title="time series", xlab=L"$t$", ylab="arbitrary units", labels=[]::Array, colorpairs=false)
	close("all")
	ser_colors(i) = colorpairs ? colors[i] : colors[2i]
	
	label(i) = i > length(labels) ? i : labels[i]
    
    # Plot records vs. time
	(tt, ser) = ttseries[1]
	p = plot(tt, ser, color=ser_colors(1), label=label(1))
	ax = gca()
	
	if length(ttseries) > 1
		
		for i in 2:length(ttseries)
			(tt, ser) = ttseries[i]
			plot(tt, ser, color=ser_colors(i),label=label(i)) 
		end
		
	end
	
    xlabel(xlab)
    ylabel(ylab)
    title(plot_title)
    legend()
    gcf()
	
end

# ╔═╡ 05b4b92a-d068-4251-a8fa-99970a27977d
plot_timeseries((tη0, xη0), (tens, xens), (tη0, yη0), (tens, yens), (tη0, zη0), (tens, zens); plot_title="ensemble average comparison", ylab="Bloch coordinates", labels=[L"$x_{η = 0}$", L"$x_{avg}$", L"$y_{η = 0}$", L"$y_{avg}$", L"$z_{η = 0}$", L"$z_{avg}$"], colorpairs=true)

# ╔═╡ 0bf61670-0339-47b7-962c-4f3b71028245
# Plotting
function plot_solutions((sol1,sol2); plot_title="Rabi Oscillation")
    close("all")
    
    tt1 = sol1[1]
    ρt1 = sol1[2]
    tt2 = sol2[1]
    ρt2 = sol2[2]
    
    # Get Bloch components
    evs1 = expects.(ρt1);
    x1,y1,z1,ρ1 = [map(x -> x[i], evs1) for i in 1:4];
    evs2 = expects.(ρt2);
    x2,y2,z2,ρ2 = [map(x -> x[i], evs2) for i in 1:4];
    
    # Plot Bloch components vs. time
	p = plot(tt1, x1, color=colors[2], linewidth=2, label=L"$x_{\eta = 0}$")
    plot(tt1, y1, color=colors[4], linewidth=2, label=L"$y_{\eta = 0}$")
    plot(tt1, z1, color=colors[6], linewidth=2, label=L"$z_{\eta = 0}$")
    plot(tt1, p1, color=colors[8], linewidth=2, label=L"(Tr $\rho^2)_{\eta = 0}$")
    plot(tt2, x2,  color=colors[1], linestyle="dashed", label=L"$x_{avg}$")
    plot(tt2, y2, color=colors[3], linestyle="dashed", label=L"$y_{avg}$")
    plot(tt2, p2, color=colors[7], linewidth=2, linestyle="dashed", label=L"(Tr $\rho^2)_{avg}$")
  
	ax = gca()
	ax.set_ylim([-1.1,1.1]) 
	xlabel(L"$t$")
	ylabel("Bloch coordinates")
	title(plot_title)
	legend()
	gcf()
end

# ╔═╡ c9f9505d-7259-45a5-9dcd-38f0ccf7f984
function plot_evals((tt1, evals); α=0.1, linewidth=1, labels=false)
    xxs,yys,zzs,ρρs = [map(x -> x[i], evals) for i in 1:4];
    if labels
        plot(tt1, xxs, color=colors[2], alpha=α, linewidth=linewidth, label=L"$x$")
        plot(tt1, yys, color=colors[4], alpha=α, linewidth=linewidth, label=L"$y$")
        plot(tt1, zzs, color=colors[6], alpha=α, linewidth=linewidth, label=L"$z$")
        plot(tt1, ρρs, color=colors[8], alpha=α, linewidth=linewidth, label=L"Tr $ρ^2$")
    else
        plot(tt1, xxs, color=colors[2], alpha=α, linewidth=linewidth)
        plot(tt1, yys, color=colors[4], alpha=α, linewidth=linewidth)
        plot(tt1, zzs, color=colors[6], alpha=α, linewidth=linewidth)
        plot(tt1, ρρs, color=colors[8], alpha=α, linewidth=linewidth)
        
    end

end

# ╔═╡ 628e171f-eeb1-4f63-9786-9318917441f8
function plot_ensemble(sol_ens; α=0.1, linewidth=1, labels=false, average=false, n=50)
    close("all")
	tt1 = sol_ens[1]
    evs = collect(map(ρs -> expects.(ρs), sol_ens[2]));

    for i in 1:n
        plot_evals((tt1, evs[i]); α=α, labels=labels, linewidth=linewidth)
    end

    if average
        plot_evals((tt1, mean(evs)), α=1, linewidth=1.5, labels=true)
        title_string = "Trajectories w/ ensemble average"
    else
        plot_evals((tt1, evs[1]), α=1, linewidth=1.5, labels=true)
        title_string = "Trajectories"
    end
    
    ax = gca()
    ax.set_ylim([-1.1,1.1]) 
	legend()

    xlabel(L"$t$")
    ylabel("Bloch coordinates")
    title(title_string)
	
	gcf()
end

# ╔═╡ 28377197-47bc-49cf-81a0-f8e3023bd640
plot_ensemble(sol_ens, average=true, n=10)

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
	tan(md"""Change the value of Δ (resonator-drive detuning) using the slider to see how it affects the coherent states.""")

		
end

# ╔═╡ 7a17cf5d-db67-499c-886f-2c56683719f9
hint(text; title="Hint") = Markdown.MD(Markdown.Admonition("hint", title, [text]))

# ╔═╡ Cell order:
# ╠═4c45fe3e-cd69-11eb-20cd-7bfb98c040cf
# ╟─377a3336-20bd-4baa-a033-af8bbc8668a8
# ╟─db847000-8785-4bcb-9f93-db82f6ead06d
# ╟─b7be6342-01f9-4afd-9250-27f0479ed5c5
# ╟─306f55d1-fd48-4830-ba0a-b54ca5dccf82
# ╟─56f69b9a-af7e-4e50-bf6a-0d4e436a7f95
# ╟─c681b1d2-a32b-4c10-bc80-6447368f13eb
# ╟─d5b3318e-a1a8-48c2-baf7-92b398902aff
# ╟─ce3f9a90-70ed-4d73-b7ac-e61791eb1aa9
# ╠═543b95ab-589b-4523-9bf5-b955d249699e
# ╠═9e816960-3731-4b98-96c8-9aa689bdf1ac
# ╟─83f4d601-33f2-4c3d-92c4-37231d55cae3
# ╟─54e9278d-7757-4af9-968b-c78dfbe04758
# ╟─91e5577b-fe11-4d8e-b37e-813a81621385
# ╟─9c860ffd-7a2d-45df-84a8-1ddcf0f46298
# ╟─897d6364-90ee-41cc-ab70-d5ef0c537bd2
# ╟─67d7fdb7-145f-42f3-8dec-cb6fc8bb3018
# ╟─1c0539b3-e6ca-4869-aa65-e268c6842d22
# ╟─e73546a2-2b9e-4062-89f2-5a0913209988
# ╟─05def7cc-6c0e-496e-9429-1fd457575b31
# ╟─38702e03-a6ff-4c7a-bbb6-4dd0bcf89ba8
# ╟─18c5a83d-cf4f-4303-842d-c800193e52b1
# ╟─c728af0d-72e3-4596-9f27-01f82fe86437
# ╟─11838257-1e11-484d-831e-10eeb28abb6d
# ╟─a86b6204-1a3a-45e8-96e6-21dd9369f059
# ╟─0f1146c2-2675-492f-883d-d6a353aab727
# ╟─57a11094-5c9b-44ad-9636-e729bc407ebe
# ╟─ca2e2708-789b-4a95-ba3f-0402f93cb84e
# ╟─f5719cfb-1acf-4f87-b12e-a0ed65f635d5
# ╟─e16b9185-418e-4846-96dc-c076ad781e70
# ╟─cf3fd908-d968-4581-a8b4-5aac46c4afb0
# ╟─c5141459-60d0-47ed-90f3-ac4e741e5707
# ╟─0d85d373-e094-4864-8a1a-83936c3362d5
# ╟─71ea6e65-46b3-4c6c-a919-87dd31a6674f
# ╟─5fe7e35f-c2d1-425d-96b4-bc3e2c377f61
# ╟─2b5dc7be-d3c8-40b8-b13a-d36b51f114d8
# ╟─e0db1833-9b68-4254-8556-c76d0162a3d1
# ╟─941b53b6-02e3-4cbc-ba8a-6e4b6b764f51
# ╟─3ca15002-9b64-416c-94e3-a7904766a87b
# ╟─a9dc195d-fe8f-4c6b-8078-597585b2ee3f
# ╟─17d6c0a1-9191-4da8-a298-e63fd4008453
# ╟─14cccb18-0187-4c91-ba53-b920ce0c3704
# ╟─e8299f27-e74b-47d1-8349-7e28c3ea4495
# ╟─ec64b4d3-7631-4a1e-af3e-cf895f2fd922
# ╟─4e368417-1107-45af-9af8-10fc16eaabe5
# ╠═a3ec8af8-9f28-468a-8d9f-286c66e69164
# ╠═61e73b8f-1ee6-4a5c-aa11-beb85bdd80d1
# ╠═0ec0004b-92dc-451b-bc74-54a8d44d85a5
# ╠═5a86ea71-0684-476f-8d5d-ee065b43c92f
# ╠═682b7883-fb96-4ebb-a196-5d85c5f91418
# ╠═3cc5e5ea-7ad4-4dc3-b478-06f5ada3cb1a
# ╠═5eb98909-4d69-4e4e-a7b0-70a1bc7e8874
# ╠═f45d3d31-decc-4517-9bfd-494604299d75
# ╠═05b4b92a-d068-4251-a8fa-99970a27977d
# ╠═28377197-47bc-49cf-81a0-f8e3023bd640
# ╟─ebfe6948-64be-4190-8470-d5e196ee43ec
# ╠═67c05b7c-4376-495a-8a21-b5d4f749c0d6
# ╟─2b8f2e1c-0f25-4be0-ab40-2bb0d96ba5e8
# ╟─e4269abf-9b63-4fbe-8001-397a449dd0e6
# ╠═e1345e44-5104-42bb-81a6-1e85db1ccf91
# ╠═7780e6ad-5241-4bde-bf43-a7410d4c1650
# ╠═202fde97-a8ed-44d9-a02e-6a7f564b348d
# ╠═07629b53-a846-4341-9333-0c442ad2fc96
# ╠═8d353b07-2dc9-4611-bb94-fd413dabf120
# ╠═64da2aa4-96ac-49b8-8028-4dad4698353e
# ╠═b464fa05-4b99-49a7-838b-93826da4f0e6
# ╠═dbaaafbb-18b0-43ff-9b20-8d337ca745a9
# ╠═c2085fde-2ec3-4d27-8aca-9b6ecc8218b1
# ╠═c0e189dc-26cc-4bbe-9b72-aa681dd57e46
# ╠═dc9478a2-0532-4b18-bd6e-2a0942eef585
# ╠═ec327412-51bd-47f8-b40e-b89b4679dc90
# ╠═cd1f1146-0b21-4b01-89de-e8bf67ed7867
# ╠═8a9c8844-0554-4992-bf92-016f46c0e2df
# ╠═79803729-7f3c-4473-a84f-1f02b9246e1e
# ╟─02d2a1f6-12af-4345-b42a-db2b11e2a883
# ╟─498bfa13-4278-4808-b1d6-279eabcd9ed8
# ╟─c7eaf3d6-73ac-4291-ab9c-5857cb77cdaf
# ╟─7f1176b8-f6d3-4fd1-a23e-31710dcfff10
# ╟─982565bc-b137-4c10-8e13-5fe37be86823
# ╠═4284173a-be05-4b58-a8d9-7189301344fd
# ╟─e51ca599-2452-4bfb-9f64-f27751d80aa0
# ╟─7e62e6c6-aa3c-4350-901c-15d017b8db42
# ╟─d7ff20d4-674d-4248-a493-f63bff94aefe
# ╟─f3c97871-bcf3-45f2-aae1-2c5d72461d4d
# ╟─3e789d78-38f8-41c0-b2ae-b230849be534
# ╟─3274ad19-1b27-4090-8aee-ef8c148b5e8e
# ╟─457bd92a-d80c-47c1-878f-3ea9c7e5511b
# ╟─4646ff94-42f3-40e0-9119-23dde83b499f
# ╟─edd523bb-7ac3-406a-b9a0-41c0967ccb69
# ╟─9e44da9d-5efb-4d9c-9b06-f648b0f37973
# ╟─1d5608cb-e976-477d-a79d-beb5f9b65722
# ╟─9a4d431a-ead6-48eb-9a18-5d56fd7b56d0
# ╟─c2d5c6c4-f321-4399-a01d-dd559867ff9e
# ╟─98e5f842-148c-4682-9e36-917b43b6e5d5
# ╟─0bf61670-0339-47b7-962c-4f3b71028245
# ╠═c9f9505d-7259-45a5-9dcd-38f0ccf7f984
# ╠═628e171f-eeb1-4f63-9786-9318917441f8
# ╟─e4306041-7bfe-4363-ba25-60556a03d4d2
# ╟─d304e60c-f47f-4b75-b34d-4ef5dc7b5750
# ╟─a7332fed-6089-4ffd-a93e-8893024d87e8
# ╟─dfb8cf89-86f5-4ab9-aac3-5867130fd6bd
# ╟─b01d4d05-7b32-4bd4-b554-e85abca61595
# ╟─7a17cf5d-db67-499c-886f-2c56683719f9
