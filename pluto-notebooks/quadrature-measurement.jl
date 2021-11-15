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
	using Random
	using Statistics
	using Distributions
	using PyPlot
	using QuantumCircuits
end

# ╔═╡ f230c876-f933-4ee4-a079-a66147dea73c
using PlutoUI

# ╔═╡ 377a3336-20bd-4baa-a033-af8bbc8668a8
md"""
# Quadrature measurements

In this interactive notebook, we'll look at how changing measurement angle affects qubit measurement.
"""

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents()

# ╔═╡ f24df928-22c4-4943-bc9c-218750fe4da7
md" ## Readout"

# ╔═╡ 3d17c4b0-e051-41b0-9f20-bb22c4db4398
md"""
Readout is a process with different steps at room temperature, in the cavity, and at the qubit. In this section, we'll walk backward through the experimental setup. We'll start from demodulating the readout at room temperature, look back at what the readout says about the resonator, then finally predict the backaction on the qubit.
"""

# ╔═╡ 3500ae58-7577-4af7-a458-5613310471d9
md"""
## Room temp: Demodulation of readout signal
"""

# ╔═╡ 4e2680a7-8035-4161-b4e3-abb3f54da99d
md" ### Background"

# ╔═╡ d330e502-e4d9-482d-949e-d9c3d0166b01
md"""
In superconducting qubit readout, the readout signal "RO" reflects from the cavity at cryogenic temperatures and travels along the transmission line. After amplification, it is demodulated at room temperature using an IQ mixer as depicted in the figure below.

"""

# ╔═╡ 8c3e20ca-4018-4448-8848-073a2bbcf3bb
md"""

$(LocalResource("/Users/sachagreenfield/Desktop/Pluto-figures/demodulation-single-quadrature.png"))"""

# ╔═╡ 54690fb9-7f41-4166-b231-d0915ef293d7
md"""
The figure depicts a readout signal $s(t) = I \cos(\omega_\text{RO} t) + Q \sin(\omega_\text{RO} t) $ that has exited the resonator along the transmission line, and has been amplified. It enters the 3-port mixer and is mixed with an LO signal $B \cos(\omega_\text{LO} t + \varphi)$, which has been phased-shifted with some phase $\varphi \in [0,2\pi)$. The LO is very close to the readout carrier frequency, compared to characteristic frequencies: $\omega_\text{LO} - \omega_\text{RO} \ll \omega_\text{LO}$. The mixer acts to:
"""

# ╔═╡ 2a548d96-451f-4bd7-b9e5-7c3fa85c2abb
md"""
**MULTIPLY the SIGNALS**, creating an output signal

$B \cos(\omega_\text{LO} t + \varphi) s(t) = B I \Big[ \cos \big( (\omega_\text{LO} + \omega_\text{RO}) t + \varphi \big) + \cos \big( (\omega_\text{LO} - \omega_\text{RO}) t - \varphi \big)\Big] + B Q \Big[ \sin \big( (\omega_\text{LO} + \omega_\text{RO}) t + \varphi \big) - \sin \big( (\omega_\text{LO} - \omega_\text{RO}) t - \varphi \big)\Big]$

according to sum-difference formulas.

**LOW-PASS FILTER the RESULT.** Since $\omega_\text{LO} \approx \omega_\text{RO}$, low-pass filtering outputs

$B \big[I \cos(\Delta_\text{RO} t + \varphi) + Q \sin(\Delta_\text{RO} t + \varphi) \big]$

along the IF port, where $\Delta_\text{RO} \equiv \omega_\text{LO} - \omega_\text{RO}$. 

"""

# ╔═╡ f1d33aff-03d3-4223-b67d-d7e9be26c168
md" ### Mixer simulation"

# ╔═╡ 019c8da1-8cc8-452b-839f-bc50772f9d79
begin
	χ = 2π * (-0.47)
	κ = 2π * (1.56)
end

# ╔═╡ 1abcbbed-e7e6-4a06-a1eb-3a0922785323
md"""
**Show:**
$(@bind show_RO CheckBox(default=true)) readout 
$(@bind show_LO CheckBox(default=true)) LO 
$(@bind show_LORO CheckBox(default=false)) LO x readout
$(@bind show_int CheckBox(default=false)) integrated signal
$(@bind show_env CheckBox(default=false)) envelope 
"""

# ╔═╡ 970d9217-e5e9-4adb-8ada-5d9a225523c3
md"""

Integrate the mixed signals (LO x readout) across `nc` bins to low-pass filter the result: 

`nc = ` 0
$(@bind nc html"<input type=range min=1 max=25 step=1 value=1>") 
25 

Decrease the detuning `ΔRO` between LO and RO to switch from heterodyne to homodyne measurement (discussed in next section):

`ΔRO =` 400 rad MHz
$(@bind ΔROm html"<input type=range min=0 max=400 step=20 value=0>")  0 rad MHz 

"""

# ╔═╡ e0eda8fc-2ff6-49b0-82b3-156954efedbe
ΔRO = 400 - ΔROm

# ╔═╡ 82eed712-4adc-4408-baaa-64f49373c405
md"""
Change the time-window displayed: `tf =` $(@bind tfstr Select(["2", "4", "8", "20", "40"])) ns
"""

# ╔═╡ fc688b0e-7de9-473f-9c15-1e8492944f28
md"""
### Types of measurement
"""

# ╔═╡ 38b793c7-75c1-4b26-a6de-cd213db0dd2a
md"""
The demodulation via mixers can be configured in several different ways, leading to distinct types of measurement.
"""

# ╔═╡ e92eca0d-854f-4257-9e54-ed7acf79ecff
md"""
###### HOMODYNE MEASUREMENT
It is possible to choose $\Delta_\text{RO} =0$. This results in an output signal

$o(t) = B \big[I \cos \varphi + Q \sin \varphi \big].$

Thus, by demodulating the readout signal using an LO of the same frequency, the output is brought brought down to DC. This is called a *homodyne* measurement, because there is one ("homo") frequency ("dyne", or power) involved.
"""

# ╔═╡ 367819f8-b540-445f-83a4-a9c7982e878c
md"""
###### HETERODYNE MEASUREMENT

Alternatively, the signal can be measured with $\Delta_\text{RO} \neq 0$, and the output signal will be slow oscillating with frequency $\Delta_\text{RO}$. The I and Q amplitudes can then be extracted via digital demodulation. Further discussion follows at the end of this section.

"""

# ╔═╡ 2f7a96c1-d043-407f-84f8-68efea7ded0e
md"""
###### SINGLE-QUADRATURE MEASUREMENT

Processing the signal on a 3-port mixer, as in the first diagram, leads to single-quadrature measurement. Notice that choosing $\varphi = 0$ leads to measurement of an output $o(t) \propto I$, while choosing $\varphi = \pi/2$ leads to $o(t) \propto Q$. Thus, treating $I$ and $Q$ as real and imaginary parts of a complex signal $\alpha = I + i Q$, we have

$o(t) \propto \text{Re} (e^{-i \varphi}\alpha).$

We measuring $\alpha$ along a new effective quadrature determined by $\varphi$. Thus we call $\varphi$ the "measurement angle." It will have important consequences for the type of information about the qubit acquired as well as the direction of measurement backaction.

It is interesting to note that this choice is made at *room temperature*, long after the readout pulse has left the fridge.

"""

# ╔═╡ c2f7d614-0d49-4283-b574-c0713e41824a
md"""
###### BALANCED (DUAL-QUADRATURE) MEASUREMENT
It is also possible to measure two orthogonal quadratures of the output by means of a 4-port mixer. This will lead to effective outputs

$o_\varphi(t) \propto \text{Re} (e^{i \varphi} \alpha), \hspace{5mm} o_{(\varphi + \pi/2)}(t) \propto \text{Re} (e^{i (\varphi + \pi/2)} \alpha).$

For example, when $\varphi = 0$, this becomes

$o_\varphi(t) \propto I, \hspace{5mm} o_{(\varphi + \pi/2)}(t) \propto Q.$

This measures both quadratures of the qubit and will lead to both informational and phase backaction.

See Ref. [1] Sec. V for more details on experimental implementation of homodyne and heterodyne measurement.


"""

# ╔═╡ b08e7d51-3e2c-431a-8471-12914bcd1384
md"""
###### NOTE ON HOMODYNE AND HETERODYNE DETECTION

In practice, homodyne measurement is practical only for readout of a single qubit. Heterodyne measurement must be used to distinguish signals of multiple qubits. In such a case, the resonators corresponding to qubits $1, 2, ..., n$ will have detuned frequencies $\omega_1, \omega_2, ..., \omega_n$, and the input signal will take the form

$s(t) = \sum_{i = 1}^n s_i(t)$

where 

$s_i(t) = I_i \cos \omega_i t + Q_i \sin \omega_i t$

and the output signal after low-pass filtering is

$B \sum_{i = 1}^n \big[I_i \cos(\Delta_i t + \varphi) + Q_i \sin(\Delta_i t + \varphi) \big],$

with $\Delta_i \equiv \omega_\text{LO} - \omega_i$.

Because $\Delta_i$ are in the intermediate frequency range, they can be demodulated digitally and the signals from multiple qubits processed, so long as the $\omega_i$ are sufficiently separated in frequency space.


"""

# ╔═╡ 3bc48dd3-363b-4dc8-8f9d-b3edc7700bce


# ╔═╡ e4f0d549-b9d2-48ef-bbe7-cd4731156c77
md"""
###### TERMINOLOGY OF HOMODYNE AND HETERODYNE

In some sources (e.g. [2], [3]), "homodyne" and "heterodyne" distinguish the number of quadratures detected for a single qubit, with "homodyne" indicating single-quadrature measurement and "heterodyne" indicating dual-quadrature measurement.

However, in the signal processing literature, "homodyne" and "heterodyne" distinguish the detuning of the LO from the signal frequency, as discussed above, with "homodyne" ("heterodyne") referring to (non-)zero detuning. Using this terminology, it is possible to make single- or dual-quadrature measurements using *either* homodyne *or* heterodyne detection. The difference lies only in whether the signal is brought down to DC analogically (by the mixer, in homodyne) or digitally (by the computer, in heterodyne).

This latter terminology will be used in these tutorials.
"""

# ╔═╡ 58d845a6-45e7-450e-b6c2-009a605e1c7c
md" ## Resonator: Cavity states and phasor diagrams"

# ╔═╡ 5702d0d1-6e34-469e-926b-c369dbde78b3
md" ### Background"

# ╔═╡ d44ea114-6e92-4790-978d-7286394e9506
md"""
More detailed discussion information is contained in `coherent-states.jl`. The phasor $\alpha$ corresponds to a coherent state amplitude that is entangled with the qubit and travels down the transmission line before being mixed with the LO.

In this discussion, the cavity amplitudes $\alpha_+$ and $\alpha_-$ are phase-shifted by $i$ relative to those in `coherent-states.jl`. I believe this is due to a phase shift on reflection of the RO tone from the cavity. I have included it to be consistent with the literature.
"""

# ╔═╡ 8a8987ef-225c-4538-92d0-ca2004266053
md" ### Cavity amplitudes simulation"

# ╔═╡ 82ffdfde-b2ff-4859-8ca6-119a4d5d95cd
begin
	ε0 = 1 # MHz
	Δs = range(-12, 12, step=0.2) # MHz
	md" 🌀 parameters"
end

# ╔═╡ bddc09b0-f90a-40e5-be77-d1afa4309fa4
begin
	αpss(Δ) = (2ε0/κ) / (1 + im * (2(Δ + χ)/κ))
	αmss(Δ) = (2ε0/κ) / (1 + im * (2(Δ - χ)/κ))
	md" 🌀 function definitions"
end

# ╔═╡ 6ae79dd2-dc81-4d00-8bce-394d788a9207
begin
	I = real(αpss(0))
	Q = imag(αpss(0))
	ωRO = 2π * (6.67e3) # rad MHz
	ωLO = ωRO + ΔRO # 2π * (6.60e3) # rad MHz
	B = 1
	ts = range(0, 0.04, step=2π/ωRO / 50)
	φ0 = 0
	
	md"🌀 define parameters"
end

# ╔═╡ 3b125df9-c86d-4e4a-949e-8701b91d58e4
begin
	s(t) = I * cos(ωRO * t) + Q * sin(ωRO * t)
	lo(t) = B * cos(ωLO * t + φ0)
	prod(t) = B * (I * cos(-ΔRO * t + φ0) + Q * sin(-ΔRO * t + φ0))
	
	md" 🌀 define signals"
end

# ╔═╡ ba89c6ca-abd3-48ce-b1da-54f282fbd94d
begin
	ωs = range(ωRO - 400, ωRO + 400, step=4)
	
	ωROs = zeros(length(ωs))
	ωROs[argmin(abs.(ωs .- ωRO))] = sqrt(I^2 + Q^2)
	
	ωLOs = zeros(length(ωs))
	ωLOs[argmin(abs.(ωs .- ωLO))] = B
	
	md" 🌀 define frequency domain series"
end

# ╔═╡ 4652677b-8602-4640-836b-07fab787e88c
let 
	tf = parse(Float64, tfstr)
	newsig = coarse_grain(s.(ts) .* lo.(ts); n=nc)
	
	close("all")
	
	# Signal time series ---------------------------------------------------------
	
	subplot(2,1,1)
	tss = ts * 1e3
	# tff = tf * 1e3
   	if show_RO 
		p = plot(tss, s.(ts), color="gray", linewidth=1, alpha=0.9, label="readout")
	end
	
	if show_LO
		plot(tss, lo.(ts) ./ 5, color="purple", linewidth=1, alpha=0.9, label="LO / 5")
	end
	
	if show_LORO
		plot(tss, s.(ts) .* lo.(ts), color="green", linewidth=1, alpha=0.9, label="LO x readout") end
	
	if show_int
		plot(tss, newsig, color="blue", linewidth=1, alpha=0.9, label="integrated")
	end
	
	if show_env
		plot(tss, 0.5* prod.(ts), color="black", linewidth=1, alpha=0.9, label="envelope", linestyle="dashed") end

	ax1 = gca()
	ax1.set_xlim([0, tf]) 
	ax1.set_ylim([-0.22, 0.22]) 
	
    xlabel("t (ns)")
    ylabel("arbitrary units")
    title("signal mixing")
	ax1.legend(loc="lower right")
	
    gcf()
	
	
	
	# Frequency domain ---------------------------------------------------------
	
	subplot(2,1,2)
   	p = plot(ωs, ωROs, color="gray") 
	p = plot(ωs, ωLOs, color="purple") 

	ax2 = gca()
	# ax2.set_xlim([0, tff]) 
	
    xlabel("ω (rad MHz)")
    ylabel("arbitrary units")
    title("frequency domain")
	ax2.legend(loc="lower right")
	tight_layout()
	
    gcf()
	
end

# ╔═╡ d6f73306-ab64-4fba-9955-55298ab566e5
md"""

Δ : -12
$(@bind Δ html"<input type=range min=-12 max=12 step=0.25 value=0>")
12

"""

# ╔═╡ aec1af64-50c5-4560-8461-e45bfc2f406d
md"""
Δ = $Δ
"""

# ╔═╡ b4cbaf81-c33c-42a9-ba79-437c32f47f70
md" ## Qubit: Measurement backaction"

# ╔═╡ de3fe955-780c-489f-9b6d-f3a8c3bf140d
md" # Weak measurement demos"

# ╔═╡ cd661d56-ebd0-4334-9ae6-47e738d11b54
md" ### Reduced-qubit description"

# ╔═╡ 88e1d88a-3565-42b3-9732-02862163a12b
md"""
In the reduced-qubit evolution, the measurement angle simply determines the quadrature of the backaction.
"""

# ╔═╡ 4284173a-be05-4b58-a8d9-7189301344fd
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

# ╔═╡ ca554621-0fe8-4d7a-bd3f-acdf795648c4
md"""

ϕ = 0 
$(@bind ϕc html"<input type=range min=0 max=16 step=1 value=0>") 
ϕ = 2π

"""

# ╔═╡ fe1986ff-00c7-4b9f-8d20-9a38c85e2a28
begin
	ϕ = ϕc * (π/8)
	md" ϕ = $(ϕc/8) π "	
end

# ╔═╡ 073c34d3-a243-4560-9862-28f15b24e685
md"""
$\ket{\psi} = cos(\theta /2) \ket{0} + \sin(\theta/2) e^{i \phi} \ket{1}$
"""

# ╔═╡ de277991-5f34-4906-a9bb-b1832131b45f
md"""

ϕ = 0 
$(@bind ϕcd html"<input type=range min=0 max=16 step=1 value=0>") 
ϕ = 2π

"""

# ╔═╡ eb1ac4e6-5bc3-4c67-85b0-e4be923362a4
begin
	ϕd = ϕcd * (π/8)
	md" ϕ = $(ϕcd/8) π "	
end

# ╔═╡ ede93004-ef9e-461f-b4ea-062f9d7879f3
md"""
### Qubit-resonator description
"""

# ╔═╡ 9d0fe6ee-aa00-43c8-b40f-937dd6024d74
Nfock = 15 # fock space dimension cutoff

# ╔═╡ 05b536fa-e464-4476-8c27-65ffe89d09a6
begin
	# Basis
	f = FockBasis(Nfock)
	If = identityoperator(f)

	# Qubit-resonator operators
	Id = Iq ⊗ If
	a = Iq ⊗ destroy(f)
	
	σx = σxq ⊗ If
	σy = σyq ⊗ If
	σz = σzq ⊗ If
	σp = σpq ⊗ If
	σm = σmq ⊗ If
	
	# projectors
	zp = dm(excited) ⊗ If
	zm = dm(ground) ⊗ If
	
	αp = a * zp
	αm = a * zm
	
	md"###### 🔶 Qubit-resonator Hilbert space operators "
end

# ╔═╡ aaf95d04-aeb2-47b4-99c6-9499fe27a33f
md"""

t : $(first(tt2)) μs
$(@bind j html"<input type=range min=1 max=12001 step=10 value=1>")
$(last(tt2)) μs

"""

# ╔═╡ fcc6db99-6ff7-41d5-b2af-ff07899a0250
md" t = $(tt2[j]) μs "

# ╔═╡ c1e0ad96-5335-4f11-abe0-5b4f5274c36e
md"""
# References

[1] P. Krantz, M. Kjaergaard, F. Yan, T. P. Orlando, S. Gustavsson, and W. D. Oliver, A Quantum Engineer’s Guide to Superconducting Qubits, Applied Physics Reviews 6, 021318 (2019).

[2] P. Campagne-Ibarcq, Measurement Back Action and Feedback in Superconducting Circuits, 223 (n.d.).

[3] J. Gambetta, A. Blais, M. Boissonneault, A. A. Houck, D. I. Schuster, and S. M. Girvin, Quantum Trajectory Approach to Circuit QED: Quantum Jumps and the Zeno Effect, Phys. Rev. A 77, 012112 (2008).



"""

# ╔═╡ 7e62e6c6-aa3c-4350-901c-15d017b8db42
md" # Utilities "

# ╔═╡ 989ba897-50b4-4d82-be56-6dcb439e8aca


# ╔═╡ 3d1472b2-30ac-4237-ace3-00ad5f5768b6
xyz(θ, ϕ) = (sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ))

# ╔═╡ 76d5307b-f58b-4ab1-8ba5-b43491c54552
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
		σs = (size(ρ[1]) == (2, 2)) ? (σxq, σyq, σzq) : (σx, σy, σz)
		x, y, z =  [real(expect(σi, ρ)) for σi in σs] 
		p = real(expect.(ρ, ρ))
		traj(t, x, y, z, p, r)
	end
	
	
end

# ╔═╡ 139565f9-0219-4939-a470-4252397dd6a0
let
	
	# Parameters -------------------------------------------------------------------
	
	ρ0 = DenseOperator(0.5*(Iq + σxq)) # dm(spindown(q)) # initial state
	dt = 1e-3  # integration time-step
	
	ΩR  = 0 #2π # Rabi frequency (rad * MHz)
	Γ = 0.15 # Measurement dephasing rate (MHz)
	τ = 1/(2Γ) #  Measurement collapse timescale (μs)
	η = 1 # collection efficiency
	
	# Kraus operators --------------------------------------------------------------
	
	H = ΩR * σyq/2
	J = [(σzq, ((1-η)*Γ))]
	C = [(exp(im * ϕ) * σzq, τ, η)]
	# Cd = [(exp(im * ϕd) * σzq, τ0, η0/2), (exp(im * (ϕd + π/2)) * σzq, τ0, η0/2)]
	
	Random.seed!(1)
	sol = bayesian((0, 4τ), ρ0, H, J, C; dt=dt)
	
	global RQ = traj(sol)
	

	
	md" ###### 🔻 Bayesian simulation (single-quadrature)"
end

# ╔═╡ b6ae1b65-aa9d-4d7e-b07e-f51ebdcd44b3
let
	
	# Parameters -------------------------------------------------------------------
	
	ρ0 = dm(spindown(q)) # initial state
	dt = 1e-3  # integration time-step
	
	ΩR  = 2π # Rabi frequency (rad * MHz)
	Γ = 0.15 # Measurement dephasing rate (MHz)
	τ = 1/(2Γ) #  Measurement collapse timescale (μs)
	η = 1 # collection efficiency
	
	# Kraus operators --------------------------------------------------------------
	
	H = ΩR * σyq/2
	J = [(σzq, ((1-η)*Γ))]
	C = [(exp(im * ϕd) * σzq, τ, η/2), (exp(im * (ϕd + π/2)) * σzq, τ, η/2)]
	
	Random.seed!(1)
	sol = bayesian((0, 4τ), ρ0, H, J, C; dt=dt)
	
	global RQ2 = traj(sol)
	

	
	md" ###### 🔻 Bayesian simulation (single-quadrature)"
end

# ╔═╡ a823cf8d-8637-44b9-af2a-8feefd22d986
begin
	
	mutable struct qr_traj
		t::Vector{Float64}
		
		x::Vector{Float64}
		y::Vector{Float64}
		z::Vector{Float64}
		
		n::Vector{Float64}
		ap::Vector{ComplexF64}
		am::Vector{ComplexF64}
		
		p::Vector{Float64}
		pr::Vector{Float64}
		r
	end
	
	
	function qr_traj(sol::QuantumCircuits.solution)
		t, ρ, r = (sol.t, sol.ρ, sol.r)
		
		# get expectation values
		x, y, z =  [real(expect(σi, ρ)) for σi in (σx, σy, σz)] 
		n, zpp, zmm = [real(expect(op, ρ)) for op in ((a' * a), zp, zm)] 
		αpp, αmm = [expect(op, ρ) for op in (αp, αm)] 
		p = real(expect.(ρ, ρ))

		# calculate functions of exp. values
		pq = 0.5 .* (1 .+ x.^2 .+ y.^2 .+ z.^2)
		ap = αpp ./ zpp
		am = αmm ./ zmm
		
		qr_traj(t, x, y, z, n, ap, am, p, pq, r)
	end
	
end

# ╔═╡ 55b946e9-20e2-4ea8-85e2-455dd9b6b8a1
let
	# Parameters -------------------------------------------------------------------
	
	Δt = 1e-3 # integration time step
	
	# carrier frequencies - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
# 	ωr = 2π * (6.67e3) # rad MHz; bare resonator frequency
# 	ωε = 2π * (6.67e3)# rad MHz; resonator readout pulse carrier frequency
# 	ωq = 2π * (5.56e3) # rad MHz; bare qubit frequency
#     ωR = 2π * (5.56e3) # 2π*(5.559871e3) # rad MHz; qubit drive carrier frequency
	
	
	# detunings - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	Δc = 0 # ωr - ωε # cavity readout pulse from bare cavity frequency
	Δq = 0 # ωq - ωR # Rabi drive from bare qubit frequency
	
	
	# envelopes - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	ΩR = 2π * (0.4) # rad MHz; Rabi drive
	ε = 2π * (0.7) # rad MHz; cavity readout drive envelope
	χ = 2π * (-0.47) # rad MHz; dispersive shift / 2
	
	
	# rates - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	T1 = 160 # μs; qubit energy decay time 
	γ1 = 1/T1 # qubit energy decay rate
	
	Tϕ = 16 # μs; cavity-induced dephasing time
	γϕ = 1/Tϕ # cavity-induced dephasing rate
	
	κ = 2π * (1.56) # rad Hz; cavity linewidth / decay rate
	τ = 1/(2κ) # measurement collapse time
	
	
	# measurement parameters - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	θR = π/2 # Rabi-drive axis angle in x-y plane, w.r.t. x-axis
	φ = 0 # measurement quadrature
	η = 1 #1 # signal collection efficiency
	
	
	# initial state  - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	θ0 = π/2 # initial qubit polar angle
	ϕ0 = 0 # initial qubit azimuthal angle
	
	ρq = dm(normalize(cos(θ0/2) * ground + exp(im * ϕ0) * sin(θ0/2) * excited))
	ρi = ρq ⊗ dm(fockstate(f, 0)) # initial qubit in ρq and cavity in vacuum
	
	
	# Kraus operators -----------------------------------------------------------
	
	# Hamiltonian - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	Hc = Δc * (a' * a) # cavity 
	Hq =  Δq * σz / 2 # qubit
	Hqc = χ * (a' * a) * σz # cavity-qubit coupling
	Hε = (ε / 2) * (a + a') # cavity (readout) drive
	HR = (ΩR / 2) * (cos(θR) * σx + sin(θR) * σy) # qubit (Rabi) drive
	
	H = Hε + Hqc  # + HR # + Hc + Hq 
	
	
	# Lindblad & jump operators - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	J =	[(a, κ)]
	# J = [(a, κ), (σm, γ1), (σz, γϕ/2)] # including T1, T2
	
	# single-quadrature measurement (phase-amplifying)
	C = [(exp(im * φ) * a, τ, η)] 	
	
	# double-quadrature measurement (phase-preserving) 
	# Cd = [(exp(im * φ) * a, τ, η/2), (exp(im * (φ + π/2)) * a, τ, η/2)] 
	
	# Simulation ----------------------------------------------------------------
	
	Random.seed!(1)
	sol = bayesian((0, 10), ρi, H, J, C; dt=Δt)
	
	global QR = qr_traj(sol)
	
	# # get expectation values
	# funcs = [σx, σy, σz, (a' * a), zp, zm]
	# funcsC = [αp, αm]
	# evs2 = expects(funcs).(ρt2)
	# evs2C = expectsC(funcsC).(ρt2)
	# xx2,yy2,zz2,nn2,zp2,zm2 = [map(x -> x[i], evs2) for i in 1:length(funcs)];
	# αpex2,αmex2 = [map(x -> x[i], evs2C) for i in 1:length(funcsC)];
	
# 	# calculate functions of exp. values
# 	ρρ2 = 0.5 * (1 .+ xx2.^2 .+ yy2.^2 .+ zz2.^2)
# 	αp2 = αpex2 ./ zp2
# 	αm2 = αmex2 ./ zm2
	
	md" ###### 🔻 Bayesian simulation (single-quadrature)"
end

# ╔═╡ 865752f3-c563-4a6c-be2f-3dbe73d38406
md"""

t : $(first(QR.t)) μs
$(@bind i html"<input type=range min=1 max=10010 step=10 value=10010>")
$(last(QR.t)) μs

"""

# ╔═╡ 9baf661a-feb9-463d-a75e-69848ad888ae
md" t = $(QR.t[i]) μs "

# ╔═╡ 3604c61d-6cf4-4d6a-94ea-678c196d336e
typeof(1.0 + im) == ComplexF64

# ╔═╡ 52e2feca-5070-444d-aef6-c85c6404d391
size(σxq) == (2,2)

# ╔═╡ b136e6e0-6e88-4008-94bb-7756716e913e
φdict = Dict("0" => 0, 
				"π/8" => π/8, 
				"π/4" => π/4,
				"3π/8" => 3π/8,
				"π/2" => π/2)

# ╔═╡ dfb8774e-bb14-4f8b-9492-ecd2220e393d
md"""
### Old utilities
"""

# ╔═╡ 71acc211-c8b8-4059-8002-be14fff7d822
expectsC(ops) = ρ -> collect(expect(ρ, s) for s in vcat(ops, ρ)) # ρ -> [<x>,<y>,<z>,<ρ>]

# ╔═╡ 4de027ac-2080-499f-9a99-7ef3007e1384
expects(ops) = ρ -> collect(real(expect(ρ, s)) for s in vcat(ops, ρ)) # ρ -> [<x>,<y>,<z>,<ρ>]

# ╔═╡ 28c3789a-303e-4c08-83d2-95066828a1d7
begin
	Random.seed!(1)
	sol3 = bayesian((0, 12), ρi, H, J, Cd; dt=Δt)
	
	# collect outputs
	tt3 = sol3[1]
    ρt3 = sol3[2]
	rs1 = collect(sol3[3][1])
	rs2 = collect(sol3[3][2])
	
	# get expectation values
	evs3 = expects(funcs).(ρt3)
	evs3C = expectsC(funcsC).(ρt3)
    xx3,yy3,zz3,nn3,zp3,zm3 = [map(x -> x[i], evs3) for i in 1:length(funcs)];
	αpex3,αmex3 = [map(x -> x[i], evs3C) for i in 1:length(funcsC)];
	
	# calculate functions of exp. values
	ρρ3 = 0.5 * (1 .+ xx3.^2 .+ yy3.^2 .+ zz3.^2)
	αp3 = αpex3 ./ zp3
	αm3 = αmex3 ./ zm3
	
	md" ###### 🔻 Bayesian simulation (dual-quadrature)"
end

# ╔═╡ b2235afc-1bc9-4cf4-ae21-6f535eb37e96
# expects(ρ) =  collect(real(expect(ρ, s)) for s in [σx,σy,σz,ρ]) # ρ -> [<x>,<y>,<z>,<ρ>]

# ╔═╡ 0a2219d8-2fbd-4b03-867f-966c4a595e78
# colorscheme
begin
	colorscheme = "Paired"
	cmap = plt.matplotlib.cm.get_cmap(colorscheme)
	colors=collect(map(x -> cmap(x), 0:11))
	md" `colorscheme`"
end

# ╔═╡ 0a6b57b5-d5bf-4598-a01a-eced9358e452
let
	close("all")
	
	αp_list = αpss.(Δs)
	αm_list = αmss.(Δs)
	
	plot(real.(αp_list), imag.(αp_list), color="gray")
	plot([real(αpss(Δ))], [imag(αpss(Δ))], color=colors[2], marker="o", label=L"\alpha_+")
	plot([real(αmss(Δ))], [imag(αmss(Δ))], color=colors[4], marker="o", label=L"\alpha_-")
	
	tight_layout()
	
	ax3 = gca()
	# ax.set_xticks(range(-2,2,step=1))
	ax3.set_xlim([-0.25, 0.25]) 
	ax3.set_ylim([-0.25, 0.25]) 
	ax3.grid()
	
    xlabel(string("Re", L"\alpha_\pm"))
    ylabel(string("Im", L"\alpha_\pm"))
    title("")
    legend()
    gcf()
end

# ╔═╡ 031be99a-eb8b-420f-a155-af7c2a25b8f0
let
	close("all")
	
	# Plot phase vs. detuning Δ ---------------------------------------------------
	
	subplot(2, 2, 1)
	
	αp_list = αpss.(Δs)
	αm_list = αmss.(Δs)
	
	p = plot(Δs, angle.(αp_list), color=colors[2], label=L"\alpha_+")
	plot(Δs, angle.(αm_list), color=colors[4], label=L"\alpha_-")
	plot([Δ,Δ], [-1.5,1.5], linestyle="dashed", color="red")
	
	ax1 = gca()
	# ax1.set_xticks(range(-8,8,step=1))
	ax1.grid()
	
    ax1.set_xlabel(string(L"$\Delta$", " (MHz)"))
    title(string("Phase: ArcTan( Re ", L"\alpha_\pm", "/ Im ", L"\alpha_\pm", ")"))
    legend()
    gcf()
	
	
	
	# Plot amplitude vs. detuning Δ ------------------------------------------------
	
	subplot(2, 2, 2)
	plot(Δs, abs.(αp_list ), color=colors[2], label=L"\alpha_+")
	plot(Δs, abs.(αm_list ), color=colors[4], label=L"\alpha_-")
	plot([Δ,Δ], [0,0.22], linestyle="dashed", color="red")
	
	ax2 = gca()
	# ax2.set_xticks(range(-4,4,step=1))
	ax2.grid()
	
	ax2.set_xlabel(string(L"$\Delta$", " (MHz)"))
	
    # xlabel(string(L"$\Delta$", " (MHz)"))
    # ylabel(string("Amplitude ", L"|\alpha_\pm|"), loc="top")
    title(string("Amplitude ", L"|\alpha_\pm|"))
    legend()
    gcf()
	
	# Plot α in complex plane ---------------------------------------------------
	
	subplot(2, 2, 3) 
	plot(real.(αp_list), imag.(αp_list), color="gray")
	plot([real(αpss(Δ))], [imag(αpss(Δ))], color=colors[2], marker="o", label=L"\alpha_+")
	plot([real(αmss(Δ))], [imag(αmss(Δ))], color=colors[4], marker="o", label=L"\alpha_-")
	
	tight_layout()
	
	ax3 = gca()
	# ax.set_xticks(range(-2,2,step=1))
	ax3.set_xlim([-0.25, 0.25]) 
	ax3.set_ylim([-0.25, 0.25]) 
	ax3.grid()
	
    xlabel(string("Re", L"\alpha_\pm"))
    ylabel(string("Im", L"\alpha_\pm"))
    title("")
    legend()
    gcf()
	

	
end

# ╔═╡ 38d58263-6411-4189-9aca-95e99e7e558a
let
	close("all")
	
	sim = RQ
	
	tt, xx, yy, zz, ρρ, r = (sim.t, sim.x, sim.y, sim.z, sim.p, sim.r[1])
	
	# Plot Bloch components vs. time -----------------------------------------
	
	subplot(2,2,1)
	
    p = plot(tt, xx, color=colors[2], label=L"$x$", linewidth=0.8)
    plot(tt, yy, color=colors[4],label=L"$y$", linewidth=0.8)
    ax1 = gca()
    ax1.set_ylim([-1.1,1.1]) 
    plot(tt, zz, color=colors[6], label=L"$z$", linewidth=0.8)
	plot(tt, ρρ, color=colors[8], label=L"Tr $\rho^2$", linewidth=0.8)
    xlabel("t (μs)")
    ylabel("Bloch coordinates")
    title("Bayesian trajectory")
    ax1.legend(loc="lower right")
    gcf()
	
	
	# Plot measurement angle -----------------------------------------------------
	
	subplot(2,2,2)
	
	phase(ϕ) = exp(im * ϕ)
	angles = range(0, 2π, step=2π/100)
	phases = phase.(angles)
	
	
	p = plot([0, real(phase(ϕ))], [0, imag(phase(ϕ))], color="gray", linestyle="dashed")
	plot(real(phases), imag(phases), color="black", linewidth=0.8)

	
	ax2 = gca()
	ax2.set_xticks(range(-1,1,step=0.5))
	ax2.set_yticks(range(-1,1,step=0.5))
	ax2.grid()
	
	xlabel("Q")
	ylabel("I")
    title("Measurement angle")
	
	tight_layout()
    gcf()
	

    

	# Plot stochastic record  ----------------------------------------------------
	
	subplot(2, 2, 4)
    
   	p = plot(tt, r, color="gray", linewidth=0.1, alpha=0.9)

	ax1 = gca()
	ax1.set_xlim([first(tt), last(tt)]) 
    ax1.set_ylim([minimum(r) - 2, maximum(r) + 2]) 
	
    xlabel("t (μs)")
    ylabel("arbitrary units")
    title("stochastic record")
	
	tight_layout()
    gcf()
	
	
end

# ╔═╡ 83ec507b-a0ea-48b9-bb9e-0e2c70e385a6
let
	close("all")
	
	# Plot Bloch components vs. time -----------------------------------------
	
	subplot(2,2,1)
	
	sim = RQ2
	
	ttd, xxd, yyd, zzd, ρρd = sim.t, sim.x, sim.y, sim.z, sim.p
	rd1, rd2 = sim.r
	
    p = plot(ttd, xxd, color=colors[2], label=L"$x$", linewidth=0.8)
    plot(ttd, yyd, color=colors[4],label=L"$y$", linewidth=0.8)
    ax1 = gca()
    ax1.set_ylim([-1.1,1.1]) 
    plot(ttd, zzd, color=colors[6], label=L"$z$", linewidth=0.8)
	plot(ttd, ρρd, color=colors[8], label=L"Tr $\rho^2$", linewidth=0.8)
    xlabel("t (μs)")
    ylabel("Bloch coordinates")
    title("Bayesian trajectory")
    ax1.legend(loc="lower right")
    gcf()
	
	
	# Plot measurement angle -----------------------------------------------------
	
	subplot(2,2,2)
	
	phase(ϕ) = exp(im * ϕ)
	angles = range(0, 2π, step=2π/100)
	phases = phase.(angles)
	
	
	p = plot([0, real(phase(ϕd))], [0, imag(phase(ϕd))], color="gray", linestyle="dashed")
	plot([0, real(phase(ϕd + π/2))], [0, imag(phase(ϕd + π/2))], color="purple", linestyle="dashed")
	plot(real(phases), imag(phases), color="black", linewidth=0.8)

	
	ax2 = gca()
	ax2.set_xticks(range(-1,1,step=0.5))
	ax2.set_yticks(range(-1,1,step=0.5))
	ax2.grid()
	
	xlabel("Q")
	ylabel("I")
    title("Measurement angle")
	
	tight_layout()
    gcf()
	

	# Plot stochastic record  ----------------------------------------------------
	
	subplot(2, 2, 4)
    
   	p = plot(ttd, rd1, color="gray", linewidth=0.1, alpha=0.9)
	plot(ttd, rd2, color="purple", linewidth=0.1, alpha=0.4)
	
	plot([last(ttd)], [last(rd1)], color="gray", label=L"$\varphi$")
	plot([last(ttd)], [last(rd2)], color="purple", label=L"$\varphi + \pi/2$")

	ax1 = gca()
	ax1.set_xlim([first(ttd), last(ttd)]) 
    ax1.set_ylim([minimum(rd1) - 2, maximum(rd1) + 2]) 
	
    xlabel("t (μs)")
    ylabel("arbitrary units")
    title("stochastic record")
	legend(loc="lower right")
	
	tight_layout()
    gcf()
	
	
end

# ╔═╡ 761b6c4f-17f7-4735-bcac-a980c559d838
let
	close("all")
	
	sim = QR
	
	tt2, xx2, yy2, zz2, ρρ2 = (sim.t, sim.x, sim.y, sim.z, sim.pr)
	nn2, αp2, αm2 = (sim.n, sim.ap, sim.am)
	rs = sim.r[1]
	
	# Plot Bloch components vs. time ------------------------------------------
	
	subplot(2, 2, 1)
    
   	p = plot(tt2[1:i], xx2[1:i], color=colors[2], label=L"$x$", linewidth=0.8)
	plot([tt2[i]], [xx2[i]], color=colors[2], marker="o")
	
    plot(tt2[1:i], yy2[1:i], color=colors[4],label=L"$y$",  linewidth=0.8)
	plot([tt2[i]], [yy2[i]], color=colors[4], marker="o")
	
    plot(tt2[1:i], zz2[1:i], color=colors[6], label=L"$z$", linewidth=0.8)
	plot([tt2[i]], [zz2[i]], color=colors[6], marker="o")
	
	plot(tt2[1:i], ρρ2[1:i], color=colors[8], label=L"Tr $\rho_q^2$", linewidth=0.8)
	plot([tt2[i]], [ρρ2[i]], color=colors[8], marker="o")
	
	ax1 = gca()
	ax1.set_xlim([first(tt2), last(tt2)]) 
    ax1.set_ylim([-1.1,1.1]) 
	
    xlabel("t (μs)")
    ylabel("Bloch coordinates")
    title("Bayesian trajectory")
    ax1.legend(loc="lower right")
    gcf()
	
	
	# Plot photon number -------------------------------------------------------
	
	subplot(2, 2, 2)
	
	p1 = plot(tt2[1:i], nn2[1:i], color="purple", linewidth=0.8)
	plot([tt2[i]], [nn2[i]], color="purple", marker="o")
	
	ax2 = gca()
	# ax2.set_yticks(range(0.0, 0.25,step=0.05))
	ax2.grid()
	
	xlabel("t (μs)")
	ylabel(L"$\langle a^\dagger a \rangle$")
    title("Photon number")
	
	ax2.set_xlim([first(tt2), last(tt2)]) 
    # ax2.set_ylim([0, 0.3]) 
	
	
	tight_layout()
    gcf()
	
	
	# Plot αp, αm -------------------------------------------------------------
	
	subplot(2, 2, 3)
	
	αps = αp2[1:i]
	αms = αm2[1:i]
	
	plot(real.(αps), imag.(αps), color=colors[2], alpha=0.9, linewidth=0.3)
	plot([real(last(αps))], [imag(last(αps))], color=colors[2], marker="o", label=L"\alpha_+")
	plot(real.(αms), imag.(αms), color=colors[4], linestyle="dashed", alpha=0.9, linewidth=0.3)
	plot([real(last(αms))], [imag(last(αms))], color=colors[4], marker="o", label=L"\alpha_-")
	# plot([real(last(αps)), real(last(αms))], [imag(last(αps)), imag(last(αms))], color="black", label=L"\Delta \alpha", linewidth=1, alpha=0.9)
	# plot([0, 0.25 * real(exp(im * φ))], [0, 0.25 * imag(exp(im * φ))], color="black", linestyle="dashed", label=L"\varphi", linewidth=1, alpha=0.9)
	
	tight_layout()
	
	ax4 = gca()
	ax4.set_xticks(range(-0.3,0.3,step=0.1))
	ax4.set_yticks(range(-0.6, 0.6, step=0.3))
	ax4.grid()
	
    xlabel(string("Re", L"\alpha_\pm"))
    ylabel(string("Im", L"\alpha_\pm"))
    title("")
    legend()
	title("cavity states")
	ax4.legend(loc="upper right")
    gcf()
	
	
		
	# Plot cavity drive -----------------------------------------------------
	
# 	subplot(2, 2, 4)
	
# 	p = plot(tt2[1:i], (t -> ε).(tt2)[1:i], color=colors[8], label="")
# 	plot([tt2[i]], [(t -> ε).(tt2[i])], color=colors[8], marker="o")

# 	ax3 = gca()
# 	ax3.set_yticks(range(0,5,step=1))
# 	ax3.grid()
	
# 	ax3.set_xlim([first(tt2), last(tt2)]) 
#     ax3.set_ylim([0,5]) 
	
	
#     xlabel(string(L"$t$", " (μs)"))
# 	ylabel(string(L"$\epsilon(t)$", " (MHz)"))
#     gcf()
	
	# Plot stochastic record  ------------------------------------------
	
	subplot(2, 2, 4)
    
   	p = plot(tt2[1:i], rs[1:i], color="gray", linewidth=0.1, alpha=0.9)
	plot([tt2[i]], [rs[i]], color="gray", marker="o")
	
	tight_layout()
	
	ax1 = gca()
	ax1.set_xlim([first(tt2), last(tt2)]) 
    ax1.set_ylim([minimum(rs) - 2, maximum(rs) + 2]) 
	
    xlabel("t (μs)")
    ylabel("arbitrary units")
    title("stochastic record")
    gcf()
	
	


	
end

# ╔═╡ 0aff657d-99b8-49dc-93f3-e3a3fe058e7c
let
	close("all")
	
	
	# Plot Bloch components vs. time ------------------------------------------
	
	subplot(2, 2, 1)
    
   	p = plot(tt3[1:j], xx3[1:j], color=colors[2], label=L"$x$", linewidth=0.8)
	plot([tt3[j]], [xx3[j]], color=colors[2], marker="o")
	
    plot(tt3[1:j], yy3[1:j], color=colors[4],label=L"$y$",  linewidth=0.8)
	plot([tt3[j]], [yy3[j]], color=colors[4], marker="o")
	
    plot(tt3[1:j], zz3[1:j], color=colors[6], label=L"$z$", linewidth=0.8)
	plot([tt3[j]], [zz3[j]], color=colors[6], marker="o")
	
	plot(tt3[1:j], ρρ3[1:j], color=colors[8], label=L"Tr $\rho_q^2$", linewidth=0.8)
	plot([tt3[j]], [ρρ3[j]], color=colors[8], marker="o")
	
	ax1 = gca()
	ax1.set_xlim([first(tt3), last(tt3)]) 
    ax1.set_ylim([-1.1,1.1]) 
	
    xlabel("t (μs)")
    ylabel("Bloch coordinates")
    title("Bayesian trajectory")
    ax1.legend(loc="lower right")
    gcf()
	
	
	# Plot photon number -------------------------------------------------------
	
	subplot(2, 2, 2)
	
	p1 = plot(tt3[1:j], nn3[1:j], color="purple", linewidth=0.8)
	plot([tt3[j]], [nn3[j]], color="purple", marker="o")
	
	ax2 = gca()
	ax2.set_yticks(range(0.0, 0.25,step=0.05))
	ax2.grid()
	
	xlabel("t (μs)")
	ylabel(L"$\langle a^\dagger a \rangle$")
    title("Photon number")
	
	ax2.set_xlim([first(tt3), last(tt3)]) 
    ax2.set_ylim([0, 0.3]) 
	
	
	tight_layout()
    gcf()
	
	
	# Plot αp, αm -------------------------------------------------------------
	
	subplot(2, 2, 3)
	
	αps = αp3[1:j]
	αms = αm3[1:j]
	
	plot(real.(αps), imag.(αps), color=colors[2], alpha=0.9, linewidth=0.3)
	plot([real(last(αps))], [imag(last(αps))], color=colors[2], marker="o", label=L"\alpha_+")
	plot(real.(αms), imag.(αms), color=colors[4], linestyle="dashed", alpha=0.9, linewidth=0.3)
	plot([real(last(αms))], [imag(last(αms))], color=colors[4], marker="o", label=L"\alpha_-")
	# plot([real(last(αps)), real(last(αms))], [imag(last(αps)), imag(last(αms))], color="black", label=L"\Delta \alpha", linewidth=1, alpha=0.9)
	# plot([0, 0.25 * real(exp(im * φ))], [0, 0.25 * imag(exp(im * φ))], color="black", linestyle="dashed", label=L"\varphi", linewidth=1, alpha=0.9)
	
	tight_layout()
	
	ax4 = gca()
	ax4.set_xticks(range(-0.3,0.3,step=0.1))
	ax4.set_yticks(range(-0.6, 0.6, step=0.3))
	ax4.grid()
	
    xlabel(string("Re", L"\alpha_\pm"))
    ylabel(string("Im", L"\alpha_\pm"))
    title("")
    legend()
	title("cavity states")
	ax4.legend(loc="upper right")
    gcf()
	
	
	# Plot stochastic record  ------------------------------------------
	
	subplot(2, 2, 4)
    
   	p = plot(tt3[1:j], rs1[1:j], color="gray", linewidth=0.1, alpha=0.9)
	plot([tt3[j]], [rs1[j]], color="gray", marker="o", label=L"$\varphi$")
	plot(tt3[1:j], rs2[1:j], color="purple", linewidth=0.1, alpha=0.4)
	plot([tt3[j]], [rs2[j]], color="purple", marker="o", label=L"$\varphi + \pi/2$")
	
	tight_layout()
	
	ax1 = gca()
	ax1.set_xlim([first(tt3), last(tt3)]) 
    ax1.set_ylim([minimum(vcat(rs1,rs2)) - 2, maximum(vcat(rs1,rs2))+ 2]) 
	
    xlabel("t (μs)")
    ylabel("arbitrary units")
    title("stochastic record")
	legend()
    gcf()
	
	


	
end

# ╔═╡ 31371407-5798-49a4-a83a-2c068d622c4b
# Plotting
function plot_solution(sol; plot_title="Rabi Oscillation", xlab=L"$t$")
	
	close("all")
    
    tt = sol[1]
    ρt = sol[2]
    
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
    xlabel(xlab)
    ylabel("Bloch coordinates")
    title(plot_title)
    legend()
    gcf()
end

# ╔═╡ 99e90182-1ee7-45ef-8656-cc5a1e01d564
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

# ╔═╡ 53f32791-a721-4d80-96ce-fae9f533dbaa
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

# ╔═╡ 7d4f4a02-f3b8-4369-aad3-03ff4bc7ec0b
function plot_ensemble(sol_ens; α=0.1, linewidth=1, labels=false, average=false)
    close("all")
	tt1 = sol_ens[1]
    evs = collect(map(ρs -> expects.(ρs), sol_ens[2]));

    for i in 1:50
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

# ╔═╡ a255f7a5-d617-4444-80a4-4e4e7762b236
green(text; title="Note") = Markdown.MD(Markdown.Admonition("correct", title, [text]))

# ╔═╡ 6e86b5ee-5b38-4245-b555-c208151c0c53
red(text; title="Note") = Markdown.MD(Markdown.Admonition("danger", title, [text]))

# ╔═╡ fb32f241-34bb-413b-99d2-c0155b363670
tan(text; title="Note") = Markdown.MD(Markdown.Admonition("warning", title, [text]))

# ╔═╡ 0d715d29-173d-4e4f-9ea2-74e0d9441b50
tan(md"""Play with the plot below to see how the mixer processes signals.""")

# ╔═╡ e730e11b-6fa5-4121-80eb-69a59d74b852
blue(text; title="Note") = Markdown.MD(Markdown.Admonition("note", title, [text]))

# ╔═╡ 30b5b2c5-1243-40c2-aa23-ec1a4e1b50da
if Δ == 0.
	green(md"""When the cavity is driven on its bare resonance, coherent states are maximally distinguishable in phase.""",title="Driving bare resonance")
	
elseif Δ == 3
	blue(md"""A detuning of Δ = 3 MHz ~ -χ corresponds to driving on resonance with the α+ coherent state.""",title="Driving α+ state.")
	
elseif Δ == -3
	blue(md"""A detuning of Δ = -3 MHz ~ χ  corresponds to driving on resonance with the α- coherent state.""", title="Driving α- state.")
	
elseif abs(Δ) > abs(3χ)
	red(md"""When the detuning is large, the coherent states are not easily distinguished in either amplitude or phase.""", title="Driving off resonance.")

else
	tan(md"""Change the value of Δ (resonator-drive detuning) using the slider to see how it affects the coherent states.""")

		
end

# ╔═╡ df85f145-8739-44aa-9f1f-2b3ec699a4a4
if ϕ == 0 || ϕ == π
	blue(md"""When ϕ = 0 or π, information is collected about the qubit state. Backaction is purely informational (towards $\ket{\pm z}$) and the state remains in the x-z plane.""",title="Informational backaction")
	
elseif ϕ == π/2 || ϕ == 3π/2
	green(md"""When ϕ = π/2 or 3π/2, information is collected about the qubit phase. Backaction is purely phasal (rotations in x-z plane). """,title="Phase backaction")
	
else
	tan(md"""For intermediate values of ϕ, there is both informational ($\ket{\pm z}$) and phase ($\sigma_z$ rotations) backaction. """, title="Informational and phase backaction")
	

		
end

# ╔═╡ c8aa6f9a-a2b7-45d7-9e39-344380a3909b
hint(text; title="Hint") = Markdown.MD(Markdown.Admonition("hint", title, [text]))

# ╔═╡ 6bb8877f-d861-492e-b816-2396fcd6d59c
md"""
!!! warning "Edit"

	We need to use simple master equation API rather than `bayesian` here. I'm not sure this currently exists in `QuantumCircuits.jl`, but should probably be the `jump-no-jump` method. Or does `bayesian` automatically do that in the absence of measurement?
"""

# ╔═╡ Cell order:
# ╠═f230c876-f933-4ee4-a079-a66147dea73c
# ╠═4c45fe3e-cd69-11eb-20cd-7bfb98c040cf
# ╠═377a3336-20bd-4baa-a033-af8bbc8668a8
# ╠═3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─f24df928-22c4-4943-bc9c-218750fe4da7
# ╟─3d17c4b0-e051-41b0-9f20-bb22c4db4398
# ╟─3500ae58-7577-4af7-a458-5613310471d9
# ╟─4e2680a7-8035-4161-b4e3-abb3f54da99d
# ╟─d330e502-e4d9-482d-949e-d9c3d0166b01
# ╟─8c3e20ca-4018-4448-8848-073a2bbcf3bb
# ╟─54690fb9-7f41-4166-b231-d0915ef293d7
# ╟─2a548d96-451f-4bd7-b9e5-7c3fa85c2abb
# ╟─f1d33aff-03d3-4223-b67d-d7e9be26c168
# ╟─3b125df9-c86d-4e4a-949e-8701b91d58e4
# ╠═019c8da1-8cc8-452b-839f-bc50772f9d79
# ╟─6ae79dd2-dc81-4d00-8bce-394d788a9207
# ╟─ba89c6ca-abd3-48ce-b1da-54f282fbd94d
# ╟─0d715d29-173d-4e4f-9ea2-74e0d9441b50
# ╟─1abcbbed-e7e6-4a06-a1eb-3a0922785323
# ╟─970d9217-e5e9-4adb-8ada-5d9a225523c3
# ╟─e0eda8fc-2ff6-49b0-82b3-156954efedbe
# ╟─82eed712-4adc-4408-baaa-64f49373c405
# ╟─4652677b-8602-4640-836b-07fab787e88c
# ╟─fc688b0e-7de9-473f-9c15-1e8492944f28
# ╟─38b793c7-75c1-4b26-a6de-cd213db0dd2a
# ╟─e92eca0d-854f-4257-9e54-ed7acf79ecff
# ╟─367819f8-b540-445f-83a4-a9c7982e878c
# ╟─2f7a96c1-d043-407f-84f8-68efea7ded0e
# ╟─c2f7d614-0d49-4283-b574-c0713e41824a
# ╟─b08e7d51-3e2c-431a-8471-12914bcd1384
# ╟─3bc48dd3-363b-4dc8-8f9d-b3edc7700bce
# ╠═e4f0d549-b9d2-48ef-bbe7-cd4731156c77
# ╟─58d845a6-45e7-450e-b6c2-009a605e1c7c
# ╠═5702d0d1-6e34-469e-926b-c369dbde78b3
# ╟─d44ea114-6e92-4790-978d-7286394e9506
# ╟─8a8987ef-225c-4538-92d0-ca2004266053
# ╟─82ffdfde-b2ff-4859-8ca6-119a4d5d95cd
# ╠═bddc09b0-f90a-40e5-be77-d1afa4309fa4
# ╟─30b5b2c5-1243-40c2-aa23-ec1a4e1b50da
# ╟─d6f73306-ab64-4fba-9955-55298ab566e5
# ╟─aec1af64-50c5-4560-8461-e45bfc2f406d
# ╠═0a6b57b5-d5bf-4598-a01a-eced9358e452
# ╟─031be99a-eb8b-420f-a155-af7c2a25b8f0
# ╟─b4cbaf81-c33c-42a9-ba79-437c32f47f70
# ╟─de3fe955-780c-489f-9b6d-f3a8c3bf140d
# ╟─cd661d56-ebd0-4334-9ae6-47e738d11b54
# ╟─88e1d88a-3565-42b3-9732-02862163a12b
# ╠═4284173a-be05-4b58-a8d9-7189301344fd
# ╠═139565f9-0219-4939-a470-4252397dd6a0
# ╟─ca554621-0fe8-4d7a-bd3f-acdf795648c4
# ╟─fe1986ff-00c7-4b9f-8d20-9a38c85e2a28
# ╠═38d58263-6411-4189-9aca-95e99e7e558a
# ╟─073c34d3-a243-4560-9862-28f15b24e685
# ╟─df85f145-8739-44aa-9f1f-2b3ec699a4a4
# ╠═b6ae1b65-aa9d-4d7e-b07e-f51ebdcd44b3
# ╟─de277991-5f34-4906-a9bb-b1832131b45f
# ╟─eb1ac4e6-5bc3-4c67-85b0-e4be923362a4
# ╠═83ec507b-a0ea-48b9-bb9e-0e2c70e385a6
# ╟─ede93004-ef9e-461f-b4ea-062f9d7879f3
# ╟─9d0fe6ee-aa00-43c8-b40f-937dd6024d74
# ╠═05b536fa-e464-4476-8c27-65ffe89d09a6
# ╠═55b946e9-20e2-4ea8-85e2-455dd9b6b8a1
# ╟─865752f3-c563-4a6c-be2f-3dbe73d38406
# ╟─9baf661a-feb9-463d-a75e-69848ad888ae
# ╠═761b6c4f-17f7-4735-bcac-a980c559d838
# ╟─28c3789a-303e-4c08-83d2-95066828a1d7
# ╟─aaf95d04-aeb2-47b4-99c6-9499fe27a33f
# ╟─fcc6db99-6ff7-41d5-b2af-ff07899a0250
# ╠═0aff657d-99b8-49dc-93f3-e3a3fe058e7c
# ╟─c1e0ad96-5335-4f11-abe0-5b4f5274c36e
# ╟─7e62e6c6-aa3c-4350-901c-15d017b8db42
# ╠═989ba897-50b4-4d82-be56-6dcb439e8aca
# ╠═3d1472b2-30ac-4237-ace3-00ad5f5768b6
# ╠═76d5307b-f58b-4ab1-8ba5-b43491c54552
# ╠═a823cf8d-8637-44b9-af2a-8feefd22d986
# ╠═3604c61d-6cf4-4d6a-94ea-678c196d336e
# ╠═52e2feca-5070-444d-aef6-c85c6404d391
# ╠═b136e6e0-6e88-4008-94bb-7756716e913e
# ╟─dfb8774e-bb14-4f8b-9492-ecd2220e393d
# ╟─71acc211-c8b8-4059-8002-be14fff7d822
# ╠═4de027ac-2080-499f-9a99-7ef3007e1384
# ╟─b2235afc-1bc9-4cf4-ae21-6f535eb37e96
# ╟─31371407-5798-49a4-a83a-2c068d622c4b
# ╟─99e90182-1ee7-45ef-8656-cc5a1e01d564
# ╟─53f32791-a721-4d80-96ce-fae9f533dbaa
# ╟─7d4f4a02-f3b8-4369-aad3-03ff4bc7ec0b
# ╟─0a2219d8-2fbd-4b03-867f-966c4a595e78
# ╟─a255f7a5-d617-4444-80a4-4e4e7762b236
# ╟─6e86b5ee-5b38-4245-b555-c208151c0c53
# ╟─fb32f241-34bb-413b-99d2-c0155b363670
# ╠═e730e11b-6fa5-4121-80eb-69a59d74b852
# ╟─c8aa6f9a-a2b7-45d7-9e39-344380a3909b
# ╟─6bb8877f-d861-492e-b816-2396fcd6d59c
