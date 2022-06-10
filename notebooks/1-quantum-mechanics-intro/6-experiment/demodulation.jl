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

# ╔═╡ 8ea33322-60b8-47b1-84e6-9069cc9c4f46
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

	include("utilities/single-qubit-operators.jl")
	include("utilities/plotting.jl")

	include("notebooks/table-of-contents.jl")
	include("resources.jl")
	
	md" # Packages and julia files"
end

# ╔═╡ 79f8fed1-b2c1-4276-9c1e-1bef826e3089
mdp("Readout is a process with different steps at room temperature, in the cavity, and at the qubit. In this section, we'll walk backward through the experimental setup. We'll start from ", demodulation📔, " of the readout signal at room temperature, look back at what the readout says about the ", resonator_states📔, ", then finally predict the ", measurement_backaction📔, " on the qubit")

# ╔═╡ 04f8e225-5cec-44c1-8711-d89bdb7da9f2
md"""
In this interactive notebook, we focus on the first part -- the role of demodulation in qubit readout.
"""

# ╔═╡ 00ae069c-3d04-4937-9e12-51637a237aaf
mdp(table_of_contents📔)

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents(title="Demodulation")

# ╔═╡ 3500ae58-7577-4af7-a458-5613310471d9
md"""
# Readout at room temperature
"""

# ╔═╡ 4e2680a7-8035-4161-b4e3-abb3f54da99d
md" ## Background"

# ╔═╡ d330e502-e4d9-482d-949e-d9c3d0166b01
md"""
In superconducting qubit readout, the readout signal "RO" reflects from the cavity at cryogenic temperatures and travels along the transmission line. After amplification, it is demodulated at room temperature using an IQ mixer as depicted in the figure below.

"""

# ╔═╡ 354d762d-c3b3-4f9f-b42a-7e9e0d78086e
md"""
$(LocalResource(string(path, "/notebooks/figures/demodulation-single-quadrature.png")))
"""

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
md" ## Mixer simulation"

# ╔═╡ 019c8da1-8cc8-452b-839f-bc50772f9d79
begin
	χ = 2π * (-0.47)
	κ = 2π * (1.56)
	ε0 = 1 # MHz
	
	αpss(Δ) = (2ε0/κ) / (1 + im * (2(Δ + χ)/κ))
	αmss(Δ) = (2ε0/κ) / (1 + im * (2(Δ - χ)/κ))
	md" 🌀 cavity parameters"
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

Integrate the mixed signals (LO x readout) across $n_c$ bins to low-pass filter the result: 

 $n_c$ = 0
$(@bind nc Slider(1:25, default=5)) 
25 

Decrease the detuning `ΔRO` between LO and RO to switch from heterodyne to homodyne measurement (discussed in next section):

 $\Delta_{RO}$ = 0 MHz
$(@bind ΔRO Slider(0:100, default=40)) 40 MHz 

"""

# ╔═╡ 6ae79dd2-dc81-4d00-8bce-394d788a9207
begin
	I = real(αpss(0))
	Q = imag(αpss(0))
	ωRO = 6.67e3 #  MHz
	ωLO = ωRO + ΔRO # (6.60e3) # MHz
	B = 1
	ts = range(0, 0.064, step=1/ωRO / 50)
	φ0 = 0
	
	md"🌀 define parameters"
end

# ╔═╡ 3b125df9-c86d-4e4a-949e-8701b91d58e4
begin
	s(t) = I * cos((2π * ωRO) * t) + Q * sin((2π * ωRO )* t)
	lo(t) = B * cos((2π * ωLO) * t + φ0)
	prod(t) = B * (I * cos(-(2π * ΔRO) * t + φ0) + Q * sin(-(2π * ΔRO) * t + φ0))
	
	md" 🌀 define signals"
end

# ╔═╡ ba89c6ca-abd3-48ce-b1da-54f282fbd94d
begin
	ωs = range(ωRO - 150, ωRO + 150, step=1)
	
	ωROs = zeros(length(ωs))
	ωROs[argmin(abs.(ωs .- ωRO))] = sqrt(I^2 + Q^2)
	
	ωLOs = zeros(length(ωs))
	ωLOs[argmin(abs.(ωs .- ωLO))] = B
	
	md" 🌀 define frequency domain series"
end

# ╔═╡ 82eed712-4adc-4408-baaa-64f49373c405
md"""
Change the time-window displayed: $t_f$ = $(@bind tfstr Select(["0.5", "1", "2", "4", "8", "16", "32", "64"], default="16")) ns
"""

# ╔═╡ 4652677b-8602-4640-836b-07fab787e88c
let 
	tf = parse(Float64, tfstr)
	newsig = coarse_grain(s.(ts) .* lo.(ts); n=nc)
	
	# Signal time series ---------------------------------------------------------
	
	ptime = plot(xlim=[0, tf], ylim=[-0.22, 0.22], xlabel = "t (ns)", ylabel = "arbitrary units", title="time domain", legend=:bottomright)
	
	tss = ts * 1e3
	# tff = tf * 1e3
   	if show_RO 
		plot!(tss, s.(ts), color="gray", linewidth=1, alpha=0.9, label="readout")
	end
	
	if show_LO
		plot!(tss, lo.(ts) ./ 5, color="purple", linewidth=1, alpha=0.9, label="LO / 5")
	end
	
	if show_LORO
		plot!(tss, s.(ts) .* lo.(ts), color="green", linewidth=1, alpha=0.9, label="LO x readout") end
	
	if show_int
		plot!(tss, newsig, color="blue", linewidth=1, alpha=0.9, label="integrated")
	end
	
	if show_env
		plot!(tss, 0.5* prod.(ts), color="black", linewidth=1, alpha=0.9, label="envelope", linestyle=:dash) end

	# Frequency domain ---------------------------------------------------------

	pfreq = plot(xlabel = "ω (GHz)", ylabel = "arbitrary units", title = "frequency domain", legend=:bottomright)

   	plot!(ωs ./1000, ωROs, color="gray", label="readout") 
	plot!(ωs ./1000, ωLOs, color="purple", label="LO / 5") 

	l = @layout [time{0.5h}; freq{0.5h}]

	plot(ptime, pfreq, layout = l)
	
end

# ╔═╡ fc688b0e-7de9-473f-9c15-1e8492944f28
md"""
## Types of measurement
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

See $Krantz_et_al Sec. V for more details on experimental implementation of homodyne and heterodyne measurement.


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

In some sources (e.g. $Campagne_Ibarcq_et_al, $Gambetta_et_al), "homodyne" and "heterodyne" distinguish the number of quadratures detected for a single qubit, with "homodyne" indicating single-quadrature measurement and "heterodyne" indicating dual-quadrature measurement.

However, in the signal processing literature, as well as in $Krantz_et_al, "homodyne" and "heterodyne" distinguish the detuning of the LO from the signal frequency, as discussed above, with "homodyne" ("heterodyne") referring to (non-)zero detuning. Using this terminology, it is possible to make single- or dual-quadrature measurements using *either* homodyne *or* heterodyne detection. The difference lies only in whether the signal is brought down to DC analogically (by the mixer -- homodyne measurement) or digitally (by the computer -- heterodyne measurement).

This latter terminology will be used in these tutorials.
"""

# ╔═╡ 7e62e6c6-aa3c-4350-901c-15d017b8db42
md" # Utilities "

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
	
	function traj(sol::Solution)
	  t, ρ, r = (sol.t, sol.ρ, sol.r)
	  x, y, z = map(op -> expectations(sol, op), [σx, σy, σz])
	  p = (typeof(ρ[1]) <: Ket) ? [1.0 for el in ρ] : real(expect.(ρ, ρ))
	  traj(t, x, y, z, p, r)
	end
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

# ╔═╡ c8aa6f9a-a2b7-45d7-9e39-344380a3909b
hint(text; title="Hint") = Markdown.MD(Markdown.Admonition("hint", title, [text]))

# ╔═╡ Cell order:
# ╟─79f8fed1-b2c1-4276-9c1e-1bef826e3089
# ╟─04f8e225-5cec-44c1-8711-d89bdb7da9f2
# ╟─00ae069c-3d04-4937-9e12-51637a237aaf
# ╟─3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─3500ae58-7577-4af7-a458-5613310471d9
# ╟─4e2680a7-8035-4161-b4e3-abb3f54da99d
# ╟─d330e502-e4d9-482d-949e-d9c3d0166b01
# ╟─354d762d-c3b3-4f9f-b42a-7e9e0d78086e
# ╟─54690fb9-7f41-4166-b231-d0915ef293d7
# ╟─2a548d96-451f-4bd7-b9e5-7c3fa85c2abb
# ╟─f1d33aff-03d3-4223-b67d-d7e9be26c168
# ╠═3b125df9-c86d-4e4a-949e-8701b91d58e4
# ╟─019c8da1-8cc8-452b-839f-bc50772f9d79
# ╟─6ae79dd2-dc81-4d00-8bce-394d788a9207
# ╟─ba89c6ca-abd3-48ce-b1da-54f282fbd94d
# ╟─0d715d29-173d-4e4f-9ea2-74e0d9441b50
# ╟─1abcbbed-e7e6-4a06-a1eb-3a0922785323
# ╟─970d9217-e5e9-4adb-8ada-5d9a225523c3
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
# ╟─e4f0d549-b9d2-48ef-bbe7-cd4731156c77
# ╟─7e62e6c6-aa3c-4350-901c-15d017b8db42
# ╟─76d5307b-f58b-4ab1-8ba5-b43491c54552
# ╟─a255f7a5-d617-4444-80a4-4e4e7762b236
# ╟─6e86b5ee-5b38-4245-b555-c208151c0c53
# ╟─fb32f241-34bb-413b-99d2-c0155b363670
# ╟─e730e11b-6fa5-4121-80eb-69a59d74b852
# ╟─c8aa6f9a-a2b7-45d7-9e39-344380a3909b
# ╠═8ea33322-60b8-47b1-84e6-9069cc9c4f46
