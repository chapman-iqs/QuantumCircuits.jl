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

# ╔═╡ 3661abf7-84e6-4e39-854b-ad7741a8ff36
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
	using ProgressMeter
	using LsqFit

	using DataFrames
	using CSV
	
	include("utilities/two-qubit-operators.jl")
	include("utilities/plotting.jl")
	include("utilities/utilities.jl")

	
	include("notebooks/table-of-contents.jl")
	include("resources.jl")

	include("notebook-utilities.jl")
	
	md" # Packages and julia files"
end

# ╔═╡ bd133c05-a98e-4b71-86d3-24c4519e9038
using FFTW

# ╔═╡ a52844ea-9a74-42b0-b37a-7cc8cd79a0a6
md"""
⏱ Tuesday June 14, 2022
"""

# ╔═╡ 01c6cbbc-1945-40bb-be1a-31abefa38bba
md"""
✍🏼 Sacha Greenfield
"""

# ╔═╡ e625f413-e93b-48c9-b14d-181c37b0fe0a
md"""
☄️ In this notebook, I discuss how to simulate fluctuators using an exponential waiting time distribution.
"""

# ╔═╡ 249f33ea-92aa-4b9c-bdfa-d35fe5254225
mdp(table_of_contents📔)

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents(title="Simulating fluctuators")

# ╔═╡ 4e92db51-ef23-4a88-ba33-92771472ff2f
md"""
# Introduction
"""

# ╔═╡ 4c838a57-4cf0-4218-980e-102359996fb9
md"""
# Modeling time between rare events
"""

# ╔═╡ 3a8257a5-6997-4f37-b4bc-170c722bbcfd
md"""
## Poisson distribution
"""

# ╔═╡ 48a30585-f3ae-4cc4-ad58-fb5b07b26795
md"""
f = γ/2π =  100 kHz $(@bind f Slider([0.1,1,10,100], default=1)) 100 MHz
"""

# ╔═╡ 25b85405-97d2-4a6e-ac77-be9823b3193f
md"""
τ = 1 ns $(@bind τ Slider(10.0 .^(-3:1), default=1)) 10 μs
"""

# ╔═╡ 49c97425-9008-4109-8544-db805396cdc4
md"""
f = $f MHz, τ = $(round(τ, digits=4)) μs
"""

# ╔═╡ a3940d5a-5b76-4224-bc4b-1cae947505e6
let
	global γ = 2π * f
	n = γ * τ
	
	Pkγ(k) = (γ * τ)^k / factorial(big(k)) * exp(-γ * τ)
	lim = minimum([Int(round(3n)), 100])
	ks = 0:lim
	pl1 = plot(ks, Pkγ.(ks), ylabel = L"P_\gamma(k)", legend=:none, size=(500,400), markers=true)

	pdist = Poisson(n)
	x = rand(pdist, 100000)
	pl2 = histogram(x, normalize=true, xlabel=string("number of events in interval τ = ", round(τ, digits=4), " μs, for rate γ = (2π) * ", f, " MHz"), ylabel="pdf", legend=:none)

	plot(pl1, pl2, layout = (2,1), link=:x)
end

# ╔═╡ 727c5a36-1f56-43be-b75a-f7cccc910a4d
md"""
## Exponential distribution
"""

# ╔═╡ 21a93094-5aac-4097-b8af-e4e974f49dcf
md"""
While the Poisson distribution accurately models the underlying physical process, it is not an efficient simulation method because it requires resampling the fluctuator sign thousands of times during the simulation. The more efficient method is to recognize the connection between the Poisson distribution and the exponential distribution. In the limit of **rare events**, i.e. $n \ll 1$, the time interval $T$ *between events* goes as

$P_\gamma(T = t) = 2\gamma e^{-2\gamma t}$

where $\gamma$ is the same parameter as before.
"""

# ╔═╡ d135dad7-8434-4c56-88f8-4feeafe3eb11
let
	Ptγ(t) = 2γ * exp(-2γ * t)
	ts = (0:0.01:5) ./ γ
	pl1 = plot(ts, Ptγ.(ts), ylabel = L"P_\gamma(t)", legend=:none, size=(500,400))

	edist = Exponential(1/γ)
	y = rand(edist, 10000)
	pl2 = histogram(y, xlabel="time between events (μs)", ylabel="counts", legend=:none) #bins=1:80, xticks=0:10:80,)

	plot(pl1, pl2, layout = (2,1), link=:x)
end

# ╔═╡ 2be13c7c-1f23-4d70-87ae-f997c8ec050d
md"""
# Modeling the fluctuators
"""

# ╔═╡ b21cfda7-a6f9-4200-8676-b089c952f0c4
md"""
Based on the above discussion, the most efficient method to simulate the fluctuators is to sample the switching times of each fluctuator in advance of the simulation. This generates a unique time-dependent Hamiltonian for each trajectory.
"""

# ╔═╡ 8ccc5fa4-493a-4e17-9860-9a2fc1500857
md"""
## Single fluctuator
"""

# ╔═╡ 495e2820-4029-4c95-a736-82e258ebedb7
md"""
As previously indicated, the number $k$ of switches that the fluctuator experiences within a time interval $T$ is Poisson distributed. Thinking in terms of fluctuators, we can understand the exponential distribution between switching events as coming from the correlation function $C(t) \equiv \langle \chi(t) \chi(0) \rangle$, where $\chi(t) = \pm 1$ is the sign of the fluctuator at time $t$. The number of switches $k$ is exactly the sign $(-1)^k$ for the $k^{\text{th}}$ contribution to $C(t)$; so letting $P_\gamma(X = k)  \equiv P_{k,\gamma}$,

$C_\gamma(t) = \sum_{k=0}^\infty \chi(t) \chi(0) P_{k,\gamma} = \sum_{k=0}^\infty (-1)^k \frac{(\gamma t)^k}{k!}e^{-\gamma t} = e^{-2 \gamma t}.$

Renormalizing yields $P_\gamma(T = t) = 2\gamma e^{-2\gamma t} \equiv P_{t,\gamma}$.
"""

# ╔═╡ 4465eafa-85be-44ca-b5ee-3f62dcf83cf8
md"""
A single fluctuator is characterized by its magnitude $|v|$ and switching rate $\gamma$. The sign $\chi(t)$ switches every interval $T \sim P_t$, so that

$v(t) = |v| \chi(t).$


"""

# ╔═╡ 6d7d40c5-1b1c-48b3-a6a7-635a628ee488
let
	v = 1
	S(ω) = v^2 * (1/π) * (2γ) / ((2γ)^2 + ω^2)
	ωs = γ .* (0:0.01:10)
	plot(ωs, S.(ωs), xlabel = "frequency ω (MHz)", ylabel = "S(ω) (units?)", legend=:none, size=(500,400))
end

# ╔═╡ 0a0a5f1c-7a2f-4a04-bee5-f368c1b12900
md"""
### Simulation
"""

# ╔═╡ b01b386c-d76b-4383-851e-9316a824b994
md"""
To simulate the fluctuating sign, we use the following approach:
1. Sample an array of switching times $\{\tau_0, \tau_1, ..., \tau_n\}$ such that $\sum_k \tau_k > t_f$, the duration of the simulation. Each time $\tau_k$ is i.i.d. with the  [exponential distribution](#727c5a36-1f56-43be-b75a-f7cccc910a4d). Define intermediate times $T_j \equiv \sum_{k = 0}^j \tau_k$.

2. Pick the initial sign to be $\chi_0 = \pm 1$ with equal probability.

3. Define the time-dependent sign as 

$\chi(t) = \chi_0 \times \begin{cases} 1 & \text{for } t < T_0 \\ 
-1 & \text{for } T_0 \le t < T_1 \\
(-1)^2 & \text{for } T_1 \le t < T_2 \\
\vdots \\
(-1)^j & \text{for } T_{j-1} \le t < T_j \\
\vdots \\
(-1)^n & \text{for } T_{n-1} \le t < t_f \\
\end{cases}$
"""

# ╔═╡ d3aa8364-5f56-442f-8990-e4ccad33ead6
md"""
An example fluctuator trajectory is plotted for a fluctuator of frequency $\gamma = 2\pi \times$ $f MHz:
"""

# ╔═╡ 72b67c33-12b0-4d04-b2aa-b18fc6faeb4a
function switching_times(γ, tf)
	dist = Exponential(1/γ)
	T = 0.0
	Ts = []
	while T < tf
		T += rand(dist)
		push!(Ts, T)
	end
	return Ts
end

# ╔═╡ cf7c7b7b-8adb-49a7-b245-83d440d2b9c5
function fluctuator_sign(γ, tf)
	
	# get list of switching times
	Ts = switching_times(γ, tf)

	# initial sign
	χ0 = rand([-1,1])

	return t -> let j = (findfirst(t .< Ts) - 1)
					return χ0 * (-1)^j
	end

end

# ╔═╡ 70648025-b722-4920-a876-c22dfa8f1d9a
begin
	tf = 10.0
	dt = 1e-3
	χ = fluctuator_sign(γ, tf)
end

# ╔═╡ 91507a5a-dc57-4b15-9fd4-309d899210d0
begin
	ts = 0.0:dt:tf
	single_fluc = χ.(ts)
	plot(ts, single_fluc, xlabel = "t (μs)", ylabel="fluctuator sign", legend=:none, title=string("γ = 2π * ",f, " MHz"), titlefontsize=13, marker=true, linealpha=0.5)
end

# ╔═╡ acf8d55f-2676-4114-8163-c6cb49a8a8a4
md"""
## Multiple fluctuators
"""

# ╔═╡ 1e6c50f7-31dd-4e61-83a9-b7a56b5c5f13
md"""
### Choosing the frequencies
"""

# ╔═╡ fdd591eb-1934-4cc3-9a7f-afcef3031972
md"""
Note that here, I have chosen the high frequency cutoff to be $1/(20 dt)$ or $(1/20dt) MHz. We could go all the way up to 1 GHz, but this is difficult to simulate, less relevant to 1/f noise, and may be physically unrealistic. I have to check this.
"""

# ╔═╡ 8631b550-01aa-438a-b659-8818bddfa4cb
N = 100

# ╔═╡ d001c6b6-5d39-4af3-b410-18d9e0a2ad57
begin
	(fmin, fmax) = 1/tf, 1/(20dt)
	flog = range(log(fmin), log(fmax), length=N)
	fi = exp.(flog)
end

# ╔═╡ b1d4e9d2-8d9f-4729-8078-564391fcd817
	pl2 = histogram(fi, xlabel="γ (MHz)", ylabel="counts", title="fluctuator distribution", legend=:none, yaxis=(:log10)) #bins=1:80, xticks=0:10:80,)

# ╔═╡ 9e067f48-df19-41f2-b76a-6222dd22a9b6
md"""
### Choosing the magnitude
"""

# ╔═╡ 49b393f8-28a9-491d-aa6a-2e62bfa29de1
md"""
When we were considering one fluctuator, the magnitude $|v|$ was not very important. However, now that there are many fluctuators, the relative magnitudes of fluctuators at different frequencies could matter. For now, I'll assume equal magnitude at all frequencies and allow the distribution in frequency space to control the effect of the fluctuators. However, I will normalize so that the strength per fluctuator is $v_i \equiv |v| / N$, with $N$ the total number of fluctuators, and I'll set $|v| = 1$ as before.
"""

# ╔═╡ 1a1e2154-b0ed-4d21-8456-92ac61c735d7
begin
	χs = []
	vi = 1/N
	for f in fi
		γ = 2π * f
		χ = fluctuator_sign(γ, tf)
		push!(χs, χ)
	end
end

# ╔═╡ da396e7c-5223-459b-b027-415f020087b9
let
	i = 10
	χ = χs[i] 
	ts = 0.0:dt:tf
	single_fluc = χ.(ts)
	plot(ts, single_fluc, xlabel = "t (μs)", ylabel="fluctuator sign", legend=:none, title=string("γ = 2π * ",fi[i], " MHz"), titlefontsize=13, marker=true, linealpha=0.5)
end

# ╔═╡ 16a6fa5f-21bb-4017-a881-a0a717079473
χtotal(t) = sum(map(χ -> χ(t), χs)) / N

# ╔═╡ 93ce97c6-9cf9-4459-8a97-3e6f2530aae8
χtotal(0.5)

# ╔═╡ c3a7a03c-a271-4c1d-be92-a14eb4b4317d
begin
	# ts = 0.0:dt:tf
	χtotal_traj = χtotal.(ts)
	plot(ts, χtotal_traj, xlabel = "t (μs)", ylabel="fluctuator sign", legend=:none, title=string("total fluctuation"), titlefontsize=13)
end

# ╔═╡ 03c48f2e-9a2f-4e8e-98f7-6a0d1da52a70
md"""
# Verifying the noise spectrum
"""

# ╔═╡ cae768c5-6417-4bfc-9937-507fe1948804
md"""
## PSD of single fluctuator time series data
"""

# ╔═╡ 4da6fe14-be74-45a8-97a0-54d170e7a9ff
md"""
$S(f) = \frac{F(f)^2}{2 \Delta f}$

where $S(f)$ is the power spectral density (PSD), $F(f)$ is the FFT result, and $\Delta f$ is the frequency spacing on the x-axis.
"""

# ╔═╡ 055cb559-2727-48e9-a858-51325f6f7c27
md"""
d = $(@bind d Slider(10 .^ (3:0.1:10)))
"""

# ╔═╡ f3734dac-3c75-4b32-8965-49861230790b
begin	
	t0 = 0.0             # Start time 
	fs = 1/dt          # Sampling rate (MHz)
	tmax = tf          # End time       
	
	t = t0:1/fs:tmax;   
	signal = single_fluc
	
	F = fftshift(fft(signal))
	freqs = fftshift(fftfreq(length(t), fs))
	positive_freqs = range(0, last(freqs), step=step(freqs))
	S(f) = (1/π) * (2γ) / ((2γ)^2 + (2π * f)^2) # Lorentzian PSD

	df = step(freqs)
	PSD = map(f -> abs(f)^2 / 2df, F)
	
	# plots 
	PSD_plot = plot(freqs, PSD, title = "power spectral density", xlim=(0, 20), legend=:none, xlabel="f (MHz)", yaxis=(:log10)) 
	plot!(positive_freqs, d * S.(positive_freqs), xlims=[1,20], yaxis=(:log10))
end

# ╔═╡ ffda5446-0e7e-4e89-962e-3397250ab879
S

# ╔═╡ 3ed66f7b-011f-4657-88b2-dc4a0817310e
md"""
d = $d
"""

# ╔═╡ 1d26e674-8e9a-4d19-9377-4ddf31bfa1f0
md"""
I expected a Lorentzian power spectrum as plotted [above](#8ccc5fa4-493a-4e17-9860-9a2fc1500857), so not sure how to interpret this.
"""

# ╔═╡ 1a874bc4-a1de-44a6-81c4-f2035693b7f0
md"""
## PSD of $N fluctuators
"""

# ╔═╡ dc88b8a0-e18b-4c50-9272-39d1176939b7
md"""
c = $(@bind c Slider(10 .^ (1:0.1:7)))
"""

# ╔═╡ bbe4c2e3-18e8-42bb-a460-23853e2f5b1d
c

# ╔═╡ 273ea734-a6c3-4c26-b233-216e3b24bb79
let	
	t0 = 0.0             # Start time 
	fs = 1/dt          # Sampling rate (MHz)
	tmax = tf          # End time       
	
	t = t0:1/fs:tmax;   
	signal = χtotal_traj
	
	F = fftshift(fft(signal))
	freqs = fftshift(fftfreq(length(t), fs))

	df = step(freqs)
	PSD = map(f -> abs(f)^2 / 2df, F)
	
	# plots 
	PSD_plot = plot(freqs, PSD, title = "power spectral density", xlim=(0, 20), legend=:none, xlabel="f (MHz)", yaxis=(:log10)) 
	plot!(positive_freqs, c ./positive_freqs, xlims=[1,20], yaxis=(:log10))
end

# ╔═╡ f9de6987-c874-4410-9d9a-d58ab72dc820
md"""
## FFTW checks
"""

# ╔═╡ 3733851e-f6ad-49cd-885c-52da3cb26b42
let
	N = 21
	xj = (0:N-1)*2*π/N
	f = 2*exp.(17*im*xj) + 3*exp.(6*im*xj) + rand(N)
	
	original_k = 1:N
	shifted_k = fftshift(fftfreq(N)*N)
	
	original_fft = fft(f)
	shifted_fft = fftshift(fft(f))
	
	p1 = plot(original_k,abs.(original_fft),title="Original FFT Coefficients", xticks=original_k[1:2:end], legend=false, ylims=(0,70));
	p1 = plot!([1,7,18],abs.(original_fft[[1,7,18]]),markershape=:circle,markersize=6,linecolor="white");
	p2 = plot(shifted_k,abs.(shifted_fft),title="Shifted FFT Coefficients",xticks=shifted_k[1:2:end], legend=false, ylims=(0,70));
	p2 = plot!([-4,0,6],abs.(shifted_fft[[7,11,17]]),markershape=:circle,markersize=6,linecolor="white");
	plot(p1,p2,layout=(2,1))
end

# ╔═╡ 7fd96fc7-09eb-4a55-a6ea-a98c634fd06c
let
		
	t0 = 0              # Start time 
	fs = 44100          # Sampling rate (Hz)
	tmax = 0.1          # End time       
	
	t = t0:1/fs:tmax;   
	signal = sin.(2π * 60 .* t)
	
	F = fftshift(fft(signal))
	freqs = fftshift(fftfreq(length(t), fs))
	
	# plots 
	time_domain = plot(t, signal, title = "Signal", label='f',legend=:top)
	freq_domain = plot(freqs, abs.(F), title = "Spectrum", xlim=(-100, +100), xticks=-100:20:100, label="abs.(F)",legend=:top) 
	plot(time_domain, freq_domain, layout = (2,1))
end

# ╔═╡ c06ee090-03e9-45a0-ba95-7baf29ad5daf
md"""
# Resources
"""

# ╔═╡ 82f94d63-8a44-4fdb-b1ba-fe85412f10bf
md"""
## Math
"""

# ╔═╡ 2b4e715a-c850-403f-8b01-48611560c1e8
md"""
## Physics
"""

# ╔═╡ 18d49165-c1b1-4a84-874a-85318515f668
md"""
## Julia
"""

# ╔═╡ 84ba819d-f07c-4d3c-9083-f969cc7587af
begin
	Modeling_time_btwn_events = html"<a href='https://medium.com/geekculture/how-to-model-time-between-events-using-the-exponential-gamma-and-poisson-distributions-4b058a357a55'>Modeling time between events, Riveroll, 2021 📘</a>"

	Bergli_et_al = html"<a href='https://arxiv.org/abs/0904.4597'> Bergli, et al. 📘</a>"

	julia_fft = html"<a href='https://www.matecdev.com/posts/julia-fft.html'> Using the FFTW Library in Julia 📘</a>"

	converting_PSD = html"<a href='https://dsp.stackexchange.com/questions/25456/conversion-of-fft-to-psd'> Conversion of FFT to PSD 📘</a>"

	
	

	md" ♦️ **URLs**"
end

# ╔═╡ 73547c2c-5b39-4329-b74e-aaafcdc9bd25
md"""
As described in $Modeling_time_btwn_events, the Poisson and exponential distributions are useful for modeling the amount of time between rare events. Here, the rare event is the fluctuator flipping sign. The Poisson distribution gives the probability of the fluctuator flipping sign in a a given time interval, in this case, the time step of simulation.

The random variable of the Poisson distribution is $X$, the number of flips in the time interval $\tau$:

$P(X = k) = \frac{(\gamma \tau)^k}{k!}e^{-\gamma \tau}$

where $n = \gamma \tau$ is the expected number of fluctuator flips, i.e.

$n = \sum_{k = 0}^{\infty} k P(X = k).$
"""

# ╔═╡ a03ec329-7017-40dd-9b6e-aa5b35d7afbe
md"""
There are many articles discussing modeling 1/f noise using fluctuators, and the relevance to quantum computing. However, the clearest reference I could find is $Bergli_et_al, so I will base my discussion here on that and include other references in the [Resources section](#c06ee090-03e9-45a0-ba95-7baf29ad5daf).
"""

# ╔═╡ 39b738ac-3f5b-4fb5-8f24-f63a755b1e41
md"""
As described in $Bergli_et_al, a single fluctuator generates a Lorentzian power spectrum

$S(\omega) = |v|^2 \frac1{\pi} \frac{2 \gamma}{(2\gamma)^2 + \omega^2}$
"""

# ╔═╡ 9e9f5d0e-a60b-4708-9f36-023a76ec1e70
md"""
To model multiple fluctuators, we sample them log-uniformly between a minimum and maximum frequency. As pointed out in $Bergli_et_al, this is important for generating a 1/f noise spectrum out of individual fluctuators, which have a Lorentzian (not 1/f) spectrum.
"""

# ╔═╡ 53f42e42-34ec-499e-a998-bb6775c55c12
md"""
Below, I compute the FFT of the single fluctuator time series data [simulated above](#0a0a5f1c-7a2f-4a04-bee5-f368c1b12900), then convert it via the following formula based on $converting_PSD.  Not sure if this is the right approach.
"""

# ╔═╡ f7d6ada7-7177-4344-8dab-2b6470e8f88e
md"""
See $julia_fft for details.
"""

# ╔═╡ 42de6802-f290-4950-9e45-b265a4221bbd
md"""
$converting_PSD
"""

# ╔═╡ 7c503194-d95b-4850-a2ae-ebcbb354cf45
md"""
 $Modeling_time_btwn_events 
"""

# ╔═╡ 217d0753-dabb-498d-a5f2-68e2b63f5922
md"""
$Bergli_et_al 
J. Bergli, Y. M. Galperin, and B. L. Altshuler, Decoherence in Qubits Due to Low-Frequency Noise, New J. Phys. 11, 025002 (2009).

"""

# ╔═╡ Cell order:
# ╟─a52844ea-9a74-42b0-b37a-7cc8cd79a0a6
# ╟─01c6cbbc-1945-40bb-be1a-31abefa38bba
# ╟─e625f413-e93b-48c9-b14d-181c37b0fe0a
# ╟─249f33ea-92aa-4b9c-bdfa-d35fe5254225
# ╠═3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─4e92db51-ef23-4a88-ba33-92771472ff2f
# ╟─4c838a57-4cf0-4218-980e-102359996fb9
# ╟─3a8257a5-6997-4f37-b4bc-170c722bbcfd
# ╟─73547c2c-5b39-4329-b74e-aaafcdc9bd25
# ╟─48a30585-f3ae-4cc4-ad58-fb5b07b26795
# ╟─25b85405-97d2-4a6e-ac77-be9823b3193f
# ╟─49c97425-9008-4109-8544-db805396cdc4
# ╟─a3940d5a-5b76-4224-bc4b-1cae947505e6
# ╟─727c5a36-1f56-43be-b75a-f7cccc910a4d
# ╟─21a93094-5aac-4097-b8af-e4e974f49dcf
# ╟─d135dad7-8434-4c56-88f8-4feeafe3eb11
# ╟─2be13c7c-1f23-4d70-87ae-f997c8ec050d
# ╟─b21cfda7-a6f9-4200-8676-b089c952f0c4
# ╟─a03ec329-7017-40dd-9b6e-aa5b35d7afbe
# ╟─8ccc5fa4-493a-4e17-9860-9a2fc1500857
# ╟─495e2820-4029-4c95-a736-82e258ebedb7
# ╟─4465eafa-85be-44ca-b5ee-3f62dcf83cf8
# ╟─39b738ac-3f5b-4fb5-8f24-f63a755b1e41
# ╟─6d7d40c5-1b1c-48b3-a6a7-635a628ee488
# ╟─0a0a5f1c-7a2f-4a04-bee5-f368c1b12900
# ╟─b01b386c-d76b-4383-851e-9316a824b994
# ╟─d3aa8364-5f56-442f-8990-e4ccad33ead6
# ╠═70648025-b722-4920-a876-c22dfa8f1d9a
# ╠═91507a5a-dc57-4b15-9fd4-309d899210d0
# ╠═72b67c33-12b0-4d04-b2aa-b18fc6faeb4a
# ╟─cf7c7b7b-8adb-49a7-b245-83d440d2b9c5
# ╟─acf8d55f-2676-4114-8163-c6cb49a8a8a4
# ╟─1e6c50f7-31dd-4e61-83a9-b7a56b5c5f13
# ╟─9e9f5d0e-a60b-4708-9f36-023a76ec1e70
# ╟─fdd591eb-1934-4cc3-9a7f-afcef3031972
# ╠═8631b550-01aa-438a-b659-8818bddfa4cb
# ╠═d001c6b6-5d39-4af3-b410-18d9e0a2ad57
# ╠═b1d4e9d2-8d9f-4729-8078-564391fcd817
# ╟─9e067f48-df19-41f2-b76a-6222dd22a9b6
# ╟─49b393f8-28a9-491d-aa6a-2e62bfa29de1
# ╠═1a1e2154-b0ed-4d21-8456-92ac61c735d7
# ╠═93ce97c6-9cf9-4459-8a97-3e6f2530aae8
# ╠═da396e7c-5223-459b-b027-415f020087b9
# ╠═16a6fa5f-21bb-4017-a881-a0a717079473
# ╠═c3a7a03c-a271-4c1d-be92-a14eb4b4317d
# ╟─03c48f2e-9a2f-4e8e-98f7-6a0d1da52a70
# ╠═bd133c05-a98e-4b71-86d3-24c4519e9038
# ╟─cae768c5-6417-4bfc-9937-507fe1948804
# ╟─53f42e42-34ec-499e-a998-bb6775c55c12
# ╟─4da6fe14-be74-45a8-97a0-54d170e7a9ff
# ╠═f3734dac-3c75-4b32-8965-49861230790b
# ╠═ffda5446-0e7e-4e89-962e-3397250ab879
# ╠═055cb559-2727-48e9-a858-51325f6f7c27
# ╠═3ed66f7b-011f-4657-88b2-dc4a0817310e
# ╟─1d26e674-8e9a-4d19-9377-4ddf31bfa1f0
# ╟─1a874bc4-a1de-44a6-81c4-f2035693b7f0
# ╠═dc88b8a0-e18b-4c50-9272-39d1176939b7
# ╠═bbe4c2e3-18e8-42bb-a460-23853e2f5b1d
# ╠═273ea734-a6c3-4c26-b233-216e3b24bb79
# ╟─f9de6987-c874-4410-9d9a-d58ab72dc820
# ╟─f7d6ada7-7177-4344-8dab-2b6470e8f88e
# ╟─3733851e-f6ad-49cd-885c-52da3cb26b42
# ╟─7fd96fc7-09eb-4a55-a6ea-a98c634fd06c
# ╟─c06ee090-03e9-45a0-ba95-7baf29ad5daf
# ╟─82f94d63-8a44-4fdb-b1ba-fe85412f10bf
# ╟─42de6802-f290-4950-9e45-b265a4221bbd
# ╟─7c503194-d95b-4850-a2ae-ebcbb354cf45
# ╟─2b4e715a-c850-403f-8b01-48611560c1e8
# ╟─217d0753-dabb-498d-a5f2-68e2b63f5922
# ╟─18d49165-c1b1-4a84-874a-85318515f668
# ╟─84ba819d-f07c-4d3c-9083-f969cc7587af
# ╠═3661abf7-84e6-4e39-854b-ad7741a8ff36
