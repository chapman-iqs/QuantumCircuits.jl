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

# ╔═╡ edf82b03-419b-41c3-9bae-a348658b57a2
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
	using Plots
	using DataFrames
	using Dates
	using QuantumCircuits

	include("utilities/single-qubit-operators.jl")
	include("utilities/utilities.jl")
	include("utilities/plotting.jl")

	include("notebooks/table-of-contents.jl")
	include("resources.jl")

	md" # Packages and julia files"
	
	
end

# ╔═╡ 5996b231-cfbd-4d92-82fb-d9eccc0898df
md"""
In this interactive notebook, we'll compare quantitative and qualitative features of the two primary integration methods implemented in `QuantumCircuits.jl` : `bayesian` and `rouchon`.
"""

# ╔═╡ c095d9bb-6360-42bc-be5d-47c8db202598
mdp(table_of_contents📔)

# ╔═╡ 539e06ce-d31d-42f1-963c-f83a3b9fdfd5
TableOfContents(title="CQED Integration methods")

# ╔═╡ 5d30c7c1-4c53-40d2-8e73-b58db5697cfa
md"""
# Description of the methods

"""

# ╔═╡ c257cccb-efbe-4302-bf7a-e8e2798d7263
md"""
## Usage
"""

# ╔═╡ e9e8b681-3e90-4f6f-b3a0-6d7f5bff918b
md"""
API is the same for both methods.

Choose the simulation method:
$(@bind solve_str Select(["bayesian", "rouchon"]))

"""

# ╔═╡ 385a0296-e3b0-4178-ab1a-212c93ca3324
solve = (solve_str == "bayesian") ? bayesian : rouchon

# ╔═╡ 639931ed-0520-499c-ba1c-90b04e854ebd
begin
	ΩR  = 2π # Rabi frequency
	
	Γ = 0.15 # Measurement dephasing rate, (MHz)
	η = 0.3
	
	H = (ΩR / 2) * σx
	J = [(σz, (1 - η) * Γ)] 
	C = [(σz, Γ, η)]
	
	ψ0 = g
	dt = 1e-3 # integration time-step
	tf = 10.0 # μs

	sol = solve((0, tf), ψ0, H, J, C; dt=dt)
end

# ╔═╡ 8ec19bbd-8ea6-44bb-bd57-b1638d519cd0
md"""
Both solvers return a `QuantumCircuits.Solution` object, a struct with fields `(t::Vector{Timescale}, ρ::Vector{State}, r::Vector{Record})`.
"""

# ╔═╡ d4e7b79b-eb8c-48b3-9a94-650481fa8cb7
md"""
##### Commonalities

The `bayesian` and `rouchon` methods have important commonalities. 

Both:
* Model quantum trajectories of weakly and continuously measured open quantum systems.

* Take `H`, `[J]`, `[C]` (Hamiltonian, Lindblad, collapse) operators as inputs.

* If simulating quantum trajectories, can generate a noisy record of a monitored variable $r = \langle x \rangle + dW/dt$, with Wiener process $dW \sim \mathcal{N}(0, \sqrt{\tau_m})$ if running a simulation. If reconstructing experimental trajectories, can process an imported experimental record.

"""

# ╔═╡ b6c73eea-d17c-49d4-a091-be10c8134718
md" ## Bayesian method "

# ╔═╡ 4fc42181-a1d4-46d9-be76-7277b9a2d6a5
md" #### Update equations"

# ╔═╡ 343161e9-e975-45b5-a444-b9f9f0b30cfb
md"""
Bayesian filtering describes the evolution of a weakly measured quantum system in fundamentally discrete chunks of time, $\Delta t$. In general, the system is described by a density matrix $\rho$. It undergoes **unitary evolution** according to a  Hamiltonian $H$, such that

$\hat \rho' = \hat U \rho {\hat U}^\dagger,$

where 

$\hat U = \exp(- i \Delta t \hat H).$

(Note that $\hbar$ has been absorbed into $\hat H$ here.) 


"""

# ╔═╡ d5af31e4-ba93-4f3f-9c55-4f2ce600d797
md"""
The **weak measurement** of an observable $\hat C$ with ensemble measurement dephasing rate $\Gamma_m$ is represented by a Kraus operator $\hat M_{r_C}$, such that

$\hat M_{r_C} = \exp\Big[\frac{\Delta t}{2 \tau_m} {r_C}^* \hat C  - \frac14 (\hat C + \hat C^\dagger)^2\Big].$

"""

# ╔═╡ 97e2560f-3dd9-45a0-85f5-aa9e12800292
md"""

Here, $r_C$ is a complex noisy record associated with measuring $\hat C$:

$r_C = \braket{\hat C} + d\zeta, \hspace{5mm} d\zeta \sim \mathcal{N}(0, \tau_m / dt)$

while $\tau_m = 1/2\Gamma_m \eta$ is the timescale of measurement collapse. Each $\hat M_{r_C}$ corresponds to a POVM (positive operator-valued measure) element 

$E_{r_C} = \hat M_{r_C}^\dagger \hat M_{r_C},$ 

while the $\{E_{r_C}\}_{r_C}$ form a POVM such that

$\sum_{r_C} E_{r_C} = \mathbb{1}.$

(There are, of course, infinitely many possible records $r_C$ corresponding to different measurement records.) 

The weak measurement updates the state according to

$\hat \rho '' = \hat M_{r_C}^\dagger \rho' \hat M_{r_C}.$

"""


# ╔═╡ d6d6e78b-15cb-4c6a-af30-8c42b35dc1b9
md"""
Finally, given Lindblad operators $\hat J$ with rates $\kappa_J$, the system **decoheres** according to the corresponding Kraus operators

$M_{\text{decay, } \hat J} = \sqrt{\kappa_J \Delta t} \hat J$

and

$M_{\text{null, } \hat J} = \sqrt{\mathbb{1} - \kappa_J \Delta t \hat J^\dagger \hat J}$

satisfying

$M_{\text{decay, } \hat J}^\dagger M_{\text{decay, } \hat J} + M_{\text{null, } \hat J}^\dagger M_{\text{null, } \hat J} = \mathbb{1},$

which is the completeness relation for the Kraus map. The state updates as

$\rho''' = M_{\text{decay, } \hat J} \rho'' M_{\text{decay, } \hat J}^\dagger + M_{\text{null, } \hat J} \rho'' M_{\text{null, } \hat J}^\dagger.$

"""

# ╔═╡ 6315198d-4962-4686-b436-80a3525e0dca
md"""
The final state after increment $\Delta t$ is obtained after renormalizing by the trace:

$\rho(t + \Delta t) = \frac{\rho'''}{\text{Tr}( \rho''')}.$

"""

# ╔═╡ 40a4c411-745d-4ac9-8cd8-ab4a66979fa9
md"""
##### Assumptions

Bayesian filtering assumes small $\Delta t$ in the interleaving of unitary, measurement, and decay evolution. However, it does not assume small $\Delta t$ in pure measurement evolution. Thus, in the absence of unitary evolution and dissipation (i.e. pure measurement evolution), the Bayesian method yields an exact result that is robust to choice of $\Delta t$.
"""

# ╔═╡ a306838f-f883-4fe9-94e6-8d5771cd9793
md" ## Rouchon method "

# ╔═╡ f66192ba-ff34-499e-ae8e-e8d08b91b9a1
md""" 
The Rouchon method was presented in [1] as a numerical integration method that provides an efficient alternative to the stochastic master equation (SME). Rouchon's method expands to second order in the Wiener increment $\Delta W_r$. As a result, the method preserves positivity of the conditioned quantum state, unlike directly integrating the SME.

"""

# ╔═╡ a541146f-9143-4a91-bb3a-be17e9bb6561
md"""
$M(t) = \mathbb{1} - \Big (i H + \frac12 \sum_j J_j^\dagger J_j ) dt + \sum_c \sqrt{\eta_c} C_c dy_c(t) +  \sum_{c,d} \frac{\sqrt{\eta_c \eta_d}}2 C_c C_d \big(dy_c(t) dy_d(t) - \delta_{c,d} dt \big)$
"""

# ╔═╡ 88179c16-f111-462d-a9d6-72e3764715eb
md"""
$dy_c(t) = \sqrt{\eta_c} \text{Tr} \big ( C_c \rho(t) + \rho(t) C_c^\dagger \big ) dt + dW_c(t)$
"""

# ╔═╡ 6e92e073-d122-4a84-8c4c-b19858c173a9
md"""
$D(t) = \Big (\sum_j J_j \rho(t) J_j^\dagger  + \sum_c C_c \rho(t) C_c^\dagger \Big ) dt$
"""

# ╔═╡ 501a6621-8054-4fa3-b9a7-f1a921e14a56
md"""
$\rho(t + dt) = \frac{M(t) \rho(t) M(t)^\dagger + D(t)}{\text{Tr} \Big ( M(t) \rho(t) M(t)^\dagger + D(t) \Big )}$
"""

# ╔═╡ cb87a45d-4ff1-4da9-8bbe-ec2c193cb7c6
md"""
The Rouchon method is similar to the Bayesian method in that it uses completely positive maps, and therefore preserves positivity of the density matrix. However, it is still a first-order approximation in $dt$. Whereas the Bayesian method will be exact when measurement is the only source of system evolution, the Rouchon method is not exact under any circumstance.
"""

# ╔═╡ 676e029d-1141-41b2-bc81-5d2e4500c684
md"""
# Method comparison
"""

# ╔═╡ df840a53-c8b0-4215-97de-ef6f1c720ed7
md" ## Coarse-graining "

# ╔═╡ 7ffe3ce0-2ce5-4012-8679-93b9154b150e
md"""
In order to compare across time-step sizes, we first must generate consistent records corresponding to various time-steps. In weak measurement, $\Delta t$ has a physical interpretation as the integration time for the measurement. Thus, when we change $\Delta t$, we must both resample *and* rescale the measurement record in order to reproduce the same evolution.

"""

# ╔═╡ e41a9229-34f2-4a42-810b-7c95dc02780b
md"""
Let's start by looking at the record generated by $solve_str in our example at the beginning of the notebook:
"""

# ╔═╡ 4f3c3104-c64d-4b77-8932-c37f4aaabaa1
begin
	tt = sol.t
	rs = sol.r[1]
end

# ╔═╡ 3e0a12d2-b4bf-4424-8b62-d9a714e393fc
md"""
Note that `bayesian` and `rouchon` both take in / output records of the form 
	
$r_i = \langle c_i \rangle + \zeta_i,$

where 

$\zeta_i \sim \mathcal N\Big(0, \sqrt{\frac{\tau_i}{dt}}\Big)$ 

is a normally distributed stochastic variable with variance $\tau_i / dt$. This record corresponds to measurement of the $i^{th}$ measurement operator $c_i$, measured at rate $\Gamma_i = 1/(2\tau_i)$. (Note that the stochastic variable can be related to a Wiener increment $dW \sim \mathcal N(0, \sqrt{dt})$ by the relation $\zeta_i = \frac{dW}{dt} \sqrt{\tau_i}$.) 

Suppose we want to reconstruct a trajectory using timestep $dt'$, but our record was sampling with time bins $dt < dt'$. Since increasing $dt$ decreases the variance of the record, we cannot simply sample the old record at the new timescale. We must also coarse-grain the old record, which naturally decreases the variance of the zero-mean Gaussian noise.

"""

# ╔═╡ 0b0f2aa0-431d-4029-b7c2-911aa4d2f69a
md"""
n = 1 $(@bind n Slider(1:1:20, default=5)) 20 
"""

# ╔═╡ b29c82ed-2c3e-4d19-a9ec-9d704b181edb
md"""

The function `subseries` (in Utilities section) implements this. (You can `command-click` a function to jump to its definition.) The code below coarse-grains the record at a scale of n = $n, which creates a new record corresponding to measurement bins of width

```math
dt' = n \hspace{1mm} dt,
```
or $dt'$ =  $(dt*n) μs.
"""

# ╔═╡ 69eaed7c-392d-40d8-9f0c-63b4bc4aa18d
md"""
You can try playing with the coarse-graining scale below:

n = $n

"""

# ╔═╡ 012e9206-f306-4776-a546-d6f7dcf00ef2
md"""
We can look at these two records as time-series or as histograms:
"""

# ╔═╡ 43232217-a1e6-4630-b3ef-d339f10f8e2f
τm = 1/(2Γ*η)

# ╔═╡ 275ffd8d-274b-4fab-946e-d9482fa340a6
md"""
Since the record is Gaussian-distributed, time-smoothing the record using `subseries` reduces the variance of the distribution. We can check that the variances reported in the plot above are close to the theoretical values at each $dt$:

fine-grained: $\sigma = \sqrt{\tau_m / dt}$ = $(round(sqrt(τm/dt),digits=3))

coarse-grained: $\sigma = \sqrt{\tau_m / (n \hspace{1mm} dt)}$ = $(round(sqrt(τm/(dt*n)),digits=3))
"""

# ╔═╡ e0cb7572-3d67-47e0-9482-70588a693d41
md" ## Time-step comparison"

# ╔═╡ f870609e-f3bc-40e4-90fb-484e7b52db4b
md"""
First we'll get the subseries of the fine-grained simulation:
"""

# ╔═╡ 3f7d824c-8a0e-4547-bb06-81717a01df40
md"""
We also need the subrecord to run a coarse-grained trajectory:
"""

# ╔═╡ 0df1c1f4-387a-4d36-a818-10e75540e762
md"""
Now we can run the trajectory reconstructions:
"""

# ╔═╡ 7c434d23-a5ed-41bb-8dc7-8b1440a0680b
md"""
Now we can compare. Change the coarse-graining scale below and see how `bayesian` reconstruction is affected. (For small $n_s$, you may not see the difference visually!)
"""

# ╔═╡ 82b18884-28d4-4ecb-b21c-d3e57dd64870
md"""
 $n_s$ = 1 $(@bind ns Slider(1:150, default=40)) 150
"""

# ╔═╡ e47ee9c6-c9d6-44b8-8fc3-e6ebb58f28e3
md"""
 $n_s$ = $ns
"""

# ╔═╡ a0bb68b2-edde-4941-b034-c1fe8c7141ae
md" ## Exactness of Bayesian method for H = 0 "

# ╔═╡ 9af5502b-3639-419a-b113-64c11e70ead7
md"""
In the absence of coherent evolution (i.e., $H = 0$), the Bayesian method is exact. Thus, trajectories should match exactly regardless of coarse-graining.
"""

# ╔═╡ a0cb40a2-aa44-405d-9fbb-bd5e59bed6b3
md"""
 $n_0$ = 1 $(@bind n0 Slider(1:150, default=40)) 150
"""

# ╔═╡ 64767210-f6d5-4d01-88bb-0e25e327fd83
md" $n_0$ = $n0 "

# ╔═╡ e7ffdfb4-622e-4248-a71b-907ea1cbb8a5
md"""

# References

[1] P. Rouchon and J. F. Ralph, Efficient Quantum Filtering for Quantum Feedback Control, Phys. Rev. A 91, 012118 (2015).

"""

# ╔═╡ 2b7aa36e-e64a-4839-875e-2c472763cb80
md" # Utilities "

# ╔═╡ 26d5f3cc-a2c0-4530-a7e0-8c4fc80e9585
md"""
## Plotting
"""

# ╔═╡ 00c42bae-6a45-4aeb-8ee9-88881b2f23e8
# Plotting
function plot_timeseries(series...; plot_title="time series", xlabel="t (μs)", ylabel="arbitrary units", labels=[], colors = palette(:lightrainbow), kwargs...)

	label(i) = i > length(labels) ? i : labels[i]

	p = plot(size=(600,300))

	# Plot records vs. time
	for (i, (tt, ser)) in enumerate(series)
		plot!(tt, ser, color=colors[i], label=label(i), xlabel=xlabel, ylabel=ylabel, title=plot_title, legend=:outerright, kwargs...)
	end
	
    p
	
end

# ╔═╡ 419ab808-4918-445c-9087-cb7b33a733d6
# Plotting
function plot_timeseries(tt::Vector{Timescale}, series...; plot_title="time series", xlabel=L"$t$", ylabel="arbitrary units", labels=[], colors = palette(:lightrainbow), kwargs...)

	ser_colors(i) = colorpairs ? colors[i] : colors[2i]
	label(i) = i > length(labels) ? i : labels[i]

	p = plot(size=(600,300))

	# Plot records vs. time
	for (i, ser) in enumerate(series)
		plot!(tt, ser, color=ser_colors(i), label=label(i), xlabel=xlabel, ylabel=ylabel, title=plot_title, legend=:outerright, kwargs)
	end
	
    p
	
end

# ╔═╡ 275e327e-1f0e-437f-9d62-7aa88b3e812b
# Plotting
function plot_records(series; plot_title="time series", xlabel="t (μs)", ylabel="arbitrary units", labels=[], colors = palette(:lightrainbow), histograms=false)

	label(i) = i > length(labels) ? i : labels[i]

	p = plot(size=(600,300))

	# Plot records vs. time
	for (i, (tt, ser)) in enumerate(series)
		plot!(tt, ser, color=colors[i], xlabel=xlabel, ylabel=ylabel, title=plot_title, legend=:none)
	end

	if !histograms
	    return p
		
	else
		
		h =	histogram(title="record histograms", xlabel="probability", orientation=:h)
		for (i, (tt, ser)) in enumerate(series)
			(μ, σ) = params(fit(Normal, ser))
			histogram!(ser, bins = :scott, normalize=true, color=colors[i], alpha=0.65, label=string(labels[i], ", σ = ", round(σ, digits=3)), orientation=:h, legend=:bottomright)
		end

		l = @layout [timeseries{0.5w} histogram{0.5w}]
		return plot(p, h, layout=l, link=:y)
		
	end
	
end

# ╔═╡ 36739d52-b59f-41a6-95e1-10676a1dc691
md"""
#### Colors
"""

# ╔═╡ 51fd2fa1-1a8f-4f29-a600-bfd37ce9b26a
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

# ╔═╡ 0d4f6084-f52f-4e07-b28c-26447f8e2368
function qubit_plot(sol::Solution; record=false, title="")

	basis = qbasis
	colors = colorsq1
	labels = qlabels

	exps = map(op -> expectations(sol, op), basis)
	p = 0.5 .* (1 .+ exps[1].^2 .+ exps[2].^2 .+ exps[3].^2)

	pl = plot(size=(600,300), legendfontsize=10, titlefontsize=12, legend=:outerright, ylabel="bloch coordinates", xlabel = record ? "" : "t (μs)", linewidth=1.2, title=title)

	for l in 1:length(basis)
		label = labels[l]
		color = colors[l]
		exp = exps[l]
		plot!(sol.t, exps[l], color=color, label=label, legend=:outerright, ylims=[-1,1])
	end

	plot!(sol.t, p, color=colors[4], label=L"Tr(\rho^2)")

	if !record
		return pl
	else
		l = @layout [blochs{0.6h}; record{0.4h}]
		pr = plot(sol.t, sol.r[1], color=colors1[1], xlabel="t (μs)", ylabel="record", label=:none, legend=:outerright, title="")
		return plot(pl, pr, layout = l, link=:y, linewidth=1.2)
	end
end

# ╔═╡ f9d3ac45-0fe1-4d0a-9df9-a03bd2e08327
function qubit_plot(sol1::Solution, sol2::Solution; record=false, title="", color1=colorsq1, color2=colorsq2, l1="", l2="")

	basis = qbasis
	labels = qlabels

	exps1 = map(op -> expectations(sol1, op), basis)
	p1 = 0.5 .* (1 .+ exps1[1].^2 .+ exps1[2].^2 .+ exps1[3].^2)

	exps2 = map(op -> expectations(sol2, op), basis)
	p2 = 0.5 .* (1 .+ exps2[1].^2 .+ exps2[2].^2 .+ exps2[3].^2)

	pl = plot(size=(600,300), legendfontsize=10, titlefontsize=12, legend=:outerright, ylabel="bloch coordinates", xlabel = record ? "" : "t (μs)", linewidth=1.2, title=title)

	for l in 1:length(basis)
		plot!(sol1.t, exps1[l], color=color1[l], label=labels[l], legend=:outerright, ylims=[-1,1], linewidth=1.2)
	end
	plot!(sol1.t, p1, color=color1[4], label=L"Tr(\rho^2)", linewidth=1.2)
	
	for l in 1:length(basis)
		plot!(sol2.t, exps2[l], color=color2[l], label=:none, legend=:outerright, ylims=[-1,1], linestyle=:dash, linewidth=2)
	end
	plot!(sol2.t, p2, color=color2[4], label=:none, linestyle=:dash, title=string(l1, " ---, ", l2, " - - -"), linewidth=2)


	if !record
		return pl
	else
		l = @layout [blochs{0.6h}; record{0.4h}]
		pr = plot(sol1.t, sol1.r[1], color=mixed_colors[1], xlabel="t (μs)", ylabel="record", label=:none, legend=:outerright, title="")
		plot(sol2.t, sol2.r[1], color=mixed_colors[2], xlabel="t (μs)", ylabel="record", label=:none, legend=:outerright, title="")
		return plot(pl, pr, layout = l, link=:y)
	end
end

# ╔═╡ 81bd33a8-8970-427b-9541-8f0e244f57f2
qubit_plot(sol, record=true, title=string(solve_str, " method solution"))

# ╔═╡ 6bfe4b1a-cf3a-48fb-9782-88735f1546f0
md"""
## Other
"""

# ╔═╡ 8eb9d84d-e658-406e-8a2a-8a28cfee9b04
function subseries(rec, T, dt; scale=2)
	ts = collect(range(first(T), last(T), step=dt))
	tts = subselect(real(coarse_grain(ts; n=scale)); n=scale)
	(ti, tf) = (tts[1], tts[end])
	dtt = dt * scale
	
	subrec = subselect(real(coarse_grain(rec; n=scale)); n=scale)
	
	(tts, subrec)
	
end

# ╔═╡ 559d34ed-fcab-46d8-9fea-1774067b2d48
(tts, Rss) = subseries(sol.r[1], (0.0, tf), dt; scale=n)

# ╔═╡ 0da08e8a-0db3-4da1-96ff-531618dab458
begin
	rec_labels = ["dt = $dt μs", "dt' = $(round(dt*n, digits=3)) μs"]
	plot_records([(tt, rs), (tts, Rss)]; plot_title="records", labels=rec_labels, colors=mixed_colors, histograms=true)
end

# ╔═╡ 7baa7a79-090e-4318-982f-6c1982f82a58
function rms(ser1::Array, ser2::Array)
	l = min(length(ser1), length(ser2))
	sqrt(sum((ser1[1:l] - ser2[1:l]).^2))/l
end


# ╔═╡ 078c6b8e-4b74-46dd-aebd-77acf392c2c9
function blochs(sol::Solution)

	t = sol.t
	x, y, z = map(op -> expectations(sol, op), qbasis)
	p = 0.5 .* (x.^2 .+ y.^2 .+ z.^2)

	(t, x, y, z, p)	
end

# ╔═╡ 43717946-c7f2-4f57-b516-9db669167587
begin
	T = (0.0, tf)
	(tb,xb,yb,zb,ρb) = map(ser -> subseries(ser, T, dt; scale=ns)[2], blochs(sol))
end

# ╔═╡ 91c7a5dd-f2e1-425b-9008-43854097f966
(ttc, Rsc) = subseries(rs, T, dt; scale=ns)

# ╔═╡ 559ce94c-8276-4716-91ae-29f4bd973c44
Rsc

# ╔═╡ 719cc6a7-946e-46e7-b8f0-fc16894c7e3e
begin
	ψ00 = normalize(g + im * e)
	H0 = 0 * Iq

	global solH0 = bayesian(T, ψ00, H0, J, C; dt=dt)	
	Rs0 = solH0.r[1] # record output from bayesian
	(tss0, Rss0) = subseries(Rs0, T, dt; scale=n0)
	
	T0 = (tss0[1], tss0[end])
	global solH0_coarse = bayesian(T0, ψ00, H0, J, C; dt=dt*n0, records=[Rss0])	
end

# ╔═╡ 4468cd55-6ac2-4ab8-af99-e521ca7a1672
let
	rec_labels = ["dt = $(1000*dt) ns", "dt' = $(1000*round(dt*n0, digits=3)) ns"]
	qubit_plot(solH0, solH0_coarse; l1=rec_labels[1], l2=rec_labels[2])
end

# ╔═╡ 7cc75874-ce1f-49da-860b-62a0f7fe800f
begin
	Tc = (ttc[1], ttc[end])
	sol_coarse = solve(Tc, ψ0, H, J, C; dt=dt*ns, records=[Rsc])	
	
	(tbc,xbc,ybc,zbc,ρbc) = blochs(sol_coarse)
end

# ╔═╡ 462bad0d-3706-435b-9918-b1b1c2b5e320
md"""
The root-mean-squared error between the two $solve_str series is $(round(rms(yb[1:end-1], ybc), digits=5)).
"""

# ╔═╡ 3ad9077b-1d06-48fb-9b6d-c98266a40bfe
begin
	rec_labels_s = ["dt = $(1000*dt) ns", "dt' = $(1000*round(dt*ns, digits=3)) ns"]
	qubit_plot(sol, sol_coarse; l1=rec_labels_s[1], l2=rec_labels_s[2])
end

# ╔═╡ 093383f1-d7a5-48c6-9927-4e42d158cc3d
green(text; title="Note") = Markdown.MD(Markdown.Admonition("correct", title, [text]))

# ╔═╡ 5ed525a9-c4b6-44d6-8ebf-6c86190a8992
red(text; title="Note") = Markdown.MD(Markdown.Admonition("danger", title, [text]))

# ╔═╡ b966b51c-e9cd-4696-a3b1-15ecb9d3808e
tann(text; title="Note") = Markdown.MD(Markdown.Admonition("warning", title, [text]))

# ╔═╡ f08bf969-d03f-47f9-9f5a-cddd18f8b109
blue(text; title="Note") = Markdown.MD(Markdown.Admonition("note", title, [text]))

# ╔═╡ d7c1fd28-5780-4efd-b085-4f1f797a3f09
hint(text; title="Hint") = Markdown.MD(Markdown.Admonition("hint", title, [text]))

# ╔═╡ Cell order:
# ╟─5996b231-cfbd-4d92-82fb-d9eccc0898df
# ╟─c095d9bb-6360-42bc-be5d-47c8db202598
# ╟─539e06ce-d31d-42f1-963c-f83a3b9fdfd5
# ╟─5d30c7c1-4c53-40d2-8e73-b58db5697cfa
# ╟─c257cccb-efbe-4302-bf7a-e8e2798d7263
# ╟─e9e8b681-3e90-4f6f-b3a0-6d7f5bff918b
# ╟─385a0296-e3b0-4178-ab1a-212c93ca3324
# ╠═639931ed-0520-499c-ba1c-90b04e854ebd
# ╠═81bd33a8-8970-427b-9541-8f0e244f57f2
# ╟─8ec19bbd-8ea6-44bb-bd57-b1638d519cd0
# ╟─d4e7b79b-eb8c-48b3-9a94-650481fa8cb7
# ╟─b6c73eea-d17c-49d4-a091-be10c8134718
# ╟─4fc42181-a1d4-46d9-be76-7277b9a2d6a5
# ╟─343161e9-e975-45b5-a444-b9f9f0b30cfb
# ╟─d5af31e4-ba93-4f3f-9c55-4f2ce600d797
# ╟─97e2560f-3dd9-45a0-85f5-aa9e12800292
# ╟─d6d6e78b-15cb-4c6a-af30-8c42b35dc1b9
# ╟─6315198d-4962-4686-b436-80a3525e0dca
# ╟─40a4c411-745d-4ac9-8cd8-ab4a66979fa9
# ╟─a306838f-f883-4fe9-94e6-8d5771cd9793
# ╟─f66192ba-ff34-499e-ae8e-e8d08b91b9a1
# ╟─a541146f-9143-4a91-bb3a-be17e9bb6561
# ╟─88179c16-f111-462d-a9d6-72e3764715eb
# ╟─6e92e073-d122-4a84-8c4c-b19858c173a9
# ╟─501a6621-8054-4fa3-b9a7-f1a921e14a56
# ╟─cb87a45d-4ff1-4da9-8bbe-ec2c193cb7c6
# ╟─676e029d-1141-41b2-bc81-5d2e4500c684
# ╟─df840a53-c8b0-4215-97de-ef6f1c720ed7
# ╟─7ffe3ce0-2ce5-4012-8679-93b9154b150e
# ╟─e41a9229-34f2-4a42-810b-7c95dc02780b
# ╠═4f3c3104-c64d-4b77-8932-c37f4aaabaa1
# ╟─3e0a12d2-b4bf-4424-8b62-d9a714e393fc
# ╟─b29c82ed-2c3e-4d19-a9ec-9d704b181edb
# ╠═559d34ed-fcab-46d8-9fea-1774067b2d48
# ╟─69eaed7c-392d-40d8-9f0c-63b4bc4aa18d
# ╟─0b0f2aa0-431d-4029-b7c2-911aa4d2f69a
# ╟─012e9206-f306-4776-a546-d6f7dcf00ef2
# ╠═0da08e8a-0db3-4da1-96ff-531618dab458
# ╠═43232217-a1e6-4630-b3ef-d339f10f8e2f
# ╟─275ffd8d-274b-4fab-946e-d9482fa340a6
# ╟─e0cb7572-3d67-47e0-9482-70588a693d41
# ╟─f870609e-f3bc-40e4-90fb-484e7b52db4b
# ╠═43717946-c7f2-4f57-b516-9db669167587
# ╟─3f7d824c-8a0e-4547-bb06-81717a01df40
# ╠═91c7a5dd-f2e1-425b-9008-43854097f966
# ╟─0df1c1f4-387a-4d36-a818-10e75540e762
# ╠═7cc75874-ce1f-49da-860b-62a0f7fe800f
# ╠═559ce94c-8276-4716-91ae-29f4bd973c44
# ╟─7c434d23-a5ed-41bb-8dc7-8b1440a0680b
# ╟─82b18884-28d4-4ecb-b21c-d3e57dd64870
# ╟─e47ee9c6-c9d6-44b8-8fc3-e6ebb58f28e3
# ╟─462bad0d-3706-435b-9918-b1b1c2b5e320
# ╠═3ad9077b-1d06-48fb-9b6d-c98266a40bfe
# ╟─a0bb68b2-edde-4941-b034-c1fe8c7141ae
# ╟─9af5502b-3639-419a-b113-64c11e70ead7
# ╠═719cc6a7-946e-46e7-b8f0-fc16894c7e3e
# ╟─a0cb40a2-aa44-405d-9fbb-bd5e59bed6b3
# ╟─64767210-f6d5-4d01-88bb-0e25e327fd83
# ╠═4468cd55-6ac2-4ab8-af99-e521ca7a1672
# ╟─e7ffdfb4-622e-4248-a71b-907ea1cbb8a5
# ╟─2b7aa36e-e64a-4839-875e-2c472763cb80
# ╟─26d5f3cc-a2c0-4530-a7e0-8c4fc80e9585
# ╟─0d4f6084-f52f-4e07-b28c-26447f8e2368
# ╟─f9d3ac45-0fe1-4d0a-9df9-a03bd2e08327
# ╟─00c42bae-6a45-4aeb-8ee9-88881b2f23e8
# ╟─419ab808-4918-445c-9087-cb7b33a733d6
# ╟─275e327e-1f0e-437f-9d62-7aa88b3e812b
# ╟─36739d52-b59f-41a6-95e1-10676a1dc691
# ╠═51fd2fa1-1a8f-4f29-a600-bfd37ce9b26a
# ╟─6bfe4b1a-cf3a-48fb-9782-88735f1546f0
# ╟─8eb9d84d-e658-406e-8a2a-8a28cfee9b04
# ╟─7baa7a79-090e-4318-982f-6c1982f82a58
# ╟─078c6b8e-4b74-46dd-aebd-77acf392c2c9
# ╟─093383f1-d7a5-48c6-9927-4e42d158cc3d
# ╟─5ed525a9-c4b6-44d6-8ebf-6c86190a8992
# ╟─b966b51c-e9cd-4696-a3b1-15ecb9d3808e
# ╟─f08bf969-d03f-47f9-9f5a-cddd18f8b109
# ╟─d7c1fd28-5780-4efd-b085-4f1f797a3f09
# ╠═edf82b03-419b-41c3-9bae-a348658b57a2
