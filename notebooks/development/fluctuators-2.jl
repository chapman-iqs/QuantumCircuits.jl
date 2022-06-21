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

# ╔═╡ df539f73-7c8a-4d65-b909-0d94720f724f
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
	using ProgressMeter
	using LsqFit

	using DataFrames
	using CSV
	

	include("notebooks/table-of-contents.jl")
	include("notebooks/resources.jl")
	
	include("utilities/two-qubit-operators.jl")
	include("utilities/plotting.jl")
	include("utilities/utilities.jl")

	include("notebooks/notebook-utilities.jl")
	
	md" # Packages and julia files"
end

# ╔═╡ e625f413-e93b-48c9-b14d-181c37b0fe0a
md"""
In this notebook, we model dephasing using a fluctuator bath.
"""

# ╔═╡ 249f33ea-92aa-4b9c-bdfa-d35fe5254225
mdp(table_of_contents📔)

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents(title="Fluctuator bath")

# ╔═╡ 7cb4592d-4b47-41c6-9cb3-cca1033b6c9c
md"""
# Setup
"""

# ╔═╡ 81394981-ea8f-4718-b4b3-138546ca8a71
md"""
## System parameters
"""

# ╔═╡ 41954632-a8c7-49d1-9c7e-559c0b0c0556
begin
	dt = 1e-3  # integration time-step
	tf = 4.0 # simulation duration (μs)
	ωz = (2π) * 0.0 # qubit z rotation

	Ω = (2π) * 0.2 # fluctuator freq. (MHz)
	ψ0 = normalize(g + e) # initial state
end

# ╔═╡ 9d1860c9-0b12-4d7b-97d6-0ea51db03a94
md"""
## Simulating the fluctuators
"""

# ╔═╡ 23bc6da9-15fc-4ee4-93a9-1abf6dd207d3
nfluctuators = 10

# ╔═╡ 9771ac72-8daf-4323-8267-b91410b150ff
md"""
Sample $nfluctuators fluctuators log-uniformly between a minimum and maximum frequency:
"""

# ╔═╡ 176571f5-291c-40dd-964b-7aa59f1ed812
begin
	ωmin = 1/tf
	ωmax = 1/dt
	ωlog = range(log(ωmin), log(ωmax), length=nfluctuators)
	ωs = exp.(ωlog)
end

# ╔═╡ ec4a876c-c4b8-4f4b-9c60-9827b8ebe337
function sgn(τ, tf)
	nflips = Int64(floor(tf / τ)) # number of times fluctuator CAN flip in simulation
	series = rand((-1, 1), (nflips + 1)) # string of random signs of fluctuator

	t -> let
		index = Int64(floor(t/tf * nflips)) + 1
		return series[index]
	end
end

# ╔═╡ 3f4a410c-a6f8-45e4-b7d0-e10d583fa91f
begin
	ωi = ωs #range(1/tf, 16, step=0.2)
	τi = map(ω -> 1/ω, ωi)
	times = range(0.0, tf, step=dt)
	si(τi, tf) = map(τ -> sgn(τ, tf), τi)
end

# ╔═╡ 181fd41d-e256-4c78-a100-899acd89a066
md"""
# Qubit simulations
"""

# ╔═╡ 7bbba21f-8d77-4d2e-8a7b-ddeb6c291427
begin
	
	# fluctuator Hamiltonian
	Hf(τi, tf) = let si = map(τ -> sgn(τ, tf), τi)
					t -> map(s -> (Ω/2) * σz * s(t), si) end
	
	# total Hamiltonian --
	# -- a function of a particular realization of fluctuator Hamiltonian Hf
	H(Hf) = t -> sum(Hf(t))
	
end

# ╔═╡ 65a619fc-ff14-4426-9c45-6d749fd8b78c
md"""
## Single trajectory
"""

# ╔═╡ 6d939bd2-2689-449d-97ca-34f7bf838858
begin
	Hi1 = Hf(τi, tf)
	sol1 = bayesian((0, tf), ψ0, H(Hi1), [], []; dt=dt)
end

# ╔═╡ 689a12f4-41cd-4efc-a88c-92f2fff864f8
begin
	qp = qubit_plot(sol1, record=false, legendpos=:outerright)
	pf = plotfluctuators(ωi, times, Hi1, Ω, labels=:none)
	plot(qp, pf, layout=grid(2,1))
end

# ╔═╡ ada56af9-023f-4364-ab67-2f4de3e15e0c
md"""
## Ensemble
"""

# ╔═╡ 567ff65b-c70d-42fc-9f31-7373f543bd9e
md"""
### No decoupling
"""

# ╔═╡ cbaac8ff-3898-4413-96c9-2d07d46ca5b5
N = 1000

# ╔═╡ 846b569e-1440-4cd1-8a5f-0db962cc7d3a
md"""
Average over N = $N realizations. The strength of fluctuations is Ω = $(Ω/(2π)) MHz.
"""

# ╔═╡ dba94b12-aac9-4cf2-8714-6cbb75e473e9
sols = map(1:N) do i
	Hi = Hf(τi, tf) # new realization of fluctuator Hamiltonian
	bayesian((0, tf), ψ0, H(Hi), [], []; dt=dt)
end

# ╔═╡ 6c0d59cd-0f9d-463b-aa97-c1617ab0aecf
ρavg = [mean(map(sol -> dm(sol.ρ[i]), sols)) for i in 1:length(sols[1].t)]

# ╔═╡ a10ad7ca-873c-41cf-a4bc-edd1e4288b02
pur = map(ρ -> real(tr(ρ * ρ)), ρavg)

# ╔═╡ 210332f6-d26b-4f5b-85f9-20ac9555a97c
md"""
Ω = $(Ω/2π) MHz
"""

# ╔═╡ fd4b2d5c-0ed3-4cf3-bf8d-dbc1543d0a71
md" $(@bind anim_plot CheckBox()) Animate bloch average "

# ╔═╡ da766962-481b-4db5-940f-662225147988
md"""
### Decoupling
"""

# ╔═╡ 547046a1-6163-41e2-b397-d41253ab0b8f
ΩR = (2π) * 25.0 # qubit Rabi frequency (for decoupling)

# ╔═╡ 1aa49c31-f055-4fee-abd3-6af4e3593d12
begin
	Hq = (ΩR/2) * σy + (ωz/2) * σz # qubit Hamiltonian, used for decoupling
	Hd(Hf) = t -> Hq + sum(Hf(t)) 
end

# ╔═╡ 490899df-edcc-4a7d-87fa-e4d0b3b103d2
sols_d = map(1:N) do i
	Hi = Hf(τi, tf) # new realization of fluctuator Hamiltonian
	bayesian((0, tf), ψ0, Hd(Hi), [], []; dt=dt)
end

# ╔═╡ 7217689a-ecda-4714-95a9-400d85f11282
md" $(@bind anim_plot_2 CheckBox()) Animate bloch average "

# ╔═╡ 6e6bb8e1-892e-4fd4-8e48-20988c6eecba
qubit_plot(sols_d[4]; record=false, title="sample trajectory")

# ╔═╡ 51ae5fda-9444-42e4-b190-a52215c39c87
md"""
## Simple Lindblad dephasing
"""

# ╔═╡ 48c21723-679d-40eb-a226-2988fa448f61
md"""
### No decoupling
"""

# ╔═╡ aaf39f5e-aeb7-4b2b-a2e7-be950c723995
md"""

γϕ = 0.1 MHz
$(@bind γϕ Slider(0.1:0.01:1.0, default=0.64))
1.0 MHz

"""

# ╔═╡ 7175e140-44f8-4961-9e83-367c20d00f26
let
	# simulation parameters
	# γϕ = 0.3 # dephasing rate (MHz)
	
	# Kraus operators
	H = Iq
	J = [(σz, γϕ)]
	C = []
	
	global sol_Lind = bayesian((0, tf), ψ0, H, J, C; dt=dt)
end

# ╔═╡ 9858f17b-7e3a-4b3c-83e0-52a45d8869db
qubit_plot(sol_Lind; record=false)

# ╔═╡ c3bd9db7-a958-47f9-98b6-4f3bb1bddbc7
md"""
γϕ = $γϕ MHz     
"""

# ╔═╡ 7232c30f-ad84-4a17-bf04-d8191244436c
γϕ / (2π)

# ╔═╡ ed3be650-e080-4928-b3d2-24445da8b9d4
exps_lind = map(op -> expectations(sol_Lind, op), qbasis)

# ╔═╡ b781fb8f-d7a0-4502-8322-e3b8b0e1ce1c
md"""
### Decoupling
"""

# ╔═╡ 6d5f64a4-84cd-4204-8169-a1103b241ae0
let	
	# Kraus operators
	H = (ΩR/2) * σy
	J = [(σz, γϕ)]
	C = []
	
	global sol_Lind_d = bayesian((0, tf), ψ0, H, J, C; dt=dt)
	global sol_Lind_d2 = bayesian((0, tf), ψ0, H, J, C; dt=dt/100)
	
end

# ╔═╡ b0f750db-6c58-4939-9f42-8388ba8b9af4
begin
	exps_lind_d = map(op -> expectations(sol_Lind_d, op), qbasis)
	exps_lind_d2 = map(op -> expectations(sol_Lind_d2, op), qbasis)
	times2 = sol_Lind_d2.t
end

# ╔═╡ b6a49a95-8e21-4e81-8f39-737c7a628a0e
exps_lind_d2

# ╔═╡ 022d6a97-2fa5-4b1d-bcac-fcb593e27cf4
qubit_plot(sol_Lind_d; record=false)

# ╔═╡ e55dcdfc-5d06-4b2a-a6f4-6271ef95a8df


# ╔═╡ c3e0da3a-9b8a-4943-bfcb-255225b8d4aa
md"""
# Two-qubit simulations
"""

# ╔═╡ 74777352-0e44-4e54-af24-97f9ed7cbaf3
begin
	
	# fluctuator Hamiltonian
	Hf1(τi, tf) = let si = map(τ -> sgn(τ, tf), τi)
					t -> map(s -> (Ω/2) * σz1 * s(t), si) end

	Hf2(τi, tf) = let si = map(τ -> sgn(τ, tf), τi)
					t -> map(s -> (Ω/2) * σz2 * s(t), si) end
	
	# total Hamiltonian --
	# -- a function of a particular realization of fluctuator Hamiltonian Hf
	H2(Hf1, Hf2) = t -> sum(Hf1(t)) + sum(Hf2(t))
	
end

# ╔═╡ b36c7a47-bf52-473b-b7b6-0ed3fefaff82
begin
	H2q = (ΩR/2) * (σy1 - σy2) #(σy1 + σy2) # qubit Hamiltonian, used for decoupling
	H2d(Hf1, Hf2) = t -> H2q + sum(Hf1(t)) + sum(Hf2(t))
end

# ╔═╡ b40c91cf-c72b-46af-8708-cd16c1d477f8
md"""
## Single trajectory
"""

# ╔═╡ 46cf17ce-10f2-4d85-b53f-ce0f4452b85f
md"""
### No decoupling
"""

# ╔═╡ 01deed2c-617d-4fbb-8e77-30f471d20f91
let
	# parameters -----------------------------------------------------------------
	ψ0 = Ψp
	dt = 1e-3  # integration time-step
	
	# Kraus operators --------------------------------------------------------------
	Hi1 = Hf1(τi, tf) # new realization of fluctuator Hamiltonian (single qubit)
	Hi2 = Hf1(τi, tf) # new realization of fluctuator Hamiltonian (single qubit)
	
	global sol2q = bayesian((0, tf), ψ0, H2(Hi1, Hi2), [], []; dt=dt)

	bell_plot(sol2q)
end

# ╔═╡ 61db37a5-f942-4342-bb90-da941c0be63c
md"""
### Decoupling
"""

# ╔═╡ 4ee9e5f7-b6ed-4869-a89a-afa1293a0851
let
	# parameters -----------------------------------------------------------------
	ψ0 = Ψp
	
	# Kraus operators --------------------------------------------------------------
	Hi1 = Hf1(τi, tf) # new realization of fluctuator Hamiltonian (single qubit)
	Hi2 = Hf1(τi, tf) # new realization of fluctuator Hamiltonian (single qubit)
	
	global sol2q_d = bayesian((0, tf), ψ0, H2d(Hi1, Hi2), [], []; dt=dt)

	# exps = map(op -> expectations(sol2q_d, dm(op)), [Ψp, Ψm])
	# plot(sol2q.t, exps[1])
	# plot!(sol2q.t, exps[2])
	
end

# ╔═╡ 085c2e52-c099-4e94-a7d9-32d574a8c9f0
bell_plot(sol2q_d)

# ╔═╡ b182f5c4-0fc6-466a-a918-7e3efd3b7ee3
md"""
## Lindblad dephasing
"""

# ╔═╡ d4a7bb23-af85-4327-9d71-b6098d54c448
let	
	ψ0 = Ψp
	J = [(σz1, γϕ), (σz2, γϕ)]

	global sol2q_lind = bayesian((0, 10.0), ψ0, Iq ⊗ Iq, J, []; dt=dt)
	global sol2q_lind_d = bayesian((0, 10.0), ψ0, H2q, J, []; dt=dt)
	
end

# ╔═╡ b0830267-af47-43a3-8e80-6dd790adc70e
bell_plot(sol2q_lind_d, sol2q_lind)

# ╔═╡ 342f3675-a184-43d1-b82a-e5798a1b2c49
md"""
## Ensemble
"""

# ╔═╡ 0224c645-8e50-4a37-acbb-22edf6ad7b58
md"""
Import from pre-simulated data:
"""

# ╔═╡ 0726ef44-9bae-417a-bd01-2d10c8129658
begin
	# load purity
	df_dd = DataFrame(CSV.File(string(pwd(), "/notebooks/data/two-qubit-fluctuators/2022-04-24T17-07-17.918/purity-trajectory.csv")))
	df_no_dd = DataFrame(CSV.File(string(pwd(), "/notebooks/data/two-qubit-fluctuators/2022-04-24T17-18-40.653/purity-trajectory.csv")))
	
	pur2q_dd = df_dd[!, 1]
	pur2q_no_dd = df_no_dd[!, 1]
end

# ╔═╡ 3d876798-2010-4227-bf49-8106e834249b
md"""
# Dephasing + measurement
"""

# ╔═╡ 16651e64-bb14-4cf2-958c-d3ad6de758c0
md"""
## Fluctuators + measurement
"""

# ╔═╡ d468b1e4-2cf0-42b7-8c9f-5c9c66108ce4
md"""
Γ = 0.1 $(@bind Γ Slider(2π.*(0.1:0.1:10.0))) 10.0
"""

# ╔═╡ 4da976ba-0424-4b73-a9d9-f9444b47157e
let
	# parameters -----------------------------------------------------------------
	ψ0 = Ψp
	
	# Kraus operators --------------------------------------------------------------
	Hi1 = Hf1(τi, tf) # new realization of fluctuator Hamiltonian (single qubit)
	Hi2 = Hf1(τi, tf) # new realization of fluctuator Hamiltonian (single qubit)

	C = [(σx1 * σx2, Γ, 1.0)]
	
	global sol2q_m = bayesian((0, tf), ψ0, H2d(Hi1, Hi2), [], C; dt=dt)
	
end

# ╔═╡ 27b96029-01e6-4323-94d8-709b1e142c57
n

# ╔═╡ 52209773-c5a8-45ff-b12b-228bf1b6856c
bell_plot(sol2q_m)

# ╔═╡ 91460147-dad7-4c5b-8b2b-c4977559d71e
md"""
## Lindblad + measurement
"""

# ╔═╡ 01978b9f-3325-479f-859f-2e77b2e5f730
Γ/(2π)

# ╔═╡ 412556d8-5072-4b36-acdc-ca80b37585b1
let	
	ψ0 = Ψp
	J = [(σz1, γϕ), (σz2, γϕ)]
	C = [(σx1 * σx2, Γ, 1.0)]

	global sol2q_lind_m = bayesian((0, 2.0), ψ0, Iq ⊗ Iq, J, C; dt=dt)
	
end

# ╔═╡ 22ed5bad-200b-4863-aca8-9b26756d44c9
bell_plot(sol2q_lind_m)

# ╔═╡ be327ea2-1f14-4f43-ba06-a278b08c1988
md"""
# Exact Lindblad solutions
"""

# ╔═╡ dfbc7eb2-ae1c-480e-90e2-37486b256c7c
begin
	# parameters
	γ = γϕ
	h = ΩR
end

# ╔═╡ b11251ea-9f0f-422d-8303-e875829adc3f
function bloch_to_w(v, h)
	u = γ / h
	
	w⃗₀ = [1, 0, 0]
	w⃗₊ = [0, 1, 0]
	w⃗₋ = [0, 0, 1]
	
	x̂ = w⃗₀
	ŷ = 1/(2√Complex(u^2 - 4)) .* (w⃗₊ .- w⃗₋)
	ẑ = (w⃗₊ .+ w⃗₋) / 2 - u * ŷ 

	x, y, z = v
	return x .* x̂ .+ y .* ŷ .+ z .* ẑ
end

# ╔═╡ 38c911b9-6371-4d63-960a-2fd2ed46c175
begin

	bloch0 = [0, 1, 0]

	# v(t, v0) = (v0[1] * exp(h * z0 * t) .* w0) + (v0[2] * exp(h * zp * t) .* wp) + (v0[3] * exp(h * zm * t) .* wm) 
	vw(t, h) = let
		u = γ / h

		w0, wp, wm = bloch_to_w(bloch0, h)
		
		# eigenvalues
		z0 = -2u
		zp = -u + √Complex(u^2 - 4)
		zm = -u - √Complex(u^2 - 4)
		
		return [w0 * exp(h * z0 * t), wp * exp(h * zp * t), wm * exp(h * zm * t)] 
	end
end

# ╔═╡ 89793870-57af-432e-85b2-ea1354d701de
function w_to_bloch(v, h)
	u = γ / h
	
	x̂ = [1, 0, 0]
	ŷ = [0, 1, 0]
	ẑ = [0, 0, 1]

	w⃗₀ = x̂
	w⃗₊ = ((-u + √Complex(u^2 - 4)) .* ŷ) + 2 .* ẑ
	w⃗₋ = ((-u - √Complex(u^2 - 4)) .* ŷ) + 2 .* ẑ

	w0, wp, wm = v
	return w0 .* w⃗₀ .+ wp .* w⃗₊ .+ wm .* w⃗₋
end

# ╔═╡ 62416024-3982-4e01-8caf-1253d7ac82e9
begin
	ts = 0:0.001:10.0
	vs = [map(t -> vw(t, h)[i], ts) for i in 1:3]
	labels = [L"a_0(t)", L"a_+(t)",  L"a_-(t)"]
	ps = []
	for i in 1:3
		push!(ps, plot(ts, real(vs[i]), color=colors1[i], label=string("Re(", labels[i], ")"), linewidth=2))
		push!(ps, plot(ts, imag(vs[i]), color=colors1[i], linestyle=:dash, label=string("Im(", labels[i], ")"), linewidth=2))
	end
end

# ╔═╡ c5b49272-8da3-4bdc-86c0-d3a2e306f516
length(ts)

# ╔═╡ 5f883c57-b754-4c44-a523-0fce0ffd262e
plot(ps..., layout=(3,2))

# ╔═╡ e5f6c8ae-3b6e-4d27-a3e1-0788e92bf914
vbloch(t, h) = real(w_to_bloch(vw(t, h), h))

# ╔═╡ 9dc47aa2-b545-472f-82a2-a4c033eb9794
let
	vs_bloch = [map(t -> vbloch(t, h)[i], ts) for i in 1:3]
	labels = [L"x(t)", L"y(t)",  L"z(t)"]
	global ps_bloch = []
	for i in 1:3
		push!(ps_bloch, plot(ts, vs_bloch[i], color=colors1[i], label=labels[i], linewidth=0.5))
	end
end

# ╔═╡ e9c32550-5983-4d44-ad9a-c41d53ee2938
plot(ps_bloch..., layout=(3,1))

# ╔═╡ fd9a9252-4d55-42bf-80e4-2d511b12dd19
function density_matrix(vbloch)
	x, y, z = vbloch
	return DenseOperator(0.5 * (Iq + x * σx + y * σy + z * σz))
end

# ╔═╡ eb6d5e68-43dc-40cc-aa6b-6fb70bc979bc
ρ0 = DenseOperator(density_matrix(bloch0))

# ╔═╡ 83956a26-8730-4867-ac3f-2cfc6e5cfd63
ρ(t, h) = density_matrix(vbloch(t, h))

# ╔═╡ f435790d-7fcc-4c56-acdb-35ca9297cbae
purity(t, h) = real(tr(ρ(t, h) * ρ(t, h)))

# ╔═╡ b027493e-f9a2-4c48-91e5-64759b456ab2
let
	p = plot(title="purity", xlabel="t (μs)")
	plot!(ts, purity.(ts, ΩR), label=string("ΩR = ", round(ΩR/(2π), digits=3), " MHz"))
	plot!(ts, purity.(ts, 0.01), label=string("ΩR = 0.001 MHz"))
	
end

# ╔═╡ 3f873dd4-166d-43f6-8e9b-79f2c676be53
md" γfit1 = 0.001 MHz $(@bind γfit1 Slider((2π) * 0.005:0.005:1.0)) 1.0 MHz"

# ╔═╡ ac09d253-e5be-4274-a6fd-99fa86f257bc
γfit1/(2π)

# ╔═╡ ddc5e3b1-49e3-46b6-9ef9-f8034dab4c32
md" γfit2 = 0.005 MHz $(@bind γfit2 Slider((2π) .* 0.005:0.005:1.0)) 1.0 MHz"

# ╔═╡ ce87ea1a-a9cf-4a6c-bf56-0ca904b901f1
γfit2 / (2π)

# ╔═╡ bd74f7dc-0ecb-433e-af87-474abf551c7e
ρexp(t, h) = let U = exp(-im * h * DenseOperator(σx) * t)
				return U * ρ0 * U'
			end

# ╔═╡ ef205170-c760-4cca-91d6-7be42d924785
fid(t, h) = fidelity(ρexp(t, h), ρ(t, h))

# ╔═╡ 70858fde-3292-4164-8537-8c56b4094755
let
	p = plot(title=string("fidelity with state ", L"U(t) \rho U^\dagger(t)"), xlabel="t (μs)")
	plot!(ts, fid.(ts, ΩR), linewidth=1.5, label=string("ΩR = ", round(ΩR/(2π), digits=3), " MHz"))
	plot!(ts, fid.(ts, 0.001 * (2π)), linewidth=1.5, label=string("ΩR = 0.0 MHz"))
	
end

# ╔═╡ ef32d7fb-5132-4fcc-abc5-40c580b37a80
fidf = last(fid.(ts, ΩR))

# ╔═╡ cfa0bd83-0aad-48bf-b7e1-abc4decebba4
let
	p = plot(title=string("fidelity with state ", L"U(t) \rho U^\dagger(t)"), xlabel="t (μs)")
	plot!(ts, fid.(ts, ΩR), linewidth=1.5, label=string("ΩR = ", round(ΩR/(2π), digits=3), " MHz"))
	plot!(ts, fid.(ts, 0.001 * (2π)), linewidth=1.5, label=string("ΩR = 0.0 MHz"))
	plot!(ts, (1 - fidf) .* exp.(-2γfit1 .* ts) .+ fidf, linewidth=1.5, label=string("ΩR = 25.0 MHz fit"))
	plot!(ts, (1 - fidf) .* exp.(-2 * 0.0965*2π .* ts) .+ fidf, linewidth=1.5, label=string("ΩR = 0.0 MHz fit"))
	
end

# ╔═╡ 468b182a-4ea5-4a5d-bc65-800d6549ec29
colors1

# ╔═╡ c6eb30e7-1975-4772-a1be-6bec1120b568
v(0.5, v0)[1]

# ╔═╡ 91128165-d015-4265-b117-afaf7cf3caf0
md"""
# Markovianity
"""

# ╔═╡ 4c342d55-83eb-4c15-8048-ed5f903317e7
md"""
Rotate view: ϕv = $(@bind ϕv Slider(2π .* 0:0.02:1))
"""

# ╔═╡ cab84d3b-89a2-44aa-b03d-e72c8272b0a7
# ψ0s = map(p -> 0.5 * (Iq + p[1] * σx + p[2] * σy + p[3] * σz), points)

# ╔═╡ 06fd35c8-bb6d-491b-86b5-a67bc16ba93c
times

# ╔═╡ dc8a4138-36a6-424e-86b1-5f18dea499f6
begin
	tm = g ⊗ g
	tp = e ⊗ e
	s = Ψm
	t0 = Ψp
end

# ╔═╡ 9220ef01-a001-4435-9ff2-907a896a7c43
dagger(t0)

# ╔═╡ 1a2879db-381a-40f5-add4-73242f776e58
projector(tm, dagger(t0))

# ╔═╡ a0610e48-fce4-4a66-bd0d-caa155435018
 tm * dagger(t0)   

# ╔═╡ 395cb4c1-22f2-444f-9a8c-5b77c1241960


# ╔═╡ ebb259a8-3234-4286-8218-94c83be5545a
function positive(deriv)
	for val in deriv
		if val > 0.01
			return true
		end
	end
	return false
end

# ╔═╡ 0f5189fb-f05c-4aa5-98ae-218b13acb4b4
md"""
# Utilities
"""

# ╔═╡ f8eb195f-e0d9-48ef-bf22-e1af09cd47a1
xyz(θ, ϕ) = (sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ))

# ╔═╡ d517818a-9f27-4b36-b473-8866a188e7f4
begin
	ϕs = [0, π/2, π, 3π/2]
	θs = [0, π/4, π/2, 3π/4, π]
	points = []
	ψ0s = []
	for ϕ in ϕs
		for θ in θs
			push!(points, xyz(θ,ϕ))
			ψ = cos(θ/2) * g + sin(θ/2) * exp(-im * ϕ) * e
			push!(ψ0s, ψ)
		end
	end
end

# ╔═╡ efd0e637-d17f-431c-bafc-8bde3f25325e
length(ϕs) * length(θs)

# ╔═╡ 77a6f475-537f-45dd-a41c-ed3a55fd58a4
ρavgs = map(ψ0s) do ψ0

	solutions = map(1:2000) do i
		Hi = Hf(τi, tf) # new realization of fluctuator Hamiltonian
		bayesian((0, tf), ψ0, H(Hi), [], []; dt=dt)
	end
	
	ρavg = 	if typeof(solutions[1].ρ[1]) <: Ket
					[mean(map(sol -> dm(sol.ρ[i]), solutions)) for i in 1:length(solutions[1].t)]
			else
					[mean(map(sol -> sol.ρ[i], solutions)) for i in 1:length(solutions[1].t)]
			end

end	

# ╔═╡ 775ec8f0-13cd-4bdd-af39-119d105b91a9
begin
	tracedistances = []
	derivs = []
	pairs = []
	for i in 1:length(ρavgs)
		ρ1 = ρavgs[i]
		for j in (i + 1):length(ρavgs)
			ρ2 = ρavgs[j]
			td = map(i -> tracedistance(ρ1[i], ρ2[i]) , 1:length(ρ1))
			td_deriv = map(i -> (td[i+1] - td[i])/dt, 1:(length(ρ1)-1))
			push!(tracedistances, td)
			push!(derivs, td_deriv)
			push!(pairs, (ρ1, ρ2))
		end
	end
end

# ╔═╡ 776fdc9a-7ab9-4594-9b30-c203305a26ca
length(tracedistances)

# ╔═╡ b4eb197e-2897-46c9-a059-6fd418f21f06
begin
	p = plot(xlabel="t (μs)", ylabel="trace distance", labelfontsize=10, tickfontsize=10)
	for td in tracedistances
		plot!(times, td, legend=:none)
	end
	p
end

# ╔═╡ eed3f4d0-fe75-4efc-bb9e-6fe5c353c6e3
begin
	pp = plot(xlabel="t (μs)", ylabel="derivative", labelfontsize=10, tickfontsize=10)
	for d in derivs
		plot!(times[2:end], d, legend=:none, linewidth=1, alpha=0.5)
	end
	pp
end

# ╔═╡ 0741932f-6653-4a32-8bc5-8b4b4716cc05
plot(p, pp, layout=grid(2,1))

# ╔═╡ e6f40df0-128e-485a-9689-ee7f6ee0a774
length(derivs)

# ╔═╡ f4747c1a-c906-4a3e-8af4-493492f2eada
begin
	pl = plot()
	for deriv in derivs[1:2]
		plot!(times[2:end], deriv, legend=:none)
	end
	pl
end

# ╔═╡ a1f3c584-7d42-43e3-9e88-3af000ca9596
derivs[1][1] > 0

# ╔═╡ deec843a-010f-46b4-8444-e6f1ee364618
begin
	maxima = map(d -> maximum(d), derivs)
	m = maximum(maxima)
	index = findall(x->x==m, maxima)
end

# ╔═╡ 73157668-1e8d-41f3-beae-263142e04883
begin
	(ρ1, ρ2) = pairs[index[1]]
	exps1 = map(op -> map(ρ -> real(expect(ρ, op)), ρ1), [σx, σy, σz])
	exps2 = map(op -> map(ρ -> real(expect(ρ, op)), ρ2), [σx, σy, σz])

	qubit_plot((times, exps1, []), (times, exps2, []); ylims=[-1,1])
end

# ╔═╡ 9ac6937c-245e-4ee5-a92b-84ae6432cc82
function purities(sol::Solution)
	ρs = (typeof(sol.ρ[1]) <: Ket) ? dm.(sol.ρ) : sol.ρ
	map(ρ -> real(tr(ρ * ρ)), ρs)
end

# ╔═╡ 0a6f05f6-6bf3-4ca9-8543-3633b6819e87
purities(exps) = 0.5 .* (1 .+ exps[1].^2 .+ exps[2].^2 .+ exps[3].^2)

# ╔═╡ fd61db11-5347-48bd-82b6-6aabc61f7959
purities(sols_d[1])

# ╔═╡ b264ef5d-eff2-406e-a4e2-0ade25d51112
let
	ts = sol2q_lind_d.t[1:2000]
	
	plot(ts, pur2q_no_dd[1:2000], label=string("fluctuators, ", L"\Omega_R = ", 0.0, " MHz"), color=:blue, linestyle=:dash, tickfontsize=:12, legendfontsize=10, title="two qubits", labelfontsize=12, linewidth=2)
		plot!(ts, pur2q_dd[1:2000], label=string("fluctuators, ", L"\Omega_R = ", (ΩR/2π), " MHz"), legend=:right, xlabel="t (μs)", ylabel="purity", color=:blue, linewidth=2)
		
		plot!(ts, purities(sol2q_lind)[1:2000], label=string("Lindblad sim, ", L"\Omega_R = ", 0.0, " MHz"), color=:black, linewidth=2, linestyle=:dash)
		plot!(ts, purities(sol2q_lind_d)[1:2000], label=string("Lindblad sim, ", L"\Omega_R = ", (ΩR/2π), " MHz"), legend=:right, xlabel="t (μs)", ylabel="purity", color=:black, linewidth=2)
end

# ╔═╡ 1fc7a91b-8303-4a24-8cd9-3430033cf316
function ensembleavg(sols::Vector{Solution}; ops = qbasis)
	exparrs = map(op -> [], ops)

	for sol in sols
		exps = map(op -> expectations(sol, op), ops)
		for (list, traj) in zip(exparrs, exps)
			push!(list, traj)
		end
	end

	return map(arr -> mean(arr), exparrs)

end

# ╔═╡ 44cf9a24-57f4-414f-902a-0c4fdac9516d
exps = ensembleavg(sols)

# ╔═╡ 998a6f49-1199-4279-92d4-5babd8978948
 qubit_plot((times, exps, []); title="", legendpos=:bottomleft)

# ╔═╡ 8c721925-1604-4cc2-a28f-f3c08a4ac2df
let
	x, y, z = exps
	t = collect(times)
	if anim_plot
		anim = @animate for i ∈ range(1, length(times), step=25)
				blochsphere(i, x, y, z, linewidth=1., linealpha=0.85, ax=true, viewϕ=0, blochmark=true) end
		gif(anim, fps = 10)
	else
		blochsphere(x, y, z, linewidth=1., linealpha=0.85, ax=true, viewϕ=0, blochmark=true)
	end
end	

# ╔═╡ e9054568-0f76-46c8-aaa1-3a3b014a432e
begin
	tt = collect(times)
	qubit_plot((tt, exps_lind, []), (tt, exps, []); ylims=[0,1], l1="lindblad", l2="fluctuators")
end

# ╔═╡ 92803ff8-1148-4061-92b5-24c7230ccd78
exps_d = ensembleavg(sols_d)

# ╔═╡ 44c75dab-6383-4070-b365-da8d8492c6f5
let
	x, y, z = exps_d
	t = collect(times)
	if anim_plot_2
		anim = @animate for i ∈ range(1, length(times), step=25)
				blochsphere(i, x, y, z, linewidth=1., linealpha=0.85, ax=true, viewϕ=0, blochmark=true) end
		gif(anim, fps = 10)
	else
		blochsphere(x, y, z, linewidth=1., linealpha=0.85, ax=true, viewϕ=0, blochmark=true)
	end
end	

# ╔═╡ 662a593a-6ce9-407b-b1db-e06ab085df82
let
	p = purities(exps)
	pd = purities(exps_d)
	plot(times, p, label="no decoupling", color=:blue, linestyle=:dash, linewidth=2)
	plot!(times, pd, label=string("decoupling at ΩR = ", (ΩR/2π), " MHz"), legend=:right, xlabel="t (μs)", ylabel="purity", color=:blue, linewidth=2)
end

# ╔═╡ 678aa2d2-1aff-467c-9219-ed1b43f1734c
let
	p = purities(exps)
	pd = purities(exps_d)
	plind = purities(exps_lind)
	plindd = purities(exps_lind_d)
	plindd2 = purities(exps_lind_d2)
	
	plot(times, p, label="no decoupling", color=:blue, linestyle=:dash)
	plot!(times, pd, label=string("decoupling at ΩR = ", (ΩR/2π), " MHz"), legend=:right, xlabel="t (μs)", ylabel="purity", color=:blue)
	
	plot!(times, plind, label="no decoupling (Lindblad)", color=:red, linestyle=:dash)
	plot!(times, plindd, label=string("decoupling (Lindblad)"), legend=:right, xlabel="t (μs)", ylabel="purity", color=:red)
	plot!(times2, plindd2, label=string("decoupling (Lindblad, dt=dt/100)"), legend=:right, xlabel="t (μs)", ylabel="purity", color=:red, linestyle=:dot)
end

# ╔═╡ 7b9a6f80-e14f-4ee9-95bb-d2aa0641343f
let
	p = purities(exps)
	pd = purities(exps_d)
	plind = purities(exps_lind)
	plindd = purities(exps_lind_d)
	plindd2 = purities(exps_lind_d2)
	
	plot(times, p, label=string("fluctuators, ", L"\Omega_R = ", 0.0, " MHz"), color=:blue, linestyle=:dash, tickfontsize=:12, legendfontsize=10, title="single qubit", labelfontsize=12)
	plot!(times, pd, label=string("fluctuators, ", L"\Omega_R = ", (ΩR/2π), " MHz"), legend=:right, xlabel="t (μs)", ylabel="purity", color=:blue)
	
	plot!(times, plind, label=string("Lindblad sim., ", L"\Omega_R = ", 0.0, " MHz"), color=:black, linestyle=:dash, linewidth=2)
	plot!(times, plindd, label=string("Lindblad sim., ", L"\Omega_R = ", (ΩR/2π), " MHz"), legend=:right, xlabel="t (μs)", ylabel="purity", color=:black, linewidth=2)

	plot!(times, purity.(times, 0.01), color=:red, linewidth=2, linestyle=:dot, label=string("Lindblad anal., ", L"\Omega_R = ", 0.0, " MHz"))
	plot!(times, purity.(times, ΩR), color=:red, linestyle=:dash, linewidth=2, label=string("Lindblad anal., ", L"\Omega_R = ", (ΩR/2π), " MHz"))
	
end

# ╔═╡ 1b7e3334-9825-4832-a0b8-d58816fd5236
begin
	@userplot PlotFluctuators
	@recipe function f(pf::PlotFluctuators)
		ωi, times, Hi, Ω = pf.args
	
		labels = permutedims(map(ω -> string(ω, " MHz"), ωi))
	
		# Plot time series --------------------------------------------------------
	
		xlabel --> "t (μs)"
		ylabel --> "fluctuator value"
		legend --> :outerright
		label --> labels
	
		palette := :tab10
		linealpha --> 1
	
		legendfontsize --> 6
		titlefontsize --> 12
		xtickfontsize --> 8
		ytickfontsize --> 8
		xguidefontsize --> 10
		yguidefontsize --> 10
		size --> (600,300)
		linewidth --> 1.5
	
		for i in reverse(1:length(ωi))
	
			fluctuators = map(t -> real(expect(Hi(t)[i], σz)) / Ω, times)
	
			@series begin
				times, fluctuators
			end
	
		end
	end
end

# ╔═╡ 2e072a88-7166-4d5d-96e6-52fcbe668505
function bloch_plots(sols::Vector{Solution}; alpha=0.1, N=50, kwargs...)
	colors = palette(:tab10) 
	
	# calculate expectation values --------------------------------------------
	t = sols[1].t
	xs, ys, zs = [], [], []
	
	for sol in sols
		x, y, z = map(op -> expectations(sol, op), qbasis)
		for (list, traj) in zip([xs, ys, zs], [x, y, z])
			push!(list, traj)
		end
		
	end
	
	
	# plot ----------------------------------------------------------------------
	function bloch(os; color=colors1q[1], xlabel="", ylabel="")
		
		po = plot(ylims = [-1,1], xlabel=xlabel, ylabel=ylabel)
		
		for o in os[1:min(N, 50)]
			plot!(t, o, alpha=alpha, label=:none, color=color)
		end
		
		oavg = [mean([os[i][j] for i in 1:N]) for j in 1:length(t)]
		plot!(t, oavg, alpha=1, color=color, linewidth=3)

		po
		
	end

	
	l = @layout [xplot{0.33h}; yplot{0.33h}; zplot{0.33h}]
	
	px = bloch(xs, color=colors[1], ylabel="x")
	py = bloch(ys, color=colors[2], ylabel="y")
	pz = bloch(zs, color=colors[3], ylabel="z", xlabel="t (μs)")
	
	plot(px, py, pz, layout = l, link=:y, size=(800,500), legendfontsize=8, titlefontsize=12, legend=:none, kwargs...)
	
end

# ╔═╡ 6e506109-b181-4be2-a243-9fc4d1d4af17
 bloch_plots(sols)

# ╔═╡ e1eb9c92-d3ef-4549-b1f1-e65c0bf77ee7
 bloch_plots(sols_d)

# ╔═╡ 06ace002-be3d-4478-860d-7dc23838c9de
function bloch_vectors(vecs...; mesh=30, ax=false, viewϕ=0, connecting=false, labels=[], size=(400,400))
	bcolors = palette(:seaborn_bright)

	# Wire frame coordinates ---------------------------------------------------

	x(θ, ϕ) = sin(θ) * cos(ϕ + viewϕ)
	y(θ, ϕ) = sin(θ) * sin(ϕ + viewϕ)
	z(θ, ϕ) = cos(θ)

	θs = range(0, 2π, length=mesh)
	ϕs = range(0, π, length=mesh) .+ viewϕ


	# Plot wireframe -----------------------------------------------------------
	wf = plot()
	
	# Longitudes
	for ϕ in ϕs
		plot!([x(θ, ϕ) for θ in θs], [y(θ, ϕ) for θ in θs], [z(θ, ϕ) for θ in θs], label=:none, linecolor="steelblue", linewidth=0.5, linealpha=1, seriestype=path3d, aspect_ratio=1.0, size=size)
			
	end

	# Latitudes
	for θ in θs
		plot!([x(θ, ϕ) for ϕ in ϕs], [y(θ, ϕ) for ϕ in ϕs], [z(θ, ϕ) for ϕ in ϕs], label=:none, linecolor="steelblue", linewidth=0.5, linealpha=1, seriestype=path3d, aspect_ratio=1.0, size=size)
	end


	# Plot reference axes ------------------------------------------------------
	
	colors = [palette(:tab10)[i] for i in 1:3]

	if ax
		plot!([0, cos(viewϕ)], [0, sin(viewϕ)], [0, 0], linecolor=colors[1], linewidth=3.0, label=:none)
		plot!([0, -sin(viewϕ)], [0, cos(viewϕ)], [0, 0], linecolor=colors[2], linewidth=3.0, label=:none)
		plot!([0, 0], [0, 0], [0, 1], linecolor=colors[3], linewidth=3.0, label=:none)
	end


	# Plot Bloch vectors input by user --------------------------------

	for (i, vec) in enumerate(vecs)

		(xvv, yvv, zvv) = vec

		xv = xvv * cos(viewϕ) - yvv * sin(viewϕ)
		yv = xvv .* sin(viewϕ) .+ yvv * cos(viewϕ)
		zv = zvv

		if connecting
			plot!([0, xv], [0, yv], [0, zv], label=:none, linewidth=2, linecolor=bcolors[i])
		end
	
		plot!([xv], [yv], [zv], legend=:outerright, legendfontsize=10, marker=(:circle, 5), markercolor=bcolors[mod(i, length(bcolors)) + 1], label=try labels[i] catch e "" end)
	end

	return wf


end

# ╔═╡ bad8f75f-a350-41cd-a960-38c26f1a14f0
bloch_vectors(points..., ax=true, size=(600,600), viewϕ = ϕv)

# ╔═╡ bb58f61d-334f-424e-80c1-4730abb59050
mod(10, 10)

# ╔═╡ Cell order:
# ╟─e625f413-e93b-48c9-b14d-181c37b0fe0a
# ╟─249f33ea-92aa-4b9c-bdfa-d35fe5254225
# ╟─3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─7cb4592d-4b47-41c6-9cb3-cca1033b6c9c
# ╟─81394981-ea8f-4718-b4b3-138546ca8a71
# ╠═41954632-a8c7-49d1-9c7e-559c0b0c0556
# ╟─9d1860c9-0b12-4d7b-97d6-0ea51db03a94
# ╟─23bc6da9-15fc-4ee4-93a9-1abf6dd207d3
# ╟─9771ac72-8daf-4323-8267-b91410b150ff
# ╠═176571f5-291c-40dd-964b-7aa59f1ed812
# ╠═3f4a410c-a6f8-45e4-b7d0-e10d583fa91f
# ╠═ec4a876c-c4b8-4f4b-9c60-9827b8ebe337
# ╟─181fd41d-e256-4c78-a100-899acd89a066
# ╠═7bbba21f-8d77-4d2e-8a7b-ddeb6c291427
# ╟─65a619fc-ff14-4426-9c45-6d749fd8b78c
# ╠═6d939bd2-2689-449d-97ca-34f7bf838858
# ╠═689a12f4-41cd-4efc-a88c-92f2fff864f8
# ╟─ada56af9-023f-4364-ab67-2f4de3e15e0c
# ╟─567ff65b-c70d-42fc-9f31-7373f543bd9e
# ╟─846b569e-1440-4cd1-8a5f-0db962cc7d3a
# ╟─cbaac8ff-3898-4413-96c9-2d07d46ca5b5
# ╠═dba94b12-aac9-4cf2-8714-6cbb75e473e9
# ╟─44cf9a24-57f4-414f-902a-0c4fdac9516d
# ╟─6c0d59cd-0f9d-463b-aa97-c1617ab0aecf
# ╟─a10ad7ca-873c-41cf-a4bc-edd1e4288b02
# ╟─210332f6-d26b-4f5b-85f9-20ac9555a97c
# ╠═6e506109-b181-4be2-a243-9fc4d1d4af17
# ╠═998a6f49-1199-4279-92d4-5babd8978948
# ╟─fd4b2d5c-0ed3-4cf3-bf8d-dbc1543d0a71
# ╟─8c721925-1604-4cc2-a28f-f3c08a4ac2df
# ╟─da766962-481b-4db5-940f-662225147988
# ╟─1aa49c31-f055-4fee-abd3-6af4e3593d12
# ╟─547046a1-6163-41e2-b397-d41253ab0b8f
# ╟─490899df-edcc-4a7d-87fa-e4d0b3b103d2
# ╠═fd61db11-5347-48bd-82b6-6aabc61f7959
# ╠═92803ff8-1148-4061-92b5-24c7230ccd78
# ╠═e1eb9c92-d3ef-4549-b1f1-e65c0bf77ee7
# ╟─7217689a-ecda-4714-95a9-400d85f11282
# ╟─44c75dab-6383-4070-b365-da8d8492c6f5
# ╟─662a593a-6ce9-407b-b1db-e06ab085df82
# ╠═6e6bb8e1-892e-4fd4-8e48-20988c6eecba
# ╟─51ae5fda-9444-42e4-b190-a52215c39c87
# ╟─48c21723-679d-40eb-a226-2988fa448f61
# ╟─7175e140-44f8-4961-9e83-367c20d00f26
# ╠═9858f17b-7e3a-4b3c-83e0-52a45d8869db
# ╠═aaf39f5e-aeb7-4b2b-a2e7-be950c723995
# ╠═c3bd9db7-a958-47f9-98b6-4f3bb1bddbc7
# ╠═7232c30f-ad84-4a17-bf04-d8191244436c
# ╠═e9054568-0f76-46c8-aaa1-3a3b014a432e
# ╟─ed3be650-e080-4928-b3d2-24445da8b9d4
# ╟─b781fb8f-d7a0-4502-8322-e3b8b0e1ce1c
# ╟─6d5f64a4-84cd-4204-8169-a1103b241ae0
# ╟─b0f750db-6c58-4939-9f42-8388ba8b9af4
# ╟─b6a49a95-8e21-4e81-8f39-737c7a628a0e
# ╠═022d6a97-2fa5-4b1d-bcac-fcb593e27cf4
# ╟─678aa2d2-1aff-467c-9219-ed1b43f1734c
# ╠═e55dcdfc-5d06-4b2a-a6f4-6271ef95a8df
# ╠═7b9a6f80-e14f-4ee9-95bb-d2aa0641343f
# ╟─c3e0da3a-9b8a-4943-bfcb-255225b8d4aa
# ╟─74777352-0e44-4e54-af24-97f9ed7cbaf3
# ╠═b36c7a47-bf52-473b-b7b6-0ed3fefaff82
# ╟─b40c91cf-c72b-46af-8708-cd16c1d477f8
# ╟─46cf17ce-10f2-4d85-b53f-ce0f4452b85f
# ╠═01deed2c-617d-4fbb-8e77-30f471d20f91
# ╟─61db37a5-f942-4342-bb90-da941c0be63c
# ╠═4ee9e5f7-b6ed-4869-a89a-afa1293a0851
# ╟─085c2e52-c099-4e94-a7d9-32d574a8c9f0
# ╟─b182f5c4-0fc6-466a-a918-7e3efd3b7ee3
# ╠═d4a7bb23-af85-4327-9d71-b6098d54c448
# ╟─b0830267-af47-43a3-8e80-6dd790adc70e
# ╟─342f3675-a184-43d1-b82a-e5798a1b2c49
# ╟─0224c645-8e50-4a37-acbb-22edf6ad7b58
# ╟─0726ef44-9bae-417a-bd01-2d10c8129658
# ╠═b264ef5d-eff2-406e-a4e2-0ade25d51112
# ╠═c5b49272-8da3-4bdc-86c0-d3a2e306f516
# ╠═3d876798-2010-4227-bf49-8106e834249b
# ╟─16651e64-bb14-4cf2-958c-d3ad6de758c0
# ╠═d468b1e4-2cf0-42b7-8c9f-5c9c66108ce4
# ╠═4da976ba-0424-4b73-a9d9-f9444b47157e
# ╠═27b96029-01e6-4323-94d8-709b1e142c57
# ╠═52209773-c5a8-45ff-b12b-228bf1b6856c
# ╟─91460147-dad7-4c5b-8b2b-c4977559d71e
# ╠═01978b9f-3325-479f-859f-2e77b2e5f730
# ╠═412556d8-5072-4b36-acdc-ca80b37585b1
# ╠═22ed5bad-200b-4863-aca8-9b26756d44c9
# ╟─be327ea2-1f14-4f43-ba06-a278b08c1988
# ╠═dfbc7eb2-ae1c-480e-90e2-37486b256c7c
# ╠═38c911b9-6371-4d63-960a-2fd2ed46c175
# ╠═b11251ea-9f0f-422d-8303-e875829adc3f
# ╠═89793870-57af-432e-85b2-ea1354d701de
# ╠═62416024-3982-4e01-8caf-1253d7ac82e9
# ╠═5f883c57-b754-4c44-a523-0fce0ffd262e
# ╠═e5f6c8ae-3b6e-4d27-a3e1-0788e92bf914
# ╠═9dc47aa2-b545-472f-82a2-a4c033eb9794
# ╠═e9c32550-5983-4d44-ad9a-c41d53ee2938
# ╠═fd9a9252-4d55-42bf-80e4-2d511b12dd19
# ╠═eb6d5e68-43dc-40cc-aa6b-6fb70bc979bc
# ╠═83956a26-8730-4867-ac3f-2cfc6e5cfd63
# ╠═ef205170-c760-4cca-91d6-7be42d924785
# ╠═f435790d-7fcc-4c56-acdb-35ca9297cbae
# ╟─b027493e-f9a2-4c48-91e5-64759b456ab2
# ╠═70858fde-3292-4164-8537-8c56b4094755
# ╠═3f873dd4-166d-43f6-8e9b-79f2c676be53
# ╠═ac09d253-e5be-4274-a6fd-99fa86f257bc
# ╠═ddc5e3b1-49e3-46b6-9ef9-f8034dab4c32
# ╠═ce87ea1a-a9cf-4a6c-bf56-0ca904b901f1
# ╠═cfa0bd83-0aad-48bf-b7e1-abc4decebba4
# ╠═ef32d7fb-5132-4fcc-abc5-40c580b37a80
# ╠═bd74f7dc-0ecb-433e-af87-474abf551c7e
# ╠═468b182a-4ea5-4a5d-bc65-800d6549ec29
# ╠═c6eb30e7-1975-4772-a1be-6bec1120b568
# ╟─91128165-d015-4265-b117-afaf7cf3caf0
# ╠═efd0e637-d17f-431c-bafc-8bde3f25325e
# ╠═d517818a-9f27-4b36-b473-8866a188e7f4
# ╠═4c342d55-83eb-4c15-8048-ed5f903317e7
# ╠═bad8f75f-a350-41cd-a960-38c26f1a14f0
# ╠═cab84d3b-89a2-44aa-b03d-e72c8272b0a7
# ╠═77a6f475-537f-45dd-a41c-ed3a55fd58a4
# ╠═775ec8f0-13cd-4bdd-af39-119d105b91a9
# ╠═776fdc9a-7ab9-4594-9b30-c203305a26ca
# ╠═06fd35c8-bb6d-491b-86b5-a67bc16ba93c
# ╠═b4eb197e-2897-46c9-a059-6fd418f21f06
# ╠═eed3f4d0-fe75-4efc-bb9e-6fe5c353c6e3
# ╠═dc8a4138-36a6-424e-86b1-5f18dea499f6
# ╠═9220ef01-a001-4435-9ff2-907a896a7c43
# ╠═1a2879db-381a-40f5-add4-73242f776e58
# ╠═a0610e48-fce4-4a66-bd0d-caa155435018
# ╠═e6f40df0-128e-485a-9689-ee7f6ee0a774
# ╠═0741932f-6653-4a32-8bc5-8b4b4716cc05
# ╠═f4747c1a-c906-4a3e-8af4-493492f2eada
# ╠═395cb4c1-22f2-444f-9a8c-5b77c1241960
# ╠═a1f3c584-7d42-43e3-9e88-3af000ca9596
# ╠═73157668-1e8d-41f3-beae-263142e04883
# ╠═deec843a-010f-46b4-8444-e6f1ee364618
# ╠═ebb259a8-3234-4286-8218-94c83be5545a
# ╟─0f5189fb-f05c-4aa5-98ae-218b13acb4b4
# ╠═f8eb195f-e0d9-48ef-bf22-e1af09cd47a1
# ╠═9ac6937c-245e-4ee5-a92b-84ae6432cc82
# ╠═0a6f05f6-6bf3-4ca9-8543-3633b6819e87
# ╠═1fc7a91b-8303-4a24-8cd9-3430033cf316
# ╠═1b7e3334-9825-4832-a0b8-d58816fd5236
# ╠═2e072a88-7166-4d5d-96e6-52fcbe668505
# ╠═06ace002-be3d-4478-860d-7dc23838c9de
# ╠═bb58f61d-334f-424e-80c1-4730abb59050
# ╠═df539f73-7c8a-4d65-b909-0d94720f724f
