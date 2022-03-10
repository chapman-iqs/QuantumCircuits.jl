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

# ╔═╡ 3d85dac8-5abb-4622-8eeb-25afb89a71a6
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
	using Plots
	using QuantumCircuits

	md" # Packages and julia files"
	
	
end

# ╔═╡ fab5ca8b-29b3-4456-8520-c00927886612
md"""
In this interactive notebook, we'll look at ensembles of reduced qubit evolution as a basic example to gain intuition for Cavity Quantum Electrodynamics (CQED).
"""

# ╔═╡ 32746f87-236e-4e09-9e19-cfc645f7296e
TableOfContents(title="Qubit ensembles")

# ╔═╡ 07821c37-73f0-4bb8-94c8-c843c23976b9
md" # Reduced qubit description"

# ╔═╡ ee061b17-e095-417d-956e-ae1fb058db65
md"""
We work in the rotating frame, so that $\omega_q = 0$ effectively. Thus the Hamiltonian is

$\hat H = \frac{\Omega_R}{2} \hat \sigma_y.$
"""

# ╔═╡ d3c789f7-6fc2-4448-85c3-c1cd51c95155
md"""
The qubit is weakly and continuously monitored in its informational basis, $|+z \rangle$ (excited) and $|-z \rangle$ (ground). This ''weak measurement'' manifests in two primary ways:

* The information *collected* from the measurement creates **measurement backaction**, characterized by measurement collapse timescale $\tau$. This is represented via the dissipation operator $\hat J = \sqrt{\Gamma} \hat \sigma_z$.

* The information *not collected* leads to **dephasing** at rate $\Gamma = 1/(2\tau)$. This is represented via the measurement collapse operator $C = \sqrt{\Gamma \eta}\hspace{1mm}\hat  \sigma_z$.

The balance of these two effects is determined by the signal collection efficiency $\eta$, where $\eta = 1$ means all information is collected, and $\eta = 0$ means no information is collected.
"""

# ╔═╡ b5cd2cf1-897b-477b-b2a6-9a0869d9f482
md"""
A *quantum trajectory* is a particular realization of the qubit evolution, based on the measurement record. The Hamiltonian evolution of the qubit is interleaved with Bayesian measurement updates based on information acquired from measurement. The resulting trajectory is a noisy perturbation of the ensemble-average behavior.

In the following code, we use `QuantumCircuits.jl`'s `bayesian` function to simulate trajectories for the reduced qubit system.
"""

# ╔═╡ 8b4d5601-1310-4773-b0d5-699d86789827
md" ## Simulation "

# ╔═╡ 9931a7b2-5062-48c2-ad8c-7e502fa06a9e
md"""
Change the Rabi rate: $\Omega_R$ = 
$(@bind ΩR_MHz Select(["0.5", "1.0", "1.5", "2.0"])) MHz
"""

# ╔═╡ 4de01d6c-70aa-4dd1-b497-df96124f8b8a
md"""
Change the measurement rate: $\Gamma$ = 
$(@bind Γstr Select(["0.1", "0.25", "0.5", "0.75", "1.0"])) rad MHz
"""

# ╔═╡ 8c211b48-722c-4baf-b9ec-fc14500d5695
md"""
Change the quantum efficiency: $\eta$ =
$(@bind ηstr Select(["1.0", "0.75", "0.5", "0.25", "0.0"]))
"""

# ╔═╡ 685f4b80-a9e5-4c9d-823a-819707f04601
md" # Utilities"

# ╔═╡ 868f6146-1e76-47c3-97ac-24f63b8a6b15
md"""
## Plotting
"""

# ╔═╡ b1a25a71-ba30-4d92-b7fc-cfb54836234e
colors1q = palette(:tab10)

# ╔═╡ d6e770d5-462d-4c16-bd1e-0d0a140c9a3b
begin
	green(text; title="Note") = Markdown.MD(Markdown.Admonition("correct", title, [text]))
	
	red(text; title="Note") = Markdown.MD(Markdown.Admonition("danger", title, [text]))
	
	tan(text; title="Note") = Markdown.MD(Markdown.Admonition("warning", title, [text]))
	
	blue(text; title="Note") = Markdown.MD(Markdown.Admonition("note", title, [text]))
	
	hint(text; title="Hint") = Markdown.MD(Markdown.Admonition("hint", title, [text]))
	
end

# ╔═╡ dbc1c08b-138b-4049-a57c-25376937735c
begin
	Solution = QuantumCircuits.Solution
	QOp = AbstractOperator
end

# ╔═╡ 032caad9-5875-4acc-b6e8-2cb8d4d86c35
function expectations(sol::Solution, op::QOp)
	map(state -> real(expect(op, state)), sol.ρ)
end

# ╔═╡ b4bee75b-7125-4242-8403-4ed4f7ee5129
begin
	"Single-qubit operators ------------------------------------------------"
	
	# Basis
	q = SpinBasis(1//2)
	Iq = identityoperator(q)
	
	# qubit operators, using convention that |-z> is ground state
	σx = sigmax(q)
	σy = sigmay(q)
	σz = sigmaz(q)
	σp = sigmap(q)
	σm = sigmam(q)
	n = σp * σm
	
	g = spindown(q)
	e = spinup(q)
	
	# operator basis
	qbasis = [σx, σy, σz]
	qlabels = ["x", "y", "z"]
	
end

# ╔═╡ 90d544e3-a262-4209-840a-cb796e113d92
let
	ψ0 = dm(g)
	dt = 1e-3
	tf = 10.0
	global N = 100 # number of realizations
	
	ΩR = parse(Float64, ΩR_MHz) * 2π # Rabi frequency
	Γ = parse(Float64, Γstr) # Measurement dephasing rate
	η = parse(Float64, ηstr)

	# Kraus operators -----------------------------------------
	H = (ΩR / 2) * σy
	J(η) = (η == 1.0) ? [] : [(σz, (1 - η) * Γ/2)]
	C(η) = (η == 0.0) ? [] : [(σz, Γ, η)]

	global sols = map(1:N) do m
		bayesian((0, tf), ψ0, H, J(η), C(η); dt=dt)
	end

	global η0_sol = bayesian((0, tf), ψ0, H, J(0.0), C(0.0); dt=dt)
end

# ╔═╡ 3eb969ef-b99c-4ca9-990e-238fde395bc9
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

# ╔═╡ de6758e1-8d5b-4bb8-b1fe-b6d07c1c0fee
bloch_plots(sols, η0_sol, alpha=0.15, N=N)

# ╔═╡ Cell order:
# ╟─fab5ca8b-29b3-4456-8520-c00927886612
# ╟─32746f87-236e-4e09-9e19-cfc645f7296e
# ╟─07821c37-73f0-4bb8-94c8-c843c23976b9
# ╟─ee061b17-e095-417d-956e-ae1fb058db65
# ╟─d3c789f7-6fc2-4448-85c3-c1cd51c95155
# ╟─b5cd2cf1-897b-477b-b2a6-9a0869d9f482
# ╟─8b4d5601-1310-4773-b0d5-699d86789827
# ╠═90d544e3-a262-4209-840a-cb796e113d92
# ╟─9931a7b2-5062-48c2-ad8c-7e502fa06a9e
# ╟─4de01d6c-70aa-4dd1-b497-df96124f8b8a
# ╟─8c211b48-722c-4baf-b9ec-fc14500d5695
# ╟─de6758e1-8d5b-4bb8-b1fe-b6d07c1c0fee
# ╟─685f4b80-a9e5-4c9d-823a-819707f04601
# ╠═032caad9-5875-4acc-b6e8-2cb8d4d86c35
# ╟─868f6146-1e76-47c3-97ac-24f63b8a6b15
# ╟─b1a25a71-ba30-4d92-b7fc-cfb54836234e
# ╟─d6e770d5-462d-4c16-bd1e-0d0a140c9a3b
# ╠═dbc1c08b-138b-4049-a57c-25376937735c
# ╟─3eb969ef-b99c-4ca9-990e-238fde395bc9
# ╠═b4bee75b-7125-4242-8403-4ed4f7ee5129
# ╠═3d85dac8-5abb-4622-8eeb-25afb89a71a6
