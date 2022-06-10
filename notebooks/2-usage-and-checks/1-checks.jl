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

# ╔═╡ 41332d40-c450-4346-93ea-cdb6d71d82eb
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

	include("utilities/single-qubit-operators.jl")
	include("utilities/utilities.jl")
	include("utilities/plotting.jl")

	include("notebooks/table-of-contents.jl")
	include("resources.jl")
	
	
end

# ╔═╡ 75236fca-cfac-11eb-2ffe-c394a1c504cf


# ╔═╡ 32746f87-236e-4e09-9e19-cfc645f7296e
TableOfContents(title="Solver checks")

# ╔═╡ 23079100-19c4-4031-b75c-c11855d5ed69
mdp(table_of_contents📔)

# ╔═╡ 337d1109-0183-4411-b269-34aeb03961b8
md"""
# Ensemble-average convergence

The ensemble-average should be independent of $\eta$, and should converge to the $\eta = 0$ solution (equivalent to the master equation solution). Below we consider the example of the reduced qubit.

"""

# ╔═╡ 8b4d5601-1310-4773-b0d5-699d86789827
md" ## Simple qubit"

# ╔═╡ b8838906-9d47-4a76-bf68-32c27b011aac
md"""
Choose the simulation method:
$(@bind solve_str Select(["bayesian", "rouchon"]))
"""

# ╔═╡ ca047367-c863-408f-9194-ce53732a8fb3
solve = (solve_str == "bayesian") ? bayesian : rouchon

# ╔═╡ 8c211b48-722c-4baf-b9ec-fc14500d5695
md"""
Change the quantum efficiency to check invariance of the ensemble average: $\eta$ =
$(@bind ηstr Select(["1.0", "0.75", "0.5", "0.25", "0.0"]))
"""

# ╔═╡ 90d544e3-a262-4209-840a-cb796e113d92
let
	ψ0 = g
	dt = 1e-3
	tf = 10.0
	global N = 100 # number of realizations
	
	ΩR = 2π # Rabi frequency
	Γ = 0.15 # Measurement dephasing rate
	η = parse(Float64, ηstr)

	# Kraus operators -----------------------------------------
	H = (ΩR / 2) * σy
	J(η) = (η == 1.0) ? [] : [(σz, (1 - η) * Γ/2)]
	C(η) = (η == 0.0) ? [] : [(σz, Γ, η)]

	global sols = map(1:N) do m
		solve((0, tf), ψ0, H, J(η), C(η); dt=dt)
	end

	global η0_sol = solve((0, tf), ψ0, H, J(0.0), C(0.0); dt=dt)
end

# ╔═╡ 549c46e5-f1d2-463b-b355-ccde14fe4249
md"""
Try changing the value of `η` in the code above to see how it affects the ensemble average.
"""

# ╔═╡ cb42e11d-85c6-425b-a48a-99087b695c5c
md" ## Time-step convergence"

# ╔═╡ a1b4f410-acee-402c-8e2c-1ccf5468e88a
md"""
Choosing the right time-step is critical to getting physically realistic simulation results. In weakly measured systems, there is a natural finite timestep of the experiment due to the detector bandwidth BW ($\Delta t_{detector} = 1/(2 \text{BW})$). Ideally, your simulation timestep $dt$ should match this experimental time-resolution. However, depending on your simulation method, choosing a large $dt$ may lead to inaccurate results if your $H$, $J$, and $C$ operators are non-commuting on that timescale.

Below are some checks to see if the integration is robust at the desired time-step.

"""

# ╔═╡ 5b49af65-c3f5-408c-9f8b-dd9a409a7ba7
md" ##### single-trajectory convergence across $\Delta t$"

# ╔═╡ bb7c273b-aabe-439e-a858-82fd2d5471f5
md"""
A single trajectory (correpsonding to a given measurement record) should be robust to changes in $dt$.
"""

# ╔═╡ 685f4b80-a9e5-4c9d-823a-819707f04601
md" # Utilities"

# ╔═╡ 6c288136-8334-455f-94c1-72b46c41d7d0
md"""
## Misc
"""

# ╔═╡ ad5b0f4f-f5f8-4669-97ef-054a82a85feb
ηdict = Dict(["η = 1.0" => 1.0, 
				"η = 0.75" => 0.75, 
				"η = 0.5" => 0.5, 
				"η = 0.25" => 0.25, 
				"η = 0.0" => 0.0])

# ╔═╡ a887b375-b80a-4d99-862c-2b9455afffba
begin
	qbasis = [σx, σy, σz]
	qlabels = ["x", "y", "z"]
end

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

# ╔═╡ d09a282e-4b9b-4718-89bf-62c42603f937
hint(md"If there seems to be a mismatch, try increasing $N$ (the number of trajectories) simulated in the ensemble. The larger $\eta$ is, the greater $N$ must be to see good convergence.")

# ╔═╡ 34a3b82f-7722-407a-b703-555ad2efb59c
green(md"If you're not getting ensemble-average convergence in your own simulations, this could be for a few reasons: (i) you need to simulate more trajectories to get convergence (especially if $\eta$ is close to 1)
	; (ii) you are having time-step convergence issues (see next section); or (iii) there is a mismatch between your measurement operators `C` and your dissipation operators `J`.
	"; title="Tip")

# ╔═╡ daf693cb-0715-46a9-98b0-c3b3d6f2fd2c
# function bloch_plots(sols::Vector{Solution}, sol_η0::Solution)
# 	colors = colors1q

# 	# calculate expectation values --------------------------------------------
# 	t = sols[1].t
# 	xs, ys, zs = [], [], []

# 	for sol in sols
# 		x, y, z = map(op -> expectations(sol, op), qbasis)
# 		for (list, traj) in zip([xs, ys, zs], [x, y, z])
# 			push!(list, traj)
# 		end

# 	end

# 	# η = 0 solution
# 	xη0, yη0, zη0 = map(op -> expectations(sol_η0, op), qbasis)



# 	# plot ----------------------------------------------------------------------
# 	function bloch(os, oη0; color=colors1q[1], xlabel="", ylabel="")

# 		po = plot(ylims = [-1,1], xlabel=xlabel, ylabel=ylabel)

# 		for o in os[1:min(N, 50)]
# 			plot!(t, o, alpha=0.1, label=:none, color=color)
# 		end

# 		oavg = [mean([os[i][j] for i in 1:N]) for j in 1:length(t)]
# 		plot!(t, oavg, alpha=1, color=color, label="average", linewidth=3)
# 		plot!(t, oη0, alpha=1, color=:black, label="η = 0", linewidth=2)

# 		po

# 	end


# 	l = @layout [xplot{0.33h}; yplot{0.33h}; zplot{0.33h}]

# 	px = bloch(xs, xη0, color=colors[1], ylabel="x")
# 	py = bloch(ys, yη0, color=colors[2], ylabel="y")
# 	pz = bloch(zs, zη0, color=colors[3], ylabel="z", xlabel="t (μs)")

# 	plot(px, py, pz, layout = l, link=:y, size=(800,500), legendfontsize=8, titlefontsize=12, legend=:outerright)
# end


# ╔═╡ 3a51d3be-248c-44b2-829a-29cdc0cda6d3
colorangle(ϕ) = RGB(200/255,abs(sin(ϕ))*200/255,200/255)

# ╔═╡ cec64cc5-ee9d-4f2b-a44d-53d52ea64a1b
function bloch_plots(sols::Vector{Solution}, sol_η0::Solution; record=false, size=(800,500))
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
			plot!(t, o, alpha=0.1, label=:none, color=color)
		end

		oavg = [mean([os[i][j] for i in 1:N]) for j in 1:length(t)]
		plot!(t, oavg, alpha=1, color=color, label="average", linewidth=3)
		plot!(t, oη0, alpha=1, color=:black, label="η = 0", linewidth=2)

		po

	end


	px = bloch(xs, xη0, color=colors[1], ylabel="x")
	py = bloch(ys, yη0, color=colors[2], ylabel="y")
	pz = bloch(zs, zη0, color=colors[3], ylabel="z", xlabel="t (μs)")

	if !record

		l = @layout [xplot{0.33h}; yplot{0.33h}; zplot{0.33h}]
		plot(px, py, pz, layout = l, link=:y, size=size, legendfontsize=8, titlefontsize=12, legend=:outerright)

	else
		rs = map(sol -> sol.r[1], sols)
		pr = plot(ylabel="r", xlabel="t (μs)")

		for r in rs[1:min(N, 50)]
			plot!(t, r, alpha=0.05, color=colorangle(0), label=:none)
		end
		ravg = [mean([rs[i][j] for i in 1:N]) for j in 1:length(t)]
		ravg = ravg .* (4/5 * maximum(rs[1]) / maximum(ravg))
		plot!(t, ravg, alpha=1, color=:white, label="average", linewidth=1)
		
		l = @layout [xplot{0.2h}; yplot{0.2h}; zplot{0.2h}; rplot{0.4h}]
		plot(px, py, pz, pr, layout = l, link=:y, size=size, legendfontsize=8, titlefontsize=12, legend=:outerright)

	end



	

	
end


# ╔═╡ de6758e1-8d5b-4bb8-b1fe-b6d07c1c0fee
bloch_plots(sols, η0_sol, record=true, size=(800, 800))

# ╔═╡ 8245a294-8615-4baf-8484-89d133a25230
# for plotting ensmebles of single-qubit trajectories
function bloch_plots_records(sols::Vector{Solution}, sol_η0::Solution)
	colors = colors1q

	# calculate expectation values --------------------------------------------
	t = sols[1].t
	xs, ys, zs = [], [], []
	rs = map(sol -> sol.r[1], sols)

	for sol in sols
		x, y, z = map(op -> expectations(sol, op), qbasis)
		for (list, traj) in zip([xs, ys, zs], [x, y, z])
			push!(list, traj)
		end

	end

	# η = 0 solution
	xη0, yη0, zη0 = map(op -> expectations(sol_η0, op), qbasis)


	# plot ----------------------------------------------------------------------
	function bloch(os, oη0, color=colors1q[1], xlabel="", ylabel="")

		po = plot(ylims = [-1,1], xlabel=xlabel, ylabel=ylabel)

		for o in os[1:min(N, 50)]
			plot!(t, o, alpha=0.1, label=:none, color=color)
		end

		oavg = [mean([os[i][j] for i in 1:N]) for j in 1:length(t)]
		plot!(t, oavg, alpha=1, color=color, label="average", linewidth=3)
		plot!(t, oη0, alpha=1, color=:black, label="η = 0", linewidth=2)

		po

	end


	l = @layout [xplot{0.25h}; yplot{0.25h}; zplot{0.25h}; rplot{0.25h}]

	px = bloch(xs, xη0, color=colors[1], ylabel="x")
	py = bloch(ys, yη0, color=colors[2], ylabel="y")
	pz = bloch(zs, zη0, color=colors[3], ylabel="z")
	pr = plot(ylabel="measurement record", xlabel="t (μs)")

	for r in rs
		plot!(t, r, alpha=0.1, color=colorangle(0))
	end

	plot(px, py, pz, pr, layout = l, link=:y, size=(800,500), legendfontsize=8, titlefontsize=12, legend=:outerright)
end

# ╔═╡ c1b40baa-0d31-4caf-8bb8-22d7ae5f97fc
bloch_plots_records(sols, η0_sol)

# ╔═╡ Cell order:
# ╠═75236fca-cfac-11eb-2ffe-c394a1c504cf
# ╟─32746f87-236e-4e09-9e19-cfc645f7296e
# ╟─23079100-19c4-4031-b75c-c11855d5ed69
# ╟─337d1109-0183-4411-b269-34aeb03961b8
# ╟─8b4d5601-1310-4773-b0d5-699d86789827
# ╟─b8838906-9d47-4a76-bf68-32c27b011aac
# ╟─ca047367-c863-408f-9194-ce53732a8fb3
# ╠═90d544e3-a262-4209-840a-cb796e113d92
# ╟─8c211b48-722c-4baf-b9ec-fc14500d5695
# ╠═de6758e1-8d5b-4bb8-b1fe-b6d07c1c0fee
# ╠═c1b40baa-0d31-4caf-8bb8-22d7ae5f97fc
# ╟─549c46e5-f1d2-463b-b355-ccde14fe4249
# ╟─d09a282e-4b9b-4718-89bf-62c42603f937
# ╟─34a3b82f-7722-407a-b703-555ad2efb59c
# ╟─cb42e11d-85c6-425b-a48a-99087b695c5c
# ╟─a1b4f410-acee-402c-8e2c-1ccf5468e88a
# ╟─5b49af65-c3f5-408c-9f8b-dd9a409a7ba7
# ╟─bb7c273b-aabe-439e-a858-82fd2d5471f5
# ╟─685f4b80-a9e5-4c9d-823a-819707f04601
# ╟─6c288136-8334-455f-94c1-72b46c41d7d0
# ╟─ad5b0f4f-f5f8-4669-97ef-054a82a85feb
# ╠═a887b375-b80a-4d99-862c-2b9455afffba
# ╟─868f6146-1e76-47c3-97ac-24f63b8a6b15
# ╟─b1a25a71-ba30-4d92-b7fc-cfb54836234e
# ╟─d6e770d5-462d-4c16-bd1e-0d0a140c9a3b
# ╟─daf693cb-0715-46a9-98b0-c3b3d6f2fd2c
# ╠═cec64cc5-ee9d-4f2b-a44d-53d52ea64a1b
# ╠═8245a294-8615-4baf-8484-89d133a25230
# ╠═3a51d3be-248c-44b2-829a-29cdc0cda6d3
# ╠═41332d40-c450-4346-93ea-cdb6d71d82eb
