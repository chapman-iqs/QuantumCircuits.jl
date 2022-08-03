### A Pluto.jl notebook ###
# v0.17.5

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

# ╔═╡ b1210b2c-21af-11ec-1019-fd2882da1e28
begin
	cd("/Users/sachagreenfield/Desktop/GitHub/QC-notebooks")
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
	using CSV
	
	include("plotting.jl")
	include("utilities.jl")
	
	md" ###### 🔶 Packages and julia files"
end

# ╔═╡ 8916212e-1d1b-4560-ab01-ebb10ff24905
md"""
### Load Data
"""

# ╔═╡ ce78f5b7-1882-4808-b87b-3fc929de8700
md"""
Start by loading in the data from CSV. This notebook is formatted to load data exported by `QuantumCircuits.jl/simulations/ensemble.jl`.
"""

# ╔═╡ a3ace32d-7fc3-4f93-8e6c-18c430266042
begin
	path_sf = string("test/trajectory_test_data-system-filter/", "2022-01-14T15-12-42.391", "/")
	path_main = string("test/trajectory_test_data-main/", "2022-01-14T14-25-52.617", "/")
end

# ╔═╡ 695757a3-e47d-4892-b211-a62acb4cc46e
function get_trajs(X, Y, Z, R, T)
	
	t = T[!, 1]
	trajs = map(1:length(X[1,:])) do i
				x, y, z, r = [o[!,i] for o in (X, Y, Z, R)]
				p = 0.5 .* (1 .+ x.^2 .+ y.^2 .+ z.^2)
				traj(t, x, y, z, p, r)
			end
	
	return trajs
end	

# ╔═╡ a666c1ab-32e7-4d51-a0a1-f5893b491a7c
function get_avg(X, Y, Z, R, T)

	t = T[!, 1]
	ens_avg = 	begin
					x, y, z = [mean.(eachrow(A)) for A in (X, Y, Z)]
					p = 0.5 .* (1 .+ x.^2 .+ y.^2 .+ z.^2)
					traj(t, x, y, z, p, ())
				end

	return ens_avg
end

# ╔═╡ 123e2852-8e45-4766-a7e6-c1bd022c3edf
function get_pars(DF::DataFrame)

	pairs = map(1:nrow(DF)) do i
			pair = copy(DF[i,:])
			k = pair[1]
			v = try 
					parse(Float64, pair[2]) 
				catch 
					pair[2] 
				end
			(k, v)
	end	
	
	return Dict(pairs)

end

# ╔═╡ a0d49386-63d1-4ac9-9780-d04e3596669a
function load(path::String)
	
	# load DataFrames
	parDF = DataFrame(CSV.File(string(path, "parameters.csv")))	
	X = DataFrame(CSV.File(string(path, "x.csv")))
	Y = DataFrame(CSV.File(string(path, "y.csv")))
	Z = DataFrame(CSV.File(string(path, "z.csv")))
	R1 = DataFrame(CSV.File(string(path, "r1.csv")))
	T = DataFrame(CSV.File(string(path, "t.csv")))

	pars = get_pars(parDF)
	trajs = get_trajs(X, Y, Z, R1, T)
	ensavg = get_avg(X, Y, Z, R1, T)

	(pars, trajs, ensavg)

end

# ╔═╡ 524ae1e4-5e3c-47b8-bc7e-1725ebc50694
md"""
###### Simulation parameters
"""

# ╔═╡ a076c21e-ee90-4619-b253-f290db2092a4
begin
	(pars_sf, trajs_sf, ens_avg_sf) = load(path_sf)
	(pars_main, trajs_main, ens_avg_main) = load(path_main)
end

# ╔═╡ 5c9940d0-7abf-45ad-ab5e-e9a326f2aa08
@bind i Slider(1:100)

# ╔═╡ 2c2d2421-3c4b-4a45-a216-628a409fc6dc
begin
	sim = trajs_sf[i]
	blochtimeseries(sim.t, sim.x, sim.y, sim.z, size=(800,300), legend=:outerright, title=string("trajectory ", i, ", main branch"))
end

# ╔═╡ d5f82d6f-8483-41c8-a33c-6f688b3a5ae4
trajs_sf[i].r

# ╔═╡ 654239f1-256e-4173-8f2f-aacd3d9b0377
trajs_main[i].r

# ╔═╡ 37524560-ebfd-409b-bd7a-734e876b0b62
md"""
# Utilities
"""

# ╔═╡ 8b7e05dd-94e6-4ca5-b82e-9c9ae350ca5f
colors = palette(:tab10)

# ╔═╡ 2955f826-8768-48f2-a6e3-49fea998eb1a
function plotresults(traj1::traj, traj2::traj)
	# plot ----------------------------------------------------------------------
	l = @layout [blochs{0.5h}; records{0.5h}]

	p1 = plot(traj1.t, traj1.x, color=colors[1], label=L"x", legend=:outerright)
	plot!(traj1.t, traj1.y, color=colors[2], label=L"y")
	plot!(traj1.t, traj1.z, color=colors[3], label=L"z")
	plot!(traj2.t, traj2.x, color=colors[1], linestyle=:dash, label=:none)
	plot!(traj2.t, traj2.y, color=colors[2], linestyle=:dash, label=:none)
	plot!(traj2.t, traj2.z, color=colors[3], linestyle=:dash, label=:none)

	p2 = plot(traj1.t, traj1.r, color=:red, label=:none)
	plot!(traj2.t, traj2.r, color=:pink, linestyle=:dash, label=:none)

	plot(p1, p2, layout = l, link=:y, size=(600,800), legendfontsize=8, titlefontsize=12, legend=:outerright, xlabel="t (μs)")
	
end

# ╔═╡ 2b02acd1-b60e-4e3c-8a67-fded1980c80a
plotresults(trajs_sf[i], trajs_main[i])

# ╔═╡ Cell order:
# ╠═b1210b2c-21af-11ec-1019-fd2882da1e28
# ╟─8916212e-1d1b-4560-ab01-ebb10ff24905
# ╟─ce78f5b7-1882-4808-b87b-3fc929de8700
# ╠═a3ace32d-7fc3-4f93-8e6c-18c430266042
# ╠═a0d49386-63d1-4ac9-9780-d04e3596669a
# ╠═695757a3-e47d-4892-b211-a62acb4cc46e
# ╠═a666c1ab-32e7-4d51-a0a1-f5893b491a7c
# ╠═123e2852-8e45-4766-a7e6-c1bd022c3edf
# ╟─524ae1e4-5e3c-47b8-bc7e-1725ebc50694
# ╠═a076c21e-ee90-4619-b253-f290db2092a4
# ╠═2c2d2421-3c4b-4a45-a216-628a409fc6dc
# ╠═5c9940d0-7abf-45ad-ab5e-e9a326f2aa08
# ╠═2b02acd1-b60e-4e3c-8a67-fded1980c80a
# ╠═d5f82d6f-8483-41c8-a33c-6f688b3a5ae4
# ╠═654239f1-256e-4173-8f2f-aacd3d9b0377
# ╟─37524560-ebfd-409b-bd7a-734e876b0b62
# ╠═8b7e05dd-94e6-4ca5-b82e-9c9ae350ca5f
# ╠═2955f826-8768-48f2-a6e3-49fea998eb1a
