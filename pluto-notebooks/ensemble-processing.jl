### A Pluto.jl notebook ###
# v0.16.0

using Markdown
using InteractiveUtils

# ╔═╡ b1210b2c-21af-11ec-1019-fd2882da1e28
begin
	cd("/Users/sachagreenfield/Desktop/GitHub/QuantumCircuits.jl")
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
loadpath = string("data/", "2021-10-01T13-09-09.342", "/")

# ╔═╡ 524ae1e4-5e3c-47b8-bc7e-1725ebc50694
md"""
###### Simulation parameters
"""

# ╔═╡ 75e1d2dd-06a2-496b-b69f-c52aaf206563
begin
	parDF = DataFrame(CSV.File(string(loadpath, "parameters.csv")))	
	
	pairs = map(1:nrow(parDF)) do i
			pair = copy(parDF[i,:])
			k = pair[1]
			v = try 
					parse(Float64, pair[2]) 
				catch 
					pair[2] 
				end
			(k, v)
	end	
	
	pars = Dict(pairs)
end

# ╔═╡ 5002b4ce-02f6-4b8e-a50e-1710bb051ec8
dualquad = (pars["dualquad"] == "true")

# ╔═╡ 64d41305-e2df-4633-8882-93ffbb9789ec
begin
	X = DataFrame(CSV.File(string(loadpath, "x.csv")))
	Y = DataFrame(CSV.File(string(loadpath, "y.csv")))
	Z = DataFrame(CSV.File(string(loadpath, "z.csv")))
	R1 = DataFrame(CSV.File(string(loadpath, "r1.csv")))
	if dualquad
		R2 = DataFrame(CSV.File(string(loadpath, "r2.csv")))
	end
	T = DataFrame(CSV.File(string(loadpath, "t.csv")))
	
	md" ###### 🟡 Load data frames"
end

# ╔═╡ c2d03585-5c34-4fe0-8f31-46c3eaebf58d
parDF

# ╔═╡ 8beea68b-c83e-4f4c-894c-40406713bf95
md"""
### Process data
"""

# ╔═╡ cafec9e4-102a-4e8c-a350-0ae626267f95
begin
	
	# process data as trajectories 
	t = T[!, 1]
	
	trajs = map(1:length(X[1,:])) do i
				if dualquad
					x, y, z, r1, r2 = [o[!,i] for o in (X, Y, Z, R1, R2)]
					r = [r1, r2]
				else
					x, y, z, r = [o[!,i] for o in (X, Y, Z, R1)]
				end
				
				p = 0.5 .* (1 .+ x.^2 .+ y.^2 .+ z.^2)
				traj(t, x, y, z, p, r)
	end
	
	# take ensemble average
	ens_avg = 	begin
					x, y, z = [mean.(eachrow(A)) for A in (X, Y, Z)]
					p = 0.5 .* (1 .+ x.^2 .+ y.^2 .+ z.^2)
					traj(t, x, y, z, p, ())
				end
end

# ╔═╡ fc62283c-6cb6-42ff-8c4d-7baeb542c19c
blochtimeseries(ens_avg.t, ens_avg.x, ens_avg.y, ens_avg.z)

# ╔═╡ d34de629-2dd8-4bf3-ba92-1c3204839f96
md"""
### Stabilization results
"""

# ╔═╡ 7a4b076d-d0fe-459b-81d8-0507652e38a8
begin
	η = pars["η"]
	T1 = pars["Τ1"]
	τm = pars["τm"]
	rmax = 1/sqrt(1/η + 2τm/T1)
	fmax  = (1 + rmax)/2
	pmax = (1 + rmax^2)/2
	r = sqrt(last(x)^2 + last(y)^2 + last(z)^2)
	f = (1+r)/2
	pur = last(ens_avg.p)
end

# ╔═╡ 1b4a4e82-ac76-4a74-ac7d-008693531f0d
md"""
Parameters: η = $(pars["η"]), τm = $(pars["τm"]) μs, T1 = $(pars["Τ1"]) μs, T2 = $(pars["Τ2"]) μs, Td = $(pars["td"]) μs, $(dualquad ? "dual quadrature" : "single quadrature"), ΔS = $(dualquad ? ΔS : "N/A")

The stabilized Bloch radius was R = $(round(r, digits=3)), compared to the maximum bloch radius for these parameters: Rmax = $(round(rmax, digits=3)).

This corresponds to fidelity F = $(round(f, digits=3)), compared to maximum fidelity Fmax = $(round(fmax, digits=3)), and purity P = $(round(pur, digits=3)), compared to maximum purity Pmax = $(round(pmax, digits=3)).

Data found in $loadpath
"""

# ╔═╡ Cell order:
# ╠═b1210b2c-21af-11ec-1019-fd2882da1e28
# ╟─8916212e-1d1b-4560-ab01-ebb10ff24905
# ╟─ce78f5b7-1882-4808-b87b-3fc929de8700
# ╠═a3ace32d-7fc3-4f93-8e6c-18c430266042
# ╟─524ae1e4-5e3c-47b8-bc7e-1725ebc50694
# ╠═75e1d2dd-06a2-496b-b69f-c52aaf206563
# ╠═5002b4ce-02f6-4b8e-a50e-1710bb051ec8
# ╠═64d41305-e2df-4633-8882-93ffbb9789ec
# ╠═c2d03585-5c34-4fe0-8f31-46c3eaebf58d
# ╟─8beea68b-c83e-4f4c-894c-40406713bf95
# ╠═cafec9e4-102a-4e8c-a350-0ae626267f95
# ╠═fc62283c-6cb6-42ff-8c4d-7baeb542c19c
# ╟─d34de629-2dd8-4bf3-ba92-1c3204839f96
# ╟─1b4a4e82-ac76-4a74-ac7d-008693531f0d
# ╠═7a4b076d-d0fe-459b-81d8-0507652e38a8
