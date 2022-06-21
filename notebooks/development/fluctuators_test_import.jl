### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

# ╔═╡ 4ff4574e-5cdb-4b06-b2e8-836ba918eb38
include("simulations/fluctuators/fluctuators.jl")

# ╔═╡ 353c68d0-f1a0-11ec-083c-2b15358477dd
function homepwd(dirname = "QuantumCircuits.jl")
	path = let 
			arr = split(pwd(), "/")
			index = findfirst(s -> s == dirname, arr)
			join(map(a -> string(a, "/"), arr[1:index]))
	end
	cd(path)
	pwd()
end

# ╔═╡ e6b0f01c-2e3f-454e-abb4-0344bd0625c7
begin
	homepwd()
	
	import Pkg
	Pkg.activate(".")
	using QuantumCircuits

	include("simulations/fluctuators/FluctuatorSims.jl")
end

# ╔═╡ 832c65a4-bf89-4525-9c70-6bed8feb15cf
import .FluctuatorSims

# ╔═╡ Cell order:
# ╠═353c68d0-f1a0-11ec-083c-2b15358477dd
# ╠═e6b0f01c-2e3f-454e-abb4-0344bd0625c7
# ╠═4ff4574e-5cdb-4b06-b2e8-836ba918eb38
# ╠═832c65a4-bf89-4525-9c70-6bed8feb15cf
