### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 8495a543-be81-4ee7-b4f1-1aac51f71235
begin
	cd("/Users/sachagreenfield/Desktop/GitHub/QuantumCircuits.jl")
	import Pkg
	Pkg.activate(".")
	using QuantumCircuits
end

# ╔═╡ e1676dd6-cd7b-11eb-29a1-6dd4f65df6ea
# using QuantumCircuits

# ╔═╡ 927b582e-acd8-4865-8f5a-9e985de45a25
function f(x, r; str="this is broadcast str")
	return str
end

# ╔═╡ a20cb006-cd0d-4a90-9781-71b86d409441
function f(x; pig="this is broadcast pig")
	return pig
end

# ╔═╡ 4484fd9c-d60b-4ca3-86ce-6109160d4907
rouchon

# ╔═╡ 1c4b8b81-3e4b-4c10-b9ed-3c2cfef03994
names(QuantumCircuits)

# ╔═╡ Cell order:
# ╠═e1676dd6-cd7b-11eb-29a1-6dd4f65df6ea
# ╠═8495a543-be81-4ee7-b4f1-1aac51f71235
# ╠═927b582e-acd8-4865-8f5a-9e985de45a25
# ╠═a20cb006-cd0d-4a90-9781-71b86d409441
# ╠═4484fd9c-d60b-4ca3-86ce-6109160d4907
# ╠═1c4b8b81-3e4b-4c10-b9ed-3c2cfef03994
