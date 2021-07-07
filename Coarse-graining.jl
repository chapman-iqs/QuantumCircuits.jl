### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ c520963b-d160-431d-a834-8e58185d0ae4
begin
	using Random
	using Statistics
	using PyPlot
	using QuantumCircuits 	# note: for some reason, QuantumCircuits has to be 
							# "used" last
end

# ╔═╡ 18ad70f0-d08e-11eb-374b-ad929f1d2a18
md" # Coarse-graining"

# ╔═╡ 45f9e3a9-aaac-414a-8852-d9959d4f4bfa
md"""
In this interactive notebook, we look at the implementation of coarse-graining in `QuantumCircuits.jl`.
"""

# ╔═╡ Cell order:
# ╠═c520963b-d160-431d-a834-8e58185d0ae4
# ╠═18ad70f0-d08e-11eb-374b-ad929f1d2a18
# ╟─45f9e3a9-aaac-414a-8852-d9959d4f4bfa
