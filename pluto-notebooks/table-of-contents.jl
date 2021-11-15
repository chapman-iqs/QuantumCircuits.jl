### A Pluto.jl notebook ###
# v0.16.4

using Markdown
using InteractiveUtils

# ╔═╡ 4c45fe3e-cd69-11eb-20cd-7bfb98c040cf
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
	
	include("plotting.jl")
	
	md" ###### 🔶 Packages and julia files"
end

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents()

# ╔═╡ ae8fc1e6-55aa-4172-ad4b-11803c1b0c29
md"""
# Quantum Circuits
"""

# ╔═╡ ce729cef-e617-4b6a-a381-620dbaa8df76
md"""
## What are Quantum Circuits?
"""

# ╔═╡ 4cd7393d-b7ac-448d-978e-322a090cee52
md"""
## What is QuantumCircuits.jl?
"""

# ╔═╡ 0ffc88ee-a8cd-4d44-9c4d-af307d77155c
md"""
## Why QuantumCircuits.jl?
"""

# ╔═╡ 96dc4699-548c-4599-a9fb-ce4aa6976b1b
md"""
# Quantum computing concepts
"""

# ╔═╡ 88bdfcc9-550f-4b88-b382-c55fb3d3f710
md"""
## Quantum states and the Bloch sphere
"""

# ╔═╡ a94f7fff-c226-4b65-adb4-d733f3fd5b5b
md"""
## What can happen to a quantum state?
"""

# ╔═╡ e2ddd951-ce50-471f-976b-968b15332e7c
md"""
### Preparation
"""

# ╔═╡ 40ebd9b5-3d3b-40b7-8afc-bb17a6307dc8
md"""
### Evolution
"""

# ╔═╡ 6042b8fd-acfe-40f3-94e9-0185dc59752b
md"""
### Measurement
"""

# ╔═╡ 4b1469e7-eff2-4f29-a3e6-74366a222866
md"""
## Types of measurement
"""

# ╔═╡ 1c7fb316-14b3-4d82-81e1-6e161801a9d7
md"""
### Strong / projective measurement
"""

# ╔═╡ a87cb141-bba1-456a-9e40-857140649abc
md"""
### Weak measurement
"""

# ╔═╡ acc7f2ad-ee29-44a1-b425-924f8daf3102
md"""
# Laboratory concepts
"""

# ╔═╡ 47fe66f6-1972-417c-b5ac-a530222a7f64
md"""
## Demodulation
"""

# ╔═╡ e5ac1d91-a2d4-4d75-a731-60f0edad1bad
md"""
# Feedback
"""

# ╔═╡ Cell order:
# ╠═4c45fe3e-cd69-11eb-20cd-7bfb98c040cf
# ╠═3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─ae8fc1e6-55aa-4172-ad4b-11803c1b0c29
# ╟─ce729cef-e617-4b6a-a381-620dbaa8df76
# ╟─4cd7393d-b7ac-448d-978e-322a090cee52
# ╟─0ffc88ee-a8cd-4d44-9c4d-af307d77155c
# ╟─96dc4699-548c-4599-a9fb-ce4aa6976b1b
# ╟─88bdfcc9-550f-4b88-b382-c55fb3d3f710
# ╟─a94f7fff-c226-4b65-adb4-d733f3fd5b5b
# ╟─e2ddd951-ce50-471f-976b-968b15332e7c
# ╟─40ebd9b5-3d3b-40b7-8afc-bb17a6307dc8
# ╟─6042b8fd-acfe-40f3-94e9-0185dc59752b
# ╟─4b1469e7-eff2-4f29-a3e6-74366a222866
# ╟─1c7fb316-14b3-4d82-81e1-6e161801a9d7
# ╟─a87cb141-bba1-456a-9e40-857140649abc
# ╟─acc7f2ad-ee29-44a1-b425-924f8daf3102
# ╟─47fe66f6-1972-417c-b5ac-a530222a7f64
# ╟─e5ac1d91-a2d4-4d75-a731-60f0edad1bad
