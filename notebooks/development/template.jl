### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

# ╔═╡ 3661abf7-84e6-4e39-854b-ad7741a8ff36
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
	
	include("utilities/two-qubit-operators.jl")
	include("utilities/plotting.jl")
	include("utilities/utilities.jl")

	
	include("notebooks/table-of-contents.jl")
	include("resources.jl")

	include("notebook-utilities.jl")
	
	md" # Packages and julia files"
end

# ╔═╡ a52844ea-9a74-42b0-b37a-7cc8cd79a0a6
md"""
⏱ Tuesday June 14, 2022
"""

# ╔═╡ 01c6cbbc-1945-40bb-be1a-31abefa38bba
md"""
✍🏼 Sacha Greenfield
"""

# ╔═╡ e625f413-e93b-48c9-b14d-181c37b0fe0a
md"""
☄️ This is an template notebook that can be copied to make other notebooks.
"""

# ╔═╡ 249f33ea-92aa-4b9c-bdfa-d35fe5254225
mdp(table_of_contents📔)

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents(title="Template")

# ╔═╡ ea282b56-7b7b-48ec-b6e8-8c14ca4c5b77
md"""
# Problem setup
"""

# ╔═╡ c4a96416-0def-479b-b7f9-a62de40d3ece
md"""
## System parameters
"""

# ╔═╡ f93517db-d1b0-4b34-91be-b0b027665f36
md"""
## Simulation
"""

# ╔═╡ b3a0454c-c016-4f46-84d6-5c3d298bd69a
md"""
# Resources
"""

# ╔═╡ 06574d54-07fa-403c-8842-799b82697548
begin
	Slichter_et_al = html"<a href='https://iopscience.iop.org/article/10.1088/1367-2630/18/5/053031/pdf'>Slichter et al., 2016 📘</a>"

	md" ♦️ **URLs**"
end

# ╔═╡ 8b05fb26-4195-49f3-84bb-0d8968c2fce8
md"""
 $Slichter_et_al D. H. Slichter, C. Müller, R. Vijay, S. J. Weber, A. Blais, and I. Siddiqi, Quantum Zeno Effect in the Strong Measurement Regime of Circuit Quantum Electrodynamics, New J. Phys. 18, 053031 (2016).
"""

# ╔═╡ Cell order:
# ╟─a52844ea-9a74-42b0-b37a-7cc8cd79a0a6
# ╟─01c6cbbc-1945-40bb-be1a-31abefa38bba
# ╟─e625f413-e93b-48c9-b14d-181c37b0fe0a
# ╟─249f33ea-92aa-4b9c-bdfa-d35fe5254225
# ╟─3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─ea282b56-7b7b-48ec-b6e8-8c14ca4c5b77
# ╟─c4a96416-0def-479b-b7f9-a62de40d3ece
# ╟─f93517db-d1b0-4b34-91be-b0b027665f36
# ╟─b3a0454c-c016-4f46-84d6-5c3d298bd69a
# ╟─06574d54-07fa-403c-8842-799b82697548
# ╠═8b05fb26-4195-49f3-84bb-0d8968c2fce8
# ╠═3661abf7-84e6-4e39-854b-ad7741a8ff36
