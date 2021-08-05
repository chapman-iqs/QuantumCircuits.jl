### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 0bcc46c4-ce48-11eb-24ce-49f8706b7525
begin
	using Random
	using Statistics
	using PyPlot
	using QuantumCircuits 	# note: for some reason, QuantumCircuits has to be 
							# "used" last
end

# ╔═╡ 5d30c7c1-4c53-40d2-8e73-b58db5697cfa
md"""
# Rouchon integration method

In this interactive notebook, we describe the Rouchon integration method as implemented in `QuantumCircuits.jl`.
"""

# ╔═╡ 2b27ab1e-f488-4565-a5b9-35eff5e37a56
md" ### Rouchon approximation "

# ╔═╡ addbae24-128e-426c-bbf6-d430e381f7e3
md""" 
The Rouchon method was presented in [1] as a numerical integration method that provides an efficient alternative to the stochastic master equation (SME). Unlike the SME method, Rouchon's method preserves positivity of the conditioned quantum state by expanding to second order in the Wiener increment $\Delta W_r$.

"""

# ╔═╡ 9e070a86-1816-4064-bbf9-acc519bc9849
md"""
$M(t) = \mathbb{1} - \Big (i H + \frac12 \sum_j J_j^\dagger J_j ) dt + \sum_c \sqrt{\eta_c} C_c dy_c(t) +  \sum_{c,d} \frac{\sqrt{\eta_c \eta_d}}2 C_c C_d \big(dy_c(t) dy_d(t) - \delta_{c,d} dt \big)$
"""

# ╔═╡ 7df8377e-9373-49d1-8bcd-2ab40848fcd9
md"""
$dy_c(t) = \sqrt{\eta_c} \text{Tr} \big ( C_c \rho(t) + \rho(t) C_c^\dagger \big ) dt + dW_c(t)$
"""

# ╔═╡ 373f2b0a-0873-4975-a8e2-c83126296aa6
md"""
$D(t) = \Big (\sum_j J_j \rho(t) J_j^\dagger  + \sum_c C_c \rho(t) C_c^\dagger \Big ) dt$
"""

# ╔═╡ f467415d-e02c-4b3b-8925-5ab0bf291d9b
md"""
$ \rho(t + dt) = \frac{M(t) \rho(t) M(t)^\dagger + D(t)}{\text{Tr} \Big ( \frac{M(t) \rho(t) M(t)^\dagger + D(t) \Big )}$
"""

# ╔═╡ 676e029d-1141-41b2-bc81-5d2e4500c684
md"""
## Timing comparison

Rouchon is faster for given $\Delta t$, but Bayesian is more robust to large $\Delta t$.
"""

# ╔═╡ e7ffdfb4-622e-4248-a71b-907ea1cbb8a5
md"""

### References

[1] P. Rouchon and J. F. Ralph, Efficient Quantum Filtering for Quantum Feedback Control, Phys. Rev. A 91, 012118 (2015).

"""

# ╔═╡ Cell order:
# ╠═0bcc46c4-ce48-11eb-24ce-49f8706b7525
# ╟─5d30c7c1-4c53-40d2-8e73-b58db5697cfa
# ╟─2b27ab1e-f488-4565-a5b9-35eff5e37a56
# ╟─addbae24-128e-426c-bbf6-d430e381f7e3
# ╟─9e070a86-1816-4064-bbf9-acc519bc9849
# ╟─7df8377e-9373-49d1-8bcd-2ab40848fcd9
# ╟─373f2b0a-0873-4975-a8e2-c83126296aa6
# ╠═f467415d-e02c-4b3b-8925-5ab0bf291d9b
# ╟─676e029d-1141-41b2-bc81-5d2e4500c684
# ╟─e7ffdfb4-622e-4248-a71b-907ea1cbb8a5
