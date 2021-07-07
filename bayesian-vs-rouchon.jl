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
# Bayesian vs. Rouchon integration methods

In this interactive notebook, we'll compare quantitative and qualitative features of the two primary integration methods -- `bayesian` and `rouchon` -- implemented in `QuantumCircuits.jl`. What follows is a brief summary of key commonalities and differences before we dive into the details.
"""

# ╔═╡ d4e7b79b-eb8c-48b3-9a94-650481fa8cb7
md"""
##### Commonalities

The `bayesian` and `rouchon` methods have important commonalities. 

Both:
* Model quantum trajectories of weakly and continuously measured open quantum systems.

* Take `H`, `[J]`, `[C]` (Hamiltonian, Lindblad, collapse) operators as inputs.

* If simulating quantum trajectories, can generate a noisy record of a monitored variable $r = \langle x \rangle + dW/dt$, with Wiener process $dW \sim \mathcal{N}(0, \sqrt{\tau_m})$ (**CHECK THIS**) if running a simulation. If reconstructing experimental trajectories, can process an imported experimental record.

"""

# ╔═╡ b6c73eea-d17c-49d4-a091-be10c8134718
md" ### Bayesian filtering "

# ╔═╡ 4fc42181-a1d4-46d9-be76-7277b9a2d6a5
md" ##### Update equations"

# ╔═╡ 343161e9-e975-45b5-a444-b9f9f0b30cfb
md"""
Bayesian filtering describes the evolution of a weakly measured quantum system in fundamentally discrete chunks of time, $\Delta t$. In general, the system is described by a density matrix $\rho$. It undergoes **unitary evolution** according to a  Hamiltonian $H$, such that

$\hat \rho' = \hat U \rho {\hat U}^\dagger,$

where 

$\hat U = \exp(- i \Delta t \hat H).$

(Note that $\hbar$ has been absorbed into $\hat H$ here.) 


"""

# ╔═╡ d5af31e4-ba93-4f3f-9c55-4f2ce600d797
md"""
The **weak measurement** of an observable $\hat C$ is represented by a Kraus operator $\hat M_{r_C}$, such that

$\hat M_{r_C} = \exp\Big[\frac{\Delta t}{2 \tau_m} {r_C}^* \hat C  - \frac14 (\hat C + \hat C^\dagger)^2\Big].$

"""

# ╔═╡ 97e2560f-3dd9-45a0-85f5-aa9e12800292
md"""

Here, $r_C$ is a complex noisy record associated with measuring $\hat C$, while $\tau_m$ is the timescale of measurement collapse. From each $\hat M_{r_C}$ we can construct a POVM (positive operator valued measure)

$E_{r_C} = \hat M_{r_C}^\dagger \hat M_{r_C}$

satisfying

$\sum_{r_C} E_{r_C} = \mathbb{1}.$

(There are, of course, infinitely many possible records $r_C$ corresponding to different noise realizations. So this might not be the right way to put it.) 

The weak measurement updates the state according to

$\hat \rho '' = \hat M_{r_C}^\dagger \rho' \hat M_{r_C}.$

"""


# ╔═╡ d6d6e78b-15cb-4c6a-af30-8c42b35dc1b9
md"""
Finally, the system **decoheres** according to decay operators

$M_{\text{decay, } \hat J} = \sqrt{\kappa_J \Delta t} \hat J$

and

$M_{\text{null, } \hat J} = \sqrt{\mathbb{1} - \kappa_J \Delta t \hat J^\dagger \hat J}$

satisfying

$M_{\text{decay, } \hat J}^\dagger M_{\text{decay, } \hat J} + M_{\text{null, } \hat J}^\dagger M_{\text{null, } \hat J} = \mathbb{1}.$

The $\hat J$ are Lindblad operators. The state updates as

$\rho''' = M_{\text{decay, } \hat J} \rho'' M_{\text{decay, } \hat J}^\dagger + M_{\text{null, } \hat J} \rho'' M_{\text{null, } \hat J}^\dagger.$

"""

# ╔═╡ 6315198d-4962-4686-b436-80a3525e0dca
md"""
The final state after increment $\Delta t$ is obtained by renormalizing by the trace:

$\rho(t + \Delta t) = \frac{\rho'''}{\text{Tr}( \rho''')}.$

"""

# ╔═╡ d1d5feec-227e-49cc-ba17-df67974218e7
md" ##### Example"

# ╔═╡ 2b27ab1e-f488-4565-a5b9-35eff5e37a56
md" ### Rouchon approximation "

# ╔═╡ addbae24-128e-426c-bbf6-d430e381f7e3
md""" 
The Rouchon method was presented in [1] as a numerical integration method that provides an efficient alternative to the stochastic master equation (SME). Unlike the SME method, Rouchon's method preserves positivity of the conditioned quantum state by expanding to second order in the Wiener increment $\Delta W_r$.

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
# ╟─d4e7b79b-eb8c-48b3-9a94-650481fa8cb7
# ╠═b6c73eea-d17c-49d4-a091-be10c8134718
# ╟─4fc42181-a1d4-46d9-be76-7277b9a2d6a5
# ╟─343161e9-e975-45b5-a444-b9f9f0b30cfb
# ╟─d5af31e4-ba93-4f3f-9c55-4f2ce600d797
# ╟─97e2560f-3dd9-45a0-85f5-aa9e12800292
# ╟─d6d6e78b-15cb-4c6a-af30-8c42b35dc1b9
# ╟─6315198d-4962-4686-b436-80a3525e0dca
# ╟─d1d5feec-227e-49cc-ba17-df67974218e7
# ╟─2b27ab1e-f488-4565-a5b9-35eff5e37a56
# ╟─addbae24-128e-426c-bbf6-d430e381f7e3
# ╟─676e029d-1141-41b2-bc81-5d2e4500c684
# ╟─e7ffdfb4-622e-4248-a71b-907ea1cbb8a5
