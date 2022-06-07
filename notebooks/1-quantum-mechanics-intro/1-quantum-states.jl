### A Pluto.jl notebook ###
# v0.17.4

using Markdown
using InteractiveUtils

# ╔═╡ c659a7cf-7922-498c-9188-648caa106225
begin

	directory_name = "QC-notebooks"
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
	using Distributions
	using QuantumCircuits
	using Plots
	using StatsPlots

	include("notebooks/table-of-contents.jl")
	include("notebooks/resources.jl")

	include("utilities/single-qubit-operators.jl")
	include("utilities/plotting.jl")
	
	md" # Packages and julia files"
end

# ╔═╡ 377a3336-20bd-4baa-a033-af8bbc8668a8
md"""
In this interactive notebook, we'll introduce the idea of the quantum state.
"""

# ╔═╡ d5e8592d-2887-4be1-acaa-88700f0c5d23
mdp(table_of_contents📔)

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents(title="Quantum states")

# ╔═╡ 3aeab324-35b2-4d53-bffe-bd40775ac908
md"""
# What / why is a "quantum state"?
"""

# ╔═╡ f10fcb92-6703-4a6a-a4aa-cc725aa5add0
md"""
The quantum state lies at the heart of quantum mechanics. To understand what it is, it helps to consider why it is necessary.

##### Classical analogy

The closest classical analogue to the quantum state is the coordinate. A coordinate is a vector containing all relevant information about the system. By "relevant", we mean all information necessary to predict system behaviors of interest, given an evolution equation. For example, a classical particle moving in one dimension has coordinate $(x, p)$, and this coordinate determines all future behavior once a potential $V(x)$ is provided. Systems with higher dimensions may be described by higher dimensional vectors $(\vec x, \vec p)$, and may also include other variables such as thermodynamic variables of entropy $S$, temperature $T$ and so forth, so that the system state is $(\vec x, \vec p, S, T)$.

"""

# ╔═╡ 9340160d-e579-4115-9e1d-a84719d584dd
md"""
More generally, a classical state may be described by a probability distribution over coordinates. For example, if we think the position and momentum are $(x_0, p_0)$ up to some uncertainties $\Delta x$, $\Delta p$, we might describe the state as a Gaussian probability density

$P(x, p) = \frac{e^{-\frac12 \big((x - x_0)/\Delta x\big)^2} e^{-\frac12 \big((p - p_0)/\Delta p\big)^2}}{2\pi (\Delta x \Delta p) }.$

This is appropriate whenever we aren't quite sure what the coordinate is. Then we can evolve the system under its Hamiltonian by evolving the whole probability density according to Liouville's equations. Wikipedia has a [nice animation of this](https://en.wikipedia.org/wiki/Liouville%27s_theorem_(Hamiltonian)):

"""

# ╔═╡ 5c4e9eb3-dcaa-4fd1-8dfd-99509a485a87
md"""
This may sound like overkill, but it is an experimental fact: all measurements are are imprecise. All numbers are reported with uncertainties, often in the form of, e.g. $x_0 \pm \Delta x$, and this always implies a distribution of the type $P(x)$ (even when the distribution is unspecified)."""

# ╔═╡ b91667a8-bea3-4df3-9aec-1ee4314428a1
md"""
What is this distribution? There are two complementary ways of understanding probability distributions, both of which are necessary to get the complete picture. 
1. *Bayesian interpretation*: The probability distribution represents our *state of knowledge* about the system.

2. *Frequentist interpretation*: The probability distribution represents *outcomes of measurements of identical experiments*.

When we think of classical probability distributions, the Bayesian interpretation seems most natural: the distribution is centered at what we *think* the actual value is, and the variance of the distribution our uncertainty about our measurement. However, we can also imagine preparing a particle at position $(6.0, 10.0) \pm 0.5$ cm on a slide, using a grid of cm precision. Then we put the slide under an optical microscope that has μm precision, and repeat the position measurement precisely. We could repeat this process 1000 times. So long as the grid is accurate, we might expect to see our precision measurements Gaussian distributed with standard deviation $0.5$ cm about $(6.0, 10.0)$ cm. So we see that the frequentist interpretation is also possible classically.
"""

# ╔═╡ b50811e6-b070-45a7-b4e5-bb1efa5bc8db
md"""
##### Quantum state
The quantum version of the classical coordinate is the pure state vector (PSV), denoted $\ket \psi$ ("ket psi"). Just as a classical state can be a coordinate OR a probability distribution over coordinates, a quantum state can be a PSV (pure state) OR a probability distribution over PSVs (mixed state). Like $\vec X$, $\ket \psi$ can be evolved under a Hamiltonian to give the state at a later time.

There are two (related) key differences between $\vec X$ and $\ket \psi$:

1. *PSVs can have probabilistic outcomes, coordinates cannot:* If my classical state is a coordinate $\vec X$, I am guaranteed the outcome $\vec X$ when I take a measurement. However, even if I *know* my quantum state is $\ket \psi$, I can still have probabilistic outcomes.

2. *PSVs can be in so-called "superposition" with complex amplitudes*, that only when squared give the probability. This superposition is itself another PSV.

Upon deeper examination we will find these two statements are equivalent, and basically amount to the statement *"The space of quantum states is larger than can be measured in a single context."* To understand what all this means, we will start with  the simplest possible quantum system: the quantum bit or "qubit". 
"""

# ╔═╡ e668d360-b7c1-4e53-ad9d-92eed2850d49
md"""
# Qubit states
"""

# ╔═╡ 781c9e75-188d-4a21-9f49-d85f507df0d6
md"""
## Pure states (coordinates)
"""

# ╔═╡ 24657e95-a146-4add-918f-afe7dbeb7f7f
md"""
A classical bit is a system with two possible coordinate values: $0$ or $1$. 

A quantum bit is a quantum system with two possible PSVs: $\ket 0$ or $\ket 1$. Because of key difference #2, these PSVs can be in superposition of the form

$\ket \psi = \alpha \ket 0 + \beta \ket 1,$

where $\alpha, \beta \in \mathbb C$ are complex numbers satisfying $|\alpha|^2 + |\beta|^2 = 1$. Because of this property, quantum physicists like to write this state as

$\ket \psi = \cos \frac{\theta}2 \ket 0 + e^{i \phi} \sin \frac{\theta}2 \ket 1.$

We can consider $\theta, \phi$ as the polar and azimuthal angles the surface of a sphere. This leads to nice visual depictions of qubit states:
"""

# ╔═╡ 0f641c3c-f5a5-4eef-85e8-e7eb7f482f4a
md"""
Clearly, in addition to $\ket 0$ and $\ket 1$, there are a continuum of states. 

What is the meaning of a state like $\ket + \equiv \frac1{\sqrt{2}} (\ket 0 + \ket 1)$? We can understand this through the frequentist and Bayesian views:

1. *frequentist* : For this state, $\alpha = \beta = \frac1{\sqrt{2}}$. One of the rules of quantum mechanics is that squaring these gives the probabilities of measuring each state, i.e. $P_0 = |\alpha|^2 = \frac12$ and $P_1 = |\beta|^2 = \frac12$. If we prepare 1000 states $\ket +$, we will measure it to be $\ket 0$ approximately 500 times and $\ket 1$ approximately 500 times.

2. *Bayesian* : $\alpha$ and $\beta$ store our *knowledge* of the quantum state. This includes the fact that we expect to measure $\ket 0$ and $\ket 1$ with equal probability, for a given measurement. However, $\alpha$ and $\beta$ are complex amplitudes, so they also capture something that this frequentist interpretation doesn't: something about the prior evolution that isn't reflected in the $P_0$, $P_1$ statistics.
"""

# ╔═╡ 432f17b9-d976-433b-8af3-6b90405a3e23
md"""
### Thought experiment
"""

# ╔═╡ 8a1158aa-7f78-4f77-8cac-6ab7b3e4b5a2
begin
	Ry(θ) = exp(im * DenseOperator(σy) * θ/2)
	sqrtY = Ry(π/2)
	ketp = (g + e)/√2
end

# ╔═╡ bf6738b3-0b80-4646-835f-80c63e40af9e
md"""
To understand why this is fundamentally different from classical bits, we can consider the following thought experiment: we have an unlimited "source" of states called $\ket +$. We don't know anything about them, other than that they are guaranteed to be identical. """

# ╔═╡ 3248edd6-01f9-479f-94ef-e37a898c958a
ketp

# ╔═╡ 39a0dbfc-d0f5-4158-abd3-2ea3da2d3a8d
md"""
We are also handed a "mystery" operation called $\sqrt{\hat Y}$. 
"""

# ╔═╡ b1219e4b-5ed9-4a59-846d-8251286c7576
sqrtY

# ╔═╡ 2492aba8-0aeb-4bae-bfed-8f26561db4b3
md"""
Finally, we are given a measuring device with two outputs: $0$ or $1$: 
"""

# ╔═╡ 475b45cc-3240-4361-8328-5cc39f198427
md"""
We are tasked with finding out as much as possible about the state $\ket +$ and how it behaves under the mystery operation and measurement. Here is one possible procedure:
"""

# ╔═╡ 8c143833-4990-4c23-9854-7ac4986418dc
md"""
###### First,
"""

# ╔═╡ d6325ff2-f95c-4836-bb10-dfa9bd48ecc7
md"We may measure $\ket +$ once.  "

# ╔═╡ 57543852-2189-43c9-b748-c390e8f2139e
md"""
###### To check this guess,
"""

# ╔═╡ 1edbf9d9-3888-45f2-a947-2789e0fe70c2
N = 5000

# ╔═╡ c5824862-04e2-476d-9bcc-f34b9e20ce01
md"To check this guess, we measure it $N more times."

# ╔═╡ a03eb3b3-3fa1-470f-91fb-32ee50b42698
md"""
###### Now,
"""

# ╔═╡ 8ab49226-ca8f-463c-8607-2b6021fcd154
md" we try applying the $\sqrt{\hat Y}$ operation before measurement, for each of 1000 measurements. "

# ╔═╡ c734517e-e628-41c6-aa52-41d8642be6d4
md"""
We find that we get $1$ $100\%$ of the time. We guess that the  $\sqrt{\hat Y}$ operation turns all the $0$s into $1$s, and leaves the $1$s alone.
"""

# ╔═╡ 0af8fa28-ebfc-482e-ac89-fd8375e010cb
md"""
###### To test this conclusion,
"""

# ╔═╡ 617e2efd-d391-4655-bc66-ffa449c748eb
md"""
we measure $N copies of

$\ket +$ and sort into groups of $0$ and $1$ to give us definite states. Then we apply $\sqrt{\hat Y}$ to both groups, and measure:
"""

# ╔═╡ c345d7f4-95e6-4be4-9328-2f768fbb9244
md"""
By now, you should be thinking there's something funny going on. $\sqrt{\hat Y}$ seems to treat $0$ and $1$ differently depending on whether or not we measured before applying  $\sqrt{\hat Y}$. Rather than turning all the $0$s to $1$s as we thought, it seems to "know" that a measurement was taken before applying $\sqrt{\hat Y}$, and turns only half of the 0s to 1s, and vice versa.

If we believe that the particles were definitely $0$ or $1$ before AND after we measured (i.e. when they were " $\ket +$"), we conclude that they must also have a "memory" of what happens that determines the action of $\sqrt{\hat Y}$.
"""

# ╔═╡ b6c0d2d9-b785-4724-917c-69093a409458
md"""
##### Quantum states are sensitive to measurement context
"""

# ╔═╡ c65fac52-d431-45b3-a788-e1a9cb30881f
labels = [L"|0\rangle", L"\frac{1}{\sqrt{2}}(|0\rangle + |1\rangle)", L"\frac{1}{\sqrt{2}}(|0\rangle + i |1\rangle)", L"\frac{1}{\sqrt{2}}(|0\rangle - |1\rangle)", L"\frac{1}{\sqrt{2}}(|0\rangle - i |1\rangle)", L"|1\rangle"]

# ╔═╡ cb3290d2-7253-4ed7-994e-2caa21315e94
counts(outcome, data) = length(data[data .== outcome])

# ╔═╡ d5c8beb5-28d9-4e21-922a-307a87e5b95a
md"""
### Superpositions


"""

# ╔═╡ 4f0b8974-2055-4a16-a65d-40601869ac34
md"""
## Mixed states (distributions)
"""

# ╔═╡ 83ac03db-fc53-40bf-b37b-4f632b562530
md"""
The most general classical bit state is a probabilistic state $(P_0, 1 - P_0)$ where $P_0 \in [0,1]$ is the probability of measuring $0$. Similarly, we can create an analogous distribution out of quantum states.
"""

# ╔═╡ 674c1c5d-88ea-4e28-9a85-1e249b8a0309
function measure(ket)
	P0 = real(expect(dm(ket), dm(g)))
	return rand() < P0 ? (0, g) : (1, e)
end

# ╔═╡ 8e376058-b33b-40e7-8c38-19df32b707fa
measure(ketp)

# ╔═╡ a4262761-ba97-4ca7-9a79-ed0d0f2b5cf5
begin
	exp1 = measure(ketp)
	v1 = exp1[1]
end

# ╔═╡ 48a3112f-c9ec-4165-81ee-c9e5eb24cdf0
md"""
We get outcome $v1. 

We guess we guess that the state $\ket +$ corresponds to measurement outcome $v1.
"""

# ╔═╡ fe9eeae4-9179-43d6-91dd-b5f65e887aed
begin
	data = [measure(ketp)[1] for i in 1:N]
	c0 = counts(0, data)
	c1 = counts(1, data)
end

# ╔═╡ edf4b27e-65a0-49fe-a6d0-5ce9c2662259
md"""
We get $c0 0s and $c1 1s. 

We guess that $\ket +$ means that the person giving us the state prepared a $0$ with $50\%$ probability and $1$ with $50\%$ probability.
"""

# ╔═╡ 6dbc1183-bc7b-42e1-8abd-a098428a3094
let
	outcomes = ["0", "1"]
	
	bar(outcomes, [counts(0, data), counts(1, data)], ylabel = "counts", 
	        title = string(N, " trials"), bar_width = 0.4, legendfontsize=11,tickfontsize=12, size=(300,400), legend=:none)
end

# ╔═╡ 34da5ac2-3682-4e16-a591-9c6faec5f7eb
let
	data = [measure(sqrtY * ketp)[1] for i in 1:N]
	global d0 = counts(0, data)
	global d1 = counts(1, data)
end

# ╔═╡ 9ba3f203-d5a8-4129-86ee-8f193660a70d
let
	outcomes = ["0", "1"]
	
	bar(outcomes, [d0, d1], ylabel = "counts", 
	        title = string(N, " trials"), bar_width = 0.4, legendfontsize=11,tickfontsize=12, size=(300,400), legend=:none)
end

# ╔═╡ d210b742-3087-48bc-9b55-c66d14dfb9ab
begin
	# measure ketp
	measured = [measure(ketp) for i in 1:N]

	# sort
	zero_states = last.(measured[first.(measured) .== 0])
	one_states = last.(measured[first.(measured) .== 1])

	# apply sqrtY and measure again
	data0 = measure.([sqrtY * state for state in zero_states])
	data1 = measure.([sqrtY * state for state in one_states])
end

# ╔═╡ 42e6069a-713e-494a-aea8-639d5c815aa2
let
	outcomes = ["0", "0", "1", "1"]

	labels = [L"\sqrt{\hat{Y}} \cdot 0", L"\sqrt{\hat{Y}} \cdot 1",
				L"\sqrt{\hat{Y}} \cdot 0", L"\sqrt{\hat{Y}} \cdot 1"]

	data01 = vcat(counts(0,first.(data0)), counts(1,first.(data0)), counts(0,first.(data1)), counts(1,first.(data1)))
	
	
	groupedbar(outcomes, data01, group = labels, ylabel = "counts", 
	        title = string(N, " trials"), bar_width = 0.4, legendfontsize=12,tickfontsize=12, size=(600,400), legend=:outerright, legendtitle="state before measurement")
end

# ╔═╡ 72f249cb-d658-43bb-b4fa-acdc008800c1
md"""
# Links
"""

# ╔═╡ e1fa4f7c-fe6d-41b8-a3ad-73b67ab96aac
Liouville_url = "https://upload.wikimedia.org/wikipedia/commons/f/f7/Hamiltonian_flow_classical.gif"

# ╔═╡ 55a3bc91-ad81-4cc9-9736-32c2e27b2baf
Resource(Liouville_url)

# ╔═╡ 7481a13f-be55-4440-a61f-bb65a5a8caca
magic_states = html"<a href='https://earltcampbell.com/research/magic-states/' title='Magic states'>Magic states</a>"

# ╔═╡ 84f6e6a7-b5ce-4d09-9c33-19b23a71b4b3
md"""
We can refer to this memory of "what happened before the measurement" as the measurement "context". This constitutes the key difference between classical and quantum states: **statistics from a single measurement context (the sequence of operations leading up to the measurement) are not sufficient to describe the state.** Or, there is more to quantum states than we can uncover in one type of measurement, even if we repeat the experiment many times.

But wait -- what if there are just hidden degrees of freedom we haven't accounted for, that help the statistics match up? E.g., maybe there is an invisible companion particle carrying the "quantum memory". People have investigated this in great depth. It turns out this is possible only for special states -- for a qubit, these are called "stabilizer states" -- and has been rigorously shown to be impossible in general. The space of quantum mechanics truly is larger than any single measurement context can capture. This feature is called "contextuality" and can be explained in as many different ways as there are physicists talking about it.

See $magic_states for a nice explanation.
"""

# ╔═╡ 133b6939-10b2-4c8e-acf8-5658ca96a0f9
md" # Utilities"

# ╔═╡ 6ccdb308-54eb-44b0-a6fc-819984c3897e
function bloch_vectors(vecs...; mesh=30, ax=false, viewϕ=0, connecting=false, labels=[], size=(400,400))
	bcolors = palette(:seaborn_bright)

	# Wire frame coordinates ---------------------------------------------------

	x(θ, ϕ) = sin(θ) * cos(ϕ + viewϕ)
	y(θ, ϕ) = sin(θ) * sin(ϕ + viewϕ)
	z(θ, ϕ) = cos(θ)

	θs = range(0, 2π, length=mesh)
	ϕs = range(0, π, length=mesh) .+ viewϕ


	# Plot wireframe -----------------------------------------------------------
	wf = plot()
	
	# Longitudes
	for ϕ in ϕs
		plot!([x(θ, ϕ) for θ in θs], [y(θ, ϕ) for θ in θs], [z(θ, ϕ) for θ in θs], label=:none, linecolor="steelblue", linewidth=0.5, linealpha=1, seriestype=path3d, aspect_ratio=1.0, size=size)
			
	end

	# Latitudes
	for θ in θs
		plot!([x(θ, ϕ) for ϕ in ϕs], [y(θ, ϕ) for ϕ in ϕs], [z(θ, ϕ) for ϕ in ϕs], label=:none, linecolor="steelblue", linewidth=0.5, linealpha=1, seriestype=path3d, aspect_ratio=1.0, size=size)
	end


	# Plot reference axes ------------------------------------------------------
	
	colors = [palette(:tab10)[i] for i in 1:3]

	if ax
		plot!([0, cos(viewϕ)], [0, sin(viewϕ)], [0, 0], linecolor=colors[1], linewidth=3.0, label=:none)
		plot!([0, -sin(viewϕ)], [0, cos(viewϕ)], [0, 0], linecolor=colors[2], linewidth=3.0, label=:none)
		plot!([0, 0], [0, 0], [0, 1], linecolor=colors[3], linewidth=3.0, label=:none)
	end


	# Plot Bloch vectors input by user --------------------------------

	for (i, vec) in enumerate(vecs)

		(xvv, yvv, zvv) = vec

		xv = xvv * cos(viewϕ) - yvv * sin(viewϕ)
		yv = xvv .* sin(viewϕ) .+ yvv * cos(viewϕ)
		zv = zvv

		if connecting
			plot!([0, xv], [0, yv], [0, zv], label=:none, linewidth=2, linecolor=bcolors[i])
		end
	
		plot!([xv], [yv], [zv], legend=:outerright, legendfontsize=10, marker=(:circle, 5), markercolor=bcolors[i], label=try labels[i] catch e "" end)
	end

	return wf


end

# ╔═╡ 832f539f-5710-450f-b054-b1c08be0ce94
colors = palette(:tab10)

# ╔═╡ 235cee23-c8af-4de2-a1c3-a2b173156703
expects(ops) = ρ -> collect(real(expect(ρ, s)) for s in vcat(ops, ρ)) # ρ -> [<x>,<y>,<z>,<ρ>]

# ╔═╡ f3ac9b41-b123-4236-8c62-8e2b988463c2
begin
	@userplot MyPlot
	
	@recipe function f(mp::MyPlot; add_marker=false)
		
		x, y = mp.args
		
		linecolor   --> :blue
		seriestype  :=  :path
		markershape --> (add_marker ? :circle : :none)
		legend := :none
		
		@series begin
			x, y
		end
	end
	
end

# ╔═╡ e27bd39c-58b7-4c5c-a677-3fe70f500ee8
xyz(θ, ϕ) = (sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ))

# ╔═╡ fd7c2d87-dd17-442f-9e5f-82801e9d2b30
bloch_vectors(xyz(0,0), xyz(π/2,0), xyz(π/2, π/2), xyz(π/2, π), xyz(π/2, 3π/2), xyz(π,0), ax=true, labels=labels, size=(600,600))

# ╔═╡ 97342f70-4ce4-480b-be81-e8e913539802
θϕ(x, y, z) = (acos(z), atan(y/x))

# ╔═╡ 28d8df47-07cf-4a0d-a447-f894d021b2bc
begin
	mutable struct traj
		t::Vector{Float64}
		x::Vector{Float64}
		y::Vector{Float64}
		z::Vector{Float64}
		p::Vector{Float64}
		r
	end
	
	function traj(t, ρ, r)
		x, y, z = [real(expect(σi, ρ)) for σi in (σx, σy, σz)]
		p = real(expect.(ρ, ρ))
		traj(t, x, y, z, p, r)
	end
	
	function traj(sol::QuantumCircuits.solution)
		t, ρ, r = (sol.t, sol.ρ, sol.r)
		x, y, z = [real(expect(σi, ρ)) for σi in (σx, σy, σz)]
		p = real(expect.(ρ, ρ))
		traj(t, x, y, z, p, r)
	end
	
end

# ╔═╡ dc3e2dcb-5bf1-492d-8337-f366dcf0170b
φdict = Dict("0" => 0, 
				"π/8" => π/8, 
				"π/4" => π/4,
				"3π/8" => 3π/8,
				"π/2" => π/2)

# ╔═╡ Cell order:
# ╟─377a3336-20bd-4baa-a033-af8bbc8668a8
# ╟─d5e8592d-2887-4be1-acaa-88700f0c5d23
# ╟─3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─3aeab324-35b2-4d53-bffe-bd40775ac908
# ╟─f10fcb92-6703-4a6a-a4aa-cc725aa5add0
# ╟─9340160d-e579-4115-9e1d-a84719d584dd
# ╠═55a3bc91-ad81-4cc9-9736-32c2e27b2baf
# ╟─5c4e9eb3-dcaa-4fd1-8dfd-99509a485a87
# ╟─b91667a8-bea3-4df3-9aec-1ee4314428a1
# ╟─b50811e6-b070-45a7-b4e5-bb1efa5bc8db
# ╟─e668d360-b7c1-4e53-ad9d-92eed2850d49
# ╟─781c9e75-188d-4a21-9f49-d85f507df0d6
# ╟─24657e95-a146-4add-918f-afe7dbeb7f7f
# ╠═fd7c2d87-dd17-442f-9e5f-82801e9d2b30
# ╟─0f641c3c-f5a5-4eef-85e8-e7eb7f482f4a
# ╟─432f17b9-d976-433b-8af3-6b90405a3e23
# ╠═8a1158aa-7f78-4f77-8cac-6ab7b3e4b5a2
# ╟─bf6738b3-0b80-4646-835f-80c63e40af9e
# ╠═3248edd6-01f9-479f-94ef-e37a898c958a
# ╟─39a0dbfc-d0f5-4158-abd3-2ea3da2d3a8d
# ╠═b1219e4b-5ed9-4a59-846d-8251286c7576
# ╟─2492aba8-0aeb-4bae-bfed-8f26561db4b3
# ╠═8e376058-b33b-40e7-8c38-19df32b707fa
# ╟─475b45cc-3240-4361-8328-5cc39f198427
# ╟─8c143833-4990-4c23-9854-7ac4986418dc
# ╟─d6325ff2-f95c-4836-bb10-dfa9bd48ecc7
# ╠═a4262761-ba97-4ca7-9a79-ed0d0f2b5cf5
# ╟─48a3112f-c9ec-4165-81ee-c9e5eb24cdf0
# ╟─57543852-2189-43c9-b748-c390e8f2139e
# ╟─c5824862-04e2-476d-9bcc-f34b9e20ce01
# ╟─1edbf9d9-3888-45f2-a947-2789e0fe70c2
# ╠═fe9eeae4-9179-43d6-91dd-b5f65e887aed
# ╟─edf4b27e-65a0-49fe-a6d0-5ce9c2662259
# ╟─6dbc1183-bc7b-42e1-8abd-a098428a3094
# ╟─a03eb3b3-3fa1-470f-91fb-32ee50b42698
# ╟─8ab49226-ca8f-463c-8607-2b6021fcd154
# ╠═34da5ac2-3682-4e16-a591-9c6faec5f7eb
# ╟─9ba3f203-d5a8-4129-86ee-8f193660a70d
# ╟─c734517e-e628-41c6-aa52-41d8642be6d4
# ╟─0af8fa28-ebfc-482e-ac89-fd8375e010cb
# ╟─617e2efd-d391-4655-bc66-ffa449c748eb
# ╠═d210b742-3087-48bc-9b55-c66d14dfb9ab
# ╟─42e6069a-713e-494a-aea8-639d5c815aa2
# ╟─c345d7f4-95e6-4be4-9328-2f768fbb9244
# ╟─b6c0d2d9-b785-4724-917c-69093a409458
# ╟─84f6e6a7-b5ce-4d09-9c33-19b23a71b4b3
# ╟─c65fac52-d431-45b3-a788-e1a9cb30881f
# ╟─cb3290d2-7253-4ed7-994e-2caa21315e94
# ╟─d5c8beb5-28d9-4e21-922a-307a87e5b95a
# ╟─4f0b8974-2055-4a16-a65d-40601869ac34
# ╟─83ac03db-fc53-40bf-b37b-4f632b562530
# ╠═674c1c5d-88ea-4e28-9a85-1e249b8a0309
# ╟─72f249cb-d658-43bb-b4fa-acdc008800c1
# ╟─e1fa4f7c-fe6d-41b8-a3ad-73b67ab96aac
# ╟─7481a13f-be55-4440-a61f-bb65a5a8caca
# ╟─133b6939-10b2-4c8e-acf8-5658ca96a0f9
# ╟─6ccdb308-54eb-44b0-a6fc-819984c3897e
# ╠═832f539f-5710-450f-b054-b1c08be0ce94
# ╟─235cee23-c8af-4de2-a1c3-a2b173156703
# ╟─f3ac9b41-b123-4236-8c62-8e2b988463c2
# ╠═e27bd39c-58b7-4c5c-a677-3fe70f500ee8
# ╠═97342f70-4ce4-480b-be81-e8e913539802
# ╠═28d8df47-07cf-4a0d-a447-f894d021b2bc
# ╠═dc3e2dcb-5bf1-492d-8337-f366dcf0170b
# ╠═c659a7cf-7922-498c-9188-648caa106225
