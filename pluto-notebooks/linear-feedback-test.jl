### A Pluto.jl notebook ###
# v0.16.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

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
	using DataFrames
	
	include("plotting.jl")
	
	md" ###### 🔶 Packages and julia files"
end

# ╔═╡ 377a3336-20bd-4baa-a033-af8bbc8668a8
md"""
# Linear feedback stabilization

In this interactive notebook, we'll look at dynamical linear feedback stabilization of quantum trajectories, based on [1].
"""

# ╔═╡ 3edd54c6-4b52-41ff-a707-a6efce05e698
TableOfContents()

# ╔═╡ 8f7c6440-eac8-49d3-bd23-c2545ec16830
md"""
## Measurement backaction principles
"""

# ╔═╡ ff721d9a-c37c-43a1-adee-015221ddd133
md"""
To understand how linear feedback stabilization works, we must first look at how measurement creates backaction in a qubit.
"""

# ╔═╡ 727ced31-7a51-4892-b1d9-eff4fcdc9767
md"""
Weak measurement of a qubit observable $\hat A$ generates a stochastic record

$\tilde r (t) = \braket{\hat A}(t) + \zeta(t)$

where $\braket{\hat A}(t)  \equiv \text{Tr} \big(\rho_q(t) \hat A \big)$ is the expectation value of the observable at time $t$, and 

$\zeta(t) \sim \mathcal{N} \Big(0, \sqrt{\tau_m/ dt} \Big)$

is a zero-mean, Gaussian distributed random variable with variance $\sigma = \tau_m / dt$. Here, $\tau_m$ is the measurement timescale determined by the strength of the measurement.
"""

# ╔═╡ 1aaa9cbe-2f98-41e1-af41-598ba6578818
md"""
In this notebook, we will consider measurement of the qubit state $\ket0$ or $\ket1$, such that $\hat A = e^{i \varphi} \hat \sigma_z$. The phase $\varphi$ is called the "measurement angle" and correpsonds to the quadrature of amplification of the signal. It determines the type of backaction experienced by the qubit.
"""

# ╔═╡ 1ee143ba-2f0d-4eee-a76c-4bf58469e0be
md"""
### Simulation
"""

# ╔═╡ e43e5329-bd96-41ce-a183-1bd206204f65
begin
	# Basis
	q = SpinBasis(1//2)

	# Operators, using convention that |-z> is ground state
	σx = sigmax(q)
	σy = sigmay(q)
	σz = sigmaz(q)
	σp = sigmap(q)
	σm = sigmam(q)
	id = identityoperator(q)
	
	ground = spindown(q)
	excited = spinup(q)

	md"###### 🔶 Qubit Hilbert space operators "
end

# ╔═╡ 875adb1b-22f7-4375-83c5-2a17c9a8e46c
md"""

Rotate view: `ϕv = 0`
$(@bind ϕvIB html"<input type=range min=0 max=32 step=1 value=0>") 
`ϕv = 2π`

"""

# ╔═╡ 3531558a-f9b7-4f6b-9800-884fe0b04712
md" $(@bind animIB_bs CheckBox()) Animate Bloch sphere "

# ╔═╡ faf0b11c-4339-4a82-a23c-9d35eb9d10b4
md" $(@bind animIB CheckBox(default=false)) Animate informational backaction"

# ╔═╡ 911d4aa1-41cd-4eea-8e0f-1af6f8024290
@bind go Button("Rerun simulation")

# ╔═╡ 3ae402b0-476a-4c11-823c-d0016793adab
md"""
Change the measurement angle to see how backaction depends on the measurement quadrature.
"""

# ╔═╡ 71604112-0704-4832-8a80-e886b6155e34
md"""
`φ =` $(@bind φIB Select(["0", "π/8", "π/4", "3π/8", "π/2"]))
"""

# ╔═╡ d3e48ac4-973f-4397-bf06-89150a0c13ed
md"""
Informational backaction occurs when $\varphi = 0$. Phase backaction occurs when $\varphi = \pi/2$. Intermediate angles result in a combination of both types of backaction.

In the following simulation, the state is initialized in the state $\ket{+x}$. By default, the measurement angle is $\varphi = 0$, leading to informational backaction. The measurement angle is currently `φ =` $φIB.

Interact with the following plots and try different parameter settings (below) to gain intuition for informational and phase backaction.
"""

# ╔═╡ 3ff551a9-3a07-4a2f-928e-880b7e3ba1fc
if 	φIB == "0"
	
	md""" **Informational backaction**  The qubit undergoes a biased random walk towards one of the poles. For $\eta = 1$, the walk is confined to the circle $x^2 + z^2 = 1$. Once it reaches a pole, it is pinned there. As an aside, this demonstrates how continuous weak monitoring can be understand as a projective measurement in the limit of measuring for a long time. The statistics of final states will follow that of projective measurement, in the absence of Hamiltonian evolution.

While the initial state creates a bias towards the pole it is closer to, it is possible to stabilize to either pole so long as the initial state is a superposition. Click the button below to see this:
	"""

elseif φIB == "π/2"
	
	md""" **Phase backaction**  The qubit undergoes a unbiased random walk in the $x$ -- $y$ plane.  For $\eta = 1$, the walk is confined to the circle $x^2 + y^2 = 1$. Note that this effect is preserved for a state initialized outside of the $x$ -- $y$ plane; the random walk will then be confined to its original plane: $x^2 + y^2 = 1 - z^2$.

	"""

else
	
	md""" **Informational and phase backaction**  The qubit undergoes a biased random walk towards one of the poles.  For $\eta = 1$, the walk is confined to the sphere $x^2 + y^2 + z^2 = 1$. 

	"""

end

# ╔═╡ 133b6939-10b2-4c8e-acf8-5658ca96a0f9
md" # Utilities"

# ╔═╡ 8e85754f-d66b-477b-8153-b162519edb7c
expects(ops) = ρ -> collect(real(expect(ρ, s)) for s in vcat(ops, ρ)) # ρ -> [<x>,<y>,<z>,<ρ>]

# ╔═╡ fdcfcea3-e2b6-4939-ac6c-eada7421f3dd
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

# ╔═╡ ccbcf668-d948-4ec6-a5f7-39a178d54c29
xyz(θ, ϕ) = (sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ))

# ╔═╡ 2b35485b-155b-4fb1-868d-b431cc867d61
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

# ╔═╡ 618168dc-53dd-4562-8061-67a0b56587aa
φdict = Dict("0" => 0, 
				"π/8" => π/8, 
				"π/4" => π/4,
				"3π/8" => 3π/8,
				"π/2" => π/2)

# ╔═╡ 49065a29-93ac-458d-ac0c-84f2f2151f9d
let
	
	# System parameters ------------------------------------------------------------
	# all times given in μs
	
	# initial state 
	x0, y0, z0 = xyz(π/2, 0)
	ρ0 = DenseOperator(0.5*(id + x0*σx + y0*σy + z0*σz))
	
	# measurement parameters
	τm = 1 					# time
	Γm = 1/(2τm) 			# rate
	η =  1 					# efficiency
	φ = get(φdict, φIB, 0) 	# angle
	
	# simulation timescales
	T = (0, 8τm) # simulation duration
	dt = 0.5e-3  # integration time-step
	
	
	# Kraus operators -------------------------------------------------------------
	
	H = 0 * id
	J = [(σz, ((1-η)*Γm))]
	C = [(exp(im * φ) * σz, τm, η)]
	
	
	
	# Bayesian simulation ---------------------------------------------------------
	go
	sol = bayesian(T, ρ0, H, J, C; dt=dt)
	
	global IB = traj(sol)
	
	md" ###### 🔻 Simulation "
	
end

# ╔═╡ f9eaaf70-4e0f-4503-8a64-2380682354ce
let 
	ϕv = ϕvIB * (π/16)
	tt, xx, yy, zz = (IB.t, IB.x, IB.y, IB.z)
	
	if animIB_bs
		
		anim = @animate for i ∈ range(1, length(tt), step=100)
			blochsphere(xx[1:i], yy[1:i], zz[1:i], linewidth=1., linealpha=0.85, ax=true, viewϕ = ϕv, blochmark=true) end
		
		gif(anim, fps = 15)
		
	else
		blochsphere(xx, yy, zz, linewidth=1., linealpha=0.85, ax=true, viewϕ = ϕv, blochmark=true)
		
	end
	
end

# ╔═╡ 6bc3e01f-b3cb-4a32-ab6b-e5dcc967b07f
let 
	tt, xx, yy, zz = (IB.t, IB.x, IB.y, IB.z)
	
	if animIB
		anim = @animate for i ∈ range(1, length(tt), step=100)
			blochprojections(xx[1:i], yy[1:i], zz[1:i], blochmark=true) end
		gif(anim, fps = 15)
	else
		blochprojections(xx, yy, zz, blochmark=true)
	end
end

# ╔═╡ df97d34b-16a7-49a0-a143-939f18248f48
blochtimeseries(IB.t, IB.x, IB.y, IB.z, title = "Monitored qubit", tf=last(IB.t), ylims=[-1.1,1.5])

# ╔═╡ 62472483-7fae-4adc-976c-9275e9d5ebfc
myplot(IB.t, coarsegrain(IB.r, n=50))

# ╔═╡ Cell order:
# ╠═4c45fe3e-cd69-11eb-20cd-7bfb98c040cf
# ╟─377a3336-20bd-4baa-a033-af8bbc8668a8
# ╟─3edd54c6-4b52-41ff-a707-a6efce05e698
# ╟─8f7c6440-eac8-49d3-bd23-c2545ec16830
# ╟─ff721d9a-c37c-43a1-adee-015221ddd133
# ╟─727ced31-7a51-4892-b1d9-eff4fcdc9767
# ╟─1aaa9cbe-2f98-41e1-af41-598ba6578818
# ╟─1ee143ba-2f0d-4eee-a76c-4bf58469e0be
# ╟─d3e48ac4-973f-4397-bf06-89150a0c13ed
# ╟─e43e5329-bd96-41ce-a183-1bd206204f65
# ╠═49065a29-93ac-458d-ac0c-84f2f2151f9d
# ╟─875adb1b-22f7-4375-83c5-2a17c9a8e46c
# ╟─3531558a-f9b7-4f6b-9800-884fe0b04712
# ╠═f9eaaf70-4e0f-4503-8a64-2380682354ce
# ╟─faf0b11c-4339-4a82-a23c-9d35eb9d10b4
# ╟─6bc3e01f-b3cb-4a32-ab6b-e5dcc967b07f
# ╠═df97d34b-16a7-49a0-a143-939f18248f48
# ╠═62472483-7fae-4adc-976c-9275e9d5ebfc
# ╟─3ff551a9-3a07-4a2f-928e-880b7e3ba1fc
# ╟─911d4aa1-41cd-4eea-8e0f-1af6f8024290
# ╟─3ae402b0-476a-4c11-823c-d0016793adab
# ╟─71604112-0704-4832-8a80-e886b6155e34
# ╟─133b6939-10b2-4c8e-acf8-5658ca96a0f9
# ╟─8e85754f-d66b-477b-8153-b162519edb7c
# ╠═fdcfcea3-e2b6-4939-ac6c-eada7421f3dd
# ╠═ccbcf668-d948-4ec6-a5f7-39a178d54c29
# ╠═2b35485b-155b-4fb1-868d-b431cc867d61
# ╠═618168dc-53dd-4562-8061-67a0b56587aa
