### A Pluto.jl notebook ###
# v0.14.8

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
	using Random
	using Statistics
	using PyPlot
	using QuantumCircuits 	# note: for some reason, QuantumCircuits has to be 
							# "used" last
end

# ╔═╡ 377a3336-20bd-4baa-a033-af8bbc8668a8
md"""
# Reduced qubit evolution

In this interactive notebook, we'll look at reduced qubit evolution as a basic example to gain intuition for Cavity Quantum Electrodynamics (CQED).
"""

# ╔═╡ 7f1176b8-f6d3-4fd1-a23e-31710dcfff10
md"""
## Problem setup

We begin by defining a two-level Hilbert space for the system. `QuantumCircuits.jl` uses `QuantumOptics.jl` as its backend: `SpinBasis`, `sigmax`, `identityoperator` and so forth are `QuantumOptics.jl` functions.
"""

# ╔═╡ 982565bc-b137-4c10-8e13-5fe37be86823
md" ##### Qubit Hilbert space"

# ╔═╡ 4284173a-be05-4b58-a8d9-7189301344fd
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
end

# ╔═╡ 85206c80-af15-4e1b-81e9-04bcb1083063
md" ##### Description of the evolution"

# ╔═╡ 629669ac-006a-4774-ae3e-67787bd0fa0a
md"""
We work in the rotating frame, so that $\omega_q = 0$ effectively. Thus the Hamiltonian is

$H_0 = \frac{\Omega_R}{2} \sigma_y.$
"""

# ╔═╡ bbbc92cd-8daa-4ed8-8bf0-d10461af4edd
md"""
The qubit is weakly and continuously monitored in its informational basis, $|+z \rangle$ (excited) and $|-z \rangle$ (ground). This ''weak measurement'' manifests in two primary ways:

* The information *collected* from the measurement creates **measurement backaction**, characterized by measurement collapse timescale $\tau$. This is represented via the dissipation operator $J_0 = \sqrt{\Gamma} \sigma_z$.

* The information *not collected* leads to **dephasing** at rate $\Gamma = 1/(2\tau)$. This is represented via the measurement collapse operator $C_0 = \sqrt{\Gamma \eta}\sigma_z$.

The balance of these two effects is determined by the signal collection efficiency $\eta$, where $\eta = 1$ means all information is collected, and $\eta = 0$ means no information is collected.
"""

# ╔═╡ 5b48435b-da1b-4547-ac7c-72c95e1d07ff
md"""
##### System parameters

The default parameters of our qubit system are chosen as
"""

# ╔═╡ 99f2b351-877a-43c8-a508-d5688f92ae0c
begin
	ΩR0  = 2π # Rabi frequency
	τ0 = 3.0 # Measurement collapse timescale
	Γ0 = 1/(2τ0) # Measurement dephasing rate
	η0 = 0.3
end

# ╔═╡ 259ed22f-b0c0-4fbe-9c6b-97a8730b6b14
H0 = ΩR0*σy/2

# ╔═╡ b290f188-47f9-4af1-a5f0-11d99c29c1e3
begin
	J0 = √Γ0*σz
	C0 = √(Γ0*η0)*σz;
end

# ╔═╡ 1aa85cc1-32d2-4fa9-a29d-b042d4d4b94a
md"""
## Quantum trajectories

A *quantum trajectory* is a particular realization of the qubit evolution, based on the measurement record. The Hamiltonian evolution of the qubit is interleaved with Bayesian measurement updates based on information acquired from measurement. The resulting trajectory is a noisy perturbation of the ensemble-average behavior.

In the following code, we use `QuantumCircuits.jl`'s `bayesian` function to simulate trajectories for the reduced qubit system.

!!! warning "Edit"
	Include description of what a quantum trajectory is and the Bayesian filtering method for simulating one.

"""

# ╔═╡ 139565f9-0219-4939-a470-4252397dd6a0
begin
	ρ0 = dm(spindown(q)) # initial state
	dt = 1e-2  # integration time-step
	
	Random.seed!(1)
	sol1 = @time bayesian((0, 4τ0), ρ0, H0, [J0], [C0]; dt=dt)
end

# ╔═╡ 62170241-f050-4634-ad6a-8842ac96096c
md" This is one of infinitely many possible trajectories for the system. Let's generate an ensemble of 500 trajectories using `ensemble`:"

# ╔═╡ d1544d60-7de6-4f54-bffb-9ea0c3ddddc1
begin
	Random.seed!(1)
	sol_ens = @time ensemble(bayesian, (0,4τ0), ρ0, H0, [J0], [C0]; dt=dt, N=500)
end

# ╔═╡ 4d42ca32-5348-48b5-8a24-eafa5b303690
md" Let's look at 50 of these trajectories, with our original highlighted: "

# ╔═╡ ee6b82a5-52b0-4467-b73a-d6534de57202
md" While there is some randomness to the trajectories, there is clearly an average behavior that is consistent with all of them. By averaging the Bloch coordinates of all 500 trajectories, for every time step, we can obtain the average qubit evolution: "

# ╔═╡ 0b2ba4be-0213-4b3b-9943-291381c2b4fc
md" η0 = $η0"

# ╔═╡ 0cb6b898-347c-4baf-988f-d1dfac1d6557
md" This is equivalent to setting $\eta = 0$, and reproduces the master equation solution. Intuitively, when no signal is collected, our ''best estimate'' of the quantum state is the Hamiltonian evolution plus additional decoherence due to the uncollected measurement. You can gain intuition for this effect in the next section."

# ╔═╡ 7961752d-faee-4799-99b1-3bd3d17d1f31
md"""
## Measurement efficiency $\eta$ and measurement collapse timescale $\tau$

Let's gain intuition for how $\eta$ and $\tau$ affect our simulation. Play with the values of $\eta$ and $\tau$ below to see what happens.
"""

# ╔═╡ dff8aec2-7216-49e9-a19d-81b7ff79eb4d
md"""
##### Parameters:

"""

# ╔═╡ 6391d48d-4332-4209-bfe1-bec2e60e30c1
md"""

ΩR = 0
$(@bind ΩR_coeff html"<input type=range min=0.1 max=6 step=0.1 value=2>")
ΩR = 10π

η = 0
$(@bind η html"<input type=range min=0. max=1. step=0.02 value=0.3>")
η = 1 | - - - - - - | τ = 0.5 
$(@bind τ html"<input type=range min=0.5 max=10 step=0.5 value=3.0>")
τ=10.0

"""


# ╔═╡ 8638af64-25f9-4793-bc8b-5d22c3ef70fa
ΩR = ΩR_coeff*π

# ╔═╡ 02ce38f3-f5e1-4fea-9e6e-b887baa548e7
md" ##### Simulate"

# ╔═╡ 5e8b047e-b70f-4478-a8b2-a2e48b4a56af
begin
	# other parameters
	Γ = 1/(2τ)
	
	# operators
	H = ΩR*σy/2
	J = √Γ*σz
	C = √(Γ*η)*σz
	md" `operators defined`"
end

# ╔═╡ 2db9b4c4-0256-4e29-990c-2ca3ee9c7ddf
md"""
ΩR = $(round(ΩR/π, digits=3)) π | - - - - - - | η = $η | - - - - - - | τ = $τ | - - - - - - | Γ = 1/(2τ) = $(round(Γ,digits=3))
"""

# ╔═╡ 09b609f6-a441-4ae2-aaff-66c1b8b6ffae
begin
	Random.seed!(1)
	sol2 = @time bayesian((0, 4τ0), ρ0, H, [J], [C]; dt=dt)
	print()
end

# ╔═╡ 7e62e6c6-aa3c-4350-901c-15d017b8db42
md" ## Utilities "

# ╔═╡ b2235afc-1bc9-4cf4-ae21-6f535eb37e96
expects = ρ -> collect(real(expect(ρ, s)) for s in [σx,σy,σz,ρ]) # ρ -> [<x>,<y>,<z>,<ρ>]

# ╔═╡ 0a2219d8-2fbd-4b03-867f-966c4a595e78
# colorscheme
begin
	colorscheme = "Paired"
	cmap = plt.matplotlib.cm.get_cmap(colorscheme)
	colors=collect(map(x -> cmap(x), 0:11))
	md" `colorscheme`"
end

# ╔═╡ 31371407-5798-49a4-a83a-2c068d622c4b
# Plotting
function plot_solution(sol; plot_title="Rabi Oscillation")
	
	close("all")
    
    tt = sol[1]
    ρt = sol[2]
    
    # Get Bloch components
    evs0 = expects.(ρt);
    xx,yy,zz,ρρ = [map(x -> x[i], evs0) for i in 1:4];
    
    # Plot Bloch components vs. time
    
    p = plot(tt, xx, color=colors[2], label=L"$x$")
    plot(tt, yy, color=colors[4],label=L"$y$")
    ax = gca()
    ax.set_ylim([-1.1,1.1]) 
    plot(tt, zz, color=colors[6], label=L"$z$")
	plot(tt, ρρ, color=colors[8], label=L"Tr $\rho^2$")
    xlabel(L"$t$")
    ylabel("Bloch coordinates")
    title(plot_title)
    legend()
    gcf()
end

# ╔═╡ 7218d350-0b14-4fed-b87e-65a9eb573780
plot_solution(sol1; plot_title="Monitored Rabi Oscillation")

# ╔═╡ cdbd897e-d689-42fa-9977-2b3169b29d08
begin
	title_string = "Monitored Rabi Oscillation; τ = $(τ), η = $(η)"
	plot_solution(sol2; plot_title=title_string)
end

# ╔═╡ 99e90182-1ee7-45ef-8656-cc5a1e01d564
# Plotting
function plot_solutions((sol1,sol2); plot_title="Rabi Oscillation")
    close("all")
    
    tt1 = sol1[1]
    ρt1 = sol1[2]
    tt2 = sol2[1]
    ρt2 = sol2[2]
    
    # Get Bloch components
    evs1 = expects.(ρt1);
    x1,y1,z1,ρ1 = [map(x -> x[i], evs1) for i in 1:4];
    evs2 = expects.(ρt2);
    x2,y2,z2,ρ2 = [map(x -> x[i], evs2) for i in 1:4];
    
    # Plot Bloch components vs. time
	p = plot(tt1, x1, color=colors[2], linewidth=2, label=L"$x_{\eta = 0}$")
    plot(tt1, y1, color=colors[4], linewidth=2, label=L"$y_{\eta = 0}$")
    plot(tt1, z1, color=colors[6], linewidth=2, label=L"$z_{\eta = 0}$")
    plot(tt1, p1, color=colors[8], linewidth=2, label=L"(Tr $\rho^2)_{\eta = 0}$")
    plot(tt2, x2,  color=colors[1], linestyle="dashed", label=L"$x_{avg}$")
    plot(tt2, y2, color=colors[3], linestyle="dashed", label=L"$y_{avg}$")
    plot(tt2, p2, color=colors[7], linewidth=2, linestyle="dashed", label=L"(Tr $\rho^2)_{avg}$")
  
	ax = gca()
	ax.set_ylim([-1.1,1.1]) 
	xlabel(L"$t$")
	ylabel("Bloch coordinates")
	title(plot_title)
	legend()
	gcf()
end

# ╔═╡ 53f32791-a721-4d80-96ce-fae9f533dbaa
function plot_evals((tt1, evals); α=0.1, linewidth=1, labels=false)
    xxs,yys,zzs,ρρs = [map(x -> x[i], evals) for i in 1:4];
    if labels
        plot(tt1, xxs, color=colors[2], alpha=α, linewidth=linewidth, label=L"$x$")
        plot(tt1, yys, color=colors[4], alpha=α, linewidth=linewidth, label=L"$y$")
        plot(tt1, zzs, color=colors[6], alpha=α, linewidth=linewidth, label=L"$z$")
        plot(tt1, ρρs, color=colors[8], alpha=α, linewidth=linewidth, label=L"Tr $ρ^2$")
    else
        plot(tt1, xxs, color=colors[2], alpha=α, linewidth=linewidth)
        plot(tt1, yys, color=colors[4], alpha=α, linewidth=linewidth)
        plot(tt1, zzs, color=colors[6], alpha=α, linewidth=linewidth)
        plot(tt1, ρρs, color=colors[8], alpha=α, linewidth=linewidth)
        
    end

end

# ╔═╡ 7d4f4a02-f3b8-4369-aad3-03ff4bc7ec0b
function plot_ensemble(sol_ens; α=0.1, linewidth=1, labels=false, average=false)
    close("all")
	tt1 = sol_ens[1]
    evs = collect(map(ρs -> expects.(ρs), sol_ens[2]));

    for i in 1:50
        plot_evals((tt1, evs[i]); α=α, labels=labels, linewidth=linewidth)
    end

    if average
        plot_evals((tt1, mean(evs)), α=1, linewidth=1.5, labels=true)
        title_string = "Trajectories w/ ensemble average"
    else
        plot_evals((tt1, evs[1]), α=1, linewidth=1.5, labels=true)
        title_string = "Trajectories"
    end
    
    ax = gca()
    ax.set_ylim([-1.1,1.1]) 
	legend()

    xlabel(L"$t$")
    ylabel("Bloch coordinates")
    title(title_string)
	
	gcf()
end

# ╔═╡ 2b5ae77e-07c0-4750-85dd-138dcf84da1f
 plot_ensemble(sol_ens)

# ╔═╡ 8696fe89-1c19-4852-afe0-2b642f5a9382
 plot_ensemble(sol_ens; average=true)

# ╔═╡ a255f7a5-d617-4444-80a4-4e4e7762b236
green(text; title="Note") = Markdown.MD(Markdown.Admonition("correct", title, [text]))

# ╔═╡ 6e86b5ee-5b38-4245-b555-c208151c0c53
red(text; title="Note") = Markdown.MD(Markdown.Admonition("danger", title, [text]))

# ╔═╡ 4585fc96-cc13-476d-99d6-de839c4ec691
if τ > τ0
	green(md"""As the measurement timescale τ increases, the dephasing rate $\Gamma = 1/(2\tau)$ decreases. Thus the system takes longer to dephase (remains high purity).""", title="τ increased from 3.0")
	
elseif τ < τ0
	red(md"""As the measurement timescale τ decreases, the dephasing rate $\Gamma = 1/(2\tau)$ increases. Thus the system dephases quickly and the purity decreases.""", title="τ decreased from 3.0")
		
end

# ╔═╡ fb32f241-34bb-413b-99d2-c0155b363670
tan(text; title="Note") = Markdown.MD(Markdown.Admonition("warning", title, [text]))

# ╔═╡ e730e11b-6fa5-4121-80eb-69a59d74b852
blue(text; title="Note") = Markdown.MD(Markdown.Admonition("note", title, [text]))

# ╔═╡ 5590bdc6-a943-4614-a24d-335512b0f73f
if η == 0.
	tan(md"""When η = 0, the estimated trajectory is exactly the master equation solution. There is no backaction.""",title="η = 0 (no signal collected)")
	
elseif η == 1
	blue(md"""When η = 1, the whole leaked signal is collected, and the state remains pure.""",title="η = 1 (all signal collected)")
	
elseif η < η0
	red(md"""As the signal collection efficiency η decreases, the backaction strength decreases. The trajectory becomes less noisy.""", title="η decreased from 0.3")

elseif η > η0
	green(md"""As the signal collection efficiency η increases, the backaction strength increases. The trajectory becomes more noisy.""",title="η increased from 0.3")
		
end

# ╔═╡ c8aa6f9a-a2b7-45d7-9e39-344380a3909b
hint(text; title="Hint") = Markdown.MD(Markdown.Admonition("hint", title, [text]))

# ╔═╡ dfcd22dc-5517-4ecb-99d1-2151d8c73005
hint(md"To see the most dramatic effects, try setting one parameter to an extreme value and then change the other parameter.")

# ╔═╡ 6bb8877f-d861-492e-b816-2396fcd6d59c
md"""
!!! warning "Edit"

	We need to use simple master equation API rather than `bayesian` here. I'm not sure this currently exists in `QuantumCircuits.jl`, but should probably be the `jump-no-jump` method. Or does `bayesian` automatically do that in the absence of measurement?
"""

# ╔═╡ Cell order:
# ╠═4c45fe3e-cd69-11eb-20cd-7bfb98c040cf
# ╟─377a3336-20bd-4baa-a033-af8bbc8668a8
# ╠═7f1176b8-f6d3-4fd1-a23e-31710dcfff10
# ╟─982565bc-b137-4c10-8e13-5fe37be86823
# ╠═4284173a-be05-4b58-a8d9-7189301344fd
# ╠═85206c80-af15-4e1b-81e9-04bcb1083063
# ╟─629669ac-006a-4774-ae3e-67787bd0fa0a
# ╠═259ed22f-b0c0-4fbe-9c6b-97a8730b6b14
# ╟─bbbc92cd-8daa-4ed8-8bf0-d10461af4edd
# ╠═b290f188-47f9-4af1-a5f0-11d99c29c1e3
# ╟─5b48435b-da1b-4547-ac7c-72c95e1d07ff
# ╠═99f2b351-877a-43c8-a508-d5688f92ae0c
# ╟─1aa85cc1-32d2-4fa9-a29d-b042d4d4b94a
# ╠═139565f9-0219-4939-a470-4252397dd6a0
# ╠═7218d350-0b14-4fed-b87e-65a9eb573780
# ╟─62170241-f050-4634-ad6a-8842ac96096c
# ╠═d1544d60-7de6-4f54-bffb-9ea0c3ddddc1
# ╟─4d42ca32-5348-48b5-8a24-eafa5b303690
# ╠═2b5ae77e-07c0-4750-85dd-138dcf84da1f
# ╟─ee6b82a5-52b0-4467-b73a-d6534de57202
# ╠═8696fe89-1c19-4852-afe0-2b642f5a9382
# ╟─0b2ba4be-0213-4b3b-9943-291381c2b4fc
# ╟─0cb6b898-347c-4baf-988f-d1dfac1d6557
# ╟─7961752d-faee-4799-99b1-3bd3d17d1f31
# ╟─dff8aec2-7216-49e9-a19d-81b7ff79eb4d
# ╟─6391d48d-4332-4209-bfe1-bec2e60e30c1
# ╟─8638af64-25f9-4793-bc8b-5d22c3ef70fa
# ╟─2db9b4c4-0256-4e29-990c-2ca3ee9c7ddf
# ╟─5590bdc6-a943-4614-a24d-335512b0f73f
# ╟─4585fc96-cc13-476d-99d6-de839c4ec691
# ╟─02ce38f3-f5e1-4fea-9e6e-b887baa548e7
# ╠═09b609f6-a441-4ae2-aaff-66c1b8b6ffae
# ╠═cdbd897e-d689-42fa-9977-2b3169b29d08
# ╟─dfcd22dc-5517-4ecb-99d1-2151d8c73005
# ╟─5e8b047e-b70f-4478-a8b2-a2e48b4a56af
# ╟─7e62e6c6-aa3c-4350-901c-15d017b8db42
# ╠═b2235afc-1bc9-4cf4-ae21-6f535eb37e96
# ╠═31371407-5798-49a4-a83a-2c068d622c4b
# ╠═99e90182-1ee7-45ef-8656-cc5a1e01d564
# ╠═53f32791-a721-4d80-96ce-fae9f533dbaa
# ╟─7d4f4a02-f3b8-4369-aad3-03ff4bc7ec0b
# ╠═0a2219d8-2fbd-4b03-867f-966c4a595e78
# ╠═a255f7a5-d617-4444-80a4-4e4e7762b236
# ╠═6e86b5ee-5b38-4245-b555-c208151c0c53
# ╠═fb32f241-34bb-413b-99d2-c0155b363670
# ╠═e730e11b-6fa5-4121-80eb-69a59d74b852
# ╠═c8aa6f9a-a2b7-45d7-9e39-344380a3909b
# ╟─6bb8877f-d861-492e-b816-2396fcd6d59c
