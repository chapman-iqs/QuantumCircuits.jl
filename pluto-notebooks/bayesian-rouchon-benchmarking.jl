### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 80ffb52c-a0dc-43d9-8282-e8acb51df4e0
begin
	cd("/Users/sachagreenfield/Desktop/GitHub/QuantumCircuits.jl")
	import Pkg
	Pkg.activate(".")
	using Random
	using Statistics
	using Distributions
	using PyPlot
	using QuantumCircuits
end

# ╔═╡ 595fc42b-943e-4a29-841a-cd36c90a2b55
md"""
##### Quantum basis
"""

# ╔═╡ 483e648c-2ac3-46f7-b49b-a3109deec27d
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

# ╔═╡ bef8c648-af4c-46c1-87b4-29f744f9aaf3
md"""
##### API

The two methods have the same API:
```
rouchon(T, ρ0, H, J, C; <keyword arguments>)

Arguments:

T :: tuple (ti,tf)
ρ0 :: initial density matrix
H :: time-dependent Hamiltonian
J :: deterministic collapse operators
C :: stochastic collapse operators

Keyword Arguments:

dt :: time step; default dt=1e-4
dy :: record; default dy=[], i.e. simulation generates record time series.
            record should be input in the shape [dy_1,...,dy_Nc] given
            length(C)=Nc collapse operators. Records should have shape
            dy_c[t] indexing the mth trajectory at time T[n]
fn : ρ → Any :: eval function (e.g. to return expectation values instead of density matrices)

Returns: (ts, ρs, dy)

ts :: list of simulation times
ρs :: fn(ρ) at each simulation time
dy :: input OR simulated record, depending on value of keyword argument dy

```
"""

# ╔═╡ 2059a835-925d-430d-aa0c-283a03f6ba2d
md"""
In most (all?) implementations of quantum circuit dynamics, inputs will take the following form:

* `H = H0` or `H = H(t)` : time-in/dependent Hamiltonian;

* `J = [J0, ..., Jn, √(1-η_0)C0, ... , √(1-η_m)Cm]` : list of $n$ dissipation operators $\{J_0, ..., J_n\}$ and $m$ fluctuation-dissipation operators $\{\sqrt{1-\eta_0} C_0, ..., \sqrt{1-\eta_m} C_m \}$ with $n, m \geq 0$;

* `C = [√η C0, ..., √Cm]` : list of $m$ fluctuation operators $\{\sqrt{\eta_0} C_0, ..., \sqrt{\eta_m} C_m \}$.
"""

# ╔═╡ 476f0b08-7d21-4010-b114-130e6dfbbae0
md"""
#### Example
###### Parameters and operators
"""

# ╔═╡ 639931ed-0520-499c-ba1c-90b04e854ebd
begin
	ΩR  = 2π # Rabi frequency
	τm = 3.0 # Measurement collapse timescale
	Γm = 1/(2τm) # Measurement dephasing rate
	η = 0.3
	
	H = ΩR*σx/2
	J = [√((1-η)*Γm)*σz] 
	C = [√(Γm*η)*σz]
	
	T = (0,4τm) # simulation duration
	ρ0 = dm(spinup(q))
	dt = 5e-4 # integration time-step
end

# ╔═╡ 5f98cd8f-2f26-4fd3-b7fb-a75ee4469875
md"""
## Systematic comparison of rouchon + bayesian
"""

# ╔═╡ a504e2af-58ab-4819-af4c-fdcecc0ad0b4
md"""
In order to systematically compare `rouchon` and `bayesian`, we'll first generate 100 simulated records using `bayesian`:
"""

# ╔═╡ 6cf47c1d-7c11-4bd9-a74d-73c6d314f691
N = 20

# ╔═╡ 25512701-e9d4-430e-8dd2-11716285d49b
ensb = ensemble(bayesian, T, ρ0, H, J, C; dt=dt, N=N)

# ╔═╡ 14391721-a7be-4055-a6ca-7664b95ff23d
md"""
then use those same records to reconstruct trajectories using rouchon:
"""

# ╔═╡ 34029d43-1946-49fa-b7b6-ab3e06406f72
ensr = ensemble(rouchon, T, ρ0, H, J, C; dt=dt, N=N)

# ╔═╡ 1aa035ea-5ced-459f-8fd9-82e37e592d12
md" #### RMS error"

# ╔═╡ f8854877-3cd7-4b4f-a534-a7269b153462
md"""
Find the root-mean-square (RMS) error between the trajectories generated from coarse-grained records at scale dt and the original trajectories.
"""

# ╔═╡ 2da3a708-9bec-4207-889f-3632a3fa940b
nlist = prepend!(collect(2:2:100), [1])

# ╔═╡ 3db7293a-ef8f-4ec5-a894-3805514c2210
dydt=[[1,2,3],[5,6,7]]

# ╔═╡ 69f145d7-758b-479a-801a-483e1f823074
dydt./dt

# ╔═╡ e7ffdfb4-622e-4248-a71b-907ea1cbb8a5
md"""

### References

[1] P. Rouchon and J. F. Ralph, Efficient Quantum Filtering for Quantum Feedback Control, Phys. Rev. A 91, 012118 (2015).

"""

# ╔═╡ 389364fb-cb2d-4a48-ba05-f7a222adba91
min(3,5)

# ╔═╡ 4070b2da-508d-467b-bacb-e3c4dacbe641


# ╔═╡ 2b7aa36e-e64a-4839-875e-2c472763cb80
md" ## Utilities "

# ╔═╡ 7baa7a79-090e-4318-982f-6c1982f82a58
function rms(ser1::Array, ser2::Array)
	l = min(length(ser1), length(ser2))
	sqrt(sum((ser1[1:l] - ser2[1:l]).^2))/l
end


# ╔═╡ 087145d6-81d4-44ef-b2a6-58c5126081ee
expects = ρ -> collect(real(expect(ρ, s)) for s in [σx,σy,σz,ρ]) # ρ -> [<x>,<y>,<z>,<ρ>]

# ╔═╡ 078c6b8e-4b74-46dd-aebd-77acf392c2c9
function blochs(sol)
	(tt, ρt, _) = sol

	# Get Bloch components
	evs0 = expects.(ρt);
	xx,yy,zz,ρρ = [map(x -> x[i], evs0) for i in 1:4];

	(collect(tt), xx, yy, zz, ρρ)
	
end

# ╔═╡ a9102fda-a602-45e2-ad7a-6fe192a7972a
purity(x,y,z) = 0.5*(1 + x^2 + y^2 + z^2)

# ╔═╡ 9c95461e-4243-447e-b412-7814ca18da43
R(x,y,z) = sqrt(x^2 + y^2 + z^2)

# ╔═╡ 81df100b-0f36-42ee-a014-e9e1f7fb9cee
test = [1,2,3,4,5]

# ╔═╡ a9b6e622-f9ff-4e05-8c6b-5d05dcdcb876
test[1]

# ╔═╡ 8eb9d84d-e658-406e-8a2a-8a28cfee9b04
function subseries(rec, T, dt; scale=2)
	ts = collect(range(first(T), last(T), step=dt))
	tts = subselect(real(coarse_grain(ts; n=scale)); n=scale)
	(ti, tf) = (tts[1], tts[end])
	dtt = dt * scale
	
	subrec = subselect(real(coarse_grain(rec; n=scale)); n=scale)
	
	(tts, subrec)
	
end

# ╔═╡ 9b2492db-94b4-470f-89f7-370ac8b5db26
function rms(solve, sols, T, ρ0, H, J, C; n=1, dt=1e-4, kwargs...)
	
	"""
	Arguments
	
	solve :: Function -- bayesian or rouchon
	dt :: time-step of reference solutions
	T :: time duration of reference solutions
	sols :: reference solutions
	n :: scale of coarse-graining
	
	Returns
	
	mean(rms_list) :: average rms between reference trajectories and corresponding
	                  coarse-grained generated trajectories
	
	"""
	
	(tt, ρlist, recs) = sols
	recs = collect(map(rec -> collect(rec[1]) , recs))
	
	rms_list = []
	
	for i in 1:length(recs)
		
		# look at one solution at a time
		Rs = recs[i]
		sol = (tt, ρlist[i], Rs)

		# first, get bloch trajectories for each of the ρ in ρs, and coarse-grain
		(ts,xs,ys,zs,rs) = map(ser -> subseries(ser, T, dt; scale=n)[2], blochs(sol))
		
		# coarse-grain the record
		(_, Rsc) = subseries(Rs, T, dt; scale=n)
		
		# generate a trajectory from coarse-grained record, and get blochs
		sol_coarse = solve(T, ρ0, H, J, C; dt=dt*n, dydt=[Rsc], kwargs...)
		(_,xc,yc,zc,rc) = blochs(sol_coarse)
		
		# calculate rms error for each bloch coordinate, then take total
		push!(rms_list, sqrt(rms(xs, xc)^2 + rms(ys, yc)^2 + rms(zs, zc)^2))
		
	end

	# return average rms error
	mean(rms_list)


end

# ╔═╡ 98be356a-7afa-4387-a9ae-57542e808b1b
begin
	rms_bayes = []
	rms_rouchon = []
	times_bayes = []
	times_rouchon = []
	
	for n in nlist
		
		rms_b = @timed rms(bayesian, ensb, T, ρ0, H, J, C; n=n, dt=dt)
		rms_r = @timed rms(rouchon, ensr, T, ρ0, H, J, C; n=n, dt=dt)
		push!(rms_bayes, rms_b.value)
		push!(rms_rouchon, rms_r.value)
		push!(times_bayes, rms_b.time)
		push!(times_rouchon, rms_r.time)
		
	end
end

# ╔═╡ b13e4f3d-df37-4a75-8b00-2d0f768bd18e
# function subseries(rec, T, dt; scale=2)
# 	smooth_rec = real(coarse_grain(rec; n=scale))
	
# 	dtt = dt * scale
# 	tts = collect(range(first(T), last(T), step=dtt))
# 	subrec = subselect(smooth_rec; n=scale)
	
# 	if length(tts) > length(subrec)
# 		push!(subrec,subrec[end])
		
# 	elseif length(tts) < length(subrec)
# 		deleteat!(subrec,length(subrec))
		
# 	end
	
# 	(tts, subrec)
	
# end

# ╔═╡ 3e2b91c7-97cb-4a59-96d5-01081f9d7a6f
ylabel

# ╔═╡ 37ac2cb1-c957-4182-be1f-29521aafbaa3
# colorscheme
begin
	colorscheme = "Paired"
	cmap = plt.matplotlib.cm.get_cmap(colorscheme)
	colors=collect(map(x -> cmap(x), 0:11))
	md" `colorscheme`"
end

# ╔═╡ 6163347e-7044-4ab7-8683-0a7186d9ac62
let
	close("all")
	p = plot(nlist, rms_bayes, color=colors[2], marker="o", label="bayesian")
	plot(nlist, rms_rouchon, color=colors[4], marker="o", label="rouchon")
	ax = gca()
	
    xlabel("n")
    ylabel("rms error")
    title(string("average rms error, w/ reference ", L"$dt = $", dt*1000, " ns"))
    legend()
    gcf()
	
end

# ╔═╡ 7aa70d12-915f-4eb6-8c1a-d4410abcc1a4
let
	close("all")
	p = plot(nlist, times_bayes, color=colors[2], marker="o", label="bayesian")
	plot(nlist, times_rouchon, color=colors[4], marker="o", label="rouchon")
	ax = gca()
	
    xlabel("n")
    ylabel("time (s)")
    title(string("reconstruction time, w/ reference ", L"$dt = $", dt*1000, " ns"))
    legend()
    gcf()
	
end

# ╔═╡ 7672f11f-c5f2-4f5e-8ffc-3f7b45b4a322
# Plotting
function plot_solution(sol; plot_title="Rabi Oscillation")
	
	close("all")
	
	(tt, ρt, dys) = sol
    
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

# ╔═╡ 08a6190e-aa1a-485f-bf30-e059e56048e5
# Plotting
function plot_blochs((tt, blochs); plot_title="Rabi Oscillation")
	close("all")
	
	(xx,yy,zz) = blochs
	pp = map(i -> purity(xx[i], yy[i], zz[i]), 1:length(tt))
    
    # Plot Bloch components vs. time
    
    p = plot(tt, xx, color=colors[2], label=L"$x$")
    plot(tt, yy, color=colors[4],label=L"$y$")
    ax = gca()
    ax.set_ylim([-1.1,1.1]) 
    plot(tt, zz, color=colors[6], label=L"$z$")
	plot(tt, pp, color=colors[8], label=L"Tr $\rho^2$")
    xlabel(L"$t$")
    ylabel("Bloch coordinates")
    title(plot_title)
    legend()
    gcf()
end

# ╔═╡ 4fd5a499-61e5-42ac-9d9d-38f461f59616
function record_histograms(records...; plot_title="record histogram", labels=[]::Array, density=false)
	close("all")
	
	μσs = []
	hist_colors = []
	hist_labels = []
	
	for i in 1:length(records)
		label = i > length(labels) ? i : labels[i]
		dys = records[i]
		
		# get mean and std dev for real part
		(μ, σ) = params(fit(Normal, dys))
		push!(μσs, map(p -> round(p, digits=4), (μ, σ)))

		# make histogram		
		n, bins, patches = hist(dys, 50, density=density, 					 									facecolor=colors[2i], alpha=1, label=label)
		push!(hist_colors, colors[2i])
		

	end
	
	# write down (μ, σ) pairs as text boxes
	μσ_strings = map(μσ -> string("(μ, σ) = (", μσ[1], ", ", μσ[2], ")\n"), μσs)
	ax = gca()
	for i in 1:length(μσ_strings)
		str = μσ_strings[i]

		ax.text(0.05, 1 - 0.05i, str, transform=ax.transAxes, fontsize=10,
			verticalalignment="top", color=hist_colors[i])
		
	end
	
	
	legend()
	xlabel("value (arbitrary units)")
	ylabel(density ? "relative frequency" : "frequency")
	title("record histograms")
	gcf()
		
end

# ╔═╡ 82f2e0a1-11ee-4aa6-a63f-e6ca170a8116
# Plotting
function plot_timeseries(tt::Array, series...; plot_title="time series", xlab=L"$t$", ylab="arbitrary units", labels=[]::Array, colorpairs=false, kwargs...)
	close("all")
	ser_colors(i) = colorpairs ? colors[i] : colors[2i]
	
	label(i) = i > length(labels) ? i : labels[i]
    
    # Plot records vs. time
	ser = series[1]
	p = plot(tt, ser, color=ser_colors(1), label=label(1), kwargs...)
	ax = gca()
	
	if length(series) > 1
		
		for i in 2:length(series)
			ser = series[i]
			plot(tt, ser, color=ser_colors(i),label=label(i), kwargs...) 
		end
		
	end
	
    xlabel(xlab)
    ylabel(ylab)
    title(plot_title)
    legend()
    gcf()
	
end

# ╔═╡ 5156e99d-7583-40fc-aa4f-73df81cbcedc
# Plotting
function plot_timeseries(ttseries...; plot_title="time series", xlab=L"$t$", ylab="arbitrary units", labels=[]::Array, colorpairs=false)
	close("all")
	ser_colors(i) = colorpairs ? colors[i] : colors[2i]
	
	label(i) = i > length(labels) ? i : labels[i]
    
    # Plot records vs. time
	(tt, ser) = ttseries[1]
	p = plot(tt, ser, color=ser_colors(1), label=label(1))
	ax = gca()
	
	if length(ttseries) > 1
		
		for i in 2:length(ttseries)
			(tt, ser) = ttseries[i]
			plot(tt, ser, color=ser_colors(i),label=label(i)) 
		end
		
	end
	
    xlabel(xlab)
    ylabel(ylab)
    title(plot_title)
    legend()
    gcf()
	
end

# ╔═╡ e0df4349-7035-4c41-995c-0c4d78069636
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

# ╔═╡ b3e6866b-6ab7-4f33-be5d-6e28da5f5fe7
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

# ╔═╡ b6577418-9cf1-4232-9fdb-fb6ebf6efd22
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

# ╔═╡ 093383f1-d7a5-48c6-9927-4e42d158cc3d
green(text; title="Note") = Markdown.MD(Markdown.Admonition("correct", title, [text]))

# ╔═╡ 5ed525a9-c4b6-44d6-8ebf-6c86190a8992
red(text; title="Note") = Markdown.MD(Markdown.Admonition("danger", title, [text]))

# ╔═╡ b966b51c-e9cd-4696-a3b1-15ecb9d3808e
tan(text; title="Note") = Markdown.MD(Markdown.Admonition("warning", title, [text]))

# ╔═╡ f08bf969-d03f-47f9-9f5a-cddd18f8b109
blue(text; title="Note") = Markdown.MD(Markdown.Admonition("note", title, [text]))

# ╔═╡ d7c1fd28-5780-4efd-b085-4f1f797a3f09
hint(text; title="Hint") = Markdown.MD(Markdown.Admonition("hint", title, [text]))

# ╔═╡ 59cfb67b-443c-431f-bade-d3f188f08e84


# ╔═╡ Cell order:
# ╠═80ffb52c-a0dc-43d9-8282-e8acb51df4e0
# ╟─595fc42b-943e-4a29-841a-cd36c90a2b55
# ╟─483e648c-2ac3-46f7-b49b-a3109deec27d
# ╟─5d30c7c1-4c53-40d2-8e73-b58db5697cfa
# ╟─d4e7b79b-eb8c-48b3-9a94-650481fa8cb7
# ╟─bef8c648-af4c-46c1-87b4-29f744f9aaf3
# ╟─2059a835-925d-430d-aa0c-283a03f6ba2d
# ╟─476f0b08-7d21-4010-b114-130e6dfbbae0
# ╠═639931ed-0520-499c-ba1c-90b04e854ebd
# ╟─5f98cd8f-2f26-4fd3-b7fb-a75ee4469875
# ╟─a504e2af-58ab-4819-af4c-fdcecc0ad0b4
# ╠═6cf47c1d-7c11-4bd9-a74d-73c6d314f691
# ╠═25512701-e9d4-430e-8dd2-11716285d49b
# ╟─14391721-a7be-4055-a6ca-7664b95ff23d
# ╠═34029d43-1946-49fa-b7b6-ab3e06406f72
# ╟─1aa035ea-5ced-459f-8fd9-82e37e592d12
# ╟─f8854877-3cd7-4b4f-a534-a7269b153462
# ╠═2da3a708-9bec-4207-889f-3632a3fa940b
# ╠═98be356a-7afa-4387-a9ae-57542e808b1b
# ╠═6163347e-7044-4ab7-8683-0a7186d9ac62
# ╠═7aa70d12-915f-4eb6-8c1a-d4410abcc1a4
# ╠═9b2492db-94b4-470f-89f7-370ac8b5db26
# ╠═3db7293a-ef8f-4ec5-a894-3805514c2210
# ╠═69f145d7-758b-479a-801a-483e1f823074
# ╟─e7ffdfb4-622e-4248-a71b-907ea1cbb8a5
# ╠═389364fb-cb2d-4a48-ba05-f7a222adba91
# ╠═4070b2da-508d-467b-bacb-e3c4dacbe641
# ╟─2b7aa36e-e64a-4839-875e-2c472763cb80
# ╠═7baa7a79-090e-4318-982f-6c1982f82a58
# ╠═078c6b8e-4b74-46dd-aebd-77acf392c2c9
# ╟─087145d6-81d4-44ef-b2a6-58c5126081ee
# ╟─a9102fda-a602-45e2-ad7a-6fe192a7972a
# ╟─9c95461e-4243-447e-b412-7814ca18da43
# ╟─7672f11f-c5f2-4f5e-8ffc-3f7b45b4a322
# ╟─08a6190e-aa1a-485f-bf30-e059e56048e5
# ╟─4fd5a499-61e5-42ac-9d9d-38f461f59616
# ╠═81df100b-0f36-42ee-a014-e9e1f7fb9cee
# ╠═a9b6e622-f9ff-4e05-8c6b-5d05dcdcb876
# ╠═8eb9d84d-e658-406e-8a2a-8a28cfee9b04
# ╠═b13e4f3d-df37-4a75-8b00-2d0f768bd18e
# ╠═82f2e0a1-11ee-4aa6-a63f-e6ca170a8116
# ╠═5156e99d-7583-40fc-aa4f-73df81cbcedc
# ╠═3e2b91c7-97cb-4a59-96d5-01081f9d7a6f
# ╟─e0df4349-7035-4c41-995c-0c4d78069636
# ╟─b3e6866b-6ab7-4f33-be5d-6e28da5f5fe7
# ╟─b6577418-9cf1-4232-9fdb-fb6ebf6efd22
# ╟─37ac2cb1-c957-4182-be1f-29521aafbaa3
# ╟─093383f1-d7a5-48c6-9927-4e42d158cc3d
# ╟─5ed525a9-c4b6-44d6-8ebf-6c86190a8992
# ╟─b966b51c-e9cd-4696-a3b1-15ecb9d3808e
# ╟─f08bf969-d03f-47f9-9f5a-cddd18f8b109
# ╟─d7c1fd28-5780-4efd-b085-4f1f797a3f09
# ╠═59cfb67b-443c-431f-bade-d3f188f08e84
