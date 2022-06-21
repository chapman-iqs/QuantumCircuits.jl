### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

# ╔═╡ 75236fca-cfac-11eb-2ffe-c394a1c504cf
begin
	path = pwd()
	using PlutoUI
end

# ╔═╡ 59173935-7989-4ca3-a0f7-5fe26cde1127
Markdown.parse("[Table of Contents 📔](./open?path=$(string(path, "/notebooks/table-of-contents.jl")))")

# ╔═╡ 5b8ed51f-0f75-494e-b3a3-745e0f17f8ed
TableOfContents(title="Notebook utilities")

# ╔═╡ 30061cc8-c84f-43b1-830d-a6dc2d9fc540
md"""
# Plotting functions
"""

# ╔═╡ 904ad4eb-20a9-41b7-8862-784aab3f04b5
md"""
## Single-qubit Hilbert space
"""

# ╔═╡ 5ac39d9e-f299-49f5-9a3f-b73ce8d8e336
md"""
### Ensembles
"""

# ╔═╡ bbf44e72-911a-41ce-954c-07eee87270c7
md"""
### Single trajectories
"""

# ╔═╡ 85f91804-5293-4caf-a72e-47543dab260c
# Plotting
function plot_timeseries(series...; plot_title="time series", xlabel="t (μs)", ylabel="arbitrary units", labels=[], colors = palette(:lightrainbow), kwargs...)

	label(i) = i > length(labels) ? i : labels[i]

	p = plot(size=(600,300))

	# Plot records vs. time
	for (i, (tt, ser)) in enumerate(series)
		plot!(tt, ser, color=colors[i], label=label(i), xlabel=xlabel, ylabel=ylabel, title=plot_title, legend=:outerright, kwargs...)
	end
	
    p
	
end

# ╔═╡ 7a1205e7-f848-4e58-b739-dd322dfc75f2
# Plotting
function plot_timeseries(tt::Vector{Timescale}, series...; plot_title="time series", xlabel=L"$t$", ylabel="arbitrary units", labels=[], colors = palette(:lightrainbow), kwargs...)

	ser_colors(i) = colorpairs ? colors[i] : colors[2i]
	label(i) = i > length(labels) ? i : labels[i]

	p = plot(size=(600,300))

	# Plot records vs. time
	for (i, ser) in enumerate(series)
		plot!(tt, ser, color=ser_colors(i), label=label(i), xlabel=xlabel, ylabel=ylabel, title=plot_title, legend=:outerright, kwargs)
	end
	
    p
	
end

# ╔═╡ 564b3cd6-800f-4f81-9a06-f419bdc2c299
# Plotting
function series_histogram(tt, series, index; title="time series", xlabel="t (μs)", ylabel="arbitrary units", labels=[], color=:blue)

	t = tt[index]
	ser = series[1]
	p = plot([t, t], [-1, 1], linestyle=:dash, color=:black, size=(600,300), linewidth=1.5)

	# Plot records vs. time
	histdata = []
	for (i, ser) in enumerate(series)
		plot!(tt, ser, color=color, xlabel=xlabel, ylabel=ylabel, title=title, legend=:none, alpha=0.5, linewidth=0.7)
		push!(histdata, ser[index])
	end

	h =	histogram(histdata, bins=(-1):0.1:1, normalize=false, color=color, alpha=0.9, orientation=:h, legend=:none, xlabel="occurences", xlims=[0, round(length(series)/2)])
		
	l = @layout [timeseries{0.7w} histogram{0.3w}]
	return plot(p, h, layout=l, link=:y)
	
end

# ╔═╡ 5b6d4ce0-c7a1-446a-a087-48a7acd673e1
# Plotting
function plot_records(series; plot_title="time series", xlabel="t (μs)", ylabel="arbitrary units", labels=[], colors = palette(:lightrainbow), histograms=false)

	label(i) = i > length(labels) ? i : labels[i]

	p = plot(size=(600,300))

	# Plot records vs. time
	for (i, (tt, ser)) in enumerate(series)
		plot!(tt, ser, color=colors[i], xlabel=xlabel, ylabel=ylabel, title=plot_title, legend=:none)
	end

	if !histograms
	    return p
		
	else
		
		h =	histogram(title="record histograms", xlabel="probability", orientation=:h)
		for (i, (tt, ser)) in enumerate(series)
			(μ, σ) = params(fit(Normal, ser))
			histogram!(ser, bins = :scott, normalize=true, color=colors[i], alpha=0.65, label=string(labels[i], ", σ = ", round(σ, digits=3)), orientation=:h, legend=:bottomright)
		end

		l = @layout [timeseries{0.5w} histogram{0.5w}]
		return plot(p, h, layout=l, link=:y)
		
	end
	
end

# ╔═╡ 5b5ed1f4-be2c-4f95-8ac0-023f27471f1b
md"""
## Two-qubit Hilbert space
"""

# ╔═╡ ee5c6288-d3b3-454c-af4a-bcdd0f5e6d87
function bell_plot(exps::Vector{Vector{Float64}}, t::Vector{Float64})

	basis = bell_basis
	colors = palette(:rainbow)
	labels = bell_basis_labels
	title = "Bell states"
	
	pl = plot(size=(600,300), legendfontsize=12, titlefontsize=12, legend=:outerright, title=title)
	
	for l in 1:length(basis)
		label = labels[l]
		color = colors[l]
		exp = exps[l]
		plot!(t, exps[l], color=color, label=label, legend=:outerright, xlabel="t (μs)", ylims=[0,1])
	end

	pl

end

# ╔═╡ bb34d9c6-6802-4129-9b3a-352364bef48b
md"""
## Colors
"""

# ╔═╡ f8bf56e6-e572-4581-a2b7-654232937148
begin
	# from where?
	colors1q = palette(:tab10)
	colors_system = palette(:rainbow)
	colors_filter = palette(:lightrainbow)	

	# 4-measurement-backaction.jl 
	colors1 = palette(:rainbow)
	colors2 = palette(:lightrainbow)
	colorsq1 = palette(:tab20)[1:2:20]
	colorsq2 = palette(:tab20)[2:2:20]

	mixed_colors = let
		cols = []
		for i in 1:length(colors1)
			push!(cols, colors1[i])
			push!(cols, colors2[i])
		end
		cols
	end
end

# ╔═╡ 05b3bc9f-81c7-4470-a887-05b000d1b084
function single_qubit_plots(sol::Solution)
	colors = colors1q 
	
	# calculate expectation values --------------------------------------------
	q1_basis = [σx1, σy1, σz1]
	q2_basis = [σx2, σy2, σz2]
	
	exps1 = map(op -> expectations(sol, op), q1_basis)
	exps2 = map(op -> expectations(sol, op), q2_basis)

	qlabels = ["x", "y", "z"]

	p1s = 0.5 * (1 .+ exps1[1].^2 .+ exps1[2].^2 .+ exps1[3].^2)
	p2s = 0.5 * (1 .+ exps2[1].^2 .+ exps2[2].^2 .+ exps2[3].^2)
	
	# plot ----------------------------------------------------------------------
	l = @layout [bloch1{0.5h}; bloch2{0.5h}]

	p1 = plot(sol.t, p1s, color=colors[4], label=L"Tr(\rho_1^2)", linestyle=:dash, title="single-qubit states")

	for l in 1:3
		label = qlabels[l]
		color = colors[l]
		exp = exps1[l]
		plot!(sol.t, exp, color=color, label=label, legend=:outerright, ylims=[-1,1])
	end

	p2 = plot(sol.t, p2s, color=colors[4], label=L"Tr(\rho_2^2)", linestyle=:dash, title="")

	for l in 1:3
		label = qlabels[l]
		color = colors[l]
		exp = exps2[l]
		plot!(sol.t, exp, color=color, label=label, legend=:outerright, xlabel="t (μs)", ylims=[-1,1])
	end
	
	plot(p1, p2, layout = l, link=:y, size=(600,300), legendfontsize=8, titlefontsize=12, legend=:outerright)
	
end

# ╔═╡ 3ddf7a52-ac73-4fa7-8ee5-30357a066a26
function bell_plot(sol::Solution)

	basis = bell_basis
	colors = colors_system
	labels = bell_basis_labels
	title = "Bell states"

	exps = map(op -> expectations(sol, dm(op)), basis)
	
	pl = plot(size=(600,300), legendfontsize=12, titlefontsize=12, legend=:outerright, title=title)
	
	for l in 1:length(basis)
		label = labels[l]
		color = colors[l]
		exp = exps[l]
		plot!(sol.t, exps[l], color=color, label=label, legend=:outerright, xlabel="t (μs)", ylims=[0,1])
	end

	pl

end

# ╔═╡ 1e8698cc-f79f-421a-9bdb-f17c6d3bfc84
function bell_plot(sys::Solution, fil::Solution)

	basis = bell_basis
	labels = bell_basis_labels
	title = "Bell states"

	exps_sys = map(op -> expectations(sys, dm(op)), basis)
	exps_fil = map(op -> expectations(fil, dm(op)), basis)
	
	pl = plot(size=(600,300), legendfontsize=12, titlefontsize=12, legend=:outerright, title=title, xlabel="t (μs)",)
	
	for l in 1:length(basis)
		plot!(sys.t, exps_sys[l], color=colors_system[l], label=labels[l], legend=:outerright, ylims=[0,1])
		plot!(fil.t, exps_fil[l], color=colors_filter[l], label=:none, linestyle=:dash)
	end

	pl

end

# ╔═╡ cad89ec7-3916-42c3-a978-48dc9437ea48
md"""
## Animation
"""

# ╔═╡ 4b31e5a9-74b0-4bdc-9b8a-9e754f4e4d3d
function animate_plot(tt, plot, args...; step=100, fps=15)
	anim = @animate for i ∈ range(1, length(tt), step=step)
			plot(i, args...) end
	gif(anim, fps = fps)
end	

# ╔═╡ f7f96c78-27bf-4c87-8e68-4a05a9e3765f
md"""
# Structs
"""

# ╔═╡ a84d8645-10dd-47d3-b2de-7a5e1054c021
begin
	mutable struct traj
	  t::Vector{Float64}
	  x::Vector{Float64}
	  y::Vector{Float64}
	  z::Vector{Float64}
	  p::Vector{Float64}
	  r
	end
	
	function traj(sol::Solution)
	  t, ρ, r = (sol.t, sol.ρ, sol.r)
	  x, y, z = map(op -> expectations(sol, op), [σx, σy, σz])
	  p = (typeof(ρ[1]) <: Ket) ? [1.0 for el in ρ] : real(expect.(ρ, ρ))
	  traj(t, x, y, z, p, r)
	end
end

# ╔═╡ d7982457-b1eb-40f7-a7ab-a9053baf4e08
md"""
# Other
"""

# ╔═╡ 4783f63a-08c4-42a3-9e74-c295e3a5fd5b
md"""
## Array manipulation
"""

# ╔═╡ 65ac7c43-bad3-4e7c-9b2c-2d15f5687307
getclosest(array, val) = argmin(abs.(val .- array))

# ╔═╡ a371616b-f652-4c30-9f89-ad658f52a89f
md"""
## Labeling
"""

# ╔═╡ 3eb3583f-7684-4706-a05f-89acdcc4652d
begin
	qbasis = [σx, σy, σz]
	qlabels = ["x", "y", "z"]
end

# ╔═╡ 0baa1d76-f0fe-44b9-96c8-5dc6befef504
function bloch_plots(sols::Vector{Solution}; alpha=0.1, N=50)
	colors = palette(:tab10) 
	
	# calculate expectation values --------------------------------------------
	t = sols[1].t
	xs, ys, zs = [], [], []
	
	for sol in sols
		x, y, z = map(op -> expectations(sol, op), qbasis)
		for (list, traj) in zip([xs, ys, zs], [x, y, z])
			push!(list, traj)
		end
		
	end

	
	
	# plot ----------------------------------------------------------------------
	function bloch(os; color=colors1q[1], xlabel="", ylabel="")
		
		po = plot(ylims = [-1,1], xlabel=xlabel, ylabel=ylabel)
		
		for o in os[1:min(N, 50)]
			plot!(t, o, alpha=alpha, label=:none, color=color)
		end
		
		oavg = [mean([os[i][j] for i in 1:N]) for j in 1:length(t)]
		plot!(t, oavg, alpha=1, color=color, label="average", linewidth=3)

		po
		
	end

	
	l = @layout [xplot{0.33h}; yplot{0.33h}; zplot{0.33h}]
	
	px = bloch(xs, color=colors[1], ylabel="x")
	py = bloch(ys, color=colors[2], ylabel="y")
	pz = bloch(zs, color=colors[3], ylabel="z", xlabel="t (μs)")
	
	plot(px, py, pz, layout = l, link=:y, size=(800,500), legendfontsize=8, titlefontsize=12, legend=:outerright)
	
end

# ╔═╡ 5d56aa00-e4c7-48fd-8e73-99501f539956
# from 3-qubit-ensembles.jl
function bloch_plots(sols::Vector{Solution}, sol_η0::Solution; alpha=0.1, N=50)
	colors = colors1q 
	
	# calculate expectation values --------------------------------------------
	t = sols[1].t
	xs, ys, zs = [], [], []
	
	for sol in sols
		x, y, z = map(op -> expectations(sol, op), qbasis)
		for (list, traj) in zip([xs, ys, zs], [x, y, z])
			push!(list, traj)
		end
		
	end

	# η = 0 solution
	xη0, yη0, zη0 = map(op -> expectations(sol_η0, op), qbasis)

	
	
	# plot ----------------------------------------------------------------------
	function bloch(os, oη0; color=colors1q[1], xlabel="", ylabel="")
		
		po = plot(ylims = [-1,1], xlabel=xlabel, ylabel=ylabel)
		
		for o in os[1:min(N, 50)]
			plot!(t, o, alpha=alpha, label=:none, color=color)
		end
		
		oavg = [mean([os[i][j] for i in 1:N]) for j in 1:length(t)]
		plot!(t, oavg, alpha=1, color=color, label="average", linewidth=3)
		plot!(t, oη0, alpha=1, color=:black, label="η = 0", linewidth=2)

		po
		
	end

	
	l = @layout [xplot{0.33h}; yplot{0.33h}; zplot{0.33h}]
	
	px = bloch(xs, xη0, color=colors[1], ylabel="x")
	py = bloch(ys, yη0, color=colors[2], ylabel="y")
	pz = bloch(zs, zη0, color=colors[3], ylabel="z", xlabel="t (μs)")
	
	plot(px, py, pz, layout = l, link=:y, size=(800,500), legendfontsize=8, titlefontsize=12, legend=:outerright)
	
end

# ╔═╡ 54670d16-0bf0-4cbb-9ff1-2947c817fc2e
function qubit_plot(sol::Solution; record=false, title="", legendpos=:bottomleft)

	basis = qbasis

	t = sol.t
	exps = map(op -> expectations(sol, op), basis)
	r = record ? sol.r[1] : []

	return qubit_plot((t, exps, r); title=title, legendpos=legendpos)
end

# ╔═╡ 65be6ef5-0e70-4c74-b73d-c0e8590fee95
function qubit_plot((t, exps, r); title="", legendpos=:bottomleft)

	record = (r != [])

	basis = qbasis
	colors = colorsq1
	labels = qlabels

	p = 0.5 .* (1 .+ exps[1].^2 .+ exps[2].^2 .+ exps[3].^2)

	pl = plot(size=(600,300), legendfontsize=10, titlefontsize=12, legend=:outerright, ylabel="bloch coordinates", xlabel = record ? "" : "t (μs)", linewidth=1.5, title=title)

	for l in 1:length(basis)
		label = labels[l]
		color = colors[l]
		exp = exps[l]
		plot!(t, exps[l], color=color, label=label, legend=legendpos, ylims=[-1,1], linewidth=1.5)
	end

	plot!(t, p, color=colors[4], label=L"Tr(\rho^2)")

	if !record
		return pl
	else
		l = @layout [blochs{0.6h}; record{0.4h}]
		pr = plot(t, r, color=colors1[1], xlabel="t (μs)", ylabel="record", label=:none, legend=legendpos, title="", linewidth=0.8)
		return plot(pl, pr, layout = l, link=:y)
	end
end

# ╔═╡ 2530efde-1ac9-4e6b-b4e6-003ef0f3d9d3
function qubit_plot(sol1::Solution, sol2::Solution; record=false, title="", color1=colorsq1, color2=colorsq2, l1="", l2="")

	basis = qbasis

	t1 = sol1.t
	exps1 = map(op -> expectations(sol1, op), basis)
	r1 = record ? sol1.r[1] : []
	
	t2 = sol2.t
	exps2 = map(op -> expectations(sol2, op), basis)
	r2 = record ? sol2.r[1] : []
	
	return qubit_plot((t1, exps1, r1), (t2, exps2, r2); title=title, color1=color1, color2=color2, l1=l1, l2=l2)
	
end

# ╔═╡ 2f33d62f-c0c1-4572-afa6-7a78dad8076e
function qubit_plot((t1, exps1, r1), (t2, exps2, r2); title="", color1=colorsq1, color2=colorsq2, l1="", l2="", ylims=[-1,1])

	basis = qbasis
	labels = qlabels

	record = ((r1 != []) || (r2 != []))

	p1 = 0.5 .* (1 .+ exps1[1].^2 .+ exps1[2].^2 .+ exps1[3].^2)
	p2 = 0.5 .* (1 .+ exps2[1].^2 .+ exps2[2].^2 .+ exps2[3].^2)

	pl = plot(size=(600,300), legendfontsize=10, titlefontsize=12, legend=:outerright, ylabel="bloch coordinates", xlabel = record ? "" : "t (μs)", linewidth=1.2, title=title)

	for l in 1:length(basis)
		plot!(t1, exps1[l], color=color1[l], label=labels[l], legend=:outerright, ylims=ylims, linewidth=1.2)
	end
	plot!(t1, p1, color=color1[4], label=L"Tr(\rho^2)", linewidth=1.2)
	
	for l in 1:length(basis)
		plot!(t2, exps2[l], color=color2[l], label=:none, legend=:outerright, ylims=ylims, linestyle=:dash, linewidth=2)
	end
	plot!(t2, p2, color=color2[4], label=:none, linestyle=:dash, title=string(l1, " ---, ", l2, " - - -"), linewidth=2)


	if !record
		return pl
	else
		l = @layout [blochs{0.6h}; record{0.4h}]
		pr = plot()
		if r1 != []	
			plot!(t1, r1, color=mixed_colors[1], xlabel="t (μs)", ylabel="record", label=:none, legend=:outerright, title="")
		end
		if r2 != []
			plot!(t2, r2, color=mixed_colors[2], xlabel="t (μs)", ylabel="record", label=:none, legend=:outerright, title="")
		end
		return plot(pl, pr, layout = l, link=:both)
	end
end

# ╔═╡ 5687f4e3-a57e-4885-bc81-95633193a17e
md"""
# Utilities
"""

# ╔═╡ 3c3d59df-1eab-4b1b-874c-bd5457da8d6f
function notebook_path(folder, index)
	notebooks = readdir(folder)
	string(path, "/", folder, "/", notebooks[index])
end

# ╔═╡ a5cbd8d7-c099-44c9-bd48-298246062c97
path

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.34"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.0-DEV.1208"
manifest_format = "2.0"
project_hash = "e766545f1b4ef968b5991d6a702e32de0b70adeb"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.0+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.73.0+4"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.9.1+2"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.24.0+2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2020.7.22"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.17+2"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "13468f237353112a01b2d6b32f3d0f80219944aa"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8979e9802b4ac3d58c503a20f2824ad67f9074dd"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.34"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "3.1.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.41.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "16.2.1+1"
"""

# ╔═╡ Cell order:
# ╠═59173935-7989-4ca3-a0f7-5fe26cde1127
# ╟─5b8ed51f-0f75-494e-b3a3-745e0f17f8ed
# ╟─30061cc8-c84f-43b1-830d-a6dc2d9fc540
# ╟─904ad4eb-20a9-41b7-8862-784aab3f04b5
# ╟─5ac39d9e-f299-49f5-9a3f-b73ce8d8e336
# ╠═0baa1d76-f0fe-44b9-96c8-5dc6befef504
# ╠═5d56aa00-e4c7-48fd-8e73-99501f539956
# ╟─bbf44e72-911a-41ce-954c-07eee87270c7
# ╠═54670d16-0bf0-4cbb-9ff1-2947c817fc2e
# ╠═65be6ef5-0e70-4c74-b73d-c0e8590fee95
# ╟─2530efde-1ac9-4e6b-b4e6-003ef0f3d9d3
# ╠═2f33d62f-c0c1-4572-afa6-7a78dad8076e
# ╠═85f91804-5293-4caf-a72e-47543dab260c
# ╠═7a1205e7-f848-4e58-b739-dd322dfc75f2
# ╠═564b3cd6-800f-4f81-9a06-f419bdc2c299
# ╠═5b6d4ce0-c7a1-446a-a087-48a7acd673e1
# ╠═5b5ed1f4-be2c-4f95-8ac0-023f27471f1b
# ╠═05b3bc9f-81c7-4470-a887-05b000d1b084
# ╠═ee5c6288-d3b3-454c-af4a-bcdd0f5e6d87
# ╠═3ddf7a52-ac73-4fa7-8ee5-30357a066a26
# ╠═1e8698cc-f79f-421a-9bdb-f17c6d3bfc84
# ╟─bb34d9c6-6802-4129-9b3a-352364bef48b
# ╠═f8bf56e6-e572-4581-a2b7-654232937148
# ╟─cad89ec7-3916-42c3-a978-48dc9437ea48
# ╠═4b31e5a9-74b0-4bdc-9b8a-9e754f4e4d3d
# ╟─f7f96c78-27bf-4c87-8e68-4a05a9e3765f
# ╠═a84d8645-10dd-47d3-b2de-7a5e1054c021
# ╟─d7982457-b1eb-40f7-a7ab-a9053baf4e08
# ╟─4783f63a-08c4-42a3-9e74-c295e3a5fd5b
# ╠═65ac7c43-bad3-4e7c-9b2c-2d15f5687307
# ╠═a371616b-f652-4c30-9f89-ad658f52a89f
# ╠═3eb3583f-7684-4706-a05f-89acdcc4652d
# ╟─5687f4e3-a57e-4885-bc81-95633193a17e
# ╟─3c3d59df-1eab-4b1b-874c-bd5457da8d6f
# ╠═75236fca-cfac-11eb-2ffe-c394a1c504cf
# ╟─a5cbd8d7-c099-44c9-bd48-298246062c97
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
