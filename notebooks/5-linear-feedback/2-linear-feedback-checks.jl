### A Pluto.jl notebook ###
# v0.17.4

using Markdown
using InteractiveUtils

# ╔═╡ 1994807c-85a7-4309-8786-c1b5d2dd6b86
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
	using Plots
	using QuantumCircuits

	include("notebooks/table-of-contents.jl")
	include("notebooks/resources.jl")

	include("utilities/single-qubit-operators.jl")
	include("utilities/utilities.jl")
	include("utilities/plotting.jl")

	md" # Packages and julia files"
	
	
end

# ╔═╡ b9bd68f2-df48-4072-8269-1898b7cf1b15
md"""
In this interactive notebook, we test the feedback matrix exponentiation code to compare to bloch equations and verify that `QuantumCircuits.jl` is working as expected.
"""

# ╔═╡ 99a066a7-7013-472e-82d8-3f4c651c3475
mdp(table_of_contents📔)

# ╔═╡ 0fb01de9-882a-46c1-8447-0018146842df
TableOfContents(title="Linear feedback checks")

# ╔═╡ 18c29abc-1f35-4b5a-bb27-a491c02cc98f
md"""
# System description

The Hamiltonian is a Rabi drive modulated by the feedback from the measurement readout:

$\hat H_c = \hbar \big[\Delta_0 + \Delta_1 \tilde{r}(t - T_d) \big] \frac{\hat \sigma_\phi}2$

where

$\hat \sigma_\phi = \cos \phi \hat \sigma_x + \sin \phi \hat \sigma_y,$

and $\hat \sigma_x$ and $\hat \sigma_y$ are the corresponding Pauli operators. Thus there is a base Rabi frequency $\Delta_0$ modulated by a linear feedback term $\Delta_1 \tilde{r} (t-T_d)$.

"""

# ╔═╡ 01f57775-b647-4fea-8e96-0b8c8ceeff05
md" ## Parameters and operators"

# ╔═╡ c5de979e-de3b-4a20-9fc4-649851a311fa
ideal = true

# ╔═╡ eae605ed-f411-4f33-8066-bd8f01fc8a2d
begin
	# all times given in μs
	# initial state
	(x0,y0,z0) = (0., 0.3, 0.91)
	ρ0 = DenseOperator(0.5 * (Iq + x0*σx + y0*σy + z0*σz))
	
	dt = 1e-3  # integration time-step
	td = 0.0 # 200e-3 # delay time
	η = 1.0 # 0.41 # measurement efficiency
	
	θs = 3π/10 # target angle on Bloch sphere
	ϕ = π # fixes plane of oscillations
	σϕ = cos(ϕ)*σx + sin(ϕ)*σy
	
	τm =  0.2 # measurement time
	Γm = 1/(2τm) # measurement rate
	tf = 4.0
	T = (0.0, tf) # simulation duration
	td = 0.0 # time delay for feedback
	
	Rs = 1.0 # 0.64 # radius of target
	
	# feedback drive parameters
	Δ0 = ideal ? -sin(2θs)/(4τm) : -sin(2θs)/(4τm*Rs^2) 
	Δ1 = ideal ? sin(θs)/τm : sin(θs)/(τm*Rs)
		
	T1 = 40 # energy decay time
	T2 = 60 # environmental dephasing time

	# corresponding rates
	Γ1 = 1/(2T1)
	Γ2 = 1/T2

	
end

# ╔═╡ ee541c05-c187-4b43-a803-2255e254efe5
begin
	# Hamiltonian defined at a certain time t, 
	# takes in readout r in some observable with time delay td
	H0(t::Timescale, r::Readout) = (Δ0 + Δ1*r[1]) * σϕ / 2
	J0 = ideal ? [(σz, (1 - η) * Γm)] : [(σz, (1 - η) * Γm), (σm, Γ1), (σz, Γ2)]
	C0 = [(σz, Γm, η)]
end

# ╔═╡ 8aa08bfb-ac91-4e5f-9fb2-dbce02a38b8a
begin
	
	title_str = ideal ? "ideal Rabi oscillation with feedback, general" : "non-ideal Rabi oscillation with feedback, general";
	title_str_hc = ideal ? "ideal Rabi oscillation with feedback, hard-coded" : "non-ideal Rabi oscillation with feedback, hard-coded";
	
end

# ╔═╡ c7f72983-771c-46e5-a1ec-94ef80391dd3
md"""
# Method comparison
"""

# ╔═╡ ea558387-68dc-4dff-ae04-f7a877379510
md" ## Bloch equations (Patti *et al.*)"

# ╔═╡ 6294f81c-ea9d-4500-9751-b45f8a348639
begin
	# set the target coordinates
	ts = range(first(T), last(T), step=dt)
	ztar = [Rs*cos(θs) for i in 1:length(ts)]
	ytar = [Rs*sin(θs) for i in 1:length(ts)]
	xtar = zeros(length(ts))
end

# ╔═╡ a5f1944a-80b7-4db3-9355-44140e05d908
md"""
### Evolution function
"""

# ╔═╡ f2892559-e623-46d4-9ff2-f04fd9253734
function blochevolve(; rs=[])
	# readout distribution
	sim = (length(rs) == 0)
	ts = range(first(T), last(T), step=dt)
	
	# first time step
	(xn, yn, zn) = (x0, y0, z0)	
	(xs, ys, zs) = ([x0], [y0], [z0])

	
	if sim
		dist = Normal(0, sqrt(τm/dt))
		push!(rs, 0) end
		

	for i in 2:length(ts)
		
		# sample from distribution		
		if sim
			r = zn + rand(dist)
			push!(rs, r)
		else
			r = rs[i] end
		
 		# update equations
		pn = cosh(r*dt/τm) + zn*sinh(r*dt/τm)
		xn1 = xn/pn
		yn1 = yn/pn
		zn1 = (zn*cosh(r*dt/τm) + sinh(r*dt/τm))/pn
		Δ = (Δ0 + Δ1*r) # the Rabi-drive depends on the readout (instantaneous here)
		
		yn2 = yn1*cos(dt*Δ) + zn1*sin(dt*Δ)
		zn2 = zn1*cos(dt*Δ) - yn1*sin(dt*Δ)
		
		# if ideal 
		xn = xn1*exp(-dt*(1-η)/(2τm*η))
		yn = yn2*exp(-dt*(1-η)/(2τm*η))
		zn = zn2

		# else 
		# 	xn = xn1*exp(-dt/(2T1) - dt/T2 - dt*(1-η)/(2τm*η))
		# 	yn = yn2*exp(-dt/(2T1) - dt/T2 - dt*(1-η)/(2τm*η))
		# 	zn = zn2*exp(-dt/T1) - (1 - exp(-dt/T1)) 
		# end

		# store values
		push!(xs, xn)
		push!(ys, yn)
		push!(zs, zn)
		

	end
	
	ρs = 0.5*(1 .+ xs.^2 + ys.^2 + zs.^2)
	
	return ts, rs, (xs, ys, zs, ρs)
	
	
end

# ╔═╡ a1358f8d-0622-49ff-bd71-4d4fb5b434f3
md"""
### Simulation
"""

# ╔═╡ 4487690f-349c-461c-a5b1-f9a6be9897bd
begin
	Random.seed!(2)
	(tt, r, (xx, yy, zz, rr)) = blochevolve()
	tt = collect(tt)
end

# ╔═╡ 15613717-87c1-45d8-b194-1c7b1085d4f7
md" ## Matrix exponentiation"

# ╔═╡ 840beeea-a731-4a38-b36b-e2ab848147d9
begin
	H = length(methods(H0)) > 0 ? H0 : t -> H0
	J = map(j -> length(methods(j)) > 0 ? j : t -> j, J0) # same, for each element of J0
	C = []
	for (c, Γ, η) in C0
		push!(C, length(methods(c)) > 0 ? (c, Γ, η) : (t -> c, Γ, η)) end
end

# ╔═╡ bc8e8a5f-6476-402f-9cc6-0af722a6c033
ρ0.data

# ╔═╡ e7353f43-899d-4899-86b5-1303676c07ec
ρ0.data[2,2]

# ╔═╡ 230e1d0b-64a5-4ef2-847c-89d344300870
md"""
### Evolution function
"""

# ╔═╡ 2c1f165f-cabd-4599-afc7-348799754f5e
function expevolve(; rs=[])
	# readout distribution
	sim = (length(rs) == 0)
	ts = range(first(T), last(T), step=dt)
	
	# first time step
	(xn, yn, zn) = (x0, y0, z0)	
	(xs, ys, zs, ρs) = ([x0], [y0], [z0], [tr(ρ0*ρ0)])
	ρ = ρ0

	
	if sim
		dist = Normal(0, sqrt(τm/dt))
		push!(rs, 0) end
		

	for i in 2:length(ts)
		
		# sample from distribution		
		if sim
			r = zn + rand(dist)
			push!(rs, r)
		else
			r = rs[i] end
		
 		# update equations
		
		# measurement backaction
		M = exp(r' * dt * DenseOperator(σz) / 2τm)
		ρ = M * ρ * M'
		ρ = ρ / tr(ρ) 
		
		# Hamiltonian evolution
		H = (Δ0 + Δ1*r) * σϕ/2
		U = exp(-im * dt * DenseOperator(H))
		ρ = U * ρ * U'
		
		# Lindblad evolution
		if ideal
			ρ.data[1,2] = ρ.data[1,2] * exp(-dt * (1 - η)/(2τm * η))
			ρ.data[2,1] = ρ.data[2,1] * exp(-dt * (1 - η)/(2τm * η))
			
		else
			ρ.data[1,2] = ρ.data[1,2] * exp(-dt/2T1 - dt/T2 - dt*(1 - η)/(2τm * η))
			ρ.data[2,1] = ρ.data[2,1] * exp(-dt/2T1 - dt/T2 - dt*(1 - η)/(2τm * η))
		end
		
		# store values
		push!(xs, real(expect(ρ, σx)))
		push!(ys, real(expect(ρ, σy)))
		push!(zs, real(expect(ρ, σz)))
		push!(ρs, tr(ρ*ρ))
		

	end
	
	return ts, rs, (xs, ys, zs, ρs)
	
	
end

# ╔═╡ bc62f9e3-3b06-4de2-8867-5b5530c441ae
md"""
### Simulation
"""

# ╔═╡ 90511e52-c343-490b-81b2-1502c546e16b
(tte, re, (xxe, yye, zze, rre)) = expevolve(rs=r)

# ╔═╡ 9e731ac9-e11c-427b-9569-4a925799042a
md" ## QuantumCircuits.jl "

# ╔═╡ ddc50208-bbde-4131-ac93-c254efbbae7a
md"""
### Simulation
"""

# ╔═╡ c012d425-be84-46c9-bb90-b8a5f9e78ca1
sol = bayesian((0, tf), ρ0, H0, J0, C0; dt=1e-3, td=0.0)

# ╔═╡ 3258df22-de38-4ab5-92b7-0a828bc32155
md" # Utilities "

# ╔═╡ 7e1e2e12-48d5-47bd-ac0f-85cebac2bf84
md"""
## Plotting
"""

# ╔═╡ 4080d8f2-6e23-45f3-b34f-d04a050ed7d0
function qubit_plot(sol::Solution; record=false, title="")

	basis = qbasis

	t = sol.t
	exps = map(op -> expectations(sol, op), basis)
	r = record ? sol.r[1] : []

	return qubit_plot((t, exps, r); title=title)
end

# ╔═╡ eb6670fc-cecb-47b7-82ce-fe7a1075f11a
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

# ╔═╡ 53ada28f-90c9-42cd-8bb4-1decd9d1d9e1
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

# ╔═╡ 5da5ba47-cef4-491a-b615-ef80f1557080
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

# ╔═╡ 14c4a025-541d-48aa-a682-7ec6e9f161c4
md"""
#### Colors
"""

# ╔═╡ 3eedcb17-4e02-4e46-9799-b438c469c74c
begin
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

# ╔═╡ 18531b2e-a721-4c12-8e83-ae4844694d8d
function qubit_plot((t, exps, r); title="")

	record = (r != [])

	basis = qbasis
	colors = colorsq1
	labels = qlabels

	p = 0.5 .* (1 .+ exps[1].^2 .+ exps[2].^2 .+ exps[3].^2)

	pl = plot(size=(600,300), legendfontsize=10, titlefontsize=12, legend=:outerright, ylabel="bloch coordinates", xlabel = record ? "" : "t (μs)", linewidth=1.2, title=title)

	for l in 1:length(basis)
		label = labels[l]
		color = colors[l]
		exp = exps[l]
		plot!(t, exps[l], color=color, label=label, legend=:outerright, ylims=[-1,1])
	end

	plot!(t, p, color=colors[4], label=L"Tr(\rho^2)")

	if !record
		return pl
	else
		l = @layout [blochs{0.6h}; record{0.4h}]
		pr = plot(t, r, color=colors1[1], xlabel="t (μs)", ylabel="record", label=:none, legend=:outerright, title="")
		return plot(pl, pr, layout = l, link=:y, linewidth=1.2)
	end
end

# ╔═╡ e1a20bd5-ec19-4c34-a084-bc2fcc22042a
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

# ╔═╡ 3e95667a-9c6a-48f0-b289-ed67b4fdc484
function qubit_plot((t1, exps1, r1), (t2, exps2, r2); title="", color1=colorsq1, color2=colorsq2, l1="", l2="")

	basis = qbasis
	labels = qlabels

	record = ((r1 != []) || (r2 != []))

	p1 = 0.5 .* (1 .+ exps1[1].^2 .+ exps1[2].^2 .+ exps1[3].^2)
	p2 = 0.5 .* (1 .+ exps2[1].^2 .+ exps2[2].^2 .+ exps2[3].^2)

	pl = plot(size=(600,300), legendfontsize=10, titlefontsize=12, legend=:outerright, ylabel="bloch coordinates", xlabel = record ? "" : "t (μs)", linewidth=1.2, title=title)

	for l in 1:length(basis)
		plot!(t1, exps1[l], color=color1[l], label=labels[l], legend=:outerright, ylims=[-1,1], linewidth=1.2)
	end
	plot!(t1, p1, color=color1[4], label=L"Tr(\rho^2)", linewidth=1.2)
	
	for l in 1:length(basis)
		plot!(t2, exps2[l], color=color2[l], label=:none, legend=:outerright, ylims=[-1,1], linestyle=:dash, linewidth=2)
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

# ╔═╡ d898bb17-b80f-4665-8fcc-e67863cb0430
qubit_plot((tt, [xtar, ytar, ztar], []), (tt, [xx, yy, zz], r), l1="target", l2="trajectory")

# ╔═╡ 2a842a56-4d47-4a4e-9c3c-7f7ce5c1ac75
qubit_plot((tt, [xtar, ytar, ztar], []), (tte, [xxe, yye, zze], re), l1="target", l2="trajectory")

# ╔═╡ 477756ac-f4ff-471e-b73a-13183585592a
qubit_plot((tt, [xx, yy, zz], r), (tte, [xxe, yye, zze], re), l1="Bloch", l2="Matrix")

# ╔═╡ 212c9d72-a828-47d3-9343-9e978e304823
let
	t2 = sol.t
	exps2 = map(op -> expectations(sol, op), qbasis)
	r2 = sol.r[1]
	qubit_plot((tt, [xx, yy, zz], r), (t2, exps2, r2), l1="Bloch", l2="QuantumCircuits")
end

# ╔═╡ 5bb19df4-b183-4a3f-b09e-789685ec5d11
md"""
## Misc
"""

# ╔═╡ 36b2bbc2-55ac-45ad-b458-a792c0acddc1
green(text; title="Note") = Markdown.MD(Markdown.Admonition("correct", title, [text]))

# ╔═╡ b5c85a10-9f32-43d4-a4ea-5ec2b5dd37da
red(text; title="Note") = Markdown.MD(Markdown.Admonition("danger", title, [text]))

# ╔═╡ 7bfcd716-0b89-4ff0-85e5-0b7fd65a305d
tan(text; title="Note") = Markdown.MD(Markdown.Admonition("warning", title, [text]))

# ╔═╡ 7fb69f36-7575-4eec-b6de-7008c347eb48
tan(md"Note that here the Bloch and Quantum Circuits stuff doesn't match because inputing a measurement record is not yet supported in QuantumCircuits.jl.")

# ╔═╡ c44d4025-7ce8-4623-bc97-7bbb555914d1
blue(text; title="Note") = Markdown.MD(Markdown.Admonition("note", title, [text]))

# ╔═╡ 3cf94694-4002-4c1c-9593-e61ea5c99b7b
hint(text; title="Hint") = Markdown.MD(Markdown.Admonition("hint", title, [text]))

# ╔═╡ Cell order:
# ╟─b9bd68f2-df48-4072-8269-1898b7cf1b15
# ╟─99a066a7-7013-472e-82d8-3f4c651c3475
# ╟─0fb01de9-882a-46c1-8447-0018146842df
# ╟─18c29abc-1f35-4b5a-bb27-a491c02cc98f
# ╟─01f57775-b647-4fea-8e96-0b8c8ceeff05
# ╠═eae605ed-f411-4f33-8066-bd8f01fc8a2d
# ╠═ee541c05-c187-4b43-a803-2255e254efe5
# ╠═c5de979e-de3b-4a20-9fc4-649851a311fa
# ╟─8aa08bfb-ac91-4e5f-9fb2-dbce02a38b8a
# ╟─c7f72983-771c-46e5-a1ec-94ef80391dd3
# ╟─ea558387-68dc-4dff-ae04-f7a877379510
# ╠═6294f81c-ea9d-4500-9751-b45f8a348639
# ╟─a5f1944a-80b7-4db3-9355-44140e05d908
# ╠═f2892559-e623-46d4-9ff2-f04fd9253734
# ╟─a1358f8d-0622-49ff-bd71-4d4fb5b434f3
# ╠═4487690f-349c-461c-a5b1-f9a6be9897bd
# ╠═d898bb17-b80f-4665-8fcc-e67863cb0430
# ╟─15613717-87c1-45d8-b194-1c7b1085d4f7
# ╠═840beeea-a731-4a38-b36b-e2ab848147d9
# ╠═bc8e8a5f-6476-402f-9cc6-0af722a6c033
# ╠═e7353f43-899d-4899-86b5-1303676c07ec
# ╟─230e1d0b-64a5-4ef2-847c-89d344300870
# ╠═2c1f165f-cabd-4599-afc7-348799754f5e
# ╟─bc62f9e3-3b06-4de2-8867-5b5530c441ae
# ╠═90511e52-c343-490b-81b2-1502c546e16b
# ╠═2a842a56-4d47-4a4e-9c3c-7f7ce5c1ac75
# ╠═477756ac-f4ff-471e-b73a-13183585592a
# ╟─9e731ac9-e11c-427b-9569-4a925799042a
# ╟─ddc50208-bbde-4131-ac93-c254efbbae7a
# ╠═c012d425-be84-46c9-bb90-b8a5f9e78ca1
# ╟─7fb69f36-7575-4eec-b6de-7008c347eb48
# ╠═212c9d72-a828-47d3-9343-9e978e304823
# ╟─3258df22-de38-4ab5-92b7-0a828bc32155
# ╟─7e1e2e12-48d5-47bd-ac0f-85cebac2bf84
# ╟─4080d8f2-6e23-45f3-b34f-d04a050ed7d0
# ╟─18531b2e-a721-4c12-8e83-ae4844694d8d
# ╟─e1a20bd5-ec19-4c34-a084-bc2fcc22042a
# ╟─3e95667a-9c6a-48f0-b289-ed67b4fdc484
# ╟─eb6670fc-cecb-47b7-82ce-fe7a1075f11a
# ╟─53ada28f-90c9-42cd-8bb4-1decd9d1d9e1
# ╟─5da5ba47-cef4-491a-b615-ef80f1557080
# ╟─14c4a025-541d-48aa-a682-7ec6e9f161c4
# ╟─3eedcb17-4e02-4e46-9799-b438c469c74c
# ╟─5bb19df4-b183-4a3f-b09e-789685ec5d11
# ╟─36b2bbc2-55ac-45ad-b458-a792c0acddc1
# ╟─b5c85a10-9f32-43d4-a4ea-5ec2b5dd37da
# ╟─7bfcd716-0b89-4ff0-85e5-0b7fd65a305d
# ╟─c44d4025-7ce8-4623-bc97-7bbb555914d1
# ╟─3cf94694-4002-4c1c-9593-e61ea5c99b7b
# ╠═1994807c-85a7-4309-8786-c1b5d2dd6b86
