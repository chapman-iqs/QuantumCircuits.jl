"-- Color palettes ----------------------------------------------------------------------------------------------------"
# single qubit
colors1q = palette(:tab10)

# system-filter
colors_system = palette(:rainbow)
colors_filter = palette(:lightrainbow)

# multiple qubits
colors2q_bell = palette(:rainbow)
colors3q_prod = palette(:okabe_ito)
colors3q_number = palette(:lightrainbow)

"-- Recipes ----------------------------------------------------------------------------------------------------"

@userplot BlochSphere
@userplot BlochTimeSeries
@userplot BlochProjections

@recipe function f(bs::BlochSphere; mesh=30, ax=false, viewϕ=0, vec=nothing, blochmark=false, blochmarkcolor="white")

	args = bs.args
	timed = length(args) == 4

	xss, yss, zss = timed ? args[2:4] : args
	i = timed ? bs.args[1] : length(xss)

	xs = xss[1:i] .* cos(viewϕ) .- yss[1:i] .* sin(viewϕ)
	ys = xss[1:i] .* sin(viewϕ) .+ yss[1:i] .* cos(viewϕ)
	zs = zss[1:i]


	# Plot trajectory ----------------------------------------------------------

	# marker := nothing
	linecolor --> "black" # makes blue only if color is not specified
	linewidth --> :2

	@series begin
		xs, ys, zs
	end


	# Wire frame coordinates ---------------------------------------------------

	x(θ, ϕ) = sin(θ) * cos(ϕ + viewϕ)
	y(θ, ϕ) = sin(θ) * sin(ϕ + viewϕ)
	z(θ, ϕ) = cos(θ)

	θs = range(0, 2π, length=mesh)
	ϕs = range(0, π, length=mesh) .+ viewϕ


	# Plot wireframe -----------------------------------------------------------


	legend := false
	linecolor := "steelblue"
	linewidth := 0.5
	linealpha := 1
	seriestype := path3d
	aspect_ratio := 1.0
	size --> (400,400)

	# Longitudes
	for ϕ in ϕs
		@series begin
			[x(θ, ϕ) for θ in θs], [y(θ, ϕ) for θ in θs], [z(θ, ϕ) for θ in θs]
		end
	end

	# Latitudes
	for θ in θs
		@series begin
			[x(θ, ϕ) for ϕ in ϕs], [y(θ, ϕ) for ϕ in ϕs], [z(θ, ϕ) for ϕ in ϕs]
		end
	end


	# Plot reference axes ------------------------------------------------------
	linewidth := 3
	colors = [palette(:tab10)[i] for i in 1:3]

	if ax

		linecolor := colors[1]
		@series begin
			[0, cos(viewϕ)], [0, sin(viewϕ)], [0, 0]
		end

		linecolor := colors[2]
		@series begin
			[0, -sin(viewϕ)], [0, cos(viewϕ)], [0, 0]
		end

		linecolor := colors[3]
		@series begin
			[0, 0], [0, 0], [0, 1]
		end
	end


	# Plot optional Bloch vector input by user --------------------------------

	if vec != nothing

		(xvv, yvv, zvv) = vec

		xv = xvv * cos(viewϕ) - yvv * sin(viewϕ)
		yv = xvv .* sin(viewϕ) .+ yvv * cos(viewϕ)
		zv = zvv

		linewidth := 2
		linecolor := "red"

		@series begin
			[0, xv], [0, yv], [0, zv]
		end

		marker := (:circle, 5)
		markercolor := "red"
		@series begin
			[xv], [yv], [zv]
		end

	end

	# Plot optional marker ------------------------------------------------

	if blochmark

		marker := (:circle, 3.5)
		markercolor := blochmarkcolor
		@series begin
			[last(xs)], [last(ys)], [last(zs)]
		end

	end
end # blochsphere

@recipe function f(bts::BlochTimeSeries; vec=nothing, tf=nothing)
	ts, xs, ys, zs = bts.args
	ps = 0.5 .* (1 .+ xs.^2 + ys.^2 + zs.^2) # purity

	# (xv, yv, zv) = vec

	if vec != nothing
		(xv, yv, zv) = vec
		ρv = 0.5 * (1 + xv^2 + yv^2 + zv^2)
	end

	# Plot time series --------------------------------------------------------

	legend --> :outerright
	label --> [L"$x$" L"$y$" L"$z$" L"$Tr ( \rho^2 )$"]
	xlabel --> "t (μs)"
	ylabel --> "bloch coordinates"

	palette := :tab10
	linealpha --> 1

	legendfontsize --> 10
	titlefontsize --> 12
	xtickfontsize --> 10
	ytickfontsize --> 10
	xguidefontsize --> 10
	yguidefontsize --> 10
	size --> (600,300)
	linewidth --> 1.5

	tf = (tf == nothing) ? last(ts) : tf
	xlims --> [first(ts), tf]
	ylims --> [-1, 1]


	for bs in (xs, ys, zs, ps)

		@series begin
			ts, bs
		end

	end

	linestyle := :dot
	label := :none
	if vec != nothing

		for (i, sv) in enumerate((xv, yv, zv, ρv))

			linecolor := i

			@series begin
				[first(ts), tf], [sv, sv]
			end

		end

	end
end

@recipe function f(bp::BlochProjections; plottitle="Bloch trajectory cross-sections", vec=nothing, blochmark=false, blochmarkcolor="white")
	bloch = bp.args # has the form bloch = (xs, ys, zs)

	# Define colors and labels -------------------------------------------------

	colors = [palette(:tab10)[i] for i in 1:3]
	labels = ["x", "y", "z"]
	palette := :tab10
	linealpha --> 1
	linewidth --> 0.85

	linecolor --> "black"

	function setaxes(i1, i2)
		# get colors
		c1 = colors[i1]
		c2 = colors[i2]

		# get labels
		l1 = labels[i1]
		l2 = labels[i2]

		# set colors
		x_guidefontcolor := c1
		y_guidefontcolor := c2
		x_foreground_color_axis := c1
		y_foreground_color_axis := c2
		x_foreground_color_text := c1
		y_foreground_color_text := c2
		x_foreground_color_border := c1
		y_foreground_color_border := c2

		# set labels
		xlabel := labels[i1]
		ylabel := labels[i2]
	end


	# General plot settings ---------------------------------------------------

	legend := false
	aspectratio := 1
	layout := @layout [ a b c ]
	titlefontsize --> 12

	xlims := [-1.1, 1.1]
	ylims := [-1.1, 1.1]
	xticks := [-1.0, 0, 1.0]
	yticks := [-1.0, 0, 1.0]

	titles = ["", plottitle, ""]


	# Plot -------------------------------------------------------------------

	for (i, (i1, i2)) in enumerate([(1, 3), (1, 2), (2, 3)])

		b1, b2 = (bloch[i1], bloch[i2])

		# Plot trajectory - - - - - - - - - - - - - - - - - - - - - - - - -

		title := titles[i]
		setaxes(i1, i2)

		@series begin

			subplot := i
			b1, b2

		end


		# Plot optional Bloch vector input by user - - - - - - - - - - - - -

		if vec != nothing

			v1, v2 = (vec[i1], vec[i2])

			@series begin

				linewidth := 2
				linestyle := :dot
				linecolor := "red"
				subplot := i

				[0, v1], [0, v2]

			end

			@series begin

				marker := (:circle, 5)
				markercolor := "red"
				subplot := i

				[v1], [v2]

			end

		end


		# Plot optional marker  - - - - - - - - - - - - - - - - - - - -

		if blochmark

			(l1, l2) = (last(b1), last(b2))

			@series begin

				marker := (:circle, 3.5)
				markercolor := blochmarkcolor
				subplot := i
				[l1], [l2]

			end

		end


	end
end


"-- Functions ----------------------------------------------------------------------------------------------------"

# for plotting ensmebles of single-qubit trajectories
function bloch_plots(sols::Vector{Solution}, sol_η0::Solution)
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
			plot!(t, o, alpha=0.1, label=:none, color=color)
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

# for plotting reduced qubit evolution in a two-qubit system
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

# for plotting bell state trajectories
function bell_plot(sol::Solution)

	basis = bell_basis
	colors = colors2q_bell
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
