@userplot BlochSphere
@userplot BlochTimeSeries
@userplot BlochProjections



@recipe function f(bs::BlochSphere; mesh=30, ax=false, viewϕ=0, vec=nothing, blochmark=false, blochmarkcolor="white")
	xss, yss, zss = bs.args

	xs = xss .* cos(viewϕ) .- yss .* sin(viewϕ)
	ys = xss .* sin(viewϕ) .+ yss .* cos(viewϕ)
	zs = zss


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


end

@recipe function f(bts::BlochTimeSeries; vec=nothing, tf=nothing)
	ts, xs, ys, zs = bts.args
	ps = 0.5 .* (1 .+ xs.^2 + ys.^2 + zs.^2) # purity

	# (xv, yv, zv) = vec

	if vec != nothing
		(xv, yv, zv) = vec
		ρv = 0.5 * (1 + xv^2 + yv^2 + zv^2)
	end

	# Plot time series --------------------------------------------------------

	legend --> :bottomright
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
	size --> (400,300)
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
