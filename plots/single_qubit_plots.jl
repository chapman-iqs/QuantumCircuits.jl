using Plots, Measures, LaTeXStrings

"-- Color palettes ----------------------------------------------------------------------------------------------------"
# single qubit
colors1q = palette(:tab10)

"Single trajectories ---------------------------------------------------------------------------------------------------------"


mutable struct BlochTimeSeries end
const blochtimeseries = BlochTimeSeries()

@recipe f(bts::BlochTimeSeries, sol::Solution) = (bts, sol.t, sol.exps...)
@recipe f(bts::BlochTimeSeries, sol::Solution, ::Records) = (bts, sol.t, sol.exps..., sol.r)
@recipe function f(::BlochTimeSeries, ts, xs, ys, zs; vec=nothing, tf=nothing)
	ps = 0.5 .* (1 .+ xs.^2 + ys.^2 + zs.^2) # purity

	# (xv, yv, zv) = vec

	if vec != nothing
		(xv, yv, zv) = vec
		ρv = 0.5 * (1 + xv^2 + yv^2 + zv^2)
	end

	# Plot time series --------------------------------------------------------

	legend --> :topright
	label --> [L"$x$" L"$y$" L"$z$" L"$Tr ( \rho^2 )$"]
	xlabel --> "t (μs)"
	ylabel --> "bloch coordinates"

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


	for (i, bs) in enumerate([xs, ys, zs, ps])

		@series begin
			color := palette(:tab10)[i]
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
@recipe function f(::BlochTimeSeries, ts, xs, ys, zs, records::Records; vec=nothing, tf=nothing, ylabels = ["bloch coordinates", "record 1", "record 2"])
	ps = 0.5 .* (1 .+ xs.^2 + ys.^2 + zs.^2) # purity

	if vec != nothing
		(xv, yv, zv) = vec
		ρv = 0.5 * (1 + xv^2 + yv^2 + zv^2)
	end

	# Plot time series --------------------------------------------------------

	layout := @layout [b{0.5h}; grid(length(records),1)]
	link := :x

	legend --> :topright
	label --> [L"$x$" L"$y$" L"$z$" L"$Tr ( \rho^2 )$"]
	xlabel --> "t (μs)"
	# ylabel --> "bloch coordinates"

	palette := :tab10
	linealpha --> 1

	legendfontsize --> 10
	titlefontsize --> 12
	xtickfontsize --> 9
	ytickfontsize --> 9
	xguidefontsize --> 10
	yguidefontsize --> 10
	size --> (600,500)
	linewidth --> 1.5

	tf = (tf == nothing) ? last(ts) : tf
	xlims --> [first(ts), tf]


	for bs in (xs, ys, zs, ps)

		@series begin
			subplot := 1
			ylims --> [-1, 1]
			xlabel := ""
			ylabel := ylabels[1]
			ts, bs
		end

	end

	if vec != nothing

		for (i, sv) in enumerate((xv, yv, zv, ρv))

			linecolor := i

			@series begin
				subplot := 1
				linestyle := :dot
				label := :none
				xlabel := ""
				[first(ts), tf], [sv, sv]
			end
		end
	end

	for (i, rs) in enumerate(records)
		@series begin
			subplot := 1 + i
			color := palette(:rainbow)[i]
			linewidth := 0.6
			label := :none
			ylims := [minimum(rs), maximum(rs)]
			if i < length(records)
				xlabel := ""
			end
			ylabel := ylabels[i + 1]
			ts, rs
		end
	end
end

@userplot BlochProjections
@recipe function f(bp::BlochProjections; plottitle="Bloch trajectory cross-sections", vec=nothing, blochmark=false, blochmarkcolor="white")
	bloch = bp.args # has the form bloch = (xs, ys, zs)

	# Define colors and labels -------------------------------------------------

	colors = palette(:tab10)[1:3]
	labels = qlabels
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
