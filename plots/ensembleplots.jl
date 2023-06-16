import Statistics: mean

"-- Color palettes ----------------------------------------------------------------------------------------------------"
# single qubit
colors1q = palette(:tab10)
qcolors = colors1q

"Ensembles ---------------------------------------------------------------------------------------------------------"
mutable struct BlochEnsemble end
const blochensemble = BlochEnsemble()

@recipe f(be::BlochEnsemble, ens::Ensemble) = (be, ens.t, ens.exps...)
@recipe f(be::BlochEnsemble, ens::Ensemble, sol::Solution) = (be, ens.t, ens.exps..., sol.exps)
@recipe function f(::BlochEnsemble, t, xs, ys, zs, η0exps = []; N=100)

	n = min(N, length(xs))

	legend --> :none
    titlefontsize --> 11
    xtickfontsize --> 10
    ytickfontsize --> 10
    xguidefontsize --> 10
	yguidefontsize --> 10
    size --> (600,400)
    linewidth --> 1
	linestyle --> :solid

	linealpha --> 0.1
	layout := @layout [ a; b; c ]
	ylims --> [-1,1]

	xlabel --> "t (μs)"

	for (i, cs) in enumerate([xs, ys, zs])

		subplot := i
		color := qcolors[i]
		ylabel := qlabels[i]

		for c in cs[1:n]
			@series begin
				label := ""
				if i < 3
					xlabel := ""
				end
				t, c
			end
		end

		@series begin
			label := ""
			if i < 3
				xlabel := ""
			end
			label := "ensemble average"
			linealpha := 1.0
			linewidth := 4
			t, mean(cs)
		end

		if length(η0exps) >= i
			@series begin
				if i < 3
					xlabel := ""
				end
				label := "η = 0"
				linealpha := 1.0
				linewidth := 1
				color := :black
				t, η0exps[i]
			end
		end
	end
end # @recipe BlochEnsemble




mutable struct SeriesHistogram end
const serieshistogram = SeriesHistogram()

@recipe function f(::SeriesHistogram, ts, cs, index = length(ts); plotline = true)

	t = ts[index]

	color --> colors1q[3]
	legend := :none
	linealpha --> 0.5
	linewidth --> 0.7
	size --> (600,300)
	xlabel := "t (μs)"
	ylabel --> "bloch coordinate"

	layout := @layout [a{0.7w} b{0.3w}]

	for (i, c) in enumerate(cs)
		@series begin
			subplot := 1
			ts, c
		end
	end

	if plotline
		@series begin
			subplot := 1
			linestyle := :dash
			color := :black
			linewidth := 1.5
			linealpha := 1
			[t, t], [-1, 1]
		end
	end

	seriestype := :histogram
	@series begin
		subplot := 2
		bins := (-1):0.1:1
		normalize := false
		alpha := 0.9
		orientation := :h
		xlabel := "occurrences"
		ylabel := ""
		xlims := [0, round(length(cs)/2)]

		map(c -> c[index], cs)
	end
end
