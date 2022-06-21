"-- Color palettes ----------------------------------------------------------------------------------------------------"
# single qubit
colors1q = palette(:tab10)

# system-filter
colors_system = palette(:rainbow)
colors_filter = palette(:lightrainbow)

"Ensembles ---------------------------------------------------------------------------------------------------------"
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

"Single trajectories ---------------------------------------------------------------------------------------------------------"
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
