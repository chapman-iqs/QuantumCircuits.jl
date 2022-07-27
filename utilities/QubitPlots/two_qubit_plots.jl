
# for plotting reduced qubit evolution in a two-qubit system
function single_qubit_plots(sol::Solution)
	colors = colors1q

	# calculate expectation values --------------------------------------------
	q1_basis = [σx1, σy1, σz1]
	q2_basis = [σx2, σy2, σz2]

	exps1 = map(op -> expectations(sol, op), q1_basis)
	exps2 = map(op -> expectations(sol, op), q2_basis)

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


@userplot BellPlot
@recipe function f(bp::BellPlot)
	sol = bp.args[1]

	basis = bell_basis
	exps = map(op -> expectations(sol, dm(op)), basis)

	# Plot time series --------------------------------------------------------

	title --> "Bell states"
	legend --> :outerright
	label --> hcat(bell_basis_labels...)
	xlabel --> "t (μs)"
	ylabel --> "Bell state populations"

	palette := :rainbow
	linealpha --> 1

	legendfontsize --> 10
	titlefontsize --> 12
	xtickfontsize --> 10
	ytickfontsize --> 10
	xguidefontsize --> 10
	yguidefontsize --> 10
	size --> (600,300)
	linewidth --> 1.5
	margin --> 5mm

	ylims --> [0, 1]


	for exp in exps

		@series begin
			sol.t, exp
		end

	end
end
