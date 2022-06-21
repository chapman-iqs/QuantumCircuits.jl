function average_purity(solutions::Vector{Solution})
	return map(ρ -> real(tr(ρ * ρ)), ensemble_average(solutions))
end

function ensemble_average(solutions::Vector{Solution})

	if typeof(solutions[1].ρ[1]) <: Ket
		return [mean(map(sol -> dm(sol.ρ[i]), solutions)) for i in 1:length(solutions[1].t)]
	else
		return [mean(map(sol -> sol.ρ[i], solutions)) for i in 1:length(solutions[1].t)]
	end

end
