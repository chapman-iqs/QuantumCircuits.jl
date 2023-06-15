"record is distributed as expected for measurement eigenstates"
function bayesian_dt_robustness(dt, scale;  ψ0 = QuantumOpticsBase.randstate(SpinBasis(1//2)), test_id="bayesian_dt_robustness", testset_id="", tolerance = 1e-2, makeplots=false, plotpath="test_result_plots", kwargs...)

    η = 1.0
    Γ = 0.2
    tf = 10.0

    C = [(σz, Γ, η)]

    sol1 = bayesian((0.0, tf), ψ0, Iq, [], C; dt=dt)
    solcoarse = coarse_sol(sol1, scale)

    sol2 = bayesian((0.0, tf), ψ0, Iq, [], C; dt=dt*scale, records=solcoarse.r)

    average_trace_distance = mean(map(zip(dm.(solcoarse.ρ), dm.(sol2.ρ))) do (ρ1, ρ2) tracedistance(ρ1, ρ2) end)
    pass = average_trace_distance < tolerance
    pass_string = pass ? "passed" : "failed"

    if makeplots
        path = mkpath(joinpath(plotpath, string(bayesian), testset_id, test_id))
        sol1 = Solution(sol1, qbasis)
        sol2 = Solution(sol2, qbasis)
        solcoarse = Solution(solcoarse, qbasis)

        p = if length(sol2.ρ) <= 100
            plot(blochtimeseries, sol1; title="dt = $dt, scale = $scale \n $pass_string d = $average_trace_distance", legend=:outerright, size=(800,600))
            plot!(blochtimeseries, sol2; seriestype=:scatter, label="reconstructed", markersize=4, markershape=:square, legend=:outerright)
            plot!(blochtimeseries, solcoarse; seriestype=:scatter, label="subselected", markersize=3, markershape=:circle, markercolor=:black, legend=:outerright, topmargin=10mm)
        else
            plot(blochtimeseries, sol1; title="dt = $dt, scale = $scale \n $pass_string d = $average_trace_distance", legend=:outerright, size=(800,600))
            plot!(blochtimeseries, sol2; linewidth=7, linealpha=0.35, label="reconstructed", legend=:outerright)
            plot!(blochtimeseries, solcoarse; linewidth=5, linestyle=:dash, linealpha=0.35, label="subselected", legend=:outerright, topmargin=10mm)
        end
        savefig(joinpath(path, "$(pass_string)_dt_$(dt)_scale_$(scale).png"))
    end
    
    return pass
end

function coarse_sol(sol::Solution, scale)
    t0, tf = first(sol.t), last(sol.t)
    dt = sol.t[2] - sol.t[1]
    return coarse_sol((t0, tf), dt, scale, sol)
end
coarse_sol((t0, tf), dt, scale, sol::Solution) = coarse_sol((t0, tf), dt, scale, sol.r, sol.ρ)
function coarse_sol((t0, tf), dt, scale, records, ρs)

    newdt = scale * dt
    newtimes = collect(t0:newdt:tf)
    newrecords::Vector{Record} = map(r -> subselect(coarse_grain(r; n=scale); n=scale), records)
    new_ρs = subselect(ρs; n=scale)

    # account for discrepancies in length from subselection
    for (old_r, new_r) in zip(records,  newrecords)
        if length(new_r) == length(newtimes) - 1
            pushfirst!(new_r, old_r[1])
        elseif length(new_r) == length(newtimes) + 1
            pop!(new_r)
        end
        @assert length(new_r) == length(newtimes)
    end

    if length(new_ρs) == length(newtimes) - 1
        pushfirst!(new_ρs, ρs[1])
    elseif length(new_ρs) == length(newtimes) + 1
        pop!(new_ρs)
    end
    @assert length(new_ρs) == length(newtimes)

    return Solution(newtimes, new_ρs, newrecords)
end

function test_bayesian(; solve=bayesian, testset_id="bayesian_specific", kwargs...)

    @testset "bayesian dt robustness" begin
        if solve==bayesian
            for scale in [5,10,25,50,100,150]
                for dt in [1e-4, 1e-3, 1e-2]
                    @test bayesian_dt_robustness(dt, scale; solve=solve, testset_id=testset_id, kwargs...)
                end
            end
        else
            println("⬇︎⬇︎⬇ No dt robustness implemented for rouchon.")
        end
    end
end
