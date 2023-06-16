"I can't remember what this actually is doing... so we may want to replace with better tests. (sacha)"
function timedelay(td; solve=bayesian, makeplots=false, test_id="timedelay", testset_id="", plotpath="test_result_plots", kwargs...)

    Ωr = 2π * 0.15
    H1(t::Timescale, ρ::State) = Ωr * σx * (1 - fidelity(ρ, normalize(g + e)))
    ψ0 = g

    (t0, tf) = 0.0, 10.0
    dt = 1e-3

    sol = solve((t0, tf), ψ0, H1, [], []; dt = dt, td = td)
    make_test_plot(makeplots, sol, solve, plotpath, testset_id, test_id)
    # if makeplots
    #     plot(blochtimeseries, sol.t, sol.exps...)
    #     savefig(joinpath(plotpath, string(solve), "timedelay_td_$td.png"))
    # end

    # solve by hand and test
    ψs = [ψ0]
    ψ = ψ0
    ts = t0:dt:tf
    td_index = round(Int, (td - t0) / dt) + 1
    for t in ts[2:end]
        ψd = t <= td ? ψ0 : reverse(ψs)[td_index]
        u = exp( -im * dt * DenseOperator(H1(t, ψd)))
        ψ = u * ψ
        push!(ψs, ψ)
    end

    exps = map(op -> expectations(ψs, op), qbasis)
    # return plot(blochtimeseries, ts, exps...)
    return maximum(map(i -> abs(1 - fidelity(ψs[i], sol.ρ[i])), 1:length(sol.ρ))) < 0.00001
end

function filter_forward_estimate(td; solve=bayesian, ψ0 = randstate(SpinBasis(1//2)), tolerance=1e-9, plotpath="test_result_plots", test_id="filter_forward_estimate", testset_id="", makeplots=true, kwargs...)
    Γ = 0.5
    η = 1.0
    tf = 5.0
    dt = 1e-2

    H(t::Timescale, ρ::State, r::Record) = Iq
    C = [(σz, Γ, η)]

    sys, fil, fow = solve(hamiltonianfe, (0.0, tf), ψ0, (H,H), ([],[]), C; td=td, dt=dt)

    delay = Int64(round(td / dt) + 1)
    fow_shifted = dm.(fow.ρ[(delay + 1):end])
    fil_cutoff = dm.(fil.ρ[1:length(fil.ρ) - (delay + 1)])

    average_trace_distance = mean(map(zip(fow_shifted, fil_cutoff)) do (ρfow, ρfil) tracedistance(ρfow, ρfil) end)

    pass = average_trace_distance < tolerance
    pass_string = pass ? "passed" : "failed"

    if makeplots   
        path = mkpath(joinpath(plotpath, string(solve), testset_id, test_id))
        sys, fil, fow = map(s -> Solution(s, qbasis), [sys,fil,fow])
        plot(blochtimeseries, sys; linewidth=7, linealpha=0.35, title="$(pass_string), td = $td, d = $(average_trace_distance)")
        plot!(blochtimeseries, fil; label="filter")
        plot!(blochtimeseries, fow; linestyle=:dash, label="forward est.")
        savefig(joinpath(path, "$(pass_string)_td_$(td).png"))
    end

    return pass
end

function test_timedelay(; solve=bayesian, testset_id="timedelay", kwargs...)

    # @testset "time delay" begin
    #     if solve==bayesian
    #         for td in 0:0.2:2.0
    #             @test timedelay(td; solve=solve, testset_id=testset_id, kwargs...)
    #         end
    #     else
    #         println("⬇︎⬇︎⬇︎ Time delay not yet implemented for rouchon method.")
    #     end
    #  end

     @testset "filter and forward estimate" begin
        if solve==bayesian
            for td in 0.0:0.5:2.0
                @test filter_forward_estimate(td; solve=solve, testset_id=testset_id, kwargs...)
            end
        else
            println("⬇︎⬇︎⬇︎ Time delay not yet implemented for rouchon method.")
        end
    end
end

