function timedelay(td; showplot=false, kwargs...)

    Ωr = 2π * 0.15
    H1(t::Timescale, ρ::State) = Ωr * σx * (1 - fidelity(ρ, normalize(g + e)))
    ψ0 = g

    (t0, tf) = 0.0, 10.0
    dt = 1e-3

    sol = bayesian((t0, tf), ψ0, H1, [], []; dt = dt, td = td)
    sol = Solution(sol, qbasis)
    if showplot
        plot(blochtimeseries, sol.t, sol.exps...)
    end

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

function test_timedelay(; kwargs...)

    @testset "time delay" begin
        for td in 0:0.2:2.0
            @test timedelay(td)
        end
    end
end
