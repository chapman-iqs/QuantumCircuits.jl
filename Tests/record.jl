"record is distributed as expected for measurement eigenstates"
function record_distribution(ψ0, τ, dt; test_id="record", testset_id="", samples=10000, solve=bayesian, tolerance = 0.25, makeplots=false, plotpath="test_result_plots", kwargs...)

    tf = dt * samples
    η = 1.0
    Γ = 1/(2τ*η)

    C = [(σz, Γ, η)]

    sol = solve((0.0, tf), ψ0, Iq, [], C; dt=dt)
    record = sol.r[1]

    μ, σ = params(fit(Normal, record))
    σ2 = σ^2

    μ_expected = real(expect(σz, ψ0))
    @assert isa(μ_expected, Float64)
    σ2_expected = τ / dt

    if makeplots
        path = mkpath(joinpath(plotpath, string(solve), testset_id, test_id, ψ0 == g ? "ground" : "excited"))
        xs = range(extrema(record)..., length=100)
        histogram(record; label="data", normalize=:pdf)
        plot!(xs, gaussian(μ, σ).(xs); label="fit")
        plot!(xs, gaussian(μ_expected, sqrt(σ2_expected)).(xs); label="expected", lineestyle=:dash, color=:black)
        plot!(title="μ = $(μ), σ^2 = $(σ2)\n μ_exp = $(μ_expected), σ2_exp = $(σ2_expected)")
        savefig(joinpath(path, "τ_$(τ)_dt_$(dt).png"))
    end

    return prod([abs(dif(μ, μ_expected) / μ_expected), abs(dif(σ2, σ2_expected) / σ2_expected)] .< tolerance)
end
gaussian(μ, σ) = x -> (1/(σ * sqrt(2π))) * exp(-(x - μ)^2 / (2σ^2) )

function test_readout(; solve=bayesian, testset_id="readout", kwargs...)

    @testset "record distribution" begin
        for ψ0 in [g, e]
            for τ in 0.002:0.002:0.01
                for dt in [1e-4, 1e-3, 1e-2]
                    @test record_distribution(ψ0, τ, dt; solve=solve, testset_id=testset_id, kwargs...)
                end
            end
        end
    end
end

# function multilindblad(; γ1 = 0.1, γ2 = 0.2, γ3 = 0.3, kwargs...)

#     H = Iq
#     Γ1 = γ1
#     Γ2(t::Timescale) = γ2
#     Γ3(t::Timescale, ρ::State) = γ3

#     J = [(σz, Γ1), (σx, Γ2), (σy, Γ3)]
#     C = []
#     tf = 10.0
#     ψ0 = normalize(g + 0.3 * im * e)

#     sol = bayesian((0.0, tf), ψ0, H, J, C; dt = 1e-3)
# end
