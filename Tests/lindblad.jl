expdecay(Γ) = t -> exp(-Γ * t)
gaussdecay(γ) = t -> exp(-(γ*t)^2)
expgauss(Γ, γ) = t -> exp(-Γ * t) * exp(-(γ*t)^2)

# checks we obtain the correct exponential solution for σz lindblad decay
function lindblad_const(γ; solve=bayesian, test_id="lindblad_const", testset_id="", tolerance = 0.001, makeplots=false, plotpath="test_result_plots", kwargs...)

    tf = 4.0
    dt = 1e-3
    ψ0 = normalize(g + e)

    H = Iq
    J = [(σz, γ/2)]
    C = []

    sol = solve((0.0, tf), ψ0, H, J, C; dt=dt)

    x, y, z = map(op -> expectations(sol, op), [σx, σy, σz])

    x_expected = expdecay(γ).(sol.t)
    y_expected = 0.0 .* sol.t
    z_expected = 0.0 .* sol.t

    pass = prod([dif(x, x_expected), dif(y, y_expected), dif(z, z_expected)] .< tolerance)
    pass_string = pass ? "passed" : "failed"

    if makeplots
        path = mkpath(joinpath(plotpath, string(solve), testset_id, test_id))
        times = collect(0.0:dt:tf)
        plot(blochtimeseries, times, x, y, z)
        plot!(blochtimeseries, times, x_expected, y_expected, z_expected; linestyle=:dash, color=:black, linewidth=0.7)
        savefig(joinpath(path, "$(pass_string)_γ_$γ.png"))
    end

    return pass
end

# tests we can reproduce exponential-gaussian decay curves using time-dependent Lindblad operators
function lindblad_timedep(Γ, γ; test_id="lindblad_timedep", testset_id="", solve=bayesian, tolerance = 0.001, makeplots=false, plotpath="test_result_plots", kwargs...)

    ψ0 = normalize(g + e)
    tf = 4.0
    dt = 1e-3

	H = Iq
	γL(t::Timescale) = 0.5 * (Γ + 2γ^2 * t)
	J = [(σz, γL)]
	C = []

    sol = solve((0.0, tf), ψ0, H, J, C; dt=dt)

    x, y, z = map(op -> expectations(sol, op), [σx, σy, σz])

    x_expected = expgauss(Γ, γ).(sol.t)
    y_expected = 0.0 .* sol.t
    z_expected = 0.0 .* sol.t

    pass = prod([dif(x, x_expected), dif(y, y_expected), dif(z, z_expected)] .< tolerance)
    pass_string = pass ? "passed" : "failed"

    if makeplots
        path = mkpath(joinpath(plotpath, string(solve), testset_id, test_id))
        times = collect(0.0:dt:tf)
        plot(blochtimeseries, times, x, y, z)
        plot!(blochtimeseries, times, x_expected, y_expected, z_expected; linestyle=:dash, color=:black, linewidth=0.7)
        savefig(joinpath(path, "$(pass_string)_γ_$(γ)_Γ_$(Γ).png"))
    end

    return pass
end

# needs to be written still
function lindblad_statedep(Γ, γ; test_id="lindblad_statedep", testset_id="", solve=bayesian, tolerance = 0.001, kwargs...)

    ψ0 = normalize(g + e)
    tf = 4.0
    dt = 1e-3

    H = Iq
	γL(t::Timescale, ρ::State) = 0.5 * (Γ + 2γ^2 * (1 - expect(σz, ρ)^2))
	J = [(σz, γL)]
	C = []

    sol = solve((0.0, tf), ψ0, H, J, C; dt=dt)
	
    x, y, z = map(op -> expectations(sol, op), [σx, σy, σz])

    println("lindblad_statedep test still needs to be written")
    return true
end


dif(x, y) = maximum(abs.(x .- y))

function test_lindblad(; testset_id="lindblad_master_eq", solve=bayesian, tolerance = 0.001, kwargs...)

        @testset "lindblad constant decay (master equation comparison)" begin
            for γ in 0:0.5:2.0
                @test lindblad_const(γ; testset_id=testset_id, solve=solve, tolerance=tolerance, kwargs...)
                # @test (x && y && z)
            end
        end

        @testset "lindblad time-dep decay (master equation comparison)" begin
            for γ in 0:0.5:1.0
                for Γ in 0:0.5:1.0
                    @test lindblad_timedep(Γ, γ; testset_id=testset_id, solve=solve, tolerance=tolerance, kwargs...)
                end
            end
        end

        # @testset "lindblad state-dep (run)" begin
        #     @test lindblad_statedep(0.2, 0.5; tolerance=tolerance, kwargs...)
        # end

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
