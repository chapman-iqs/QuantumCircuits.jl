expdecay(Γ) = t -> exp(-Γ * t)
gaussdecay(γ) = t -> exp(-(γ*t)^2)
expgauss(Γ, γ) = t -> exp(-Γ * t) * exp(-(γ*t)^2)

# checks we obtain the correct exponential solution for σz lindblad decay
function lindblad(γ; tolerance = 0.001)

    tf = 10.0
    ψ0 = normalize(g + e)

    H = Iq
    J = [(σz, γ/2)]
    C = []

    sol = bayesian((0.0, tf), ψ0, H, J, C; dt = 1e-3)

    x, y, z = map(op -> expectations(sol, op), [σx, σy, σz])

    x_expected = expdecay(γ).(sol.t)
    y_expected = 0.0 .* sol.t
    z_expected = 0.0 .* sol.t

    return [dif(x, x_expected), dif(y, y_expected), dif(z, z_expected)] .< tolerance
end

# tests we can reproduce exponential-gaussian decay curves using time-dependent Lindblad operators
function lindblad_timedep(Γ, γ; tolerance = 0.001)

    ψ0 = normalize(g + e)
    tf = 10.0

	H = Iq
	γL(t::Timescale) = 0.5 * (Γ + 2γ^2 * t)
	J = [(σz, γL)]
	C = []

    sol = bayesian((0.0, tf), ψ0, H, J, C; dt = 1e-3)

    x, y, z = map(op -> expectations(sol, op), [σx, σy, σz])

    x_expected = expgauss(Γ, γ).(sol.t)
    y_expected = 0.0 .* sol.t
    z_expected = 0.0 .* sol.t

    return [dif(x, x_expected), dif(y, y_expected), dif(z, z_expected)] .< tolerance
end


function lindblad_statedep(Γ, γ; tolerance = 0.001)

    ψ0 = normalize(g + e)
    tf = 10.0

    H = Iq
	γL(t::Timescale, ρ::State) = 0.5 * (Γ + 2γ^2 * (1 - expect(σz, ρ)^2))
	J = [(σz, γL)]
	C = []

    sol = bayesian((0.0, tf), ψ0, H, J, C; dt = 1e-3)
	
    x, y, z = map(op -> expectations(sol, op), [σx, σy, σz])

    return true
end


dif(x, y) = maximum(abs.(x .- y))

function test_lindblad(; tolerance = 0.001)

        @testset "lindblad constant decay" begin
            for γ in 0:0.5:2.0
                x, y, z = lindblad(γ; tolerance=tolerance)
                @test (x && y && z)
            end
        end

        @testset "lindblad time-dep decay" begin
            for γ in 0:0.5:1.0
                for Γ in 0:0.5:1.0
                    x, y, z = lindblad_timedep(Γ, γ; tolerance=tolerance)
                    @test (x && y && z)
                end
            end
        end

        @testset "lindblad state-dep (run)" begin
            @test lindblad_statedep(0.2, 0.5; tolerance=tolerance)
        end

end

function multilindblad(; γ1 = 0.1, γ2 = 0.2, γ3 = 0.3)

    H = Iq
    Γ1 = γ1
    Γ2(t::Timescale) = γ2
    Γ3(t::Timescale, ρ::State) = γ3

    J = [(σz, Γ1), (σx, Γ2), (σy, Γ3)]
    C = []
    tf = 10.0
    ψ0 = normalize(g + 0.3 * im * e)

    sol = bayesian((0.0, tf), ψ0, H, J, C; dt = 1e-3)
end