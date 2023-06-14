function bayesian_update(Γ, η;  solve = bayesian,
                                r=0.2, 
                                ρ = QuantumOpticsBase.randstate(SpinBasis(1//2)),
                                dt=1e-2, 
                                tolerance=1e-6,
                                kwargs...)

    # first, use solver...
    C = [(σz, Γ, η)]
    sol = solve((0.0, dt), ρ, Iq, [], C; dt=dt, records=[[r, r]])

    # do by hand for single-qubit σz
    τ = 1 / (2Γ * η)
    Pg = 1/sqrt(2π * τ * dt) * exp(-dt * (r + 1)^2 / 2τ)
    Pe = 1/sqrt(2π * τ * dt) * exp(-dt * (r - 1)^2 / 2τ)
    M = sqrt(Pg) * dm(g) + sqrt(Pe) * dm(e)

    # ρ1 = M * ρ
    # ρ2 = ρ1 / sqrt(ρ1' * ρ1)
    ρ2 =    if typeof(ρ) <: Ket
                ρ1 = M * ρ
                ρ1 / sqrt(ρ1' * ρ1)
            else
                ρ1 = M * ρ * M' 
                ρ1 / tr(ρ1)
            end

    # trace distance
    ρf = dm(last(sol.ρ))
    ρ2 = dm(ρ2)
    # ρf = last(sol.ρ)
    # ρf = isa(ρf, Ket) ? dm(ρf) : ρf
    # ρ2 = isa(ρ2, Ket) ? dm(ρ2) : ρ2

    d = tracedistance(ρf, ρ2)
    return d < tolerance
end

function lindblad_Γ2_decay(Γ;
                                solve = bayesian,
                                ψ = QuantumOpticsBase.randstate(SpinBasis(1//2)), 
                                dt=1e-2, 
                                tolerance=1e-4,
                                kwargs...)

    # first, use solver...
    J = [(σz, Γ/2)]
    sol = solve((0.0, dt), ψ, Iq, J, []; dt=dt)

    # do by hand for single-qubit σz dephasing
    ρ = dm(ψ)
    ρ.data[1,2] = ρ.data[1,2] * exp(-Γ * dt)
    ρ.data[2,1] = ρ.data[2,1] * exp(-Γ * dt)

    # compare via trace distance
    d = tracedistance(dm(last(sol.ρ)), dm(ρ))
    return d < tolerance
end

function ham_update(Ω;  solve = bayesian,
                        ψ = QuantumOpticsBase.randstate(SpinBasis(1//2)), 
                        dt=1e-2, 
                        tolerance=1e-4,
                        kwargs...)

    # first, use solver...
    H = (Ω/2)σy
    sol = solve((0.0, dt), ψ, H, [], []; dt=dt)

    # do by hand for small dt
    U = Iq - im * H * dt
    ρ = U * dm(ψ) * U'
    ρ = ρ / tr(ρ)

    # compare via trace distance
    d = tracedistance(dm(last(sol.ρ)), ρ)
    return d < tolerance
end

function test_single_timestep(; solve=bayesian, kwargs...)
    @testset "bayesian update (single time step)" begin
        if solve == bayesian
            for η in 0.25:0.25:1.0
                for Γ in 0.5:0.5:1.5
                    for r in rand(5)
                        @test bayesian_update(Γ, η; solve=solve, kwargs...)
                    end
                end
            end   
        else
            println("⬇︎⬇︎⬇︎ Measurement unit test not yet written for rouchon method.")
        end  
    end



    @testset "lindblad decay (single time step)" begin
        for Γ in 0.5:0.5:1.5
            @test lindblad_Γ2_decay(Γ; solve=solve, kwargs...)
        end
    end

    @testset "hamiltonian rabi (single time step)" begin
        for Ω in 2π .* (0.5:0.5:1.5)
            @test ham_update(Ω; solve=solve, kwargs...)
        end
    end
end