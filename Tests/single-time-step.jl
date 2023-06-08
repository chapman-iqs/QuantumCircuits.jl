function bayesian_update(Γ, η;  
                                basis=SpinBasis(1//2), 
                                r=0.2, 
                                ρ = QuantumOpticsBase.randstate(basis),
                                dt=1e-2, 
                                tolerance=1e-6,
                                kwargs...)

    # first, use solver...
    C = [(σz, Γ, η)]
    sol = bayesian((0.0, dt), ρ, Iq, [], C; dt=dt, records=[[r, r]])

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

    # return sol.ρ, ρ2

    d = if typeof(ρ) <: Ket
                tracedistance(dm(last(sol.ρ)), dm(ρ2))
            else
                tracedistance(last(sol.ρ), ρ2)
            end
    return d < tolerance
end

function test_single_timestep(; kwargs...)

    @testset "bayesian update" begin
        for η in 0.25:0.25:1.0
            for Γ in 0.5:0.5:1.5
                for r in rand(10)
                    @test bayesian_update(Γ, η; kwargs...)
                end
            end
        end     
    end
end