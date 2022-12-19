import .QuantumCircuits: supdate, update, ham, sham, lind, slind, meas, smeas

"""
A series of checks that the superoperator functions match what we get from sandwiching operations.
"""

function superops(; basis=SpinBasis(1//2), tolerance = 1e-9, verbose = true, kwargs...)

    # generate a random state and random operator for testing
    ρ = dm(QuantumOpticsBase.randstate(basis))
    A = QuantumOpticsBase.randoperator(basis)

    superprod = supdate(A, ρ)
    sandprod = update(A, ρ)

    d = tracedistance(superprod, sandprod)
    bool = d < tolerance
    if verbose
        println("Superoperator update trace distance $d $(bool ? "<" : ">") $tolerance")
    end

    return bool
end

function superham(; basis=SpinBasis(1//2), tolerance = 1e-9, verbose = true, dt=1e-3, kwargs...)

        # generate a random state and random Hermitian operator for testing
        ρ = dm(QuantumOpticsBase.randstate(basis))
        A = QuantumOpticsBase.randoperator(basis)
        H = A + A'
    
        ham1::SuperOperator = sham(dt, H)
        ham2::Function = ham(dt, H) 

        superprod = supdate(ham1, ρ)
        sandprod = ham2(0.0, ρ)
    
        d = tracedistance(superprod, sandprod)
        bool = d < tolerance
        if verbose
            println("Superoperator Hamiltonian update trace distance $d $(bool ? "<" : ">") $tolerance")
        end
    
        return bool
end

function superlind(; basis=SpinBasis(1//2), tolerance = 1e-9, verbose = true, dt=1e-3, kwargs...)

        # generate a random state and random operator for testing
        ρ = dm(QuantumOpticsBase.randstate(basis))
        A = QuantumOpticsBase.randoperator(basis)
        γ = rand()
        J = [(A, γ)]
    
        lind1::Function = slind(dt, J)
        lind2::Function = lind(dt, J) 

        superprod = supdate(lind1(0.0, ρ), ρ)
        sandprod = lind2(0.0, ρ)
    
        d = tracedistance(superprod, sandprod)
        bool = d < tolerance
        if verbose
            println("Superoperator Lindblad update trace distance $d $(bool ? "<" : ">") $tolerance")
        end
    
        return bool
end


function supermeas(; basis=SpinBasis(1//2), tolerance = 1e-9, verbose = true, dt=1e-3, kwargs...)

        ρ = dm(QuantumOpticsBase.randstate(basis))
        A = QuantumOpticsBase.randoperator(basis)
        γ = rand()
        η = rand()
        C = [(A, γ, η)]
    
        meas1::Function = smeas(dt, C)
        meas2::Function = meas(dt, C; sample=false)
        
        M, ro = meas1(0.0, ρ)

        superprod = normalize(supdate(M, ρ))
        sandprod = meas2(0.0, ρ, ro)
    
        d = tracedistance(superprod, sandprod)
        bool = d < tolerance
        if verbose
            println("Superoperator Measurement update trace distance $d $(bool ? "<" : ">") $tolerance")
        end
    
        return bool
end


function test_superops(; n = 10, kwargs...)

    @testset "Superoperator Measurement update" begin
        for i in 1:n
            @test supermeas(; kwargs...)
        end
    end


    @testset "Superoperator update" begin
        for i in 1:n
            @test superops(; kwargs...)
        end
    end

    @testset "Superoperator Hamiltonian update" begin
        for i in 1:n
            @test superham(; kwargs...)
        end
    end

    @testset "Superoperator Lindblad update" begin
        for i in 1:n
            @test superlind(; kwargs...)
        end
    end


end
