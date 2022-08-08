using ProgressMeter
using Test
using .QuantumCircuits
using .QuantumCircuits.SingleQubitOperators
using Statistics

Ω = 2π * 1.0
Γ = 0.2
dt = 1e-3
tf = 10.0
ρ0 = normalize(g + e)
N = 5000
η1 = 1.0
η2 = 0.5

H = Ω * σy / 2
J(η) = η == 1.0 ? [] : [(σz, (1 - η) * Γ)]
C(η) = η > 0.0 ? [(σz, 2Γ, η)] : []

# get solutions and ensembles
sol0 = bayesian((0, tf), ρ0, H, J(0.0), C(0.0); dt=dt)
ens1 = ensemble(bayesian, (0, tf), ρ0, H, J(η1), C(η1); dt=dt, N = N)
ens2 = ensemble(bayesian, (0, tf), ρ0, H, J(η2), C(η2); dt=dt, N = N)

# get expectation values in pauli basis
sol0 = Solution(sol0, qbasis)
ens1 = Ensemble(ens1, qbasis)
ens2 = Ensemble(ens2, qbasis)

t, x, y, z = sol0.t, sol0.exps...
t1, x1, y1, z1 = ens1.t, mean.(ens1.exps)...
t2, x2, y2, z2 = ens2.t, mean.(ens2.exps)...

# the most a bloch coordinate value can differ from another is 2.0, since it is
# within the range [-1,1]. The convergence accurancy scales as 1/√N with the
# number of trajectories N. So a reasonable error to aim for is 2.0 / √N .
# Note that for small values of N (< 1000), statistical fluctuations may cause the target
# error to be exceeded even if it is functioning correctly. So it is recommended to
# test for N >= 1000.

const targeterror = 2.0 / √N
maxerr(a, a1) = maximum(abs.(a .- a1))

@testset "bloch_convergence" begin
    @test maxerr(x, x1) < targeterror * η1
    @test maxerr(y, y1) < targeterror * η1
    @test maxerr(z, z1) < targeterror * η1

    @test maxerr(x, x2) < targeterror * η2 # not sure if this is the right scaling factor with η2. Maybe it should be √η ?
    @test maxerr(y, y2) < targeterror * η2
    @test maxerr(z, z2) < targeterror * η2
end
