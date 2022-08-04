using ProgressMeter
using .QuantumCircuits
using .QuantumCircuits.SingleQubitOperators

Ω = 2π * 1.0
Γ = 0.2
dt = 1e-3
tf = 10.0
ρ0 = normalize(g + e)

H = Ω * σy / 2
J(η = 0) = η == 1.0 ? [] : [(σz, (1 - η) * Γ)]
C(η = 0) = η > 0.0 ? [(σz, Γ, η)] : []

# get solutions and ensembles
sol0 = bayesian((0, tf), ρ0, H, J(η = 0), C(η = 0); dt=dt)
ens1 = @showprogress ensemble(bayesian, (0, tf), ρ0, H, J(η = 1.0), C(η = 1.0); dt=dt, N = 1000)
ens05 = @showprogress ensemble(bayesian, (0, tf), ρ0, H, J(η = 0.5), C(η = 0.5); dt=dt, N = 1000)

# get expectation values in pauli basis
sol0 = Solution(sol0, qbasis)
ens1 = Ensemble(ens1, qbasis)
ens05 = Ensemble(ens05, qbasis)
