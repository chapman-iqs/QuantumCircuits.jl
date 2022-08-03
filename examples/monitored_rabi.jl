using QuantumCircuits
using QuantumCircuits.SingleQubitOperators
using QubitPlots
# using QuantumCircuits.QubitPlots

Ω = 2π * 1.0
Γ = 2.0
η = 1.0

ψ0 = normalize(g + e)
tf = 10.0
dt = 1e-3

H = Ω/2 * σy
J = η < 1.0 ? [(σz, (1 - η) * Γ)] : []
C = [(σz, Γ, η)]

sol = Solution(bayesian((0.0, tf), ψ0, H, J, C; dt=dt), qbasis, :qbasis)
plot(blochtimeseries, sol, legend=:outerright)
