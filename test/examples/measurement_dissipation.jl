using QuantumCircuits
using QuantumCircuits.SingleQubitOperators
using QubitPlots

# define parameters of the simulation
Ω = 2π * 0.5
Γ = 0.2
η = 1.0

ψ0 = normalize(g + e) # initialize in the +x state
tf = 10.0
dt = 1e-3
N = 300

# set up Kraus operators
H = Ω/2 * σy
J = η < 1.0 ? [(σz, (1 - η) * Γ)] : []
C = [(σz, 2Γ, η)]

# solve an ensemble of trajectories, and find the corresponding η = 0 solution
ens = Ensemble(ensemble(bayesian, (0.0, tf), ψ0, H, J, C; dt=dt, N = N), qbasis, :qbasis)
sol0 = Solution(bayesian((0.0, tf), ψ0, H, [(σz, Γ)], []; dt=dt), qbasis, :qbasis)

# plot the ensemble of trajectories and  η = 0 trajectory along with the master equation solution
t = ens.t
x(t) = exp(-2Γ * t) # master equation solution for x(t)
plot(blochensemble, ens, sol0)
plot!(t, x.(t), subplot=1, color=:black, linestyle=:dash, linewidth=2)

# you can also just plot the ensemble average
m = mean.(ens.exps)
plot(blochtimeseries, t, m..., title="ensemble average of $N trajectories")
plot!(blochtimeseries, sol0, linestyle=:dash, label = :none)


"""
There's a mismatch between the master equation solution and the measurement
unless we include the 2Γ in front of the measurement rate. Not sure where
the discrepancy is from.
"""
