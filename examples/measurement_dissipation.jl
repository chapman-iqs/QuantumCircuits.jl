include("../src/QuantumCircuits.jl")

using .QuantumCircuits
using .QuantumCircuits.SingleQubitOperators
using Plots, Measures
include("../plots/single_qubit_plots.jl")
include("../plots/ensembleplots.jl")

# define parameters of the simulation
Ω = 2π * 1.0
Γ = 0.2
η = 1.0

ψ0 = normalize(g + e) # initialize in the +x state
tf = 10.0
dt = 1e-3
N = 300

# set up Kraus operators
H = Ω/2 * σy
J = η < 1.0 ? [(σz, (1 - η) * Γ/2)] : [] # here, we use Γ/2 because of the relation between bloch coordinate dephasing and density matrix dephasing for a single qubit
C = [(σz, Γ, η)]

# solve an ensemble of trajectories, and find the corresponding η = 0 solution
ens = Ensemble(ensemble(bayesian, (0.0, tf), ψ0, H, J, C; dt=dt, N = N), qbasis, :qbasis);
sol0 = Solution(bayesian((0.0, tf), ψ0, H, [(σz, Γ/2)], []; dt=dt), qbasis, :qbasis);

# plot the ensemble of trajectories and  η = 0 trajectory along with the master equation solution
t = ens.t
x(t) = exp(-Γ * t) # master equation solution for x(t)

plot(blochensemble, ens, sol0; legend=:outerright, size=(1000,500), bottommargin=5mm)
plot!(t, x.(t), subplot=1, color=:black, linestyle=:dash, linewidth=2)

# you can also just plot the ensemble average
m = mean.(ens.exps)
plot(blochtimeseries, t, m..., title="ensemble average of $N trajectories")
plot!(blochtimeseries, sol0, linestyle=:dash, label = :none, legend=:outerright)
