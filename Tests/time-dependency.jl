using .QuantumCircuits.SingleQubitOperators
using QubitPlots

Ω(t) = floor(t) % 4 == 0 ? (2π) * 1.0 : 0.0
H(t::Timescale) = Ω(t) * σx
ψ0 = g

sol = bayesian((0.0, 10.0), ψ0, H, [], []; dt = 1e-3)
sol = Solution(sol, qbasis)


plot(blochtimeseries, sol.t, sol.exps...)
