include("../src/QuantumCircuits.jl")

using .QuantumCircuits
using .QuantumCircuits.SingleQubitOperators

Hs(t::Timescale, ρ::State) = (2π) * 0.1 * σx
Hf(t::Timescale, ρ::State) = Iq

Js = [(σz, 0.5)]
Jf = []

C = [(σz, 0.1, 1.0)]

ψ0 = normalize(g + e)

sys, fil = bayesian((0.0, 5.0), ψ0, (Hs, Hf), (Js, Jf), C; dt = 1e-3)
sys = Solution(sys, qbasis)
fil = Solution(fil, qbasis)

using QubitPlots

plot(blochtimeseries, sys.t, sys.exps...)
