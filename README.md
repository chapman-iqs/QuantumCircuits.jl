# QuantumCircuits.jl

**QuantumCircuits.jl** is a numerical framework written in Julia for simulating quantum evolution, with a focus on superconducting circuits, and builds on [QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl).

## Usage

```jl
using QuantumOptics
q = SpinBasis(1//2)
ρ0 = dm(spindown(q))
J(t) = [√Γ*σz]
C(t) = [√(Γ*η)*σz]
H(t) = [Ω*σy/2]

using QuantumCircuits
tt, ρs, dy = rouchon(tspan, ρ0, H, J, C; dt=1e-4, dy=dy)
tt, ρs, dy = bayesian(tspan, ρ0, H, J, C; dt=1e-4, dy=dy)
tt, ρs, dy = ensemble(method, tspan, ρ0, H, J, C; N=10, dt=1e-4, dy=dy)
```

[Efficient Quantum Filtering for Quantum Feedback Control](https://arxiv.org/abs/1410.5345) (Pierre Rouchon, Jason F. Ralph)