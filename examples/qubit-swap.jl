include("../src/QuantumCircuits.jl")

using .QuantumCircuits
using .QuantumCircuits.SingleQubitOperators
using Plots, Measures
include("../plots/single_qubit_plots.jl")


####### operators
# cavity
N = 4
const f = FockBasis(N) 
const If = identityoperator(f) # identity
const af = destroy(f)       # annihilation operator
const nf = af' * af         # number operator
const v = fockstate(f, 0) # vacuum state
fbasis = map(n -> fockstate(f,n), [0,1,2,3,4])

const σzq = σz ⊗ If ⊗ Iq  # qubit z
const σxq = σx ⊗ If ⊗ Iq  # qubit x
const σzd = Iq ⊗ If ⊗ σz # dissipator z
const σxd = Iq ⊗ If ⊗ σx # dissipator x
const a = Iq ⊗ af ⊗ Iq # resonator a
const n = Iq ⊗ nf ⊗ Iq # resonator n
const I = Iq ⊗ If ⊗ Iq # identity

####### parameters (reference scale: 1 GHz)
ωq = 2π * 4
ωd = 2π * 9
ωc = 2π * 6
ωm = ωc - ωq
gqc = 2π * 0.05
gqd = 2π * 0.1
Δ = 2π * 0.5

dt = (2π/ωd) / 8
tf = 5000 * dt

ψq = e # qubit is excited
ψc = fockstate(f, 0)# cavity is vacuum
ψd = g # dissipator is ground
ψ0 = ψq ⊗ ψc ⊗ ψd

####### Hamiltonian
Hqubit = -(ωq/2)σzq
Hdiss(t) = (-ωd/2 + (Δ/2)cos(ωm * t)) * σzd 
Hcav = ωc * n
Hint = gqd*(a' + a)*σxd + gqc*(a' + a)*σxq

H(t) =  Hcav #+ Hdiss(t) + Hint + Hqubit


####### solve

sol = bayesian((0.0, tf), ψ0, H, [], []; dt=dt)

ρq = map(ρ -> ptrace(ρ, [2,3]), sol.ρ) # qubit solution
solq = Solution(Solution(sol.t, ρq), qbasis) 
plot(solq.t, solq.exps[3])

ρd = map(ρ -> ptrace(ρ, [1,2]), sol.ρ) # dissipator solution
sold = Solution(Solution(sol.t, ρd), qbasis) 
plot(sold.t, solq.exps[3])

ρc = map(ρ -> ptrace(ρ, [1,3]), sol.ρ) # cavity solution
solc = Solution(Solution(sol.t, ρc), fbasis)
plot(sol.t, solc.exps[1])
###### let's see if we got swaps
# plot(blochtimeseries, sold, legend=:outerright, xlabel="t (ns)")

