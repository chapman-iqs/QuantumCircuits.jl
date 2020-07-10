using Test
using QuantumCircuits
using Statistics

@testset "rouchon" begin

# Test Params
target_max = 0.1

# Parameters
n = 100 # ensemble size
Ω  = 2π # Rabi frequency
τ = 3.0 # Measurement collapse timescale
Γ = 1/(2τ) # Measurement dephasing rate
dt = 1e-4
T = (0.0, 6τ) # Time duration of simulation
η = 0.3

# Basis
q = SpinBasis(1//2)

# Operators
σx = sigmax(q)
σy = sigmay(q)
σz = sigmaz(q)

# Initial condition
ρ0 = dm(spindown(q))

# Time-dependent Hamiltonian
H(t) = 2exp(-(t-3τ)^2/2)/sqrt(2π) * (Ω/2) * σy

# Measurement dephasing
# J = [sqrt(Γ/2)σz]
J = [sqrt(Γ)σz]

# Stochastic monitoring (quantum-limited efficiency)
C = [sqrt(η*Γ)σz]

# Solve using Rouchon method

expects = ρ -> collect(real(expect(ρ, s)) for s in [σx,σy,σz]) # ρ -> [<x>,<y>,<z>]

# Deterministic expectation values
t, devs = @time rouchon(T, ρ0, H, J, []; fn=expects, dt=dt);
devs = hcat(devs...)'

# Stochastic
t, evals = @time ensemble(rouchon, T, ρ0, H, J, C; fns=[expects], dt=dt, N=n);

sevsx = []
sevsy = []
sevsz = []

for i in 1:n
  evs = hcat(evals[1,i,:]...)'
  
  # Compute average stochastic expectation values
  if (i == 1)
    sevsx = evs[:,1]
    sevsy = evs[:,2]
    sevsz = evs[:,3]
  else
    sevsx += evs[:,1]
    sevsy += evs[:,2]
    sevsz += evs[:,3]
  end
end

sevsx /= n
sevsy /= n
sevsz /= n

# Difference between deterministic and average stochastic expectation values
diffx = devs[:,1]-sevsx
diffy = devs[:,2]-sevsy
diffz = devs[:,3]-sevsz

# Naive convergence test
maxx = maximum(diffx)
maxy = maximum(diffy)
maxz = maximum(diffz)

@test maxx < target_max
@test maxy < target_max
@test maxz < target_max
end