using Distributed
using Test
using Statistics
using Distributions
using LinearAlgebra
using Random

procs = 10
addprocs(procs)

@everywhere begin
  using QuantumCircuits

  # Parameters
  n = 150 # ensemble size
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
  id = identityoperator(q)

  # Initial condition
  ρ0 = dm(spindown(q))

  # Helper
  expects = ρ -> collect(real(expect(ρ, s)) for s in [σx,σy,σz]) # ρ -> [<x>,<y>,<z>]
  ispositive = ρ -> begin
      mineig = eigmin(real(ρ))
      if mineig > 0 || abs(mineig) < 1e-10
          return true
      else 
          return false
      end
  end
end

# Single Step Test
  @testset "rouchon_step" begin
  Random.seed!(1)

  atol = 1e-6

  J0 = √Γ*σz
  C0 = √(Γ*η)*σz
  H0 = Ω*σy/2
  tt, ρt = @time rouchon((0, dt), ρ0, H0, [J0], [C0]; dt=dt)
  Δρ = ρt[1] - ρ0;

  # compute first step by hand
  dist = Normal(0, sqrt(dt))
  dW = rand(dist,1)[1]
  Δρ_expected = -im*(H0*ρ0 - ρ0*H0)*dt + (J0*ρ0*J0' - 0.5(J0'*J0*ρ0 + ρ0*J0'*J0))*dt + (C0*ρ0 + ρ0*C0' - tr(C0*ρ0 + ρ0*C0')*ρ0)*dW

  # compute error
  err = (Δρ - Δρ_expected).data

  # test
  @test isapprox.(err, 0.0, atol=atol, rtol=0) == BitArray([1 1; 1 1])
end

# Positivity Test
@testset "rouchon_positivity" begin
  tt, ρt = @time rouchon((0, π), ρ0, σy, [σz], [√η*σz]);
  @test sum(ispositive.(map(e->e.data,ρt))) == length(tt)
end

# Ensemble Test
@testset "rouchon_ensemble" begin

# Test Params
target_max = 0.1

# Time-dependent Hamiltonian
@everywhere H = (Ω/2) * σy

# Measurement dephasing
@everywhere J = [√(Γ)σz]

# Stochastic monitoring (quantum-limited efficiency)
@everywhere C = [√(η*Γ)σz]

# Deterministic expectation values
t, evs = @time rouchon(T, ρ0, H, J, []; fn=expects, dt=dt);
devs = hcat(evs...)'

# Stochastic
println("Sampling $n trajectories across $procs processes (estimate: ~3 minutes)...")
t, ens = @time ensemble(rouchon, T, ρ0, H, J, C; fn=expects, dt=dt, N=n);
ens = [ens[m,n][ax] for n in 1:length(t), m in 1:n, ax in 1:3]
avg = hcat(mean(ens[:,:,1], dims=2), mean(ens[:,:,2], dims=2), mean(ens[:,:,3], dims=2));

# Difference between deterministic and average stochastic expectation values
diffx = devs[:,1]-avg[:,1]
diffy = devs[:,2]-avg[:,2]
diffz = devs[:,3]-avg[:,3]

# Naive convergence test
maxx = maximum(diffx)
maxy = maximum(diffy)
maxz = maximum(diffz)

@test maxx < target_max
@test maxy < target_max
@test maxz < target_max
end