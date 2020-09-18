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
  dt = 1e-3
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
@testset "bayesian_step" begin
  Random.seed!(1)

  atol = 1e-5

  J0 = √Γ*σz
  C0 = √(Γ*η)*σz
  H0 = Ω*σy/2
  tt, ρt = @time bayesian((0, dt), ρ0, H0, [J0], [C0]; dt=dt)

  # compute first step by hand
  r = √(1/dt)*randn() + real(expect(ρ0, (C0 + C0')/2))
  MC = exp(DenseOperator(r*dt*C0/2 - dt*(C0 + C0')^2/16))
  MJ = J0
  U = exp( -im * dt * DenseOperator(H0))

  ρ_expected = U*MJ*MC*ρ0*MC'*MJ'*U'
  ρ_expected = ρ_expected/tr(ρ_expected)

  # compute error
  err = (ρt[1] - ρ_expected).data
  @test isapprox.(err, 0.0, atol=atol, rtol=0) == BitArray([1 1; 1 1])
end

# Positivity Test
@testset "bayesian_positivity" begin
  tt, ρt = @time bayesian((0, π), ρ0, σy, [σz], [√η*σz]);
  @test sum(ispositive.(map(e->e.data,ρt))) == length(tt)
end

# Ensemble Test
@testset "bayesian_ensemble" begin
  # Test Params
  target_max = 0.1

  # Time-dependent Hamiltonian
  @everywhere H = (Ω/2) * σy

  # Measurement dephasing
  @everywhere J = [√(Γ)σz]

  # Stochastic monitoring (quantum-limited efficiency)
  @everywhere C = [√(η*Γ)σz]

  # Deterministic expectation values
  t, evs = @time bayesian(T, ρ0, H, J, []; fn=expects, dt=dt);
  devs = hcat(evs...)'

  # Stochastic
  println("Sampling $n trajectories across $procs processes (estimate: ~20 minutes)...")
  t, ens = @time ensemble(bayesian, T, ρ0, H, J, C; fn=expects, dt=dt, N=n);
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

  println(maxx)
  println(maxy)
  println(maxz)

  @test maxx < target_max
  @test maxy < target_max
  @test maxz < target_max
end