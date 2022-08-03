# # Single Step Test
# @test "bayesian_step" begin
#   Random.seed!(1)
#
#   using ..QuantumCircuits.SingleQubitOperators
#
#   atol = 1e-5
#
#   H = Ω * σy / 2
#   J = ()
#
#   J0 = √Γ*σz
#   C0 = √(Γ*η)*σz
#   H0 = Ω*σy/2
#   tt, ρt = @time bayesian((0, dt), ρ0, H0, [J0], [C0]; dt=dt)
#
#   # compute first step by hand
#   r = √(1/dt)*randn() + real(expect(ρ0, (C0 + C0')/2))
#   MC = exp(DenseOperator(r*dt*C0/2 - dt*(C0 + C0')^2/16))
#   MJ = J0
#   U = exp( -im * dt * DenseOperator(H0))
#
#   ρ_expected = U*MJ*MC*ρ0*MC'*MJ'*U'
#   ρ_expected = ρ_expected/tr(ρ_expected)
#
#   # compute error
#   err = (ρt[1] - ρ_expected).data
#   @test isapprox.(err, 0.0, atol=atol, rtol=0) == BitArray([1 1; 1 1])
# end
