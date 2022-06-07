cd("/Users/sachagreenfield/Desktop/GitHub/QuantumCircuits.jl")
import Pkg
Pkg.activate(".")

using QuantumCircuits
using Random
using Statistics
using Distributions
using DataFrames
using CSV


"Qubit Hilbert space operators ------------------------------------------------------------------------"

# Basis
q = SpinBasis(1//2)

# Operators, using convention that |-z> is ground state
σx = sigmax(q)
σy = sigmay(q)
σz = sigmaz(q)
σp = sigmap(q)
σm = sigmam(q)
id = identityoperator(q)

ground = spindown(q)
excited = spinup(q)


# System parameters ------------------------------------------------------------
# all times given in μs

# initial state
(x0,y0,z0) = (0., 0.3, 0.91)
ρ0 = DenseOperator(0.5*(id + x0*σx + y0*σy + z0*σz))
#
#
# # measurement parameters
#
# τm = parse(Float64, τmstr)			# time
# Γm = 1/(2τm) 						# rate
# η =  parse(Float64, ηstr)			# efficiency
# φ = get(φdict, φIB, 0) 				# angle
# td = parse(Float64, Tdstr) * 1e-3 	# time delay for feedback
#
#
# # simulation timescales
#
# T = (0, 8τm) # simulation duration
# dt = 0.5e-3  # integration time-step
#
#
# # feedback target parameters
#
# θs = 3π/10 		# target angle on Bloch sphere
# Rs = 1 # 0.64 			# radius of target
# ϕ = π 					# fixes plane of oscillations
# vec = (sin(θs) * cos(π/2), sin(θs) * sin(π/2), cos(θs))
# σϕ = cos(ϕ)*σx + sin(ϕ)*σy
#
#
# # feedback drive parameters
#
# Δ0 = ideal ? -sin(2θs)/(4τm) : -sin(2θs)/(4τm*Rs^2)
# Δ1 = ideal ? sin(θs)/τm : sin(θs)/(τm*Rs)
#
#
# # decay times and rates
#
# T1 = 40 	# energy decay time
# Γ1 = 1/(2T1)# energy decay rate
# T2 = 60 	# environmental dephasing time
# Γ2 = 1/T2 	# environemntal dephasing rate
#
#
# # Kraus operators -------------------------------------------------------------
#
# H = Hf(t, r) = (Δ0 + Δ1*r[1]) * σϕ/2
# J = ideal ? [(σz, ((1-η)*Γm))] : [(σz, ((1-η)*Γm)), (σm, Γ1), (σz, Γ2)]
# C = [(exp(im * φ) * σz, τm, η)]
#
#
# # Bayesian simulation ---------------------------------------------------------
# trajs = []
#
# # loop over N simulations
# for i in 1:N
# 	sol = bayesian(T, ρ0, H, J, C; dt=dt, td=td)
# 	tr = traj(sol)
# 	push!(trajs, tr)
# end
#
# # 	xm, ym, zm = (mean(xs), mean(ys), mean(zs))
#
# # 	p = 0.5 .* (1 .+ xm.^2 + ym.^2 + zm.^2) # state purity
# t = collect(T[1]:dt:T[2])   # list of times
#
# global dfx = DataFrame([tr.x for tr in trajs], [string("x", i) for i in 1:N])
# global dfy = DataFrame([tr.y for tr in trajs], [string("y", i) for i in 1:N])
# global dfz = DataFrame([tr.z for tr in trajs], [string("z", i) for i in 1:N])
# global dfr = DataFrame([tr.r[1] for tr in trajs], [string("r", i) for i in 1:N])
# global dft = DataFrame([t], ["t"])
#
# # global ILFens = traj(t, xm, ym, zm, p, ())
