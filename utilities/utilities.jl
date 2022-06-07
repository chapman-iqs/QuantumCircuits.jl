"Trajectory utilities ------------------------------------------------"
mutable struct traj
  t::Vector{Float64}
  x::Vector{Float64}
  y::Vector{Float64}
  z::Vector{Float64}
  p::Vector{Float64}
  r
end

function traj(sol::Solution)
  t, ρ, r = (sol.t, sol.ρ, sol.r)
  x, y, z = map(op -> expectations(sol, op), [σx, σy, σz])
  p = (typeof(ρ[1]) <: Ket) ? [1.0 for el in ρ] : real(expect.(ρ, ρ))
  traj(t, x, y, z, p, r)
end
