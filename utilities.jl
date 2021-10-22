mutable struct traj
  t::Vector{Float64}
  x::Vector{Float64}
  y::Vector{Float64}
  z::Vector{Float64}
  p::Vector{Float64}
  r
end

function traj(t, ρ, r)
  x, y, z = [real(expect(σi, ρ)) for σi in (σx, σy, σz)]
  p = real(expect.(ρ, ρ))
  traj(t, x, y, z, p, r)
end

function traj(sol::QuantumCircuits.solution)
  t, ρ, r = (sol.t, sol.ρ, sol.r)
  x, y, z = [real(expect(σi, ρ)) for σi in (σx, σy, σz)]
  p = real(expect.(ρ, ρ))
  traj(t, x, y, z, p, r)
end
