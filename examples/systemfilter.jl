include("../src/QuantumCircuits.jl")
include("../plots/single_qubit_plots.jl")

using .QuantumCircuits
using .QuantumCircuits.SingleQubitOperators

#### Increase γ, ε, or td to see the system/filter discrepancy increase

using Distributions

function whitenoise(tf, dt)
    times = 0.0:dt:tf
    χtimeseries = rand(Normal(0,1), length(times))
    return t -> let i = findfirst(times .>= t)
                    χtimeseries[i]
                end
end

tf = 5.0
dt = 1e-3
Ω = 2π * 2
γ = 2π * 0.5
Γ2 = 0.2
ε = 0.1
Γm = 0.1
η = 1.0
td = 0.5

χ = whitenoise(tf, dt)

Hf(t::Timescale, ρ::State) = (Ω/2) * σx
Hs(t::Timescale, ρ::State) = Hf(t, ρ) + (γ/2) * χ(t) * σz


Js = [(σz, Γ2)]
Jf = [(σz, Γ2 - ε)]

C = [(σz, Γm, η)]

ψ0 = normalize(g + e)

sys, fil = bayesian(hamiltonianfe, (0.0, tf), ψ0, (Hs, Hf), (Js, Jf), C; dt = 1e-3, td = td, forwardestimate=false)
sys = Solution(sys, qbasis)
fil = Solution(fil, qbasis)

function plot_system_filter(sys::Solution, fil::Solution; layout=(2,1), size=(800,500), kwargs...)

    p1 = plot(blochtimeseries, sys.t, sys.exps...; title="system", kwargs...)
    p2 = plot(blochtimeseries, fil.t, fil.exps...; title = "filter", kwargs...)

    plot(p1, p2; layout=layout, size=size, kwargs...)
end


plot_system_filter(sys, fil; size=(600,500), legend=:outerright)

plot(blochtimeseries, sys.t, sys.exps...; legend=:outerright)
plot!(blochtimeseries, fil.t, fil.exps..., linewidth=4, linealpha=0.5, label=:none, legend=:outerright)