"""
This brief tutorial assumes you want to save data and simulation parameters as part of your workflow.

Start by using QuantumCircuits, as well as operators for the system you'll want to simulation.
For commonly simulated systems, QuantumCircuits exports submodules (e.g. SingleQubitOperators,
TwoQubitOperators) with operators for those dimensions with standard names.

You can also import the submodule QubitPlots if you plan to plot anything.
"""

using QuantumCircuits
using QuantumCircuits.SingleQubitOperators
using QuantumCircuits.QubitPlots

using Parameters
using ProgressMeter

"""
Next, define parameters of the simulation. For later interpretation of your results, I recommend
creating a parameter dictionary which is saved in the same folder as your data.
"""

@with_kw struct Pars
    Ω = 2π * 0.0
    Γ = 0.2
    η = 1.0

    ψ0 = normalize(g + e) # initialize in the +x state
    tf = 10.0
    dt = 1e-3
    N = 300

end

"Wrap your simulation in a function to make it more easily reproducible."
function sim(pars::Pars)
    @unpack_Pars pars

    # Kraus operators
    H = Ω/2 * σy
    J = η < 1.0 ? [(σz, (1 - η) * Γ)] : []
    C = [(σz, 2Γ, η)]

    ens = Ensemble(ensemble(bayesian, (0.0, tf), ψ0, H, J, C; dt=dt, N = N, showprogress=false), qbasis, :qbasis)
    return ens.t, mean.(ens.exps)
end

"Now, it is relatively easy and clear to run simulations while varying only one parameter..."
Ωs = 2π .* (0.0:0.25:1.0)
means = []
times = []
@showprogress for Ω in Ωs
    t, m = sim(Pars(N = 250, Ω = Ω))
    times = t
    push!(means, m)
end

p = plot(xlabel="t (μs)", ylabel="x bloch coordinate", legendtitle = "Ω (MHz)")
for (m, Ω) in zip(means, Ωs)
    plot!(times, m[1], label=round(Ω/2π, digits=3))
end
p
