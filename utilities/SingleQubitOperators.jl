module SingleQubitOperators

" ---- Export symbols ---- "

export q, Iq, σx, σy, σz, σp, σm, g, e, qbasis, qlabels, qcolors


" ---- Dependencies ---- "

using QuantumCircuits
import Plots.palette


" ---- Definitions of constants ---- "

const q = SpinBasis(1//2)
const Iq = identityoperator(q)

# qubit operators, using convention that |-z> (spin down) is ground state
const σx = sigmax(q)
const σy = sigmay(q)
const σz = sigmaz(q)
const σp = sigmap(q)
const σm = sigmam(q)

# qubit states
const g = spindown(q)
const e = spinup(q)

# operator basis
const qbasis = [σx, σy, σz]
const qlabels = ["x", "y", "z"]
const qcolors =  palette(:tab10)[1:3]


end # module SingleQubitOperators
