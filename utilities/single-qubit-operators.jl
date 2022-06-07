"Single-qubit operators ------------------------------------------------"

# Basis
q = SpinBasis(1//2)
Iq = identityoperator(q)

# qubit operators, using convention that |-z> is ground state
σx = sigmax(q)
σy = sigmay(q)
σz = sigmaz(q)
σp = sigmap(q)
σm = sigmam(q)
n = σp * σm

g = spindown(q)
e = spinup(q)

# operator basis
qbasis = [σx, σy, σz]
qlabels = ["x", "y", "z"]
