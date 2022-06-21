# Basis
if !@isdefined q

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
else
    println("Qubit basis q already defined; not re-including single-qubit operators.")
end
