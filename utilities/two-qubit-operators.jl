"Single-qubit operators ------------------------------------------------"

# Basis
q = SpinBasis(1//2)
Iq = identityoperator(q)

# qubit operators, using convention that |-z> is ground state (aligned with QuantumOpticsBase convention)
σx = sigmax(q)
σy = sigmay(q)
σz = sigmaz(q)
σp = sigmap(q)
σm = sigmam(q)
nq = σp * σm

g = spindown(q)
e = spinup(q)

"Two-qubit operators ------------------------------------------------"
I = Iq ⊗ Iq

σx1 = σx ⊗ Iq
σx2 = Iq ⊗ σx
σy1 = σy ⊗ Iq
σy2 = Iq ⊗ σy
σz1 = σz ⊗ Iq
σz2 = Iq ⊗ σz
σp1 = σp ⊗ Iq
σm1 = σm ⊗ Iq
σp2 = Iq ⊗ σp
σm2 = Iq ⊗ σm

# number operators
n1 = nq ⊗ Iq
n2 = Iq ⊗ nq
n = n1 + n2

# basis states
ket00 = g ⊗ g
ket01 = g ⊗ e
ket10 = e ⊗ g
ket11 = e ⊗ e

Φp = normalize(g ⊗ g + e ⊗ e)
Φm = normalize(g ⊗ g - e ⊗ e)
Ψp = normalize(g ⊗ e + e ⊗ g)
Ψm = normalize(g ⊗ e - e ⊗ g)


"Two-qubit labels ------------------------------------------------"

const bell_basis = [Ψp, Ψm, Φp, Φm]
const bell_basis_strings = [ "Ψp", "Ψm", "Φp", "Φm"]
const bell_basis_labels = [L"| \Psi_+ \rangle",
				L"| \Psi_- \rangle",
				L"| \Phi_+ \rangle",
				L"| \Phi_- \rangle"]
const number_basis = [ket00, Ψp, Ψm, ket11]
const number_basis_strings = [ "00", "Ψp", "Ψm", "11"]
const number_basis_labels = [L"| 00 \rangle",
								 L"| \Psi_+ \rangle",
								 L"| \Psi_- \rangle",
								 L"| 11 \rangle"]
