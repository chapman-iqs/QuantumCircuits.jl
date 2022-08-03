"""
Still needs to be updated to module format, as in SingleQubitOperators.jl and TwoQubitOperators.jl.
"""

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

"Three-qubit operators ------------------------------------------------"
I = Iq ⊗ Iq ⊗ Iq

σx1 = σx ⊗ Iq ⊗ Iq
σx2 = Iq ⊗ σx ⊗ Iq
σx3 = Iq ⊗ Iq ⊗ σx
σy1 = σy ⊗ Iq ⊗ Iq
σy2 = Iq ⊗ σy ⊗ Iq
σy3 = Iq ⊗ Iq ⊗ σy
σz1 = σz ⊗ Iq ⊗ Iq
σz2 = Iq ⊗ σz ⊗ Iq
σz3 = Iq ⊗ Iq ⊗ σz
n1 = nq ⊗ Iq ⊗ Iq
n2 = Iq ⊗ nq ⊗ Iq
n3 = Iq ⊗ Iq ⊗ nq
n = n1 + n2 + n3

"Basis states and labels ------------------------------------------------"
# product states
ket000 = g ⊗ g ⊗ g
ket001 = g ⊗ g ⊗ e
ket010 = g ⊗ e ⊗ g
ket100 = e ⊗ g ⊗ g
ket110 = e ⊗ e ⊗ g
ket101 = e ⊗ g ⊗ e
ket011 = g ⊗ e ⊗ e
ket111 = e ⊗ e ⊗ e
basis_3q = [ket000, ket001, ket010, ket100, ket110, ket101, ket011, ket111]

basis_3q_labels = [L"|000\rangle",
					L"|001\rangle",
					L"|010\rangle",
					L"|100\rangle",
					L"|110\rangle",
					L"|101\rangle",
					L"|011\rangle",
					L"|111\rangle"]

# number eigenstates
ket0 = ket000
ket1(α, β, γ) = normalize(α * ket001 + β * ket010 + γ * ket100)
ket2(α, β, γ) = normalize(α * ket100 + β * ket101 + γ * ket011)
ket3 = ket111

# number POVM
Λ0 = dm(ket000)
Λ1 = dm(ket001) + dm(ket010) + dm(ket100)
Λ2 = dm(ket110) + dm(ket101) + dm(ket011)
Λ3 = dm(ket111)

number_POVM = [Λ0, Λ1, Λ2, Λ3]
number_labels = [L"Tr (\rho \Lambda_0)",
					L"Tr (\rho \Lambda_1)",
					L"Tr (\rho \Lambda_2)",
					L"Tr (\rho \Lambda_3)"]
