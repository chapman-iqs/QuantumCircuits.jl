"Three-qubit operators ------------------------------------------------"
const I3 = Iq ⊗ Iq ⊗ Iq

const σx1_3q = σx ⊗ Iq ⊗ Iq
const σx2_3q = Iq ⊗ σx ⊗ Iq
const σx3_3q = Iq ⊗ Iq ⊗ σx
const σy1_3q = σy ⊗ Iq ⊗ Iq
const σy2_3q = Iq ⊗ σy ⊗ Iq
const σy3_3q = Iq ⊗ Iq ⊗ σy
const σz1_3q = σz ⊗ Iq ⊗ Iq
const σz2_3q = Iq ⊗ σz ⊗ Iq
const σz3_3q = Iq ⊗ Iq ⊗ σz
const n1_3q = nq ⊗ Iq ⊗ Iq
const n2_3q = Iq ⊗ nq ⊗ Iq
const n3_3q = Iq ⊗ Iq ⊗ nq
const n_3q = n1 + n2 + n3

"Basis states and labels ------------------------------------------------"
# product states
const ket000 = g ⊗ g ⊗ g
const ket001 = g ⊗ g ⊗ e
const ket010 = g ⊗ e ⊗ g
const ket100 = e ⊗ g ⊗ g
const ket110 = e ⊗ e ⊗ g
const ket101 = e ⊗ g ⊗ e
const ket011 = g ⊗ e ⊗ e
const ket111 = e ⊗ e ⊗ e
const basis_3q = [ket000, ket001, ket010, ket100, ket110, ket101, ket011, ket111]

const basis_3q_labels = [L"|000\rangle",
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
const Λ0 = dm(ket000)
const Λ1 = dm(ket001) + dm(ket010) + dm(ket100)
const Λ2 = dm(ket110) + dm(ket101) + dm(ket011)
const Λ3 = dm(ket111)

const number_POVM = [Λ0, Λ1, Λ2, Λ3]
const number_labels = [L"Tr (\rho \Lambda_0)",
					L"Tr (\rho \Lambda_1)",
					L"Tr (\rho \Lambda_2)",
					L"Tr (\rho \Lambda_3)"]
