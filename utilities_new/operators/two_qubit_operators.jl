using LaTeXStrings

if !@isdefined I2
	"Two-qubit operators ------------------------------------------------"
	const I2 = Iq ⊗ Iq

	const σx1 = σx ⊗ Iq
	const σx2 = Iq ⊗ σx
	const σy1 = σy ⊗ Iq
	const σy2 = Iq ⊗ σy
	const σz1 = σz ⊗ Iq
	const σz2 = Iq ⊗ σz
	const σp1 = σp ⊗ Iq
	const σm1 = σm ⊗ Iq
	const σp2 = Iq ⊗ σp
	const σm2 = Iq ⊗ σm

	# number operators
	const n1 = σp1 * σm1
	const n2 = σp2 * σm2
	const n = n1 + n2

	# basis states
	const ket00 = g ⊗ g
	const ket01 = g ⊗ e
	const ket10 = e ⊗ g
	const ket11 = e ⊗ e

	const Φp = normalize(g ⊗ g + e ⊗ e)
	const Φm = normalize(g ⊗ g - e ⊗ e)
	const Ψp = normalize(g ⊗ e + e ⊗ g)
	const Ψm = normalize(g ⊗ e - e ⊗ g)


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
else
	println("Two-qubit operator I2 already defined; not re-including two-qubit operators.")
end
