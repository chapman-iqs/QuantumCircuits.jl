"Single-qubit operators ------------------------------------------------"

# Basis
q = SpinBasis(1//2)
Iq = identityoperator(q)

# qubit operators, using convention that |-z> is ground state (aligned with QuantumOpticsBase convention)
X = sigmax(q)
Y = sigmay(q)
Z = sigmaz(q)
P = sigmap(q)
M = sigmam(q)

g = spindown(q)
e = spinup(q)


"Resonator operators ------------------------------------------------------------------------------------"

# Basis
Nfock = 30 # fock space dimension cutoff
f = FockBasis(Nfock)
If = identityoperator(f)

# Qubit-resonator operators
Id = Iq ⊗ If
a = Iq ⊗ destroy(f)

σx = X ⊗ If
σy = Y ⊗ If
σz = Z ⊗ If
σp = P ⊗ If
σm = M ⊗ If

# projectors
zp = dm(e) ⊗ If
zm = dm(g) ⊗ If

αp = a * zp
αm = a * zm


"bases and labels ------------------------------------------------------------------------------------"
qbasis = [σx, σy, σz]
qlabels = ["x", "y", "z"]
