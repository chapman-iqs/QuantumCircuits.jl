### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 45b3db58-09c7-11ec-26e4-6fdfa55e989f
begin
	cd("/Users/sachagreenfield/Desktop/GitHub/QuantumCircuits.jl")
	import Pkg
	Pkg.activate(".")
	using Random
	using Statistics
	using Distributions
	using PyPlot
	using QuantumCircuits
end

# ╔═╡ 395cf98f-af1f-4696-84eb-5db637215f14
md"""
Using convention $|target\rangle \otimes |control\rangle$:
"""

# ╔═╡ e499d688-1c53-426e-acd6-7c903d9843df
begin
	# Basis
	q = SpinBasis(1//2)

	# single qubit operators
	σx = sigmax(q)
	σy = sigmay(q)
	σz = sigmaz(q)
	σp = sigmap(q)
	σm = sigmam(q)
	id = identityoperator(q)
	zp = 0.5(id + σz)
	zm = 0.5(id - σz)
	p01 = 0*id
	p01.data[1,2] = 1
	p10 = 0*id
	p10.data[2,1] = 1
	
	Rx(θ) = cos(θ/2)*id - im*sin(θ/2)*σx
	Ry(θ) = cos(θ/2)*id - im*sin(θ/2)*σy
	Rz(θ) = cos(θ/2)*id - im*sin(θ/2)*σz
	
	# two-qubit operators
	CNOT = (σx ⊗ zm) + (id ⊗ zp)
	CNOT2 = (zm ⊗ σx) + (zp ⊗ id)
	CZ = (σz ⊗ zm) + (id ⊗ zp)	
	
end

# ╔═╡ 0b67337f-3bda-444f-badb-dca8cb5a8540
begin
	iSWAP = 0*(id ⊗ id)
	iSWAP.data[1,1] = 1
	iSWAP.data[4,4] = 1
	iSWAP.data[2,3] = 1*im
	iSWAP.data[3,2] = 1*im
end

# ╔═╡ d0ac9630-8e3c-42eb-aebc-88d9671479c4
begin
	H = 0*id
	H.data[1,1]=1/sqrt(2)
	H.data[1,2] = 1/sqrt(2)
	H.data[2,1] =1/sqrt(2)
	H.data[2,2] = -1/sqrt(2)
end

# ╔═╡ 2e6fed05-ec81-4cfc-be38-15431249d0fa
(id ⊗ H) * CNOT * CNOT2 * (H ⊗ id)

# ╔═╡ bdaa2be7-04d2-4b6f-9b40-777b4b332fe0
(id ⊗ H) * CNOT

# ╔═╡ 3ac09d40-c499-4b49-a6fc-64989dd7717f
H

# ╔═╡ 7b39081e-3581-4100-ad7c-7c0fb13265ac
H

# ╔═╡ 0fc46c39-bcf0-4ae1-868e-a2ee2c709fcb
roundop(op) = round.(op.data, digits=5)

# ╔═╡ 9b207662-25ab-4f1c-bbf5-8321c9271b1a
im*roundop(Rx(π) * Ry(π/2))

# ╔═╡ 8d8fdaa2-8ca9-4bff-8f1a-7cada9c18bbd
roundop((id ⊗ Rz(π/2) )* iSWAP * (Rx(π/2) ⊗ id) * iSWAP * (Rz(-π/2) ⊗ Rz(π/2)) * (id ⊗ Rx(π/2)))

# ╔═╡ 21f80c93-cca4-4383-9c47-98d3eef1bf73
roundop((id ⊗ H) * CNOT *  CNOT2 * (H ⊗ id))

# ╔═╡ 1111261b-37d0-448e-bfe3-521b31aadd26
cz=roundop((Ry(-π/2) ⊗ id) * CNOT * (Ry(π/2) ⊗ id))

# ╔═╡ 93c8bdf3-ed7f-4a40-afdf-ddd48777e28d
roundop((id ⊗ Ry(-π/2) ) * CNOT * (id⊗ Ry(π/2) ))

# ╔═╡ 74d278dc-f794-433c-a4ab-0ed2e494871b
begin
	Id = id ⊗ id
	
end

# ╔═╡ Cell order:
# ╠═45b3db58-09c7-11ec-26e4-6fdfa55e989f
# ╟─395cf98f-af1f-4696-84eb-5db637215f14
# ╠═e499d688-1c53-426e-acd6-7c903d9843df
# ╠═2e6fed05-ec81-4cfc-be38-15431249d0fa
# ╠═bdaa2be7-04d2-4b6f-9b40-777b4b332fe0
# ╠═9b207662-25ab-4f1c-bbf5-8321c9271b1a
# ╠═3ac09d40-c499-4b49-a6fc-64989dd7717f
# ╠═0b67337f-3bda-444f-badb-dca8cb5a8540
# ╠═d0ac9630-8e3c-42eb-aebc-88d9671479c4
# ╠═7b39081e-3581-4100-ad7c-7c0fb13265ac
# ╠═8d8fdaa2-8ca9-4bff-8f1a-7cada9c18bbd
# ╠═21f80c93-cca4-4383-9c47-98d3eef1bf73
# ╠═1111261b-37d0-448e-bfe3-521b31aadd26
# ╠═0fc46c39-bcf0-4ae1-868e-a2ee2c709fcb
# ╠═93c8bdf3-ed7f-4a40-afdf-ddd48777e28d
# ╠═74d278dc-f794-433c-a4ab-0ed2e494871b
