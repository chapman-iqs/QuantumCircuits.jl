### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

# ╔═╡ 75236fca-cfac-11eb-2ffe-c394a1c504cf
begin
	directory_name = "QC-notebooks"
	path = let 
			arr = split(pwd(), "/")
			index = findfirst(s -> s == directory_name, arr)
			s = join(map(a -> string(a, "/"), arr[1:index]))
			string(s, "notebooks")
			
	end
	cd(path)
	
	using PlutoUI
end

# ╔═╡ ca0ba638-d91e-41c5-bf99-92d31fe016c0
TableOfContents()

# ╔═╡ e6781a16-226c-4d78-91a2-c50b40c7720e
md"""
# Notebooks
"""

# ╔═╡ 7d387428-bc2f-4486-865f-814ce68702c7
md" ## 1. Quantum mechanics intro "

# ╔═╡ e5300f3d-b000-40f6-ae00-9973a678494c
md"""
### Basics
"""

# ╔═╡ 8be32eea-7335-493c-963e-c5f68f2e8096
md"""
### Measurement
"""

# ╔═╡ 6846fd9e-29cb-44cb-9e16-b5421000a52e
md"""
### Experiment
Introduces experimental concepts in superconducting qubit dynamics and measurement.
"""

# ╔═╡ 543db20b-af94-46ba-90a9-274769f2d0a9
md" ## 2. Usage and checks"

# ╔═╡ 4db78c6f-92ec-4d7d-9dda-db61789dfc1b
md" ## 3. Demos"

# ╔═╡ 14b2f760-bb74-43c9-bfb5-14821d6222d7
md" ## 4. Excitation feedback"

# ╔═╡ 1bf03b51-5349-46b3-b904-9782d57a9f9f
md" ## 5. Linear feedback"

# ╔═╡ 3247d5f7-60f6-4201-9d7b-5554a69fd6f8
md"""
## 6. Other
"""

# ╔═╡ 099b4fbb-77a2-4827-8466-5f68d9638641
md"""
## In development
"""

# ╔═╡ d7982457-b1eb-40f7-a7ab-a9053baf4e08
md"""
# Utilities
"""

# ╔═╡ 10707e74-03ef-45cc-88bc-8a81fbcca162
begin
	quantum_states📔 = "[Quantum States 📔](./open?path=$(string(path, "/1-quantum-mechanics-intro/1-quantum-states.jl")))"
	
	bloch_sphere📔 = "[Bloch Sphere 📔](./open?path=$(string(path, "/1-quantum-mechanics-intro/2-bloch-sphere.jl")))"

	qubit_ensembles📔 = "[Qubit Ensembles 📔](./open?path=$(string(path, "/1-quantum-mechanics-intro/3-qubit-ensembles.jl")))"

	measurement_backaction📔 = "[Measurement Backaction 📔](./open?path=$(string(pwd(), "/1-quantum-mechanics-intro/4-measurement-backaction.jl")))"

	quantum_zeno📔 = "[Quantum Zeno Effect 📔](./open?path=$(string(pwd(), "/1-quantum-mechanics-intro/5-quantum-zeno.jl")))"

	coherent_states📔 = "[Coherent State Dynamics 📔](./open?path=$(string(pwd(), "/1-quantum-mechanics-intro/measurement/coherent-states.jl")))"

	demodulation📔= "[Demodulation 📔](./open?path=$(string(pwd(), "/1-quantum-mechanics-intro/6-experiment/demodulation.jl")))"

	resonator_states📔= "[Resonator States 📔](./open?path=$(string(pwd(), "/1-quantum-mechanics-intro/6-experiment/resonator-states.jl")))"

	checks📔 = "[Checks 📔](./open?path=$(string(path, "/2-usage-and-checks/1-checks.jl")))"

	bayesian_rouchon📔 = "[Bayesian and Rouchon 📔](./open?path=$(string(path, "/2-usage-and-checks/2-bayesian-rouchon.jl")))"

	rabi_dragging📔 = "[Rabi Dragging 📔](./open?path=$(string(path, "/3-demos/1-rabi-dragging.jl")))"

	excitation_feedback📔 = "[Excitation Feedback 📔](./open?path=$(string(path, "/4-two-qubit-feedback/1-excitation-feedback.jl")))"

	two_qubit_decoherence📔 = "[Two-qubit Decoherence 📔](./open?path=$(string(path, "/4-two-qubit-feedback/2-two-qubit-decoherence.jl")))"

	system_filter📔 = "[System-filter 📔](./open?path=$(string(path, "/4-two-qubit-feedback/3-system-filter.jl")))"

	excitation_feedback_3q📔 = "[Three-qubit Excitation Feedback 📔](./open?path=$(string(path, "/4-two-qubit-feedback/4-excitation-feedback-3q.jl")))"

	linear_feedback📔 = "[Linear Feedback 📔](./open?path=$(string(path, "/5-linear-feedback/1-linear-feedback.jl")))"

	linear_feedback_checks📔 = "[Linear Feedback Checks 📔](./open?path=$(string(path, "/5-linear-feedback/2-linear-feedback-checks.jl")))"

	table_of_contents📔 = "[Table of Contents 📔](./open?path=$(string(path, "/table-of-contents.jl")))"

	resources📔 = "[Resources 📔](./open?path=$(string(path, "/resources.jl")))"

	utilities📔 = "[Utilities 📔](./open?path=$(string(path, "/notebook-utilities.jl")))"

	quadrature_measurement📔 = "[Quadrature measurement 📔](./open?path=$(string(path, "/development/quadrature-measurement.jl")))"

	


end

# ╔═╡ f03c972c-8917-4bb3-8436-604a93889bc1
mdp(strings...) = let
					s = ""
					for str in strings
						s = string(s, str)
					end
					Markdown.parse(s)
end

# ╔═╡ c830573d-c042-480c-90f2-671598580cbd
mdp("", quantum_states📔, " Introduces the idea of the quantum state.")

# ╔═╡ d4b07505-53de-4c53-8ccb-82ae02f86dcd
mdp("", bloch_sphere📔, " Looks at quantum trajectories on the Bloch sphere.")

# ╔═╡ 94b1dcfb-340b-4281-9af8-2328293e7ef9
mdp("", qubit_ensembles📔, " Looks at qubit ensemble behavior under weak measurement.")

# ╔═╡ 3ddb3f0f-47a5-4932-afca-39f7b0a71d6f
mdp("", measurement_backaction📔, " Compares informational and phase backaction for a qubit.")

# ╔═╡ 8627ef8f-99fa-4b9a-b99d-d0e1912f73ec
mdp("", quantum_zeno📔, " Introduces the Quantum Zeno effect to give intuition for readout.")

# ╔═╡ 13b8f2bb-f72e-4258-81a1-f4a7c35ff1e0
mdp("", coherent_states📔, " Discusses the coherent state approximation for the readout resonator and resulting reduced qubit evolution.")

# ╔═╡ 3f99dd1a-20a5-4715-8bda-411ba1790440
mdp("", demodulation📔, " Introduces the concept of demodulation in qubit readout.")

# ╔═╡ 260bd30f-6c30-4960-89e2-e9474ef5a537
mdp("", checks📔, " Outlines various checks on bayesian and rouchon solvers. Currently only checks that ensemble average converges to η = 0 solution.")

# ╔═╡ 82e9aca0-8ba8-4638-a7bf-5d112aaceef0
mdp("", bayesian_rouchon📔, " Compares features of `bayesian` and `rouchon` solvers, including checks on record normalization, consistency across different choices of `dt`, and reproducing results between the two solvers.")

# ╔═╡ 731a165d-3470-43ec-baba-06023ee9dd73
mdp("", rabi_dragging📔, " Explores adiabatically changing the Rabi axis during evolution to drag the quantum state.")

# ╔═╡ 29d57063-367a-408c-b784-58a08a31eecf
mdp("", excitation_feedback📔, " Explores generating two-qubit entangled states (Bell states) using different weak measurement feedback protocols. Focuses on measurement of the total excitation number.")

# ╔═╡ 9c2f2aa9-9310-4fb6-a77b-bab7d6787916
mdp("", two_qubit_decoherence📔, " Checks what happens to entangled qubit states under different decoherence mechanisms. Currently only implemented for random σz rotations.")

# ╔═╡ 1684ca47-791d-4975-94d1-9675a668a440
mdp("", system_filter📔, " Tests and explains the system-filter functionality of QuantumCircuits.jl, using excitation feedback as an example.")

# ╔═╡ bf24618d-5f34-41c9-b2cc-29da7d0599a7
mdp("", excitation_feedback_3q📔, " Implements excitation feedback on three-qubit entangled states.")

# ╔═╡ 46e30a5b-16b6-446a-914d-c10110970aaa
mdp("", linear_feedback📔, "Introduces linear feedback stabilization.")

# ╔═╡ f3982cad-3f00-446e-8869-26cfd18bf557
mdp("", linear_feedback_checks📔, " Compares bloch equation solutions to the linear feedback problem, and compares to matrix exponentation and `bayesian` solutions.")

# ╔═╡ ab1baadd-963e-40ef-b9b3-5b81635d36dc
mdp("", resources📔, " Includes citations and links for articles cited in notebooks.")

# ╔═╡ 1777e7d8-7338-4a27-8b9e-b354defd760e
mdp("", utilities📔, " Includes plotting utilities commonly used in notebooks.")

# ╔═╡ ad4bb00b-a172-4e50-9ef4-0cd3cca26188
mdp("", quadrature_measurement📔, "Resonator fock space simulations for qubit-resonator dynamics. Has an unresolved bug for dual quadrature measurement.")

# ╔═╡ a5cbd8d7-c099-44c9-bd48-298246062c97
path

# ╔═╡ 3c3d59df-1eab-4b1b-874c-bd5457da8d6f
function notebook_path(folder, index)
	notebooks = readdir(folder)
	string(path, "/", folder, "/", notebooks[index])
end

# ╔═╡ 0d8972f3-6a8b-402d-bb5a-6ca084881881
# begin
# 	folders = readdir()
# 	deleteat!(folders, findall(x->x==".DS_Store",folders))
# 	deleteat!(folders, findall(x->x=="processing.jl",folders))
# 	deleteat!(folders, findall(x->x=="100-refactoring",folders))
# 	deleteat!(folders, findall(x->x=="table-of-contents.jl",folders))
# 	folder_names = map(folder -> replace(folder, "-" => " "), folders)
# end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.38"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.0-DEV.1208"
manifest_format = "2.0"
project_hash = "9189d135f87bb42ffdf8c0fdd7f2d0c021e78ebf"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.0+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.73.0+4"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.9.1+2"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.24.0+2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2020.7.22"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.17+2"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "621f4f3b4977325b9128d5fae7a8b4829a0c2222"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.4"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "670e559e5c8e191ded66fa9ea89c97f10376bb4c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.38"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "3.1.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.41.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "16.2.1+1"
"""

# ╔═╡ Cell order:
# ╟─ca0ba638-d91e-41c5-bf99-92d31fe016c0
# ╟─e6781a16-226c-4d78-91a2-c50b40c7720e
# ╟─7d387428-bc2f-4486-865f-814ce68702c7
# ╟─e5300f3d-b000-40f6-ae00-9973a678494c
# ╟─c830573d-c042-480c-90f2-671598580cbd
# ╟─d4b07505-53de-4c53-8ccb-82ae02f86dcd
# ╟─94b1dcfb-340b-4281-9af8-2328293e7ef9
# ╟─8be32eea-7335-493c-963e-c5f68f2e8096
# ╟─3ddb3f0f-47a5-4932-afca-39f7b0a71d6f
# ╟─8627ef8f-99fa-4b9a-b99d-d0e1912f73ec
# ╟─13b8f2bb-f72e-4258-81a1-f4a7c35ff1e0
# ╟─6846fd9e-29cb-44cb-9e16-b5421000a52e
# ╟─3f99dd1a-20a5-4715-8bda-411ba1790440
# ╟─543db20b-af94-46ba-90a9-274769f2d0a9
# ╟─260bd30f-6c30-4960-89e2-e9474ef5a537
# ╟─82e9aca0-8ba8-4638-a7bf-5d112aaceef0
# ╟─4db78c6f-92ec-4d7d-9dda-db61789dfc1b
# ╟─731a165d-3470-43ec-baba-06023ee9dd73
# ╟─14b2f760-bb74-43c9-bfb5-14821d6222d7
# ╟─29d57063-367a-408c-b784-58a08a31eecf
# ╟─9c2f2aa9-9310-4fb6-a77b-bab7d6787916
# ╟─1684ca47-791d-4975-94d1-9675a668a440
# ╟─bf24618d-5f34-41c9-b2cc-29da7d0599a7
# ╟─1bf03b51-5349-46b3-b904-9782d57a9f9f
# ╟─46e30a5b-16b6-446a-914d-c10110970aaa
# ╟─f3982cad-3f00-446e-8869-26cfd18bf557
# ╟─3247d5f7-60f6-4201-9d7b-5554a69fd6f8
# ╟─ab1baadd-963e-40ef-b9b3-5b81635d36dc
# ╟─1777e7d8-7338-4a27-8b9e-b354defd760e
# ╟─099b4fbb-77a2-4827-8466-5f68d9638641
# ╟─ad4bb00b-a172-4e50-9ef4-0cd3cca26188
# ╟─d7982457-b1eb-40f7-a7ab-a9053baf4e08
# ╠═10707e74-03ef-45cc-88bc-8a81fbcca162
# ╠═f03c972c-8917-4bb3-8436-604a93889bc1
# ╠═75236fca-cfac-11eb-2ffe-c394a1c504cf
# ╠═a5cbd8d7-c099-44c9-bd48-298246062c97
# ╠═3c3d59df-1eab-4b1b-874c-bd5457da8d6f
# ╠═0d8972f3-6a8b-402d-bb5a-6ca084881881
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
