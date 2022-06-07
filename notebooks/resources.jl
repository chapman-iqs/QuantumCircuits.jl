### A Pluto.jl notebook ###
# v0.17.4

using Markdown
using InteractiveUtils

# ╔═╡ 75236fca-cfac-11eb-2ffe-c394a1c504cf
begin
	path = pwd()
	using PlutoUI
end

# ╔═╡ 59173935-7989-4ca3-a0f7-5fe26cde1127
Markdown.parse("[Table of Contents 📔](./open?path=$(string(path, "/notebooks/table-of-contents.jl")))")

# ╔═╡ 30061cc8-c84f-43b1-830d-a6dc2d9fc540
md"""
# Resources
"""

# ╔═╡ 6754eae6-ea22-4316-8682-c47356f05e5b
begin
	Slichter_et_al = html"<a href='https://iopscience.iop.org/article/10.1088/1367-2630/18/5/053031/pdf'>Slichter et al., 2016 📘</a>"

	Szigeti_et_al = html"<a href='https://arxiv.org/abs/1211.2490'>Szigeti et al., 2013 📘</a>"

	Korotkov_2011 = html"<a href='https://arxiv.org/pdf/1111.4016.pdf'>Korotkov, 2011 📘</a>"

	Jacobs_Steck_2006 = html"<a href='https://arxiv.org/abs/quant-ph/0611067'>Jacobs and Steck, 2006 📘</a>"

	Patti_et_al = html"<a href='https://arxiv.org/abs/1705.03878/'>Patti et al., 2017 📘</a>"

	Itano_et_al = html"<a href='https://journals.aps.org/pra/abstract/10.1103/PhysRevA.41.2295'> Itano et al., 1990 📘</a>"

	Krantz_et_al = html"<a href='https://arxiv.org/abs/1904.06560'> Krantz et al., 2019 📘</a>"

	Campagne_Ibarcq_et_al = html"<a href='https://hal.archives-ouvertes.fr/tel-01248789/file/thesisCampagneibarcq.pdf'> Campagne-Ibarcq et al., 2019 📘</a>"

	Gambetta_et_al = html"<a href='https://arxiv.org/abs/0709.4264'> Gambetta et al., 2019 📘</a>"

	Budini_2014 = html"<a href='https://journals.aps.org/pre/abstract/10.1103/PhysRevE.89.012147'> Budini, 2014 📘</a>"

	Dressel_notes_2019 = "Dressel, 2019"
	

	md" ♦️ **URLs**"
end

# ╔═╡ 31319a7d-dcfd-48ce-a7fd-98ac69145879
md"""
 $Slichter_et_al D. H. Slichter, C. Müller, R. Vijay, S. J. Weber, A. Blais, and I. Siddiqi, Quantum Zeno Effect in the Strong Measurement Regime of Circuit Quantum Electrodynamics, New J. Phys. 18, 053031 (2016).
"""

# ╔═╡ 9db3ea0a-a39d-4c72-b2db-ed080e51f9b2
md" $Korotkov_2011 A. N. Korotkov, Quantum Bayesian Approach to Circuit QED Measurement, ArXiv:1111.4016 \[Cond-Mat, Physics:Quant-Ph\] (2011)."

# ╔═╡ 092fe67d-2623-43fb-abea-2fad5fa93dc7
md"""
 $Jacobs_Steck_2006 K. Jacobs and D. A. Steck, A Straightforward Introduction to Continuous Quantum Measurement, Contemporary Physics 47, 279 (2006).

"""

# ╔═╡ 75bd289c-f727-4cf1-8980-d8785ba7cade
md"""
$Patti_et_al T. L. Patti, A. Chantasri, L. P. García-Pintos, A. N. Jordan, and J. Dressel, Linear Feedback Stabilization of a Dispersively Monitored Qubit, Phys. Rev. A 96, 022311 (2017).

"""

# ╔═╡ dd6ecb57-7de7-40cd-b015-fde3b504376e
md"""
$Itano_et_al Wayne M. Itano, D. J. Heinzen, J. J. Bollinger, and D. J. Wineland, Quantum Zeno effect, Phys. Rev. A 41, 2295 (1990).
"""

# ╔═╡ 76f443c4-f7d9-4b23-ba17-c061d51ac7f8
md"""
$Krantz_et_al
P. Krantz, M. Kjaergaard, F. Yan, T. P. Orlando, S. Gustavsson, and W. D. Oliver, A Quantum Engineer’s Guide to Superconducting Qubits, Applied Physics Reviews 6, 021318 (2019).
"""

# ╔═╡ fa27466c-4bc6-47ca-bb50-07ad365327f4
md"""
$Campagne_Ibarcq_et_al
P. Campagne-Ibarcq, Measurement Back Action and Feedback in Superconducting Circuits, 223 (n.d.).
"""

# ╔═╡ 643b3d2c-ee3e-43b0-a022-9d7f5275ee32
md"""
$Gambetta_et_al
Gambetta, A. Blais, M. Boissonneault, A. A. Houck, D. I. Schuster, and S. M. Girvin, Quantum Trajectory Approach to Circuit QED: Quantum Jumps and the Zeno Effect, Phys. Rev. A 77, 012112 (2008).
"""

# ╔═╡ c7fc0ba1-78b9-4515-b354-30c07b0735cb
md"""
$Szigeti_et_al
S. S. Szigeti, S. J. Adlong, M. R. Hush, A. R. R. Carvalho, and J. J. Hope, Robustness of System-Filter Separation for the Feedback Control of a Quantum Harmonic Oscillator Undergoing Continuous Position Measurement, Phys. Rev. A 87, 013626 (2013).

"""



# ╔═╡ fc25bb7f-b056-41e7-ad8c-af8ced0e1a45
md"""
$Budini_2014 A. A. Budini, Post-Markovian Quantum Master Equations from Classical Environment Fluctuations, Phys. Rev. E 89, 012147 (2014).
"""

# ╔═╡ ca4aa0b9-a0f9-4763-ae00-20a1c1e72ceb
md"""
$Dressel_notes_2019 Dressel, J. Phase-sensitive Qubit Monitoring with Resonator Transients. Unpublished notes, 2019.
"""

# ╔═╡ d7982457-b1eb-40f7-a7ab-a9053baf4e08
md"""
# Utilities
"""

# ╔═╡ 3c3d59df-1eab-4b1b-874c-bd5457da8d6f
function notebook_path(folder, index)
	notebooks = readdir(folder)
	string(path, "/", folder, "/", notebooks[index])
end

# ╔═╡ a5cbd8d7-c099-44c9-bd48-298246062c97
path

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.34"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.0-DEV.1208"
manifest_format = "2.0"
project_hash = "e766545f1b4ef968b5991d6a702e32de0b70adeb"

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
git-tree-sha1 = "13468f237353112a01b2d6b32f3d0f80219944aa"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8979e9802b4ac3d58c503a20f2824ad67f9074dd"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.34"

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
# ╟─59173935-7989-4ca3-a0f7-5fe26cde1127
# ╟─30061cc8-c84f-43b1-830d-a6dc2d9fc540
# ╟─31319a7d-dcfd-48ce-a7fd-98ac69145879
# ╟─9db3ea0a-a39d-4c72-b2db-ed080e51f9b2
# ╟─092fe67d-2623-43fb-abea-2fad5fa93dc7
# ╟─75bd289c-f727-4cf1-8980-d8785ba7cade
# ╠═dd6ecb57-7de7-40cd-b015-fde3b504376e
# ╟─76f443c4-f7d9-4b23-ba17-c061d51ac7f8
# ╟─fa27466c-4bc6-47ca-bb50-07ad365327f4
# ╟─643b3d2c-ee3e-43b0-a022-9d7f5275ee32
# ╟─c7fc0ba1-78b9-4515-b354-30c07b0735cb
# ╟─fc25bb7f-b056-41e7-ad8c-af8ced0e1a45
# ╟─ca4aa0b9-a0f9-4763-ae00-20a1c1e72ceb
# ╟─6754eae6-ea22-4316-8682-c47356f05e5b
# ╟─d7982457-b1eb-40f7-a7ab-a9053baf4e08
# ╟─3c3d59df-1eab-4b1b-874c-bd5457da8d6f
# ╠═75236fca-cfac-11eb-2ffe-c394a1c504cf
# ╟─a5cbd8d7-c099-44c9-bd48-298246062c97
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
