### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 80ffb52c-a0dc-43d9-8282-e8acb51df4e0
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

# ╔═╡ 595fc42b-943e-4a29-841a-cd36c90a2b55
md"""
##### Quantum basis
"""

# ╔═╡ 483e648c-2ac3-46f7-b49b-a3109deec27d
begin
	# Basis
	q = SpinBasis(1//2)

	# Operators, using convention that |-z> is ground state
	σx = sigmax(q)
	σy = sigmay(q)
	σz = sigmaz(q)
	σp = sigmap(q)
	σm = sigmam(q)
	id = identityoperator(q)
end

# ╔═╡ 52a501dd-c011-41cb-8b47-ce5ff894761c
md"""
#### Parameters and operators
"""

# ╔═╡ 5d30c7c1-4c53-40d2-8e73-b58db5697cfa
md"""
# Bayesian vs. Rouchon integration methods

In this interactive notebook, we'll compare quantitative and qualitative features of the two primary integration methods -- `bayesian` and `rouchon` -- implemented in `QuantumCircuits.jl`. What follows is a brief summary of key commonalities and differences before we dive into the details.
"""

# ╔═╡ d4e7b79b-eb8c-48b3-9a94-650481fa8cb7
md"""
##### Commonalities

The `bayesian` and `rouchon` methods have important commonalities. 

Both:
* Model quantum trajectories of weakly and continuously measured open quantum systems.

* Take `H`, `[J]`, `[C]` (Hamiltonian, Lindblad, collapse) operators as inputs.

* If simulating quantum trajectories, can generate a noisy record of a monitored variable $r = \langle x \rangle + dW/dt$, with Wiener process $dW \sim \mathcal{N}(0, \sqrt{\tau_m})$ (**CHECK THIS**) if running a simulation. If reconstructing experimental trajectories, can process an imported experimental record.

"""

# ╔═╡ bef8c648-af4c-46c1-87b4-29f744f9aaf3
md"""
##### API

The two methods have the same API:
```
rouchon(T, ρ0, H, J, C; <keyword arguments>)

Arguments:

T :: tuple (ti,tf)
ρ0 :: initial density matrix
H :: time-dependent Hamiltonian
J :: deterministic collapse operators
C :: stochastic collapse operators

Keyword Arguments:

dt :: time step; default dt=1e-4
dy :: record; default dy=[], i.e. simulation generates record time series.
            record should be input in the shape [dy_1,...,dy_Nc] given
            length(C)=Nc collapse operators. Records should have shape
            dy_c[t] indexing the mth trajectory at time T[n]
fn : ρ → Any :: eval function (e.g. to return expectation values instead of density matrices)

Returns: (ts, ρs, dy)

ts :: list of simulation times
ρs :: fn(ρ) at each simulation time
dy :: input OR simulated record, depending on value of keyword argument dy

```
"""

# ╔═╡ 2059a835-925d-430d-aa0c-283a03f6ba2d
md"""
In most (all?) implementations of quantum circuit dynamics, inputs will take the following form:

* `H = H0` or `H = H(t)` : time-in/dependent Hamiltonian;

* `J = [J0, ..., Jn, √(1-η_0)C0, ... , √(1-η_m)Cm]` : list of $n$ dissipation operators $\{J_0, ..., J_n\}$ and $m$ fluctuation-dissipation operators $\{\sqrt{1-\eta_0} C_0, ..., \sqrt{1-\eta_m} C_m \}$ with $n, m \geq 0$;

* `C = [√η C0, ..., √Cm]` : list of $m$ fluctuation operators $\{\sqrt{\eta_0} C_0, ..., \sqrt{\eta_m} C_m \}$.
"""

# ╔═╡ 476f0b08-7d21-4010-b114-130e6dfbbae0
md"""
#### Example
###### Parameters and operators
"""

# ╔═╡ 639931ed-0520-499c-ba1c-90b04e854ebd
begin
	ΩR  = 2π # Rabi frequency
	τm = 3.0 # Measurement collapse timescale
	Γm = 1/(2τm) # Measurement dephasing rate
	η = 0.3
	
	H = ΩR*σx/2
	J = [√((1-η)*Γm)*σz] 
	C = [√(Γm*η)*σz]
	
	T = (0,4τm) # simulation duration
	ρ0 = dm(spinup(q))
	dt = 1e-3 # integration time-step
end

# ╔═╡ 68627f37-30c6-47c2-85e3-a031cf4a2e05
md" ###### Simulation: Bayesian"

# ╔═╡ 0c3770b6-02b2-46af-a875-cea82999f88f
begin
	Random.seed!(10)
	solb = bayesian(T, ρ0, H, J, C; dt=dt, heterodyne=false)
end

# ╔═╡ eb250fd3-7c43-4f75-ab9e-6f7ae0035ad3
begin
	(tt, ρs, dys) = solb
	Rs = collect(dys[1]) # record output from bayesian
	dy = Rs*dt # rescale record for input into rouchon
end

# ╔═╡ 75f1cc94-5222-42d5-8c0e-ceedb3c53d3b
md" ###### Simulation: Rouchon"

# ╔═╡ bdf4c000-2845-4e2f-baf6-0d1161e477e3
begin
	Random.seed!(1)
	solr = rouchon(T, ρ0, H, J, C; dt=dt, dy=[dy])
end

# ╔═╡ b6c73eea-d17c-49d4-a091-be10c8134718
md" ### Bayesian filtering "

# ╔═╡ 4fc42181-a1d4-46d9-be76-7277b9a2d6a5
md" ##### Update equations"

# ╔═╡ 343161e9-e975-45b5-a444-b9f9f0b30cfb
md"""
Bayesian filtering describes the evolution of a weakly measured quantum system in fundamentally discrete chunks of time, $\Delta t$. In general, the system is described by a density matrix $\rho$. It undergoes **unitary evolution** according to a  Hamiltonian $H$, such that

$\hat \rho' = \hat U \rho {\hat U}^\dagger,$

where 

$\hat U = \exp(- i \Delta t \hat H).$

(Note that $\hbar$ has been absorbed into $\hat H$ here.) 


"""

# ╔═╡ d5af31e4-ba93-4f3f-9c55-4f2ce600d797
md"""
The **weak measurement** of an observable $\hat C$ is represented by a Kraus operator $\hat M_{r_C}$, such that

$\hat M_{r_C} = \exp\Big[\frac{\Delta t}{2 \tau_m} {r_C}^* \hat C  - \frac14 (\hat C + \hat C^\dagger)^2\Big].$

"""

# ╔═╡ 97e2560f-3dd9-45a0-85f5-aa9e12800292
md"""

Here, $r_C$ is a complex noisy record associated with measuring $\hat C$, while $\tau_m$ is the timescale of measurement collapse. From each $\hat M_{r_C}$ we can construct a POVM (positive operator valued measure)

$E_{r_C} = \hat M_{r_C}^\dagger \hat M_{r_C}$

satisfying

$\sum_{r_C} E_{r_C} = \mathbb{1}.$

(There are, of course, infinitely many possible records $r_C$ corresponding to different noise realizations. So this might not be the right way to put it.) 

The weak measurement updates the state according to

$\hat \rho '' = \hat M_{r_C}^\dagger \rho' \hat M_{r_C}.$

"""


# ╔═╡ d6d6e78b-15cb-4c6a-af30-8c42b35dc1b9
md"""
Finally, the system **decoheres** according to decay operators

$M_{\text{decay, } \hat J} = \sqrt{\kappa_J \Delta t} \hat J$

and

$M_{\text{null, } \hat J} = \sqrt{\mathbb{1} - \kappa_J \Delta t \hat J^\dagger \hat J}$

satisfying

$M_{\text{decay, } \hat J}^\dagger M_{\text{decay, } \hat J} + M_{\text{null, } \hat J}^\dagger M_{\text{null, } \hat J} = \mathbb{1}.$

The $\hat J$ are Lindblad operators. The state updates as

$\rho''' = M_{\text{decay, } \hat J} \rho'' M_{\text{decay, } \hat J}^\dagger + M_{\text{null, } \hat J} \rho'' M_{\text{null, } \hat J}^\dagger.$

"""

# ╔═╡ 6315198d-4962-4686-b436-80a3525e0dca
md"""
The final state after increment $\Delta t$ is obtained by renormalizing by the trace:

$\rho(t + \Delta t) = \frac{\rho'''}{\text{Tr}( \rho''')}.$

"""

# ╔═╡ 40a4c411-745d-4ac9-8cd8-ab4a66979fa9
md"""
##### Assumptions

Bayesian filtering assumes small $\Delta t$ in the interleaving of unitary, measurement, and decay evolution. However, it does not assume small $\Delta t$ in pure measurement evolution. Thus, in the absence of unitary evolution and dissipation (i.e. pure measurement evolution), the Bayesian method yields an exact result that is robust to choice of $\Delta t$.
"""

# ╔═╡ a306838f-f883-4fe9-94e6-8d5771cd9793
md" ### Rouchon approximation "

# ╔═╡ f66192ba-ff34-499e-ae8e-e8d08b91b9a1
md""" 
The Rouchon method was presented in [1] as a numerical integration method that provides an efficient alternative to the stochastic master equation (SME). Rouchon's method expands to second order in the Wiener increment $\Delta W_r$. As a result, the method preserves positivity of the conditioned quantum state, unlike directly integrating the SME.

"""

# ╔═╡ a541146f-9143-4a91-bb3a-be17e9bb6561
md"""
$M(t) = \mathbb{1} - \Big (i H + \frac12 \sum_j J_j^\dagger J_j ) dt + \sum_c \sqrt{\eta_c} C_c dy_c(t) +  \sum_{c,d} \frac{\sqrt{\eta_c \eta_d}}2 C_c C_d \big(dy_c(t) dy_d(t) - \delta_{c,d} dt \big)$
"""

# ╔═╡ 88179c16-f111-462d-a9d6-72e3764715eb
md"""
$dy_c(t) = \sqrt{\eta_c} \text{Tr} \big ( C_c \rho(t) + \rho(t) C_c^\dagger \big ) dt + dW_c(t)$
"""

# ╔═╡ 6e92e073-d122-4a84-8c4c-b19858c173a9
md"""
$D(t) = \Big (\sum_j J_j \rho(t) J_j^\dagger  + \sum_c C_c \rho(t) C_c^\dagger \Big ) dt$
"""

# ╔═╡ 501a6621-8054-4fa3-b9a7-f1a921e14a56
md"""
$\rho(t + dt) = \frac{M(t) \rho(t) M(t)^\dagger + D(t)}{\text{Tr} \Big ( M(t) \rho(t) M(t)^\dagger + D(t) \Big )}$
"""

# ╔═╡ cb87a45d-4ff1-4da9-8bbe-ec2c193cb7c6
md"""
The Rouchon method is similar to the Bayesian method in that it uses completely positive maps, and therefore preserves positivity of the density matrix. However, it is still a first-order approximation in $dt$. Whereas the Bayesian method will be exact when measurement is the only source of system evolution, the Rouchon method is not exact under any circumstance.
"""

# ╔═╡ 676e029d-1141-41b2-bc81-5d2e4500c684
md"""
## Timing comparison

Rouchon is faster for given $\Delta t$, but Bayesian is more robust to large $\Delta t$.
"""

# ╔═╡ df840a53-c8b0-4215-97de-ef6f1c720ed7
md" ##### Coarse-graining: bayesian "

# ╔═╡ 7ffe3ce0-2ce5-4012-8679-93b9154b150e
md"""
In order to compare across time-step sizes, we first must generate consistent records corresponding to various time-steps. In weak measurement, $\Delta t$ has a physical interpretation as the integration time for the measurement. Thus, when we change $\Delta t$, we must both resample *and* rescale the measurement record in order to reproduce the same evolution.

"""

# ╔═╡ e41a9229-34f2-4a42-810b-7c95dc02780b
md"""
Let's start by looking at the record generated by Bayesian in our example at the beginning of the notebook:
"""

# ╔═╡ 4f3c3104-c64d-4b77-8932-c37f4aaabaa1
Rs

# ╔═╡ 0b0f2aa0-431d-4029-b7c2-911aa4d2f69a
md"""
1
$(@bind n html"<input type=range min=1. max=20. step=1 value=10>")
 20 
"""

# ╔═╡ 3e0a12d2-b4bf-4424-8b62-d9a714e393fc
md"""
Note that `bayesian` takes in / outputs a record of the form 
	
$R \equiv \frac{r}{\sqrt{\tau_m}} = \frac{z}{\sqrt{\tau_m}} + \frac{dW}{dt},$
	
where $dW \sim \mathcal N(0,\sqrt{dt})$ is a Wiener increment. Because of the $dt$ in the denominator of $dW/dt$, the record is distributed with standard deviation
$\sigma = 1/\sqrt{dt}$. So changing $dt$ changes the variance of our record. In order to get an "equivalent" record as we increase $dt$, we coarse-grain the record.

`QuantumCircuits` has a built-in function `subseries` **(not built-in yet)** that does this automatically. The code below coarse-grains the record at a scale of n = $n, which creates a new record corresponding to measurement bins of width

```math
dt' = 1/\sqrt{n \hspace{1mm}dt}
```
or $dt'$ =  $(dt*n) μs.

"""

# ╔═╡ 69eaed7c-392d-40d8-9f0c-63b4bc4aa18d
md"""
You can try playing with the coarse-graining scale below:

n = $n

"""

# ╔═╡ 012e9206-f306-4776-a546-d6f7dcf00ef2
md"""
We can look at these two records as time-series or as histograms:
"""

# ╔═╡ 275ffd8d-274b-4fab-946e-d9482fa340a6
md"""
Since the record is Gaussian-distributed, time-smoothing the record using `subseries` reduces the variance of the distribution. We can check that the variances reported in the plot above are close to the theoretical values at each $dt$:

fine-grained: $\sigma = 1/\sqrt{dt}$ = $(round(sqrt(1/dt),digits=3))

coarse-grained: $\sigma = 1/\sqrt{n \hspace{1mm} dt}$ = $(round(sqrt(1/(dt*n)),digits=3))
"""

# ╔═╡ e0cb7572-3d67-47e0-9482-70588a693d41
md" ##### Time-step comparison: bayesian"

# ╔═╡ f870609e-f3bc-40e4-90fb-484e7b52db4b
md"""
First we'll get the subseries of the fine-grained simulation:
"""

# ╔═╡ 3f7d824c-8a0e-4547-bb06-81717a01df40
md"""
We also need the subrecord to run a coarse-grained trajectory:
"""

# ╔═╡ 7c434d23-a5ed-41bb-8dc7-8b1440a0680b
md"""
Now we can compare:
"""

# ╔═╡ 82b18884-28d4-4ecb-b21c-d3e57dd64870
md"""
1
$(@bind ns html"<input type=range min=1. max=100. step=1 value=10>")
 100
"""

# ╔═╡ e47ee9c6-c9d6-44b8-8fc3-e6ebb58f28e3
md"""
ns = $ns
"""

# ╔═╡ 82d43f38-9cf6-4fee-bdc8-502234781678
pairs=((L"$y$",round(dt*ns,digits=3)), (L"$y$",dt), (L"$z$",round(dt*ns,digits=3)),(L"$z$",dt))

# ╔═╡ 499291df-d136-4dab-8a6d-0d72cf812435
md" ##### Coarse-graining: rouchon"

# ╔═╡ 3d3d4838-2c53-44b5-a103-5a7656d58611
md"""
Let's do the same for `rouchon`. Note that `rouchon` takes in / outputs a record of the form 
	
```math
dy \equiv r \frac{dt}{\sqrt{\tau_m}} = z \frac{dt}{\sqrt{\tau_m}}+ dW,
```
	
where $dW \sim \mathcal N(0,\sqrt{dt})$ is a Wiener increment. As a result, the record is distributed with standard deviation
$\sigma = \sqrt{dt}$. So multiplying $dt' = n \hspace{1mm} dt$ changes the standard deviation of the record by $\sigma' = \sqrt{n} \sigma$.

In order to get the equivalent record, we cannot simply coarse-grain $dy$, as $dy'$ should actually have a *larger* standard deviation. Instead, we coarse grain the equivalent bayesian record $R$, which is related to $dy$ by $dy = R \hspace{1mm} dt$ to obtain $R'$, then get $dy' = R' \hspace{1mm} dt'$.

"""

# ╔═╡ b9209585-d439-488e-b65f-d49a61d53bd8
md" ##### Timing comparison: rouchon"

# ╔═╡ 7e1f7d9e-df60-410e-8ece-0f46afea968b
md"""
1
$(@bind nr html"<input type=range min=1. max=100. step=1 value=10>")
 100
"""

# ╔═╡ 814f26e1-9e46-42be-8897-b443585e9e98
md"""
As desired, the larger $dt$ record has a larger standard deviation that agrees with the theory:

fine-grained: $\sigma = \sqrt{dt}$ = $(round(sqrt(dt),digits=3))

coarse-grained: $\sigma' = \sqrt{n \hspace{1mm} dt}$ = $(round(sqrt(dt*nr),digits=3))
"""

# ╔═╡ a91a41ba-6f89-43d5-a4d2-7fc0256d3cb6
md"""
nr = $nr
"""

# ╔═╡ 5f98cd8f-2f26-4fd3-b7fb-a75ee4469875
md"""
## Systematic comparison of rouchon + bayesian
"""

# ╔═╡ a504e2af-58ab-4819-af4c-fdcecc0ad0b4
md"""
In order to systematically compare `rouchon` and `bayesian`, we'll first generate 100 simulated records using `bayesian`:
"""

# ╔═╡ 25512701-e9d4-430e-8dd2-11716285d49b
ensb = ensemble(bayesian, T, ρ0, H, J, C; dt=dt, N=100, heterodyne=false)

# ╔═╡ 589ba5a1-2418-4ae5-b9c9-cc073f90c5bb
records = collect(map(rec -> collect(rec[1]) , ensb[3]))

# ╔═╡ 852deadb-3c98-4c61-b9d0-6e4bab29bb92
collect(ensb[2])

# ╔═╡ 41936367-19c9-48ef-bbec-2bed2bdd059f
sols

# ╔═╡ 1aa035ea-5ced-459f-8fd9-82e37e592d12
md" #### RMS error"

# ╔═╡ f8854877-3cd7-4b4f-a534-a7269b153462
md"""
Find the root-mean-square (RMS) error between the trajectories generated from coarse-grained records at scale dt and the original trajectories.
"""

# ╔═╡ 2da3a708-9bec-4207-889f-3632a3fa940b
nlist = [1,2,5,10,15,20,30,40,50,80,100]

# ╔═╡ 3db7293a-ef8f-4ec5-a894-3805514c2210
dydt=[[1,2,3],[5,6,7]]

# ╔═╡ 69f145d7-758b-479a-801a-483e1f823074
dydt./dt

# ╔═╡ e7ffdfb4-622e-4248-a71b-907ea1cbb8a5
md"""

### References

[1] P. Rouchon and J. F. Ralph, Efficient Quantum Filtering for Quantum Feedback Control, Phys. Rev. A 91, 012118 (2015).

"""

# ╔═╡ 2b7aa36e-e64a-4839-875e-2c472763cb80
md" ## Utilities "

# ╔═╡ 7baa7a79-090e-4318-982f-6c1982f82a58
rms(ser1::Array, ser2::Array) = sqrt(sum((ser1 - ser2).^2))/length(ser1)

# ╔═╡ 087145d6-81d4-44ef-b2a6-58c5126081ee
expects = ρ -> collect(real(expect(ρ, s)) for s in [σx,σy,σz,ρ]) # ρ -> [<x>,<y>,<z>,<ρ>]

# ╔═╡ 078c6b8e-4b74-46dd-aebd-77acf392c2c9
function blochs(sol)
	(tt, ρt, _) = sol

	# Get Bloch components
	evs0 = expects.(ρt);
	xx,yy,zz,ρρ = [map(x -> x[i], evs0) for i in 1:4];

	(collect(tt), xx, yy, zz, ρρ)
	
end

# ╔═╡ a9102fda-a602-45e2-ad7a-6fe192a7972a
purity(x,y,z) = 0.5*(1 + x^2 + y^2 + z^2)

# ╔═╡ 9c95461e-4243-447e-b412-7814ca18da43
R(x,y,z) = sqrt(x^2 + y^2 + z^2)

# ╔═╡ b13e4f3d-df37-4a75-8b00-2d0f768bd18e
function subseries(rec, T, dt; scale=2)
	smooth_rec = real(coarse_grain(rec; n=scale))
	
	dtt = dt * scale
	tts = collect(range(first(T), last(T), step=dtt))
	subrec = subselect(smooth_rec; n=scale)
	
	if length(tts) > length(subrec)
		push!(subrec,subrec[end])
		
	elseif length(tts) < length(subrec)
		deleteat!(subrec,length(subrec))
		
	end
	
	(tts, subrec)
	
end

# ╔═╡ 559d34ed-fcab-46d8-9fea-1774067b2d48
(tts, Rss) = subseries(Rs, T, dt; scale=n)

# ╔═╡ 43717946-c7f2-4f57-b516-9db669167587
(tb,xb,yb,zb,ρb) = map(ser -> subseries(ser, T, dt; scale=ns)[2], blochs(solb))

# ╔═╡ 91c7a5dd-f2e1-425b-9008-43854097f966
(ttc, Rsc) = subseries(Rs, T, dt; scale=ns)

# ╔═╡ 7cc75874-ce1f-49da-860b-62a0f7fe800f
solb_coarse = bayesian(T, ρ0, H, J, C; dt=dt*ns, heterodyne=false, dy=[Rsc])

# ╔═╡ aad7d7ba-b3c5-4ce2-99f5-60c13db14969
(tc,xc,yc,zc,ρc) = blochs(solb_coarse)

# ╔═╡ d4f586f6-0e6d-4e96-9426-e00f3eee17da
begin
	(ttr, Rsr) = subseries(Rs, T, dt; scale=nr)
	dyc = Rsr*(dt*nr)
end

# ╔═╡ e4beb554-bc47-477d-ab3a-48432b6b7677
solr_coarse = rouchon(T, ρ0, H, J, C; dt=dt*nr, dy=[dyc])

# ╔═╡ 64710b6d-0238-4f31-a3d0-dbdcbcd97507
(trc,xrc,yrc,zrc,ρrc) = blochs(solr_coarse)

# ╔═╡ 57ff41a2-15d3-4c31-8ce1-ec2bca82a504
(tr,xr,yr,zr,ρr) = map(ser -> subseries(ser, T, dt; scale=nr)[2], blochs(solr))

# ╔═╡ 9b2492db-94b4-470f-89f7-370ac8b5db26
function rms(solve, sols, T, ρ0, H, J, C; n=1, dt=1e-4, kwargs...)
	
	"""
	Arguments
	
	solve :: Function -- bayesian or rouchon
	dt :: time-step of reference solutions
	T :: time duration of reference solutions
	sols :: reference solutions
	n :: scale of coarse-graining
	
	Returns
	
	mean(rms_list) :: average rms between reference trajectories and corresponding
	                  coarse-grained generated trajectories
	
	"""
	
	(tt, ρlist, recs) = sols
	recs = collect(map(rec -> collect(rec[1]) , recs))
	
	rms_list = []
	
	for i in 1:length(recs)
		
		# look at one solution at a time
		Rs = recs[i]
		sol = (tt, ρlist[i], Rs)

		# first, get bloch trajectories for each of the ρ in ρs, and coarse-grain
		(ts,xs,ys,zs,rs) = map(ser -> subseries(ser, T, dt; scale=n)[2], blochs(sol))
		
		# coarse-grain the record
		(_, Rsc) = subseries(Rs, T, dt; scale=n)
		
		# generate a trajectory from coarse-grained record, and get blochs
		sol_coarse = solve(T, ρ0, H, J, C; dt=dt*n, dy=[Rsc], kwargs...)
		(_,xc,yc,zc,rc) = blochs(sol_coarse)
		
		# calculate rms error for each bloch coordinate, then take total
		push!(rms_list, sqrt(rms(xs, xc)^2 + rms(ys, yc)^2 + rms(zs, zc)^2))
		
	end

	# return average rms error
	mean(rms_list)


end

# ╔═╡ 462bad0d-3706-435b-9918-b1b1c2b5e320
md"""
The root-mean-squared error between the two series is $(round(rms(yb, yc), digits=5)).
"""

# ╔═╡ abce17a3-cc09-43ea-ad5e-7606e05d9d69
md"""
The root-mean-squared error between the two series is $(round(rms(yr, yrc), digits=5)).
"""

# ╔═╡ 98be356a-7afa-4387-a9ae-57542e808b1b
begin
	rms_series = []
	
	for n in nlist
		
		rms_n = rms(bayesian, ensb, T, ρ0, H, J, C; n=n, dt=dt, heterodyne=false)
		push!(rms_series, rms_n)
		
	end
end

# ╔═╡ 87664eb9-2a30-4340-a733-a8392bf3dbab
rms_bayes = rms_series

# ╔═╡ 3e2b91c7-97cb-4a59-96d5-01081f9d7a6f
ylabel

# ╔═╡ 37ac2cb1-c957-4182-be1f-29521aafbaa3
# colorscheme
begin
	colorscheme = "Paired"
	cmap = plt.matplotlib.cm.get_cmap(colorscheme)
	colors=collect(map(x -> cmap(x), 0:11))
	md" `colorscheme`"
end

# ╔═╡ 6163347e-7044-4ab7-8683-0a7186d9ac62
begin
	close("all")
	p = plot(nlist, rms_series, color=colors[2], marker="o")
	ax = gca()
	
    xlabel("n")
    ylabel("rms error")
    title(string("rms error for bayesian method, compared to ", L"$dt = 1$", " ns"))
    legend()
    gcf()
	
end

# ╔═╡ 7672f11f-c5f2-4f5e-8ffc-3f7b45b4a322
# Plotting
function plot_solution(sol; plot_title="Rabi Oscillation")
	
	close("all")
	
	(tt, ρt, dys) = sol
    
    # Get Bloch components
    evs0 = expects.(ρt);
    xx,yy,zz,ρρ = [map(x -> x[i], evs0) for i in 1:4];
    
    # Plot Bloch components vs. time
    
    p = plot(tt, xx, color=colors[2], label=L"$x$")
    plot(tt, yy, color=colors[4],label=L"$y$")
    ax = gca()
    ax.set_ylim([-1.1,1.1]) 
    plot(tt, zz, color=colors[6], label=L"$z$")
	plot(tt, ρρ, color=colors[8], label=L"Tr $\rho^2$")
    xlabel(L"$t$")
    ylabel("Bloch coordinates")
    title(plot_title)
	legend()
	gcf()

end

# ╔═╡ af608dc7-55fc-4bdd-bef7-52e6a7eefcb8
plot_solution(solb; plot_title="monitored Rabi oscillations; bayesian method")

# ╔═╡ 9c78b702-3734-4341-92dd-b475e4dd717c
plot_solution(solr; plot_title="monitored Rabi oscillations; rouchon method")

# ╔═╡ 08a6190e-aa1a-485f-bf30-e059e56048e5
# Plotting
function plot_blochs((tt, blochs); plot_title="Rabi Oscillation")
	close("all")
	
	(xx,yy,zz) = blochs
	pp = map(i -> purity(xx[i], yy[i], zz[i]), 1:length(tt))
    
    # Plot Bloch components vs. time
    
    p = plot(tt, xx, color=colors[2], label=L"$x$")
    plot(tt, yy, color=colors[4],label=L"$y$")
    ax = gca()
    ax.set_ylim([-1.1,1.1]) 
    plot(tt, zz, color=colors[6], label=L"$z$")
	plot(tt, pp, color=colors[8], label=L"Tr $\rho^2$")
    xlabel(L"$t$")
    ylabel("Bloch coordinates")
    title(plot_title)
    legend()
    gcf()
end

# ╔═╡ 4fd5a499-61e5-42ac-9d9d-38f461f59616
function record_histograms(records...; plot_title="record histogram", labels=[]::Array, density=false)
	close("all")
	
	μσs = []
	hist_colors = []
	hist_labels = []
	
	for i in 1:length(records)
		label = i > length(labels) ? i : labels[i]
		dys = records[i]
		
		# get mean and std dev for real part
		(μ, σ) = params(fit(Normal, dys))
		push!(μσs, map(p -> round(p, digits=4), (μ, σ)))

		# make histogram		
		n, bins, patches = hist(dys, 50, density=density, 					 									facecolor=colors[2i], alpha=1, label=label)
		push!(hist_colors, colors[2i])
		

	end
	
	# write down (μ, σ) pairs as text boxes
	μσ_strings = map(μσ -> string("(μ, σ) = (", μσ[1], ", ", μσ[2], ")\n"), μσs)
	ax = gca()
	for i in 1:length(μσ_strings)
		str = μσ_strings[i]

		ax.text(0.05, 1 - 0.05i, str, transform=ax.transAxes, fontsize=10,
			verticalalignment="top", color=hist_colors[i])
		
	end
	
	
	legend()
	xlabel("value (arbitrary units)")
	ylabel(density ? "relative frequency" : "frequency")
	title("record histograms")
	gcf()
		
end

# ╔═╡ ca72e503-c34c-4627-aa6f-ddc6aa3292aa
record_histograms(Rs, Rss; labels=["dt = $dt μs", "dt = $(dt*n) μs"], density=true)

# ╔═╡ d6d8de9a-6f61-4fdc-ab81-c433303ae452
record_histograms(dy, dyc; labels=["dt = $dt μs", "dt = $(dt*nr) μs"], density=true)

# ╔═╡ 82f2e0a1-11ee-4aa6-a63f-e6ca170a8116
# Plotting
function plot_timeseries(tt::Array, series...; plot_title="time series", xlab=L"$t$", ylab="arbitrary units", labels=[]::Array, colorpairs=false, kwargs...)
	close("all")
	ser_colors(i) = colorpairs ? colors[i] : colors[2i]
	
	label(i) = i > length(labels) ? i : labels[i]
    
    # Plot records vs. time
	ser = series[1]
	p = plot(tt, ser, color=ser_colors(1), label=label(1), kwargs...)
	ax = gca()
	
	if length(series) > 1
		
		for i in 2:length(series)
			ser = series[i]
			plot(tt, ser, color=ser_colors(i),label=label(i), kwargs...) 
		end
		
	end
	
    xlabel(xlab)
    ylabel(ylab)
    title(plot_title)
    legend()
    gcf()
	
end

# ╔═╡ 5156e99d-7583-40fc-aa4f-73df81cbcedc
# Plotting
function plot_timeseries(ttseries...; plot_title="time series", xlab=L"$t$", ylab="arbitrary units", labels=[]::Array, colorpairs=false)
	close("all")
	ser_colors(i) = colorpairs ? colors[i] : colors[2i]
	
	label(i) = i > length(labels) ? i : labels[i]
    
    # Plot records vs. time
	(tt, ser) = ttseries[1]
	p = plot(tt, ser, color=ser_colors(1), label=label(1))
	ax = gca()
	
	if length(ttseries) > 1
		
		for i in 2:length(ttseries)
			(tt, ser) = ttseries[i]
			plot(tt, ser, color=ser_colors(i),label=label(i)) 
		end
		
	end
	
    xlabel(xlab)
    ylabel(ylab)
    title(plot_title)
    legend()
    gcf()
	
end

# ╔═╡ a55e5440-8a5c-484d-83db-2ee32495f89b
plot_timeseries((tt,Rs),(tts,Rss); plot_title="records", labels=["dt = $dt μs", "dt = $(dt*n) μs"])

# ╔═╡ 576b3970-b370-4798-9e04-13321d97d587
plot_timeseries(tc, yc, yb, zc, zb; plot_title="bayesian comparison", ylab="Bloch coordinates", labels=pairs, colorpairs=true)

# ╔═╡ 0eb75fe9-6123-4a2b-b8c9-196cc46cf0de
plot_timeseries((tt,dy),(ttr,dyc); plot_title="records", labels=["dt = $dt μs", "dt = $(round(dt*nr, digits=3)) μs"])

# ╔═╡ 10ac1593-b4d1-4774-a158-488a0cf65448
plot_timeseries(trc, yrc, yr, zrc, zr; plot_title="rouchon comparison", ylab="Bloch coordinates", labels=pairs, colorpairs=true)

# ╔═╡ e0df4349-7035-4c41-995c-0c4d78069636
# Plotting
function plot_solutions((sol1,sol2); plot_title="Rabi Oscillation")
    close("all")
    
    tt1 = sol1[1]
    ρt1 = sol1[2]
    tt2 = sol2[1]
    ρt2 = sol2[2]
    
    # Get Bloch components
    evs1 = expects.(ρt1);
    x1,y1,z1,ρ1 = [map(x -> x[i], evs1) for i in 1:4];
    evs2 = expects.(ρt2);
    x2,y2,z2,ρ2 = [map(x -> x[i], evs2) for i in 1:4];
    
    # Plot Bloch components vs. time
	p = plot(tt1, x1, color=colors[2], linewidth=2, label=L"$x_{\eta = 0}$")
    plot(tt1, y1, color=colors[4], linewidth=2, label=L"$y_{\eta = 0}$")
    plot(tt1, z1, color=colors[6], linewidth=2, label=L"$z_{\eta = 0}$")
    plot(tt1, p1, color=colors[8], linewidth=2, label=L"(Tr $\rho^2)_{\eta = 0}$")
    plot(tt2, x2,  color=colors[1], linestyle="dashed", label=L"$x_{avg}$")
    plot(tt2, y2, color=colors[3], linestyle="dashed", label=L"$y_{avg}$")
    plot(tt2, p2, color=colors[7], linewidth=2, linestyle="dashed", label=L"(Tr $\rho^2)_{avg}$")
  
	ax = gca()
	ax.set_ylim([-1.1,1.1]) 
	xlabel(L"$t$")
	ylabel("Bloch coordinates")
	title(plot_title)
	legend()
	gcf()
end

# ╔═╡ b3e6866b-6ab7-4f33-be5d-6e28da5f5fe7
function plot_evals((tt1, evals); α=0.1, linewidth=1, labels=false)
    xxs,yys,zzs,ρρs = [map(x -> x[i], evals) for i in 1:4];
    if labels
        plot(tt1, xxs, color=colors[2], alpha=α, linewidth=linewidth, label=L"$x$")
        plot(tt1, yys, color=colors[4], alpha=α, linewidth=linewidth, label=L"$y$")
        plot(tt1, zzs, color=colors[6], alpha=α, linewidth=linewidth, label=L"$z$")
        plot(tt1, ρρs, color=colors[8], alpha=α, linewidth=linewidth, label=L"Tr $ρ^2$")
    else
        plot(tt1, xxs, color=colors[2], alpha=α, linewidth=linewidth)
        plot(tt1, yys, color=colors[4], alpha=α, linewidth=linewidth)
        plot(tt1, zzs, color=colors[6], alpha=α, linewidth=linewidth)
        plot(tt1, ρρs, color=colors[8], alpha=α, linewidth=linewidth)
        
    end

end

# ╔═╡ b6577418-9cf1-4232-9fdb-fb6ebf6efd22
function plot_ensemble(sol_ens; α=0.1, linewidth=1, labels=false, average=false)
    close("all")
	tt1 = sol_ens[1]
    evs = collect(map(ρs -> expects.(ρs), sol_ens[2]));

    for i in 1:50
        plot_evals((tt1, evs[i]); α=α, labels=labels, linewidth=linewidth)
    end

    if average
        plot_evals((tt1, mean(evs)), α=1, linewidth=1.5, labels=true)
        title_string = "Trajectories w/ ensemble average"
    else
        plot_evals((tt1, evs[1]), α=1, linewidth=1.5, labels=true)
        title_string = "Trajectories"
    end
    
    ax = gca()
    ax.set_ylim([-1.1,1.1]) 
	legend()

    xlabel(L"$t$")
    ylabel("Bloch coordinates")
    title(title_string)
	
	gcf()
end

# ╔═╡ 093383f1-d7a5-48c6-9927-4e42d158cc3d
green(text; title="Note") = Markdown.MD(Markdown.Admonition("correct", title, [text]))

# ╔═╡ 5ed525a9-c4b6-44d6-8ebf-6c86190a8992
red(text; title="Note") = Markdown.MD(Markdown.Admonition("danger", title, [text]))

# ╔═╡ b966b51c-e9cd-4696-a3b1-15ecb9d3808e
tan(text; title="Note") = Markdown.MD(Markdown.Admonition("warning", title, [text]))

# ╔═╡ f08bf969-d03f-47f9-9f5a-cddd18f8b109
blue(text; title="Note") = Markdown.MD(Markdown.Admonition("note", title, [text]))

# ╔═╡ d7c1fd28-5780-4efd-b085-4f1f797a3f09
hint(text; title="Hint") = Markdown.MD(Markdown.Admonition("hint", title, [text]))

# ╔═╡ Cell order:
# ╠═80ffb52c-a0dc-43d9-8282-e8acb51df4e0
# ╟─595fc42b-943e-4a29-841a-cd36c90a2b55
# ╟─483e648c-2ac3-46f7-b49b-a3109deec27d
# ╟─52a501dd-c011-41cb-8b47-ce5ff894761c
# ╟─5d30c7c1-4c53-40d2-8e73-b58db5697cfa
# ╟─d4e7b79b-eb8c-48b3-9a94-650481fa8cb7
# ╟─bef8c648-af4c-46c1-87b4-29f744f9aaf3
# ╟─2059a835-925d-430d-aa0c-283a03f6ba2d
# ╟─476f0b08-7d21-4010-b114-130e6dfbbae0
# ╠═639931ed-0520-499c-ba1c-90b04e854ebd
# ╟─68627f37-30c6-47c2-85e3-a031cf4a2e05
# ╠═0c3770b6-02b2-46af-a875-cea82999f88f
# ╠═af608dc7-55fc-4bdd-bef7-52e6a7eefcb8
# ╠═eb250fd3-7c43-4f75-ab9e-6f7ae0035ad3
# ╟─75f1cc94-5222-42d5-8c0e-ceedb3c53d3b
# ╠═bdf4c000-2845-4e2f-baf6-0d1161e477e3
# ╠═9c78b702-3734-4341-92dd-b475e4dd717c
# ╟─b6c73eea-d17c-49d4-a091-be10c8134718
# ╟─4fc42181-a1d4-46d9-be76-7277b9a2d6a5
# ╟─343161e9-e975-45b5-a444-b9f9f0b30cfb
# ╟─d5af31e4-ba93-4f3f-9c55-4f2ce600d797
# ╟─97e2560f-3dd9-45a0-85f5-aa9e12800292
# ╟─d6d6e78b-15cb-4c6a-af30-8c42b35dc1b9
# ╟─6315198d-4962-4686-b436-80a3525e0dca
# ╟─40a4c411-745d-4ac9-8cd8-ab4a66979fa9
# ╟─a306838f-f883-4fe9-94e6-8d5771cd9793
# ╟─f66192ba-ff34-499e-ae8e-e8d08b91b9a1
# ╟─a541146f-9143-4a91-bb3a-be17e9bb6561
# ╟─88179c16-f111-462d-a9d6-72e3764715eb
# ╟─6e92e073-d122-4a84-8c4c-b19858c173a9
# ╟─501a6621-8054-4fa3-b9a7-f1a921e14a56
# ╟─cb87a45d-4ff1-4da9-8bbe-ec2c193cb7c6
# ╟─676e029d-1141-41b2-bc81-5d2e4500c684
# ╟─df840a53-c8b0-4215-97de-ef6f1c720ed7
# ╟─7ffe3ce0-2ce5-4012-8679-93b9154b150e
# ╟─e41a9229-34f2-4a42-810b-7c95dc02780b
# ╠═4f3c3104-c64d-4b77-8932-c37f4aaabaa1
# ╟─3e0a12d2-b4bf-4424-8b62-d9a714e393fc
# ╠═559d34ed-fcab-46d8-9fea-1774067b2d48
# ╟─69eaed7c-392d-40d8-9f0c-63b4bc4aa18d
# ╟─0b0f2aa0-431d-4029-b7c2-911aa4d2f69a
# ╟─012e9206-f306-4776-a546-d6f7dcf00ef2
# ╠═a55e5440-8a5c-484d-83db-2ee32495f89b
# ╠═ca72e503-c34c-4627-aa6f-ddc6aa3292aa
# ╟─275ffd8d-274b-4fab-946e-d9482fa340a6
# ╟─e0cb7572-3d67-47e0-9482-70588a693d41
# ╟─f870609e-f3bc-40e4-90fb-484e7b52db4b
# ╠═43717946-c7f2-4f57-b516-9db669167587
# ╟─3f7d824c-8a0e-4547-bb06-81717a01df40
# ╠═91c7a5dd-f2e1-425b-9008-43854097f966
# ╠═7cc75874-ce1f-49da-860b-62a0f7fe800f
# ╠═aad7d7ba-b3c5-4ce2-99f5-60c13db14969
# ╟─7c434d23-a5ed-41bb-8dc7-8b1440a0680b
# ╟─e47ee9c6-c9d6-44b8-8fc3-e6ebb58f28e3
# ╟─82b18884-28d4-4ecb-b21c-d3e57dd64870
# ╟─462bad0d-3706-435b-9918-b1b1c2b5e320
# ╠═576b3970-b370-4798-9e04-13321d97d587
# ╠═82d43f38-9cf6-4fee-bdc8-502234781678
# ╟─499291df-d136-4dab-8a6d-0d72cf812435
# ╟─3d3d4838-2c53-44b5-a103-5a7656d58611
# ╠═d4f586f6-0e6d-4e96-9426-e00f3eee17da
# ╠═0eb75fe9-6123-4a2b-b8c9-196cc46cf0de
# ╠═d6d8de9a-6f61-4fdc-ab81-c433303ae452
# ╟─814f26e1-9e46-42be-8897-b443585e9e98
# ╟─b9209585-d439-488e-b65f-d49a61d53bd8
# ╠═57ff41a2-15d3-4c31-8ce1-ec2bca82a504
# ╠═e4beb554-bc47-477d-ab3a-48432b6b7677
# ╠═64710b6d-0238-4f31-a3d0-dbdcbcd97507
# ╟─a91a41ba-6f89-43d5-a4d2-7fc0256d3cb6
# ╟─7e1f7d9e-df60-410e-8ece-0f46afea968b
# ╠═abce17a3-cc09-43ea-ad5e-7606e05d9d69
# ╠═10ac1593-b4d1-4774-a158-488a0cf65448
# ╟─5f98cd8f-2f26-4fd3-b7fb-a75ee4469875
# ╟─a504e2af-58ab-4819-af4c-fdcecc0ad0b4
# ╠═25512701-e9d4-430e-8dd2-11716285d49b
# ╠═589ba5a1-2418-4ae5-b9c9-cc073f90c5bb
# ╠═852deadb-3c98-4c61-b9d0-6e4bab29bb92
# ╠═41936367-19c9-48ef-bbec-2bed2bdd059f
# ╟─1aa035ea-5ced-459f-8fd9-82e37e592d12
# ╠═f8854877-3cd7-4b4f-a534-a7269b153462
# ╠═2da3a708-9bec-4207-889f-3632a3fa940b
# ╠═98be356a-7afa-4387-a9ae-57542e808b1b
# ╠═87664eb9-2a30-4340-a733-a8392bf3dbab
# ╠═6163347e-7044-4ab7-8683-0a7186d9ac62
# ╠═9b2492db-94b4-470f-89f7-370ac8b5db26
# ╠═3db7293a-ef8f-4ec5-a894-3805514c2210
# ╠═69f145d7-758b-479a-801a-483e1f823074
# ╟─e7ffdfb4-622e-4248-a71b-907ea1cbb8a5
# ╟─2b7aa36e-e64a-4839-875e-2c472763cb80
# ╠═7baa7a79-090e-4318-982f-6c1982f82a58
# ╠═078c6b8e-4b74-46dd-aebd-77acf392c2c9
# ╟─087145d6-81d4-44ef-b2a6-58c5126081ee
# ╟─a9102fda-a602-45e2-ad7a-6fe192a7972a
# ╟─9c95461e-4243-447e-b412-7814ca18da43
# ╟─7672f11f-c5f2-4f5e-8ffc-3f7b45b4a322
# ╟─08a6190e-aa1a-485f-bf30-e059e56048e5
# ╟─4fd5a499-61e5-42ac-9d9d-38f461f59616
# ╟─b13e4f3d-df37-4a75-8b00-2d0f768bd18e
# ╠═82f2e0a1-11ee-4aa6-a63f-e6ca170a8116
# ╠═5156e99d-7583-40fc-aa4f-73df81cbcedc
# ╠═3e2b91c7-97cb-4a59-96d5-01081f9d7a6f
# ╟─e0df4349-7035-4c41-995c-0c4d78069636
# ╟─b3e6866b-6ab7-4f33-be5d-6e28da5f5fe7
# ╟─b6577418-9cf1-4232-9fdb-fb6ebf6efd22
# ╟─37ac2cb1-c957-4182-be1f-29521aafbaa3
# ╟─093383f1-d7a5-48c6-9927-4e42d158cc3d
# ╟─5ed525a9-c4b6-44d6-8ebf-6c86190a8992
# ╟─b966b51c-e9cd-4696-a3b1-15ecb9d3808e
# ╟─f08bf969-d03f-47f9-9f5a-cddd18f8b109
# ╟─d7c1fd28-5780-4efd-b085-4f1f797a3f09
