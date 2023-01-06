abstract type Feedback end
abstract type SystemFilter <: Feedback end
abstract type ForwardEstimation <: SystemFilter end

struct NoFeedback <: Feedback end
struct SimpleTimeDelay <: Feedback end
struct SystemFilter2H <: SystemFilter end
struct SystemFilter2H2L <: SystemFilter end
struct HamiltonianFE <: ForwardEstimation end
struct HamLindFE <: ForwardEstimation end

const nofeedback = NoFeedback()
const simpletimedelay = SimpleTimeDelay()
const hamiltonianfe = HamiltonianFE()
const hamlindfe = HamLindFE()
const systemfilter2h = SystemFilter2H()
const systemfilter2h2l = SystemFilter2H2L()

# these cases are included for integration with QuantumSystems
bayesian(::NoFeedback, args...; kwargs...) = bayesian(args...; kwargs...)
trajectory(::NoFeedback, args...; kwargs...) = trajectory(args...; kwargs...)

function trajectory(::HamiltonianFE, ts, (ρs0, ρf0), (Us, Uf), (ls, lf), (ms, mf); dt=1e-3, r0=[0.0], td = 0.0)

	# initialize system
    ρs = ρs0
    ρs_list = [ρs]

    # initialize filter
    ρf = ρf0
    ρf_list = [ρf]

    # initialize forward-estimated state
    ρe = ρf0
    ρe_list = [ρe]

    # delay and buffer maps
    delay = Int64(round(td / dt)) + 1
    Id = identityoperator(Us(0.0, ρs0))
    buffermaps = [Id for i in 1:(delay-1)]
    ubuffer = Id

    # initialize readout
	readouts = Readouts([r0])

    for i in 2:length(ts)
        t = ts[i]

        # delayed filter state
        # before td is reached, we just use the initial filter state
        ρf_delayed = i <= delay ? ρf0 : ρf_list[i - delay] 

        # old version
        # foldr does "last applied first" -- so the most recent unitary is first in the list, and is applied last
        # get the forward-estimated state as a function of the delayed filter state and the intervening maps
        # ρe = foldr(update, buffermaps; init=ρf_delayed)

        # get the forward-estimated state as a function of the delayed filter state and the intervening maps
        ρe = update(ubuffer, ρf_delayed)

        # get the unitaries for this timestep, and update buffer
        us = Us(t, ρe)
        uf = Uf(t, ρe)
        pushfirst!(buffermaps, uf)
        uf_last = pop!(buffermaps)

        # update buffer; uf gets tacked on at end, and (uf_last)' at beginning to cancel it out
        ubuffer = uf_last' * ubuffer * uf

        # state updates: system
        ρs, ro = ms(t, ρs)
        ρs = update(us, ρs)
        ρs = ls(t, ρs)

        # state updates: filter
        ρf = mf(t, ρf, ro) # takes system readout to get measurement operator; does not depend on ρf, only updates it
        ρf = update(uf, ρf)
        ρf = lf(t, ρf)
        
		push!(ρs_list, ρs)
        push!(ρf_list, ρf)
        push!(ρe_list, ρe)
        push!(readouts, ro)
	end
    
    records = trans(readouts)
	return Solution(collect(ts), ρs_list, records), Solution(collect(ts), ρf_list), Solution(collect(ts), ρe_list)
end

function bayesian(feedback::HamiltonianFE, (t0, tf), ρ, (Hs, Hf), (Js, Jf), C; dt=1e-3, td = 0.0)
	# Hs, Hf = convertham.(Hpair)
    # currently, Hs and Hf are both assumed to be functions of (t::Timescale, ρ::State)

	ts = range(t0, tf, step=dt)
	r0::Readout = [0.0 for i in 1:length(C)] # initial, empty Readout object

    # define unitary operators (NOT superoperators) -- since we'll need to track these in a buffer
    # system
    # Us() = exp( -im * dt * DenseOperator(Hs)
    # Us(t::Timescale) = exp( -im * dt * DenseOperator(Hs(t)))
    Us(t::Timescale, ρe::State) = exp( -im * dt * DenseOperator(Hs(t, ρe)))
    
    # filter
    # Uf() = exp( -im * dt * DenseOperator(Hf)) 
    # Uf(t::Timescale) = exp( -im * dt * DenseOperator(Hf(t)))
    Uf(t::Timescale, ρe::State) = exp( -im * dt * DenseOperator(Hf(t, ρe)))

    # define lindblad / measurement superoperator functions (not matrices, for efficiency)
	ms = meas(dt, C; sample=true) 	# "system" measurement function (generates a readout based on system state)
    mf = meas(dt, C; sample=false) 	# "filter" measurement function (ALSO generates a readout based on system state)
	ls = lind(dt, Js) # decay function for "system"
	lf = lind(dt, Jf) # decay function for "filter"

    # type conversion
    ρs0 = length(Js) > 0 && typeof(ρ) <: Ket ? dm(ρ) : ρ
    ρf0 = length(Jf) > 0 && typeof(ρ) <: Ket ? dm(ρ) : ρ

    trajectory(feedback, ts, (ρs0, ρf0), (Us, Uf), (ls, lf), (ms, mf); dt=dt, r0=r0, td=td)
end # bayesian

# OLD TRAJECTORY STUFF THAT NEEDS TO BE REFACTORED
#
# currently, feeding in a pre-determined measurement record is not supported for feedback,
# since it's not clear when this would be appropriate.
# NOTE: chose to apply unitary, then lindblad, then measurement, so that feedback is applied directly after measurement obtained
# but should check if this is correct.
function trajectory(::SimpleTimeDelay, m::Function, h::Function, l::Function, ts, ρ; fn::Function=ρ->ρ, dt=1e-4, r0=[0.0], td=0.0)

	# initialize
	ρ0 = ρ
	ρs = [ρ0]
	fnρs = [fn(ρ0)]
	readouts = [r0]

	# if there is no time delay, delay = 1, and we use the record / state from most recent update in the
	# feedback Hamiltonian
	delay = Int64(round(td / dt))

	# before feedback delay buffer is filled, we must feed back on "fake" record or state
	for t in ts[2:(delay + 1)]

		# measurement evolution
		ρ1, ro = m(t, ρ)
		push!(readouts, ro)
		push!(ρs, ρ1)
		push!(fnρs, fn(ρ1))

		# unitary evolution
		ρ2 = h(t, ρ1, ρ0, r0)
		# ρ2 = apply(h, t, ρ1, ρ0, r0)

		# lindblad evolution
		ρ = l(t, ρ2)
	end

	# after feedback delay buffer is filled, we can feed back on acquired readout / filtered state
	for i in (delay + 2):length(ts)
		t = ts[i]

		# measurement evolution
		ρ1, ro = m(t, ρ)
		push!(readouts, ro)
		push!(ρs, ρ1)
		push!(fnρs, fn(ρ1))

		# unitary evolution
		ρd = ρs[i - delay]
		rd = readouts[i - delay]
		ρ2 = h(t, ρ1, ρd, rd)

		# rd = readouts[i - delay]
		# ρd = ρs[i - delay]
		# ρ2 = apply(h, t, ρ1, ρd, rd)

		# lindblad evolution
		ρ = l(t, ρ2)
    end

	records = trans(readouts)

	return Solution(collect(ts), fnρs, records)
end

function trajectory(::SystemFilter2H, (ms, mf), (hs, hf), l::Function, ts, ρ0; fn::Function=ρ->ρ, dt=1e-4, r0=[0.0], td=0.0, readout_integration_bins=1)

	# initialize system
	ρs = ρ0
	ρss = [ρs]
	fnρs = [fn(ρs)]

	# initialize filter
	ρf = ρ0
	ρfs = [ρf]
	fnρf = [fn(ρ0)]

	# shared readout for system and filter
	readouts = [r0]

	# if there is no time delay, delay = 1, and we use the record / state from most recent update in the
	# feedback Hamiltonian
	delay = Int64(round(td / dt)) + 1


	# feed back on filter state; obtain record using system state
	# before feedback delay buffer is filled, we must feed back on "fake" record or state
	for t in ts[2:delay]

		### system evolution
		ρs1 = hs(t, ρs, ρ0, r0)# unitary
		ρs2 = l(t, ρs1) # lindblad
		ρs, ro = ms(t, ρs2) # measurement

		push!(ρss, ρs)
		push!(readouts, ro)
		push!(fnρs, fn(ρs))

		### filter evolution
		ρf1 = hf(t, ρf, ρ0, r0) # unitary
		ρf2 = l(t, ρf1) # lindblad
		ρf = mf(t, ρf2, ro) # measurement (based on filtered state expectation values, but with system readout)

		push!(ρfs, ρf)
		push!(fnρf, fn(ρf))
	end

	# after feedback delay buffer is filled, we can feed back on acquired readout / filtered state
	for i in (delay + 1):length(ts)
		t = ts[i]
		j = i - delay
		ρd = ρfs[j] # state used in feedback is the filtered state from time td ago
		startsmooth = max(1, j - (readout_integration_bins - 1)) # don't index readouts to smooth from an index less than 1
		rd = mean(readouts[startsmooth:j]) # mean(trans(readouts[startsmooth:j])) 	# readout for feedback is an average of the readout_integration_bins readouts
													# before the delayed time index j

		### system evolution
		ρs1 = hs(t, ρs, ρd, rd) # unitary
		ρs2 = l(t, ρs1) # lindblad
		ρs, ro = ms(t, ρs2) # measurement

		push!(ρss, ρs)
		push!(readouts, ro)
		push!(fnρs, fn(ρs))

		### filter evolution
		ρf1 = hf(t, ρf, ρd, rd) # unitary
		ρf2 = l(t, ρf1) # lindblad
		ρf = mf(t, ρf2, ro) # measurement (based on filtered state expectation values, but with system readout)

		push!(ρfs, ρf)
		push!(fnρf, fn(ρf))
    end

	records = trans(readouts)

	return (Solution(collect(ts), fnρs, records), Solution(collect(ts), fnρf, records))
end
function bayesian(feedback::SystemFilter2H, (t0, tf), ρ, (Hs0, Hf0), J, C; fn=ρ->ρ, dt=1e-4, td::Timescale = 0.0, readout_integration_bins=1)

	if !(typeof(Hs0) <: Function) || !(typeof(Hf0) <: Function) || applicable(Hs0, tt) || applicable(Hf0, tt)
		str = string("Non-feedback system-filter Hamiltonians are not yet supported. Please input Hamiltonians ",
						"as a tuple of functions (Hs(args1...), Hf(args2...)), where args1, args2 = (t::Timescale, ρ::State),",
						"(t::Timescale, r::Record), or (t::Timescale, ρ::State, r::Record), and Hs and Hf ",
						"are the system and the filter Hamiltonians, respectively. This feature should ",
						"be added in the future...."
						)
		error(str)
	end

	Hs = convertham(Hs0)
	Hf = convertham(Hf0)

	ts = range(t0, tf, step=dt)
	r0::Readout = [0.0 for i in 1:length(C)] # initial, empty Readout object

	# define superoperator functions for state update
	ms = meas(dt, C; sample=true) 	# "system" measurement superoperator (generates a readout based on system state)
	mf = meas(dt, C; sample=false) 	# "filter" measurement superoperator (ALSO generates a readout based on system state
									# -- same superoperator, but needs to be able to take readout as argument)
	hs = ham(dt, Hs) # unitary superoperator for experimentally "known" Hamiltonian (applied to both "system" and "filter")
	hf = ham(dt, Hf) # unitary superoperator for experimentally "unknown" Hamiltonian (applied to only "system")
	l = lind(dt, J) # decay superoperator

	if length(J) > 0 && typeof(ρ) <: Ket
		ρ = dm(ρ)
	end

	trajectory(feedback::SystemFilter2H, (ms, mf), (hs, hf), l, ts, ρ; fn=fn, r0=r0, dt=dt, td=td, readout_integration_bins=readout_integration_bins)
end # bayesian

function trajectory(::SystemFilter2H2L, (ms, mf), (hs, hf), (ls, lf), ts, ρ0; fn::Function=ρ->ρ, dt=1e-4, r0=[0.0], td=0.0, readout_integration_bins=1)

	# initialize system
	ρs = ρ0
	ρss = [ρs]
	fnρs = [fn(ρs)]

	# initialize filter
	ρf = ρ0
	ρfs = [ρf]
	fnρf = [fn(ρ0)]

	# shared readout for system and filter
	readouts = [r0]

	# if there is no time delay, delay = 1, and we use the record / state from most recent update in the
	# feedback Hamiltonian
	delay = Int64(round(td / dt)) + 1

	# feed back on filter state; obtain record using system state
	# before feedback delay buffer is filled, we must feed back on "fake" record or state
	for t in ts[2:delay]

		### system evolution
		ρs1 = hs(t, ρs, ρ0, r0)# unitary
		ρs2 = ls(t, ρs1) # lindblad
		ρs, ro = ms(t, ρs2) # measurement

		push!(ρss, ρs)
		push!(readouts, ro)
		push!(fnρs, fn(ρs))

		### filter evolution
		ρf1 = hf(t, ρf, ρ0, r0) # unitary
		ρf2 = lf(t, ρf1) # lindblad
		ρf = mf(t, ρf2, ro) # measurement (based on filtered state expectation values, but with system readout)

		push!(ρfs, ρf)
		push!(fnρf, fn(ρf))
	end

	# after feedback delay buffer is filled, we can feed back on acquired readout / filtered state
	for i in (delay + 1):length(ts)
		t = ts[i]
		j = i - delay
		ρd = ρfs[j] # state used in feedback is the filtered state from time td ago
		startsmooth = max(1, j - (readout_integration_bins - 1)) # don't index readouts to smooth from an index less than 1
		rd = mean(readouts[startsmooth:j]) # mean(trans(readouts[startsmooth:j])) 	# readout for feedback is an average of the readout_integration_bins readouts
													# before the delayed time index j

		### system evolution
		ρs1 = hs(t, ρs, ρd, rd) # unitary
		ρs2 = ls(t, ρs1) # lindblad
		ρs, ro = ms(t, ρs2) # measurement

		push!(ρss, ρs)
		push!(readouts, ro)
		push!(fnρs, fn(ρs))

		### filter evolution
		ρf1 = hf(t, ρf, ρd, rd) # unitary
		ρf2 = lf(t, ρf1) # lindblad
		ρf = mf(t, ρf2, ro) # measurement (based on filtered state expectation values, but with system readout)

		push!(ρfs, ρf)
		push!(fnρf, fn(ρf))
    end

	records = trans(readouts)

	return (Solution(collect(ts), fnρs, records), Solution(collect(ts), fnρf, records))
end
function bayesian(feedback::SystemFilter2H2L, (t0, tf), ρ, Hpair::Tuple, Jpair::Tuple, C; fn=ρ->ρ, dt=1e-4, td::Timescale = 0.0, readout_integration_bins=1)
	Hs0, Hf0 = Hpair
	Js, Jf = Jpair

	if !(typeof(Hs0) <: Function) || !(typeof(Hf0) <: Function) || applicable(Hs0, tt) || applicable(Hf0, tt)
		str = string("Non-feedback system-filter Hamiltonians are not yet supported. Please input Hamiltonians ",
						"as a tuple of functions (Hs(args1...), Hf(args2...)), where args1, args2 = (t::Timescale, ρ::State),",
						"(t::Timescale, r::Record), or (t::Timescale, ρ::State, r::Record), and Hs and Hf ",
						"are the system and the filter Hamiltonians, respectively. This feature should ",
						"be added in the future...."
						)
		error(str)
	end

	Hs, Hf = convertham.(Hs0, Hf0)

	ts = range(t0, tf, step=dt)
	r0::Readout = [0.0 for i in 1:length(C)] # initial, empty Readout object

	# define superoperator functions for state update
	ms = meas(dt, C; sample=true) 	# "system" measurement superoperator (generates a readout based on system state)
	mf = meas(dt, C; sample=false) 	# "filter" measurement superoperator (ALSO generates a readout based on system state
									# -- same superoperator, but needs to be able to take readout as argument)
	hs = ham(dt, Hs) # unitary superoperator for experimentally "known" Hamiltonian (applied to both "system" and "filter")
	hf = ham(dt, Hf) # unitary superoperator for experimentally "unknown" Hamiltonian (applied to only "system")
	ls = lind(dt, Js) # decay superoperator for "system"
	lf = lind(dt, Jf) # decay superoperator for "filter"

	if length(Js) > 0 || length(Jf) > 0 && typeof(ρ) <: Ket
		ρ = dm(ρ)
	end

	trajectory(feedback::SystemFilter2H2L, (ms, mf), (hs, hf), (ls, lf), ts, ρ; fn=fn, r0=r0, dt=dt, td=td, readout_integration_bins=readout_integration_bins)
end # bayesian
