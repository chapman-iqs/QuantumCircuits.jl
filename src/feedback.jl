abstract type Feedback end

struct NoFeedback <: Feedback end
const nofeedback = NoFeedback()

struct SimpleTimeDelay <: Feedback end
const simpletimedelay = SimpleTimeDelay()

struct HamiltonianFE <: Feedback end
const hamiltonianfe = HamiltonianFE()


# these cases are included for integration with QuantumSystems
bayesian(::NoFeedback, args...; kwargs...) = bayesian(args...; kwargs...)
trajectory(::NoFeedback, args...; kwargs...) = trajectory(args...; kwargs...)

function trajectory(::HamiltonianFE, ts, (ρs0, ρf0), (Us, Uf), (ls, lf), (ms, mf); dt=1e-3, r0=[0.0], td = 0.0, forwardestimate=true)

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
    Id = identityoperator(Us(0.0, ρs0, r0))
    buffermaps = [Id for i in 1:(delay-1)]
    ubuffer = Id

    # initialize readout
	readouts = Readouts([r0])

    for i in 2:length(ts)
        t = ts[i]

        # delayed filter state
        # before td is reached, we just use the initial filter state
        ρf_delayed = i <= delay ? ρf0 : ρf_list[i - delay]
		r_delayed = i <= delay ? r0 : readouts[i - delay]

        # get the forward-estimated state as a function of the delayed filter state and the intervening maps
        ρe = update(ubuffer, ρf_delayed)

        # get the unitaries for this timestep, and update buffer
		ρeff = forwardestimate ? ρe : ρf_delayed
        us = Us(t, ρeff, r_delayed)
        uf = Uf(t, ρeff, r_delayed)
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

function bayesian(feedback::HamiltonianFE, (t0, tf), ρ, (Hs, Hf), (Js, Jf), C; dt=1e-3, td = 0.0, forwardestimate=true)
	Hs, Hf = convertham.((Hs, Hf))
    # currently, Hs and Hf are both assumed to be functions of (t::Timescale, ρ::State)

	ts = range(t0, tf, step=dt)
	r0::Readout = [0.0 for i in 1:length(C)] # initial, empty Readout object

    # define unitary operators (NOT superoperators) -- since we'll need to track these in a buffer
    Us(t::Timescale, ρe::State, r::Readout) = exp( -im * dt * DenseOperator(Hs(t, ρe, r)))
    
    # filter
    Uf(t::Timescale, ρe::State, r::Readout) = exp( -im * dt * DenseOperator(Hf(t, ρe, r)))

    # define lindblad / measurement superoperator functions (not matrices, for efficiency)
	ms = meas(dt, C; sample=true) 	# "system" measurement function (generates a readout based on system state)
    mf = meas(dt, C; sample=false) 	# "filter" measurement function (ALSO generates a readout based on system state)
	ls = lind(dt, Js) # decay function for "system"
	lf = lind(dt, Jf) # decay function for "filter"

    # type conversion
    ρs0 = length(Js) > 0 && typeof(ρ) <: Ket ? dm(ρ) : ρ
    ρf0 = length(Jf) > 0 && typeof(ρ) <: Ket ? dm(ρ) : ρ

    trajectory(feedback, ts, (ρs0, ρf0), (Us, Uf), (ls, lf), (ms, mf); dt=dt, r0=r0, td=td, forwardestimate=forwardestimate)
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