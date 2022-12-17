abstract type Feedback end
abstract type ForwardEstimation <: Feedback end
struct NoFeedback <: Feedback end
struct HamiltonianFE <: ForwardEstimation end
struct HamLindFE <: ForwardEstimation end
const nofeedback = NoFeedback()
const hamiltonianfe = HamiltonianFE()
const hamlindfe = HamLindFE()

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

function bayesian(feedback::HamiltonianFE, (t0, tf), ρ, (Hs, Hf)::Tuple, (Js, Jf), C; dt=1e-3, td::Timescale = 0.0)
	# Hs, Hf = convertham.(Hpair)

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

    trajectory(feedback, ts, (ρs0, ρf0), (Us, Uf), (ls, lf), (ms, mf); dt=1e-3, r0=r0, td=td)
end # bayesian


