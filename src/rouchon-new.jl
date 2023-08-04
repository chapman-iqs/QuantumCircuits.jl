"""
Arguments:

H :: Function -- Hamiltonian as a function of time
J :: Vector{Tuple} -- vector of tuples of dissipation operators [(ji, Γi)]
C :: Vector{Tuple} -- vector of tuples of measurement operators [(mi, Γi, ηi)]
"""
rouchon_M(H::Operator, args...; kwargs...) = rouchon_M(t -> H, args...; kwargs...)
function rouchon_M(H::Function, J, C, dt; dys::Function = rouchon_dy(C, dt))
    
    # initialization and conversion
    Id = identityoperator(H(0))

    # deterministic part
    Jsum(t) = mapreduce(+, J; init=0Id) do (j, Γ)        Γ(t) * j(t)' * j(t) end
    msum(t) = mapreduce(+, C; init=0Id) do (m, Γ, η)        η * Γ(t)/2 * m(t)' * m(t)  end

    # stochastic part
    mdy(t, ρ, dWs) = mapreduce(+, zip(C, dys(t, ρ, dWs)); init=0Id) do ((m, Γ, η), dy)
                        sqrt(η * Γ(t)/2) * m(t) * dy
                  end
    mm(t, dWs) = mapreduce(+, 1:length(C), 1:length(C); init=0Id) do r, s
            (mr, Γr, ηr) = C[r]
            (ms, Γs, ηs) = C[s]
            sqrt(ηr * ηs * Γr(t) * Γs(t)) / 4 * mr(t) * ms(t) * (dWs[r] * dWs[s] - δ(r,s) * dt)
    end

    # total   
    Mn(t, ρ, dWs) = Id - (im * H(t) + Jsum(t)/2 + msum(t)/2)*dt + mdy(t, ρ, dWs) + mm(t, dWs)
    return Mn
end

function rouchon_D(J, C, dt)
    Id(ρ) = identityoperator(ρ.basis_l)
    Jsum(t, ρ) = mapreduce(+, J; init=0Id(ρ)) do (j, Γ)        Γ(t) * j(t) * ρ * j(t)'           end
    # msum(ρ) = mapreduce(+, C; init=0Id(ρ)) do (m, Γ, η)     Γ/2 * (1-η) * m * ρ * m'  end
    return (t, ρ) -> Jsum(t, ρ)*dt # + msum(ρ)*dt
end

function rouchon_dy(C, dt)
   return (t, ρ, dWs) -> map(zip(C, dWs)) do ((m, Γ, η), dW)
                            c(t) = sqrt(Γ(t) * η/2) * m(t)
                            real(expect(c(t) + c(t)', ρ)) * dt + dW
                end
end

function record_to_dy(record::Record, (c, Γ, η), dt, ts)
    τms = 1/(2Γ.(ts)*η)
    return map(zip(record, τms)) do (r, τm)     r * dt / sqrt(τm) end
end
function dy_to_record(dys, (m, Γ, η), dt, ts)::Record
    τms = 1/(2Γ.(ts)*η)
    return map(zip(dys, τms)) do (dy, τm)     dy * sqrt(τm) / dt end
end
function records_to_dys(records::Vector{Record}, C, dt, ts)
    map(zip(records, C)) do (record, ctup)
        record_to_dy(record, ctup, dt, ts)
    end
end
function dys_to_records(dys, C, dt, ts)::Vector{Record}
    map(zip(dys, C)) do (dy, ctup)
        dy_to_record(dy, ctup, dt, ts)
    end
end

function convertfunction!(A::Vector{Tuple}, args::Tuple)
    for (i,a) in A
        A[i] = convertfunction!(a, args)
    end
    return A
end
function convertfunction!(a::Tuple, args::Tuple)
    a = map(el -> convertfunction(el), a)
    return a
end
function convertfunction(el, args::Tuple)
    if applicable(el, args...)
        return el
    elseif !isa(el, Function)
        return (anyargs...) -> el
    else
        error("convertfunction not yet implemented for arbitrary argument types (and maybe cannot be implemented arbitrarily!)")
    end
end
function convertfunction!(el, args::Tuple)
    el = convertfunction(el, args)
end

function rouchon((t0, tf), ρ, H, J, C; fn=ρ->ρ, dt=1e-4, records=Record[])

    ### initialization and conversion
    ts = range(t0, tf, step=dt)                         # time series
    J = map(J) do (j, Γ)
        jj = applicable(j, tt) ? j : (t -> j)
        ΓΓ = applicable(Γ, tt) ? Γ : (t -> Γ)
        (jj, ΓΓ)
    end
    C = map(C) do (m, Γ, η)
            mm = applicable(m, tt) ? m : (t -> m)
            ΓΓ = applicable(Γ, tt) ? Γ : (t -> Γ)
            (mm, ΓΓ, η)
    end


    ### check if a pure-state simulation will suffice; otherwise, convert initial state to density matrix
    pure = (typeof(ρ) <: Ket) && isempty(J)
    if !pure
        ρ = dm(ρ)
    end

    # loop elements
    ρ0 = ρ
    ρs = [ρ0]
    dys::Function = rouchon_dy(C, dt)
    dWs = map(t -> rand(Normal(0, √dt), length(C)), ts)

    # initalize rouchon
    Id = pure ? identityoperator(ρ.basis) : identityoperator(ρ.basis_l)
    M = rouchon_M(H, J, C, dt; dys=dys)
    D = rouchon_D(J, C, dt)

    # loop over times
    for (t, dW) in zip(ts[2:end], dWs[2:end])
        ρ = pure ? normalize(M(t, ρ, dW) * ρ) : normalize(M(t, ρ, dW) * ρ * M(t, ρ, dW)' + D(t, ρ))
        push!(ρs, ρ)
    end
    dy_list::Readouts = map(zip(ts, ρs, dWs)) do (t, ρ, dW)
                   dys(t, ρ, dW)
               end
    records::Records = dys_to_records(trans(dy_list)::Records, C, dt, ts)

	# revert to r (normalized) version of record
    # records = sim ? dys_to_records(dys, C0, dt) : records
    # for (i, (c, Γ, η)) in enumerate(C0)
    #     τm = 1/(2Γ*η)
    #     push!(records, dy[i] .* sqrt(τm) ./ dt)
    # end

    return Solution(collect(ts), ρs, records)
end # rouchon


