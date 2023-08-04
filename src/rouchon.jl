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
    Jsum = mapreduce(+, J; init=0Id) do (j, Γ)        Γ * j' * j end
    msum = mapreduce(+, C; init=0Id) do (m, Γ, η)     η * Γ/2 * m' * m  end

    # stochastic part
    mdy(ρ, dWs) = mapreduce(+, zip(C, dys(ρ, dWs)); init=0Id) do ((m, Γ, η), dy)
                        sqrt(η * Γ/2) * m * dy
                  end
    mm(dWs) = mapreduce(+, 1:length(C), 1:length(C); init=0Id) do r, s
            (mr, Γr, ηr) = C[r]
            (ms, Γs, ηs) = C[s]
            sqrt(ηr * ηs * Γr * Γs) / 4 * mr * ms * (dWs[r] * dWs[s] - δ(r,s) * dt)
    end

    # total   
    Mn(t, ρ, dWs) = Id - (im * H(t) + Jsum/2 + msum/2)*dt + mdy(ρ, dWs) + mm(dWs)
    return Mn
end

function rouchon_D(J)
    Id(ρ) = identityoperator(ρ.basis_l)
    return ρ -> mapreduce(+, J; init=0Id(ρ)) do (j, Γ)        Γ * j * ρ * j'           end
end

function rouchon_dy(C, dt)
   return (ρ, dWs) -> map(zip(C, dWs)) do ((m, Γ, η), dW)
                            c = sqrt(Γ * η/2) * m
                            real(expect(c + c', ρ)) * dt + dW
                end
end

function record_to_dy(record::Record, (c, Γ, η), dt, ts)
    τm = 1/(2Γ*η)
    return record .* dt ./ sqrt(τm)
end
function dy_to_record(dys, (m, Γ, η), dt)::Record
    τm = 1/(2Γ*η)
    return dys .* sqrt(τm) ./ dt
end
function records_to_dys(records::Vector{Record}, C, dt)
    map(zip(records, C)) do (record, ctup)
        record_to_dy(record, ctup, dt)
    end
end
function dys_to_records(dys, C, dt)::Vector{Record}
    map(zip(dys, C)) do (dy, ctup)
        dy_to_record(dy, ctup, dt)
    end
end

function rouchon((t0, tf), ρ, H, J, C; fn=ρ->ρ, dt=1e-4, records=Record[])

    ### initialization and conversion
    ts = range(t0, tf, step=dt)                         # time series

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
    M = rouchon_M(H, J, C, dt; dys=dys)
    D = rouchon_D(J)

    # loop over times
    for (t, dW) in zip(ts[2:end], dWs[2:end])
        ρ = pure ? normalize(M(t, ρ, dW) * ρ) : normalize(M(t, ρ, dW) * ρ * M(t, ρ, dW)' + D(ρ)*dt)
        push!(ρs, ρ)
    end
    dy_list::Readouts = map(zip(ρs, dWs)) do (ρ, dW)
                   dys(ρ, dW)
               end
    records::Records = dys_to_records(trans(dy_list)::Records, C, dt)

	# revert to r (normalized) version of record
    # records = sim ? dys_to_records(dys, C0, dt) : records
    # for (i, (c, Γ, η)) in enumerate(C0)
    #     τm = 1/(2Γ*η)
    #     push!(records, dy[i] .* sqrt(τm) ./ dt)
    # end

    return Solution(collect(ts), ρs, records)
end # rouchon


