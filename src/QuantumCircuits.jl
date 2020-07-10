module QuantumCircuits

using Reexport
@reexport using QuantumOpticsBase
using Distributions

δ(i,j) = Int(i == j);

"""
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
            dy_c[t] indexing the mth trajectory at time T[t]
fn : ρ → Any :: eval function (e.g. to return expectation values instead of density matrices)
"""

function rouchon(T, ρ, H, J, C; fn=ρ->ρ, dt=1e-4, dy=[])
    T = range(T[1], T[2], step=dt)
    Id = identityoperator(ρ.basis_l)
    Nt = length(T)
    Nj = length(J)
    Nc = length(C)
    
    ρs = Any[]
    
    # randomly sample noise time series for EACH stochastic collapse operator C
    if length(dy) == 0
        sim = true
        dist = Normal(0, sqrt(dt))
        dW=[]
        for c in 1:Nc
            push!(dW,rand(dist,Nt))
            push!(dy,zeros(Nt))
        end
    else
        sim = false
    end


    for t in 1:Nt
        M = Id - im*H(T[t])*dt

        # iterate over deterministic collapse operators J
        D = 0Id
        for j in 1:Nj
          M += -0.5J[j]'*J[j]*dt
          D += J[j]*ρ*J[j]'*dt
        end
            
        # initialize dy value, if simulation
        if sim==true
            for c in 1:Nc
                dy[c][t] = real(tr(C[c]*ρ + ρ*C[c]')*dt) + dW[c][t]
            end
        end

        # iterate over stochastic collapse operators C
        for c in 1:Nc
                
            M += C[c]*dy[c][t]
            D += -C[c]*ρ*C[c]'*dt
                
            # nested sum    
            for s in 1:Nc
                M += 0.5C[c]*C[s]*(dy[c][t]*dy[s][t] - δ(c,s)*dt)
            end
                
        end


        # update ρ according to Rouchon method
        ρ = M*ρ*M' + D
        ρ = ρ/tr(ρ)
        push!(ρs, fn(ρ))
            
    end
        
    return (T,ρs,dy)
        
end

"""
    ensemble(solve, T, ρ0, H, J, C; <keyword arguments>)

solve :: solver function
T :: tuple (ti,tf)
ρ0 :: initial density matrix
H :: time-dependent Hamiltonian
J :: deterministic collapse operators
C :: stochastic collapse operators

Keyword Arguments:

dt :: time step (default dt=1e-4)
record :: (default record=[]), i.e. simulation generates record time series. 
            record should be input in the shape [rec_1,...,rec_Nc] given 
            length(C)=Nc collapse operators. Records should have shape
            rec_c[m,t] indexing the mth trajectory at time T[t]
N :: number of trajectories (default N=10)
"""

function ensemble(solve, T, ρ0, H, J, C; fns=[ρ->ρ], dt=1e-4, record=[], N=10)
    ti = T[1]; tf = T[2];
    T = range(ti, tf, step=dt)
    Nfns = length(fns)
    Nt = length(T)
    Nc = length(C)
    
    evals = Array{Any, 3}(undef, Nfns, N, Nt)
    
    # if input record is empty, "record" will store readout after it is generated in simulation
    # initialize with zeros
    if length(record) == 0
        for c in 1:Nc
            push!(record,zeros(N,Nt))
        end
        sim = true
    else
        sim = false
    end
    
    
    for m in 1:N
        dy = []
        # if input record is NOT empty, dy will take mth trajectory's record as input
        if sim==false
            for c in 1:Nc
                push!(dy,record[c][m,:])
            end
        end
        T,ρs,dy = solve((ti,tf), ρ0, H, J, C; dt=dt, dy=dy)

        for t in 1:Nt
            # readout record
            for c in 1:Nc
                record[c][m,t] = dy[c][t]
            end

            # fn evaluations
            if Nfns > 0
                for nfn in 1:Nfns
                    evals[nfn,m,t] = fns[nfn](ρs[t])
                end
            else 

            end
        end
    end
    
    return (T,evals,record)
    
end

export δ, rouchon, ensemble

end # module
