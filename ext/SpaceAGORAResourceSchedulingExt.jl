module SpaceAGORAResourceSchedulingExt

# Time-coupled battery/thermal/data MILP for the constellation terminal-gain
# scheduler, ported from Constellation-Controllability's
# src/resource_constrained_scheduling.jl. Its coefficients come from
# `ConstellationScheduling.resource_model`, which reads the mounted hardware;
# replay and warm start live there too.
using SpaceAGORA
using JuMP, HiGHS
const CS = SpaceAGORA.ConstellationScheduling

function CS.optimize_resources(scores, available, edges, times, env, p;
        battery=true, thermal=true, data=true, time_limit=120.,initial_u=nothing)
    ne,n=size(scores);ns=length(p.capacity_kWh);dt=diff(times)
    @assert length(dt)==n && size(available)==size(scores)
    @assert size(env.sun)==(ns,n) && size(env.sink_K)==(ns,n)
    @assert all(t==1 && s!=1 for (s,t) in edges)
    outgoing=[[e for (e,(s,_)) in enumerate(edges) if s==i] for i in 1:ns]
    incoming=[[e for (e,(_,t)) in enumerate(edges) if t==i] for i in 1:ns]
    model=Model(HiGHS.Optimizer);set_silent(model)
    set_optimizer_attribute(model,"threads",1)
    set_optimizer_attribute(model,"mip_rel_gap",1e-4)
    set_time_limit_sec(model,time_limit)
    @variable(model,u[1:ne,1:n],Bin)
    @constraint(model,[e=1:ne,k=1:n],u[e,k]<=available[e,k])
    @constraint(model,[k=1:n],sum(u[:,k])<=1)
    @variable(model,v[1:n],Bin)
    @variable(model,f[1:ne,1:n]>=0)
    @variable(model,d[1:n]>=0)
    if data
        @variable(model,q[1:ns,1:n+1]>=0)
        @constraint(model,[i=1:ns,k=1:n+1],q[i,k]<=p.buffer_Gbit[i])
        @constraint(model,[i=1:ns],q[i,1]==p.initial_Gbit[i])
        @constraint(model,[i=1:ns],q[i,n+1]<=p.terminal_Gbit[i])
        @constraint(model,[k=1:n],v[k]<=env.ground[k])
        @constraint(model,[k=1:n],d[k]<=p.downlink_Gbps*v[k])
        @constraint(model,[e=1:ne,k=1:n],f[e,k]<=p.link_Gbps[e]*u[e,k])
        for i in 1:ns,k in 1:n
            out=sum((f[e,k] for e in outgoing[i]);init=0.)+(i==1 ? d[k] : 0.)
            inflow=sum((f[e,k] for e in incoming[i]);init=0.)
            @constraint(model,dt[k]*out<=q[i,k])
            @constraint(model,q[i,k+1]==q[i,k]+dt[k]*(p.generation_Gbps[i]+inflow-out))
        end
    else
        fix.(v,0.;force=true);fix.(f,0.;force=true);fix.(d,0.;force=true)
    end
    @expression(model,load[i=1:ns,k=1:n],p.bus_kW[i]+
        sum(p.laser_kW[e]*u[e,k] for e in outgoing[i])+(i==1 ? p.radio_kW*v[k] : 0.))
    @expression(model,heat[i=1:ns,k=1:n],p.bus_kW[i]+
        sum(p.laser_heat_emitter_kW[e]*u[e,k] for e in outgoing[i])+
        sum(p.laser_heat_receiver_kW[e]*u[e,k] for e in incoming[i])+
        (i==1 ? p.radio_kW*v[k] : 0.))
    if battery
        @variable(model,E[1:ns,1:n+1])
        @variable(model,solar[1:ns,1:n]>=0)
        # Battery: a surplus charges at charge_efficiency, a shortfall discharges
        # over discharge_efficiency, each within its power limit.
        @variable(model,charge[1:ns,1:n]>=0)
        @variable(model,discharge[1:ns,1:n]>=0)
        @constraint(model,[i=1:ns,k=1:n+1],E[i,k]>=p.reserve_kWh[i])
        @constraint(model,[i=1:ns,k=1:n+1],E[i,k]<=p.capacity_kWh[i])
        @constraint(model,[i=1:ns],E[i,1]==p.initial_kWh[i])
        @constraint(model,[i=1:ns,k=1:n],solar[i,k]<=p.solar_kW[i]*env.sun[i,k])
        @constraint(model,[i=1:ns,k=1:n],solar[i,k]-load[i,k]==charge[i,k]-discharge[i,k])
        @constraint(model,[i=1:ns,k=1:n],charge[i,k]<=p.max_charge_kW[i])
        @constraint(model,[i=1:ns,k=1:n],discharge[i,k]<=p.max_discharge_kW[i])
        @constraint(model,[i=1:ns,k=1:n],E[i,k+1]==E[i,k]+dt[k]/3600*
            (p.charge_efficiency[i]*charge[i,k]-discharge[i,k]/p.discharge_efficiency[i]))
    end
    if thermal
        @variable(model,T[1:ns,1:n+1])
        @constraint(model,[i=1:ns,k=1:n+1],T[i,k]>=p.minimum_K[i])
        @constraint(model,[i=1:ns,k=1:n+1],T[i,k]<=p.maximum_K[i])
        @constraint(model,[i=1:ns],T[i,1]==p.initial_K[i])
        for i in 1:ns,k in 1:n
            alpha=exp(-p.conductance_kW_K[i]*dt[k]/p.heat_capacity_kJ_K[i])
            @constraint(model,T[i,k+1]==alpha*T[i,k]+(1-alpha)*
                (env.sink_K[i,k]+heat[i,k]/p.conductance_kW_K[i]))
        end
    end
    scale=max(maximum(abs,scores),eps())
    @objective(model,Max,sum(scores.*u)/scale)
    seed=CS.resource_seed(available,edges,times,env,p;battery,thermal,data,initial_u)
    initial_u!==nothing && seed===nothing && throw(ArgumentError("Supplied seed violates enabled resource bounds."))
    if seed!==nothing
        set_start_value.(u,seed.u);set_start_value.(v,seed.v)
        set_start_value.(f,seed.flow);set_start_value.(d,seed.downlink)
        data && set_start_value.(q,seed.states.queue)
        if battery
            set_start_value.(E,seed.states.energy);set_start_value.(solar,seed.states.solar)
            net=seed.states.solar.-CS._loads(seed.u,seed.v,edges,p)
            set_start_value.(charge,max.(net,0.));set_start_value.(discharge,max.(.-net,0.))
        end
        thermal && set_start_value.(T,seed.states.temperature)
    end
    optimize!(model)
    status=string(termination_status(model))
    has_values(model) || return (;status,feasible=false)
    isapprox(maximum(abs,value.(u)-round.(value.(u))),0.;atol=1e-5) || error("Nonbinary incumbent")
    return (status=status,feasible=true,u=round.(Int,value.(u)),v=round.(Int,value.(v)),
        flow=value.(f),downlink=value.(d),energy=battery ? value.(E) : nothing,
        solar=battery ? value.(solar) : nothing,temperature=thermal ? value.(T) : nothing,
        queue=data ? value.(q) : nothing,objective=objective_value(model)*scale,
        bound=objective_bound(model)*scale,gap=relative_gap(model),seconds=solve_time(model))
end

end
