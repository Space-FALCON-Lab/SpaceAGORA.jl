# Battery, thermal and data constraints for the fixed-reference terminal-gain
# scheduler. The constraint formulation is ported from Constellation-Controllability
# (src/resource_constrained_scheduling.jl, docs/RESOURCE_CONSTRAINED_SCHEDULING.md
# there); every hardware coefficient is read from the parts mounted on the
# spacecraft: `Battery`, `SolarArray`, `PowerLoad`, `LaserTerminal`,
# `RadioTerminal`, `DataStore`, `DataSource`, `ThermalNode`, and one `DataLink`
# (`FixedRateBudget`) beside each laser link. The time-coupled MILP itself,
# `optimize_resources`, lives in the JuMP/HiGHS package extension.

using ..SimulationModel: SolarArray, Battery, PowerLoad, add_power!, battery, initial_battery_wh, bus_load_w,
    solar_flux_w_m2, AU_M, DataStore, DataSource, add_data!, data_store, ThermalNode, add_thermal!, thermal_node,
    LaserTerminal, RadioTerminal, add_actuator!, mounted_actuators

"""Ground stations (latitude, longitude in degrees) used by `resource_environment`."""
const RESOURCE_STATIONS = [(0.0, 0.0), (35.0, 120.0), (-35.0, -120.0)]

"""
    walker27_resource_constellation(config=LaserSchedulingConfig(); laser_efficiency=0.35,
        laser_data_rate_bps=1e7)

`walker27_constellation` with illustrative resource hardware mounted on every root
link: a Sun-tracking `SolarArray` (1.5 kW), a 3 kWh `Battery` at 80 % charge, a
300 W bus `PowerLoad`, a `LaserTerminal` of the configured power, a `DataStore`
and a `ThermalNode`. Helpers carry a 3 kbit/s `DataSource`; the target carries a
50 W laser-receiver `PowerLoad` and a `RadioTerminal` for its ground downlink
(150 W draw, 200 Mbit/s). Laser links stay on the virtual terminals, so the force
model is unchanged. `data_interlinks` holds one `DataLink` per laser edge
(`FixedRateBudget(laser_data_rate_bps)`), keyed by `data_keys` in edge order: an
active laser can also carry data from helper to target.
"""
function walker27_resource_constellation(config=LaserSchedulingConfig();laser_efficiency=.35,laser_data_rate_bps=1e7)
    case=walker27_constellation(config)
    for (id,sc) in enumerate(case.spacecraft)
        root=sc.root
        add_power!(sc,root,SolarArray(normal=[1.,0,0],area_m2=5.,efficiency=.25,max_power_w=1500.))
        add_power!(sc,root,Battery(capacity_wh=3000.,initial_soc=.8,max_charge_w=35e3,max_discharge_w=35e3))
        add_power!(sc,root,PowerLoad(power_w=300.,name="bus"))
        add_actuator!(sc,root,LaserTerminal(boresight=[1.,0,0],max_power_w=config.power_W,
            magnification=config.reflection_coefficient,range_m=config.max_range_m,efficiency=laser_efficiency))
        add_thermal!(sc,root,ThermalNode(heat_capacity_j_k=150e3,conductance_w_k=150.,initial_k=293.,
            minimum_k=270.,maximum_k=310.,optical_absorptivity=.05,optical_reflectivity=.95))
        if id==1
            add_power!(sc,root,PowerLoad(power_w=50.,name="laser receiver electronics"))
            add_actuator!(sc,root,RadioTerminal(boresight=[0.,0,1],max_power_w=60.,gain_dbi=10.,efficiency=.4,max_rate_bps=2e8))
            add_data!(sc,root,DataStore(capacity_bits=15e9,initial_bits=0.))
        else
            add_data!(sc,root,DataStore(capacity_bits=1e9,initial_bits=.6e9))
            add_data!(sc,root,DataSource(rate_bps=3e3,name="payload"))
        end
    end
    data_interlinks=InterlinkModel()
    data_link=DataLink(transmit_power_w=0.,budget=FixedRateBudget(rate_bps=laser_data_rate_bps),range_m=config.max_range_m)
    data_keys=[register_link!(data_interlinks,case.spacecraft,(s,1),(t,1),data_link) for (s,t) in case.edges]
    return merge(case,(;data_interlinks,data_keys))
end

_only_part(sc,T,what)=begin
    parts=[entry.actuator for entry in mounted_actuators(sc) if entry.actuator isa T]
    length(parts)==1 || throw(ArgumentError("Spacecraft $(sc.id) needs exactly one $what, has $(length(parts))."))
    only(parts)
end
function _required(x,what,sc)
    x===nothing && throw(ArgumentError("Spacecraft $(sc.id) needs a $what for resource scheduling."))
    x
end

"""
    resource_model(case; reserve_soc=0.2, terminal_bits=..., return_capture_fraction=1.0)

The scheduling coefficients, read from the parts mounted on `case.spacecraft`
(spacecraft 1 is the target; see `walker27_resource_constellation`), in the units
the MILP uses (kW, kWh, kJ/K, kW/K, K, Gbit, Gbit/s):

- battery: `Battery` capacity, initial energy, charge/discharge efficiencies and
  limits; `reserve_soc` (a planning choice, not hardware) is the lowest planned
  state of charge.
- solar: Sun-tracking `SolarArray`s, `min(max_power_w, area * efficiency * flux)`
  at 1 AU, scaled by the sampled sunlight; constant `PowerLoad`s are the bus.
- laser: the emitter's `LaserTerminal` draws `max_power_w / efficiency`; its waste
  heat is that draw less the optical power, plus the return the emitter absorbs
  (`return_capture_fraction` times the target's reflectivity times the emitter's
  absorptivity); the target absorbs its absorptivity times the optical power.
- radio: the target's `RadioTerminal` draws `max_power_w / efficiency` and
  downlinks at `max_rate_bps`.
- data: `DataStore` capacity and initial bits, constant `DataSource` production,
  and each edge's `DataLink` `FixedRateBudget` rate; `terminal_bits` (a delivery
  requirement, not hardware) caps each store's data at the end.
- thermal: each `ThermalNode`; all electrical draw except the laser's optical
  output becomes internal heat.
"""
function resource_model(case;reserve_soc=.2,
        terminal_bits=vcat(.1e9,fill(.15e9,length(case.spacecraft)-1)),return_capture_fraction=1.)
    sc=case.spacecraft;ns=length(sc);edges=case.edges
    0<=reserve_soc<=1 || throw(ArgumentError("reserve_soc must lie in [0, 1]."))
    0<=return_capture_fraction<=1 || throw(ArgumentError("return_capture_fraction must lie in [0, 1]."))
    length(terminal_bits)==ns || throw(DimensionMismatch("terminal_bits needs one entry per spacecraft."))
    all(t==1 && s!=1 for (s,t) in edges) || throw(ArgumentError("Single-target star topology requires target ID 1."))
    bat=[_required(battery(m),"Battery",m) for m in sc]
    thermal=[_required(thermal_node(m),"ThermalNode",m) for m in sc]
    store=[_required(data_store(m),"DataStore",m) for m in sc]
    flux=solar_flux_w_m2(AU_M)
    solar=[sum((min(a.max_power_w,a.area_m2*a.efficiency*flux) for l in m.links for a in l.solar_arrays);init=0.) for m in sc]
    production=map(sc) do m
        sources=[s for l in m.links for s in l.data_sources]
        all(s->s isa DataSource,sources) || throw(ArgumentError("Spacecraft $(m.id): only constant DataSource production is supported."))
        sum((s.rate_bps for s in sources);init=0.)
    end
    lasers=[_only_part(sc[s],LaserTerminal,"LaserTerminal") for (s,_) in edges]
    all(l->l.max_power_w==case.config.power_W,lasers) || throw(ArgumentError(
        "Each emitter's LaserTerminal max_power_w must equal the laser force power, $(case.config.power_W) W."))
    radio=_only_part(sc[1],RadioTerminal,"RadioTerminal")
    graph=case.data_interlinks.graph
    links=[graph.edge_properties[graph.edge_index_lookup[key]] for key in case.data_keys]
    all(l->l isa DataLink && l.budget isa FixedRateBudget,links) || throw(ArgumentError(
        "Each edge needs a DataLink with a FixedRateBudget in case.data_interlinks."))
    target=thermal[1]
    # Optical power hitting the target and the return reaching each emitter, kW.
    optical=[l.max_power_w/1e3 for l in lasers]
    laser_heat_emitter=[optical[e]/lasers[e].efficiency-optical[e]+
        thermal[s].optical_absorptivity*target.optical_reflectivity*return_capture_fraction*optical[e]
        for (e,(s,_)) in enumerate(edges)]
    return (capacity_kWh=[b.capacity_wh/1e3 for b in bat],initial_kWh=[initial_battery_wh(m)/1e3 for m in sc],
        reserve_kWh=[reserve_soc*b.capacity_wh/1e3 for b in bat],
        charge_efficiency=[b.charge_efficiency for b in bat],discharge_efficiency=[b.discharge_efficiency for b in bat],
        max_charge_kW=[b.max_charge_w/1e3 for b in bat],max_discharge_kW=[b.max_discharge_w/1e3 for b in bat],
        solar_kW=solar./1e3,bus_kW=[bus_load_w(m)/1e3 for m in sc],
        laser_kW=optical./[l.efficiency for l in lasers],laser_heat_emitter_kW=laser_heat_emitter,
        laser_heat_receiver_kW=target.optical_absorptivity.*optical,
        radio_kW=radio.max_power_w/radio.efficiency/1e3,downlink_Gbps=radio.max_rate_bps/1e9,
        heat_capacity_kJ_K=[t.heat_capacity_j_k/1e3 for t in thermal],conductance_kW_K=[t.conductance_w_k/1e3 for t in thermal],
        initial_K=[t.initial_k for t in thermal],minimum_K=[t.minimum_k for t in thermal],maximum_K=[t.maximum_k for t in thermal],
        buffer_Gbit=[d.capacity_bits/1e9 for d in store],initial_Gbit=[d.initial_bits/1e9 for d in store],
        generation_Gbps=production./1e9,terminal_Gbit=terminal_bits./1e9,
        link_Gbps=[l.budget.rate_bps/1e9 for l in links])
end

"""
    resource_environment(case, times, reference; stations=RESOURCE_STATIONS)

Conservative per-interval coefficients from 21 samples of `reference` per interval:
`sun[i,k]` is 1 only if spacecraft i is lit throughout (cylindrical shadow, Sun on +x),
`sink_K[i,k]` is 290 K if lit at any sample and 280 K otherwise, and `ground[k]` is 1
only if the target sees a station above 10 degrees elevation throughout.
"""
function resource_environment(case,times,reference;stations=RESOURCE_STATIONS)
    ns=length(case.spacecraft);n=length(times)-1;Re=case.config.planet.Re
    sun=ones(Int,ns,n);hot=zeros(Int,ns,n);ground=ones(Int,n)
    for k in 1:n,t in range(times[k],times[k+1];length=21)
        state=reference(t)
        for i in 1:ns
            r=ks_state_to_cartesian(view(state,slice(i))).position_ii_m
            lit=!(r[1]<0 && hypot(r[2],r[3])<Re)
            sun[i,k]=min(sun[i,k],lit);hot[i,k]=max(hot[i,k],lit)
            if i==1
                visible=false
                for (lat,lon) in stations
                    phi=deg2rad(lat);theta=deg2rad(lon)+7.2921159e-5*t
                    normal=[cos(phi)*cos(theta),cos(phi)*sin(theta),sin(phi)]
                    los=r-Re*normal
                    visible |= dot(los,normal)/norm(los)>=sin(deg2rad(10.))
                end
                ground[k]=min(ground[k],visible)
            end
        end
    end
    (;sun,sink_K=280 .+10hot,ground)
end

# The battery's energy rate (kW into storage) for a net bus power (collected solar
# minus load), as Battery charges and discharges it.
_stored_kW(net,p,i)=net>=0 ? net*p.charge_efficiency[i] : net/p.discharge_efficiency[i]

"""Closed-form resource replay; infeasible states are reported, not clipped.
Without `solar`, collect as much sunlight as the battery can absorb; with it,
cap the scheduled collection at available sunlight."""
function resource_replay(u,v,flow,downlink,edges,times,env,p;solar=nothing)
    ns=length(p.initial_kWh);n=length(times)-1
    E=zeros(ns,n+1);T=copy(E);q=copy(E);used=zeros(ns,n)
    E[:,1]=p.initial_kWh;T[:,1]=p.initial_K;q[:,1]=p.initial_Gbit
    for k in 1:n
        dt=times[k+1]-times[k];load=copy(p.bus_kW);heat=copy(p.bus_kW)
        net=copy(p.generation_Gbps)
        for (e,(s,t)) in enumerate(edges)
            load[s]+=p.laser_kW[e]*u[e,k]
            heat[s]+=p.laser_heat_emitter_kW[e]*u[e,k]
            heat[t]+=p.laser_heat_receiver_kW[e]*u[e,k]
            net[s]-=flow[e,k];net[t]+=flow[e,k]
        end
        load[1]+=p.radio_kW*v[k];heat[1]+=p.radio_kW*v[k];net[1]-=downlink[k]
        for i in 1:ns
            available=p.solar_kW[i]*env.sun[i,k]
            room=(p.capacity_kWh[i]-E[i,k])*3600/dt/p.charge_efficiency[i]
            used[i,k]=solar===nothing ? min(available,max(0.,load[i]+room)) : min(available,solar[i,k])
            E[i,k+1]=E[i,k]+dt/3600*_stored_kW(used[i,k]-load[i],p,i)
            alpha=exp(-p.conductance_kW_K[i]*dt/p.heat_capacity_kJ_K[i])
            T[i,k+1]=alpha*T[i,k]+(1-alpha)*(env.sink_K[i,k]+heat[i]/p.conductance_kW_K[i])
            q[i,k+1]=q[i,k]+dt*net[i]
        end
    end
    (;energy=E,temperature=T,queue=q,solar=used)
end

"""Conservative service warm start, returned only if every enabled resource bound holds."""
function resource_seed(available,edges,times,env,p;battery,thermal,data,initial_u=nothing)
    ne,n=size(available);dt=diff(times)
    u=zeros(Int,ne,n);v=zeros(Int,n);f=zeros(ne,n);d=zeros(n)
    if initial_u!==nothing
        data && throw(ArgumentError("A supplied laser-only seed requires data=false."))
        size(initial_u)==size(available) && all(x->x==0 || x==1,initial_u) &&
            all(initial_u.<=available) && all(sum(initial_u;dims=1).<=1) ||
            throw(ArgumentError("Seed must be binary, available, and use at most one link per block."))
        u.=initial_u
    end
    q=copy(p.initial_Gbit)
    last_ground=findlast(!iszero,env.ground)
    data && last_ground===nothing && return nothing
    if data
        for k in 1:n
            if env.ground[k]!=0 && q[1]>1e-12
                v[k]=1;d[k]=min(p.downlink_Gbps,q[1]/dt[k])
            end
            if k<last_ground
                for (e,(s,t)) in enumerate(edges)
                    needed=k>=last_ground-100 ?
                        q[s]+p.generation_Gbps[s]*(last(times)-times[k])>p.terminal_Gbit[s]+1e-10 :
                        q[s]>.5p.buffer_Gbit[s]
                    if available[e,k] && needed
                        rate=min(p.link_Gbps[e],q[s]/dt[k],max(0.,p.buffer_Gbit[t]-q[t])/dt[k])
                        if rate>1e-12
                            u[e,k]=1;f[e,k]=rate
                            break
                        end
                    end
                end
            end
            q .+=dt[k].*p.generation_Gbps
            q[1]-=dt[k]*d[k]
            for (e,(s,t)) in enumerate(edges)
                q[s]-=dt[k]*f[e,k];q[t]+=dt[k]*f[e,k]
            end
        end
    end
    states=resource_replay(u,v,f,d,edges,times,env,p)
    if battery
        all(p.reserve_kWh.-1e-7 .<= states.energy .<= p.capacity_kWh.+1e-7) || return nothing
        net=states.solar.-_loads(u,v,edges,p)
        all(-p.max_discharge_kW.-1e-7 .<= net .<= p.max_charge_kW.+1e-7) || return nothing
    end
    thermal && !all(p.minimum_K.-1e-7 .<= states.temperature .<= p.maximum_K.+1e-7) && return nothing
    if data
        minimum(states.queue)>=-1e-7 && all(states.queue.<=p.buffer_Gbit.+1e-7) &&
            all(states.queue[:,end].<=p.terminal_Gbit.+1e-7) || return nothing
    end
    return (;u,v,flow=f,downlink=d,states)
end

# Electrical load (kW) per spacecraft and interval for a laser and radio schedule.
function _loads(u,v,edges,p)
    load=repeat(p.bus_kW,1,size(u,2))
    for (e,(s,_)) in enumerate(edges)
        load[s,:].+=p.laser_kW[e].*u[e,:]
    end
    load[1,:].+=p.radio_kW.*v
    load
end

"""
    optimize_resources(scores, available, edges, times, env, p;
        battery=true, thermal=true, data=true, time_limit=120., initial_u=nothing)

Signed terminal-gain MILP with optional time-coupled battery, lumped-thermal and
data-queue constraints (HiGHS); `p` comes from `resource_model`. All edges must
end at spacecraft 1, whose `RadioTerminal` downlinks its queue. Requires
`using JuMP, HiGHS`.

The battery splits the net bus power into charging and discharging. With unit
efficiencies its energy is exact; below them the program may charge and
discharge in the same interval, which wastes energy but never creates it, so the
battery (`resource_replay` without `solar`) holds at least the planned energy and
the reserve still holds. The returned `energy` is then a lower bound.
"""
function optimize_resources end

"""
    plan_resource_schedule(case, fine_times; objective=:a, stride=5, battery=true,
        thermal=true, data=true, parameters=resource_model(case), reference=nothing,
        time_limit=120., linearization=:ks)

`plan_laser_schedule` terminal gains on `fine_times`, summed over command blocks of
`stride` fine intervals (a block is available only if every fine interval is), then
`optimize_resources` on those blocks with the hardware of `case`
(`walker27_resource_constellation`). `command` is the block schedule expanded to
`fine_times`; `predicted` is the linear terminal change (a in metres, angles in
radians). Requires `using JuMP, HiGHS`.
"""
function plan_resource_schedule(case,fine;objective=:a,stride=5,battery=true,thermal=true,data=true,
        parameters=nothing,reference=nothing,time_limit=120.,linearization=:ks)
    isempty(methods(optimize_resources)) && error("Resource scheduling needs `using JuMP, HiGHS`.")
    nf=length(fine)-1
    nf%stride==0 || throw(ArgumentError("The fine interval count must be a multiple of stride."))
    parameters===nothing && (parameters=resource_model(case))
    plan=plan_laser_schedule(case,fine;objective,reference,linearization)
    n=nf÷stride;ne=length(case.edges);times=fine[1:stride:end]
    gains=zeros(5,ne,n);available=trues(ne,n)
    for b in 1:n
        ids=stride*(b-1)+1:stride*b
        gains[:,:,b]=dropdims(sum(plan.gains[:,:,ids];dims=3);dims=3)
        available[:,b]=vec(all(plan.available[:,ids];dims=2))
    end
    scores=dropdims(sum(gains.*reshape(plan.weights,5,1,1);dims=1);dims=1)
    env=resource_environment(case,times,plan.reference)
    result=optimize_resources(scores,available,case.edges,times,env,parameters;battery,thermal,data,time_limit)
    context=(;times,fine_times=collect(fine),gains,available,scores,weights=plan.weights,env,
        reference=plan.reference,parameters)
    result.feasible || return merge(result,context)
    return merge(result,context,(;command=repeat(result.u;inner=(1,stride)),
        predicted=vec(sum(gains.*reshape(result.u,1,ne,n);dims=(2,3)))))
end
