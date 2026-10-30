module ConstellationScheduling
using LinearAlgebra, StaticArrays, OrdinaryDiffEq
using ..SimulationModel.DynamicsKS
using ..SimulationModel: SpacecraftModel, Link, CartesianInitialCondition
using ..SimulationModel.InterLinkModels
using ..SimulationModel.InterLinkModels: link_endpoint
using ..SimulationModel.SchedulingPolicies
import ..SimulationModel.SchedulingPolicies: candidate_score
using ..SimulationModel.FrameTransforms: orbitalelemtorv, rvtoorbitalelement, rtn_dcm_from_inertial
using ..OrbitalElementSensitivity: orbital_element_jacobian

export LaserSchedulingConfig, paper_walker_constellation, plan_laser_schedule,
    replay_laser_schedule, laser_element_history, laser_geometry, signed_link_schedule,
    TerminalGainPolicy

"""Fixed-interval terminal-gain policy for SpaceAGORA's native interlink matcher.

Scores are terminal contributions, not instantaneous rates. A plan's policy
contains only its selected activations so replay does not re-optimize or select
a substitute when a commanded link becomes unavailable.
"""
struct TerminalGainPolicy <: AbstractSchedulingPolicy
    times::Vector{Float64}
    keys::Vector{InterLinkKey}
    scores::Matrix{Float64}
    function TerminalGainPolicy(times, keys, scores)
        length(times)>=2 && all(isfinite,times) && all(diff(times).>0) ||
            throw(ArgumentError("Policy times must be finite and strictly increasing."))
        size(scores)==(length(keys),length(times)-1) || throw(DimensionMismatch("Policy scores/time grid mismatch."))
        all(isfinite,scores) && allunique(keys) || throw(ArgumentError("Invalid scores or duplicate link keys."))
        new(collect(Float64,times),collect(InterLinkKey,keys),Matrix{Float64}(scores))
    end
end

function candidate_score(policy::TerminalGainPolicy, candidate::LinkCandidate, ctx::SchedulingContext, t::Float64)
    isfinite(t) || throw(ArgumentError("A terminal-gain policy needs physical time in LinkContext."))
    first(policy.times)<=t<last(policy.times) || return 0.0
    e=findfirst(==(candidate.key),policy.keys)
    e===nothing && return 0.0
    return policy.scores[e,searchsortedlast(policy.times,t)]
end

Base.@kwdef struct LaserSchedulingConfig
    planet::KSPropagationParams=KSPropagationParams(Re=6378136.3,μ=3.986004418e14,J2=1.08262668e-3)
    power_W::Float64=10000.
    reflection_coefficient::Float64=100.
    mass_kg::Float64=300.
    max_range_m::Float64=200000.
    scheduling_range_m::Float64=199000.
    atmosphere_height_m::Float64=100000.
    tolerance::Float64=2e-12
end

function paper_walker_constellation(config=LaserSchedulingConfig())
    slots=vcat([(0,0,0)],[(l,p,s) for l in -1:1 for p in -1:1 for s in -1:1 if (l,p,s)!=(0,0,0)])
    x=zeros(6,length(slots));elements=zeros(length(slots),10)
    for (id,(l,p,s)) in enumerate(slots)
        a=config.planet.Re+1e6+5000l;e=.01;i=deg2rad(55);node=deg2rad(40+.6p);omega=deg2rad(30)
        M=deg2rad(.6s+.2p);E=M
        for _ in 1:8
            E-=(E-e*sin(E)-M)/(1-e*cos(E))
        end
        nu=2atan(sqrt(1+e)*sin(E/2),sqrt(1-e)*cos(E/2))
        r,v=orbitalelemtorv(SVector(a,e,i,node,omega,nu,0.),config.planet)
        x[:,id].=vcat(r,v)
        elements[id,:].=[id,l,p,s,a,e,i,node,omega,M]
    end
    spacecraft=map(axes(x,2)) do id
        bus=Link(root=true,m=config.mass_kg)
        SpacecraftModel(root=bus,links=[bus],inertia_tensor=bus.inertia,
            initial_condition=CartesianInitialCondition(SVector{3}(x[1:3,id]),SVector{3}(x[4:6,id])),
            id=id,n_terminal=1)
    end
    edges=[(s,1) for s in 2:length(slots)]
    interlinks=InterlinkModel(record_available=true)
    link=LaserBeamLink(power_w=config.power_W,magnification=config.reflection_coefficient,
        range_m=config.max_range_m,occultation=(;min_ray_altitude_m=config.atmosphere_height_m))
    keys=[register_link!(interlinks,spacecraft,(s,1),(t,1),link) for (s,t) in edges]
    return (;initial_cartesian=x,initial_elements=elements,edges,spacecraft,interlinks,keys,
        period_s=2pi*sqrt((config.planet.Re+1e6)^3/config.planet.μ),config)
end

slice(s)=10s-9:10s

function interlink_state(state,case)
    sc=map(1:length(state)÷10) do s
        x=ks_state_to_cartesian(view(state,slice(s)))
        (;pos=x.position_ii_m,vel=x.velocity_ii_m,mass=case.config.mass_kg)
    end
    return (;sc)
end
link_context(case,t)=LinkContext(time_s=t,planet=(;Rp_e=case.config.planet.Re))

"""Physical range plus opaque Earth/atmosphere segment-clearance gate."""
function laser_geometry(state,case;model=case.interlinks,t=state[10])
    native=interlink_state(state,case);ctx=link_context(case,t)
    available=update_availability!(model,case.spacecraft,native;context=ctx)
    return map(eachindex(case.edges)) do e
        source,target=case.edges[e];first,second=case.keys[e]
        a=link_endpoint(native.sc[first[1]],case.spacecraft[first[1]],first[1])
        b=link_endpoint(native.sc[second[1]],case.spacecraft[second[1]],second[1])
        fa,_,fb,_=interlink_wrench(model.graph.edge_properties[e],nothing,nothing,a,b,ctx)
        d=native.sc[target].pos-native.sc[source].pos;distance=norm(d)
        direction=distance>0 ? d/distance : zeros(3)
        (;distance,direction,available=available[e],target_force=first[1]==target ? fa : fb)
    end
end

function constellation_rhs!(du,state,p,t)
    case,times,command,model=p;c=case.config;n=length(state)÷10
    acceleration=zeros(3,n);k=clamp(searchsortedlast(times,t),1,size(command,2))
    if any(!iszero,view(command,:,k))
        native=interlink_state(state,case);ctx=link_context(case,t)
        update_availability!(model,case.spacecraft,native;context=ctx)
        model.graph.active.=(view(command,:,k).==1) .& model.graph.available
        for s in 1:n
            acceleration[:,s].=interlink_force_on_spacecraft(model,native,s;
                spacecraft=case.spacecraft,context=ctx)/c.mass_kg
        end
    else
        fill!(model.graph.active,false)
    end
    for s in 1:n
        rows=slice(s);y=view(state,rows)
        # Reuse SpaceAGORA's h=-epsilon dynamics; lift only the new laser input.
        f=nonlinear_ks_rhs(y,c.planet)+ks_acceleration_input_matrix(y)*acceleration[:,s]
        du[rows].=f/dot(view(y,1:4),view(y,1:4))
    end
end

function initial_ks(case,t0)
    x=case.initial_cartesian
    vcat([Vector(cartesian_to_ks_state(x[1:3,s],x[4:6,s],case.config.planet;elapsed_time_s=t0)) for s in axes(x,2)]...)
end

function replay_laser_schedule(case,times,command)
    length(times)>=2 && all(diff(times).>0) || throw(ArgumentError("Times must increase."))
    size(command)==(length(case.edges),length(times)-1) || throw(DimensionMismatch("Schedule/time grid mismatch."))
    all(u->u==0 || u==1,command) || throw(ArgumentError("Laser commands must be binary."))
    all(sum(command;dims=1).<=1) || throw(ArgumentError("At most one directed link is supported per interval."))
    c=case.config
    c.mass_kg>0 && c.power_W>=0 && c.reflection_coefficient>=0 || throw(ArgumentError("Invalid laser parameters."))
    u0=initial_ks(case,first(times))
    sol=solve(ODEProblem(constellation_rhs!,u0,(first(times),last(times)),(case,times,command,deepcopy(case.interlinks))),Tsit5();
        reltol=c.tolerance,abstol=c.tolerance,dtmax=minimum(diff(times))/4,tstops=times[2:end-1])
    occursin("Success",string(sol.retcode)) || error("Nonlinear KS replay failed: $(sol.retcode)")
    sol
end

function zoh(A,B,dt)
    n=size(A,1);m=size(B,2);Z=zeros(n+m,n+m);Z[1:n,1:n].=A;Z[1:n,n+1:end].=B
    E=exp(dt*Z);return E[1:n,1:n],E[1:n,n+1:end]
end

function rtn_model(y,c)
    physical=ks_physical_tangent(y,c.planet)
    x=ks_state_to_cartesian(y);r=x.position_ii_m;v=x.velocity_ii_m
    M=Matrix(rtn_dcm_from_inertial(r,v));I3=Matrix{Float64}(I,3,3)
    rd=(I3-M[:,1]*M[:,1]')*v/norm(r)
    nd=(I3-M[:,3]*M[:,3]')*cross(r,physical.acceleration)/norm(cross(r,v))
    td=cross(nd,M[:,1])+cross(M[:,3],rd);Md=hcat(rd,td,nd)
    R=[M' zeros(3,3);zeros(3,3) M'];Rd=[Md' zeros(3,3);zeros(3,3) Md']
    return (;A=Rd*R'+R*physical.A_cartesian*R',B=R*physical.B_cartesian,M)
end

function elements(y,c)
    x=ks_state_to_cartesian(y)
    collect(rvtoorbitalelement(x.position_ii_m,x.velocity_ii_m,c.planet)[1:5])
end

function laser_element_history(case,times,sol)
    reduce(hcat,[elements(view(sol(t),slice(1)),case.config) for t in times])
end

"""Exact interval-separable solution of the paper's signed-gain binary MILP."""
function signed_link_schedule(scores,available)
    size(scores)==size(available) || throw(DimensionMismatch("Scores and availability differ."))
    all(isfinite,scores) || throw(ArgumentError("Scores must be finite."))
    u=zeros(Int,size(scores))
    for k in axes(scores,2)
        best=0.;selected=0
        for e in axes(scores,1)
            if available[e,k] && scores[e,k]>best
                best=scores[e,k];selected=e
            end
        end
        selected>0 && (u[selected,k]=1)
    end
    u
end

"""Original all-off fixed-reference paper scheduler, not active-reference MPC.

Units: a in metres, angles in radians. Objectives :a, :i, :omega, :combined.
Joint scales are 1000 m, one degree, one degree; weights are each 1/3.
"""
function plan_laser_schedule(case,times;objective=:a,reference=nothing)
    objective in (:a,:i,:omega,:combined) || throw(ArgumentError("Unknown objective."))
    all(s!=1 && t==1 for (s,t) in case.edges) || throw(ArgumentError("Paper star topology requires target ID 1."))
    c=case.config;n=length(times)-1;ne=length(case.edges)
    reference===nothing && (reference=replay_laser_schedule(case,times,zeros(Int,ne,n)))
    Fs=Matrix{Float64}[];Vs=Matrix{Float64}[];available=trues(ne,n)
    for k in 1:n
        state=reference((times[k]+times[k+1])/2);model=rtn_model(view(state,slice(1)),c)
        pairs=laser_geometry(state,case)
        B=reduce(hcat,[model.B*(p.target_force/c.mass_kg) for p in pairs])
        F,V=zoh(model.A,B,times[k+1]-times[k]);push!(Fs,F);push!(Vs,V)
        for t in range(times[k],times[k+1];length=11)
            geom=laser_geometry(reference(t),case;t)
            for e in 1:ne
                available[e,k]&=geom[e].available && geom[e].distance<=c.scheduling_range_m
            end
        end
    end
    terminal=ks_state_to_cartesian(view(reference(last(times)),slice(1)))
    M=Matrix(rtn_dcm_from_inertial(terminal.position_ii_m,terminal.velocity_ii_m))
    C=orbital_element_jacobian(terminal.position_ii_m,terminal.velocity_ii_m,c.planet.μ)*[M zeros(3,3);zeros(3,3) M]
    S=copy(C);gains=zeros(5,ne,n)
    for k in n:-1:1
        gains[:,:,k].=S*Vs[k];S=S*Fs[k]
    end
    # Original native-unit roundoff cutoffs, converted from km to SI for a.
    for (j,threshold) in enumerate([1e-9,1e-16,1e-14,1e-14,1e-14])
        gains[j,:,:].=map(x->abs(x)<threshold ? 0. : x,gains[j,:,:])
    end
    weights=objective==:combined ? [1/3000,0.,1/(3deg2rad(1)),0.,1/(3deg2rad(1))] :
        Float64.((1:5).==(objective==:a ? 1 : objective==:i ? 3 : 5))
    scores=dropdims(sum(gains.*reshape(weights,5,1,1);dims=1);dims=1)
    # The shared target terminal makes native maximum-weight matching exactly
    # the paper's at-most-one-positive-link problem, not an instantaneous policy.
    policy=TerminalGainPolicy(times,case.keys,scores)
    scheduler=SchedulingPolicyModel(policy);graph=deepcopy(case.interlinks)
    command=zeros(Int,ne,n)
    for k in 1:n
        t=(times[k]+times[k+1])/2
        graph.graph.available.=available[:,k]
        score_candidates!(graph,scheduler,interlink_state(reference(t),case),c.planet.μ;
            spacecraft=case.spacecraft,context=link_context(case,t))
        select_interlinks!(graph,scheduler)
        command[:,k].=graph.graph.active
    end
    predicted=zeros(5,n+1);dx=zeros(6)
    for k in 1:n
        dx=Fs[k]*dx+Vs[k]*command[:,k]
        x=ks_state_to_cartesian(view(reference(times[k+1]),slice(1)))
        M=Matrix(rtn_dcm_from_inertial(x.position_ii_m,x.velocity_ii_m))
        predicted[:,k+1].=orbital_element_jacobian(x.position_ii_m,x.velocity_ii_m,c.planet.μ)*[M zeros(3,3);zeros(3,3) M]*dx
    end
    replay_policy=TerminalGainPolicy(times,case.keys,scores.*command)
    return (;command,reference,predicted,gains,available,scores,weights,C,policy=replay_policy,times=collect(times))
end
end
