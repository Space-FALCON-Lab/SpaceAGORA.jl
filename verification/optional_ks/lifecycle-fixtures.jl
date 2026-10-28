include("baseline-fixtures.jl")
using Serialization
const L=SpaceAGORA.SimulationLifecycle
const SC=S.SimulationCallbacks
mutable struct PrescribedLifecycleControl <: S.AbstractTypes.AbstractControlEffectorModel
    root::S.Link
    held_force_N::Float64
    ticks::Vector{Float64}
    actions::Vector{Tuple{Float64,Symbol}}
    roots::Vector{Vector{Float64}}
    clocks::Vector{Tuple{Float64,Float64}}
    trace::Vector{Tuple{Float64,Symbol}}
    initial_audits::Vector{NamedTuple}
end
PrescribedLifecycleControl(root)=PrescribedLifecycleControl(root,0.,Float64[],Tuple{Float64,Symbol}[],Vector{Float64}[],Tuple{Float64,Float64}[],Tuple{Float64,Symbol}[],NamedTuple[])
SpaceAGORA.trial_rhs_safe(::PrescribedLifecycleControl)=true
function SpaceAGORA.update_control!(c::PrescribedLifecycleControl,u,p,t::Float64,i::Int64)
    c.held_force_N=.02*sin(t/1000.)
    push!(c.ticks,t)
    push!(c.trace,(t,:control))
    return nothing
end
SpaceAGORA.control_wrench(c::PrescribedLifecycleControl,x,p,t)=(S.SVector(c.held_force_N,0.,0.),S.SVector(0.,0.,0.))
SpaceAGORA.control_propellant_rate(::PrescribedLifecycleControl,x,p,t)=0.
SpaceAGORA.gnc_clock_started!(c::PrescribedLifecycleControl,period,origin)=(push!(c.clocks,(period,origin));nothing)
L.checkpoint_state(c::PrescribedLifecycleControl)=(held_force_N=c.held_force_N,area=c.root.ref_area,
    ticks=copy(c.ticks),actions=copy(c.actions),roots=deepcopy(c.roots),trace=copy(c.trace))
function L.restore_checkpoint_state!(c::PrescribedLifecycleControl,data)
    c.held_force_N=data.held_force_N;c.root.ref_area=data.area
    empty!(c.ticks);append!(c.ticks,data.ticks)
    empty!(c.actions);append!(c.actions,data.actions)
    empty!(c.roots);append!(c.roots,deepcopy(data.roots))
    empty!(c.trace);append!(c.trace,data.trace)
    return nothing
end
struct LifecycleStateTimeForce <: S.AbstractTypes.AbstractForceTorqueModel end
SpaceAGORA.trial_rhs_safe(::LifecycleStateTimeForce)=true
function S.EffectorSampling.wrench(::LifecycleStateTimeForce,x::S.EffectorSampling.StateSample,
        env::S.EffectorSampling.EnvironmentSample,t::Float64)
    acceleration=1e-7*sin(t/700.)*(x.pos_ii/7e6)+S.SVector(0.,2e-8*cos(t/900.),0.)
    return x.mass_kg*acceleration,S.SVector(0.,0.,0.)
end
function lifecycle_config(mode;directory=mktempdir(),checkpoint=false,resume=false,coincident=false)
    duration=coincident ? 442. : 6000.
    cfg=baseline_config(:J2_plus_drag,mode;duration,rate=17.,rtol=mode===:ks ? 1e-12 : 2.5e-13,
        atol=mode===:ks ? 1e-13 : 2.5e-13,dtmax=5.)
    root=S.Link(root=true,m=480.,ref_area=10.,α=pi/2,β=0.,θ=0.)
    original=only(cfg.dynamics_model.spacecraft)
    spacecraft=S.SpacecraftModel(joints=S.Joint[],links=[root],root=root,prop_mass=20.,
        inertia_tensor=root.inertia,initial_condition=original.initial_condition,id=1)
    c=PrescribedLifecycleControl(root)
    burn=S.ThrusterModels.BaseThrusterModel(thrust=[coincident ? 0. : .2],direction=[0.],Δv=[0.],
        start_burn_time=[611.25],stop_burn_time=[1879.125],Isp=[300.],actuation=:ideal)
    return S.SimConfig._with_configuration(cfg;
        dynamics_model=S.DynamicsModel([spacecraft],(S.InverseSquaredJ2GravityModel(),
            S.AerodynamicCoefficientConstant(),LifecycleStateTimeForce())),
        control_model=coincident ? S.ControlModel(control_effectors=(c,),control_rates=[13.]) :
            S.ControlModel(control_effectors=(c,burn),control_rates=[13.,13.]),
        simulation_settings=S.SimulationSettings(results=false,verbose=false,generate_plots=false,normalize=false,
            results_directory=directory,save_csv=false,checkpoint_enabled=checkpoint,
            checkpoint_interval_s=coincident ? 221. : 997.3,checkpoint_directory=joinpath(directory,"checkpoints"),
            resume_from_checkpoint=resume))
end
function lifecycle_callbacks(;interrupt=false,coincident=false)
    function factory(args)
        c=args.control_model.control_effectors[1]
        epochs=coincident ? [221.] : [611.25,1879.125,3599.75]
        labels=coincident ? [:impulse] : [:area15,:impulse,:area8]
        interruption_time=coincident ? 300. : 1500.
        applied(label)=any(last(event)===label for event in c.actions)
        function initialize(cb,u,t,int)
            controls=int.p.args.control_model.control_effectors
            push!(c.initial_audits,(epoch=Float64(t),area=c.root.ref_area,
                held_force_N=c.held_force_N,last_tick=isempty(c.ticks) ? NaN : last(c.ticks),
                burn_active=length(controls)>1 ? copy(controls[2].burning) : Bool[],
                root_alias=c.root===only(int.p.args.dynamics_model.spacecraft).root,
                committed_actions=copy(c.actions)))
            segment_end=int.p.shared_buffers.solve_segment_end_time[]
            for epoch in epochs
                t<epoch<=segment_end && SE.DiffEqBase.add_tstop!(int,epoch)
            end
            interrupt && t<interruption_time<=segment_end && SE.DiffEqBase.add_tstop!(int,interruption_time)
            return nothing
        end
        function condition(u,t,int)
            interrupt && t==interruption_time && return true
            return any(abs(t-ep)<=1e-9 && !applied(label) for (ep,label) in zip(epochs,labels))
        end
        function affect!(int)
            interrupt && int.t==interruption_time && error("deliberate lifecycle interruption at$(int.t)s")
            for (ep,label) in zip(epochs,labels)
                if abs(int.t-ep)<=1e-9 && !applied(label)
                    if label===:area15
                        c.root.ref_area=15.
                    elseif label===:area8
                        c.root.ref_area=8.
                    else
                        int.u.sc[1].vel[2]+=.05
                    end
                    push!(c.actions,(int.t,label))
                    push!(c.trace,(int.t,label))
                end
            end
            return nothing
        end
        discrete=SE.DiscreteCallback(condition,affect!;initialize,save_positions=(true,true))
        function down!(int)
            push!(c.roots,vcat(int.t,Vector(int.u)))
            push!(c.trace,(int.t,:descending_root))
            SE.DiffEqBase.derivative_discontinuity!(int,false)
            return nothing
        end
        continuous=SC.ContinuousCallback((u,t,int)->u.sc[1].pos[3],nothing,down!;
            save_positions=(true,true),abstol=1e-8,reltol=0.)
        return SE.CallbackSet(continuous,discrete)
    end
    return factory
end
function lifecycle_table_state(table,row)
    return vcat([table[row,"sc1_pos_$j"] for j in 1:3],
        [table[row,"sc1_vel_$j"] for j in 1:3],table[row,"sc1_mass"],table[row,"sc1_heat_load"])
end

function lifecycle_save_fields(args)
    fields=SC.default_save_fields(args)
    function observe_output(u,t,int)
        c=int.p.args.control_model.control_effectors[1]
        push!(c.trace,(Float64(t),:output))
        return c.held_force_N
    end
    push!(fields,SC.SaveField(:lifecycle_held_force_N,observe_output))
    return fields
end
