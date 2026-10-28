include("baseline-fixtures.jl")
using KSPropagator, Test
mutable struct ManufacturedHeldForce <: S.AbstractTypes.AbstractForceTorqueModel
    acceleration::Float64
end
SpaceAGORA.trial_rhs_safe(::ManufacturedHeldForce)=true
function S.EffectorSampling.wrench(model::ManufacturedHeldForce,x::S.EffectorSampling.StateSample,
        env::S.EffectorSampling.EnvironmentSample,t::Float64)
    acceleration=env.planet.μ*x.pos_ii/norm(x.pos_ii)^3+S.SVector(model.acceleration,0.,0.)
    return x.mass_kg*acceleration,S.SVector(0.,0.,0.)
end
function manufactured_truth(t;impulse=0.,side=:right)
    r,v=orbit_initial();r=r+v*t+[.02*t^2/2,0.,0.];v=v+[.02*t,0.,0.]
    for (epoch,change,jump) in ((1.25,-.05,impulse),(2.75,.07,0.),(3.25,-.03,0.))
        if t>epoch || (t==epoch && side==:right)
            h=t-epoch
            r .+= [change*h^2/2+jump*h,0.,0.]
            v .+= [change*h+jump,0.,0.]
        end
    end
    return vcat(r,v)
end
function run_manufactured(impulse)
    cfg=baseline_config(:two_body,:ks;duration=4.,rate=.05,rtol=1e-11,atol=1e-12,dtmax=.02)
    force=ManufacturedHeldForce(.02)
    cfg=S.SimConfig._with_configuration(cfg;dynamics_model=S.DynamicsModel(cfg.dynamics_model.spacecraft,
        (S.InverseSquaredGravityModel(),force)))
    events=Tuple{Float64,Symbol}[];snapshot=Ref{Any}(nothing)
    function factory(args)
        command=args.dynamics_model.dynamic_effectors[2]
        first_affect=int->begin
            command.acceleration=-.03
            int.u.sc[1].vel[1]+=impulse
            push!(events,(int.t,:first))
        end
        second_affect=int->begin
            snapshot[]=copy(int.sol(.625))
            command.acceleration=.04
            push!(events,(int.t,:second))
        end
        discrete_affect=int->begin
            command.acceleration=.01
            push!(events,(int.t,:discrete))
        end
        discrete_init=(cb,u,t,int)->SE.DiffEqBase.add_tstop!(int,3.25)
        cb1=S.SimulationCallbacks.ContinuousCallback((u,t,int)->t-1.25,first_affect,nothing;
            save_positions=(true,true),abstol=1e-13,reltol=0.)
        cb2=S.SimulationCallbacks.ContinuousCallback((u,t,int)->t-2.75,second_affect,nothing;
            save_positions=(true,true),abstol=1e-13,reltol=0.)
        cb3=SE.DiscreteCallback((u,t,int)->t==3.25,discrete_affect;
            initialize=discrete_init,save_positions=(true,true))
        return SE.CallbackSet(cb1,cb2,cb3)
    end
    sol=run_simulation(cfg;return_solution=true,extra_callback_factories=(factory,))
    @test length(events)==3
    @test [last(x) for x in events]==[:first,:second,:discrete]
    @test maximum(abs.([first(x) for x in events]-[1.25,2.75,3.25]))<=1e-10
    maxpos=0.;maxvel=0.
    # Query all old intervals AFTER every command has run. At each event, test
    # both public knot continuities and strict interior samples on both sides.
    for t in (.1,.625,1.249,1.25,1.251,2.,2.749,2.75,2.751,3.249,3.25,3.251,3.999)
        for side in (:left,:right)
            got=sol(t;continuity=side);truth=manufactured_truth(t;impulse,side)
            pe=norm(got[1:3]-truth[1:3]);ve=norm(got[4:6]-truth[4:6])
            maxpos=max(maxpos,pe);maxvel=max(maxvel,ve)
            @test pe<=1e-3
            @test ve<=1e-5
            # Allocating and in-place/component-selection dense paths agree.
            out=zero(got);sol(out,t;continuity=side)
            @test out≈got atol=1e-12 rtol=0.
            @test sol(t;idxs=1,continuity=side)≈got[1] atol=1e-12 rtol=0.
        end
    end
    @test snapshot[]≈sol(.625) atol=1e-12 rtol=0.
    @test sol.t[end]==4.
    println("DENSE_HISTORY_CASE impulse=",impulse," max_position_error_m=",maxpos,
        " max_velocity_error_m_s=",maxvel," attempts=",sol.alg.diagnostics.attempts,
        " root_events=",events);flush(stdout)
end
@testset "Retained physical interpolation across command changes and impulses" begin
    run_manufactured(0.)
    run_manufactured(.7)
end
