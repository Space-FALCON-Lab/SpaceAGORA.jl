module GravityGradientTorqueTests
using Test, LinearAlgebra, StaticArrays, ComponentArrays, SpaceAGORA
using SpaceAGORA.SimulationModel
const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const ZERO3 = SVector(0.0,0.0,0.0)
const INERTIA = SMatrix{3,3,Float64}(2.0,0.15,-0.08, 0.15,3.0,0.11, -0.08,0.11,4.0)
const Q0 = normalize(SVector(0.13,-0.21,0.17,1.0))
const W0 = SVector(0.0002,-0.0001,0.0003)
const DURATION = 20.0

# Independent vector form of the scalar-last inertial-to-body rotation.
function reference_torque(r, q, inertia, mu)
    v, s = SVector(q[1],q[2],q[3]), q[4]
    rb = (s*s-dot(v,v))*r + 2v*dot(v,r) - 2s*cross(v,r)
    n = rb/norm(rb)
    return (3mu/norm(r)^3)*cross(n,inertia*n)
end

function configuration(mode; orientation=true)
    planet = Earth()
    spacecraft = SpacecraftModel[]
    for (i,id) in enumerate((41,97))
        inertia = i*INERTIA
        root = Link(root=true,m=500.0,ref_area=12.0,inertia=inertia)
        ic = InitialCondition(ra=planet.Rp_e+500e3+100i,rp=planet.Rp_e+500e3+100i,
            i=53.0,ω=24.0,Ω=10.0,ν=70.0(i-1),q=Q0,ang_vel=W0)
        push!(spacecraft,SpacecraftModel(Joint[],[root],root,true,500.0,0.0,inertia,0,0,ic,id))
    end
    effectors = mode == :separate ? (InverseSquaredGravityModel(),SM.GravityGradientTorqueModel()) :
        mode == :builtin ? (InverseSquaredGravityModel(gravity_gradient=true),) :
        mode == :disabled ? (InverseSquaredGravityModel(),SM.GravityGradientTorqueModel(gravity_gradient=false)) :
        (InverseSquaredGravityModel(),)
    SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false,verbose=false,generate_plots=false,normalize=false,save_csv=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime,keplerian=true,number_of_orbits=1,
            mission_time=DURATION,orientation_sim=orientation,num_steps_to_save=30,data_rate=1.0),
        environment_model=EnvironmentModel(planet=planet,EI=300.0,density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0,planet=planet),
            topography=false,wind=false,ephemerides_model=SimpleEphemeridesModel()),
        dynamics_model=DynamicsModel(spacecraft,effectors),
        guidance_model=GuidanceModel(guidance_effectors=(),guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(),navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(),control_rates=Float64[]),
        initial_time=InitialTime(year=2020,month=1,day=1,hour=0,minute=0,second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-10,abstol_orbit=1e-10,
            reltol_quaternion=1e-11,abstol_quaternion=1e-12,reltol_angular_rate=1e-11,
            abstol_angular_rate=1e-12,dt_max_orbit=0.25),solver_config=SolverConfig(solver_mode=:tsit5))
end

function rhs(args)
    u = SE.build_initial_conditions(args)
    p = ODEParams(n_sats=2,args=deepcopy(args))
    du = zero(u)
    p.shared_buffers.rhs_env_config[] = withenv("SPACEAGORA_RHS_EXECUTION_MODE"=>"serial") do
        SE._snapshot_rhs_plan_env_config()
    end
    SE.spacecraft_dynamics!(du,u,p,0.0)
    return u,du,p
end

@testset "torque-only gravity uses the existing owner and independent body-frame law" begin
    model = SM.GravityGradientTorqueModel()
    args = configuration(:separate)
    sc = args.dynamics_model.spacecraft[1]
    env = EnvironmentSample(args.environment_model.planet)
    @test SM.DynamicEffectors.GravityGradientTorqueModel === SM.GravityEffectors.GravityGradientTorqueModel
    @test model.gravity_gradient
    @test SE._dynamic_effector_threadsafe(model)
    @test SM.solver_partition(model) == :explicit
    for r in (SVector(7e6,2e6,-1e6),SVector(-2e6,6e6,3e6)), q in (Q0,-Q0,normalize(SVector(-0.3,0.1,0.4,0.8)))
        x = StateSample(r,ZERO3,500.0;q_ib=q,ω_body=W0,spacecraft=sc)
        expected = reference_torque(r,q,INERTIA,env.planet.μ)
        before = deepcopy((sc.inertia_tensor,sc.root.r,sc.root.q,sc.root.ω))
        f,tau = SM.wrench(model,x,env,0.0)
        @test f == ZERO3
        @test norm(expected)>1e-8
        @test tau ≈ expected rtol=2e-14 atol=1e-18
        @test (sc.inertia_tensor,sc.root.r,sc.root.q,sc.root.ω)==before
        @test SM.wrench(SM.GravityGradientTorqueModel(gravity_gradient=false),x,env,0.0)==(ZERO3,ZERO3)
    end
    for x in (StateSample(SVector(7e6,0.0,0.0),ZERO3,500.0;spacecraft=sc),
              StateSample(SVector(7e6,0.0,0.0),ZERO3,500.0;q_ib=Q0),
              StateSample(ZERO3,ZERO3,500.0;q_ib=Q0,spacecraft=sc),
              StateSample(SVector(NaN,0.0,0.0),ZERO3,500.0;q_ib=Q0,spacecraft=sc))
        @test SM.wrench(model,x,env,0.0)==(ZERO3,ZERO3)
    end
end

@testset "actual RHS applies one gravity force and one torque per spacecraft" begin
    args = configuration(:separate)
    u,du,p = rhs(args)
    _,builtin,_ = rhs(configuration(:builtin))
    _,bare,_ = rhs(configuration(:none))
    _,disabled,_ = rhs(configuration(:disabled))
    @test [sc.id for sc in args.dynamics_model.spacecraft] == [41,97]
    for i in eachindex(u.sc)
        sc = args.dynamics_model.spacecraft[i]
        x = u.sc[i]
        r = SVector{3,Float64}(x.pos)
        tau = reference_torque(r,SVector{4,Float64}(x.q),sc.inertia_tensor,args.environment_model.planet.μ)
        expected_rate = sc.inertia_tensor \ (tau-cross(W0,sc.inertia_tensor*W0))
        @test SVector{3,Float64}(du.sc[i].ω) ≈ expected_rate rtol=1e-12 atol=1e-16
        @test norm(SVector{3,Float64}(du.sc[i].ω)-SVector{3,Float64}(bare.sc[i].ω))>1e-8
        @test du.sc[i].vel ≈ bare.sc[i].vel rtol=2e-14 atol=1e-13
        @test SM.calcForceTorque(SM.GravityGradientTorqueModel(),x,p,Int64(i))[1] == ZERO3
        @test SM.calcForceTorque(SM.GravityGradientTorqueModel(),x,p,Int64(i))[2] ≈ tau rtol=2e-14 atol=1e-18
    end
    @test collect(ComponentArrays.getdata(du)) ≈ collect(ComponentArrays.getdata(builtin)) rtol=1e-12 atol=1e-13
    @test collect(ComponentArrays.getdata(disabled)) ≈ collect(ComponentArrays.getdata(bare)) rtol=1e-12 atol=1e-13
    @test SM.calcForceTorque(SM.GravityGradientTorqueModel(),u.sc[1],p,Int64(0))==(ZERO3,ZERO3)
    @test SM.calcForceTorque(SM.GravityGradientTorqueModel(),u.sc[1],p,Int64(3))==(ZERO3,ZERO3)
    noatt,du_noatt,p_noatt = rhs(configuration(:separate;orientation=false))
    _,baseline_noatt,_ = rhs(configuration(:none;orientation=false))
    @test collect(ComponentArrays.getdata(du_noatt)) ≈ collect(ComponentArrays.getdata(baseline_noatt)) rtol=1e-12 atol=1e-13
    @test SM.calcForceTorque(SM.GravityGradientTorqueModel(),noatt.sc[1],p_noatt,Int64(1))==(ZERO3,ZERO3)
end

function run_recorded(args; execution_mode="serial", return_runtime=false)
    fields=SaveField[
        SaveField(:q,(u,t,int)->[SVector{4,Float64}(sc.q) for sc in u.sc];per_satellite=true),
        SaveField(:omega,(u,t,int)->[SVector{3,Float64}(sc.ω) for sc in u.sc];per_satellite=true),
        SaveField(:position,(u,t,int)->[SVector{3,Float64}(sc.pos) for sc in u.sc];per_satellite=true)]
    recorder=TrajectoryRecorder(args;save_fields=fields)
    result=withenv("SPACEAGORA_RHS_EXECUTION_MODE"=>execution_mode,"SPACEAGORA_RHS_CALIBRATE"=>"off",
        "SPACEAGORA_RHS_IDENTIFY"=>"0","SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS"=>"0",
        "SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST"=>"0") do
        run_simulation(args;return_solver_metadata=true,return_solution=return_runtime,visualization=false,
            extra_callbacks=(get_trajectory_recorder_callback(recorder),))
    end
    @test result.retcode=="Success"
    @test recorder.count>=15
    @test first(trajectory_times(recorder))≈0.0 atol=1e-12
    @test last(trajectory_times(recorder))≈DURATION atol=1e-12
    times, rows = collect(trajectory_times(recorder)),trajectory_save_data(recorder)
    return return_runtime ? (times=times, rows=rows, runtime=result.solution.prob.p) : (times,rows)
end

@testset "separate and built-in torque agree in real attitude propagation" begin
    args=configuration(:separate)
    initial=copy(ComponentArrays.getdata(SE.build_initial_conditions(args)))
    roots=deepcopy([(sc.root.r,sc.root.q,sc.root.ω) for sc in args.dynamics_model.spacecraft])
    times,separate=run_recorded(args)
    builtin_times,builtin=run_recorded(configuration(:builtin))
    bare_times,bare=run_recorded(configuration(:none))
    @test times==builtin_times==bare_times
    @test ComponentArrays.getdata(SE.build_initial_conditions(args))==initial
    @test [(sc.root.r,sc.root.q,sc.root.ω) for sc in args.dynamics_model.spacecraft]==roots
    for (a,b) in zip(separate,builtin), i in 1:2
        @test a[:q][i] ≈ b[:q][i] rtol=0.0 atol=1e-10
        @test a[:omega][i] ≈ b[:omega][i] rtol=0.0 atol=1e-11
        @test a[:position][i] ≈ b[:position][i] rtol=0.0 atol=1e-6
    end
    for i in 1:2
        @test norm(separate[end][:omega][i]-bare[end][:omega][i])>1e-7
    end
end

@testset "stateless gravity torque supports concurrent calls and the real flat route" begin
    args=configuration(:separate)
    model=SM.GravityGradientTorqueModel()
    craft=args.dynamics_model.spacecraft[1]
    sample=StateSample(SVector(7e6,2e6,-1e6),ZERO3,500.0;
        q_ib=Q0,ω_body=W0,spacecraft=craft)
    environment=EnvironmentSample(args.environment_model.planet)
    expected=SM.wrench(model,sample,environment,0.0)
    before=deepcopy((craft.inertia_tensor,craft.root.r,craft.root.q,craft.root.ω))
    # Reuse the same model and read-only physical inputs on all tasks.
    outputs=Vector{Vector{typeof(expected)}}(undef,16)
    @sync for k in eachindex(outputs)
        Threads.@spawn outputs[k]=[SM.wrench(model,sample,environment,Float64(j)) for j in 1:64]
    end
    @test all(values->all(==(expected),values),outputs)
    @test (craft.inertia_tensor,craft.root.r,craft.root.q,craft.root.ω)==before
    if Threads.nthreads() < 2
        @test_skip Threads.nthreads() >= 2
    else
        planet=args.environment_model.planet
        coefficients=joinpath(@__DIR__,"..","..","..","data","Gravity_harmonics_data","EarthGGM05C.csv")
        harmonics=GravitationalHarmonicsModel(4,4,coefficients,planet)
        args=SM.SimConfig._with_configuration(args;dynamics_model=
            DynamicsModel(args.dynamics_model.spacecraft,(harmonics,model)))
        withenv("SPACEAGORA_INNER_THREAD_BUDGET"=>"2",
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE"=>"0",
            "SPACEAGORA_PARALLEL_PROFILE"=>nothing,
            "SPACEAGORA_EFFECTOR_PARALLEL"=>"on",
            "SPACEAGORA_EFFECTOR_PARALLEL_HEAVY_ONLY"=>"0",
            "SPACEAGORA_EFFECTOR_THREAD_THRESHOLD"=>"1",
            "SPACEAGORA_EFFECTOR_MAX_THREADS"=>"2",
            # Keep auto-routing deterministic for this small regression fixture.
            "SPACEAGORA_EFFECTOR_FLAT_MIN_SATS"=>"2",
            "SPACEAGORA_EFFECTOR_FLAT_MIN_THREAD_BUDGET"=>"2",
            "SPACEAGORA_EFFECTOR_FLAT_WORK_NS_THRESHOLD"=>"1",
            "SPACEAGORA_RHS_FLAT_WORK_PER_WORKER_NS_THRESHOLD"=>"1") do
            serial=run_recorded(args;execution_mode="serial",return_runtime=true)
            flat=run_recorded(args;execution_mode="flat",return_runtime=true)
            automatic=run_recorded(args;execution_mode="auto",return_runtime=true)
            effects=flat.runtime.args.dynamics_model.dynamic_effectors
            serial_plan=SE._rhs_execution_plan(serial.runtime.args,serial.runtime,effects,2)
            flat_plan=SE._rhs_execution_plan(flat.runtime.args,flat.runtime,effects,2)
            @test serial_plan.mode===:serial
            @test flat_plan.mode===:flat_constellation_effector_queue
            @test flat_plan.allotment==2
            auto_plan=SE._rhs_execution_plan(automatic.runtime.args,automatic.runtime,effects,2)
            @test auto_plan.mode===:flat_constellation_effector_queue
            @test auto_plan.allotment==2
            # Only the actual flat dispatcher fills these per-effector slots.
            slots=flat.runtime.shared_buffers.rhs_flat_effector_partials[]
            @test size(slots)==(6,2,2)
            if size(slots)==(6,2,2)
                @test all(iszero,@view slots[1:3,2,:])
                @test any(!iszero,@view slots[4:6,2,:])
            end
            @test isempty(serial.runtime.shared_buffers.rhs_flat_effector_partials[])
            auto_slots=automatic.runtime.shared_buffers.rhs_flat_effector_partials[]
            @test size(auto_slots)==(6,2,2)
            if size(auto_slots)==(6,2,2)
                @test any(!iszero,@view auto_slots[4:6,2,:])
            end
            bits(xs)=reinterpret(UInt64,Float64.(collect(xs)))
            @test bits(serial.times)==bits(flat.times)
            @test length(serial.rows)==length(flat.rows)
            @test all(zip(serial.rows,flat.rows)) do (left,right)
                all(key->all(i->bits(left[key][i])==bits(right[key][i]),1:2),(:q,:omega,:position))
            end
            @test bits(serial.times)==bits(automatic.times)
            @test length(serial.rows)==length(automatic.rows)
            @test all(zip(serial.rows,automatic.rows)) do (left,right)
                all(key->all(i->bits(left[key][i])==bits(right[key][i]),1:2),(:q,:omega,:position))
            end
        end
    end
end
end
