module MeshAeroThreadingTests
using Test, LinearAlgebra, StaticArrays, SpaceAGORA
using SpaceAGORA.SimulationModel
const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const Z3 = SVector(0.0, 0.0, 0.0)
const Q0 = normalize(SVector(0.13, -0.21, 0.17, 1.0))

function configuration(; orientation=true)
    planet = make_no_gram_planet(:earth)
    inertia = SMatrix{3,3,Float64}(20.0, 0.0, 0.0, 0.0, 30.0, 0.0, 0.0, 0.0, 40.0)
    # Direction and speed dependent coefficients exercise the per-call SH buffer.
    a = reshape([0.01sin(i) for i in 1:108], 6, 18)
    a[1,1] = -1.0
    sur = MeshAeroSurrogate(2, 1, a, 0.1a, 2.0, 1.0,
        SVector(0.0, 0.2, -0.1), 1.0, 1.0, 3.0, 20.0)
    mesh = AerodynamicCoefficientMeshSurrogate(Dict(1=>sur, 2=>sur))
    spacecraft = SpacecraftModel[]
    for (i,id) in enumerate((41,97))
        root = Link(root=true, m=450.0, ref_area=2.0, inertia=inertia)
        child = Link(root=false, m=50.0, r=[0.0, 1.5i, 0.2],
            q=collect(normalize(SVector(0.0, 0.12i, 0.0, 1.0))), inertia=inertia)
        ic = InitialCondition(ra=planet.Rp_e+180e3+100i, rp=planet.Rp_e+180e3+100i,
            i=53.0, ω=24.0, Ω=10.0, ν=70.0(i-1), q=Q0,
            ang_vel=SVector(0.0002i, -0.0001, 0.0003))
        push!(spacecraft, SpacecraftModel(Joint[], [root,child], root, true,
            500.0, 0.0, inertia, 0, 0, ic, id))
    end
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false,
            generate_plots=false, normalize=false, save_csv=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime,
            keplerian=true, number_of_orbits=1, mission_time=4.0,
            orientation_sim=orientation, num_steps_to_save=30, data_rate=0.5),
        environment_model=EnvironmentModel(planet=planet, EI=300.0,
            density_model=ExponentialAtmosphereModel(1e-10,180e3,50e3;temperature_k=250.0),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0,planet=planet),
            topography=false,wind=false,ephemerides_model=SimpleEphemeridesModel()),
        dynamics_model=DynamicsModel(spacecraft,(InverseSquaredGravityModel(),mesh)),
        guidance_model=GuidanceModel(guidance_effectors=(),guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(),navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(),control_rates=Float64[]),
        initial_time=InitialTime(year=2020,month=1,day=1,hour=0,minute=0,second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-10,abstol_orbit=1e-10,
            reltol_quaternion=1e-11,abstol_quaternion=1e-12,reltol_angular_rate=1e-11,
            abstol_angular_rate=1e-12,dt_max_orbit=0.25),
        solver_config=SolverConfig(solver_mode=:tsit5))
end

function run_recorded(args, mode)
    fields = SaveField[
        SaveField(:pos,(u,t,int)->[SVector{3,Float64}(sc.pos) for sc in u.sc];per_satellite=true),
        SaveField(:vel,(u,t,int)->[SVector{3,Float64}(sc.vel) for sc in u.sc];per_satellite=true)]
    if args.mission_configuration.orientation_sim
        append!(fields, [
            SaveField(:q,(u,t,int)->[SVector{4,Float64}(sc.q) for sc in u.sc];per_satellite=true),
            SaveField(:omega,(u,t,int)->[SVector{3,Float64}(sc.ω) for sc in u.sc];per_satellite=true)])
    end
    for (name,cache) in ((:drag,:drag_cache),(:lift,:lift_cache),(:cross,:cross_cache))
        push!(fields,SaveField(name,(u,t,int)->[SVector{3,Float64}(v)
            for v in getproperty(int.p.save_cache,cache)];per_satellite=true))
    end
    recorder=TrajectoryRecorder(args;save_fields=fields)
    result=withenv("SPACEAGORA_RHS_EXECUTION_MODE"=>mode,
        "SPACEAGORA_RHS_CALIBRATE"=>"off","SPACEAGORA_RHS_IDENTIFY"=>"0",
        "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS"=>"0",
        "SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST"=>"0") do
        run_simulation(args;return_solver_metadata=true,return_solution=true,visualization=false,
            extra_callbacks=(get_trajectory_recorder_callback(recorder),))
    end
    @test result.retcode=="Success"
    @test last(trajectory_times(recorder)) ≈ 4.0 atol=1e-12
    return (times=collect(trajectory_times(recorder)),rows=trajectory_save_data(recorder),
        runtime=result.solution.prob.p)
end

@testset "mesh coefficients are shared read-only and aerodynamic diagnostics have one writer" begin
    args=configuration(); mesh=args.dynamics_model.dynamic_effectors[2]
    @test SE._dynamic_effector_threadsafe(mesh)
    @test SE._dynamic_effectors_parallel_supported(args.dynamics_model.dynamic_effectors)
    for effectors in ((mesh,mesh),(mesh,AerodynamicCoefficientfM()),
        (AerodynamicCoefficientfM(),mesh))
        @test !SE._dynamic_effectors_parallel_supported(effectors)
        @test !SE._rhs_flat_supported(effectors)
    end
    before=deepcopy(mesh.surrogates[1].coeff_a)
    inputs=[(normalize(SVector(0.2+k,0.3,-0.1k)),3.0+0.1k) for k in 1:32]
    reference=[mesh_aero_coefficients(mesh.surrogates[1],v,s) for (v,s) in inputs]
    outputs=Vector{typeof(reference)}(undef,16)
    @sync for k in eachindex(outputs)
        Threads.@spawn outputs[k]=[mesh_aero_coefficients(mesh.surrogates[1],v,s) for (v,s) in inputs]
    end
    @test all(==(reference),outputs)
    @test mesh.surrogates[1].coeff_a==before
end

@testset "real mesh serial, forced-flat and automatic routes agree" begin
    if Threads.nthreads()<2
        @test_skip Threads.nthreads()>=2
    else
        withenv("SPACEAGORA_INNER_THREAD_BUDGET"=>"2",
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE"=>"0","SPACEAGORA_PARALLEL_PROFILE"=>nothing,
            "SPACEAGORA_EFFECTOR_PARALLEL"=>"on","SPACEAGORA_EFFECTOR_PARALLEL_HEAVY_ONLY"=>"0",
            "SPACEAGORA_EFFECTOR_THREAD_THRESHOLD"=>"1","SPACEAGORA_EFFECTOR_MAX_THREADS"=>"2",
            "SPACEAGORA_EFFECTOR_FLAT_MIN_SATS"=>"2","SPACEAGORA_EFFECTOR_FLAT_MIN_THREAD_BUDGET"=>"2",
            "SPACEAGORA_EFFECTOR_FLAT_WORK_NS_THRESHOLD"=>"1",
            "SPACEAGORA_RHS_FLAT_WORK_PER_WORKER_NS_THRESHOLD"=>"1") do
            for orientation in (false,true)
                args=configuration(;orientation)
                before=deepcopy(args.dynamics_model.dynamic_effectors[2].surrogates[1].coeff_a)
                serial=run_recorded(args,"serial")
                bits(x)=reinterpret(UInt64,Float64.(collect(x)))
                for mode in ("flat","auto")
                    run=run_recorded(args,mode)
                    effects=run.runtime.args.dynamics_model.dynamic_effectors
                    plan=SE._rhs_execution_plan(run.runtime.args,run.runtime,effects,2)
                    @test plan.mode===:flat_constellation_effector_queue
                    @test plan.allotment==2
                    slots=run.runtime.shared_buffers.rhs_flat_effector_partials[]
                    @test size(slots)==(6,2,2)
                    if size(slots)==(6,2,2)
                        @test any(!iszero,@view slots[1:3,2,:])
                        @test orientation ? any(!iszero,@view slots[4:6,2,:]) : all(iszero,@view slots[4:6,2,:])
                    end
                    @test bits(run.times)==bits(serial.times)
                    @test length(run.rows)==length(serial.rows)
                    @test all(zip(serial.rows,run.rows)) do (a,b)
                        all(key->all(i->bits(a[key][i])==bits(b[key][i]),1:2),keys(a))
                    end
                    @test any(row->norm(row[:drag][1])>0.0,run.rows)
                    @test any(row->row[:drag][1]!=row[:drag][2],run.rows)
                end
                @test isempty(serial.runtime.shared_buffers.rhs_flat_effector_partials[])
                @test args.dynamics_model.dynamic_effectors[2].surrogates[1].coeff_a==before
            end
        end
    end
end
end
