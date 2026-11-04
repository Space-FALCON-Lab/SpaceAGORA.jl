module CygnssCommandReplayTests
using Test, LinearAlgebra, StaticArrays, SpaceAGORA
include(joinpath(@__DIR__, "../../../scripts/dev/viewer_demos/cygnss_command_replay.jl"))
using .CygnssCommandReplay
const SM = SpaceAGORA.SimulationModel
@testset "Recorded wheel commands: momentum, sign and sampling" begin
    A = [0.0 -1.0; 1.0 0.0; 0.0 0.0]
    t = [0.0, 0.2, 1.3, 2.0]
    u = hcat(2 .* t, fill(3.0,length(t)))
    m = CommandedWheelReplay(t,u,A,[4.0,-2.0];time_offset_s=0.1)
    for tt in [0.0,0.1,0.55,1.2,1.9]
        x=tt+0.1; v=command_values(m,tt)
        @test v.rate ≈ A*[2x,3]
        @test v.momentum ≈ A*[4+x^2,-2+3x]
        @test v.axial_momentum ≈ [4+x^2,-2+3x]
        w=SVector(0.01,0.02,0.03)
        @test SM.wheel_reaction_torque(v.momentum,v.rate,w) ≈ -A*[2x,3]-cross(w,A*[4+x^2,-2+3x])
    end
    @test_throws DomainError command_values(m,-0.2)
    @test_throws DomainError command_values(m,2.0)
    @test_throws ArgumentError CommandedWheelReplay([0.,0.],zeros(2,2),A,[0.,0.])
    @test_throws ArgumentError CommandedWheelReplay([0.,1.],fill(NaN,2,2),A,[0.,0.])
    @test_throws ArgumentError CommandedWheelReplay([0.,1.],zeros(2,2),A,[0.])
end
using SpaceAGORA.SimulationModel
@testset "Command replay in the engine conserves total angular momentum" begin
    planet=Earth(); inertia=SMatrix{3,3,Float64}(2,0,0,0,3,0,0,0,4)
    q0=normalize(SVector(0.1,0.2,-0.1,1.0)); w0=SVector(0.01,-0.02,0.03)
    root=Link(root=true,m=500.0,ref_area=1.0,inertia=inertia)
    ic=InitialCondition(ra=planet.Rp_e+550000,rp=planet.Rp_e+550000,
        i=53.0,ω=0.0,Ω=10.0,ν=0.0,q=q0,ang_vel=w0)
    sc=SpacecraftModel(Joint[],[root],root,true,500.0,0.0,inertia,0,0,ic,97)
    A=[1.0 0.2; -0.1 0.9; 0.15 -0.2]
    cmd=CommandedWheelReplay([0.,1.,2.],[0.01 -0.02;0.02 -0.01;0.01 -0.02],A,[0.2,-0.1])
    args=SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false,verbose=false,generate_plots=false,normalize=false,save_csv=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime,keplerian=true,number_of_orbits=1,
            mission_time=2.0,orientation_sim=true,num_steps_to_save=30,data_rate=0.1),
        environment_model=EnvironmentModel(planet=planet,EI=300.0,density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0,planet=planet),
            topography=false,wind=false,ephemerides_model=SimpleEphemeridesModel()),
        dynamics_model=DynamicsModel([sc],(InverseSquaredGravityModel(),cmd)),
        guidance_model=GuidanceModel(guidance_effectors=(),guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(),navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(),control_rates=Float64[]),
        initial_time=InitialTime(year=2020,month=1,day=1,hour=0,minute=0,second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-10,abstol_orbit=1e-10,
            reltol_quaternion=1e-11,abstol_quaternion=1e-12,reltol_angular_rate=1e-11,
            abstol_angular_rate=1e-12,dt_max_orbit=0.02),solver_config=SolverConfig(solver_mode=:tsit5))
    sample=StateSample(SVector(7e6,0.,0.),SVector(0.,7500.,0.),500.;q_ib=q0,ω_body=w0,spacecraft=sc)
    env=EnvironmentSample(planet); v=command_values(cmd,0.3)
    @test SM.wrench(cmd,sample,env,0.3)[2] ≈ -v.rate-cross(w0,v.momentum)
    @test SM.wrench_caching!(cmd,sample,env,0.3,nothing,2)==(zero(w0),zero(w0))
    @test SM.wrench_caching!(cmd,sample,env,0.3,nothing,1)==SM.wrench(cmd,sample,env,0.3)
    fields=[SaveField(:q,(u,t,it)->[SVector{4,Float64}(u.sc[1].q)];per_satellite=true),
            SaveField(:omega,(u,t,it)->[SVector{3,Float64}(u.sc[1].ω)];per_satellite=true)]
    recorder=TrajectoryRecorder(args;save_fields=fields)
    meta=withenv("SPACEAGORA_RHS_EXECUTION_MODE"=>"serial","SPACEAGORA_RHS_CALIBRATE"=>"off",
                 "SPACEAGORA_RHS_IDENTIFY"=>"0","SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST"=>"0") do
        run_simulation(args;return_solver_metadata=true,visualization=false,
            extra_callbacks=(get_trajectory_recorder_callback(recorder),))
    end
    @test meta.retcode=="Success"
    times=collect(trajectory_times(recorder)); values=trajectory_save_data(recorder)
    @test length(times)>=20
    h0=SM.rot(q0)'*(inertia*w0+command_values(cmd,0.).momentum)
    for (t,snap) in zip(times,values)
        total=SM.rot(snap[:q][1])'*(inertia*snap[:omega][1]+command_values(cmd,t).momentum)
        @test norm(total-h0)<1e-7
    end
    @test norm(values[end][:omega][1]-w0)>1e-3
end

end
