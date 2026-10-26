module CampaignRequestStabilityTests
using Test, SpaceAGORA, CSV, DataFrames
using OrdinaryDiffEq: DiscreteCallback
const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const L = SpaceAGORA.SimulationLifecycle
const CH = SM.ControlHooks
const GH = SM.GuidanceHooks
const PLANET = make_no_gram_planet(:earth)

function configuration(; duration=18000.0, drag=true, guidance_period=1.0, control_period=1.0,
        directory="", resume=false, interval=2200.0, flight_scaling=true)
    sats = [(r = Link(root=true, m=100.0, ref_area=2.0);
        SpacecraftModel(joints=Joint[], links=[r], root=r, prop_mass=50.0,
            inertia_tensor=r.inertia, id=k, initial_condition=InitialCondition(
                ra=PLANET.Rp_e+2000e3, rp=PLANET.Rp_e+(drag ? 140e3 : 400e3),
                i=60.0, ω=20.0k, Ω=40.0k, ν=10.0+30.0k))) for k in 1:3]
    campaign = SM.AerobrakingCampaignPropulsiveManeuverGuidanceModel(
        maneuver_orbit_number=[1,2,3], maneuver_Δv=[2.0,2.0,-1.0],
        maneuver_flight_apoapsis_radius_m=flight_scaling ? fill(PLANET.Rp_e+2010e3,3) : Float64[])
    scheduler = BaseThrusterModel(thrust=fill(2000.0,3), direction=zeros(3), Δv=zeros(3),
        start_burn_time=fill(-1.0,3), stop_burn_time=fill(-1.0,3), Isp=fill(300.0,3))
    density = drag ? SM.ExponentialAtmosphereModel(1e-9, 1.4e5, 7e3;
        temperature_k=800.0, valid_min_altitude_m=0.0, valid_max_altitude_m=2e6) : NoAtmosphereModel()
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false,
            normalize=false, checkpoint_enabled=!isempty(directory), checkpoint_directory=directory,
            checkpoint_interval_s=interval, resume_from_checkpoint=resume),
        mission_configuration=MissionConfiguration(mission_type=MissionTime, keplerian=true,
            number_of_orbits=1, mission_time=duration, orientation_sim=false, num_steps_to_save=1000),
        environment_model=EnvironmentModel(planet=PLANET, EI=200.0, density_model=density,
            ephemerides_model=SimpleEphemeridesModel(), topography=false, wind=false,
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=PLANET)),
        dynamics_model=DynamicsModel(SpacecraftModel[sats...],
            (InverseSquaredJ2GravityModel(), SM.AerodynamicCoefficientConstant())),
        guidance_model=GuidanceModel((campaign,), [guidance_period]),
        navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((scheduler,), [control_period]),
        initial_time=InitialTime(year=2020,month=1,day=1,hour=0,minute=0,second=0.0))
end
function traced(args; path=tempname()*".csv", extra_callbacks=())
    sol=withenv("SPACEAGORA_MANEUVER_TRACE_CSV"=>path) do
        run_simulation(args; return_solution=true, extra_callbacks)
    end
    return sol, CSV.read(path,DataFrame)
end
function check_campaign(sol,trace)
    plans=trace[trace.event .== "schedule_set",:]
    starts=trace[trace.event .== "burn_start",:]
    stops=trace[trace.event .== "burn_end",:]
    for i in 1:3
        p=plans[plans.spacecraft_idx .== i,:]
        @test p.orbit_counter == [1,2,3]
        @test count(==(i),starts.spacecraft_idx) == 3
        @test count(==(i),stops.spacecraft_idx) == 3
        @test sol.u[end].sc[i].mass ≈ 150*exp(-sum(p.dv_cmd_mps)/(300*9.80665)) rtol=1e-9
    end
end

@testset "Flight-scaled campaign consumes each maneuver once" begin
    for (drag,guidance_period,control_period) in ((true,1.0,1.0),(false,30.0,10.0))
        sol,trace=traced(configuration(;drag,guidance_period,control_period))
        check_campaign(sol,trace)
    end
    @testset "completed short burn survives checkpoint and resume" begin
        mktempdir() do whole_dir
            whole,whole_trace=traced(configuration(directory=whole_dir))
            mktempdir() do directory
                path=tempname()*".csv"
                stop=DiscreteCallback((u,t,integrator)->t>2205.0,
                    integrator->error("campaign interruption");save_positions=(false,false))
                @test_throws Exception traced(configuration(;directory);path,extra_callbacks=(stop,))
                resumed,trace=traced(configuration(;directory,resume=true);path)
                check_campaign(resumed,trace)
                @test isequal(resumed.u[end],whole.u[end])
                @test isequal(trace,whole_trace)
            end
        end
    end
end

@testset "Campaign latch preserves request and checkpoint contracts" begin
    args=configuration(duration=100.0)
    p=SM.ODEParams(n_sats=3,args=args)
    u=SE.build_initial_conditions(args)
    campaign=only(args.guidance_model.guidance_effectors)
    scheduler=only(args.control_model.control_effectors)
    L.initialize_control!(scheduler,u,p,0.0)
    p.orbit_counter .= 1
    update_guidance!(campaign,u,p,1.0,2)
    before=p.shared_buffers.maneuver_commands[2]
    # Before acceptance, state changes must still update the flight scaling.
    u.sc[2].vel .*= 1.0001
    update_guidance!(campaign,u,p,1.5,2)
    refreshed=p.shared_buffers.maneuver_commands[2]
    @test refreshed.delta_v_mps != before.delta_v_mps
    before=refreshed
    update_control!(scheduler,u,p,2.0,2)
    accepted=scheduler.consumed_command[2]
    @test accepted !== nothing
    @test accepted == (before.delta_v_mps,before.direction_rad,1)
    snapshot=L.checkpoint_state(scheduler)
    @test snapshot.version == 2
    @test L.checkpoint_state(campaign) === nothing
    # Restore a completed request with no live plan, as in an existing checkpoint.
    L.initialize_control!(scheduler,u,p,0.0)
    L.restore_checkpoint_state!(scheduler,snapshot)
    p.shared_buffers.maneuver_burn_plans[2]=PropulsiveBurnPlan()
    u.sc[2].vel .*= 1.0001
    bound=L.bind_spacecraft(campaign,2)
    update_guidance!(bound,u,p,2.0,2)
    after=p.shared_buffers.maneuver_commands[2]
    @test (after.delta_v_mps,after.direction_rad,after.source_orbit) == accepted
    first=p.shared_buffers.maneuver_commands[1]
    update_guidance!(bound,u,p,2.0,1)
    @test p.shared_buffers.maneuver_commands[1] === first
    p.orbit_counter[2]=2
    update_guidance!(bound,u,p,3.0,2)
    @test p.shared_buffers.maneuver_commands[2].source_orbit == 2
    # A custom producer may deliberately change its request on the same orbit.
    p.shared_buffers.maneuver_commands[2]=SM.PropulsiveManeuverCommand(valid=true,
        delta_v_mps=accepted[1]+1.0,direction_rad=accepted[2],source_orbit=1)
    @test CH._commanded_maneuver(scheduler,p,2).delta_v_mps == accepted[1]+1.0
    L.initialize_control!(scheduler,u,p,0.0)
    @test GH._accepted_campaign_request(scheduler,2,1) === nothing
    @test GH._accepted_campaign_request(nothing,2,1) === nothing
    # The first accepting scheduler freezes the shared campaign request for a
    # slower scheduler too; changing guidance values must not split the command.
    slow=deepcopy(scheduler)
    two=SM.SimConfig._with_configuration(args; control_model=ControlModel((scheduler,slow),[1.0,10.0]))
    p2=SM.ODEParams(n_sats=3,args=two)
    p2.orbit_counter .= 1
    L.initialize_control!(slow,u,p2,0.0)
    update_guidance!(campaign,u,p2,1.0,2)
    update_control!(scheduler,u,p2,1.0,2)
    first=scheduler.consumed_command[2]
    @test first !== nothing
    @test slow.consumed_command[2] === nothing
    u.sc[2].vel .*= 1.0001
    update_guidance!(campaign,u,p2,10.0,2)
    update_control!(slow,u,p2,10.0,2)
    @test slow.consumed_command[2] == first
end

@testset "Idle maneuver ticks do not save duplicate solver states" begin
    args=configuration(duration=100.0,drag=false)
    args=SM.SimConfig._with_configuration(args;guidance_model=GuidanceModel((),Float64[]))
    sol=run_simulation(args;return_solution=true)
    @test length(sol.t) == length(unique(sol.t))
    @test sol.t[end] == 100.0
    @test all(i->sol.u[end].sc[i].mass==150.0,1:3)
end
end
