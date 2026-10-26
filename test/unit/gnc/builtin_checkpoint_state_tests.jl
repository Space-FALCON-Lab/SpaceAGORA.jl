# The built-in guidance and control models' run state in checkpoint files
# (checkpoint_state / restore_checkpoint_state!), the burn scheduler's commands and
# plans and the orbit counters: a run resumed from a checkpoint ends as an
# uninterrupted run with the same checkpoints, bit for bit, in configurations where
# each model's state changes (the kit's
# `checkpoint_round_trip`, on runs it cannot make with its own fixture); the RPO
# planner binding resumes without planning again and refuses a planner whose own
# state a checkpoint cannot hold; the RPO MPC solver continues from its saved
# iterate and step size exactly.
module BuiltinCheckpointStateTests

using Test, Random, Serialization, StaticArrays, LinearAlgebra, SpaceAGORA
import SpaceAGORA.TelemetryVerification: make_example_config
const SM = SpaceAGORA.SimulationModel
const C = SpaceAGORA.Conformance
const L = SpaceAGORA.SimulationLifecycle

round_trip(model, kind, configuration) = C._gnc_checkpoint_round_trip(model, kind, configuration)
# The RPO MPC controller's run state but its solver workspace and last result, which
# a copy rebuilds (and whose result holds the solve's wall-clock times): the held
# command with its write record, the warm start, the mounted levels and the command
# log, which records every solve's command and status.
mpc_state(c) = (held=c.held, run_written=c.run_written, U_prev=c.controller.U_prev, mounted_levels=c.mounted_levels,
    command_log=c.command_log)
function mpc_round_trip(configuration)
    run = (cfg, factories) -> C._model_solve(cfg, :control_model, 1; factories)
    difference = function (a, b)
        gap = C._run_difference(a, b)
        gap === nothing || return gap
        return C._same_state(mpc_state(a.solution.model), mpc_state(b.solution.model)) ? nothing :
            "the controller's held command, warm start, mounted levels or command log differ"
    end
    return C._checkpoint_round_trip_check(run, difference, only(configuration.control_model.control_effectors), configuration)
end
# The data as a resumed run reads it from the file.
from_file(data) = (io = IOBuffer(); serialize(io, data); seekstart(io); deserialize(io))
const Z3 = SVector(0.0, 0.0, 0.0)

# A two-spacecraft RPO run through the planner binding and the MPC controller,
# with a replan after the halfway checkpoint.
rpo_configuration(planner=DirectRPOPlanner(segment_clearance_at=SpaceAGORA.RPOPlannerLifecycle.station_segment_clearance);
                  kwargs...) =
    make_rpo_configuration(; planner=planner, mission_time_s=4.0, planning_events=(RPOPlanningEvent(3.0; reason=:replan),),
        kwargs...)

# User planners holding the chaser at its goal: one that counts its proposals in
# its own state, with or without declaring restart, and one that fails an initial
# request made after the run's start (a resumed run must not plan again).
function hold_result(r::SpaceAGORA.RPOPlanningRequest, diagnostics)
    reference = SpaceAGORA.RPOReference(t_ref_s=[0.0, r.reference_dt_s],
        r_ref_rtn_m=hcat(collect(r.x_rtn[1:3]), collect(r.goal_rtn_m)), v_ref_rtn_mps=zeros(3, 2),
        origin_time_s=r.time_s, valid_until_s=r.valid_until_s, chaser_id=r.chaser_id, target_id=r.target_id,
        geometry_revision=r.geometry_revision)
    return SpaceAGORA.RPOPlanningResult(request_id=r.request_id, status=:candidate, termination=:completed,
        reference=reference, diagnostics=diagnostics)
end
struct CountingPlanner{R} <: SpaceAGORA.AbstractRPOPlanner end
SpaceAGORA.planner_capabilities(::CountingPlanner{R}) where {R} =
    SpaceAGORA.RPOPlannerCapabilities(state_sources=(:truth,), frames=(:target_rtn,), restart=R)
SpaceAGORA.initialize_planner(::CountingPlanner, context) = Ref(0)
SpaceAGORA.plan_rpo!(state, ::CountingPlanner, r::SpaceAGORA.RPOPlanningRequest, rng::Random.AbstractRNG) =
    (state[] += 1; hold_result(r, (calls=state[], draw=rand(rng))))
struct StartOnlyPlanner <: SpaceAGORA.AbstractRPOPlanner end
SpaceAGORA.planner_capabilities(::StartOnlyPlanner) =
    SpaceAGORA.RPOPlannerCapabilities(state_sources=(:truth,), frames=(:target_rtn,))
function SpaceAGORA.plan_rpo!(::Nothing, ::StartOnlyPlanner, r::SpaceAGORA.RPOPlanningRequest, rng::Random.AbstractRNG)
    r.reason === :initial && r.time_s > 0 && error("an initial request after the run's start")
    return hold_result(r, (draw=rand(rng),))
end

# Guidance models for the round trip's other outcomes: one that refuses checkpointed
# runs at setup, one whose run state is not in a checkpoint, and one that returns
# run state it cannot restore.
mutable struct RefusesCheckpoints <: SpaceAGORA.AbstractGuidanceModel
    ticks::Int
end
SpaceAGORA.update_guidance!(m::RefusesCheckpoints, u, p, t::Float64, i::Int) = (m.ticks += 1; nothing)
L.initialize_guidance!(m::RefusesCheckpoints, u, p, t) = (m.ticks = 0; nothing)
function L.preflight_guidance(m::RefusesCheckpoints, args; isolate_state)
    args.simulation_settings.checkpoint_enabled && throw(ArgumentError("RefusesCheckpoints does not run with checkpoints."))
    return nothing
end
mutable struct NoCheckpoint <: SpaceAGORA.AbstractGuidanceModel
    ticks::Int
end
SpaceAGORA.update_guidance!(m::NoCheckpoint, u, p, t::Float64, i::Int) = (m.ticks += 1; nothing)
L.initialize_guidance!(m::NoCheckpoint, u, p, t) = (m.ticks = 0; nothing)
# A model that fails during a checkpointed run, after its start.
mutable struct FailsCheckpointedTick <: SpaceAGORA.AbstractGuidanceModel
    ticks::Int
end
function SpaceAGORA.update_guidance!(m::FailsCheckpointedTick, u, p, t::Float64, i::Int)
    p.args.simulation_settings.checkpoint_enabled && t > 5.0 && throw(ArgumentError("a tick of a checkpointed run failed"))
    m.ticks += 1
    return nothing
end
L.initialize_guidance!(m::FailsCheckpointedTick, u, p, t) = (m.ticks = 0; nothing)
mutable struct CannotRestore <: SpaceAGORA.AbstractGuidanceModel
    ticks::Int
end
SpaceAGORA.update_guidance!(m::CannotRestore, u, p, t::Float64, i::Int) = (m.ticks += 1; nothing)
L.initialize_guidance!(m::CannotRestore, u, p, t) = (m.ticks = 0; nothing)
L.checkpoint_state(m::CannotRestore) = (ticks=m.ticks,)

function with_settings(configuration; kwargs...)
    st = configuration.simulation_settings
    fields = NamedTuple{fieldnames(typeof(st))}(map(n -> getfield(st, n), fieldnames(typeof(st))))
    return SM.SimConfig._with_configuration(configuration; simulation_settings=SM.SimulationSettings(; merge(fields, (; kwargs...))...))
end

# The ideal-actuation lunar descent of the Apollo unit tests, for 20 s.
function apollo_configuration()
    planet = SM.Moon()
    ephem = SM.SimpleEphemeridesModel()
    initial_time = SM.InitialTime(year=1969, month=7, day=20, hour=20, minute=10, second=0.0)
    lat, lon, h0 = 0.674, 23.473, 60.0
    et = SM.ephemerides_epoch_seconds(ephem, initial_time)
    l_pi = SM.ephemerides_frame_lpi(ephem, planet, et)
    upv = SVector(cosd(lat) * cosd(lon), cosd(lat) * sind(lon), sind(lat))
    r_p = (planet.Rp_e + h0) * upv
    r_i = SVector{3, Float64}(l_pi' * r_p)
    v_i = SVector{3, Float64}(l_pi' * cross(SVector{3, Float64}(planet.ω), r_p))
    up_i = normalize(r_i)
    tilt = SVector(0.0, sin(deg2rad(8.0) / 2), 0.0)
    fwd_i = normalize(SVector{3, Float64}(l_pi' * SVector(sind(lon), -cosd(lon), 0.0)))
    q0 = SM.quat_mult(descent_attitude_command(up_i, up_i, fwd_i), SVector(tilt[1], tilt[2], tilt[3], cos(deg2rad(8.0) / 2)))
    q0 = q0 / norm(q0)
    ic = SM.CartesianInitialCondition(r_i, v_i; q=SVector{4, Float64}(q0))
    root = SM.Link(root=true, m=6_900.0, dims=MVector{3, Float64}(4.2, 4.2, 4.0), ref_area=1.0, q=MVector{4, Float64}(q0...))
    sc = SM.SpacecraftModel(links=[root], root=root, prop_mass=1_000.0,
        inertia_tensor=SMatrix{3, 3, Float64}(2.2e4, 0, 0, 0, 2.4e4, 0, 0, 0, 2.0e4), initial_condition=ic, id=97)
    ph = DescentPhaseTargets(name=:braking, r_T=Z3, v_T=Z3, a_T=Z3, t_go_initial_s=10.0, t_go_min_s=1.0, altitude_switch_m=Inf)
    gcfg = ApolloDescentConfig(reference_radius_m=planet.Rp_e, site_lat_deg=lat, site_lon_deg=lon, site_height_m=0.0,
        approach_azimuth_deg=270.0, braking=ph, approach=ph, ignition_trim_s=0.0, engine_cutoff_altitude_m=0.0)
    state = ApolloDescentState(1)
    guidance = ApolloDescentGuidanceModel(gcfg, state, NoTerrainModel())
    control = ApolloDescentControlModel(ApolloDescentControlConfig(touchdown_height_m=0.5), gcfg, state, NoTerrainModel())
    base = make_example_config(planet=planet, spacecraft=sc, mission_time=20.0, initial_time=initial_time,
        dynamic_effectors=(SM.InverseSquaredGravityModel(),), density_model=SM.NoAtmosphereModel(), ephemerides_model=ephem,
        orientation_sim=true, keplerian=false, EI_km=1.0, verbose=false, results=false)
    args = SM.SimConfig._with_configuration(base;
        guidance_model=SM.GuidanceModel(guidance_effectors=(guidance,), guidance_rates=[1.0]),
        control_model=SM.ControlModel(control_effectors=(control,), control_rates=[0.05]),
        integration_tolerances=SM.IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9,
            dt_max_orbit=0.25, reltol_atmosphere=1e-9, abstol_atmosphere=1e-9, dt_max_atmosphere=0.25),
        solver_config=SM.SolverConfig(solver_mode=:tsit5))
    return args, guidance, control
end

# Three spacecraft on 400 km by 2000 km orbits with a maneuver guidance model and one
# burn scheduler for all of them, for 14000 s (about two orbits). The periapsis-raise
# guidance issues its commands before apoapsis, and the scheduler plans them at the
# next scheduling pass, the start of the segment after the halfway checkpoint: the
# commands and plans waiting there must survive the checkpoint. The campaign
# guidance, with an atmosphere interface above the periapsis (no atmosphere, so no
# drag) whose exits are scheduling passes, plans each orbit's maneuver from the
# orbit counter, which must survive it too.
function maneuver_configuration(guidance=nothing; EI=120.0, mission=(;))
    planet = SM.make_no_gram_planet(:earth)
    sats = map(1:3) do k
        root = SM.Link(root=true, m=100.0, ref_area=2.0)
        ic = SM.InitialCondition(ra=planet.Rp_e + 2000e3, rp=planet.Rp_e + 400e3, i=30.0, ω=0.0, Ω=40.0k, ν=10.0 + 30.0k)
        SM.SpacecraftModel(joints=SM.Joint[], links=[root], root=root, prop_mass=50.0, inertia_tensor=root.inertia,
            initial_condition=ic, id=k)
    end
    guidance = guidance !== nothing ? guidance : SM.ApoapsisTargetPeriapsisRaiseGuidanceModel(
        target_apoapsis_radius_m=planet.Rp_e + 3000e3, target_periapsis_altitude_m=500e3)
    scheduler = SM.ThrusterModels.BaseThrusterModel(thrust=fill(50.0, 3), direction=zeros(3), Δv=zeros(3),
        start_burn_time=fill(-1.0, 3), stop_burn_time=fill(-1.0, 3), Isp=fill(300.0, 3))
    env = SM.EnvironmentModel(planet=planet, EI=EI, density_model=SM.NoAtmosphereModel(),
        ephemerides_model=SM.SimpleEphemeridesModel(),
        thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet), topography=false, wind=false)
    args = SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=SM.MissionConfiguration(; mission_type=SM.MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=14000.0, orientation_sim=false, num_steps_to_save=1000, mission...),
        environment_model=env,
        dynamics_model=SM.DynamicsModel(SM.SpacecraftModel[sats...], (SM.InverseSquaredGravityModel(),)),
        guidance_model=SM.GuidanceModel((guidance,), [30.0]),
        control_model=SM.ControlModel(control_effectors=(scheduler,), control_rates=[1.0]),
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
    return args, guidance, scheduler
end

@testset "Built-in models' run state survives a checkpoint" begin
    @testset "a resumed run ends as an uninterrupted one" begin
        # The momentum manager with its command and field callables: it accumulates
        # wheel momentum and holds a rod torque that turns the spacecraft.
        manager = SM.ControlHooks.MagneticMomentumManagerModel(h_wheels_0=SVector(0.05, -0.02, 0.01),
            commanded_torque=(t, r, v, q, w) -> SVector(1.0e-5, -2.0e-5, 5.0e-6), b_field_ii=(t, r) -> SVector(2.0e-5, 1.0e-5, -3.0e-5))
        cases = Any[(:control_model, manager, C.conformance_fixture(control=(manager,)).configuration)]
        # A burn that is on at the halfway checkpoint.
        burn = SM.ThrusterModels.BaseThrusterModel(thrust=[0.5], direction=[0.0], Δv=[0.1], start_burn_time=[25.0],
            stop_burn_time=[40.0], Isp=[300.0])
        push!(cases, (:control_model, burn, C.conformance_fixture(control=(burn,)).configuration))
        # The RPO planner, which replans after the checkpoint (its controller below).
        rpo = rpo_configuration()
        push!(cases, (:guidance_model, only(rpo.guidance_model.guidance_effectors), rpo))
        # The Apollo descent guidance and control.
        apollo, guidance, control = apollo_configuration()
        push!(cases, (:guidance_model, guidance, apollo), (:control_model, control, apollo))
        # The periapsis-raise guidance and the burn scheduler, with commands and plans
        # waiting at the halfway checkpoint.
        raise, raise_guidance, raise_scheduler = maneuver_configuration()
        push!(cases, (:guidance_model, raise_guidance, raise), (:control_model, raise_scheduler, raise))
        # The campaign guidance and the burn scheduler, past the first orbit at the
        # halfway checkpoint.
        campaign, campaign_guidance, campaign_scheduler = maneuver_configuration(
            SM.AerobrakingCampaignPropulsiveManeuverGuidanceModel(maneuver_orbit_number=[1, 2], maneuver_Δv=[2.0, -1.0]);
            EI=500.0)
        push!(cases, (:guidance_model, campaign_guidance, campaign), (:control_model, campaign_scheduler, campaign))
        for (kind, model, configuration) in cases
            check = round_trip(model, kind, configuration)
            @testset let model = nameof(typeof(model)), status = check.status, detail = check.detail
                @test status === :pass
            end
        end
        check = mpc_round_trip(rpo)
        @testset let model = :RPOMPCControlModel, status = check.status, detail = check.detail
            @test status === :pass
        end
        # Their checkpoint holds the commands and plans only when one is set.
        buffers(p) = SpaceAGORA.SimulationEngine._maneuver_buffer_checkpoint_state(p)
        p = SM.ODEParams(n_sats=3, args=raise)
        @test buffers(p) === nothing
        p.shared_buffers.maneuver_commands[2] = SM.PropulsiveManeuverCommand(valid=true, delta_v_mps=1.5)
        @test buffers(p).commands[2].delta_v_mps == 1.5
        p.shared_buffers.maneuver_burn_plans[3] = SM.PropulsiveBurnPlan(valid=true, delta_v_mps=1.5, start_burn_s=8.0,
            stop_burn_s=9.0)
        entry = buffers(p)
        @test entry.version == 1 && entry.plans[3].valid && !entry.plans[1].valid
        q = SM.ODEParams(n_sats=3, args=raise)
        SpaceAGORA.SimulationEngine._restore_maneuver_buffer_checkpoint_state!(
            (runtime_state=(maneuver_buffers=from_file(entry),),), q)
        @test q.shared_buffers.maneuver_commands == p.shared_buffers.maneuver_commands
        @test q.shared_buffers.maneuver_burn_plans == p.shared_buffers.maneuver_burn_plans
        @test_throws "spacecraft count" SpaceAGORA.SimulationEngine._restore_maneuver_buffer_checkpoint_state!(
            (runtime_state=(maneuver_buffers=entry,),), SM.ODEParams(n_sats=2, args=raise))
        # The orbit counters, once some spacecraft is past its first orbit.
        orbits(p) = SpaceAGORA.SimulationEngine._orbit_counter_checkpoint_state(p)
        p = SM.ODEParams(n_sats=3, args=raise)
        @test orbits(p) === nothing
        p.orbit_counter[2] = 3
        q = SM.ODEParams(n_sats=3, args=raise)
        SpaceAGORA.SimulationEngine._restore_orbit_counter_checkpoint_state!((runtime_state=(orbits=from_file(orbits(p)),),), q)
        @test q.orbit_counter == [1, 3, 1]
        @test_throws "spacecraft count" SpaceAGORA.SimulationEngine._restore_orbit_counter_checkpoint_state!(
            (runtime_state=(orbits=orbits(p),),), SM.ODEParams(n_sats=2, args=raise))
        # A run that ends after two orbits, resumed after its first checkpoint, ends
        # where it ends without the break, not an orbit later.
        checkpointed(dir; resume=false) = (checkpoint_enabled=true, checkpoint_interval_s=4000.0,
            checkpoint_directory=dir, resume_from_checkpoint=resume)
        two_orbits(mission_time, settings) = with_settings(first(maneuver_configuration(nothing; EI=500.0,
            mission=(mission_type=SM.MissionOrbits, number_of_orbits=2, mission_time=mission_time))); settings...)
        whole = run_simulation(two_orbits(30000.0, checkpointed(mktempdir())); return_solution=true)
        split = mktempdir()
        run_simulation(two_orbits(4000.0, checkpointed(split)))
        resumed = run_simulation(two_orbits(30000.0, checkpointed(split; resume=true)); return_solution=true)
        @test whole.t[end] < 30000.0 && resumed.t[end] == whole.t[end]
        @test isequal(collect(resumed.u[end]), collect(whole.u[end]))
        # The same checkpoint without the new entries, as an earlier SpaceAGORA wrote
        # it, resumes with no error and counts from the first orbit, as it did then:
        # the run ends late.
        earlier = mktempdir()
        run_simulation(two_orbits(4000.0, checkpointed(earlier)))
        paths = SpaceAGORA.SimulationEngine._checkpoint_paths(two_orbits(4000.0, checkpointed(earlier)))
        payload = deserialize(paths.data)
        @test hasproperty(payload.runtime_state, :orbits)
        runtime = Base.structdiff(payload.runtime_state, NamedTuple{(:orbits, :maneuver_buffers)})
        payload = isempty(runtime) ? Base.structdiff(payload, NamedTuple{(:runtime_state,)}) :
            merge(payload, (runtime_state=runtime,))
        open(io -> serialize(io, payload), paths.data, "w")
        late = run_simulation(two_orbits(30000.0, checkpointed(earlier; resume=true)); return_solution=true)
        @test late.t[end] > whole.t[end]
        # The entry-end counter is owned by the run, not a fresh callback closure.
        entries(p) = SpaceAGORA.SimulationEngine._entry_counter_checkpoint_state(p)
        restore_entries!(entry, p) = SpaceAGORA.SimulationEngine._restore_entry_counter_checkpoint_state!(
            (runtime_state=(entries=entry,),), p)
        p = SM.ODEParams(n_sats=3, args=raise)
        @test entries(p) === nothing
        p.shared_buffers.entry_counter[2] = 1
        saved_entries = from_file(entries(p))
        p.shared_buffers.entry_counter[2] = 2
        q = SM.ODEParams(n_sats=3, args=raise)
        restore_entries!(saved_entries, q)
        @test q.shared_buffers.entry_counter == [0, 1, 0]
        @test_throws "spacecraft count" restore_entries!(saved_entries, SM.ODEParams(n_sats=2, args=raise))
        @test_throws "version 1" restore_entries!(merge(saved_entries, (version=2,)), q)
        @test_throws "nonnegative" restore_entries!(merge(saved_entries, (entry_counter=[0, -1, 0],)), q)
        @test_throws "spacecraft count" restore_entries!((version=1,), q)
        fresh = SM.ODEParams(n_sats=3, args=raise)
        SpaceAGORA.SimulationEngine._restore_entry_counter_checkpoint_state!((;), fresh)
        @test fresh.shared_buffers.entry_counter == [0, 0, 0]

        # Three unpowered spacecraft cross a 500 km interface on eccentric
        # orbits. Drag is present to enable the entry-end callback, but density
        # is zero. Resume after each has entered once, before the second entry.
        function two_entries(mission_time, settings)
            base = first(maneuver_configuration(nothing; EI=500.0, mission=(mission_time=mission_time,)))
            config = SM.SimConfig._with_configuration(base;
                dynamics_model=SM.DynamicsModel(base.dynamics_model.spacecraft,
                    (SM.InverseSquaredGravityModel(), SM.AerodynamicCoefficientConstant())),
                guidance_model=SM.GuidanceModel(), control_model=SM.ControlModel(),
                solver_config=SM.SolverConfig(solver_mode=:tsit5))
            return with_settings(config; settings...)
        end
        entry_settings(dir; resume=false) = (checkpoint_enabled=true, checkpoint_interval_s=7000.0,
            checkpoint_directory=dir, resume_from_checkpoint=resume)
        withenv("SPACEAGORA_ENTRY_TARGET_COUNT" => "2") do
            whole_entries = run_simulation(two_entries(30000.0, entry_settings(mktempdir())); return_solution=true)
            split_entries = mktempdir()
            part = two_entries(7000.0, entry_settings(split_entries))
            run_simulation(part)
            entry_paths = SpaceAGORA.SimulationEngine._checkpoint_paths(part)
            payload = deserialize(entry_paths.data)
            @test payload.runtime_state.entries.entry_counter == [1, 1, 1]
            # Retain the old-format checkpoint for a compatibility comparison.
            legacy_payload = merge(payload, (runtime_state=Base.structdiff(payload.runtime_state, NamedTuple{(:entries,)}),))
            continued = run_simulation(two_entries(30000.0, entry_settings(split_entries; resume=true)); return_solution=true)
            @test whole_entries.t[end] < 30000.0 && continued.t[end] == whole_entries.t[end]
            @test isequal(collect(continued.u[end]), collect(whole_entries.u[end]))
            # No count in old checkpoints: preserve the old zero-count behavior.
            # Run past the new checkpoint's second-entry end and inspect the
            # saved counters, rather than assuming a complete third orbit.
            serialize(entry_paths.data, legacy_payload)
            legacy_horizon = whole_entries.t[end] + 1.0
            legacy = run_simulation(two_entries(legacy_horizon, entry_settings(split_entries; resume=true)); return_solution=true)
            @test legacy.t[end] == legacy_horizon
            @test deserialize(entry_paths.data).runtime_state.entries.entry_counter == [1, 1, 1]
        end
        # A run's random stream compares by its state.
        @test C._same_state(MersenneTwister(3), from_file(MersenneTwister(3)))
        @test !C._same_state(MersenneTwister(3), MersenneTwister(4))
    end

    @testset "the round trip skips a refusal and fails a model whose state does not survive" begin
        # Added in 0.4.0: advisory in 0.4, required from 0.5.0.
        @test C.CHECKPOINT_ROUND_TRIP_SINCE == v"0.4.0"
        failing = C._failure_status(v"0.4.0")
        @test C._failure_status(v"0.4.0", v"0.4.3") === :advisory_fail && C._failure_status(v"0.4.0", v"0.5.0") === :fail
        refusing = RefusesCheckpoints(0)
        check = round_trip(refusing, :guidance_model, C.conformance_fixture(guidance=(refusing,)).configuration)
        @test check.status === :skip && occursin("refuses to run with checkpoints", check.detail)
        # An ArgumentError once the solve has started is a failure, not a refusal.
        failing_tick = FailsCheckpointedTick(0)
        check = round_trip(failing_tick, :guidance_model, C.conformance_fixture(guidance=(failing_tick,)).configuration)
        @test check.status === failing && occursin("the fixture failed with checkpoints", check.detail) &&
            occursin("a tick of a checkpointed run failed", check.detail)
        forgetful = NoCheckpoint(0)
        check = round_trip(forgetful, :guidance_model, C.conformance_fixture(guidance=(forgetful,)).configuration)
        @test check.status === failing && check.since == v"0.4.0" &&
            occursin("and an uninterrupted run differ after it: the model's own state", check.detail)
        unrestorable = CannotRestore(0)
        check = round_trip(unrestorable, :guidance_model, C.conformance_fixture(guidance=(unrestorable,)).configuration)
        @test check.status === failing && occursin("has no restore_checkpoint_state! method", check.detail)
    end

    @testset "the RPO planner binding" begin
        # A resumed run takes the binding's state from the file: its initializer makes
        # no initial request when the engine says a restore follows.
        args = rpo_configuration()
        g = deepcopy(only(args.guidance_model.guidance_effectors))
        L.initialize_guidance!(g, nothing, (args=args, shared_buffers=(model_resume=Ref(true),)), 0.0)
        @test g.runtime.request_id == 0 && isempty(g.runtime.records) && g.runtime.product === nothing
        @test g.chaser_index == 1 && g.target_index == 2 && !g.plan_buffer.valid
        # In a run: the resumed run's initializer would fail an initial request.
        start_only = rpo_configuration(StartOnlyPlanner(); goal_rtn_m=(3.0, 0.0, 0.0))
        check = round_trip(only(start_only.guidance_model.guidance_effectors), :guidance_model, start_only)
        @test check.status === :pass
        # A planner with state of its own runs with checkpoints when it declares
        # restart, and its state then resumes with the run.
        restartable = rpo_configuration(CountingPlanner{true}(); goal_rtn_m=(3.0, 0.0, 0.0))
        check = round_trip(only(restartable.guidance_model.guidance_effectors), :guidance_model, restartable)
        @test check.status === :pass
        refused = with_settings(rpo_configuration(CountingPlanner{false}(); goal_rtn_m=(3.0, 0.0, 0.0));
            checkpoint_enabled=true, checkpoint_directory=mktempdir())
        err = try
            run_simulation(refused)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin("does not declare restart", sprint(showerror, err))
        # Without checkpoints it runs.
        @test string(run_simulation(rpo_configuration(CountingPlanner{false}(); goal_rtn_m=(3.0, 0.0, 0.0));
            return_solution=true).retcode) == "Success"
    end

    @testset "the RPO MPC controller continues from its saved warm start, iterate and step size" begin
        model = only(rpo_configuration().control_model.control_effectors)
        controller = model.controller
        reference = zeros(6, controller.horizon + 1)
        reference[1, :] .= 5.0
        state(k) = SVector(3.0 + 0.1k, 0.01k, 0.0, 0.002, 0.0, 0.0)
        for k in 1:30
            SM.ControlHooks.rpo_lqmpc_control(controller, state(k), reference)
        end
        data = from_file(L.checkpoint_state(model))
        # A copy has a fresh workspace, whose iterate and step size start over.
        fresh, unrestored = deepcopy(model), deepcopy(model)
        L.restore_checkpoint_state!(fresh, data)
        copyto!(unrestored.controller.U_prev, data.controller.U_prev)
        @test SM.ControlHooks._rpo_lqmpc_solver_state(fresh.controller.qp_model) == data.controller.osqp
        commands(c) = [SM.ControlHooks.rpo_lqmpc_control(c, state(k), reference) for k in 31:40]
        expected = commands(controller)
        @test commands(fresh.controller) == expected                  # bit for bit
        @test commands(unrestored.controller) != expected             # the warm start alone does not repeat them
        @test_throws ArgumentError SM.ControlHooks._restore_rpo_lqmpc_solver_state!(fresh.controller.qp_model,
            (y=zeros(2), rho=0.1))
        # A checkpoint of a model without a controller is refused.
        @test_throws ArgumentError L.restore_checkpoint_state!(fresh,
            from_file(L.checkpoint_state(SM.ControlHooks.RPOMPCControlModel())))
    end

    @testset "each built-in model's data restores into a fresh copy" begin
        # The momentum manager, with and without a write of its state recorded (G-8).
        manager() = SM.ControlHooks.MagneticMomentumManagerModel(h_wheels_0=SVector(0.1, 0.0, 0.0))
        written = manager()
        written.h_wheels, written.held_torque_body, written.held_dipole_am2 = SVector(1.0, 2.0, 3.0), SVector(0.1, 0.0, 0.0),
            SVector(0.0, 2.0, 0.0)
        written.ticks, written.last_update_s, written.initialized = 7, 6.5, true
        restored = manager()
        L.restore_checkpoint_state!(restored, from_file(L.checkpoint_state(written)))
        @test C._same_state(restored, written) && restored.run_written === nothing
        SM.ControlHooks._note_manager_writes!(written)
        L.restore_checkpoint_state!(restored, from_file(L.checkpoint_state(written)))
        @test C._same_state(restored, written)
        # The robot-arm controller: state a run wrote, saved, read back and put into
        # a fresh copy after its reset.
        arm = SM.default_robot_arm_model()
        model = SM.RobotArmControlEffector(controller=SM.init_robot_arm_joint_mpc(arm; dt_s=0.1, horizon=4))
        model.held = SM.RobotArmHeldActuation(joint_torque_nm=[0.5, -0.25, 0.0])
        model.run_written = deepcopy(model.held)
        model.controller.U_prev .= 1.5
        restored_arm = SM.RobotArmControlEffector(controller=SM.init_robot_arm_joint_mpc(arm; dt_s=0.1, horizon=4))
        L.restore_checkpoint_state!(restored_arm, from_file(L.checkpoint_state(model)))
        @test C._same_state(restored_arm, model)
        # A checkpoint from another controller is refused.
        plain = SM.RobotArmControlEffector()
        @test_throws ArgumentError L.restore_checkpoint_state!(plain, from_file(L.checkpoint_state(model)))
        # The periapsis-raise guidance and the burn scheduler.
        raise = SpaceAGORA.ApoapsisTargetPeriapsisRaiseGuidanceModel(target_apoapsis_radius_m=7.0e6, command_state=[2, 0])
        other = SpaceAGORA.ApoapsisTargetPeriapsisRaiseGuidanceModel(target_apoapsis_radius_m=7.0e6)
        L.restore_checkpoint_state!(other, from_file(L.checkpoint_state(raise)))
        @test other.command_state == [2, 0]
        two(Δv) = SM.ThrusterModels.BaseThrusterModel(thrust=[50.0, 50.0], direction=[0.0, 0.0], Δv=Δv,
            start_burn_time=[5.0, 6.0], stop_burn_time=[20.0, 21.0], Isp=[300.0, 300.0])
        scheduler = two([0.1, 0.2])
        scheduler.start_burn_time[2], scheduler.burning[2] = 7.5, true
        SM.ControlHooks._note_schedule!(scheduler, 2)
        restored = two([0.3, 0.4])
        L.restore_checkpoint_state!(restored, from_file(L.checkpoint_state(scheduler)))
        @test restored.Δv == [0.1, 0.2] && restored.start_burn_time == [5.0, 7.5] && restored.burning == [false, true]
        @test restored.written_schedule == scheduler.written_schedule
        @test_throws ArgumentError L.restore_checkpoint_state!(SM.ThrusterModels.BaseThrusterModel(thrust=[50.0],
            direction=[0.0], Δv=[0.1], start_burn_time=[5.0], stop_burn_time=[20.0], Isp=[300.0]), from_file(L.checkpoint_state(scheduler)))
    end
end

end
