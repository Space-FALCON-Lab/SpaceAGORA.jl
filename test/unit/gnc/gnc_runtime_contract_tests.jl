# The GNC runtime contract (GNC design 3.8): every built-in guidance, navigation and
# control model survives a deep copy (R-S3), gives each copy `bind_spacecraft` makes
# its own run state (R-S4), and resets its run state when a run starts, so that two
# runs of one configuration with isolate_state=false start equal (R-S2, G-8).
#
# The models that run in the contributor kit's fixture get the kit's checks in
# test/unit/conformance/builtin_conformance_tests.jl: `deepcopy_equivalent`,
# `initialize_resets` and `bind_fresh_state` (the kit's walker). This file covers
# the other built-in models with the same walker, and the engine's reset of
# navigation models, of custom sensor state and of every built-in sensor kind.
module GNCRuntimeContractTests

using Test, Random, StaticArrays, LinearAlgebra, SpaceAGORA
import SpaceAGORA.TelemetryVerification: make_example_config
const SM = SpaceAGORA.SimulationModel
const C = SpaceAGORA.Conformance
const L = SpaceAGORA.SimulationLifecycle

include(joinpath(@__DIR__, "..", "..", "fixtures", "conformance_instances.jl"))
using .ConformanceInstances: gnc_instances, gnc_exemptions, sensor_instances

const PLANET = SpaceAGORA.make_no_gram_planet(:earth)
const Z3 = SVector(0.0, 0.0, 0.0)

apollo_config() = SM.ApolloDescentConfig(reference_radius_m=PLANET.Rp_e - 2_000.0, site_lat_deg=30.0, site_lon_deg=0.0,
    braking=SM.apollo11_descent_targets()[1], approach=SM.apollo11_descent_targets()[2])

# A planner run's configuration: the planner guidance and an RPO MPC controller with
# its LQ-MPC controller (an OSQP workspace).
rpo_configuration() = make_rpo_configuration(;
    planner=DirectRPOPlanner(segment_clearance_at=SpaceAGORA.RPOPlannerLifecycle.station_segment_clearance),
    mission_time_s=2.0)

# The built-in models that cannot run in the kit's fixture, one per type, with how
# `bind_spacecraft` treats them: `:indexed` binds copies that share only the state
# the model keeps per spacecraft; `:refused` is the scope fallback's refusal (the
# model is declared at configuration level).
function contract_models()
    edg = SM.AerobrakingEnergyDepletionConfig()
    edg_state = SM.AerobrakingEnergyDepletionState(num_sats=1)
    apollo_state = SM.ApolloDescentState(2)
    rpo = rpo_configuration()
    burn = AlgorithmGuidance(algorithm=ApoapsisBurnPlanner(delta_v_mps=10.0), spacecraft_index=1,
        events=[AlgorithmEvent(0.0)])
    joint = AlgorithmGuidance(algorithm=QuinticJointPlanner(), spacecraft_index=1, link=RobotArmPlanBuffer(),
        settings=JointReferenceSettings(model=SM.default_robot_arm_model(), base_pose=SM.RobotArmBasePose([0.0, 0.0, 0.0]),
            target=(0.4, 0.1, 0.1)))
    return [
        (model=SM.ApolloDescentGuidanceModel(apollo_config(), apollo_state; spacecraft_indices=(1,)), bind=:indexed,
            shared=m -> (m.state,)),
        (model=SM.ApolloDescentControlModel(SM.ApolloDescentControlConfig(), apollo_config(), apollo_state;
            spacecraft_indices=(1,)), bind=:indexed, shared=m -> (m.state, m.actuators)),
        (model=SM.AerobrakingEnergyDepletionGuidanceModel(edg, edg_state), bind=:refused, shared=nothing),
        (model=SM.AerobrakingEnergyDepletionControlModel(edg, edg_state), bind=:refused, shared=nothing),
        (model=SM.GuidanceModels.RPOGuidanceModel(), bind=:refused, shared=nothing),
        (model=only(rpo.guidance_model.guidance_effectors), bind=:refused, shared=nothing),
        (model=only(rpo.control_model.control_effectors), bind=:refused, shared=nothing),
        (model=SM.SolarPanelAngleOfAttackControlModel(), bind=:refused, shared=nothing),
        (model=burn, bind=:run, shared=nothing),
        (model=joint, bind=:run, shared=nothing),
    ]
end

# A burn-plan run with one slot per spacecraft and the binding in the guidance model.
function burn_configuration(binding; spacecraft_count=1)
    n = spacecraft_count
    scheduler = SM.ThrusterModels.BaseThrusterModel(thrust=fill(50.0, n), direction=zeros(n), Δv=zeros(n), start_burn_time=fill(-1.0, n),
        stop_burn_time=fill(-1.0, n), Isp=fill(300.0, n))
    spacecraft = map(1:n) do id
        root = Link(root=true, m=100.0, ref_area=2.0)
        ic = InitialCondition(ra=PLANET.Rp_e + 2000e3, rp=PLANET.Rp_e + 400e3, i=30.0, ω=0.0, Ω=0.0, ν=10.0 + 5.0 * id)
        SpacecraftModel(joints=Joint[], links=[root], root=root, prop_mass=50.0, inertia_tensor=root.inertia,
            initial_condition=ic, id=id)
    end
    env = EnvironmentModel(planet=PLANET, EI=120.0, density_model=SM.NoAtmosphereModel(),
        ephemerides_model=SM.SimpleEphemeridesModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=PLANET), topography=false, wind=false)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=20.0, orientation_sim=false, num_steps_to_save=50),
        environment_model=env,
        dynamics_model=DynamicsModel(SpacecraftModel[spacecraft...], (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel((binding,), [1.0e6]),
        navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel(control_effectors=(scheduler,), control_rates=[1.0]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
end

# A short lunar descent of one spacecraft from 60 m, as in apollo_descent_tests.jl.
function apollo_descent_run()
    planet = SM.Moon()
    ephem = SM.SimpleEphemeridesModel()
    initial_time = SM.InitialTime(year=1969, month=7, day=20, hour=20, minute=10, second=0.0)
    lat, lon, h0 = 0.674, 23.473, 60.0
    et = SM.ephemerides_epoch_seconds(ephem, initial_time)
    l_pi = SM.ephemerides_frame_lpi(ephem, planet, et)
    up = SVector(cosd(lat) * cosd(lon), cosd(lat) * sind(lon), sind(lat))
    r_p = (planet.Rp_e + h0) * up
    r_i = SVector{3, Float64}(l_pi' * r_p)
    v_i = SVector{3, Float64}(l_pi' * cross(SVector{3, Float64}(planet.ω), r_p))
    up_i = normalize(r_i)
    fwd_i = normalize(SVector{3, Float64}(l_pi' * SVector(sind(lon), -cosd(lon), 0.0)))
    q0 = SM.GuidanceHooks.descent_attitude_command(up_i, up_i, fwd_i)
    root = SM.Link(root=true, m=6_900.0, dims=MVector{3, Float64}(4.2, 4.2, 4.0), ref_area=1.0, q=MVector{4, Float64}(q0...))
    sc = SM.SpacecraftModel(links=[root], root=root, prop_mass=1_000.0,
        inertia_tensor=SMatrix{3, 3, Float64}(2.2e4, 0, 0, 0, 2.4e4, 0, 0, 0, 2.0e4),
        initial_condition=SM.CartesianInitialCondition(r_i, v_i; q=SVector{4, Float64}(q0)), id=97)
    phase = SM.DescentPhaseTargets(name=:braking, r_T=Z3, v_T=Z3, a_T=Z3, t_go_initial_s=10.0, t_go_min_s=1.0,
        altitude_switch_m=Inf)
    config = SM.ApolloDescentConfig(reference_radius_m=planet.Rp_e, site_lat_deg=lat, site_lon_deg=lon, site_height_m=0.0,
        approach_azimuth_deg=270.0, braking=phase, approach=phase, ignition_trim_s=0.0, engine_cutoff_altitude_m=0.0)
    state = SM.ApolloDescentState(1)
    guidance = SM.ApolloDescentGuidanceModel(config, state, SM.NoTerrainModel())
    control = SM.ApolloDescentControlModel(SM.ApolloDescentControlConfig(touchdown_height_m=0.5), config, state, SM.NoTerrainModel())
    base = make_example_config(planet=planet, spacecraft=sc, mission_time=3.0, initial_time=initial_time,
        dynamic_effectors=(SM.InverseSquaredGravityModel(),), density_model=SM.NoAtmosphereModel(), ephemerides_model=ephem,
        orientation_sim=true, keplerian=false, EI_km=1.0, verbose=false, results=false)
    return SM.SimConfig._with_configuration(base;
        guidance_model=SM.GuidanceModel(guidance_effectors=(guidance,), guidance_rates=[1.0]),
        control_model=SM.ControlModel(control_effectors=(control,), control_rates=[0.05]),
        solver_config=SM.SolverConfig(solver_mode=:tsit5))
end

const ABSTRACT_TYPES = (guidance_model=SpaceAGORA.AbstractGuidanceModel, navigation_model=SpaceAGORA.AbstractNavigationModel,
    control_model=SpaceAGORA.AbstractControlEffectorModel)

@testset "Every built-in GNC model has the runtime contract checked" begin
    runnable = [entry.model for entry in gnc_instances() if entry.runs]
    others = [entry.model for entry in contract_models()]
    for (kind, T) in pairs(ABSTRACT_TYPES), S in C.package_types(T)
        @testset let kind = kind, type = S
            # In the kit's runs (deepcopy_equivalent, initialize_resets, bind_fresh_state) or here.
            @test count(m -> m isa S, runnable) + count(m -> m isa S, others) >= 1
        end
    end
    # The burn scheduler runs in the control model without subtyping the control effector type.
    @test any(m -> m isa SM.ThrusterModels.BaseThrusterModel, runnable)
    # The kit's reports of the runnable models carry the three checks.
    @test all(id -> id in C.CHECK_IDS.guidance_model, (:deepcopy_equivalent, :initialize_resets, :bind_fresh_state))
    @test all(id -> id in C.CHECK_IDS.control_model, (:deepcopy_equivalent, :initialize_resets, :bind_fresh_state))
end

@testset "Deep copies share no mutable state (R-S3)" begin
    models = [[entry.model for entry in gnc_instances()]; [entry.model for entry in contract_models()]]
    for model in models
        copied = deepcopy(model)
        @testset let type = typeof(model)
            @test isempty(C._shared_mutables(model, copied))
            # A native solver workspace is rebuilt, so the RPO controller is compared by its command.
            model isa SM.RPOMPCControlModel || @test C._same_state(model, copied)
        end
    end
    controller = only(rpo_configuration().control_model.control_effectors).controller
    copied = deepcopy(controller)
    @test copied.qp_model !== controller.qp_model
    x, ref = SVector(1.0, 0.5, -0.2, 0.0, 0.01, 0.0), zeros(6, controller.horizon + 1)
    @test SM.ControlHooks.rpo_lqmpc_control(copied, x, ref) == SM.ControlHooks.rpo_lqmpc_control(controller, x, ref)
end

@testset "Bound copies get their own run state (R-S4)" begin
    bind = SpaceAGORA.bind_spacecraft
    for entry in contract_models()
        model = entry.model
        @testset let type = typeof(model), rule = entry.bind
            if entry.bind === :refused
                err = try bind(model, 2) catch e e end
                @test err isa ArgumentError && occursin("Declare it at configuration level", sprint(showerror, err))
            elseif entry.bind === :indexed
                # The copies for spacecraft 1 and 2 share only the state the model keeps
                # per spacecraft, and each acts on its own spacecraft.
                one, two = bind(model, 1), bind(model, 2)
                declared = IdDict{Any, Nothing}()
                foreach(x -> C._reachable_mutables(x, declared), entry.shared(model))
                @test all(x -> haskey(declared, x), C._shared_mutables(one, two))
                @test one.spacecraft_indices == (1,) && two.spacecraft_indices == (2,)
            end
        end
    end
    # The generic binding of a family that finds its own install target runs bound on
    # two spacecraft: the kit's walker over a run with one copy per spacecraft.
    burn = only(entry.model for entry in contract_models()
        if entry.model isa AlgorithmGuidance && entry.model.algorithm isa ApoapsisBurnPlanner)
    two = burn_configuration(burn; spacecraft_count=2)
    report = C._bind_fresh_state_check(burn, :guidance_model, C.FIRST_RELEASE; configuration=two, index=1)
    @test report.status === :pass
    # A binding with a configured link shares it, so a second bound copy is refused
    # before the run: one plan buffer has one producer.
    joint = only(entry.model for entry in contract_models()
        if entry.model isa AlgorithmGuidance && entry.model.algorithm isa QuinticJointPlanner)
    bound = (bind(joint, 1), bind(joint, 2))
    @test bound[1].link === bound[2].link === joint.link
    @test bound[1].runtime === nothing && bound[2].runtime === nothing
end

# ---- initialize_*! resets run state (G-8) -------------------------------------

initialize!(m::SpaceAGORA.AbstractGuidanceModel) = L.initialize_guidance!(m, nothing, nothing, 0.0)
initialize!(m) = L.initialize_control!(m, nothing, nothing, 0.0)

# `model` at a run's start, after a run that `run!` stands for, equals `model` at the
# start of its first run.
function starts_equal(model, run!)
    initialize!(model)
    first_start = deepcopy(model)
    run!(model)
    initialize!(model)
    return C._same_state(model, first_start)
end

@testset "initialize_*! resets the run state of each built-in model (G-8)" begin
    @testset "Apollo descent guidance and control: their spacecraft's entries" begin
        state = SM.ApolloDescentState(2)
        guidance = SM.ApolloDescentGuidanceModel(apollo_config(), state; spacecraft_indices=(1,))
        control = SM.ApolloDescentControlModel(SM.ApolloDescentControlConfig(), apollo_config(), state; spacecraft_indices=(1,))
        fresh = SM.ApolloDescentState(2)
        # A finished descent of spacecraft 1, and spacecraft 2's entries, which belong to its own models.
        state.phase .= :landed
        state.phase_start_s .= 3.0
        state.touchdown_s .= 120.0
        state.site_frame[1] = (site=Z3,)
        control.actuators.thrust_n .= 900.0
        control.actuators.torque_nm[1] = SVector(1.0, 0.0, 0.0)
        push!(control.actuators.thruster_level[1], 0.5)
        initialize!(guidance)
        initialize!(control)
        for name in fieldnames(SM.ApolloDescentState)
            data, expected = getfield(state, name), getfield(fresh, name)
            @testset let field = name
                @test name === :phase_start_s ? isequal(data[:, 1], expected[:, 1]) : isequal(data[1], expected[1])
            end
        end
        @test state.phase[2] === :landed && state.touchdown_s[2] == 120.0
        @test control.actuators.thrust_n == [0.0, 900.0] && control.actuators.torque_nm[1] == Z3
        @test isempty(control.actuators.thruster_level[1]) && isnan(control.actuators.last_update_s[1])
    end
    @testset "Apollo descent: two runs on one configuration end equal" begin
        args = apollo_descent_run()
        first_run = run_simulation(args; isolate_state=false, return_solution=true)
        guidance = only(args.guidance_model.guidance_effectors)
        control = only(args.control_model.control_effectors)
        ended = deepcopy((guidance.state, control.actuators))
        @test isfinite(guidance.state.phase_start_s[1, 1])
        second_run = run_simulation(args; isolate_state=false, return_solution=true)
        @test first_run.t == second_run.t && all(isequal(collect(a), collect(b)) for (a, b) in zip(first_run.u, second_run.u))
        @test C._same_state((guidance.state, control.actuators), ended)
    end
    @testset "Energy-depletion guidance and control: the shared state, through its link" begin
        edg = SM.AerobrakingEnergyDepletionConfig()
        state = SM.AerobrakingEnergyDepletionState(num_sats=1)
        guidance = SM.AerobrakingEnergyDepletionGuidanceModel(edg, state)
        control = SM.AerobrakingEnergyDepletionControlModel(edg, state)
        link! = s -> SpaceAGORA.link_initialize!(s, (guidance,), (control,), nothing, nothing, 0.0)
        link!(state)
        first_start = deepcopy(state)
        state.control_clock_origin_s = 50.0
        state.energy_bracketing_count[1] = 4
        state.last_alpha_rad[1] = 0.3
        push!(state.heat_load_entry_j_cm2[1], 12.0)
        link!(state)
        @test C._same_state(state, first_start)
    end
    @testset "The deprecated RPO guidance: goal, plan buffer and replanning record" begin
        model = SM.GuidanceModels.RPOGuidanceModel(goal_rtn=SVector(5.0, 0.0, 0.0))
        buffer = model.plan_buffer
        @test starts_equal(model, m -> begin
            m.goal_rtn = SVector(1.0, 2.0, 3.0)
            SM.update_rpo_plan_buffer!(m.plan_buffer, SM.RPOPlan(valid=true, t_ref_s=[0.0]), 4.0)
            push!(m.replanning_events, (time_s=4.0, action=:brake))
            m.replan_count, m.safe_hold_count, m.replanning_phase, m.force_replan = 2, 1, :braking, true
        end)
        @test model.plan_buffer === buffer && !buffer.valid   # the controller's link, restored in place
    end
    @testset "The RPO MPC controller: held command, warm start, solver workspace, log" begin
        controller = only(rpo_configuration().control_model.control_effectors)
        controller.command_log = SM.RPOControlCommandLog()
        log = controller.command_log
        initialize!(controller)
        workspace = controller.controller.qp_model
        controller.held = SM.RPOHeldActuation(force_ii=SVector(1.0, 0.0, 0.0))
        controller.controller.U_prev .= 1.0
        push!(log.t_s, 1.0)
        initialize!(controller)
        @test controller.held.force_ii == Z3 && all(iszero, controller.controller.U_prev)
        @test controller.controller.qp_model !== workspace
        @test controller.command_log === log && isempty(log.t_s)
    end
    @testset "The robot-arm controller: held command from before the first run, zero warm start" begin
        arm = SM.default_robot_arm_model()
        configured = SM.RobotArmHeldActuation(joint_torque_nm=[0.5, -0.25, 0.0], base_force_ii=SVector(1.0, 2.0, 3.0))
        model = SM.RobotArmControlEffector(held=configured, controller=SM.init_robot_arm_joint_mpc(arm; dt_s=0.1, horizon=4))
        @test starts_equal(model, m -> (m.held = SM.RobotArmHeldActuation(joint_torque_nm=[9.0, 9.0, 9.0]);
            m.controller.U_prev .= 1.0))
        @test model.held.base_force_ii == SVector(1.0, 2.0, 3.0) && all(iszero, model.controller.U_prev)
    end
    @testset "The momentum manager: accumulator, held command and ticks" begin
        model = SM.ControlHooks.MagneticMomentumManagerModel(h_wheels_0=SVector(0.1, 0.0, 0.0))
        @test starts_equal(model, m -> (m.h_wheels = SVector(5.0, 5.0, 5.0); m.held_torque_body = SVector(1.0, 0.0, 0.0);
            m.initialized = true; m.ticks = 30; m.last_update_s = 29.0; push!(m.mounted_levels, 1.0)))
        @test !model.initialized && model.ticks == 0 && isempty(model.mounted_levels)
    end
    @testset "The burn scheduler: what it wrote returns, a change made since stays" begin
        model = SM.ThrusterModels.BaseThrusterModel(thrust=[50.0, 50.0], direction=[0.0, 0.0], Δv=[0.1, 0.2], start_burn_time=[5.0, 6.0],
            stop_burn_time=[20.0, 21.0], Isp=[300.0, 300.0])
        initialize!(model)
        # What a run leaves: both schedules cleared and replanned by the scheduler, a
        # burn on, a refusal recorded.
        for i in 1:2
            model.start_burn_time[i], model.stop_burn_time[i] = -1.0, -1.0
            SM.ControlHooks._note_schedule!(model, i)
        end
        model.burn_windows[1] = (5.0, 20.0, 50.0)
        model.burning[1] = true
        push!(model.plan_rejections[1], (time_s=1.0, rejection=:insufficient_propellant, delta_v_mps=0.1,
            direction_rad=0.0, source_orbit=-1))
        model.Δv[2] = 0.5   # the user's change between the runs
        initialize!(model)
        @test model.start_burn_time == [5.0, 6.0] && model.stop_burn_time == [20.0, 21.0]
        @test model.Δv == [0.1, 0.5]
        @test isequal(model.burn_windows, [(NaN, NaN, 0.0), (NaN, NaN, 0.0)]) && !any(model.burning)
        @test all(isempty, model.plan_rejections)
        # The changed value is the configuration the next run starts from.
        model.Δv[2] = 0.7
        SM.ControlHooks._note_schedule!(model, 2)
        initialize!(model)
        @test model.Δv == [0.1, 0.5]
    end
    @testset "The periapsis-raise guidance: its command states" begin
        model = SM.ApoapsisTargetPeriapsisRaiseGuidanceModel(target_apoapsis_radius_m=PLANET.Rp_e + 800e3)
        @test starts_equal(model, m -> append!(m.command_state, [1, 2]))
        @test isempty(model.command_state)
    end
    @testset "The planner guidance refuses reuse: every run gets a fresh copy" begin
        err = try run_simulation(rpo_configuration(); isolate_state=false) catch e e end
        @test err isa ArgumentError && occursin("isolate_state=true", sprint(showerror, err))
    end
end

# ---- Navigation models and custom sensor state --------------------------------

# A gyro whose bias walks at each sample: custom state its initialize_sensor! resets
# (or, for `Reset=false`, does not).
mutable struct WalkingGyro{Reset} <: SpaceAGORA.AbstractSensor
    period_s::Float64
    seed::Int
    rng::Any
    latest::Any
    latest_time_s::Float64
    bias::SVector{3, Float64}
    samples::Int
end
WalkingGyro{R}() where {R} = WalkingGyro{R}(1.0, 7, nothing, SVector(NaN, NaN, NaN), NaN, Z3, 0)
SpaceAGORA.sensor_dimension(::WalkingGyro) = 3
SpaceAGORA.initialize_sensor!(s::WalkingGyro{true}, truth, ctx, rng) = (s.bias = Z3; s.samples = 0; nothing)
function SpaceAGORA.sensor_measure(s::WalkingGyro, truth, ctx, rng)
    s.samples += 1
    reading = truth.ω_body + s.bias
    s.bias += 1.0e-5 * SVector(randn(rng), randn(rng), randn(rng))
    return reading
end

# An estimator that averages the gyro's readings: run state its
# initialize_navigation! resets (or, for `Reset=false`, does not).
mutable struct AveragingEstimator{Reset} <: SpaceAGORA.AbstractNavigationModel
    sum::SVector{3, Float64}
    count::Int
end
AveragingEstimator{R}() where {R} = AveragingEstimator{R}(Z3, 0)
function SpaceAGORA.update_navigation!(m::AveragingEstimator, u, p, t::Float64, i::Int)
    reading = SpaceAGORA.sensor_measurement(p.args.dynamics_model.spacecraft[i], 1)
    all(isfinite, reading) || return nothing
    m.sum += reading
    m.count += 1
    return nothing
end
SpaceAGORA.initialize_navigation!(m::AveragingEstimator{true}, u, p, t) = (m.sum = Z3; m.count = 0; nothing)

# Two runs with isolate_state=false on one configuration: the estimator and the
# sensor as each run ended them, and the runs' saved times and states.
function reused_runs(estimator, sensor)
    args = C.conformance_fixture(navigation=(estimator,), sensors=(sensor,), navigation_rate_s=2.0,
        mission_time_s=20.0).configuration
    ends = map(1:2) do _
        solution = run_simulation(args; isolate_state=false, return_solution=true)
        (estimator=deepcopy(only(args.navigation_model.navigation_effectors)),
            sensor=deepcopy(only(SpaceAGORA.mounted_sensors(only(args.dynamics_model.spacecraft))).sensor),
            t=solution.t, u=collect.(solution.u))
    end
    return ends
end

# An estimator that sums the finite readings of the first mounted sensor, of any
# dimension, and resets its sums when a run starts.
mutable struct ReadingSum <: SpaceAGORA.AbstractNavigationModel
    total::Float64
    count::Int
end
ReadingSum() = ReadingSum(0.0, 0)
function SpaceAGORA.update_navigation!(m::ReadingSum, u, p, t::Float64, i::Int)
    reading = SpaceAGORA.sensor_measurement(p.args.dynamics_model.spacecraft[i], 1)
    all(isfinite, reading) || return nothing
    m.total += sum(reading)
    m.count += 1
    return nothing
end
SpaceAGORA.initialize_navigation!(m::ReadingSum, u, p, t) = (m.total = 0.0; m.count = 0; nothing)

@testset "Navigation models and custom sensor state reset when a run starts (G-8)" begin
    first_end, second_end = reused_runs(AveragingEstimator{true}(), WalkingGyro{true}())
    @test first_end.estimator.count > 0 && first_end.sensor.samples > 0
    @test C._same_state(first_end.estimator, second_end.estimator)
    @test C._same_state(first_end.sensor, second_end.sensor)
    @test first_end.t == second_end.t && first_end.u == second_end.u
    # Negative controls: without the estimator's or the sensor's reset, the second run
    # starts from the first run's state.
    stale = reused_runs(AveragingEstimator{false}(), WalkingGyro{true}())
    @test stale[2].estimator.count == 2 * stale[1].estimator.count
    stale = reused_runs(AveragingEstimator{true}(), WalkingGyro{false}())
    @test stale[2].sensor.samples == 2 * stale[1].sensor.samples
    @test stale[2].estimator.sum != stale[1].estimator.sum
    # The kit's checks report the same: the estimator and the gyro pass, and each
    # model without its reset fails its reset check.
    @test C.passed(C.check_navigation_model(AveragingEstimator{true}(); sensors=(WalkingGyro{true}(),)))
    @test C.passed(C.check_sensor(WalkingGyro{true}()))
    status(report, id) = only(check.status for check in report.checks if check.id === id)
    @test status(C.check_navigation_model(AveragingEstimator{false}(); sensors=(WalkingGyro{true}(),)),
        :initialize_resets) === :fail
    @test status(C.check_sensor(WalkingGyro{false}()), :custom_state_reset) === :fail
end

@testset "Every built-in sensor kind starts each run afresh (G-8)" begin
    # The engine resets the fields it owns (stream, latest reading and its time)
    # when a run starts; the built-in kinds keep no other run state.
    for sensor in sensor_instances()
        first_end, second_end = reused_runs(ReadingSum(), sensor)
        @testset let kind = typeof(sensor)
            @test first_end.estimator.count > 0 && isfinite(first_end.sensor.latest_time_s)
            @test C._same_state(first_end.sensor, second_end.sensor)
            @test C._same_state(first_end.estimator, second_end.estimator)
            @test first_end.t == second_end.t && first_end.u == second_end.u
        end
    end
    @test Set(typeof.(sensor_instances())) ⊇ Set((SpaceAGORA.GNSSReceiver, SpaceAGORA.SunSensor, SpaceAGORA.Accelerometer))
end

# One log shared by a link and a model of each kind, which record the order the
# engine initializes them in.
mutable struct InitLog
    events::Vector{Symbol}
end
struct LoggedGuidance <: SpaceAGORA.AbstractGuidanceModel
    log::InitLog
end
struct LoggedNavigation <: SpaceAGORA.AbstractNavigationModel
    log::InitLog
end
struct LoggedControl <: SpaceAGORA.AbstractControlEffectorModel
    log::InitLog
end
SpaceAGORA.gnc_links(g::LoggedGuidance) = (Produces(g.log),)
SpaceAGORA.gnc_links(c::LoggedControl) = (Consumes(c.log),)
SpaceAGORA.link_initialize!(log::InitLog, producers, consumers, u, p, t) = (push!(log.events, :link); nothing)
SpaceAGORA.initialize_navigation!(m::LoggedNavigation, u, p, t) = (push!(m.log.events, :navigation); nothing)
SpaceAGORA.initialize_guidance!(m::LoggedGuidance, u, p, t) = (push!(m.log.events, :guidance); nothing)
SpaceAGORA.initialize_control!(m::LoggedControl, u, p, t) = (push!(m.log.events, :control); nothing)
SpaceAGORA.update_guidance!(::LoggedGuidance, u, p, t::Float64, i::Int) = nothing
SpaceAGORA.update_navigation!(::LoggedNavigation, u, p, t::Float64, i::Int) = nothing
SpaceAGORA.update_control!(::LoggedControl, u, p, t::Float64, i::Int) = nothing
SpaceAGORA.control_wrench(::LoggedControl, x::SpaceAGORA.StateSample, p, t::Float64) = (Z3, Z3)

@testset "The engine initializes links, then navigation, guidance and control" begin
    log = InitLog(Symbol[])
    args = C.conformance_fixture(guidance=(LoggedGuidance(log),), navigation=(LoggedNavigation(log),),
        control=(LoggedControl(log),), mission_time_s=5.0).configuration
    run_simulation(args)
    @test isempty(log.events)   # isolate_state=true: the run initialized its own copies
    run_simulation(args; isolate_state=false)
    run_simulation(args; isolate_state=false)
    @test log.events == repeat([:link, :navigation, :guidance, :control], 2)
end

end # module GNCRuntimeContractTests
