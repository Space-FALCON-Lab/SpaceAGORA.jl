module OperationalSegmentsTests

using Test
using SpaceAGORA
using StaticArrays
using JSON
using OrdinaryDiffEq: DiscreteCallback
using SciMLBase: terminate!

const SM = SpaceAGORA.SimulationModel
const MW = SpaceAGORA.MissionWorkflow
const EP = SpaceAGORA.ExternalPropagationServices
const ZERO3 = SVector(0.0, 0.0, 0.0)
const R0 = SVector(7.0e6, 0.0, 0.0)
const V0 = SVector(10.0, 7500.0, 0.0)
const INITIAL_MASS = 15.0
const THRUST = 0.8
const CONSUMPTION = 0.1

# A command held by an ordinary control effector. Segment actions change only
# this command, while the same state and effector remain in the same solve.
mutable struct CommandedThrust <: SM.AbstractControlEffectorModel
    throttle::Float64
    spacecraft_idx::Int
end
SM.calcControlEffect!(::CommandedThrust, u, p, t::Float64, i::Int) = nothing
function SM.calcControlForceTorque(model::CommandedThrust, u::AbstractVector,
        p::SM.ODEParams, i::Int, t::Float64)
    level = i == model.spacecraft_idx ? model.throttle : 0.0
    return SVector(THRUST * level, 0.0, 0.0), ZERO3
end
function SM.calcControlMassFlowRate(model::CommandedThrust, u::AbstractVector,
        p::SM.ODEParams, i::Int, t::Float64)
    return i == model.spacecraft_idx ? -CONSUMPTION * model.throttle : 0.0
end

function configuration(; duration_s=8.0, dtmax=2.0, n_spacecraft=1, solver=:tsit5)
    planet = make_no_gram_planet(:earth)
    spacecraft = map(1:n_spacecraft) do i
        body = Link(root=true, m=10.0, ref_area=1.0)
        SpacecraftModel(; joints=Joint[], links=[body], root=body,
            prop_mass=5.0, inertia_tensor=body.inertia, id=i,
            initial_condition=SM.CartesianInitialCondition(R0 + SVector(0.0, 0.0, 100.0 * (i - 1)), V0))
    end
    config = make_example_config(; planet, spacecraft=first(spacecraft),
        mission_time=Float64(duration_s), initial_time=InitialTime(year=2024, month=1, day=1),
        dynamic_effectors=(), density_model=NoAtmosphereModel(),
        ephemerides_model=SimpleEphemeridesModel(), orientation_sim=true,
        keplerian=true, verbose=false, results=false,
        solver_config=SolverConfig(solver_mode=solver, parallel=false))
    return SM.SimConfig._with_configuration(config;
        dynamics_model=DynamicsModel(spacecraft, ()),
        control_model=ControlModel(control_effectors=(CommandedThrust(0.0, 1),), control_rates=[1.0]),
        mission_configuration=MissionConfiguration(mission_type=SM.MissionTime,
            mission_time=duration_s, orientation_sim=true, data_rate=0.2),
        integration_tolerances=SM.IntegrationTolerances(
            reltol_orbit=1e-11, abstol_orbit=1e-11,
            reltol_atmosphere=1e-11, abstol_atmosphere=1e-11,
            reltol_mass=1e-12, abstol_mass=1e-12,
            dt_max_orbit=dtmax, dt_max_atmosphere=dtmax))
end

function execute(plan, root, name; kwargs...)
    return run_mission(plan; output_dir=joinpath(root, name), save_scene=false,
        save_payload=false, return_solution=true, kwargs...)
end

function command_actions(name, level, observations, configurations)
    return config -> begin
        push!(configurations, config)
        command = only(config.control_model.control_effectors)
        function record(ctx, kind)
            push!(observations, (; segment=name, kind, time_s=ctx.time_s,
                state=copy(ctx.state), spacecraft=ctx.spacecraft,
                configuration=ctx.configuration))
            # Context state is a snapshot, not a way to reset physical state.
            ctx.state .= -9.0e20
            nothing
        end
        return (on_enter=ctx -> begin
                    record(ctx, :enter)
                    command.throttle = level
                    nothing
                end,
                on_exit=ctx -> begin
                    record(ctx, :exit)
                    command.throttle = 0.0
                    nothing
                end)
    end
end

function thrust_plan(builder=configuration; observations=Any[], configurations=Any[])
    return MissionPlan("coast, maneuver, coast", builder; segments=(
        MissionSegment(:initial_coast; until=MissionDeadline(2.0),
            build_actions=command_actions(:initial_coast, 0.0, observations, configurations)),
        MissionSegment(:maneuver; until=MissionDuration(3.0),
            build_actions=command_actions(:maneuver, 1.0, observations, configurations)),
        MissionSegment(:final_coast;
            build_actions=command_actions(:final_coast, 0.0, observations, configurations))))
end

# Closed-form solution for constant force and constant propellant consumption.
# No gravity model, external package, or reference simulation defines this oracle.
function analytic_state(t)
    burn_time = clamp(t - 2.0, 0.0, 3.0)
    mass = INITIAL_MASS - CONSUMPTION * burn_time
    delta_v = THRUST / CONSUMPTION * log(INITIAL_MASS / mass)
    burn_displacement = THRUST / CONSUMPTION *
        (burn_time - mass / CONSUMPTION * log(INITIAL_MASS / mass))
    displacement = burn_displacement + max(t - 5.0, 0.0) * delta_v
    return (pos=R0 + V0 * t + SVector(displacement, 0.0, 0.0),
        vel=V0 + SVector(delta_v, 0.0, 0.0), mass=mass)
end

@testset "Operational segment definitions reject ambiguous boundaries" begin
    @test SM.MissionTime isa SM.MissionType # existing configuration enum is unchanged
    for value in (0.0, -1.0, Inf, NaN)
        @test_throws ArgumentError MissionDeadline(value)
        @test_throws ArgumentError MissionDuration(value)
    end
    @test_throws ArgumentError MissionCondition(identity; direction=:either)
    @test_throws ArgumentError MissionSegment(" ")
    @test_throws ArgumentError MissionSegment(:bad; until=2.0)
    @test_throws ArgumentError MissionPlan("bad", configuration; segments=(nothing,))
    @test_throws ArgumentError MissionPlan("bad", configuration;
        segments=(MissionSegment(:same; until=MissionDeadline(2.0)), MissionSegment(:same)))
    @test_throws ArgumentError MissionPlan("bad", configuration;
        segments=(MissionSegment(:unbounded), MissionSegment(:later)))
    segment = MissionSegment(:coast; until=MissionDuration(1.0))
    @test MissionPlan("vector input", configuration; segments=[segment]).segments === (segment,)
end

@testset "One continuous solve executes coast, thrust, coast" begin
    mktempdir() do root
        observations, captured_configs = Any[], Any[]
        builder_calls = Ref(0)
        original = configuration(; n_spacecraft=2)
        plan = thrust_plan(() -> (builder_calls[] += 1; original);
            observations, configurations=captured_configs)
        mission = execute(plan, root, "maneuver")
        solution = mission.results.solution
        @test builder_calls[] == 1
        @test solution.t[1] == 0.0
        @test solution.t[end] == 8.0
        @test 2.0 in solution.t
        @test 5.0 in solution.t
        @test all(isfinite, reduce(vcat, solution.u))
        @test length(captured_configs) == 3
        @test all(c -> c === mission.results.configuration, captured_configs)
        @test only(original.control_model.control_effectors).throttle == 0.0
        @test mission.results.configuration !== original
        @test getproperty.(observations, :segment) ==
            [:initial_coast, :initial_coast, :maneuver, :maneuver, :final_coast, :final_coast]
        @test getproperty.(observations, :kind) == [:enter, :exit, :enter, :exit, :enter, :exit]
        @test getproperty.(observations, :time_s) == [0.0, 2.0, 2.0, 5.0, 5.0, 8.0]
        @test observations[2].state == observations[3].state
        @test observations[4].state == observations[5].state
        @test observations[1].state == first(solution.u)
        @test observations[end].state == last(solution.u)
        for observation in observations
            expected = analytic_state(observation.time_s)
            state = observation.state.sc[1]
            @test state.pos ≈ expected.pos atol=2e-6 rtol=0
            @test state.vel ≈ expected.vel atol=2e-8 rtol=0
            @test state.mass ≈ expected.mass atol=2e-9 rtol=0
            @test state.q == [0.0, 0.0, 0.0, 1.0]
            @test state.ω == zeros(3)
            @test observation.spacecraft[1].pos_ii ≈ state.pos
            @test observation.spacecraft[1].vel_ii ≈ state.vel
            @test observation.spacecraft[1].mass_kg ≈ state.mass
            @test observation.configuration === mission.results.configuration
        end
        untouched = last(solution.u).sc[2]
        @test untouched.pos ≈ R0 + SVector(0.0, 0.0, 100.0) + 8.0 * V0 atol=2e-6 rtol=0
        @test untouched.vel == V0
        @test untouched.mass == INITIAL_MASS
        @test getproperty.(mission.history.segments, :name) == ["initial_coast", "maneuver", "final_coast"]
        @test getproperty.(mission.history.segments, :start_s) == [0.0, 2.0, 5.0]
        @test getproperty.(mission.history.segments, :end_s) == [2.0, 5.0, 8.0]
        @test getproperty.(mission.history.transitions, :time_s) == [2.0, 5.0, 8.0]
        @test getproperty.(mission.history.transitions, :from) == ["initial_coast", "maneuver", "final_coast"]
        @test getproperty.(mission.history.transitions, :to) == ["maneuver", "final_coast", nothing]
        @test mission.history.transitions[end].reason == "horizon"
        @test mission.history_path == mission.prefix * "_mission.json"
        @test mission.history_path in mission.files
        saved = JSON.parsefile(mission.history_path)
        @test saved["schema"] == 1
        @test saved["name"] == mission.name
        @test [s["name"] for s in saved["segments"]] == ["initial_coast", "maneuver", "final_coast"]
        @test [t["time_s"] for t in saved["transitions"]] == [2.0, 5.0, 8.0]
        @test saved["transitions"][end]["to"] === nothing
    end
end

@testset "Physical boundaries locate roots independently of maximum step" begin
    mktempdir() do root
        for solver in (:tsit5, :auto_stiff, :rodas5p, :dp8), direction in (:rising, :falling), dtmax in (0.13, 2.0)
            residual = direction == :rising ?
                (ctx -> ctx.spacecraft[1].pos_ii[1] - (R0[1] + 23.5)) :
                (ctx -> (R0[1] + 23.5) - ctx.spacecraft[1].pos_ii[1])
            plan = MissionPlan("position crossing", () -> configuration(; dtmax, solver); segments=(
                MissionSegment(:approach; until=MissionCondition(residual; direction)),
                MissionSegment(:after_crossing; until=MissionDuration(0.37))))
            run = execute(plan, root, "$(solver)-$(direction)-$(dtmax)")
            @test run.history.segments[1].end_s ≈ 2.35 atol=2e-7 rtol=0
            @test run.history.segments[2].start_s == run.history.segments[1].end_s
            @test run.history.segments[2].end_s ≈ 2.72 atol=2e-7 rtol=0
            @test run.results.solution.t[end] ≈ 2.72 atol=2e-7 rtol=0
            # Regular telemetry cadence is unchanged. The actual transition
            # remains in mission metadata beyond the final saved pose.
            @test run.results.table.time[end] < run.history.segments[end].end_s
            @test run.results.table.time[end] ≈ 2.6 atol=2e-7 rtol=0
            @test run.history.transitions[1].reason == "condition"
            @test run.history.transitions[end].reason == "duration"
            @test run.history.transitions[end].to === nothing
            @test last(run.results.solution.u).sc[1].mass == INITIAL_MASS
        end
        no_crossing = MissionPlan("no crossing", configuration; segments=(
            MissionSegment(:approach; until=MissionCondition(ctx -> ctx.spacecraft[1].pos_ii[1] - (R0[1] + 1e6))),
            MissionSegment(:unreached)))
        run = execute(no_crossing, root, "no-crossing")
        @test length(run.history.segments) == 1
        @test run.history.segments[1].name == "approach"
        @test run.history.segments[1].end_s == 8.0
        @test run.history.transitions[end].reason == "horizon"
    end
end

@testset "Durations meet the horizon despite floating-point addition" begin
    mktempdir() do root
        plan = MissionPlan("fractional duration", () -> configuration(; duration_s=0.3); segments=(
            MissionSegment(:first; until=MissionDuration(0.1)),
            MissionSegment(:last; until=MissionDuration(0.2))))
        run = execute(plan, root, "fractional")
        @test last(run.results.solution.t) == 0.3
        @test getproperty.(run.history.segments, :end_s) == [0.1, 0.3]
        @test run.history.transitions[end].reason == "duration"
        too_long = MissionPlan("past horizon", () -> configuration(; duration_s=0.3); segments=(
            MissionSegment(:first; until=MissionDuration(0.1)),
            MissionSegment(:last; until=MissionDuration(0.200001))))
        @test_throws ArgumentError execute(too_long, root, "past-horizon")
        explicit_late = MissionPlan("explicit deadline", () -> configuration(; duration_s=0.3); segments=(
            MissionSegment(:last; until=MissionDeadline(nextfloat(0.3))),))
        @test_throws ArgumentError execute(explicit_late, root, "explicit-late")
        @test !ispath(joinpath(root, "explicit-late"))
    end
end

@testset "Consecutive physical roots and deadline-to-condition transitions" begin
    mktempdir() do root
        for first_boundary in (:condition, :deadline)
            entry = first_boundary == :condition ?
                MissionCondition(ctx -> ctx.spacecraft[1].pos_ii[1] - (R0[1] + 23.5)) : MissionDeadline(2.35)
            function burn_actions(config)
                controller = only(config.control_model.control_effectors)
                return (on_enter=_ -> (controller.throttle = 1.0),
                    on_exit=_ -> (controller.throttle = 0.0))
            end
            plan = MissionPlan("consecutive roots", configuration; segments=(
                MissionSegment(:approach; until=entry),
                MissionSegment(:burn; until=MissionCondition(ctx -> ctx.spacecraft[1].mass_kg - 14.8;
                    direction=:falling), build_actions=burn_actions),
                MissionSegment(:coast; until=MissionDuration(0.37))))
            run = execute(plan, root, String(first_boundary))
            @test getproperty.(run.history.segments, :name) == ["approach", "burn", "coast"]
            @test getproperty.(run.history.transitions, :time_s) ≈ [2.35, 4.35, 4.72] atol=3e-7 rtol=0
            @test getproperty.(run.history.transitions, :reason) ==
                [first_boundary == :condition ? "condition" : "deadline", "condition", "duration"]
            final = last(run.results.solution.u).sc[1]
            @test final.mass ≈ 14.8 atol=2e-10 rtol=0
            @test final.vel ≈ V0 + SVector(THRUST / CONSUMPTION * log(15.0 / 14.8), 0.0, 0.0) atol=3e-8 rtol=0
            @test last(run.results.solution.t) ≈ 4.72 atol=3e-7 rtol=0
        end
    end
end

@testset "Supported solvers refresh the commanded dynamics" begin
    mktempdir() do root
        expected = analytic_state(8.0)
        for solver in (:auto_stiff, :rodas5p, :dp8)
            run = execute(thrust_plan(() -> configuration(; solver)), root, String(solver))
            final = last(run.results.solution.u).sc[1]
            @test final.mass ≈ expected.mass atol=2e-8 rtol=0
            @test final.vel ≈ expected.vel atol=3e-7 rtol=0
            @test final.pos ≈ expected.pos atol=2e-5 rtol=0
            @test getproperty.(run.history.segments, :end_s) == [2.0, 5.0, 8.0]
        end
    end
end

@testset "Saved operational history survives moving the viewer bundle" begin
    mktempdir() do root
        calls = Ref(0)
        plan = MissionPlan("portable segments", () -> (calls[] += 1; configuration()); segments=(
            MissionSegment(:first; until=MissionDeadline(2.35)),
            MissionSegment(:second; until=MissionDuration(0.37))))
        run = run_mission(plan; output_dir=joinpath(root, "original"))
        payload = JSON.parsefile(run.payload_path)
        @test payload["mission"] == JSON.parsefile(run.history_path)
        @test payload["mission"]["segments"][end]["end_s"] ≈ 2.72 atol=2e-7 rtol=0
        moved = joinpath(root, "moved")
        cp(run.directory, moved)
        rm(run.directory; recursive=true)
        prefix = joinpath(moved, basename(run.prefix))
        reconstructed = SM.SceneVisualization.export_viewer_payload(prefix;
            out=joinpath(moved, "rebuilt.json"), textures=false)
        @test JSON.parsefile(reconstructed) == payload
        @test calls[] == 1
    end
end

@testset "Repeated missions isolate commands and action factories" begin
    mktempdir() do root
        original = configuration()
        observations, captured = Any[], Any[]
        plan = thrust_plan(() -> original; observations, configurations=captured)
        first_run = execute(plan, root, "first")
        second_run = execute(plan, root, "second")
        @test length(captured) == 6
        @test captured[1] !== captured[4]
        @test only(captured[1].control_model.control_effectors) !== only(captured[4].control_model.control_effectors)
        @test only(original.control_model.control_effectors).throttle == 0.0
        @test first_run.results.solution.t == second_run.results.solution.t
        @test first_run.results.solution.u == second_run.results.solution.u
        @test first_run.history == second_run.history
        @test JSON.parsefile(first_run.history_path) == JSON.parsefile(second_run.history_path)
    end
end

@testset "Earlier engine termination closes the active segment once" begin
    mktempdir() do root
        exits = Float64[]
        stop = MissionEvent(:stop, _ -> DiscreteCallback((u, t, i) -> t >= 2.1,
            i -> terminate!(i); save_positions=(false, false)))
        segment = MissionSegment(:operation;
            build_actions=_ -> (on_exit=ctx -> push!(exits, ctx.time_s),))
        mission = execute(MissionPlan("early stop", configuration;
            segments=(segment,), events=(stop,)), root, "early-stop")
        final_time = last(mission.results.solution.t)
        @test 2.1 <= final_time < 8.0
        @test exits == [final_time]
        @test only(mission.history.segments).end_s == final_time
        @test only(mission.history.transitions).time_s == final_time
        @test only(mission.history.transitions).to === nothing
        @test only(mission.history.transitions).reason == "simulation_stop"
    end
end

struct UnsupportedOwner <: EP.AbstractExternalPropagator end

@testset "Invalid operations fail explicitly" begin
    mktempdir() do root
        cases = (
            ("beyond horizon", (MissionSegment(:one; until=MissionDeadline(9.0)),)),
            ("past deadline", (MissionSegment(:one; until=MissionDeadline(2.0)), MissionSegment(:two; until=MissionDeadline(1.0)))),
            ("long duration", (MissionSegment(:one; until=MissionDuration(9.0)),)),
            ("zero at entry", (MissionSegment(:one; until=MissionCondition(ctx -> ctx.spacecraft[1].pos_ii[1] - R0[1])),)),
            ("already rising", (MissionSegment(:one; until=MissionCondition(_ -> 1.0; direction=:rising)),)),
            ("already falling", (MissionSegment(:one; until=MissionCondition(_ -> -1.0; direction=:falling)),)),
            ("nonfinite", (MissionSegment(:one; until=MissionCondition(_ -> NaN)),)),
            ("infinite", (MissionSegment(:one; until=MissionCondition(_ -> Inf)),)),
            ("nonreal", (MissionSegment(:one; until=MissionCondition(_ -> 1 + im)),)),
            ("bad actions", (MissionSegment(:one; build_actions=_ -> 42),)),
            ("bad handler", (MissionSegment(:one; build_actions=_ -> (on_enter=42,)),)),
        )
        for (name, segments) in cases
            @test_throws ArgumentError execute(MissionPlan(name, configuration; segments), root, name)
        end
        segments = (MissionSegment(:operation),)
        external = () -> SM.SimConfig._with_configuration(configuration(); external_propagators=(UnsupportedOwner(),))
        @test_throws ArgumentError execute(MissionPlan("external owner", external; segments), root, "external")
        for mode in (:gravity_backbone_split, :multirate, :split_imex, :symplectic)
            builder = () -> configuration(; solver=mode)
            destination = joinpath(root, String(mode))
            @test_throws ArgumentError execute(MissionPlan("unsupported solve", builder; segments), root, String(mode))
            @test !ispath(destination)
        end
    end
end

end # module OperationalSegmentsTests
