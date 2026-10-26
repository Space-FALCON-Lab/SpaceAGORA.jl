# The run-time declarations of control models (GNCRuntime), the run-start
# initialization of navigation effectors and of a sensor's own state, and the
# models' own run state in checkpoint files.
module GNCRuntimeDeclarationsTests
using Test, StaticArrays, LinearAlgebra, Logging, SpaceAGORA
using DiffEqCallbacks: PresetTimeCallback
using OrdinaryDiffEq: DiscreteCallback
import SpaceAGORA.TelemetryVerification: make_example_config
const SM = SpaceAGORA.SimulationModel
const L = SpaceAGORA.SimulationLifecycle
const ZERO3 = SVector(0.0, 0.0, 0.0)

# The order in which the run-start initializers ran.
const LOG = Symbol[]

struct LoggedNavigation <: AbstractNavigationModel end
SpaceAGORA.update_navigation!(::LoggedNavigation, u, p, t::Float64, i::Int) = nothing
L.initialize_navigation!(::LoggedNavigation, u, p, t) = (push!(LOG, :navigation); nothing)

struct LoggedGuidance <: AbstractGuidanceModel end
SpaceAGORA.update_guidance!(::LoggedGuidance, u, p, t::Float64, i::Int) = nothing
L.initialize_guidance!(::LoggedGuidance, u, p, t) = (push!(LOG, :guidance); nothing)

# A periodic controller that records its clock starts, ticks and stop requests.
mutable struct ClockedControl <: SM.AbstractControlEffectorModel
    clock::Vector{Tuple{Float64, Float64}}
    ticks::Int
    stops::Int
end
ClockedControl() = ClockedControl(Tuple{Float64, Float64}[], 0, 0)
SpaceAGORA.update_control!(m::ClockedControl, u, p, t, i) = (m.ticks += 1; nothing)
SpaceAGORA.control_wrench(::ClockedControl, x, p, t) = (ZERO3, ZERO3)
L.initialize_control!(::ClockedControl, u, p, t) = (push!(LOG, :control); nothing)
SpaceAGORA.gnc_clock_started!(m::ClockedControl, period_s, origin_s) = (push!(m.clock, (period_s, origin_s)); nothing)
SpaceAGORA.register_tstops!(integrator, m::ClockedControl, spacecraft_index) = (m.stops += 1; nothing)

# An event-driven controller: no periodic tick, one callback of its own.
mutable struct ScheduledControl <: SM.AbstractControlEffectorModel
    fired::Vector{Float64}
    ticks::Int
    period_s::Float64
end
ScheduledControl() = ScheduledControl(Float64[], 0, NaN)
SpaceAGORA.update_control!(m::ScheduledControl, u, p, t, i) = (m.ticks += 1; nothing)
SpaceAGORA.control_wrench(::ScheduledControl, x, p, t) = (ZERO3, ZERO3)
SpaceAGORA.tick_schedule(::ScheduledControl) = :event_driven
function SpaceAGORA.gnc_callbacks(m::ScheduledControl, num_sats, args, period_s)
    m.period_s = period_s
    return (PresetTimeCallback([1.5], integrator -> push!(m.fired, integrator.t); save_positions=(false, false)),)
end

# A schedule the engine does not know.
struct OddControl <: SM.AbstractControlEffectorModel end
SpaceAGORA.update_control!(::OddControl, u, p, t, i) = nothing
SpaceAGORA.control_wrench(::OddControl, x, p, t) = (ZERO3, ZERO3)
SpaceAGORA.tick_schedule(::OddControl) = :sometimes

# A thread-safe periodic law that records, per spacecraft, the thread of its
# last tick and of its last stop registration.
struct ThreadedLaw <: SM.AbstractControlEffectorModel
    tick_thread::Vector{Int}
    stop_thread::Vector{Int}
end
ThreadedLaw(n) = ThreadedLaw(zeros(Int, n), zeros(Int, n))
SpaceAGORA.update_control!(m::ThreadedLaw, u, p, t, i) = (m.tick_thread[i] = Threads.threadid(); nothing)
SpaceAGORA.control_wrench(::ThreadedLaw, x, p, t) = (ZERO3, ZERO3)
SM.SimulationCallbacks.control_model_threadsafe(::ThreadedLaw) = true
SpaceAGORA.register_tstops!(integrator, m::ThreadedLaw, i) = (m.stop_thread[i] = Threads.threadid(); nothing)

function constellation_args(control; n=4, mission_time=5.0)
    planet = SM.make_no_gram_planet(:earth)
    sats = map(1:n) do id
        root = SM.Link(root=true, m=100.0, ref_area=2.0)
        ic = SM.InitialCondition(ra=planet.Rp_e + 700e3, rp=planet.Rp_e + 650e3, i=45.0, ω=0.0, Ω=30.0id, ν=10.0)
        SM.SpacecraftModel(joints=SM.Joint[], links=[root], root=root, prop_mass=0.0, inertia_tensor=root.inertia,
            initial_condition=ic, id=id)
    end
    env = SM.EnvironmentModel(planet=planet, EI=120.0, density_model=SM.NoAtmosphereModel(),
        ephemerides_model=SM.SimpleEphemeridesModel(),
        thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet), topography=false, wind=false)
    return SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=mission_time, orientation_sim=false, num_steps_to_save=100),
        environment_model=env,
        dynamics_model=SM.DynamicsModel(SM.SpacecraftModel[sats...], (SM.InverseSquaredJ2GravityModel(),)),
        guidance_model=SM.GuidanceModel((), Float64[]), navigation_model=SM.NavigationModel((), Float64[]),
        control_model=SM.ControlModel(control_effectors=(control,), control_rates=[1.0]),
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
end

# A periodic law that declares the period it assumes.
struct PeriodicLaw <: SM.AbstractControlEffectorModel
    period_s::Float64
end
SpaceAGORA.update_control!(::PeriodicLaw, u, p, t, i) = nothing
SpaceAGORA.control_wrench(::PeriodicLaw, x, p, t) = (ZERO3, ZERO3)
SpaceAGORA.gnc_period_s(m::PeriodicLaw) = m.period_s

# The message of the `ArgumentError` that `f()` refuses with, or `nothing` when it
# returns.
function refusal(f)
    try
        f()
        return nothing
    catch e
        e isa ArgumentError || rethrow()
        return sprint(showerror, e)
    end
end
period_refusal(args) = refusal(() -> SM.SimulationCallbacks.check_control_periods(args))

# A sensor with state of its own: the number of samples since the run started.
mutable struct CountingSensor <: AbstractSensor
    period_s::Float64
    seed::Int
    rng::Any
    latest::SVector{1, Float64}
    latest_time_s::Float64
    count::Int
    initialized_at::Vector{Float64}
end
CountingSensor() = CountingSensor(1.0, 1, nothing, SVector(NaN), NaN, 0, Float64[])
SpaceAGORA.sensor_dimension(::CountingSensor) = 1
SpaceAGORA.sensor_measure(s::CountingSensor, truth::SensorTruth, ctx, rng) = SVector(Float64(s.count += 1))
function SpaceAGORA.initialize_sensor!(s::CountingSensor, truth::SensorTruth, ctx, rng)
    s.count = 0
    push!(s.initialized_at, truth.t_s)
    return nothing
end

function runtime_args(control; sensor=CountingSensor(), mission_time=4.0, guidance=(LoggedGuidance(),),
                      navigation=(LoggedNavigation(),), navigation_period=1.0, settings=(;))
    root = SM.Link(root=true, m=50.0, dims=MVector(1.0, 1.0, 1.0))
    r = SM.Earth().Rp_e + 500e3
    ic = SM.CartesianInitialCondition(SVector(r, 0.0, 0.0), SVector(0.0, sqrt(SM.Earth().μ / r), 0.0);
        q=SVector(0.0, 0.0, 0.0, 1.0), ang_vel=ZERO3)
    sc = SM.SpacecraftModel(links=[root], root=root, inertia_tensor=root.inertia, initial_condition=ic, id=3)
    add_sensor!(sc, root, sensor)
    base = make_example_config(planet=SM.Earth(), spacecraft=sc, mission_time=mission_time,
        initial_time=SM.InitialTime(year=2024, month=3, day=20, hour=0, minute=0, second=0.0),
        dynamic_effectors=(), ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=true,
        keplerian=false, EI_km=100.0, verbose=false, results=false)
    extra = (;)
    if !isempty(settings)
        st = base.simulation_settings
        fields = NamedTuple{fieldnames(typeof(st))}(map(n -> getfield(st, n), fieldnames(typeof(st))))
        extra = (; simulation_settings=SM.SimulationSettings(; merge(fields, settings)...))
    end
    return SM.SimConfig._with_configuration(base; extra...,
        navigation_model=SM.NavigationModel(navigation_effectors=navigation, navigation_rates=fill(navigation_period, length(navigation))),
        guidance_model=SM.GuidanceModel(guidance_effectors=guidance, guidance_rates=fill(1.0, length(guidance))),
        control_model=SM.ControlModel(control_effectors=control, control_rates=fill(0.5, length(control))))
end

# Run state kept through checkpoint files: a navigation model counting its ticks
# (and, apart, this run's, which the file does not hold), a sensor counting its
# samples, and a model that saves state it cannot restore.
mutable struct TickCounter <: AbstractNavigationModel
    ticks::Int
    this_run::Int
    calls::Vector{Symbol}
end
TickCounter() = TickCounter(0, 0, Symbol[])
SpaceAGORA.update_navigation!(m::TickCounter, u, p, t::Float64, i::Int) = (m.ticks += 1; m.this_run += 1; nothing)
L.initialize_navigation!(m::TickCounter, u, p, t) = (m.ticks = 0; m.this_run = 0; push!(m.calls, :initialize); nothing)
L.checkpoint_state(m::TickCounter) = (ticks=m.ticks,)
L.restore_checkpoint_state!(m::TickCounter, data) = (m.ticks = data.ticks; push!(m.calls, :restore); nothing)
mutable struct KeptCountSensor <: AbstractSensor
    period_s::Float64
    seed::Int
    rng::Any
    latest::SVector{1, Float64}
    latest_time_s::Float64
    count::Int
    initialized_at::Vector{Float64}
end
KeptCountSensor() = KeptCountSensor(1.0, 1, nothing, SVector(NaN), NaN, 0, Float64[])
SpaceAGORA.sensor_dimension(::KeptCountSensor) = 1
SpaceAGORA.sensor_measure(s::KeptCountSensor, truth::SensorTruth, ctx, rng) = SVector(Float64(s.count += 1))
SpaceAGORA.initialize_sensor!(s::KeptCountSensor, truth::SensorTruth, ctx, rng) = (s.count = 0; push!(s.initialized_at, truth.t_s); nothing)
L.checkpoint_state(s::KeptCountSensor) = (count=s.count,)
L.restore_checkpoint_state!(s::KeptCountSensor, data) = (s.count = data.count; nothing)
# A kind that keeps state only when its parameter says so.
mutable struct MaybeKept <: AbstractNavigationModel
    keeps::Bool
    ticks::Int
end
SpaceAGORA.update_navigation!(m::MaybeKept, u, p, t::Float64, i::Int) = (m.ticks += 1; nothing)
L.checkpoint_state(m::MaybeKept) = m.keeps ? (ticks=m.ticks,) : nothing
L.restore_checkpoint_state!(m::MaybeKept, data) = (m.ticks = data.ticks; nothing)
struct HalfKept <: AbstractNavigationModel end
SpaceAGORA.update_navigation!(::HalfKept, u, p, t::Float64, i::Int) = nothing
L.checkpoint_state(::HalfKept) = (ticks=0,)

# What a resumed run's start shows, in order: a navigation model records what its
# initializer sees (the callbacks' and the models' resume flags) and when its own
# run state comes back, and a callback records what the callbacks see when they
# initialize (the model's ticks and the sensor's latest sample time).
const RESUME_LOG = Any[]
mutable struct ResumeProbe <: AbstractNavigationModel
    ticks::Int
end
SpaceAGORA.update_navigation!(m::ResumeProbe, u, p, t::Float64, i::Int) = (m.ticks += 1; nothing)
function L.initialize_navigation!(m::ResumeProbe, u, p, t)
    m.ticks = 0
    push!(RESUME_LOG, (:initialize, p.shared_buffers.callback_resume[], p.shared_buffers.model_resume[]))
    return nothing
end
L.checkpoint_state(m::ResumeProbe) = (ticks=m.ticks,)
L.restore_checkpoint_state!(m::ResumeProbe, data) = (m.ticks = data.ticks; push!(RESUME_LOG, (:restore, data.ticks)); nothing)
resume_probe_callback(args) = DiscreteCallback((u, t, integrator) -> false, integrator -> nothing;
    initialize=(cb, u, t, integrator) -> begin
        a = integrator.p.args
        sensor = only(SM.mounted_sensors(only(a.dynamics_model.spacecraft))).sensor
        push!(RESUME_LOG, (:callbacks, only(a.navigation_model.navigation_effectors).ticks, sensor.latest_time_s))
    end, save_positions=(false, false))

refusal(f) = try
    f()
    nothing
catch e
    e isa ArgumentError || rethrow()
    sprint(showerror, e)
end

# A contributor's link: a slot one guidance model writes and one control model reads.
mutable struct SharedSlot
    initialized::Int
    checked::Vector{Tuple{Vector{Any}, Vector{Any}}}
end
SharedSlot() = SharedSlot(0, Tuple{Vector{Any}, Vector{Any}}[])
SpaceAGORA.link_cardinality(::SharedSlot) = (producers = 1:1, consumers = 0:1)
SpaceAGORA.check_link(s::SharedSlot, producers, consumers, args) = (push!(s.checked, (producers, consumers)); nothing)
SpaceAGORA.link_initialize!(s::SharedSlot, producers, consumers, u, p, t) = (s.initialized += 1; push!(LOG, :link); nothing)
mutable struct SlotWriter <: AbstractGuidanceModel
    slot::SharedSlot
    prepared::Int
end
SpaceAGORA.update_guidance!(::SlotWriter, u, p, t::Float64, i::Int) = nothing
SpaceAGORA.gnc_links(g::SlotWriter) = (Produces(g.slot),)
L.initialize_guidance!(::SlotWriter, u, p, t) = (push!(LOG, :guidance); nothing)
mutable struct Bystander <: AbstractGuidanceModel
    prepared::Int
end
SpaceAGORA.update_guidance!(::Bystander, u, p, t::Float64, i::Int) = nothing
mutable struct SlotReader <: SM.AbstractControlEffectorModel
    slot::SharedSlot
    ticks::Int
end
SpaceAGORA.update_control!(c::SlotReader, u, p, t, i) = (c.ticks += 1; nothing)
SpaceAGORA.control_wrench(::SlotReader, x, p, t) = (ZERO3, ZERO3)
SpaceAGORA.gnc_links(c::SlotReader) = (Consumes(c.slot),)
L.before_reference_control!(g::SlotWriter, c::SlotReader, u, p, t) = (g.prepared += 1; nothing)
L.before_reference_control!(g::Bystander, c::SlotReader, u, p, t) = (g.prepared += 1; nothing)
struct BadLinks <: AbstractGuidanceModel end
SpaceAGORA.update_guidance!(::BadLinks, u, p, t::Float64, i::Int) = nothing
SpaceAGORA.gnc_links(::BadLinks) = (:not_a_link,)

# What ran, and when: a sensor, navigation, a guidance model that produces the
# link a controller consumes, plain guidance, periodic controllers, an
# event-driven controller's own callback, a saved field, and priority and extra
# callbacks each record themselves, to check the order of callbacks that apply at
# one time ("Callback order" in docs/contributing/architecture.md).
const ORDER = Tuple{Symbol, Float64}[]
mutable struct OrderSensor <: AbstractSensor
    period_s::Float64
    seed::Int
    rng::Any
    latest::SVector{1, Float64}
    latest_time_s::Float64
end
OrderSensor() = OrderSensor(1.0, 1, nothing, SVector(NaN), NaN)
SpaceAGORA.sensor_dimension(::OrderSensor) = 1
SpaceAGORA.sensor_measure(::OrderSensor, truth::SensorTruth, ctx, rng) = (push!(ORDER, (:sensor, truth.t_s)); SVector(0.0))
SpaceAGORA.initialize_sensor!(::OrderSensor, truth::SensorTruth, ctx, rng) = nothing
struct OrderNavigation <: AbstractNavigationModel end
SpaceAGORA.update_navigation!(::OrderNavigation, u, p, t::Float64, i::Int) = (push!(ORDER, (:navigation, t)); nothing)
mutable struct OrderSlot end
struct OrderWriter <: AbstractGuidanceModel
    slot::OrderSlot
end
SpaceAGORA.gnc_links(g::OrderWriter) = (Produces(g.slot),)
SpaceAGORA.update_guidance!(::OrderWriter, u, p, t::Float64, i::Int) = (push!(ORDER, (:writer, t)); nothing)
struct OrderGuidance <: AbstractGuidanceModel end
SpaceAGORA.update_guidance!(::OrderGuidance, u, p, t::Float64, i::Int) = (push!(ORDER, (:guidance, t)); nothing)
struct OrderControl <: SM.AbstractControlEffectorModel
    name::Symbol
    slot::Union{Nothing, OrderSlot}
end
SpaceAGORA.update_control!(c::OrderControl, u, p, t, i) = (push!(ORDER, (c.name, t)); nothing)
SpaceAGORA.control_wrench(::OrderControl, x, p, t) = (ZERO3, ZERO3)
SpaceAGORA.gnc_links(c::OrderControl) = c.slot === nothing ? () : (Consumes(c.slot),)
L.before_reference_control!(::OrderWriter, c::OrderControl, u, p, t) = (push!(ORDER, (:prepare, t)); nothing)
# A callback that records `name` at whole seconds, where the 0.5 s and 1 s ticks meet.
order_callback(name) = DiscreteCallback((u, t, integrator) -> t > 0 && isinteger(t),
    integrator -> (push!(ORDER, (name, integrator.t)); nothing); save_positions=(false, false))
struct OrderEvents <: SM.AbstractControlEffectorModel end
SpaceAGORA.update_control!(::OrderEvents, u, p, t, i) = nothing
SpaceAGORA.control_wrench(::OrderEvents, x, p, t) = (ZERO3, ZERO3)
SpaceAGORA.tick_schedule(::OrderEvents) = :event_driven
SpaceAGORA.gnc_callbacks(::OrderEvents, num_sats, args, period_s) = (order_callback(:event),)

# Where the hardware's run state starts a run: a controller holds a hinged panel
# at 0.5 rad, and a priority callback's initialization and a sensor's samples
# record the hinge's angle when they read it.
const SEEN = Tuple{Symbol, Float64, Float64}[]
struct HoldHinge <: SM.AbstractControlEffectorModel end
SpaceAGORA.actuation_mode(::HoldHinge) = :mounted
SpaceAGORA.actuator_levels(::HoldHinge, i::Int) = [0.5]
SpaceAGORA.update_control!(::HoldHinge, u, p, t, i) = nothing
mutable struct HingeSensor <: AbstractSensor
    period_s::Float64
    seed::Int
    rng::Any
    latest::SVector{1, Float64}
    latest_time_s::Float64
    hinge::PanelHinge
end
HingeSensor(hinge) = HingeSensor(1.0, 1, nothing, SVector(NaN), NaN, hinge)
SpaceAGORA.sensor_dimension(::HingeSensor) = 1
SpaceAGORA.sensor_measure(s::HingeSensor, truth::SensorTruth, ctx, rng) =
    (push!(SEEN, (:sensor, truth.t_s, s.hinge.angle_rad)); SVector(0.0))
SpaceAGORA.initialize_sensor!(::HingeSensor, truth::SensorTruth, ctx, rng) = nothing
hinge_seer(hinge) = DiscreteCallback((u, t, integrator) -> false, integrator -> nothing;
    initialize=(cb, u, t, integrator) -> (push!(SEEN, (:priority, t, hinge.angle_rad)); nothing), save_positions=(false, false))

@testset "GNC runtime declarations" begin
    @testset "defaults keep the periodic tick" begin
        @test tick_schedule(:any_model) === :periodic
        @test gnc_callbacks(:any_model, 2, nothing, 0.1) === ()
        @test register_tstops!(nothing, :any_model, 1) === nothing
        @test gnc_clock_started!(:any_model, 0.1, 0.0) === nothing
        @test L.initialize_navigation!(:any_model, nothing, nothing, 0.0) === nothing
        @test checkpoint_state(:any_model) === nothing && checkpoint_state === L.checkpoint_state
        @test !applicable(restore_checkpoint_state!, :any_model, (ticks=0,))   # no default: a model that saves state restores it
        @test tick_schedule === SpaceAGORA.GNCRuntime.tick_schedule
        @test initialize_sensor! === SM.Sensors.initialize_sensor!
        # The sensor fallback types nothing, so a sensor's method that types only
        # some arguments is not ambiguous with it.
        s = CountingSensor(); s.count = 3
        truth = SensorTruth(ZERO3, ZERO3, SVector(0.0, 0.0, 0.0, 1.0), ZERO3, 1.0, 0.0, ZERO3)
        @test initialize_sensor!(s, truth, actuator_context(SM.Link(root=true, m=1.0)), nothing) === nothing
        @test s.count == 0 && s.initialized_at == [0.0]
    end

    @testset "built-in models declare what the engine used to test by type" begin
        thruster = SM.BaseThrusterModel(thrust=[0.0], direction=[0.0], Δv=[0.0], start_burn_time=[0.0],
            stop_burn_time=[0.0], Isp=[300.0])
        @test tick_schedule(thruster) === :event_driven
        @test tick_schedule(ClockedControl()) === :periodic
        @test !isdefined(SM.SimulationCallbacks, :control_requires_periodic_callback)   # removed in 0.5.0
        # Initialization only; planning and burn edges run in the late phase after guidance.
        @test length(gnc_callbacks(thruster, 1, runtime_args((thruster,)), 0.5)) == 1
        own(f, T) = which(f, T).module !== SpaceAGORA.GNCRuntime
        @test own(register_tstops!, Tuple{Any, SM.BaseThrusterModel, Int})
        @test own(gnc_clock_started!, Tuple{SM.AerobrakingEnergyDepletionControlModel, Float64, Float64})
    end

    @testset "a run follows the declarations and initializes in order" begin
        empty!(LOG)
        sensor, clocked, scheduled = CountingSensor(), ClockedControl(), ScheduledControl()
        args = runtime_args((clocked, scheduled); sensor=sensor)
        ticks_per_run = Ref(0)
        for run in 1:2
            run_simulation(args; isolate_state=false)
            run == 1 && (ticks_per_run[] = clocked.ticks)
            @test LOG == repeat([:navigation, :guidance, :control], run)
            # Each run reports the clock when it starts, before the initializers, and
            # again when its solve starts.
            @test clocked.clock == fill((0.5, 0.0), 2run)
            @test ticks_per_run[] > 0 && clocked.ticks == run * ticks_per_run[]
            @test clocked.stops == clocked.ticks                   # one stop request per tick
            @test scheduled.ticks == 0 && scheduled.fired == fill(1.5, run) && scheduled.period_s == 0.5
            @test sensor.initialized_at == fill(0.0, run)         # once per run, before the first sample
            @test sensor.count == 5                               # this run's samples only: 0, 1, 2, 3, 4 s
        end
    end

    @testset "stop registration runs on one thread when the ticks are threaded" begin
        law = ThreadedLaw(4)
        withenv("SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "on", "SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD" => "2") do
            run_simulation(constellation_args(law); isolate_state=false)
        end
        @test all(>(0), law.tick_thread) && all(>(0), law.stop_thread)
        @test length(unique(law.stop_thread)) == 1                 # one thread registers every stop
        Threads.nthreads() > 1 && @test length(unique(law.tick_thread)) > 1   # while the ticks spread
    end

    @testset "declared periods are checked before the run" begin
        @test gnc_period_s(:any_model) === nothing
        # The MPC controllers keep the step they were discretized with.
        ctrl = SM.init_rpo_lqmpc(0.001, 0.2, Matrix(1.0I, 6, 6), Matrix(1.0I, 3, 3), Matrix(1.0I, 6, 6), 4)
        @test ctrl.dt_s == 0.2 && deepcopy(ctrl).dt_s == 0.2
        # A contributor's declared period against its callback's rate (0.5 s).
        @test period_refusal(runtime_args((PeriodicLaw(0.5),))) === nothing
        refused = period_refusal(runtime_args((PeriodicLaw(0.25),)))
        @test occursin("PeriodicLaw declares a 0.25 s period", refused) && occursin("every 0.5 s (control_rates[1])", refused)
        # A refusal, not a warning, whatever --depwarn says: nothing is logged.
        logger = Test.TestLogger(min_level=Logging.Warn)
        refused = with_logger(() -> period_refusal(runtime_args((PeriodicLaw(0.25),))), logger)
        @test refused isa String && isempty(logger.logs)
        # The RPO pilot builds its controller with the model's step; a different step is refused.
        rpo = make_rpo_configuration(planner=DirectRPOPlanner(), mission_time_s=0.2)
        mpc = only(rpo.control_model.control_effectors)
        @test gnc_period_s(mpc) == mpc.control_dt_s == mpc.controller.dt_s
        @test period_refusal(rpo) === nothing
        mpc.controller.dt_s = 2mpc.control_dt_s
        refused = period_refusal(rpo)
        @test occursin("RPOMPCControlModel ticks every 0.1 s but its controller was discretized with a 0.2 s step", refused)
        # The engine runs the check before every run.
        refused = refusal(() -> run_simulation(runtime_args((PeriodicLaw(0.25),))))
        @test occursin("0.25 s period", refused)
    end

    @testset "links pair models, and the engine checks and initializes each once" begin
        @test gnc_links(:any_model) === ()
        @test link_cardinality(:any_link) == (producers = 0:typemax(Int), consumers = 0:typemax(Int))
        @test check_link(:any_link, [], [], nothing) === nothing
        @test link_initialize!(:any_link, [], [], nothing, nothing, 0.0) === nothing
        empty!(LOG)
        slot = SharedSlot()
        writer, bystander, reader = SlotWriter(slot, 0), Bystander(0), SlotReader(slot, 0)
        args = runtime_args((reader,); guidance=(writer, bystander))
        groups = SM.SimulationCallbacks.gnc_link_groups(args)
        @test length(groups) == 1 && only(groups).link === slot
        @test only(groups).producers == [writer] && only(groups).consumers == [reader]
        run_simulation(args; isolate_state=false)
        @test slot.initialized == 1 && length(slot.checked) == 1
        @test slot.checked[1][1] == [writer] && slot.checked[1][2] == [reader]
        @test LOG[1:3] == [:link, :navigation, :guidance]               # the link before the models
        # The linked producer prepares the reader once per tick; the bystander never.
        @test reader.ticks > 0 && writer.prepared == reader.ticks && bystander.prepared == 0
        # A second writer breaks the link's cardinality; a malformed declaration is refused.
        err = try
            run_simulation(runtime_args((SlotReader(slot, 0),); guidance=(SlotWriter(slot, 0), SlotWriter(slot, 0))))
            nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin("2 models produce one SharedSlot", sprint(showerror, err))
        @test_throws ArgumentError SM.SimulationCallbacks.gnc_link_groups(runtime_args((ClockedControl(),); guidance=(BadLinks(),)))
        # SolverConfig(parallel=true) makes run_simulation re-enter itself; each link
        # is still checked once. The flag's persisted state goes to a scratch folder.
        mktempdir() do dir
            constants = joinpath(dir, "cost_constants.toml")
            SM.ParallelCost.ensure_machine_constants!(path=constants,
                calibrate=() -> SM.ParallelCost.calibrate_machine(k=1))
            withenv("SPACEAGORA_COST_CONSTANTS_PATH" => constants,
                    "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => joinpath(dir, "policy_state.toml"),
                    "SPACEAGORA_RHS_CALIBRATION_PATH" => joinpath(dir, "rhs_calibration.toml"),
                    "SPACEAGORA_OUTER_ROUTE_STATE_PATH" => joinpath(dir, "outer_route_state.toml"),
                    "SPACEAGORA_CAMPAIGN_CORRECTIONS_PATH" => joinpath(dir, "campaign_corrections.toml"),
                    "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing) do
                parallel_slot = SharedSlot()
                parallel_args = SM.SimConfig._with_configuration(
                    runtime_args((SlotReader(parallel_slot, 0),); guidance=(SlotWriter(parallel_slot, 0),));
                    solver_config=SM.SolverConfig(parallel=true))
                run_simulation(parallel_args; isolate_state=false)
                @test length(parallel_slot.checked) == 1 && parallel_slot.initialized == 1
            end
        end
    end

    @testset "callbacks that apply at one time run in the documented order" begin
        empty!(ORDER)
        slot = OrderSlot()
        args = runtime_args((OrderControl(:reader, slot), OrderEvents(), OrderControl(:control, nothing));
            sensor=OrderSensor(), guidance=(OrderWriter(slot), OrderGuidance()), navigation=(OrderNavigation(),))
        m = args.mission_configuration                       # saved every second
        args = SM.SimConfig._with_configuration(args; mission_configuration=SM.MissionConfiguration(
            m.mission_type, m.keplerian, m.number_of_orbits, m.mission_time, m.orientation_sim, m.num_steps_to_save, 1.0))
        saved = vcat(SM.default_save_fields(args), [SM.SaveField(:order, (u, t, integrator) -> (push!(ORDER, (:save, t)); 0.0))])
        run_simulation(args; return_results=true, save_fields=saved,
            priority_callback_factories=(a -> order_callback(:priority),),
            extra_callbacks=(order_callback(:extra),), extra_callback_factories=(a -> order_callback(:extra_factory),))
        # At 1, 2 and 3 s (before the run's 4 s end) every recorder applies. The
        # periodic control models tick in configuration order around the
        # event-driven model's callback; the reader's linked producer prepares it
        # first. Guidance that does not publish before control comes after it.
        expected = [:priority, :sensor, :navigation, :prepare, :reader, :event, :control, :writer, :guidance,
                    :save, :extra, :extra_factory]
        for t in (1.0, 2.0, 3.0)
            @test [name for (name, at) in ORDER if at == t] == expected
        end
        # Between them, at 0.5 s steps, only the 0.5 s control ticks apply.
        @test [name for (name, at) in ORDER if at == 1.5] == [:prepare, :reader, :control]
    end

    @testset "the hardware's run state starts before the callbacks that read it at the run's start" begin
        # Two runs of one configuration with isolate_state=false. At the second
        # run's start, the priority callbacks and the sensors' first sample see
        # the hinge at its configured angle, as at the first run's start.
        hinge = PanelHinge(axis=SVector(0.0, 1.0, 0.0), max_rate_rad_s=0.1)
        args = runtime_args((HoldHinge(),); sensor=HingeSensor(hinge), guidance=(), navigation=())
        sc = only(args.dynamics_model.spacecraft)
        panel = SM.Link(m=2.0, r=MVector(0.0, 1.5, 0.0), dims=MVector(0.02, 2.0, 1.0))
        push!(sc.links, panel)
        add_actuator!(sc, panel, hinge)
        for run in 1:2
            empty!(SEEN)
            run_simulation(args; isolate_state=false, priority_callback_factories=(a -> hinge_seer(hinge),))
            @testset let run = run
                @test [(who, angle) for (who, at, angle) in SEEN if at == 0.0] == [(:priority, 0.0), (:sensor, 0.0)]
                @test hinge.angle_rad > 0.2                                  # the run moved it
            end
        end
    end

    @testset "a resumed run restores the run state models keep in checkpoint files" begin
        # Navigation ticks every 0.3 s, off the 1 s checkpoint boundaries (the
        # next testset ticks on them).
        checkpointed(dir; kw...) = (checkpoint_enabled=true, checkpoint_interval_s=1.0, checkpoint_directory=dir, kw...)
        counted(mission_time, settings; navigation=TickCounter(), sensor=KeptCountSensor()) =
            runtime_args((); sensor=sensor, navigation=(navigation,), navigation_period=0.3, mission_time=mission_time,
                         settings=settings)
        run!(args) = (run_simulation(args; isolate_state=false); args)
        straight = run!(counted(6.0, (;)))
        dir = mktempdir()
        early = run!(counted(3.0, checkpointed(dir)))
        early_nav, early_sensor = only(early.navigation_model.navigation_effectors), only(SM.mounted_sensors(only(early.dynamics_model.spacecraft))).sensor
        saved = SM.IOSerialization._load_checkpoint(early).runtime_state.models
        @test saved.navigation == Any[(ticks=early_nav.ticks,)] && saved.sensors == [Any[(count=early_sensor.count,)]]
        @test saved.guidance == Any[nothing]                                  # LoggedGuidance keeps none
        resumed = run!(counted(6.0, checkpointed(dir; resume_from_checkpoint=true)))
        nav = only(resumed.navigation_model.navigation_effectors)
        sensor = only(SM.mounted_sensors(only(resumed.dynamics_model.spacecraft))).sensor
        straight_sensor = only(SM.mounted_sensors(only(straight.dynamics_model.spacecraft))).sensor
        @test nav.calls == [:initialize, :restore]                            # the reset first, then the saved state
        @test early_nav.ticks > 0 && nav.this_run > 0
        @test nav.ticks == early_nav.ticks + nav.this_run                     # counting on from the checkpoint
        @test early_sensor.count == 4 && straight_sensor.count == 7           # samples at 0..3 s and at 0..6 s
        @test sensor.count == straight_sensor.count                           # each sample of the run counted once
        @test sensor_measurement(sensor) == sensor_measurement(straight_sensor)
        @test isempty(sensor.initialized_at)                                  # no initialize_sensor! on a resume
        # Another model at the position is refused, and so is a checkpoint
        # holding no state for a model that keeps some.
        msg = refusal(() -> run!(counted(6.0, checkpointed(dir; resume_from_checkpoint=true); navigation=LoggedNavigation())))
        @test msg !== nothing && occursin("differ from the checkpoint", msg)
        plain = mktempdir()
        plain_args = run!(counted(3.0, checkpointed(plain); navigation=LoggedNavigation(), sensor=CountingSensor()))
        @test !hasproperty(SM.IOSerialization._load_checkpoint(plain_args).runtime_state, :models)   # nothing kept: no entry
        msg = refusal(() -> run!(counted(6.0, checkpointed(plain; resume_from_checkpoint=true); sensor=CountingSensor())))
        @test msg !== nothing && occursin("The run state of TickCounter", msg)
        # A snapshot written by another model (the sensor) holds nothing for a
        # model of the same kind that keeps state only on the resume.
        mixed = mktempdir()
        run!(counted(3.0, checkpointed(mixed); navigation=MaybeKept(false, 0)))
        @test SM.IOSerialization._load_checkpoint(counted(3.0, checkpointed(mixed))).runtime_state.models.navigation == Any[nothing]
        msg = refusal(() -> run!(counted(6.0, checkpointed(mixed; resume_from_checkpoint=true); navigation=MaybeKept(true, 0))))
        @test msg !== nothing && occursin("The run state of MaybeKept", msg)
        @test refusal(() -> run!(counted(6.0, checkpointed(mixed; resume_from_checkpoint=true); navigation=MaybeKept(false, 0)))) === nothing
        # A model that saves state it cannot restore is refused when the run
        # checkpoints, and runs otherwise.
        msg = refusal(() -> run!(counted(3.0, checkpointed(mktempdir()); navigation=HalfKept())))
        @test msg !== nothing && occursin("has no restore_checkpoint_state!", msg)
        @test refusal(() -> run!(counted(3.0, (;); navigation=HalfKept()))) === nothing
    end

    @testset "a resumed run restores the engine's state before the models start, and theirs after their initializers" begin
        checkpointed(dir; kw...) = (checkpoint_enabled=true, checkpoint_interval_s=1.0, checkpoint_directory=dir, kw...)
        # Navigation every 0.4 s has no tick at the 3 s checkpoint, so none is
        # pending there for the resumed run to take before its callbacks initialize.
        probed(mission_time, settings) = runtime_args((); sensor=KeptCountSensor(), navigation=(ResumeProbe(0),),
            navigation_period=0.4, mission_time=mission_time, settings=settings)
        dir = mktempdir()
        early = probed(3.0, checkpointed(dir))
        run_simulation(early; isolate_state=false)
        saved = SM.IOSerialization._load_checkpoint(early).runtime_state
        ticks = only(saved.models.navigation).ticks
        latest = only(only(saved.callbacks.spacecraft).sensors).latest_time_s
        @test ticks > 0 && latest == 3.0
        empty!(RESUME_LOG)
        run_simulation(probed(6.0, checkpointed(dir; resume_from_checkpoint=true)); isolate_state=false,
            extra_callback_factories=(resume_probe_callback,))
        # The initializer already sees both resume flags: the engine's and the
        # callbacks' state came back first. The model's own state comes back after
        # its initializer, and the callbacks then initialize with both.
        @test RESUME_LOG[1:3] == [(:initialize, true, true), (:restore, ticks), (:callbacks, ticks, latest)]
    end

    @testset "a resumed run ticks as a run in one piece, at checkpoint boundaries too" begin
        # Navigation ticks every 0.25 s: the 1 s checkpoint boundaries and the
        # 3 s resume time fall on its grid.
        checkpointed(dir; kw...) = (checkpoint_enabled=true, checkpoint_interval_s=1.0, checkpoint_directory=dir, kw...)
        counted(mission_time, settings) = runtime_args((); sensor=KeptCountSensor(), navigation=(TickCounter(),),
            navigation_period=0.25, mission_time=mission_time, settings=settings)
        run!(args) = (run_simulation(args; isolate_state=false); args)
        nav(args) = only(args.navigation_model.navigation_effectors)
        straight = run!(counted(6.0, (;)))
        split = run!(counted(6.0, checkpointed(mktempdir())))
        dir = mktempdir()
        early = run!(counted(3.0, checkpointed(dir)))
        resumed = run!(counted(6.0, checkpointed(dir; resume_from_checkpoint=true)))
        @test nav(straight).ticks == 23                                    # 0.25 s to 5.75 s
        @test nav(split).ticks == nav(straight).ticks                      # with the ticks at 1 s to 5 s
        @test nav(early).ticks == 11                                       # none at the run's end
        @test nav(resumed).this_run == 12 && nav(resumed).ticks == nav(straight).ticks   # 3 s to 5.75 s
    end

    @testset "an unknown schedule is refused" begin
        err = try
            run_simulation(runtime_args((OddControl(),)))
            nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin("tick_schedule", sprint(showerror, err))
    end
end
end
