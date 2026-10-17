module ConformanceNegativeControlsTests

# Each check of SpaceAGORA.Conformance gets at least one broken input that it must
# fail, so that a check passing a defective model is noticed. The file fails when a
# check id has no control.

using Test, SpaceAGORA, Logging
using StaticArrays: SVector
using SpaceAGORA.Conformance: CHECK_IDS, check_force_torque_model, check_service_facade, check_sensor, conformance_fixture,
    passed, check_guidance_model, check_navigation_model, check_control_model
using Random: Random

const SM = SpaceAGORA.SimulationModel
const C = SpaceAGORA.Conformance
const Z3 = SVector(0.0, 0.0, 0.0)
const ForceModel = SpaceAGORA.AbstractForceTorqueModel
const Sample = SpaceAGORA.StateSample
const Env = SpaceAGORA.EnvironmentSample

along_velocity(x, newtons) = (newtons * x.vel_ii / sqrt(sum(abs2, x.vel_ii)), Z3)

# A clean model: a small constant drag-like force along the velocity.
struct NcClean <: ForceModel end
SpaceAGORA.wrench(::NcClean, x::Sample, env::Env, t::Float64) = along_velocity(x, -1.0e-3)

# hooks_present: neither wrench nor the deprecated calcForceTorque.
struct NcNoHook <: ForceModel end

# supported_sample_signature: a method on a concrete StateSample{Nothing}, which
# only samples built outside a run reach.
struct NcConcreteSample <: ForceModel end
SpaceAGORA.wrench(::NcConcreteSample, x::Sample, env::Env, t::Float64) = along_velocity(x, -1.0e-3)
SpaceAGORA.wrench(::NcConcreteSample, x::Sample{Nothing}, env::Env, t::Float64) = (Z3, Z3)

# repeatable_run: state kept outside the model, advanced at every run start.
const NC_RUN_EPOCH = Ref(0)
struct NcRunEpoch <: ForceModel end
function SpaceAGORA.wrench(::NcRunEpoch, x::Sample, env::Env, t::Float64)
    t == 0.0 && (NC_RUN_EPOCH[] += 1; return (Z3, Z3))
    return along_velocity(x, 1.0e-3 * NC_RUN_EPOCH[])
end

# deepcopy_equivalent: a gain registered under the identity of the model's own
# field, which a deep copy does not carry.
const NC_GAINS = IdDict{Any, Float64}()
struct NcIdentityGain <: ForceModel
    key::Base.RefValue{Int}
end
NcIdentityGain() = (model = NcIdentityGain(Ref(0)); NC_GAINS[model.key] = 1.0e-3; model)
SpaceAGORA.wrench(m::NcIdentityGain, x::Sample, env::Env, t::Float64) = along_velocity(x, get(NC_GAINS, m.key, 0.0))

# global_rng_independent: a draw from the global generator.
struct NcGlobalRandom <: ForceModel end
SpaceAGORA.wrench(::NcGlobalRandom, x::Sample, env::Env, t::Float64) = along_velocity(x, 1.0e-3 * rand())

# concurrent_equals_sequential: scratch shared by every instance, reset at the
# start of a run, so that a run alone is repeatable but two overlapping runs mix.
const NC_SCRATCH = Ref(0.0)
struct NcSharedScratch <: ForceModel end
function SpaceAGORA.wrench(::NcSharedScratch, x::Sample, env::Env, t::Float64)
    t == 0.0 && (NC_SCRATCH[] = 0.0)
    NC_SCRATCH[] += 1.0
    return along_velocity(x, 1.0e-6 * NC_SCRATCH[])
end

# output_hook_pure: a call counter in the model, and a log message.
mutable struct NcCountsCalls <: ForceModel
    calls::Int
end
SpaceAGORA.wrench(m::NcCountsCalls, x::Sample, env::Env, t::Float64) = (m.calls += 1; along_velocity(x, -1.0e-3))
struct NcLogs <: ForceModel end
SpaceAGORA.wrench(::NcLogs, x::Sample, env::Env, t::Float64) = (@info "NcLogs wrench" maxlog = 1; along_velocity(x, -1.0e-3))

# output_hook_pure: two fields share one vector, and each call replaces one of them
# with an equal copy, so the values stay the same but the aliasing does not.
mutable struct NcBreaksAliasing <: ForceModel
    a::Vector{Float64}
    b::Vector{Float64}
end
NcBreaksAliasing() = (v = [1.0e-3]; NcBreaksAliasing(v, v))
SpaceAGORA.wrench(m::NcBreaksAliasing, x::Sample, env::Env, t::Float64) = (m.b = copy(m.b); along_velocity(x, -m.a[1]))

# concurrent_equals_sequential: the first call after the token is set sleeps, so one
# of two concurrent runs stalls while the other waits for it at the rendezvous.
const NC_STALL = Threads.Atomic{Bool}(false)
struct NcStallsOnce <: ForceModel end
function SpaceAGORA.wrench(::NcStallsOnce, x::Sample, env::Env, t::Float64)
    t > 0.0 && Threads.atomic_cas!(NC_STALL, true, false) && sleep(4.0)
    return along_velocity(x, -1.0e-3)
end

# environment_declared: reads the atmosphere without declaring it, and returns a
# non-finite force when the undeclared solar sample is missing.
struct NcUndeclaredAtmosphere <: ForceModel end
SpaceAGORA.wrench(::NcUndeclaredAtmosphere, x::Sample, env::Env, t::Float64) =
    along_velocity(x, -env.atmosphere.rho_kg_m3)
struct NcUndeclaredSolar <: ForceModel end
SpaceAGORA.wrench(::NcUndeclaredSolar, x::Sample, env::Env, t::Float64) =
    env.solar === nothing ? (SVector(NaN, NaN, NaN), Z3) : along_velocity(x, 1.0e-6)

# environment_declared, third bodies: declares the Moon and reads a second body.
# Needs SPICE ephemerides, so it runs when the kernels are present.
struct NcUndeclaredThirdBody <: ForceModel end
SpaceAGORA.environment_requirements(::NcUndeclaredThirdBody) =
    SpaceAGORA.EffectorEnvironmentRequirements(third_body_names=("Moon",))
SpaceAGORA.wrench(::NcUndeclaredThirdBody, x::Sample, env::Env, t::Float64) =
    along_velocity(x, 1.0e-20 * sum(abs2, env.third_bodies.positions_ii[2]))

include(joinpath(@__DIR__, "..", "..", "fixtures", "conformance_instances.jl"))
using .ConformanceInstances: on_spice, spice_available

# Sensors: one record shape, a 3-component measurement every second, and a counter
# some controls keep. The clean one adds noise from its own stream to the body rate.
for name in (:NcCleanSensor, :NcShortSensor, :NcNaNSensor, :NcRunEpochSensor, :NcIdentitySensor, :NcGlobalRandSensor,
        :NcSharedSensor, :NcFixedSeedSensor, :NcCountingSensor, :NcAtmosphereSensor)
    @eval begin
        mutable struct $name <: SpaceAGORA.AbstractSensor
            period_s::Float64
            seed::Int
            rng::Any
            latest::Any
            latest_time_s::Float64
            count::Int
        end
        $name() = $name(1.0, 1, nothing, SVector(NaN, NaN, NaN), NaN, 0)
        SpaceAGORA.sensor_dimension(::$name) = 3
    end
end
noisy(truth, rng) = truth.ω_body + 1.0e-4 * SVector(randn(rng), randn(rng), randn(rng))
SpaceAGORA.sensor_measure(::NcCleanSensor, truth, ctx, rng) = noisy(truth, rng)
# measurement_shape: two components where it declares three; a NaN component.
SpaceAGORA.sensor_measure(::NcShortSensor, truth, ctx, rng) = SVector(truth.ω_body[1], truth.ω_body[2])
SpaceAGORA.sensor_measure(::NcNaNSensor, truth, ctx, rng) = SVector(NaN, truth.ω_body[2], truth.ω_body[3])
# repeatable_run: state outside the sensor, advanced at every run start.
const NC_SENSOR_EPOCH = Ref(0)
SpaceAGORA.initialize_sensor!(::NcRunEpochSensor, truth, ctx, rng) = (NC_SENSOR_EPOCH[] += 1; nothing)
SpaceAGORA.sensor_measure(::NcRunEpochSensor, truth, ctx, rng) = noisy(truth, rng) .+ 1.0e-3 * NC_SENSOR_EPOCH[]
# deepcopy_equivalent: an offset registered under the identity of the sensor.
const NC_SENSOR_OFFSETS = IdDict{Any, Float64}()
NcRegisteredSensor() = (s = NcIdentitySensor(); NC_SENSOR_OFFSETS[s] = 1.0e-3; s)
SpaceAGORA.sensor_measure(s::NcIdentitySensor, truth, ctx, rng) = noisy(truth, rng) .+ get(NC_SENSOR_OFFSETS, s, 0.0)
# global_rng_independent: noise from the global generator.
SpaceAGORA.sensor_measure(::NcGlobalRandSensor, truth, ctx, rng) = truth.ω_body .+ 1.0e-4 * rand()
# concurrent_equals_sequential: scratch shared by every instance, reset at a run's start.
const NC_SENSOR_SCRATCH = Ref(0.0)
function SpaceAGORA.sensor_measure(::NcSharedSensor, truth, ctx, rng)
    truth.t_s == 0.0 && (NC_SENSOR_SCRATCH[] = 0.0)
    NC_SENSOR_SCRATCH[] += 1.0
    return noisy(truth, rng) .+ 1.0e-6 * NC_SENSOR_SCRATCH[]
end
# seed_repeatable: noise from its stream, but the stream is reset to a fixed seed.
SpaceAGORA.initialize_sensor!(::NcFixedSeedSensor, truth, ctx, rng) = (Random.seed!(rng, 1234); nothing)
SpaceAGORA.sensor_measure(::NcFixedSeedSensor, truth, ctx, rng) = noisy(truth, rng)
# custom_state_reset: a sample counter that initialize_sensor! does not reset.
SpaceAGORA.sensor_measure(s::NcCountingSensor, truth, ctx, rng) = (s.count += 1; noisy(truth, rng) .+ 1.0e-6 * s.count)
# requirements_supported: asks for the atmosphere.
SpaceAGORA.sensor_requirements(::NcAtmosphereSensor) = SpaceAGORA.ActuatorRequirements(atmosphere=true)
SpaceAGORA.sensor_measure(::NcAtmosphereSensor, truth, ctx, rng) = noisy(truth, rng)
# run_state_fields: not mutable; without latest_time_s.
struct NcImmutableSensor <: SpaceAGORA.AbstractSensor
    period_s::Float64
    seed::Int
    rng::Any
    latest::Any
    latest_time_s::Float64
end
NcImmutableSensor() = NcImmutableSensor(1.0, 1, nothing, SVector(NaN, NaN, NaN), NaN)
SpaceAGORA.sensor_dimension(::NcImmutableSensor) = 3
SpaceAGORA.sensor_measure(::NcImmutableSensor, truth, ctx, rng) = noisy(truth, rng)
mutable struct NcNoTimeSensor <: SpaceAGORA.AbstractSensor
    period_s::Float64
    seed::Int
    rng::Any
    latest::Any
end
NcNoTimeSensor() = NcNoTimeSensor(1.0, 1, nothing, SVector(NaN, NaN, NaN))
SpaceAGORA.sensor_dimension(::NcNoTimeSensor) = 3
SpaceAGORA.sensor_measure(::NcNoTimeSensor, truth, ctx, rng) = noisy(truth, rng)
# Clean, and it needs initialize_sensor!: its offset is NaN until that sets it.
mutable struct NcNeedsInitSensor <: SpaceAGORA.AbstractSensor
    period_s::Float64
    seed::Int
    rng::Any
    latest::Any
    latest_time_s::Float64
    offset::Float64
end
NcNeedsInitSensor() = NcNeedsInitSensor(1.0, 1, nothing, SVector(NaN, NaN, NaN), NaN, NaN)
SpaceAGORA.sensor_dimension(::NcNeedsInitSensor) = 3
SpaceAGORA.initialize_sensor!(s::NcNeedsInitSensor, truth, ctx, rng) = (s.offset = 0.0; nothing)
SpaceAGORA.sensor_measure(s::NcNeedsInitSensor, truth, ctx, rng) = noisy(truth, rng) .+ s.offset
# seed_repeatable: like a gyro whose bias walks, its first reading draws nothing (the
# current bias) and later ones do; initialize_sensor! resets the stream to a fixed seed.
mutable struct NcFixedSeedWalkSensor <: SpaceAGORA.AbstractSensor
    period_s::Float64
    seed::Int
    rng::Any
    latest::Any
    latest_time_s::Float64
    bias::SVector{3, Float64}
end
NcFixedSeedWalkSensor() = NcFixedSeedWalkSensor(1.0, 1, nothing, SVector(NaN, NaN, NaN), NaN, SVector(0.0, 0.0, 0.0))
SpaceAGORA.sensor_dimension(::NcFixedSeedWalkSensor) = 3
SpaceAGORA.initialize_sensor!(s::NcFixedSeedWalkSensor, truth, ctx, rng) =
    (Random.seed!(rng, 1234); s.bias = SVector(0.0, 0.0, 0.0); nothing)
function SpaceAGORA.sensor_measure(s::NcFixedSeedWalkSensor, truth, ctx, rng)
    reading = truth.ω_body + s.bias
    s.bias += 1.0e-4 * SVector(randn(rng), randn(rng), randn(rng))
    return reading
end

const SENSOR_RUN_LEVEL = (:repeatable_run, :deepcopy_equivalent, :global_rng_independent, :concurrent_equals_sequential,
    :seed_repeatable, :custom_state_reset)
sensor_except(ids...) = filter(id -> !(id in ids), CHECK_IDS.sensor)

# no_deprecated_hooks (force/torque): a calcForceTorque method of its own beside wrench.
struct NcDeprecatedForce <: ForceModel end
SpaceAGORA.wrench(::NcDeprecatedForce, x::Sample, env::Env, t::Float64) = along_velocity(x, -1.0e-3)
SpaceAGORA.calcForceTorque(::NcDeprecatedForce, x, p, i) = (Z3, Z3)

# Guidance and navigation models: one record shape and one behaviour per control.
# `acc` is the model's own state, which the run-level checks compare at the end of a
# run; `ticked` marks a run's first tick; `buffer` is what bind_spacecraft shares.
const NC_GNC_EPOCH = Ref(0)
const NC_GNC_SCRATCH = Ref(0.0)
const NC_GNC_GAINS = IdDict{Any, Float64}()
gnc_step(::Val{:Clean}, m, u, i) = 1.0e-9 * u.sc[i].pos[1]
gnc_step(::Val{:RunEpoch}, m, u, i) = (m.ticked || (m.ticked = true; NC_GNC_EPOCH[] += 1); Float64(NC_GNC_EPOCH[]))
gnc_step(::Val{:Identity}, m, u, i) = get(NC_GNC_GAINS, m, 0.0) + 1.0e-9 * u.sc[i].pos[1]
gnc_step(::Val{:Rand}, m, u, i) = rand()
gnc_step(::Val{:Shared}, m, u, i) = (m.ticked || (m.ticked = true; NC_GNC_SCRATCH[] = 0.0); NC_GNC_SCRATCH[] += 1.0)
gnc_step(::Val{:SharedBind}, m, u, i) = (m.buffer[1] += 1.0; m.buffer[1])
gnc_step(::Val{:NoReset}, m, u, i) = 1.0
for (prefix, Abstract, tick, initializer) in ((:NcG, SpaceAGORA.AbstractGuidanceModel, :update_guidance!, :initialize_guidance!),
        (:NcN, SpaceAGORA.AbstractNavigationModel, :update_navigation!, :initialize_navigation!))
    for variant in (:Clean, :NoTick, :RunEpoch, :Identity, :Rand, :Shared, :SharedBind, :NoReset)
        name = Symbol(prefix, variant)
        @eval begin
            mutable struct $name <: $Abstract
                acc::Float64
                ticked::Bool
                buffer::Vector{Float64}
            end
            $name() = $name(0.0, false, [0.0])
        end
        variant === :NoTick && continue
        @eval SpaceAGORA.$tick(m::$name, u, p, t::Float64, i::Int) = (m.acc += gnc_step(Val($(QuoteNode(variant))), m, u, i); nothing)
        # The models reset their own state at a run's start, except NoReset.
        variant !== :NoReset &&
            @eval SpaceAGORA.SimulationLifecycle.$initializer(m::$name, u, p, t) = (m.acc = 0.0; m.ticked = false; nothing)
        variant === :SharedBind && @eval SpaceAGORA.bind_spacecraft(m::$name, i::Int) = $name(m.acc, m.ticked, m.buffer)
    end
end
registered(model) = (NC_GNC_GAINS[model] = 1.0e-3; model)
# Clean, and its bound copies need the fixture's mounted sensor: each tick reads it.
mutable struct NcNSensorReader <: SpaceAGORA.AbstractNavigationModel
    acc::Float64
end
SpaceAGORA.update_navigation!(m::NcNSensorReader, u, p, t::Float64, i::Int) =
    (m.acc += length(SpaceAGORA.sensor_measurement(p.args.dynamics_model.spacecraft[i], 1)); nothing)
SpaceAGORA.SimulationLifecycle.initialize_navigation!(m::NcNSensorReader, u, p, t) = (m.acc = 0.0; nothing)
SpaceAGORA.bind_spacecraft(m::NcNSensorReader, i::Int) = NcNSensorReader(m.acc)

# Control models: the tick sets `cmd`, and the output hooks read it.
mutable struct NcCtl{V} <: SpaceAGORA.AbstractControlEffectorModel
    acc::Float64
    ticked::Bool
    buffer::Vector{Float64}
    calls::Int
end
NcCtl{V}() where {V} = NcCtl{V}(0.0, false, [0.0], 0)
const NC_CTL_STEPS = (:Clean, :RunEpoch, :Identity, :Rand, :Shared, :SharedBind, :NoReset)
for variant in (:Clean, :RunEpoch, :Identity, :Rand, :Shared, :SharedBind, :NoReset, :NoWrench, :Deprecated, :Narrow,
        :Impure, :Negative)
    step = variant in NC_CTL_STEPS ? variant : :Clean
    @eval SpaceAGORA.update_control!(m::NcCtl{$(QuoteNode(variant))}, u, p, t::Float64, i::Int) =
        (m.acc += gnc_step(Val($(QuoteNode(step))), m, u, i); nothing)
    variant === :NoWrench && continue
    @eval SpaceAGORA.control_wrench(m::NcCtl{$(QuoteNode(variant))}, x::Sample, p, t::Float64) =
        (SVector(1.0e-4 * m.acc / (1 + abs(m.acc)), 0.0, 0.0), Z3)
    variant === :Negative || @eval SpaceAGORA.control_propellant_rate(m::NcCtl{$(QuoteNode(variant))}, x::Sample, p, t::Float64) =
        1.0e-6 * abs(m.acc) / (1 + abs(m.acc))
end
# The models reset their own state at a run's start, except NoReset.
SpaceAGORA.SimulationLifecycle.initialize_control!(m::NcCtl, u, p, t) = (m.acc = 0.0; m.ticked = false; m.calls = 0; nothing)
SpaceAGORA.SimulationLifecycle.initialize_control!(m::NcCtl{:NoReset}, u, p, t) = nothing
SpaceAGORA.bind_spacecraft(m::NcCtl{:SharedBind}, i::Int) = NcCtl{:SharedBind}(m.acc, m.ticked, m.buffer, m.calls)
SpaceAGORA.calcControlForceTorque(::NcCtl{:Deprecated}, u, p, i, t) = (Z3, Z3)
SpaceAGORA.control_wrench(m::NcCtl{:Narrow}, x::Sample{Nothing}, p, t::Float64) = (Z3, Z3)
SpaceAGORA.control_wrench(m::NcCtl{:Impure}, x::Sample, p, t::Float64) = (m.calls += 1; (Z3, Z3))
SpaceAGORA.control_propellant_rate(m::NcCtl{:Negative}, x::Sample, p, t::Float64) = -1.0e-5

gnc_except(kind, ids...) = filter(id -> !(id in ids), CHECK_IDS[kind])
const GNC_RUN_LEVEL = (:repeatable_run, :deepcopy_equivalent, :global_rng_independent, :concurrent_equals_sequential,
    :initialize_resets)
function gnc_controls()
    controls = Any[]
    for (kind, prefix, check) in ((:guidance_model, :NcG, check_guidance_model), (:navigation_model, :NcN, check_navigation_model))
        model(variant) = getfield(@__MODULE__, Symbol(prefix, variant))()
        initializer = kind === :navigation_model ? "initialize_navigation!" : "initialize_guidance!"
        push!(controls,
            (kind=kind, id=:hooks_present, broken="no tick method", expect="has no method",
                report=() -> check(model(:NoTick)), passing=(:bind_fresh_state,)),
            (kind=kind, id=:repeatable_run, broken="state outside the model advanced at every run start",
                expect="two runs of the fixture with isolate_state=true differ: ", report=() -> check(model(:RunEpoch)),
                passing=(:hooks_present, :bind_fresh_state)),
            (kind=kind, id=:deepcopy_equivalent, broken="a gain keyed by the model's identity",
                expect="and a run of its deep copy (isolate_state=true) differ: ",
                report=() -> check(registered(model(:Identity))), passing=gnc_except(kind, :deepcopy_equivalent)),
            (kind=kind, id=:global_rng_independent, broken="rand() in the tick", expect="seeded 1 and 2 differ: ",
                report=() -> check(model(:Rand)), passing=(:hooks_present, :bind_fresh_state)),
            (kind=kind, id=:concurrent_equals_sequential, broken="scratch shared by every instance",
                expect="and sequential run ", report=() -> check(model(:Shared)),
                passing=gnc_except(kind, :concurrent_equals_sequential)),
            # Its shared buffer is never reset either.
            (kind=kind, id=:bind_fresh_state, broken="bind_spacecraft copies share a buffer the tick changes",
                expect="changed during the run", report=() -> check(model(:SharedBind)),
                passing=gnc_except(kind, :bind_fresh_state, :initialize_resets)),
            (kind=kind, id=:initialize_resets, broken="its own state never reset in $(initializer)",
                expect="reset the model's run state in $(initializer)",
                report=() -> check(model(:NoReset)), passing=gnc_except(kind, :initialize_resets)))
    end
    ctl(variant) = NcCtl{variant}()
    kind = :control_model
    push!(controls,
        (kind=kind, id=:hooks_present, broken="a tick but no control_wrench", expect="has no method",
            report=() -> check_control_model(ctl(:NoWrench)), passing=(:no_deprecated_hooks, :supported_sample_signature)),
        (kind=kind, id=:no_deprecated_hooks, broken="a calcControlForceTorque method beside control_wrench",
            expect="implements both control_wrench and calcControlForceTorque",
            report=() -> check_control_model(ctl(:Deprecated)), passing=(:hooks_present, :supported_sample_signature)),
        (kind=kind, id=:supported_sample_signature, broken="a control_wrench method on StateSample{Nothing}",
            expect="StateSample{Nothing}, narrower than StateSample",
            report=() -> check_control_model(ctl(:Narrow)), passing=gnc_except(kind, :supported_sample_signature)),
        (kind=kind, id=:repeatable_run, broken="state outside the model advanced at every run start",
            expect="two runs of the fixture with isolate_state=true differ: ",
            report=() -> check_control_model(ctl(:RunEpoch)), passing=(:hooks_present, :no_deprecated_hooks)),
        (kind=kind, id=:deepcopy_equivalent, broken="a gain keyed by the model's identity",
            expect="and a run of its deep copy (isolate_state=true) differ: ",
            report=() -> check_control_model(registered(ctl(:Identity))), passing=gnc_except(kind, :deepcopy_equivalent)),
        (kind=kind, id=:global_rng_independent, broken="rand() in the tick", expect="seeded 1 and 2 differ: ",
            report=() -> check_control_model(ctl(:Rand)), passing=(:hooks_present, :no_deprecated_hooks)),
        (kind=kind, id=:concurrent_equals_sequential, broken="scratch shared by every instance", expect="and sequential run ",
            report=() -> check_control_model(ctl(:Shared)), passing=gnc_except(kind, :concurrent_equals_sequential)),
        (kind=kind, id=:output_hook_pure, broken="control_wrench counts its calls in the model",
            expect="calls of control_wrench at t = ", report=() -> check_control_model(ctl(:Impure)),
            passing=gnc_except(kind, :output_hook_pure)),
        (kind=kind, id=:positive_consumption, broken="a negative propellant rate",
            expect="control_propellant_rate returned -1.0e-5",
            report=() -> check_control_model(ctl(:Negative)), passing=gnc_except(kind, :positive_consumption)),
        (kind=kind, id=:initialize_resets, broken="its own state never reset in initialize_control!",
            expect="reset the model's run state in initialize_control!",
            report=() -> check_control_model(ctl(:NoReset)), passing=gnc_except(kind, :initialize_resets)),
        # Its shared buffer is never reset either.
        (kind=kind, id=:bind_fresh_state, broken="bind_spacecraft copies share a buffer the tick changes",
            expect="changed during the run", report=() -> check_control_model(ctl(:SharedBind)),
            passing=gnc_except(kind, :bind_fresh_state, :initialize_resets)),
        (kind=:force_torque_model, id=:no_deprecated_hooks, broken="a calcForceTorque method beside wrench",
            expect="implements both wrench and calcForceTorque",
            report=() -> check_force_torque_model(NcDeprecatedForce()), passing=(:hooks_present, :supported_sample_signature)))
    return controls
end

# Planning algorithms: RPO planners (a core family) with one defect each, valid for
# an unchanged goal like the example hold planner; `state` counts proposals.
const RP = SpaceAGORA.RPOPlannerInterfaces
const GA = SpaceAGORA.GNCAlgorithms
struct NcPlanner{D} <: RP.AbstractRPOPlanner end
const NC_PLANNER_SHARED = Ref(0)
const NC_PLANNER_STATES = IdDict{Any, Nothing}()
SpaceAGORA.planner_capabilities(::NcPlanner{D}) where {D} =
    RP.RPOPlannerCapabilities(state_sources=(:truth,), frames=(:target_rtn,), retiming=D !== :SilentReplan)
SpaceAGORA.initialize_planner(::NcPlanner, context) = Ref(0)
SpaceAGORA.initialize_planner(::NcPlanner{:SharedState}, context) = NC_PLANNER_SHARED
SpaceAGORA.initialize_planner(::NcPlanner{:IdentityState}, context) = (state = Ref(0); NC_PLANNER_STATES[state] = nothing; state)
GA.consumes_rng(::NcPlanner{:DrawsDespiteDeclaring}) = false
RP.require_planner_support(::NcPlanner{:NeedsPackage}) = error("NcPlanner{:NeedsPackage} needs a package that is not installed")
function nc_plan(state, ::NcPlanner{D}, r::RP.RPOPlanningRequest, rng) where {D}
    state[] += 1
    D === :Mutates && (r.geometry.station.points_body .= 100.0)
    D === :IdentityState && !haskey(NC_PLANNER_STATES, state) &&
        return RP.RPOPlanningResult(request_id=r.request_id, status=:infeasible, termination=:copied_state)
    draw = D === :GlobalRand ? rand() : rand(rng)
    position = collect(r.x_rtn[1:3])
    reference = RP.RPOReference(t_ref_s=[0.0, r.reference_dt_s], r_ref_rtn_m=hcat(position, collect(r.goal_rtn_m)),
        v_ref_rtn_mps=zeros(3, 2), origin_time_s=r.time_s, valid_until_s=r.valid_until_s, chaser_id=r.chaser_id,
        target_id=r.target_id, geometry_revision=r.geometry_revision)
    return RP.RPOPlanningResult(request_id=D === :WrongId ? r.request_id + 100 : r.request_id, status=:candidate,
        termination=:completed, reference=reference, diagnostics=(draw=draw,))
end
SpaceAGORA.plan_rpo!(state, planner::NcPlanner, r::RP.RPOPlanningRequest, rng::Random.AbstractRNG) =
    nc_plan(state, planner, r, rng)
SpaceAGORA.retime_rpo!(state, planner::NcPlanner, r::RP.RPOPlanningRequest, active::RP.RPOReference,
    rng::Random.AbstractRNG) = nc_plan(state, planner, r, rng)
nc_hold_requests() = C.rpo_conformance_requests(goal_rtn_m=nothing)
nc_algorithm_report(planner; simulation=nothing) =
    C.check_algorithm(planner, nc_hold_requests(); trusted=C.rpo_trusted_queries(), simulation)
nc_simulation() = C.rpo_conformance_simulation(goal_rtn_m=(3.0, 0.0, 0.0), mission_time_s=0.1)

# An algorithm with no family.
struct NcFamilyless <: GA.AbstractGNCAlgorithm end

# Families of this file: the variant names the defect, the policy is `on_rejection`.
# A request above the trusted limit is refused.
struct NcRequest <: GA.AbstractAlgorithmRequest
    request_id::Int
    reason::Symbol
    time_s::Float64
    value::Float64
    scratch::Vector{Float64}
end
NcRequest(id, time_s, value) = NcRequest(id, :forced, time_s, value, [0.0])
struct NcFam{Variant, Policy} <: GA.AbstractAlgorithmFamily end
struct NcOtherFamily <: GA.AbstractAlgorithmFamily end
GA.on_rejection(::NcFam{V, Policy}) where {V, Policy} = Policy
function GA.validate_proposal(::NcFam{V}, r::NcRequest, result, trusted) where {V}
    V === :RandomValidation && rand()
    V === :MutatingValidation && (r.scratch[1] += 1.0)
    return result.product <= trusted.limit ? GA.ValidationResult(true, :within_limit) :
        GA.ValidationResult(false, :above_limit)
end
function GA.install!(link::Base.RefValue{Any}, ::NcFam{V}, product, r::NcRequest, t; validation, diagnostics) where {V}
    V === :SilentInstall || (link[] = (product=product, t=t))
    return nothing
end
GA.is_current(::NcFam, r::NcRequest, product, t) = true
GA.fallback_product(::NcFam{:Clean, :install_fallback}, r::NcRequest) = 0.0
C.conformance_link(::NcFam) = Ref{Any}(nothing)
struct NcDoubler{F} <: GA.AbstractGNCAlgorithm end
GA.algorithm_family(::NcDoubler{F}) where {F} = F()
GA.consumes_rng(::NcDoubler) = false
GA.propose!(state, ::NcDoubler, r::NcRequest, rng) =
    GA.AlgorithmResult(request_id=r.request_id, status=:candidate, termination=:completed, product=2r.value)
nc_requests() = [NcRequest(1, 0.0, 2.0), NcRequest(2, 1.0, 9.0), NcRequest(3, 2.0, 1.0)]
# Families with an obligation defect the walk itself cannot get past, or never reaches.
GA.algorithm_result_type(::NcFam{:NotAType}) = :not_a_type
GA.validate_proposal(::NcFam{:RejectAll}, r::NcRequest, result, trusted) = GA.ValidationResult(false, :never)
struct NcBareFamily <: GA.AbstractAlgorithmFamily end
GA.validate_proposal(::NcBareFamily, r::NcRequest, result, trusted) = GA.ValidationResult(false, :never)
# A family whose revisions the capability check refuses, and an algorithm that answers
# a refused revision `:unsupported` yet with a product.
struct NcLooseResult <: GA.AbstractAlgorithmResult
    request_id::Int
    status::Symbol
    termination::Symbol
    product::Any
    diagnostics::NamedTuple
end
GA.proposed_product(r::NcLooseResult) = r.product
GA.algorithm_result_type(::NcFam{:Loose}) = NcLooseResult
GA.is_revision_request(::NcFam{:Loose}, r::NcRequest) = r.reason === :revise
GA.validate_capabilities(::NcFam{:Loose}, algorithm, r::NcRequest; checkpoint::Bool=false) =
    r.reason === :revise ? GA.ValidationResult(false, :revision_unsupported) : GA.ValidationResult(true, :supported)
struct NcLooseAlgorithm <: GA.AbstractGNCAlgorithm end
GA.algorithm_family(::NcLooseAlgorithm) = NcFam{:Loose, :fail}()
GA.consumes_rng(::NcLooseAlgorithm) = false
GA.propose!(state, ::NcLooseAlgorithm, r::NcRequest, rng) =
    NcLooseResult(r.request_id, :candidate, :completed, 2r.value, NamedTuple())
GA.propose_revision!(state, ::NcLooseAlgorithm, r::NcRequest, active, rng) =
    NcLooseResult(r.request_id, :unsupported, :revision_not_supported, active, NamedTuple())
nc_loose_requests() = [NcRequest(1, 0.0, 2.0), NcRequest(2, :revise, 1.0, 1.0, [0.0])]
nc_family_report(F; family=F()) =
    C.check_algorithm_family(family, NcDoubler{F}(), nc_requests(); link=Ref{Any}(nothing), trusted=(limit=10.0,))

algorithm_except(ids...) = filter(id -> !(id in ids), CHECK_IDS.algorithm)
function algorithm_controls()
    kind = :algorithm
    return [
        (kind=kind, id=:family_resolves, broken="no algorithm_family method", expect="declares no algorithm family",
            report=() -> C.check_algorithm(NcFamilyless(), nc_requests()), passing=()),
        (kind=kind, id=:family_trusted, broken="a family defined outside core",
            expect="families are trusted and core-owned",
            report=() -> C.check_algorithm(NcDoubler{NcFam{:Clean, :fail}}(), nc_requests(); trusted=(limit=10.0,)),
            passing=algorithm_except(:family_trusted)),
        (kind=kind, id=:result_contract, broken="a result for another request id",
            expect="the result's request_id is 101, the request's 1",
            report=() -> nc_algorithm_report(NcPlanner{:WrongId}()), passing=algorithm_except(:result_contract)),
        (kind=kind, id=:stream_untouched_if_declared, broken="consumes_rng = false, yet it draws",
            expect="advanced the stream", report=() -> nc_algorithm_report(NcPlanner{:DrawsDespiteDeclaring}()),
            passing=algorithm_except(:stream_untouched_if_declared)),
        (kind=kind, id=:deterministic_given_stream, broken="rand() in its diagnostics",
            expect="with the global random number generator seeded 1 and 2",
            report=() -> nc_algorithm_report(NcPlanner{:GlobalRand}()), passing=algorithm_except(:deterministic_given_stream)),
        (kind=kind, id=:request_not_mutated, broken="it writes into the request's station points",
            expect="changed its copy of request(s) 1, 2, 3", report=() -> nc_algorithm_report(NcPlanner{:Mutates}()),
            passing=algorithm_except(:request_not_mutated)),
        (kind=kind, id=:initialize_independent, broken="initialize_planner returns one shared counter",
            expect="share 1 mutable object(s)", report=() -> nc_algorithm_report(NcPlanner{:SharedState}()),
            passing=algorithm_except(:initialize_independent)),
        (kind=kind, id=:refuses_explicitly, broken="no retiming declared, yet retime_rpo! plans in full",
            expect="yet propose_revision! returned :candidate", report=() -> nc_algorithm_report(NcPlanner{:SilentReplan}()),
            passing=algorithm_except(:refuses_explicitly)),
        (kind=kind, id=:refuses_explicitly, broken="a refused revision answered :unsupported, yet with a product",
            expect="yet propose_revision! returned :unsupported with a product",
            report=() -> C.check_algorithm(NcLooseAlgorithm(), nc_loose_requests(); trusted=(limit=10.0,)),
            passing=algorithm_except(:family_trusted, :refuses_explicitly)),
        (kind=kind, id=:driver_round_trip, broken="it answers differently for a copy of its state",
            expect="yet the driver's run gave an installed product",
            report=() -> nc_algorithm_report(NcPlanner{:IdentityState}()), passing=algorithm_except(:driver_round_trip)),
        (kind=kind, id=:simulation_run, broken="require_planner_support throws", expect="needs a package that is not installed",
            report=() -> nc_algorithm_report(NcPlanner{:NeedsPackage}(); simulation=nc_simulation()),
            passing=algorithm_except(:simulation_run)),
        (kind=:algorithm_family, id=:family_resolves, broken="the algorithm declares another family",
            expect="the request driver refuses the pair",
            report=() -> nc_family_report(NcFam{:Clean, :fail}; family=NcOtherFamily()), passing=()),
        (kind=:algorithm_family, id=:family_obligations, broken=":install_fallback without fallback_product",
            expect="fallback_product has no method", report=() -> nc_family_report(NcFam{:NoFallback, :install_fallback}),
            passing=(:family_resolves,)),
        (kind=:algorithm_family, id=:family_obligations, broken="validation draws from the global stream",
            expect="drew from the global random number generator",
            report=() -> nc_family_report(NcFam{:RandomValidation, :fail}), passing=(:family_resolves, :driver_round_trip)),
        (kind=:algorithm_family, id=:family_obligations, broken="validation writes into the request",
            expect="validate_proposal changed the request", report=() -> nc_family_report(NcFam{:MutatingValidation, :fail}),
            passing=(:family_resolves, :driver_round_trip)),
        (kind=:algorithm_family, id=:family_obligations, broken="algorithm_result_type returns a symbol, not a type",
            expect="algorithm_result_type returns :not_a_type, not a subtype of AbstractAlgorithmResult",
            report=() -> nc_family_report(NcFam{:NotAType, :fail}), passing=(:family_resolves,)),
        (kind=:algorithm_family, id=:family_obligations, broken="validation rejects every proposal",
            expect="validate_proposal accepted none of the proposals",
            report=() -> nc_family_report(NcFam{:RejectAll, :fail}), passing=(:family_resolves, :driver_round_trip)),
        (kind=:algorithm_family, id=:family_obligations, broken="no install! or is_current, and every proposal rejected",
            expect="install! has no method for a Base.RefValue{Any} link",
            report=() -> nc_family_report(NcBareFamily), passing=(:family_resolves, :driver_round_trip)),
        (kind=:algorithm_family, id=:driver_round_trip, broken="install! leaves the link unchanged",
            expect="install! left the link unchanged", report=() -> nc_family_report(NcFam{:SilentInstall, :fail}),
            passing=(:family_resolves, :family_obligations)),
    ]
end

# A facade with the shape the adapter standard requires; each control breaks one rule.
function control_facade(; version=v"1.0.0", require=:exact, implemented=(:serve,), consumed=(:Value,),
        omit::Union{Nothing, Symbol}=nothing)
    facade = Module(:NcFacade)
    check = require === :exact ? :(v == CONTRACT_VERSION ? nothing : throw(CompatibilityError(string(v)))) :
        require === :major ? :(v.major == CONTRACT_VERSION.major ? nothing : throw(CompatibilityError(string(v)))) :
        :(v == CONTRACT_VERSION ? nothing : throw(ArgumentError(string(v))))
    Core.eval(facade, quote
        struct CompatibilityError <: Exception
            msg::String
        end
        const CONTRACT_VERSION = $version
        require_version(v::VersionNumber) = $check
        serve() = 1
        const Value = 2.0
        const Empty = nothing
        const Nested = (; serve)
    end)
    omit === :IMPLEMENTED_FUNCTIONS || Core.eval(facade, :(const IMPLEMENTED_FUNCTIONS = $(QuoteNode(implemented))))
    omit === :CONSUMED_SERVICES || Core.eval(facade, :(const CONSUMED_SERVICES = $(QuoteNode(consumed))))
    return facade
end

# kind, check id, what is broken, a fragment of the failure detail it must give, the
# report, and the checks that must still pass.
const ALL_FORCE = CHECK_IDS.force_torque_model
except(ids...) = filter(id -> !(id in ids), ALL_FORCE)
const RUN_LEVEL = (:repeatable_run, :deepcopy_equivalent, :global_rng_independent, :concurrent_equals_sequential)

controls() = [
    (kind=:force_torque_model, id=:hooks_present, broken="no wrench and no calcForceTorque", expect="has no method wrench(",
        report=() -> check_force_torque_model(NcNoHook()), passing=(:supported_sample_signature,)),
    (kind=:force_torque_model, id=:supported_sample_signature, broken="a method on StateSample{Nothing}", expect="StateSample{Nothing}, narrower than StateSample",
        report=() -> check_force_torque_model(NcConcreteSample()), passing=except(:supported_sample_signature)),
    (kind=:force_torque_model, id=:repeatable_run, broken="state outside the model advanced at each run start", expect="two runs of the fixture with isolate_state=true differ: ",
        report=() -> check_force_torque_model(NcRunEpoch()),
        passing=(:hooks_present, :supported_sample_signature, :output_hook_pure, :environment_declared)),
    (kind=:force_torque_model, id=:deepcopy_equivalent, broken="a gain keyed by the identity of a field", expect="and a run of its deep copy (isolate_state=true) differ: ",
        report=() -> check_force_torque_model(NcIdentityGain()), passing=except(:deepcopy_equivalent)),
    (kind=:force_torque_model, id=:global_rng_independent, broken="rand() from the global generator", expect="seeded 1 and 2 differ: ",
        report=() -> check_force_torque_model(NcGlobalRandom()),
        passing=(:hooks_present, :supported_sample_signature, :environment_declared)),
    (kind=:force_torque_model, id=:concurrent_equals_sequential, broken="scratch shared by every instance", expect="and sequential run ",
        report=() -> check_force_torque_model(NcSharedScratch()),
        passing=(:hooks_present, :supported_sample_signature, :repeatable_run, :deepcopy_equivalent,
            :global_rng_independent, :environment_declared)),
    (kind=:force_torque_model, id=:output_hook_pure, broken="a call counter in the model", expect="the calls changed the model",
        report=() -> check_force_torque_model(NcCountsCalls(0)), passing=except(:output_hook_pure)),
    (kind=:force_torque_model, id=:output_hook_pure, broken="a log message", expect="the first: Info: NcLogs wrench",
        report=() -> check_force_torque_model(NcLogs()), passing=except(:output_hook_pure)),
    (kind=:force_torque_model, id=:output_hook_pure, broken="replaces one of two aliased fields with an equal copy",
        expect="the calls changed the model",
        report=() -> check_force_torque_model(NcBreaksAliasing()), passing=except(:output_hook_pure)),
    (kind=:force_torque_model, id=:environment_declared, broken="reads the atmosphere without declaring it", expect="wrench threw at t = 0.0 s with only the declared environment fields (planet_frame = nothing, atmosphere = nothing",
        report=() -> check_force_torque_model(NcUndeclaredAtmosphere()),
        passing=(:hooks_present, :supported_sample_signature)),
    (kind=:force_torque_model, id=:environment_declared, broken="NaN without the undeclared solar sample", expect="solar = nothing, third_bodies = nothing); it must return a finite force and torque",
        report=() -> check_force_torque_model(NcUndeclaredSolar()),
        passing=(:hooks_present, :supported_sample_signature, :output_hook_pure)),
    (spice_available() ? [(kind=:force_torque_model, id=:environment_declared,
        broken="reads a third body it did not declare", expect="BoundsError",
        report=() -> (entry = on_spice(planet -> NcUndeclaredThirdBody());
            check_force_torque_model(entry.model; fixture=entry.fixture)),
        passing=(:hooks_present, :supported_sample_signature))] : [])...,
    (kind=:sensor, id=:run_state_fields, broken="an immutable record", expect="is not a mutable struct",
        report=() -> check_sensor(NcImmutableSensor()), passing=()),
    (kind=:sensor, id=:run_state_fields, broken="no latest_time_s field", expect="lacks the field(s) latest_time_s",
        report=() -> check_sensor(NcNoTimeSensor()), passing=()),
    (kind=:sensor, id=:requirements_supported, broken="asks for the atmosphere", expect="add_sensor! refuses the sensor",
        report=() -> check_sensor(NcAtmosphereSensor()), passing=(:run_state_fields,)),
    (kind=:sensor, id=:measurement_shape, broken="two components where it declares three",
        expect="not 3 finite component(s)",
        report=() -> check_sensor(NcShortSensor()), passing=sensor_except(:measurement_shape)),
    (kind=:sensor, id=:measurement_shape, broken="a NaN component", expect="not 3 finite component(s)",
        report=() -> check_sensor(NcNaNSensor()), passing=sensor_except(:measurement_shape)),
    (kind=:sensor, id=:repeatable_run, broken="state outside the sensor advanced at every run start",
        expect="two runs of the fixture with isolate_state=true differ: ",
        report=() -> check_sensor(NcRunEpochSensor()), passing=(:run_state_fields, :requirements_supported, :measurement_shape)),
    (kind=:sensor, id=:deepcopy_equivalent, broken="an offset keyed by the sensor's identity",
        expect="and a run of its deep copy (isolate_state=true) differ: ",
        report=() -> check_sensor(NcRegisteredSensor()), passing=sensor_except(:deepcopy_equivalent)),
    (kind=:sensor, id=:global_rng_independent, broken="noise from the global generator", expect="seeded 1 and 2 differ: ",
        report=() -> check_sensor(NcGlobalRandSensor()), passing=(:run_state_fields, :requirements_supported, :measurement_shape)),
    (kind=:sensor, id=:concurrent_equals_sequential, broken="scratch shared by every instance", expect="and sequential run ",
        report=() -> check_sensor(NcSharedSensor()), passing=sensor_except(:concurrent_equals_sequential)),
    (kind=:sensor, id=:seed_repeatable, broken="its stream reset to a fixed seed",
        expect="yet runs with those seeds give the same measurements",
        report=() -> check_sensor(NcFixedSeedSensor()), passing=sensor_except(:seed_repeatable)),
    (kind=:sensor, id=:seed_repeatable, broken="its stream reset to a fixed seed, and its first reading draws nothing",
        expect="yet runs with those seeds give the same measurements",
        report=() -> check_sensor(NcFixedSeedWalkSensor()), passing=sensor_except(:seed_repeatable)),
    (kind=:sensor, id=:custom_state_reset, broken="a sample counter it never resets",
        expect="two runs with isolate_state=false on one configuration differ: ",
        report=() -> check_sensor(NcCountingSensor()), passing=sensor_except(:custom_state_reset)),
    gnc_controls()...,
    algorithm_controls()...,
    (kind=:service_facade, id=:facade_shape, broken="CONTRACT_VERSION other than the reviewed one", expect="CONTRACT_VERSION is v\"1.0.1\", not the reviewed 1.0.0",
        report=() -> facade_report(control_facade(version=v"1.0.1"); version=v"1.0.0"), passing=()),
    (kind=:service_facade, id=:facade_shape, broken="require_version accepts any minor or patch", expect="accepted a version other than 1.0.0",
        report=() -> facade_report(control_facade(require=:major); version=v"1.0.0"), passing=()),
    (kind=:service_facade, id=:facade_shape, broken="require_version throws another error", expect="threw ArgumentError, not CompatibilityError",
        report=() -> facade_report(control_facade(require=:argument); version=v"1.0.0"), passing=()),
    (kind=:service_facade, id=:facade_shape, broken="an implemented name listed twice", expect="implemented entries listed more than once: serve",
        report=() -> facade_report(control_facade(implemented=(:serve, :serve)); version=v"1.0.0"), passing=()),
    (kind=:service_facade, id=:facade_shape, broken="a name both implemented and consumed", expect="listed as both implemented and consumed: serve",
        report=() -> facade_report(control_facade(consumed=(:Value, :serve)); version=v"1.0.0"), passing=()),
    (kind=:service_facade, id=:facade_shape, broken="an implemented name that is not a function", expect="implemented Value is not a function",
        report=() -> facade_report(control_facade(implemented=(:Value,), consumed=()); version=v"1.0.0"), passing=()),
    (kind=:service_facade, id=:facade_shape, broken="an implemented path that does not resolve", expect="implemented Nested.absent does not resolve",
        report=() -> facade_report(control_facade(implemented=((:Nested, :absent),)); version=v"1.0.0"), passing=()),
    (kind=:service_facade, id=:facade_shape, broken="a consumed name bound to nothing", expect="consumed Empty is nothing",
        report=() -> facade_report(control_facade(consumed=(:Empty,)); version=v"1.0.0"), passing=()),
    (kind=:service_facade, id=:facade_shape, broken="IMPLEMENTED_FUNCTIONS a bare name, not a collection",
        expect="IMPLEMENTED_FUNCTIONS is :serve, not a collection of names or name paths",
        report=() -> facade_report(control_facade(implemented=:serve); version=v"1.0.0"), passing=()),
    (kind=:service_facade, id=:facade_shape, broken="CONSUMED_SERVICES nothing",
        expect="CONSUMED_SERVICES is nothing, not a collection of names or name paths",
        report=() -> facade_report(control_facade(consumed=nothing); version=v"1.0.0"), passing=()),
    (kind=:service_facade, id=:facade_shape, broken="no CONSUMED_SERVICES", expect="does not define CONSUMED_SERVICES",
        report=() -> facade_report(control_facade(omit=:CONSUMED_SERVICES); version=v"1.0.0"), passing=()),
]

# The facades are defined at run time; the check runs in the world that has them.
facade_report(facade; version) = Base.invokelatest(check_service_facade, facade; version)

status_of(report, id) = only(check.status for check in report.checks if check.id === id)
const MAY_SKIP = (:concurrent_equals_sequential, :bind_fresh_state, :stream_untouched_if_declared, :refuses_explicitly,
    :simulation_run)

@testset "Conformance negative controls" begin
    all_controls = controls()
    @testset "every check id has a negative control" begin
        for (kind, ids) in pairs(CHECK_IDS), id in ids
            @testset let kind = kind, id = id
                @test any(control -> control.kind === kind && control.id === id, all_controls)
            end
        end
    end
    @testset "clean inputs pass" begin
        report = check_force_torque_model(NcClean())
        @test [check.id for check in report.checks] == collect(CHECK_IDS.force_torque_model)
        @test all(check -> check.status === :pass, report.checks) || Threads.nthreads() < 2
        @test passed(report)
        facade = facade_report(control_facade(); version=v"1.0.0")
        @test [check.id for check in facade.checks] == collect(CHECK_IDS.service_facade)
        @test passed(facade)
        nested = facade_report(control_facade(implemented=(:serve, (:Nested, :serve))); version=v"1.0.0")
        @test passed(nested)
        for (report, kind) in ((check_guidance_model(NcGClean()), :guidance_model),
                (check_navigation_model(NcNClean()), :navigation_model), (check_control_model(NcCtl{:Clean}()), :control_model),
                (nc_algorithm_report(NcPlanner{:Clean}(); simulation=nc_simulation()), :algorithm),
                (nc_family_report(NcFam{:Clean, :fail}), :algorithm_family),
                (nc_family_report(NcFam{:Clean, :keep_previous}), :algorithm_family),
                (nc_family_report(NcFam{:Clean, :install_fallback}), :algorithm_family))
            @test [check.id for check in report.checks] == collect(CHECK_IDS[kind])
            @test all(check -> check.status in (:pass, :skip), report.checks)
            @test passed(report)
        end
        # A clean model is not answerable for another model of its fixture that uses a
        # deprecated hook.
        target = NcClean()
        mixed = check_force_torque_model(target; fixture=conformance_fixture(effectors=(target, NcDeprecatedForce())))
        @test status_of(mixed, :no_deprecated_hooks) === :pass
        controller = NcCtl{:Clean}()
        mixed = check_control_model(controller; fixture=conformance_fixture(control=(controller, NcCtl{:Deprecated}())))
        @test status_of(mixed, :no_deprecated_hooks) === :pass
        # The two-spacecraft bind check keeps the fixture's sensor, which the model reads.
        reader = check_navigation_model(NcNSensorReader(0.0); sensors=(NcCleanSensor(),))
        @test status_of(reader, :bind_fresh_state) === :pass
        @test passed(reader)
        # The second needs initialize_sensor! before it can measure, as the engine calls it.
        for sensor in (check_sensor(NcCleanSensor()), check_sensor(NcNeedsInitSensor()))
            @test [check.id for check in sensor.checks] == collect(CHECK_IDS.sensor)
            @test all(check -> check.status === :pass, sensor.checks) || Threads.nthreads() < 2
            @test passed(sensor)
        end
    end
    for control in all_controls
        # A check reports a broken input as a failing check; it never throws.
        report = try
            with_logger(NullLogger()) do
                control.report()
            end
        catch err
            err
        end
        expected = control.id === :concurrent_equals_sequential && Threads.nthreads() < 2 ? :skip : :fail
        @testset "$(control.id): $(control.broken)" begin
            @testset let id = control.id, thrown = report isa Exception ? sprint(showerror, report) : ""
                @test !(report isa Exception)
            end
            if !(report isa Exception)
                detail = only(check.detail for check in report.checks if check.id === control.id)
                @testset let id = control.id, status = status_of(report, control.id), detail = detail
                    @test status === expected
                    expected === :fail && @test occursin(control.expect, detail)
                    # A run-level control fails by its behaviour, not by a run that stopped.
                    control.id in (control.kind === :sensor ? SENSOR_RUN_LEVEL :
                        control.kind in (:guidance_model, :navigation_model, :control_model) ? GNC_RUN_LEVEL : RUN_LEVEL) &&
                        @test !occursin("run failed", detail)
                end
                expected === :fail && @test !passed(report)
                for id in control.passing
                    @testset let id = id, status = status_of(report, id)
                        # Skipped is not failed: concurrency needs two threads, a model without
                        # bind_spacecraft has no bound copies to check, and an algorithm check
                        # skips what the algorithm does not declare or the control does not give.
                        @test status === :pass || (id in MAY_SKIP && status === :skip)
                    end
                end
            end
        end
    end
    @testset "algorithm checks leave the global random number generator as they found it" begin
        # A validation that calls rand() is reported, and the generator is restored.
        for report in (() -> nc_family_report(NcFam{:RandomValidation, :fail}),
                       () -> C.check_algorithm(NcDoubler{NcFam{:RandomValidation, :fail}}(), nc_requests();
                           trusted=(limit=10.0,)))
            before = copy(Random.default_rng())
            report()
            @test rand(copy(before), UInt64) == rand(copy(Random.default_rng()), UInt64)
        end
    end
    @testset "algorithm checks report a defective family instead of throwing" begin
        # A family whose result type is not a type gets a failing report from either
        # check, not an exception from the comparisons that need the type.
        algorithm_report = C.check_algorithm(NcDoubler{NcFam{:NotAType, :fail}}(), nc_requests(); trusted=(limit=10.0,))
        @test !passed(algorithm_report)
        @test occursin("algorithm_result_type returns :not_a_type",
            only(check.detail for check in algorithm_report.checks if check.id === :result_contract))
        bare = nc_family_report(NcBareFamily)
        @test occursin("is_current has no method for this family",
            only(check.detail for check in bare.checks if check.id === :family_obligations))
    end
    @testset "concurrent_equals_sequential: a run stalls past the rendezvous limit" begin
        model = NcStallsOnce()
        fixture = conformance_fixture(effectors=(model,))
        base = C._solve(fixture.configuration)
        again = C._solve(fixture.configuration)
        NC_STALL[] = true
        check = C._concurrency_check(fixture.configuration, base, again, v"0.3.0"; timeout_s=1.0)
        @test !NC_STALL[]   # one concurrent run took the token and stalled
        if Threads.nthreads() >= 2
            @testset let status = check.status, detail = check.detail
                @test status === :fail
                @test occursin("the concurrent runs fell out of step", detail)
                @test !occursin("differ", detail)   # equal results do not rescue it
            end
        else
            @test check.status === :skip
        end
    end
    @testset "the state comparison matches mutable objects one to one" begin
        shared = [1.0, 2.0]
        aliased = NcBreaksAliasing(shared, shared)
        separate = NcBreaksAliasing([1.0, 2.0], [1.0, 2.0])
        @test C._same_state(aliased, deepcopy(aliased))
        @test C._same_state(separate, deepcopy(separate))
        @test !C._same_state(aliased, separate)
        @test !C._same_state(separate, aliased)
    end
    # The checks run copies; only the deepcopy_equivalent run mutates the model given.
    counter = NcCountsCalls(0)
    with_logger(NullLogger()) do
        check_force_torque_model(counter)
    end
    @test counter.calls > 0
end

end # module ConformanceNegativeControlsTests
