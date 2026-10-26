# GNC models: guidance, navigation and control effectors. The engine's pre-run
# check, run-level reproducibility (the state and the model's own state at the end
# of each run), the reset at each run's start, the control output hooks, and the
# state shared by the copies `bind_spacecraft` makes.

# The control output hooks the dynamics call on `(model, x::StateSample, p, t)`.
const _CONTROL_OUTPUT_HOOKS = (SimulationModel.control_wrench, SimulationModel.control_propellant_rate,
    SimulationModel.control_wheel_torque)

function _container_effectors(configuration, kind::Symbol)
    kind === :guidance_model && return configuration.guidance_model.guidance_effectors
    kind === :navigation_model && return configuration.navigation_model.navigation_effectors
    kind === :control_model && return configuration.control_model.control_effectors
    throw(ArgumentError("unknown model kind $(repr(kind))"))
end

function _container_rates(configuration, kind::Symbol)
    kind === :guidance_model && return configuration.guidance_model.guidance_rates
    kind === :navigation_model && return configuration.navigation_model.navigation_rates
    return configuration.control_model.control_rates
end

const _CONTAINER_NAMES = (guidance_model="guidance_model.guidance_effectors",
    navigation_model="navigation_model.navigation_effectors", control_model="control_model.control_effectors")

const _PER_SPACECRAFT_KEYS = (guidance_model=:guidance, navigation_model=:navigation, control_model=:control)

const _INITIALIZERS = (guidance_model="initialize_guidance!", navigation_model="initialize_navigation!",
    control_model="initialize_control!")

"""
    check_guidance_model(model; fixture=conformance_fixture(guidance=(model,)), spacecraft_index=1)
        -> ConformanceReport

Check a guidance model. `fixture` must hold `model` itself (found by identity) among
its configuration's guidance effectors. The checks, in report order:

- `hooks_present`: the engine's pre-run hook check, on the fixture's concrete argument
  types, reports nothing.
- `repeatable_run`, `deepcopy_equivalent`, `global_rng_independent`,
  `concurrent_equals_sequential`: as for force/torque models, comparing the saved
  states and the model's own state at the end of each run, bit for bit.
- `initialize_resets`: two runs with `isolate_state=false` on one configuration end
  with equal states and an equal model, so `initialize_guidance!` resets the model's
  run state.
- `bind_fresh_state`: for a model with a `bind_spacecraft` method, no mutable object
  reachable from the copies bound to spacecraft 1 and 2 changes during a run of the
  fixture with the model declared on its first two spacecraft instead of in the
  shared container (a fixture with one spacecraft gets a copy of it as the second,
  its mounted sensors included). Skipped without a method.
- `checkpoint_round_trip` (added in 0.4.0, advisory in 0.4): the fixture run with a
  checkpoint halfway, interrupted after it and resumed from the file, ends its second
  half with the states and the model's own state of an uninterrupted run with the same
  checkpoints, bit for bit: `checkpoint_state` and `restore_checkpoint_state!` carry
  the run state the model keeps from the checkpoint into the fresh copy the resumed
  run starts from. Skipped when the fixture refuses checkpointed runs at setup.

A pass certifies interface behaviour, not the guidance law.
"""
check_guidance_model(model; fixture::ConformanceFixture=conformance_fixture(guidance=(model,)),
    spacecraft_index::Integer=1) = _check_gnc_model(model, :guidance_model, fixture, Int(spacecraft_index))

"""
    check_navigation_model(model; sensors=(), fixture=conformance_fixture(navigation=(model,), sensors=sensors),
        spacecraft_index=1) -> ConformanceReport

Check a navigation model, with `sensors` mounted for it to read. The checks are those
of [`check_guidance_model`](@ref), with `initialize_resets` checking that
`initialize_navigation!` resets the model's run state.
"""
check_navigation_model(model; sensors::Tuple=(),
    fixture::ConformanceFixture=conformance_fixture(navigation=(model,), sensors=sensors), spacecraft_index::Integer=1) =
    _check_gnc_model(model, :navigation_model, fixture, Int(spacecraft_index))

"""
    check_control_model(model; fixture=conformance_fixture(control=(model,)), spacecraft_index=1)
        -> ConformanceReport

Check a control model. The checks, in report order: `hooks_present`;
`supported_sample_signature` (no method of `control_wrench`,
`control_propellant_rate` or `control_wheel_torque` that applies to the model types its
sample argument narrower than `StateSample`); the four run-level checks of
[`check_guidance_model`](@ref); `initialize_resets` (`initialize_control!` resets the
model's run state, as for guidance); `output_hook_pure` (each output hook the model has,
called twice on the same arguments, returns equal values, leaves the model unchanged
and logs nothing at `Info` level or above); `positive_consumption`
(`control_propellant_rate` is finite and not negative); `bind_fresh_state`;
`checkpoint_round_trip`, as for guidance; `outputs_in_scope` (added in 0.7.0,
advisory in 0.7): after the model's ticks, its output hooks and levels
(`actuator_levels`, `control_thruster_levels`) return zeros or `nothing` for a
fixture spacecraft outside the `spacecraft_scope` the run took, once the model's
run-start initializers ran. From 0.7.0 the engine calls them only inside the scope,
but direct callers may ask about any spacecraft. A fixture with one spacecraft is
checked with a copy of it as the second; skipped when the scope names every
spacecraft of the fixture.

The output hooks are called on the states, parameters and copies of the model the
engine holds at the start, quarters and end of a fixture run, on spacecraft
`spacecraft_index`. A model with no output hook skips those checks.
"""
check_control_model(model; fixture::ConformanceFixture=conformance_fixture(control=(model,)),
    spacecraft_index::Integer=1) = _check_gnc_model(model, :control_model, fixture, Int(spacecraft_index))

"""
    test_guidance_model(model; kwargs...)

Run [`check_guidance_model`](@ref) and record its report with [`test_report`](@ref).
Needs the `Test` standard library: `using Test` loads the method.
"""
function test_guidance_model end

"""
    test_navigation_model(model; kwargs...)

Run [`check_navigation_model`](@ref) and record its report with [`test_report`](@ref).
Needs the `Test` standard library: `using Test` loads the method.
"""
function test_navigation_model end

"""
    test_control_model(model; kwargs...)

Run [`check_control_model`](@ref) and record its report with [`test_report`](@ref).
Needs the `Test` standard library: `using Test` loads the method.
"""
function test_control_model end

function _check_gnc_model(model, kind::Symbol, fixture::ConformanceFixture, spacecraft_index::Int)
    configuration = fixture.configuration
    index = findfirst(effector -> effector === model, _container_effectors(configuration, kind))
    index === nothing && throw(ArgumentError(
        "check_$(kind): the fixture's $(_CONTAINER_NAMES[kind]) does not contain this model (it is found by " *
        "identity); use conformance_fixture($(_PER_SPACECRAFT_KEYS[kind])=(model,)) or a configuration that holds it."))
    since = FIRST_RELEASE
    hooks = _hooks_present_check(fixture, since)
    run(cfg=configuration; kwargs...) = _model_solve(cfg, kind, index; kwargs...)
    base = run()
    again = run()
    seeded = map(seed -> _with_global_seed(() -> run(), seed), (1, 2))
    concurrency = _concurrency_check(factories -> run(; factories), base, again, since;
        difference=_model_run_difference)
    if kind === :control_model
        samples = _capture_control_samples(configuration, index, spacecraft_index)
        control_checks = [
            _output_hooks_pure_check(model, samples, since),
            _positive_consumption_check(model, samples, since),
        ]
    end
    reused = deepcopy(configuration)
    reset_gap = _model_run_difference(run(reused; isolate_state=false), run(reused; isolate_state=false))
    resets = _check(:initialize_resets, since, reset_gap === nothing,
        "two runs with isolate_state=false on one configuration differ: $(reset_gap); reset the model's run " *
        "state in $(_INITIALIZERS[kind])")
    bound = _bind_fresh_state_check(model, kind, since; configuration, index)
    round_trip = _gnc_checkpoint_round_trip(model, kind, configuration)
    kind === :control_model && (in_scope = _outputs_in_scope_check(model, configuration, index))
    # Last: this run uses the caller's objects.
    original = run(; isolate_state=false)

    repeat_gap = _model_run_difference(base, again)
    copy_gap = _model_run_difference(original, base)
    seed_gap = _model_run_difference(seeded...)
    runs = [
        _check(:repeatable_run, since, repeat_gap === nothing,
            "two runs of the fixture with isolate_state=true differ: $(repeat_gap)"),
        _check(:deepcopy_equivalent, since, copy_gap === nothing,
            "a run of the model itself (isolate_state=false) and a run of its deep copy (isolate_state=true) " *
            "differ: $(copy_gap)"),
        _check(:global_rng_independent, since, seed_gap === nothing,
            "runs with the global random number generator seeded 1 and 2 differ: $(seed_gap); draw random " *
            "numbers only from a stream the engine supplies, never from rand() or Random.default_rng()"),
        concurrency,
    ]
    checks = if kind === :control_model
        [hooks, _sample_signature_check(model, since; hooks=_CONTROL_OUTPUT_HOOKS), runs..., resets, control_checks...,
            bound, round_trip, in_scope]
    else
        [hooks, runs..., resets, bound, round_trip]
    end
    return ConformanceReport(string(nameof(typeof(model))), kind, checks)
end

# The checks that need no run: the pre-run check and, for a control model, its
# sample signature. Core's tests use them on built-in models that cannot run in the
# fixture.
function _static_gnc_report(model, kind::Symbol;
        fixture::ConformanceFixture=conformance_fixture(; NamedTuple{(_PER_SPACECRAFT_KEYS[kind],)}(((model,),))...))
    since = FIRST_RELEASE
    checks = [_hooks_present_check(fixture, since)]
    kind === :control_model && push!(checks, _sample_signature_check(model, since; hooks=_CONTROL_OUTPUT_HOOKS))
    return ConformanceReport(string(nameof(typeof(model))), kind, checks)
end

# A run of the configuration that also returns the model as it ended the run.
function _model_solve(configuration, kind::Symbol, index::Int; isolate_state::Bool=true, factories::Tuple=())
    try
        result = SimulationEngine.run_simulation(configuration; isolate_state, return_solution=true,
            return_results=true, extra_callback_factories=factories)
        ran = _container_effectors(result.configuration, kind)[index]
        return _RunOutcome((t=result.solution.t, u=result.solution.u, model=deepcopy(ran)), nothing)
    catch err
        return _RunOutcome(nothing, err)
    end
end

function _model_run_difference(a::_RunOutcome, b::_RunOutcome)
    gap = _run_difference(a, b)
    gap === nothing || return gap
    _same_state(a.solution.model, b.solution.model) ||
        return "the model's own state at the end of the runs differs"
    return nothing
end

# The model, the spacecraft's state sample, the run's parameters and the time at the
# start, quarters and end of a fixture run: what the dynamics pass to the output hooks.
# `spacecraft_index` is a run position, or a function of the run's parameters that
# gives one (or `nothing`, for no samples), asked at the first sample.
function _capture_control_samples(configuration, index::Int, spacecraft_index)
    targets = [fraction * configuration.mission_configuration.mission_time for fraction in _SAMPLE_FRACTIONS]
    arguments = Any[]
    failure = Ref{Any}(nothing)
    next = Ref(1)
    resolved = Ref{Any}(spacecraft_index isa Integer ? Int(spacecraft_index) : missing)
    function capture!(integrator)
        t = Float64(integrator.t)
        if failure[] === nothing && next[] <= length(targets) && t >= targets[next[]]
            while next[] <= length(targets) && t >= targets[next[]]
                next[] += 1
            end
            try
                p = integrator.p
                args = p.args
                resolved[] === missing && (resolved[] = spacecraft_index(p))
                k = resolved[]
                if k !== nothing
                    spacecraft = args.dynamics_model.spacecraft[k]
                    x = SimulationEngine.build_state_sample(integrator.u.sc[k], spacecraft,
                        args.mission_configuration.orientation_sim, k)
                    push!(arguments, (deepcopy(args.control_model.control_effectors[index]), x, p, t))
                end
            catch err
                failure[] = err
            end
        end
        # Sampling reads the run and changes no force input.
        SciMLBase.derivative_discontinuity!(integrator, false)
        return nothing
    end
    factory = configuration -> SciMLBase.DiscreteCallback((u, t, integrator) -> true, capture!;
        initialize=(callback, u, t, integrator) -> capture!(integrator), save_positions=(false, false))
    run = _solve(configuration; extra_callback_factories=(factory,))
    failure[] === nothing && _failed(run) && isempty(arguments) && (failure[] = run.error)
    return _Samples(arguments, failure[])
end

# Whether the model has an output hook of its own: the untyped fallbacks of
# control_propellant_rate and control_wheel_torque, which serve every model, do not
# count.
function _has_output_hook(f, model, p)
    tt = Tuple{typeof(model), SimulationModel.StateSample, typeof(p), Float64}
    hasmethod(f, tt) || return false
    fallback = Tuple{Any, Any, Any, Any}
    return !hasmethod(f, fallback) || which(f, tt) !== which(f, fallback)
end

function _output_hooks_pure_check(model, samples::_Samples, since::VersionNumber)
    samples.error === nothing || return _check(:output_hook_pure, since, false, _sampling_failure(samples))
    isempty(samples.arguments) && return _skip(:output_hook_pure, since, "no samples to call the output hooks with")
    problems = String[]
    called = false
    log = _MessageLog(String[])
    Logging.with_logger(log) do
        for (subject, x, p, t) in samples.arguments, hook in _CONTROL_OUTPUT_HOOKS
            _has_output_hook(hook, subject, p) || continue
            called = true
            snapshot = deepcopy(subject)
            outputs = try
                (hook(subject, x, p, t), hook(subject, x, p, t))
            catch err
                push!(problems, "$(nameof(hook)) threw at t = $(t) s: " * _summary(err))
                continue
            end
            isequal(outputs...) || push!(problems,
                "two calls of $(nameof(hook)) with the same arguments at t = $(t) s returned $(outputs[1]) and then $(outputs[2])")
            _same_state(subject, snapshot) ||
                push!(problems, "calls of $(nameof(hook)) at t = $(t) s changed the model")
        end
    end
    called || return _skip(:output_hook_pure, since, "the model has none of the control output hooks " *
        "(control_wrench, control_propellant_rate, control_wheel_torque)")
    isempty(log.messages) ||
        push!(problems, "the calls logged $(length(log.messages)) message(s), the first: $(first(log.messages))")
    return _check(:output_hook_pure, since, isempty(problems), join(unique(problems), "\n"))
end

function _positive_consumption_check(model, samples::_Samples, since::VersionNumber)
    samples.error === nothing || return _check(:positive_consumption, since, false, _sampling_failure(samples))
    hook = SimulationModel.control_propellant_rate
    relevant = [sample for sample in samples.arguments if _has_output_hook(hook, sample[1], sample[3])]
    isempty(relevant) && return _skip(:positive_consumption, since,
        "the model has no control_propellant_rate method, so the dynamics take its consumption as 0")
    for (subject, x, p, t) in relevant
        rate = try
            hook(deepcopy(subject), x, p, t)
        catch err
            return _check(:positive_consumption, since, false, "control_propellant_rate threw at t = $(t) s: " * _summary(err))
        end
        rate isa Real && isfinite(rate) && rate >= 0 && continue
        return _check(:positive_consumption, since, false, "control_propellant_rate returned $(repr(rate)) at " *
            "t = $(t) s; it is the propellant used, kg/s, finite and not negative")
    end
    return _check(:positive_consumption, since, true)
end

# `outputs_in_scope`: a fixture spacecraft outside the scope the run took for the
# model, once its run-start initializers ran, sampled as the output hooks see it,
# after the model's ticks; each output hook and level the model has returns zeros
# or `nothing` there. A fixture with one spacecraft gets a copy of it as the
# second, so a model bound to the first has a spacecraft outside its scope.
const _OUTPUTS_IN_SCOPE_SINCE = v"0.7.0"

function _outputs_in_scope_check(model, configuration, index::Int)
    since = _OUTPUTS_IN_SCOPE_SINCE
    added = length(configuration.dynamics_model.spacecraft) == 1
    added && GNCRuntime.spacecraft_scope(model) isa GNCRuntime.AllSpacecraft && return _skip(:outputs_in_scope, since,
        "the model's spacecraft_scope names every spacecraft and the fixture has one; a fixture with a spacecraft " *
        "outside its scope checks this")
    fixture = if added
        try
            _with_second_spacecraft(configuration)
        catch err
            return _check(:outputs_in_scope, since, false,
                "the two-spacecraft fixture could not be built: " * _summary(err))
        end
    else
        configuration
    end
    # With the spacecraft this check adds, a run or a hook that fails for it says more
    # about the fixture than about the model, as for a model sized for one spacecraft.
    unchecked(why) = _skip(:outputs_in_scope, since, "the check runs the fixture with a copy of its spacecraft as " *
        "the second, and $(why); a fixture with a spacecraft outside the model's scope checks this")
    samples = _capture_control_samples(fixture, index, p -> _outside_run_scope(p, index))
    if samples.error !== nothing
        added && return unchecked("that run failed: " * _summary(samples.error))
        return _check(:outputs_in_scope, since, false, _sampling_failure(samples))
    end
    isempty(samples.arguments) && return _skip(:outputs_in_scope, since,
        "the scope the run took for the model names every spacecraft of the fixture; a fixture with a spacecraft " *
        "outside its scope checks this")
    for (subject, x, p, t) in samples.arguments
        problem = _output_outside_scope(subject, x, p, t, x.spacecraft_index)
        problem === nothing && continue
        added && problem.threw && return unchecked("an output hook threw for the copy: " * problem.message)
        return _check(:outputs_in_scope, since, false, problem.message)
    end
    return _check(:outputs_in_scope, since, true)
end

# The first spacecraft of the run outside the positions the run took for control
# model `index`, or `nothing`.
function _outside_run_scope(p, index::Int)
    scopes = p.shared_buffers.gnc_scopes[]
    scopes === nothing && return nothing
    positions = scopes.control[index].positions
    return findfirst(i -> !(i in positions), eachindex(p.args.dynamics_model.spacecraft))
end

# A deep copy of the configuration, with a copy of its spacecraft as the second (the
# next free id) when it has one.
function _with_second_spacecraft(configuration)
    variant = deepcopy(configuration)
    spacecraft = collect(SimulationModel.SpacecraftModel, variant.dynamics_model.spacecraft)
    length(spacecraft) == 1 || return variant
    second = deepcopy(only(spacecraft))
    second.id = only(spacecraft).id + 1
    push!(spacecraft, second)
    return SimulationModel.SimConfig._with_configuration(variant;
        dynamics_model=SimulationModel.DynamicsModel(spacecraft, variant.dynamics_model.dynamic_effectors))
end

# What an output hook or level of `subject` returns for spacecraft `outside` that is
# not zero or `nothing`, as `(message, threw)`, or `nothing` when every one is.
function _output_outside_scope(subject, x, p, t, outside::Int)
    where_ = "for spacecraft $(outside), outside its spacecraft_scope, at t = $(t) s"
    returned(message) = (message=message, threw=false)
    try
        if _has_output_hook(SimulationModel.control_wrench, subject, p)
            force, torque = SimulationModel.control_wrench(deepcopy(subject), x, p, t)
            all(iszero, force) && all(iszero, torque) || return returned(
                "control_wrench returned force $(force) and torque $(torque) $(where_); return zeros there")
        end
        if _has_output_hook(SimulationModel.control_propellant_rate, subject, p)
            rate = SimulationModel.control_propellant_rate(deepcopy(subject), x, p, t)
            iszero(rate) || return returned("control_propellant_rate returned $(repr(rate)) $(where_); return 0.0 there")
        end
        if _has_output_hook(SimulationModel.control_wheel_torque, subject, p)
            wheel = SimulationModel.control_wheel_torque(deepcopy(subject), x, p, t)
            wheel === nothing || all(iszero, wheel) ||
                return returned("control_wheel_torque returned $(repr(wheel)) $(where_); return nothing there")
        end
        levels = SimulationModel.actuator_levels(subject, outside)
        levels === nothing || return returned("actuator_levels returned $(repr(levels)) $(where_); return nothing there")
        thrusters = SimulationModel.control_thruster_levels(subject, outside)
        thrusters === nothing ||
            return returned("control_thruster_levels returned $(repr(thrusters)) $(where_); return nothing there")
    catch err
        return (message="an output hook threw $(where_): " * _summary(err), threw=true)
    end
    return nothing
end

# The check's own configuration with the model, at `index` of its container, declared
# on the first two spacecraft instead: every other model, sensor, environment and
# setting is kept, so a model whose bound copies need them runs as in the main runs.
# A configuration with one spacecraft gets a copy of it as the second, with the next
# free id.
function _bind_configuration(configuration, kind::Symbol, index::Int, subject)
    variant = deepcopy(configuration)
    effectors = collect(Any, _container_effectors(variant, kind))
    rates = collect(Float64, _container_rates(variant, kind))
    rate = rates[index]
    deleteat!(effectors, index)
    deleteat!(rates, index)
    spacecraft = collect(SimulationModel.SpacecraftModel, variant.dynamics_model.spacecraft)
    if length(spacecraft) == 1
        second = deepcopy(only(spacecraft))
        second.id = only(spacecraft).id + 1
        push!(spacecraft, second)
    end
    for sc in spacecraft[1:2]
        if kind === :guidance_model
            sc.guidance = SimulationModel.GuidanceModel((sc.guidance.guidance_effectors..., subject),
                [sc.guidance.guidance_rates; rate])
        elseif kind === :navigation_model
            sc.navigation = SimulationModel.NavigationModel((sc.navigation.navigation_effectors..., subject),
                [sc.navigation.navigation_rates; rate])
        else
            sc.control = SimulationModel.ControlModel((sc.control.control_effectors..., subject),
                [sc.control.control_rates; rate])
        end
    end
    container = kind === :guidance_model ? (guidance_model=SimulationModel.GuidanceModel(Tuple(effectors), rates),) :
        kind === :navigation_model ? (navigation_model=SimulationModel.NavigationModel(Tuple(effectors), rates),) :
        (control_model=SimulationModel.ControlModel(Tuple(effectors), rates),)
    return SimulationModel.SimConfig._with_configuration(variant; container...,
        dynamics_model=SimulationModel.DynamicsModel(spacecraft, variant.dynamics_model.dynamic_effectors))
end

function _bind_fresh_state_check(model, kind::Symbol, since::VersionNumber; configuration, index::Int)
    bind = SimulationLifecycle.bind_spacecraft
    T = typeof(model)
    hasmethod(bind, Tuple{T, Int}) && which(bind, Tuple{T, Int}) !== which(bind, Tuple{Any, Int}) ||
        return _skip(:bind_fresh_state, since, "the model has no bind_spacecraft method, so it cannot be declared " *
            "on spacecraft")
    subject = deepcopy(model)
    shared = try
        _shared_mutables(bind(subject, 1), bind(subject, 2))
    catch err
        return _check(:bind_fresh_state, since, false, "bind_spacecraft threw: " * _summary(err))
    end
    snapshots = IdDict{Any, Any}(object => deepcopy(object) for object in shared)
    two = try
        _bind_configuration(configuration, kind, index, subject)
    catch err
        return _check(:bind_fresh_state, since, false, "the two-spacecraft fixture could not be built: " * _summary(err))
    end
    # The engine binds the declared model to each spacecraft; with isolate_state=false
    # the bound copies are the ones the run uses.
    run = _solve(two; isolate_state=false)
    _failed(run) && return _check(:bind_fresh_state, since, false,
        "the two-spacecraft run failed: " * _summary(run.error))
    changed = [object for object in shared if !_same_state(object, snapshots[object])]
    return _check(:bind_fresh_state, since, isempty(changed),
        "$(length(changed)) mutable object(s) reachable from both the copy bound to spacecraft 1 and the copy " *
        "bound to spacecraft 2 changed during the run (the first is a $(isempty(changed) ? "" : typeof(first(changed))));" *
        " give each bound copy its own run state in bind_spacecraft")
end
