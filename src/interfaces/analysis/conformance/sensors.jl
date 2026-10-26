# Sensors: the run state the engine fills, the requirements it supports, the
# measurement's shape, run-level reproducibility of the measurements, the sensor's
# own random stream and the reset of its own state at the start of a run.

"""
    check_sensor(sensor; fixture=conformance_fixture(sensors=(sensor,)), seeds=(7, 8))
        -> ConformanceReport

Check a sensor kind against the engine's sensor contract. `fixture` must have
`sensor` itself (found by identity) mounted on one of its spacecraft. The checks, in
report order:

- `run_state_fields`: the sensor is a mutable struct with the fields the engine fills
  (`period_s`, `seed`, `rng`, `latest`, `latest_time_s`), a finite positive period, a
  positive `sensor_dimension` and a `latest` with that many components.
- `requirements_supported`: the engine's `add_sensor!` accepts a copy of the sensor; it
  refuses the environment requests (`sensor_requirements`) the engine does not sample
  for sensors, and any record it cannot mount.
- `measurement_shape`: `sensor_measure` returns `sensor_dimension` finite components
  (a number or a vector for one component, a vector otherwise), in the fixture run and
  when called on the truths the engine builds at the start, quarters and end of it.
- `repeatable_run`, `deepcopy_equivalent`, `global_rng_independent` and
  `concurrent_equals_sequential`: as for force/torque models, comparing the sensor's
  measurements and their times, bit for bit, instead of the state.
- `seed_repeatable`: two runs with the first of `seeds` give the same measurements, and
  a sensor whose measurement depends on its random stream (two streams give different
  measurements of the same truth) measures differently in a run with the second seed.
- `custom_state_reset`: two runs with `isolate_state=false` on one configuration give
  the same measurements, so any state the sensor keeps beyond the engine's fields is
  reset by `initialize_sensor!`.
- `checkpoint_round_trip` (added in 0.4.0, advisory in 0.4): the fixture run with a
  checkpoint halfway, interrupted after it and resumed from the file, gives the
  measurements of an uninterrupted run with the same checkpoints from halfway on, so
  any state the sensor keeps beyond the engine's fields is in `checkpoint_state`.
  Skipped when the fixture refuses checkpointed runs at setup.

When the fixture cannot be built (`add_sensor!` refuses the sensor), the checks that
need a run fail with its error. Every run deep-copies the fixture except the
`deepcopy_equivalent` run, which runs `sensor` itself after all the others. A pass
certifies interface behaviour, not the sensor's model of the instrument.
"""
function check_sensor(sensor; fixture::Union{Nothing, ConformanceFixture}=nothing, seeds::Tuple{Integer, Integer}=(7, 8))
    since = FIRST_RELEASE
    fields = _sensor_fields_check(sensor, since)
    requirements = _sensor_requirements_check(sensor, since)
    built = fixture === nothing ? _built_fixture(() -> conformance_fixture(sensors=(sensor,))) : fixture
    if built isa Exception
        why = "the fixture could not be built: " * _summary(built)
        rest = [_check(id, id === :checkpoint_round_trip ? CHECKPOINT_ROUND_TRIP_SINCE : since, false, why)
                for id in CHECK_IDS.sensor[3:end]]
        return ConformanceReport(string(nameof(typeof(sensor))), :sensor, [fields, requirements, rest...])
    end
    configuration = built.configuration
    position = _sensor_position(configuration, sensor)
    position === nothing && throw(ArgumentError(
        "check_sensor: no spacecraft of the fixture has this sensor mounted (it is found by identity); use " *
        "conformance_fixture(sensors=(sensor,)) or a configuration on which it is mounted with add_sensor!."))
    run(cfg=configuration; kwargs...) = _sensor_run(cfg, position; kwargs...)
    base = run()
    again = run()
    seeded = map(seed -> _with_global_seed(() -> run().outcome, seed), (1, 2))
    concurrency = _concurrency_check(factories -> run(; factories).outcome, base.outcome, again.outcome, since;
        noun="measurement")
    shape = _measurement_shape_check(sensor, base, first(seeds), since)
    seeding = _seed_check(sensor, configuration, position, base, seeds, since)
    reused = deepcopy(configuration)
    first_reuse = run(reused; isolate_state=false)
    second_reuse = run(reused; isolate_state=false)
    half = configuration.mission_configuration.mission_time / 2
    round_trip = _checkpoint_round_trip_check((cfg, factories) -> run(cfg; factories).outcome,
        (a, b) -> _run_difference(_measurements_from(a, half), _measurements_from(b, half); noun="measurement"),
        sensor, configuration)
    # Last: this run uses the caller's objects.
    original = run(; isolate_state=false)

    measure(a, b) = _run_difference(a.outcome, b.outcome; noun="measurement")
    repeat_gap = measure(base, again)
    copy_gap = measure(original, base)
    seed_gap = _run_difference(seeded...; noun="measurement")
    reset_gap = measure(first_reuse, second_reuse)
    checks = [
        fields,
        requirements,
        shape,
        _check(:repeatable_run, since, repeat_gap === nothing,
            "two runs of the fixture with isolate_state=true differ: $(repeat_gap)"),
        _check(:deepcopy_equivalent, since, copy_gap === nothing,
            "a run of the sensor itself (isolate_state=false) and a run of its deep copy (isolate_state=true) " *
            "differ: $(copy_gap)"),
        _check(:global_rng_independent, since, seed_gap === nothing,
            "runs with the global random number generator seeded 1 and 2 differ: $(seed_gap); draw noise only " *
            "from the `rng` argument of sensor_measure, never from rand() or Random.default_rng()"),
        concurrency,
        seeding,
        _check(:custom_state_reset, since, reset_gap === nothing,
            "two runs with isolate_state=false on one configuration differ: $(reset_gap); reset the state the " *
            "sensor keeps beyond the engine's fields in initialize_sensor!"),
        round_trip,
    ]
    return ConformanceReport(string(nameof(typeof(sensor))), :sensor, checks)
end

"""
    test_sensor(sensor; kwargs...)

Run [`check_sensor`](@ref) and record its report with [`test_report`](@ref). Needs
the `Test` standard library: `using Test` loads the method.
"""
function test_sensor end

# A run's measurements taken at `t_s` or later.
function _measurements_from(outcome::_RunOutcome, t_s::Float64)
    _failed(outcome) && return outcome
    kept = findall(>=(t_s), outcome.solution.t)
    return _RunOutcome((t=outcome.solution.t[kept], u=outcome.solution.u[kept]), nothing)
end

_built_fixture(build) = try
    build()
catch err
    err
end

# The spacecraft and mounting position of `sensor` in `configuration`.
function _sensor_position(configuration, sensor)
    for (i, spacecraft) in enumerate(configuration.dynamics_model.spacecraft)
        k = findfirst(entry -> entry.sensor === sensor, SimulationModel.mounted_sensors(spacecraft))
        k === nothing || return (i, k)
    end
    return nothing
end

_mounted_entry(configuration, (i, k)) = SimulationModel.mounted_sensors(configuration.dynamics_model.spacecraft[i])[k]

function _sensor_fields_check(sensor, since::VersionNumber)
    T = typeof(sensor)
    problems = String[]
    ismutabletype(T) || push!(problems, "$(nameof(T)) is not a mutable struct; the engine fills its run state")
    absent = [f for f in SimulationModel.SensorMounting._RUN_STATE_FIELDS if !hasfield(T, f)]
    isempty(absent) || push!(problems, "$(nameof(T)) lacks the field(s) $(join(absent, ", "))")
    hasfield(T, :period_s) && try
        SimulationModel.SensorMounting._check_period(sensor)
    catch err
        push!(problems, _summary(err))
    end
    hasfield(T, :seed) && !(getfield(sensor, :seed) isa Integer) &&
        push!(problems, "its seed is $(repr(getfield(sensor, :seed))), not an integer")
    dimension = try
        SimulationModel.sensor_dimension(sensor)
    catch err
        push!(problems, "sensor_dimension threw: " * _summary(err))
        nothing
    end
    if dimension !== nothing
        dimension isa Integer && dimension >= 1 || push!(problems, "sensor_dimension is $(repr(dimension)), not a positive integer")
        if hasfield(T, :latest) && dimension isa Integer
            latest = getfield(sensor, :latest)
            length(latest) == dimension || push!(problems,
                "its latest measurement has $(length(latest)) component(s), not sensor_dimension = $(dimension)")
        end
    end
    return _check(:run_state_fields, since, isempty(problems), join(problems, "\n"))
end

# The engine decides which environment requests it samples for sensors: its
# add_sensor! refuses the others. A copy of the sensor is mounted on a spacecraft of
# a scratch configuration, so the rule is the engine's own.
function _sensor_requirements_check(sensor, since::VersionNumber)
    spacecraft = first(_fixture_configuration().dynamics_model.spacecraft)
    try
        SimulationModel.add_sensor!(spacecraft, spacecraft.root, deepcopy(sensor))
    catch err
        return _check(:requirements_supported, since, false, "add_sensor! refuses the sensor: " * _summary(err))
    end
    return _check(:requirements_supported, since, true)
end

# A sensor run: its measurement history (`outcome.solution` holds `t` and `u`, one
# measurement per sample) and the engine's truths and mounting contexts at the
# capture times.
struct _SensorRun
    outcome::_RunOutcome
    samples::Vector{Any}
    failure::Any
end

function _sensor_run(configuration, position::Tuple{Int, Int}; isolate_state::Bool=true, factories::Tuple=())
    times = Float64[]
    values = Vector{Float64}[]
    samples = Any[]
    failure = Ref{Any}(nothing)
    targets = [fraction * configuration.mission_configuration.mission_time for fraction in _SAMPLE_FRACTIONS]
    next = Ref(1)
    function record!(integrator)
        try
            args = integrator.p.args
            entry = _mounted_entry(args, position)
            s = entry.sensor
            t = Float64(s.latest_time_s)
            if !isnan(t) && (isempty(times) || t != times[end])
                push!(times, t)
                push!(values, Float64[x for x in s.latest])
            end
            now = Float64(integrator.t)
            if failure[] === nothing && next[] <= length(targets) && now >= targets[next[]]
                while next[] <= length(targets) && now >= targets[next[]]
                    next[] += 1
                end
                truth = SimulationModel.SimulationCallbacks._sensor_truth(integrator, first(position), s)
                push!(samples, (truth, SimulationModel.SimulationCallbacks.actuator_context(entry.link)))
            end
        catch err
            failure[] === nothing && (failure[] = err)
        end
        # Recording reads the sensor and changes no force input.
        SciMLBase.derivative_discontinuity!(integrator, false)
        return nothing
    end
    recorder = configuration -> SciMLBase.DiscreteCallback((u, t, integrator) -> true, record!;
        initialize=(callback, u, t, integrator) -> record!(integrator), save_positions=(false, false))
    run = _solve(configuration; isolate_state, extra_callback_factories=(factories..., recorder))
    outcome = _failed(run) ? run : _RunOutcome((t=times, u=values), nothing)
    return _SensorRun(outcome, samples, failure[])
end

_has_shape(m, dimension) = (m isa Real && dimension == 1) ||
    (m isa AbstractVector{<:Real} && length(m) == dimension)

function _measurement_shape_check(sensor, run::_SensorRun, seed::Integer, since::VersionNumber)
    dimension = try
        SimulationModel.sensor_dimension(sensor)
    catch err
        return _check(:measurement_shape, since, false, "sensor_dimension threw: " * _summary(err))
    end
    problems = String[]
    if _failed(run.outcome)
        push!(problems, "the fixture run failed: " * _summary(run.outcome.error))
    else
        history = run.outcome.solution
        isempty(history.u) && push!(problems, "the fixture run recorded no measurement")
        k = findfirst(m -> length(m) != dimension || !all(isfinite, m), history.u)
        k === nothing || push!(problems, "the measurement at t = $(history.t[k]) s is $(history.u[k]), not " *
            "$(dimension) finite component(s)")
    end
    run.failure === nothing || push!(problems, "the fixture run could not be sampled: " * _summary(run.failure))
    for (truth, ctx) in run.samples
        m = try
            only(_fresh_measurements(sensor, truth, ctx, seed, 1))
        catch err
            push!(problems, "a copy of the sensor, reset and initialized as the engine does, could not measure at " *
                "t = $(truth.t_s) s: " * _summary(err))
            break
        end
        _has_shape(m, dimension) && all(isfinite, m) && continue
        push!(problems, "sensor_measure returned $(repr(m)) at t = $(truth.t_s) s; it must return " *
            "sensor_dimension = $(dimension) finite component(s)")
        break
    end
    return _check(:measurement_shape, since, isempty(problems), join(problems, "\n"))
end

# A copy of `sensor` made ready as the engine makes a sensor ready for a run: mounted
# on a scratch spacecraft and reset by the engine's `reset_sensors!` (its stream from
# `seed`, its latest measurement unset, a gyro's bias at its initial value), then
# `initialize_sensor!` with the first sample's truth and context. `count` successive
# measurements of `truth` follow, each drawing from `stream` (default: the engine's
# stream for `seed`). Initialization gets the engine's stream, so a sensor that
# reseeds the stream it is given there cannot hide its draws from `stream`.
function _fresh_measurements(sensor, truth, ctx, seed::Integer, count::Integer; stream=nothing)
    spacecraft = first(_fixture_configuration().dynamics_model.spacecraft)
    s = deepcopy(sensor)
    hasfield(typeof(s), :seed) && (s.seed = seed)
    SimulationModel.add_sensor!(spacecraft, spacecraft.root, s)
    SimulationModel.reset_sensors!(spacecraft)
    SimulationModel.initialize_sensor!(s, truth, ctx, s.rng)
    rng = stream === nothing ? s.rng : stream
    return [SimulationModel.sensor_measure(s, truth, ctx, rng) for _ in 1:count]
end

# Whether the sensor's measurements depend on its random stream: a short sequence of
# measurements of one truth, from engine-ready copies of the sensor, with streams of
# two seeds. A sequence, not one measurement: a gyro whose bias walks returns the
# current bias before drawing its next step.
const _NOISE_PROBE_LENGTH = 4
function _draws_noise(sensor, run::_SensorRun, seeds)
    isempty(run.samples) && return false
    truth, ctx = first(run.samples)
    try
        a = _fresh_measurements(sensor, truth, ctx, seeds[1], _NOISE_PROBE_LENGTH; stream=Random.Xoshiro(seeds[1]))
        b = _fresh_measurements(sensor, truth, ctx, seeds[1], _NOISE_PROBE_LENGTH; stream=Random.Xoshiro(seeds[2]))
        return !isequal(a, b)
    catch
        return false
    end
end

function _seed_check(sensor, configuration, position, base::_SensorRun, seeds, since::VersionNumber)
    hasfield(typeof(sensor), :seed) || return _skip(:seed_repeatable, since, "the sensor has no seed field")
    function seeded(seed)
        variant = deepcopy(configuration)
        _mounted_entry(variant, position).sensor.seed = seed
        return _sensor_run(variant, position).outcome
    end
    first_run, repeat_run, other_run = seeded(seeds[1]), seeded(seeds[1]), seeded(seeds[2])
    problems = String[]
    gap = _run_difference(first_run, repeat_run; noun="measurement")
    gap === nothing || push!(problems, "two runs with seed $(seeds[1]) differ: $(gap)")
    if _draws_noise(sensor, base, seeds) && _run_difference(first_run, other_run; noun="measurement") === nothing
        push!(problems, "the sensor's measurement depends on its random stream (streams of seeds $(seeds[1]) and " *
            "$(seeds[2]) give different measurements of one truth), yet runs with those seeds give the same " *
            "measurements: the run does not draw from the stream its seed selects")
    end
    return _check(:seed_repeatable, since, isempty(problems), join(problems, "\n"))
end
