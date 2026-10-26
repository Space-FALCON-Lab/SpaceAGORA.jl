# The checkpoint round trip: a model's run state saved in a checkpoint file and put
# back into the fresh copy a resumed run starts from gives the rest of the run that
# an uninterrupted run gives. Added in 0.4.0, so advisory in 0.4 (K-4).

const CHECKPOINT_ROUND_TRIP_SINCE = v"0.4.0"

# Thrown by the interrupting callback, after the checkpoint written halfway.
struct _Interrupted <: Exception end
Base.showerror(io::IO, ::_Interrupted) = print(io, "the run was interrupted after its checkpoint")

_unwrapped(err) = err isa TaskFailedException ? _unwrapped(err.task.exception) : err

# `configuration` checkpointed halfway through its mission into `directory`, or
# resumed from there.
function _checkpointed(configuration, directory::AbstractString; resume::Bool=false)
    settings = SimulationModel.SimConfig._with_fields(configuration.simulation_settings;
        checkpoint_enabled=true, checkpoint_interval_s=configuration.mission_configuration.mission_time / 2,
        checkpoint_directory=String(directory), resume_from_checkpoint=resume)
    return SimulationModel.SimConfig._with_configuration(configuration; simulation_settings=settings)
end

# Warnings of these runs repeat the fixture's other runs' (and the engine's note that
# a checkpointed run returns its last segment's solution).
_without_warnings(f) = Logging.with_logger(f, Logging.ConsoleLogger(stderr, Logging.Error))

# The round trip of `subject` (a model or a sensor) in `configuration`, checkpointed
# halfway: an uninterrupted run, and a run interrupted after its checkpoint and then
# resumed from it. `run(configuration, factories)` runs a configuration with extra
# callback factories and returns its outcome; `difference(a, b)` compares the resumed
# and the uninterrupted runs' outcomes over the second half, or returns `nothing`.
# The interrupting callback also marks when the solve starts (the solver initializes
# it then), so only an `ArgumentError` before that counts as a refusal at setup.
function _checkpoint_round_trip_check(run::Function, difference::Function, subject, configuration,
        since::VersionNumber=CHECKPOINT_ROUND_TRIP_SINCE)
    id = :checkpoint_round_trip
    half = configuration.mission_configuration.mission_time / 2
    armed, solving = Ref(false), Ref(false)
    started! = function (callback, u, t, integrator)
        solving[] = true
        SciMLBase.derivative_discontinuity!(integrator, false)
        return nothing
    end
    interrupt = _ -> SciMLBase.DiscreteCallback((u, t, integrator) -> armed[] && t > half,
        integrator -> throw(_Interrupted()); initialize=started!, save_positions=(false, false))
    directory = mktempdir()
    try
        uninterrupted = _without_warnings(() -> run(_checkpointed(configuration, joinpath(directory, "a")), (interrupt,)))
        if _failed(uninterrupted)
            err = _unwrapped(uninterrupted.error)
            unrestorable = sprint(showerror, SimulationEngine._missing_restore(subject))
            sprint(showerror, err) == unrestorable && return _check(id, since, false, unrestorable)
            # A refusal of checkpointed runs at setup, by the model or the engine,
            # before the solve started; any error after that fails the check.
            err isa ArgumentError && !solving[] && return _skip(id, since,
                "the fixture refuses to run with checkpoints: " * _summary(err))
            return _check(id, since, false, "the fixture failed with checkpoints every $(half) s: " * _summary(err))
        end
        armed[] = true
        interrupted = try
            _without_warnings(() -> run(_checkpointed(configuration, joinpath(directory, "b")), (interrupt,)))
        finally
            armed[] = false
        end
        _failed(interrupted) && _unwrapped(interrupted.error) isa _Interrupted || return _check(id, since, false,
            "a run interrupted after its checkpoint at $(half) s did not stop there: " *
            (_failed(interrupted) ? _summary(interrupted.error) : "it finished"))
        resumed = _without_warnings(() -> run(_checkpointed(configuration, joinpath(directory, "b"); resume=true), (interrupt,)))
        gap = difference(resumed, uninterrupted)
        return _check(id, since, gap === nothing,
            "a run resumed from its checkpoint at $(half) s and an uninterrupted run differ after it: $(gap); " *
            "save the run state with checkpoint_state and put it back with restore_checkpoint_state!")
    finally
        rm(directory; recursive=true, force=true)
    end
end

# The round trip of the guidance, navigation or control model `model`, found by
# identity among `configuration`'s models of `kind`: the states and the model's own
# state at the end of the second half.
function _gnc_checkpoint_round_trip(model, kind::Symbol, configuration)
    index = findfirst(effector -> effector === model, _container_effectors(configuration, kind))
    index === nothing && throw(ArgumentError("the configuration's $(_CONTAINER_NAMES[kind]) does not contain this model"))
    return _checkpoint_round_trip_check((cfg, factories) -> _model_solve(cfg, kind, index; factories),
        _model_run_difference, model, configuration)
end
