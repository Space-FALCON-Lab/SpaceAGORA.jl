import Distributed
import SHA

# This restart context belongs only to the disposable forecast process. Ordinary
# runs keep the existing file-checkpoint path and clock behavior.
const _CONTROL_FORECAST_RESTART = Base.ScopedValues.ScopedValue{Any}(nothing)
const _CONTROL_FORECAST_LOCK = ReentrantLock()

"""
    ControlForecastSnapshot

A copy-owned restart of one spacecraft immediately before its sole periodic
controller applies a command. Created by [`capture_control_forecast`](@ref).
The snapshot contains configuration and plain checkpoint state, not an
integrator, a borrowed parameter bundle, a worker, or previous output rows.
"""
struct ControlForecastSnapshot{A,C,K}
    configuration::A
    checkpoint::C
    clock_origin_s::Float64
    owner_pid::Int
    kernels::K
    environment::Dict{String,String}
end

_forecast_environment() = Dict(k=>v for (k,v) in ENV if startswith(k,"SPACEAGORA_"))

# Controller templates/state must be data, not borrowed services. This catches
# ordinary accidental retention of closures, native pointers or process handles;
# it is not a sandbox for untrusted user-defined Julia methods.
function _check_forecast_data(value,seen=IdDict{Any,Nothing}())
    value isa Union{Function,Task,IO,Ptr,Module,Type,Distributed.AbstractWorkerPool,Distributed.AbstractRemoteRef,ReentrantLock,AbstractChannel} &&
        throw(ArgumentError("Forecast controllers and checkpoint state cannot contain functions or runtime handles."))
    value === nothing || value isa Union{Number,AbstractString,Symbol,Enum} ? (return nothing) : nothing
    haskey(seen,value) && return nothing
    seen[value]=nothing
    if value isa AbstractArray
        value isa Union{Array,BitArray,SArray,MArray} || throw(ArgumentError(
            "Forecast controller arrays must own their storage (Array, BitArray, SArray or MArray)."))
        foreach(v->_check_forecast_data(v,seen),value)
    elseif value isa Union{Tuple,NamedTuple}
        foreach(v->_check_forecast_data(v,seen),value)
    elseif value isa AbstractDict
        for (k,v) in value
            _check_forecast_data(k,seen);_check_forecast_data(v,seen)
        end
    elseif isstructtype(typeof(value))
        for name in fieldnames(typeof(value))
            isdefined(value,name) && _check_forecast_data(getfield(value,name),seen)
        end
    else
        throw(ArgumentError("Unsupported forecast controller data type $(typeof(value))."))
    end
    nothing
end

function _control_forecast_kernels(args)
    needs_kernels = args.environment_model.ephemerides_model isa SimulationModel.SpiceEphemeridesModel ||
        any(e->e isa Union{SimulationModel.NBodyGravityModel,SimulationModel.SolarRadiationPressureModel},
            args.dynamics_model.dynamic_effectors)
    needs_kernels || return nothing
    return lock(RuntimeServices.SPICE_LOCK) do
        entries=NamedTuple{(:path,:kind,:sha256),Tuple{String,String,String}}[]
        for i in 1:SPICE.ktotal("ALL")
            file,kind,source,_=SPICE.kdata(i,"ALL",8192,8192)
            isempty(source) && kind!="META" || throw(ArgumentError("Control forecasts currently require directly loaded SPICE kernels, not meta-kernels."))
            path=abspath(file)
            isfile(path) || throw(ArgumentError("A loaded SPICE kernel cannot be located at $path."))
            digest=open(path,"r") do io;bytes2hex(SHA.sha256(io));end
            push!(entries,(path=path,kind=kind,sha256=digest))
        end
        isempty(entries) && throw(ArgumentError("The SPICE forecast has no loaded kernels."))
        entries
    end
end

function _initialize_control_forecast_kernels!(kernels)
    kernels===nothing && return nothing
    lock(RuntimeServices.SPICE_LOCK) do
        # Verify all assets before changing the disposable worker's kernel pool.
        for entry in kernels
            isfile(entry.path) || throw(ArgumentError("Missing forecast kernel $(entry.path)."))
            digest=open(entry.path,"r") do io;bytes2hex(SHA.sha256(io));end
            digest==entry.sha256 || throw(ArgumentError("Forecast kernel contents changed: $(entry.path)."))
        end
        SPICE.kclear()
        SimulationModel.Planets._reset_furnished_kernels!()
        # Reproduce order and duplicate loads: both affect SPICE precedence.
        # Managed furnishing updates planet caches and pool generation.
        for entry in kernels
            SimulationModel.Planets._furnsh_once(entry.path;force=true)
        end
    end
    # This engine revision's reuse keys do not include kernel generation.
    # A forecast must not reuse tables built from an earlier worker pool.
    _clear_ephemeris_reuse_cache!()
    nothing
end

function _check_control_forecast_configuration(args)
    length(args.dynamics_model.spacecraft) == 1 || throw(ArgumentError("Control forecasts currently require one spacecraft."))
    !args.mission_configuration.orientation_sim && !args.mission_configuration.keplerian ||
        throw(ArgumentError("Control forecasts require Cartesian orbit dynamics without attitude simulation."))
    args.interlink_model === nothing && isempty(args.external_propagators) ||
        throw(ArgumentError("Interlinks and external propagation are not supported by control forecasts."))
    isempty(args.guidance_model.guidance_effectors) && isempty(args.navigation_model.navigation_effectors) ||
        throw(ArgumentError("Control forecasts currently require no navigation or guidance effectors."))
    length(args.control_model.control_effectors) == length(args.control_model.control_rates) == 1 ||
        throw(ArgumentError("Control forecasts require exactly one periodic controller."))
    control = only(args.control_model.control_effectors)
    isempty(GNCRuntime.gnc_links(control)) || throw(ArgumentError("Linked control producers cannot be replayed by this forecast boundary."))
    sc = only(args.dynamics_model.spacecraft)
    sc.instant_actuation && isempty(sc.joints) && isempty(sc.attachments) &&
        isempty(SimulationModel.mounted_sensors(sc)) && isempty(SimulationModel.mounted_actuators(sc)) ||
        throw(ArgumentError("Control forecasts currently require instantaneous geometry without mounted devices or compliant attachments."))
    all(link->isempty(link.solar_arrays) && isempty(link.batteries) && isempty(link.power_loads) &&
        isempty(link.faults.faults) && isempty(link.data_stores) && isempty(link.data_sources),sc.links) ||
        throw(ArgumentError("Power devices, faults and onboard data models are not yet supported by control forecasts."))
    args.environment_model.density_model isa Union{SimulationModel.NoAtmosphereModel,SimulationModel.PolynomialFitAtmosphereModel} ||
        throw(ArgumentError("Control forecasts currently support only no-atmosphere and polynomial-fit atmospheres."))
    args.environment_model.ephemerides_model isa Union{SimulationModel.SimpleEphemeridesModel,SimulationModel.SpiceEphemeridesModel} ||
        throw(ArgumentError("Unsupported forecast ephemerides model."))
    args.environment_model.thermal_model isa SimulationModel.MaxwellianHeat ||
        throw(ArgumentError("Unsupported forecast thermal model."))
    allowed = Union{SimulationModel.InverseSquaredGravityModel,SimulationModel.InverseSquaredJ2GravityModel,
        SimulationModel.NBodyGravityModel,SimulationModel.GravitationalHarmonicsModel,
        SimulationModel.SolarRadiationPressureModel,SimulationModel.AerodynamicCoefficientfM}
    all(e -> e isa allowed, args.dynamics_model.dynamic_effectors) ||
        throw(ArgumentError("An unsupported dynamic effector is present in the control forecast."))
    return control
end

"""
    capture_control_forecast(u, p, t, controller) -> ControlForecastSnapshot

Capture a forecast at the beginning of the sole controller's `update_control!`,
before it changes its command, plan, phase, panel geometry or histories. This
initial interface supports one spacecraft and one periodic controller, with no
linked producers, navigation, guidance, mounted devices, or attitude dynamics.
Unsupported models fail before copying state. The controller must expose its
previous command time through `held_command_time` and implement checkpoint hooks.

The current control tick is still pending in the forecast, although core has
already advanced its clock index before entering `update_control!`. This is a
pre-command forecast boundary, not a general way to write a checkpoint midway
through arbitrary callback code. The caller must not have partially updated its
controller. Both the absolute simulation time and original ephemeris epoch are
preserved. Evaluate the snapshot only in a separate process.

SPICE forecasts require directly loaded, unchanged kernel files. Meta-kernels
and programmatic kernel-pool assignments are unsupported. Scoped engine overrides
are rejected; resolved solver settings and SPACEAGORA environment values are
frozen. Controller state must contain ordinary data, without closures, remote
references or other live runtime handles.
"""
function capture_control_forecast(u, p, t::Real, controller)
    _CONTROL_FORECAST_RESTART[]===nothing || throw(ArgumentError("A candidate cannot capture a nested forecast."))
    _engine_active_overrides_ref[]===nothing || throw(ArgumentError("Scoped engine overrides are not yet supported by control forecasts."))
    args = p.args
    _check_control_forecast_configuration(args) === controller ||
        throw(ArgumentError("The controller must be the live run's sole controller."))
    time = Float64(t)
    isfinite(time) && time >= 0 && all(isfinite,u) || throw(ArgumentError("Forecast state and time must be finite."))
    command_time = GNCRuntime.held_command_time(controller,1)
    command_time isa Real && !isnan(command_time) && command_time < time ||
        throw(ArgumentError("Capture must precede this tick's command; the previous held-command time is required."))
    length(p.shared_buffers.periodic_clocks) == 1 || throw(ArgumentError("Expected one initialized control clock."))
    clock = only(p.shared_buffers.periodic_clocks)
    period = only(args.control_model.control_rates)
    isfinite(clock.origin[]) && clock.period_s == period &&
        SimulationModel.SimulationCallbacks._same_run_time(clock.origin[]+(clock.index[]-1)*period,time) ||
        throw(ArgumentError("Capture must occur inside the active periodic control tick."))
    mode = _run_solver_config(args).solver_mode
    mode in (:split_imex,:auto_stiff) || throw(ArgumentError("Control forecasts currently support split_imex and auto_stiff."))
    _check_model_checkpoint_hooks(args)
    _check_forecast_data(controller)
    SimulationLifecycle.checkpoint_state(controller) === nothing &&
        throw(ArgumentError("The forecast controller must preserve its state through checkpoint hooks."))
    runtime = (activity=_activity_checkpoint_state(p.is_active),
        callbacks=_callback_checkpoint_state(args), models=_model_checkpoint_state(args),
        periodic_ticks=(version=2,time_s=time,pending=[(clock=1,period_s=Float64(period))]),
        gnc_scopes=_gnc_scope_record(args,p.shared_buffers.gnc_scopes[]),
        effector_scopes=_effector_scope_checkpoint_state(args),
        link_attitudes=_link_attitude_checkpoint_state(args),
        maneuver_buffers=_maneuver_buffer_checkpoint_state(p),
        orbits=_orbit_counter_checkpoint_state(p),entries=_entry_counter_checkpoint_state(p))
    checkpoint = (t=time,u=deepcopy(_checkpoint_state(args,u)),solver_mode=string(mode),runtime_state=deepcopy(runtime))
    _check_forecast_data(runtime.models)
    owned = SimulationModel.SimConfig._with_configuration(deepcopy(args); solver_config=deepcopy(_run_solver_config(args)),
        simulation_settings=SimulationModel.SimulationSettings(results=false,verbose=false,generate_plots=false))
    return ControlForecastSnapshot(owned,checkpoint,Float64(clock.origin[]),Int(getpid()),
        _control_forecast_kernels(args),_forecast_environment())
end

"""
    run_control_forecast(snapshot, candidate; duration_s, sample_period_s)

Evaluate a copied candidate controller from a pre-command snapshot in a dedicated
forecast process, preferably created by [`with_control_forecast_worker`](@ref).
Do not run unrelated simulations in that process. Forecast calls are serialized
there because SPICE owns process-wide state. Nested forecasts are rejected.
The candidate replaces the snapshot's sole controller and must preserve
its own state through checkpoint hooks. It must keep the intended plan fixed
rather than recursively invoking this service. Its spacecraft scope and control
period must match the snapshot. Candidate constraints and failure policies remain
the caller's responsibility.

Returns copied saved rows, the final state/time and solver return code. It returns
no integrator or live runtime handles and writes no results or checkpoint files.
The saved time coordinate remains the original run's absolute elapsed time.
"""
function run_control_forecast(snapshot::ControlForecastSnapshot,candidate;
        duration_s::Real,sample_period_s::Real)
    _CONTROL_FORECAST_RESTART[]===nothing || throw(ArgumentError("A candidate cannot recursively invoke a forecast."))
    return lock(_CONTROL_FORECAST_LOCK) do
        _run_control_forecast(snapshot,candidate;duration_s,sample_period_s)
    end
end

function _run_control_forecast(snapshot::ControlForecastSnapshot,candidate;
        duration_s::Real,sample_period_s::Real)
    getpid() != snapshot.owner_pid || throw(ArgumentError("A control forecast must run in a separate process."))
    duration, sampling = Float64(duration_s),Float64(sample_period_s)
    isfinite(duration) && duration>0 && isfinite(sampling) && sampling>0 ||
        throw(ArgumentError("Forecast duration and sample period must be finite and positive."))
    isfinite(snapshot.checkpoint.t+duration) || throw(ArgumentError("Forecast end time must be finite."))
    _forecast_environment()==snapshot.environment || throw(ArgumentError("SPACEAGORA environment settings differ from the source forecast snapshot."))
    _engine_active_overrides_ref[]===nothing || throw(ArgumentError("Forecast evaluation cannot use scoped engine overrides."))
    _check_forecast_data(candidate)
    args = deepcopy(snapshot.configuration)
    model = deepcopy(candidate)
    controls = SimulationModel.ControlModel(control_effectors=(model,),control_rates=copy(args.control_model.control_rates))
    mission = SimulationModel.SimConfig._with_fields(args.mission_configuration;
        mission_type=SimulationModel.MissionTime,mission_time=snapshot.checkpoint.t+duration,data_rate=sampling)
    args = SimulationModel.SimConfig._with_configuration(args;control_model=controls,mission_configuration=mission)
    _check_control_forecast_configuration(args)
    _check_model_checkpoint_hooks(args)
    SimulationLifecycle.checkpoint_state(model) === nothing && throw(ArgumentError("Candidate checkpoint state is required."))
    # Only the copied controller's state changes. Continuous state, geometry,
    # activity, counters and pending ticks come from the common source snapshot.
    runtime = merge(deepcopy(snapshot.checkpoint.runtime_state),(models=_model_checkpoint_state(args),))
    _check_forecast_data(runtime.models)
    checkpoint = merge(deepcopy(snapshot.checkpoint),(runtime_state=runtime,))
    restart = (checkpoint=checkpoint,clock_origin_s=snapshot.clock_origin_s)
    _initialize_control_forecast_kernels!(snapshot.kernels)
    out = Base.ScopedValues.with(_CONTROL_FORECAST_RESTART=>restart) do
        run_simulation(args;return_results=true,return_solution=true,visualization=false)
    end
    out.solution === nothing && error("The forecast did not propagate")
    return (table=deepcopy(out.table),final_state=copy(out.solution.u[end]),
        final_time_s=Float64(out.solution.t[end]),retcode=string(out.solution.retcode),
        start_time_s=snapshot.checkpoint.t,clock_origin_s=snapshot.clock_origin_s)
end

"""
    with_control_forecast_worker(f, snapshot; packages=())

Run `f(evaluate)` with a dedicated local Julia worker, where
`evaluate(candidate; duration_s, sample_period_s)` calls [`run_control_forecast`](@ref).
The worker uses the current Julia/project and is removed on success or failure.
`packages` contains the `Base.PkgId`s of companion packages defining candidate
controllers; core is loaded automatically. No existing worker is adopted.
"""
function with_control_forecast_worker(f,snapshot::ControlForecastSnapshot;packages=())
    _CONTROL_FORECAST_RESTART[]===nothing || throw(ArgumentError("A candidate cannot create a nested forecast worker."))
    project=Base.active_project()
    project===nothing && throw(ArgumentError("An active project is required for a forecast worker."))
    workers=Distributed.addprocs(1;exeflags=`--startup-file=no --project=$(dirname(project)) --threads=1`,
        env=collect(snapshot.environment))
    worker=only(workers)
    try
        Distributed.remotecall_fetch(Base.require,worker,Base.PkgId(parentmodule(@__MODULE__)))
        for package in packages
            package isa Base.PkgId || throw(ArgumentError("Forecast packages must be Base.PkgId values."))
            Distributed.remotecall_fetch(Base.require,worker,package)
        end
        evaluate=(candidate;kwargs...)->Distributed.remotecall_fetch(
            run_control_forecast,worker,snapshot,candidate;kwargs...)
        return f(evaluate)
    finally
        Distributed.rmprocs(workers)
    end
end
