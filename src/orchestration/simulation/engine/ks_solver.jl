# Optional numerical backend; physical state, force models and lifecycle stay here.
"""
    trial_rhs_safe(model) -> Bool

Declare that a model's RHS evaluation is repeatable for perturbed trial states
and times without advancing physical, random or controller history. This is an
interface contract, not scientific validation. Unknown models default to false.
Accepted-step callbacks may still update discrete commands normally.
"""
trial_rhs_safe(::Any) = false
trial_rhs_safe(::Union{SimulationModel.InverseSquaredGravityModel,
    SimulationModel.InverseSquaredJ2GravityModel,SimulationModel.AerodynamicCoefficientConstant,
    SimulationModel.MaxwellianHeat}) = true
trial_rhs_safe(model::SimulationModel.ThrusterModels.BaseThrusterModel) =
    SimulationModel.actuation_mode(model) === :ideal

function _ks_extension()
    root = parentmodule(@__MODULE__)
    extension = Base.get_extension(root, :SpaceAGORAKSPropagatorExt)
    if extension === nothing
        pkg = Base.PkgId(Base.UUID("d8aef2b7-6f1f-4b9e-9f52-5b5a2773f0d2"), "KSPropagator")
        try
            Base.require(pkg)
        catch err
            err isa ArgumentError || rethrow()
            throw(ArgumentError("solver_mode=:ks requires the optional KSPropagator package with coupled_cartesian_step_doubling. Install the documented compatible source; the ordinary solver remains available. Original error: $(sprint(showerror,err))"))
        end
        extension = Base.get_extension(root, :SpaceAGORAKSPropagatorExt)
    end
    extension === nothing && throw(ArgumentError("The optional KS solver extension could not load; check KSPropagator and OrdinaryDiffEqCore compatibility."))
    return extension
end

# Not checked on purpose: guidance, navigation and scheduling models act only
# through accepted-step callbacks, never inside the right-hand side, and SRP
# facets are read only by SRP models, which the trial-safety walk refuses.
function _validate_ks_config!(args)
    length(args.dynamics_model.spacecraft)==1 || throw(ArgumentError("The initial KS route supports exactly one spacecraft."))
    args.mission_configuration.orientation_sim && throw(ArgumentError("The initial KS route does not yet qualify propagated attitude; use orientation_sim=false."))
    isempty(args.external_propagators) || throw(ArgumentError("KS cannot be combined with shadow external propagators."))
    args.interlink_model === nothing || throw(ArgumentError("The initial KS route does not support interlinks."))
    density = args.environment_model.density_model
    density isa Union{SimulationModel.NoAtmosphereModel,SimulationModel.ExponentialAtmosphereModel} ||
        throw(ArgumentError("The initial KS route supports NoAtmosphereModel or deterministic ExponentialAtmosphereModel."))
    isfinite(args.environment_model.planet.μ) && args.environment_model.planet.μ>0 ||
        throw(ArgumentError("KS requires a positive finite central gravitational parameter."))
    trial_rhs_safe(args.environment_model.thermal_model) || throw(ArgumentError("KS requires a repeatable thermal model declaring trial_rhs_safe."))
    args.environment_model.ephemerides_model isa SimulationModel.SimpleEphemeridesModel ||
        throw(ArgumentError("The initial KS qualification supports SimpleEphemeridesModel only."))
    args.environment_model.topography && throw(ArgumentError("Topography is outside the initial KS qualification."))
    sc = only(args.dynamics_model.spacecraft)
    length(sc.links)==1 && isempty(sc.joints) && isempty(sc.attachments) ||
        throw(ArgumentError("The initial KS route supports a single rigid spacecraft link."))
    link = sc.root
    for field in (:thrusters,:magnets,:magnetorquers,:laser_terminals,:hinges,:gimbals,
                  :custom_actuators,:solar_arrays,:batteries,:power_loads,:sensors,
                  :data_stores,:data_sources)
        isempty(getproperty(link,field)) || throw(ArgumentError("Mounted $field are outside the initial KS qualification."))
    end
    link.rw_assembly.n_wheels==0 || throw(ArgumentError("Reaction wheels are outside the initial KS qualification."))
    isempty(link.rw_assembly.cmgs) && isempty(link.rw_assembly.cmg_inner) ||
        throw(ArgumentError("Control moment gyros are outside the initial KS qualification."))
    isempty(link.faults.faults) || throw(ArgumentError("Actuator faults are outside the initial KS qualification."))
    effects = args.dynamics_model.dynamic_effectors
    scopes = _effector_scopes_of(args)
    central = count(eachindex(effects)) do index
        effects[index] isa Union{SimulationModel.InverseSquaredGravityModel,
            SimulationModel.InverseSquaredJ2GravityModel} && _effector_applies(scopes,index,1)
    end
    central==1 || throw(ArgumentError("KS requires exactly one configured inverse-square gravity model (with or without J2) applying to spacecraft 1, avoiding missing or doubled central gravity."))
    for model in (effects...,args.control_model.control_effectors...)
        trial_rhs_safe(model) || throw(ArgumentError("KS requires repeatable trial RHS evaluations; $(typeof(model)) has not declared trial_rhs_safe."))
    end
    extension = _ks_extension()
    Base.invokelatest(extension.require_capability)
    return nothing
end

function _validate_ks_runtime!(p)
    cb = SimulationModel.SimulationCallbacks._callback_env_config(p)
    isempty(p.shared_buffers.density_models) || throw(ArgumentError("KS trial evaluation admits only analytic atmospheres; this run created per-run density model instances."))
    cb.density_history_dependent && throw(ArgumentError("KS trial evaluation does not support history-dependent density."))
    density = p.args.environment_model.density_model
    if !(density isa SimulationModel.NoAtmosphereModel)
        cb.vacuum_gram_cache_enabled && throw(ArgumentError("KS initial analytic-atmosphere qualification requires vacuum GRAM caching off."))
        cb.density_freeze_per_step && throw(ArgumentError("KS with a state-dependent atmosphere requires SPACEAGORA_DENSITY_FREEZE_PER_STEP=false; frozen samples have different semantics."))
        cb.gram_track_cache.mode === :off || throw(ArgumentError("KS initial analytic-atmosphere qualification requires GRAM track caching off."))
    end
    return nothing
end
