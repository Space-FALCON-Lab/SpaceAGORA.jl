@inline function _uses_atmospheric_dynamic_effector(effectors::Tuple)::Bool
    @inbounds for effector in effectors
        if effector isa AerodynamicCoefficientConstant || effector isa AerodynamicCoefficientfM || effector isa AerodynamicCoefficientNoBallisticFlight
            return true
        end
    end
    return false
end

@inline function _uses_j2_gravity_effector(effectors::Tuple)::Bool
    @inbounds for effector in effectors
        if effector isa InverseSquaredJ2GravityModel
            return true
        end
    end
    return false
end

@inline function _requires_density_callback(effectors::Tuple, args::SimulationConfiguration)::Bool
    return _uses_atmospheric_dynamic_effector(effectors) || !(args.environment_model.density_model isa NoAtmosphereModel)
end

@inline function _requires_guidance_orbit_counter(args::SimulationConfiguration)::Bool
    @inbounds for guidance_model in args.guidance_model.guidance_effectors
        hasproperty(guidance_model, :maneuver_orbit_number) && return true
    end
    return false
end

@inline function _requires_orbit_end_callback(args::SimulationConfiguration)::Bool
    return args.mission_configuration.mission_type == MissionOrbits ||
           _requires_guidance_orbit_counter(args)
end

@inline function _entry_target_count()::Int
    raw = strip(get(ENV, "SPACEAGORA_ENTRY_TARGET_COUNT", "0"))
    parsed = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_ENTRY_TARGET_COUNT must be an integer value, got '$raw'"))
    end
    parsed >= 0 || throw(ArgumentError("SPACEAGORA_ENTRY_TARGET_COUNT must be >= 0, got $parsed"))
    return parsed
end

@inline function _requires_entry_end_callback(effectors::Tuple, args::SimulationConfiguration)::Bool
    return _entry_target_count() > 0 && _requires_density_callback(effectors, args)
end

@inline function _requires_drag_state_callback(effectors::Tuple, args::SimulationConfiguration)::Bool
    if !_requires_density_callback(effectors, args)
        return false
    end
    tol = args.integration_tolerances
    return tol.dt_max_atmosphere != tol.dt_max_orbit ||
           tol.reltol_atmosphere != tol.reltol_orbit ||
           tol.abstol_atmosphere != tol.abstol_orbit
end

@inline function _requires_thermal_callback(effectors::Tuple, args::SimulationConfiguration)::Bool
    # Thermal-rate evaluation requires atmospheric state (ρ, T, wind), populated by density callback.
    return _requires_density_callback(effectors, args)
end

@inline function _requires_quaternion_projection_callback(args::SimulationConfiguration)::Bool
    return args.mission_configuration.orientation_sim
end

@inline function _resolved_component_tolerance(component_tol::Float64, baseline_tol::Float64)::Float64
    return component_tol == 0.0 ? baseline_tol : component_tol
end

function _callback_tolerances_for_phase(template_reltol, template_abstol, args::SimulationConfiguration, in_atmosphere::Bool)
    tol = args.integration_tolerances
    baseline_reltol = in_atmosphere ? tol.reltol_atmosphere : tol.reltol_orbit
    baseline_abstol = in_atmosphere ? tol.abstol_atmosphere : tol.abstol_orbit
    if template_reltol isa Number && template_abstol isa Number
        return baseline_reltol, baseline_abstol
    end

    reltol_mass = _resolved_component_tolerance(tol.reltol_mass, baseline_reltol)
    abstol_mass = _resolved_component_tolerance(tol.abstol_mass, baseline_abstol)
    reltol_heat = _resolved_component_tolerance(tol.reltol_heat_load, baseline_reltol)
    abstol_heat = _resolved_component_tolerance(tol.abstol_heat_load, baseline_abstol)
    reltol_ω = _resolved_component_tolerance(tol.reltol_angular_rate, baseline_reltol)
    abstol_ω = _resolved_component_tolerance(tol.abstol_angular_rate, baseline_abstol)

    reltol_new = copy(template_reltol)
    abstol_new = copy(template_abstol)
    reltol_new .= baseline_reltol
    abstol_new .= baseline_abstol

    @inbounds for i in eachindex(reltol_new.sc)
        reltol_new.sc[i].mass = reltol_mass
        abstol_new.sc[i].mass = abstol_mass
        reltol_new.sc[i].heat_loads .= reltol_heat
        abstol_new.sc[i].heat_loads .= abstol_heat
        if hasproperty(reltol_new.sc[i], :q)
            reltol_new.sc[i].q .= tol.reltol_quaternion
            abstol_new.sc[i].q .= tol.abstol_quaternion
        end
        if hasproperty(reltol_new.sc[i], :ω)
            reltol_new.sc[i].ω .= reltol_ω
            abstol_new.sc[i].ω .= abstol_ω
        end
    end
    return reltol_new, abstol_new
end

@inline _append_callback(callbacks::Tuple, callback) = (callbacks..., callback)
@inline _append_callback(callbacks::Tuple, ::Nothing) = callbacks
@inline _append_callbacks(callbacks::Tuple, extra::Tuple) = (callbacks..., extra...)
@inline _append_callbacks(callbacks::Tuple, extra::AbstractVector) = (callbacks..., extra...)

function get_callbacks(
    num_sats::Int,
    effectors::Tuple,
    args::SimulationConfiguration;
    saved_values=nothing,
    save_fields=nothing
)::CallbackSet
    save_fields_resolved = _resolve_save_fields(save_fields, args)
    callbacks = (
        get_impact_callback(num_sats),
        update_planet_frame_callback(),
    )

    if _requires_density_callback(effectors, args)
        callbacks = _append_callback(callbacks, get_density_callback(num_sats, effectors, args))
    end

    if _requires_thermal_callback(effectors, args)
        callbacks = _append_callback(callbacks, get_thermal_callback(num_sats, args))
    end

    if _requires_orbit_end_callback(args)
        callbacks = _append_callback(callbacks, get_orbit_end_callback(num_sats))
    end

    if _requires_entry_end_callback(effectors, args)
        callbacks = _append_callback(callbacks, get_entry_end_callback(num_sats, args))
    end

    if _requires_drag_state_callback(effectors, args)
        callbacks = _append_callback(callbacks, get_drag_state_callback(num_sats))
    end

    callbacks = _append_callbacks(callbacks, get_navigation_callbacks(num_sats, args))
    callbacks = _append_callbacks(callbacks, get_control_callbacks(num_sats, args))
    callbacks = _append_callbacks(callbacks, get_guidance_callbacks(num_sats, args))
    if _requires_quaternion_projection_callback(args)
        callbacks = _append_callback(callbacks, get_quaternion_projection_callback(num_sats, args))
    end
    if args.simulation_settings.verbose && _progress_interval_s() > 0.0
        callbacks = _append_callback(callbacks, get_progress_callback(num_sats, args))
    end
    callbacks = _append_callback(callbacks, get_data_saving_callback(num_sats, args, save_fields_resolved, saved_values))

    return CallbackSet(callbacks...)
end
