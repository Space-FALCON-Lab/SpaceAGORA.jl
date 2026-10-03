# Purpose: Check whether any effector is one of the supported aerodynamic models.
# Input: Tuple of dynamic effectors.
# Output: true if an aerodynamic model is present; otherwise false.
@inline function _uses_atmospheric_dynamic_effector(effectors::Tuple)::Bool
    @inbounds for effector in effectors
        if effector isa AerodynamicCoefficientConstant || effector isa AerodynamicCoefficientfM || effector isa AerodynamicCoefficientNoBallisticFlight
            return true
        end
    end
    return false
end

# Purpose: Check whether the effectors include the inverse-square-plus-J2 gravity model.
# Input: Tuple of dynamic effectors.
# Output: true if that gravity model is present; otherwise false.
@inline function _uses_j2_gravity_effector(effectors::Tuple)::Bool
    @inbounds for effector in effectors
        if effector isa InverseSquaredJ2GravityModel
            return true
        end
    end
    return false
end

# Purpose: Decide whether the dynamics setup calls for atmospheric density.
# Input: Dynamic effectors and the simulation configuration.
# Output: true if an aerodynamic effector or a non-NoAtmosphere density model is present.
@inline function _requires_density_for_rhs(effectors::Tuple, args::SimulationConfiguration)::Bool
    return _uses_atmospheric_dynamic_effector(effectors) || !(args.environment_model.density_model isa NoAtmosphereModel)
end

# Backward-compat alias — existing callers that check whether density computation
# is possible at all (thermal, drag-state, entry-end, tests) keep working.
# Purpose: Keep the older density-check name working by forwarding to the current check.
# Input: Dynamic effectors and the simulation configuration.
# Output: The same Boolean result as _requires_density_for_rhs.
@inline function _requires_density_callback(effectors::Tuple, args::SimulationConfiguration)::Bool
    return _requires_density_for_rhs(effectors, args)
end

# Returns true only when a non-RHS consumer needs density pre-staged in
# shared_buffers before the step.  The RHS computes density inline via
# sample_buffered_atmosphere → sample_atmosphere, so it does not require the
# staged callback.  Callers: get_callbacks density callback installation only.
# Purpose: Decide whether callbacks need density prepared in shared buffers.
# Input: Effectors, simulation configuration, and the SPACEAGORA_FORCE_DENSITY_CALLBACK override.
# Output: true if density is needed and a consumer or override requires the staging callback.
@inline function _requires_staged_density_callback(effectors::Tuple, args::SimulationConfiguration)::Bool
    _requires_density_for_rhs(effectors, args) || return false
    # Explicit compatibility/debug override.
    ParallelPolicy.parse_bool_env("SPACEAGORA_FORCE_DENSITY_CALLBACK", false) && return true
    # Thermal callback fires at each step and reads shared_buffers.densities.
    _requires_thermal_callback(effectors, args) && return true
    # Drag-state switching uses staged density to pick the next tolerance set.
    _requires_drag_state_callback(effectors, args) && return true
    # Entry-end detection compares altitude against EI using staged density.
    _requires_entry_end_callback(effectors, args) && return true
    return false
end

# Purpose: Check whether any guidance model needs to count completed orbits.
# Input: Simulation configuration containing the guidance effectors.
# Output: true if a guidance effector has a maneuver_orbit_number property.
@inline function _requires_guidance_orbit_counter(args::SimulationConfiguration)::Bool
    @inbounds for guidance_model in args.guidance_model.guidance_effectors
        hasproperty(guidance_model, :maneuver_orbit_number) && return true
    end
    return false
end

# Purpose: Decide whether to install the callback that tracks completed orbits.
# Input: Simulation configuration with mission and guidance settings.
# Output: true for an orbit-count mission or guidance that needs an orbit counter.
@inline function _requires_orbit_end_callback(args::SimulationConfiguration)::Bool
    return args.mission_configuration.mission_type == MissionOrbits ||
           _requires_guidance_orbit_counter(args)
end

# Purpose: Read and validate the requested number of atmospheric-entry targets.
# Input: No arguments; reads SPACEAGORA_ENTRY_TARGET_COUNT, defaulting to "0".
# Output: A nonnegative Int; throws ArgumentError for an invalid or negative value.
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

# Purpose: Decide whether to install the atmospheric-entry completion callback.
# Input: Effectors, simulation configuration, and the entry-target environment setting.
# Output: true if the target count is positive and the density check passes.
@inline function _requires_entry_end_callback(effectors::Tuple, args::SimulationConfiguration)::Bool
    return _entry_target_count() > 0 && _requires_density_callback(effectors, args)
end

# Purpose: Decide whether atmospheric flight needs a callback to switch solver limits.
# Input: Effectors and configuration containing atmospheric and orbital tolerances.
# Output: true if density is needed and the phase step limits or error tolerances differ.
@inline function _requires_drag_state_callback(effectors::Tuple, args::SimulationConfiguration)::Bool
    if !_requires_density_callback(effectors, args)
        return false
    end
    tol = args.integration_tolerances
    return tol.dt_max_atmosphere != tol.dt_max_orbit ||
           tol.reltol_atmosphere != tol.reltol_orbit ||
           tol.abstol_atmosphere != tol.abstol_orbit
end

# Purpose: Decide whether to install a callback for spacecraft heating.
# Input: Effectors and configuration containing the density model and spacecraft.
# Output: true if the density check passes and at least one spacecraft has links.
@inline function _requires_thermal_callback(effectors::Tuple, args::SimulationConfiguration)::Bool
    _requires_density_callback(effectors, args) || return false
    @inbounds for spacecraft in args.dynamics_model.spacecraft
        !isempty(spacecraft.links) && return true
    end
    return false
end

# Purpose: Decide whether attitude quaternions need a normalization callback.
# Input: Simulation configuration containing the orientation_sim setting.
# Output: true when attitude simulation is enabled; otherwise false.
@inline function _requires_quaternion_projection_callback(args::SimulationConfiguration)::Bool
    return args.mission_configuration.orientation_sim
end

# Purpose: Use the baseline tolerance when a component tolerance is left at zero.
# Input: Component-specific tolerance and baseline tolerance, both Float64.
# Output: The baseline for a zero component tolerance; otherwise the component tolerance.
@inline function _resolved_component_tolerance(component_tol::Float64, baseline_tol::Float64)::Float64
    return component_tol == 0.0 ? baseline_tol : component_tol
end

# Purpose: Build error tolerances for atmospheric or orbital flight without changing the templates.
# Input: Relative/absolute tolerance templates, simulation configuration, and the flight-phase flag.
# Output: A (reltol, abstol) tuple of scalars or copied state-shaped tolerances with component overrides.
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

# Purpose: Add one callback to the end of a callback tuple.
# Input: Existing callback tuple and one additional callback.
# Output: A tuple containing the original callbacks followed by the new callback.
@inline _append_callback(callbacks::Tuple, callback) = (callbacks..., callback)

# Purpose: Leave the callback list unchanged when there is no callback to add.
# Input: Existing callback tuple and nothing.
# Output: The original callback tuple.
@inline _append_callback(callbacks::Tuple, ::Nothing) = callbacks

# Purpose: Add callbacks supplied as a tuple to the existing callback list.
# Input: Existing callback tuple and a tuple of additional callbacks.
# Output: A combined tuple preserving the order of both inputs.
@inline _append_callbacks(callbacks::Tuple, extra::Tuple) = (callbacks..., extra...)

# Purpose: Add callbacks supplied as a vector to the existing callback list.
# Input: Existing callback tuple and a vector of additional callbacks.
# Output: A combined tuple preserving the order of both inputs.
@inline _append_callbacks(callbacks::Tuple, extra::AbstractVector) = (callbacks..., extra...)

# Purpose: Assemble the event handlers needed by the configured simulation.
# Input: Spacecraft count, effectors, configuration, optional saving inputs and extra callbacks; reads solver policy.
# Output: A CallbackSet with applicable built-in callbacks followed by the supplied extra callbacks.
function get_callbacks(
    num_sats::Int,
    effectors::Tuple,
    args::SimulationConfiguration;
    saved_values=nothing,
    save_fields=nothing,
    extra_callbacks=()
)::CallbackSet
    save_fields_resolved = _resolve_save_fields(save_fields, args)
    backbone_mode = _simulation_engine_module()._solver_policy_mode() == :gravity_backbone_split
    callbacks = if backbone_mode
        (get_impact_callback(num_sats),)
    else
        (
            get_impact_callback(num_sats),
            update_planet_frame_callback(),
        )
    end

    # If we are not using the gravity-backbone solver, and we need a callback to prepare atmospheric density, add that callback to the list.
    if !backbone_mode && _requires_staged_density_callback(effectors, args)
        callbacks = _append_callback(callbacks, get_density_callback(num_sats, effectors, args))
    end 

    # If we are not using the gravity-backbone solver, and we need a callback to handle thermal effects, add that callback to the list.
    if !backbone_mode && _requires_thermal_callback(effectors, args)
        callbacks = _append_callback(callbacks, get_thermal_callback(num_sats, args))
    end

    # If we need a callback to detect the end of the orbit, add that callback to the list.
    if _requires_orbit_end_callback(args)
        callbacks = _append_callback(callbacks, get_orbit_end_callback(num_sats))
    end

    # If we are not using the gravity-backbone solver, and we need a callback to detect the end of the entry phase, add that callback to the list.
    if !backbone_mode && _requires_entry_end_callback(effectors, args)
        callbacks = _append_callback(callbacks, get_entry_end_callback(num_sats, args))
    end

    # If we are not using the gravity-backbone solver, and we need a callback to handle drag state, add that callback to the list.
    if !backbone_mode && _requires_drag_state_callback(effectors, args)
        callbacks = _append_callback(callbacks, get_drag_state_callback(num_sats))
    end

    # If we are not using the gravity-backbone solver, add navigation, control, and guidance callbacks.
    if !backbone_mode
        callbacks = _append_callbacks(callbacks, get_navigation_callbacks(num_sats, args))
        callbacks = _append_callbacks(callbacks, get_control_callbacks(num_sats, args))
        callbacks = _append_callbacks(callbacks, get_guidance_callbacks(num_sats, args))
    end

    # If we are not using the gravity-backbone solver, and we need a callback to project quaternions, add that callback to the list.
    if !backbone_mode && _requires_quaternion_projection_callback(args)
        callbacks = _append_callback(callbacks, get_quaternion_projection_callback(num_sats, args))
    end
    # If we are not using the gravity-backbone solver, and result output is enabled, add the data saving callback to the list.
    if !backbone_mode && args.simulation_settings.results
        callbacks = _append_callback(callbacks, get_data_saving_callback(num_sats, args, save_fields_resolved, saved_values))
    end
    callbacks = _append_callbacks(callbacks, extra_callbacks)

    return CallbackSet(callbacks...)
end
