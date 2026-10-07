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

@inline function _requires_density_for_rhs(effectors::Tuple, args::SimulationConfiguration)::Bool
    return _uses_atmospheric_dynamic_effector(effectors) || !(args.environment_model.density_model isa NoAtmosphereModel)
end

# Backward-compat alias — existing callers that check whether density computation
# is possible at all (thermal, drag-state, entry-end, tests) keep working.
@inline function _requires_density_callback(effectors::Tuple, args::SimulationConfiguration)::Bool
    return _requires_density_for_rhs(effectors, args)
end

# Returns true only when a non-RHS consumer needs density pre-staged in
# shared_buffers before the step.  The RHS computes density inline via
# sample_buffered_atmosphere → sample_atmosphere, so it does not require the
# staged callback.  Callers: get_callbacks density callback installation only.
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
    # Stateful guidance/control needs crossings even when both solver phases
    # use equal tolerances or the configured density is identically zero.
    any(SimulationLifecycle.requires_atmosphere_events, args.guidance_model.guidance_effectors) && return true
    any(SimulationLifecycle.requires_atmosphere_events, args.control_model.control_effectors) && return true
    if !_requires_density_callback(effectors, args)
        return false
    end
    tol = args.integration_tolerances
    return tol.dt_max_atmosphere != tol.dt_max_orbit ||
           tol.reltol_atmosphere != tol.reltol_orbit ||
           tol.abstol_atmosphere != tol.abstol_orbit
end

@inline function _requires_thermal_callback(effectors::Tuple, args::SimulationConfiguration)::Bool
    _requires_density_callback(effectors, args) || return false
    @inbounds for spacecraft in args.dynamics_model.spacecraft
        !isempty(spacecraft.links) && return true
    end
    return false
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
        if hasproperty(reltol_new.sc[i], :joint_q)
            reltol_new.sc[i].joint_q .= tol.reltol_quaternion
            abstol_new.sc[i].joint_q .= tol.abstol_quaternion
            reltol_new.sc[i].joint_qd .= reltol_ω
            abstol_new.sc[i].joint_qd .= abstol_ω
        end
        if hasproperty(reltol_new.sc[i], :att_q)
            reltol_new.sc[i].att_q .= tol.reltol_quaternion
            abstol_new.sc[i].att_q .= tol.abstol_quaternion
            reltol_new.sc[i].att_ω .= reltol_ω
            abstol_new.sc[i].att_ω .= abstol_ω
        end
    end
    return reltol_new, abstol_new
end

# An integrator advances every active spacecraft with one cap/tolerance set.
# Keep the atmospheric phase until the last active member leaves it.
@inline function _active_atmospheric_phase(p)::Bool
    return any(i -> p.is_active[i] && p.shared_buffers.in_atmosphere[i], eachindex(p.is_active))
end

function _active_phase_solver_settings(p, reltol, abstol)
    inside = _active_atmospheric_phase(p)
    tol = p.args.integration_tolerances
    cap = inside ? tol.dt_max_atmosphere : tol.dt_max_orbit
    reltol, abstol = _callback_tolerances_for_phase(reltol, abstol, p.args, inside)
    return cap, reltol, abstol
end

function _apply_active_phase_solver_settings!(integrator)
    p = integrator.p
    _requires_density_callback(p.args.dynamics_model.dynamic_effectors, p.args) || return nothing
    # Fixed-step symplectic/backbone drivers retain their prescribed step.
    hasproperty(integrator.opts, :adaptive) && !integrator.opts.adaptive && return nothing
    cap, reltol, abstol = _active_phase_solver_settings(p, integrator.opts.reltol, integrator.opts.abstol)
    integrator.opts.dtmax = cap
    integrator.opts.reltol = reltol
    integrator.opts.abstol = abstol
    return nothing
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
    save_fields=nothing,
    extra_callbacks=(),
    record_saved_values::Bool=false
)::CallbackSet
    save_fields_resolved = _resolve_save_fields(save_fields, args)
    backbone_mode = _simulation_engine_module()._solver_policy_mode() == :gravity_backbone_split
    touchdown_specs = _touchdown_specs(args, num_sats)
    has_touchdown = any(spec -> spec !== nothing, touchdown_specs)
    # Touchdown spacecraft and externally propagated shadow entries have no impact event.
    impact_excluded = has_touchdown ? map(spec -> spec !== nothing, touchdown_specs) : nothing
    owned_mask = _external_owned_mask(args, num_sats)
    owned_mask === nothing || (impact_excluded = impact_excluded === nothing ? owned_mask : impact_excluded .| owned_mask)
    impact_callback = impact_excluded === nothing ? get_impact_callback(num_sats) :
        get_impact_callback(num_sats; excluded_spacecraft=impact_excluded)
    callbacks = if backbone_mode
        (impact_callback,)
    else
        (
            impact_callback,
            update_planet_frame_callback(),
        )
    end
    # First among the discrete callbacks that follow: later ones read the synced shadow entries.
    owned_mask === nothing || (callbacks = _append_callback(callbacks, get_external_propagation_callback(num_sats)))

    has_touchdown && (callbacks = _append_callback(callbacks, get_touchdown_callback(touchdown_specs)))

    # Opt-in GRAM perturbed density (off by default -> nothing appended). Placed
    # before every callback that stages or reads density, so a factor updated at
    # an accepted step is the one those callbacks see at that step.
    if !backbone_mode
        callbacks = _append_callback(callbacks, get_gram_density_perturbation_callback(num_sats, args))
    end

    if !backbone_mode && _requires_staged_density_callback(effectors, args)
        callbacks = _append_callback(callbacks, get_density_callback(num_sats, effectors, args))
    end

    if !backbone_mode && _requires_thermal_callback(effectors, args)
        callbacks = _append_callback(callbacks, get_thermal_callback(num_sats, args))
    end

    if _requires_orbit_end_callback(args)
        callbacks = _append_callback(callbacks, get_orbit_end_callback(num_sats))
    end

    if !backbone_mode && _requires_entry_end_callback(effectors, args)
        callbacks = _append_callback(callbacks, get_entry_end_callback(num_sats, args))
    end

    if !backbone_mode && _requires_drag_state_callback(effectors, args)
        callbacks = _append_callback(callbacks, get_drag_state_callback(num_sats))
    end

    if !backbone_mode
        callbacks = _append_callbacks(callbacks, get_navigation_callbacks(num_sats, args))
        callbacks = _append_callbacks(callbacks, get_control_callbacks(num_sats, args))
        callbacks = _append_callbacks(callbacks, get_guidance_callbacks(num_sats, args))
    end
    if !backbone_mode && _requires_quaternion_projection_callback(args)
        callbacks = _append_callback(callbacks, get_quaternion_projection_callback(num_sats, args))
    end
    callbacks = _append_callback(callbacks, get_plume_callback(args))
    engine = _simulation_engine_module()
    output_solver_mode = args.solver_config === nothing ? engine._solver_policy_mode() :
        engine._solver_policy_mode(args.solver_config)
    callbacks = _append_callback(callbacks,
        get_initial_force_output_callback(effectors; solver_mode=output_solver_mode))
    if !backbone_mode && (args.simulation_settings.results || record_saved_values)
        callbacks = _append_callback(callbacks, get_data_saving_callback(num_sats, args, save_fields_resolved, saved_values))
    end
    callbacks = _append_callbacks(callbacks, extra_callbacks)

    return CallbackSet(callbacks...)
end
