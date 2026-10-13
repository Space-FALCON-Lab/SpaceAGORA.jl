# Journal VI.A-D, evaluated on the executable command grid. This certificate
# covers the sampled family, not the journal's continuous-family theorem.
@inline function _edg_replan_due(t::Float64, last_t::Float64, interval::Float64)
    return !isfinite(last_t) || _edg_switch_due(t, last_t + interval)
end

function _edg_candidate_constraint_sample(config, p, spacecraft, r, v, t, alpha;
    heat_rate_control::Bool, structural_control::Bool)
    env = _edg_targeting_prediction_environment(p, r, v, t)
    rate = _edg_maxwellian_heat_rate(p, env, alpha)
    minimum_rate = _edg_maxwellian_heat_rate(p, env, config.min_alpha_rad)
    reference_area = _energy_depletion_struct_drag_area(spacecraft, env.temperature,
        env.molecular_speed_ratio, config.controlled_panel_links, config.max_alpha_rad, config)
    minimum_area = _energy_depletion_struct_drag_area(spacecraft, env.temperature,
        env.molecular_speed_ratio, config.controlled_panel_links, config.min_alpha_rad, config)
    minimum_load = env.dynamic_pressure * minimum_area / max(reference_area, eps(Float64))
    finite = all(isfinite, (rate, minimum_rate, minimum_load))
    admissible = finite && (!heat_rate_control || minimum_rate <= config.heat_rate_limit_w_cm2) &&
        (!structural_control || minimum_load <= config.structural_load_limit_pa)
    return (;rate, minimum_rate, minimum_load, admissible)
end

function _edg_targeting_candidate_certificate(config, p, spacecraft, mass, t, outcome;
    heat_load_j_cm2::Float64, heat_rate_control::Bool, structural_control::Bool)
    track = outcome.track
    heat = 0.0
    admissible = true
    max_minimum_rate = 0.0
    max_minimum_load = 0.0
    samples = 0
    # Simpson quadrature holds each command over its own interval, including the
    # fractional exit interval. It never interpolates across a command jump.
    for j in 1:length(track.time)-1
        tau = track.time[j]
        tau >= outcome.duration_s && break
        stop = min(track.time[j + 1], outcome.duration_s)
        dt = stop - tau
        r, v, alpha = track.positions[j], track.velocities[j], track.alpha_profile[j]
        accel(rr, vv, tt) = _edg_prediction_gravity_acceleration(p, rr, vv, mass, t + tt) +
            _edg_targeting_aero_acceleration(config, p, spacecraft, rr, vv, mass, t + tt, alpha)
        rm, vm = _edg_prediction_advance(accel, r, v, tau, dt / 2)
        re, ve = stop == outcome.duration_s ? (outcome.exit_position, outcome.exit_velocity) :
            (track.positions[j + 1], track.velocities[j + 1])
        a = _edg_candidate_constraint_sample(config, p, spacecraft, r, v, t + tau, alpha;
            heat_rate_control, structural_control)
        b = _edg_candidate_constraint_sample(config, p, spacecraft, rm, vm, t + tau + dt/2, alpha;
            heat_rate_control, structural_control)
        c = _edg_candidate_constraint_sample(config, p, spacecraft, re, ve, t + stop, alpha;
            heat_rate_control, structural_control)
        heat += dt * (a.rate + 4b.rate + c.rate) / 6
        admissible &= a.admissible && b.admissible && c.admissible
        max_minimum_rate = max(max_minimum_rate, a.minimum_rate, b.minimum_rate, c.minimum_rate)
        max_minimum_load = max(max_minimum_load, a.minimum_load, b.minimum_load, c.minimum_load)
        samples += 3
    end
    total_heat = heat_load_j_cm2 + heat
    heat_feasible = isfinite(total_heat) &&
        (!(:heat_load in config.max_energy_submodes) || total_heat <= config.heat_load_limit_j_cm2)
    return (heat_load_j_cm2=total_heat, instantaneous_admissible=admissible,
        heat_feasible=heat_feasible, max_minimum_heat_rate_w_cm2=max_minimum_rate,
        max_minimum_structural_load_pa=max_minimum_load, samples=samples)
end

function _edg_ordered_feasible_prefix(switches, evaluate, reference_energy::Float64)
    accepted = NamedTuple[]
    status = :horizon_end
    rejected = nothing
    for switch in switches
        candidate = evaluate(switch)
        certificate = candidate.certificate
        if !(certificate.instantaneous_admissible && certificate.heat_feasible)
            status = certificate.instantaneous_admissible ? :heat_boundary : :instantaneous_infeasibility
            rejected = candidate
            break
        elseif !isempty(accepted) && !(candidate.energy_jkg < last(accepted).energy_jkg)
            status = :ordering_not_resolved
            rejected = candidate
            break
        elseif !isempty(accepted) && candidate.energy_jkg < reference_energy
            status = :depletion_reference_boundary
            rejected = candidate
            break
        end
        push!(accepted, candidate)
        if candidate.energy_jkg == reference_energy
            status = :depletion_reference_reached
            break
        end
    end
    return (;accepted, status, rejected)
end

# With the atmospheric perturbations already propagated, Eq. (32) needs only
# the candidate-specific coast correction. Adding passage work again is wrong.
function _edg_candidate_energy_mapping(config, p, spacecraft, mass, t, candidate)
    coast = _edg_vacuum_apoapsis_correction(p, candidate.exit_position,
        candidate.exit_velocity, mass, t + candidate.duration_s)
    coast.event_reached || error("EDG candidate coast did not reach apoapsis.")
    target = _edg_corrected_target_energy_from_apoapsis(p.args.environment_model.planet,
        config.target_apoapsis_radius_m, coast.periapsis_radius_m, coast.energy_change_jkg)
    residual = candidate.energy_jkg - target
    isfinite(residual) || error("EDG candidate apoapsis mapping has a nonfinite residual.")
    return (residual_jkg=residual, target_energy_jkg=target,
        coast_energy_change_jkg=coast.energy_change_jkg,
        apoapsis_radius_m=norm(coast.position),
        apoapsis_error_m=norm(coast.position)-config.target_apoapsis_radius_m,
        periapsis_radius_m=coast.periapsis_radius_m)
end

function _edg_select_mapped_candidate(candidates, mapping)
    isempty(candidates) && throw(ArgumentError("No admissible one-switch candidates."))
    cache = Dict{Int, NamedTuple}()
    mapped(k) = get!(() -> mapping(candidates[k]), cache, k)
    left, right = 1, length(candidates)
    fl, fr = mapped(left).residual_jkg, mapped(right).residual_jkg
    # Energy decreases with later switches. A negative upper residual means the
    # target requires *less* depletion than even minimum drag can provide.
    if fl < 0.0 && fr < 0.0
        return (index=left, mode=:safe_low_drag, status=:above_upper_energy,
            left=left, right=right, left_residual_jkg=fl, right_residual_jkg=fr, mapping=mapped(left))
    elseif fl > 0.0 && fr > 0.0
        return (index=right, mode=:max_energy_depletion, status=:below_lower_energy,
            left=left, right=right, left_residual_jkg=fl, right_residual_jkg=fr, mapping=mapped(right))
    elseif fl < fr
        return (index=left, mode=:safe_low_drag, status=:mapping_ordering_not_resolved,
            left=left, right=right, left_residual_jkg=fl, right_residual_jkg=fr, mapping=mapped(left))
    end
    while right-left > 1 && fl != 0.0 && fr != 0.0
        middle = (left+right) ÷ 2
        fm = mapped(middle).residual_jkg
        if !(fr <= fm <= fl)
            return (index=left, mode=:safe_low_drag, status=:mapping_ordering_not_resolved,
                left=left, right=right, left_residual_jkg=fl, right_residual_jkg=fr, mapping=mapped(left))
        elseif fm > 0.0
            left, fl = middle, fm
        else
            right, fr = middle, fm
        end
    end
    chosen = abs(fl) <= abs(fr) ? left : right
    return (index=chosen, mode=:targeting,
        status=mapped(chosen).residual_jkg == 0.0 ? :exact_sample : :quantized_bracket,
        left=left, right=right, left_residual_jkg=fl, right_residual_jkg=fr, mapping=mapped(chosen))
end

function _edg_journal_targeting_plan(config, p, spacecraft, pos, vel, mass_state, t;
    heat_load_j_cm2::Float64=0.0)
    mass = _edg_predict_mass(spacecraft, mass_state)
    heat_rate_control = :heat_rate in config.max_energy_submodes
    structural_control = :structural_load in config.max_energy_submodes
    function evaluate(switch)
        outcome = _edg_predict_targeting_outcome(config, p, spacecraft, pos, vel, mass, t, switch;
            heat_rate_control, structural_control)
        certificate = _edg_targeting_candidate_certificate(config, p, spacecraft, mass, t, outcome;
            heat_load_j_cm2, heat_rate_control, structural_control)
        # Keep compact candidate data, not hundreds of duplicate trajectory arrays.
        return (switch_s=switch, energy_jkg=outcome.energy_jkg,
            exit_position=outcome.exit_position, exit_velocity=outcome.exit_velocity,
            duration_s=outcome.duration_s, certificate=certificate)
    end
    clock = _edg_targeting_control_clock(config, p)
    first_switch = _edg_at_control_tick(t, clock.period, clock.origin) ? t :
        clock.origin + (floor(Int, (t-clock.origin)/clock.period)+1)*clock.period
    low = evaluate(first_switch)
    if !(low.certificate.instantaneous_admissible && low.certificate.heat_feasible)
        # Preserve the depletion fallback when no EDG-T family is certified.
        # The upper-energy direction still requires minimum drag (journal VI.A).
        # Neither fallback is labeled constraint-feasible by this failed certificate.
        mapping = _edg_candidate_energy_mapping(config, p, spacecraft, mass, t, low)
        mode = mapping.residual_jkg < 0.0 || !(:max_energy_depletion in config.guidance_modes) ?
            :safe_low_drag : :max_energy_depletion
        return (mode=mode, switch_s=Inf, status=:no_admissible_family,
            family_status=low.certificate.instantaneous_admissible ? :heat_boundary : :instantaneous_infeasibility,
            target_energy_jkg=mapping.target_energy_jkg, coast_energy_change_jkg=mapping.coast_energy_change_jkg,
            lower_energy_jkg=NaN, upper_energy_jkg=low.energy_jkg,
            candidate_count=1, certificate=low.certificate,
            minimum_drag_apoapsis_radius_m=mapping.apoapsis_radius_m,
            minimum_drag_residual_jkg=mapping.residual_jkg)
    end
    env = _edg_targeting_prediction_environment(p, pos, vel, t)
    reference = _edg_predict_max_energy_depletion_outcome(config, p, spacecraft, pos, vel, mass, t, env,
        heat_load_j_cm2; heat_rate_control, structural_control)
    switches = first_switch .+ clock.period .* (0:ceil(Int, low.duration_s/clock.period))
    prefix = _edg_ordered_feasible_prefix(switches,
        sw -> sw == first_switch ? low : evaluate(sw), reference.energy_jkg)
    candidates = prefix.accepted
    isempty(candidates) && error("EDG minimum-drag candidate lost its certificate.")
    selected = _edg_select_mapped_candidate(candidates,
        c -> _edg_candidate_energy_mapping(config, p, spacecraft, mass, t, c))
    candidate = candidates[selected.index]
    mode = selected.mode == :max_energy_depletion && !(:max_energy_depletion in config.guidance_modes) ?
        :safe_low_drag : selected.mode
    return (mode=mode, switch_s=mode == :targeting ? candidate.switch_s : Inf,
        status=selected.status, family_status=prefix.status,
        target_energy_jkg=selected.mapping.target_energy_jkg,
        coast_energy_change_jkg=selected.mapping.coast_energy_change_jkg,
        lower_energy_jkg=last(candidates).energy_jkg, upper_energy_jkg=low.energy_jkg,
        candidate_count=length(candidates), certificate=candidate.certificate,
        residual_jkg=selected.mapping.residual_jkg,
        predicted_apoapsis_radius_m=selected.mapping.apoapsis_radius_m,
        predicted_apoapsis_error_m=selected.mapping.apoapsis_error_m,
        left_switch_s=candidates[selected.left].switch_s,
        right_switch_s=candidates[selected.right].switch_s,
        left_residual_jkg=selected.left_residual_jkg,
        right_residual_jkg=selected.right_residual_jkg)
end
