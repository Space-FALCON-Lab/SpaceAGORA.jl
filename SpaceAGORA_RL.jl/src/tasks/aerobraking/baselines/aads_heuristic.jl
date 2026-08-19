Base.@kwdef struct AADSHeuristicPolicy <: AbstractPolicy
    target_heat_rate_w_cm2::Float64 = 0.15
    bisection_xtol_mps::Float64 = 0.05
end

const _AADS_INITIAL_BISECTION_BOUND_MPS = 1.0
const _AADS_EXPANDED_BISECTION_BOUNDS_MPS = (6.0, 3.0)
const _AADS_PREDICTION_SEED_OFFSET = 1_000_003
const _AADS_MAX_GRAM_SEED = 2_147_483_646

function _aads_prediction_seed(actual_seed::Integer)
    return mod(Int(actual_seed) - 1 + _AADS_PREDICTION_SEED_OFFSET, _AADS_MAX_GRAM_SEED) + 1
end

function _aads_prediction_state(state::AerobrakingDecisionState)
    return AerobrakingDecisionState(
        pass_index=state.pass_index,
        apoapsis_radius_m=state.apoapsis_radius_m,
        periapsis_altitude_m=state.periapsis_altitude_m,
        inclination_rad=state.inclination_rad,
        raan_rad=state.raan_rad,
        argument_of_periapsis_rad=state.argument_of_periapsis_rad,
        true_anomaly_rad=state.true_anomaly_rad,
        epoch=state.epoch,
        mission_elapsed_s=state.mission_elapsed_s,
        total_delta_v_mps=state.total_delta_v_mps,
        maneuver_count=state.maneuver_count,
        previous_density_kg_m3=state.previous_density_kg_m3,
        previous_heat_rate_w_cm2=state.previous_heat_rate_w_cm2,
        last_drag_passage_time_s=state.last_drag_passage_time_s,
        previous_pass_minimum_altitude_m=state.previous_pass_minimum_altitude_m,
        gram_seed=_aads_prediction_seed(state.gram_seed),
        aerodynamic_cd_scale=state.aerodynamic_cd_scale,
        aerodynamic_cl_scale=state.aerodynamic_cl_scale,
    )
end

function predicted_heat_rate(config, state, action::AerobrakingAction)
    prediction_state = _aads_prediction_state(state)
    prediction_seed = prediction_state.gram_seed + state.pass_index
    prediction_rng = MersenneTwister(prediction_seed)
    next_state = propagate_pass(
        SpaceAGORACoreAdapter(config.backend_mode),
        config,
        prediction_state,
        action,
        prediction_rng,
    )
    return next_state.previous_heat_rate_w_cm2
end

predicted_heat_rate(config, state, action_index::Integer) =
    predicted_heat_rate(config, state, action_from_index(action_index))

function _aads_bisect_action(policy::AADSHeuristicPolicy, config, state,
                             direction::Float64, zero_heat::Float64)
    xtol = policy.bisection_xtol_mps
    isfinite(xtol) && xtol > 0.0 ||
        throw(ArgumentError("AADS bisection tolerance must be finite and positive"))

    heat_cache = Dict{Float64,Float64}(0.0 => zero_heat)
    function target_error(magnitude::Float64)
        heat = get!(heat_cache, magnitude) do
            action = continuous_action_from_delta_v(direction * magnitude)
            predicted_heat_rate(config, state, action)
        end
        isfinite(heat) || throw(ErrorException("AADS predictor returned a non-finite heat rate"))
        return heat - policy.target_heat_rate_w_cm2
    end

    f_zero = target_error(0.0)
    f_zero == 0.0 && return continuous_action_from_delta_v(0.0)
    bounds = (_AADS_INITIAL_BISECTION_BOUND_MPS, _AADS_EXPANDED_BISECTION_BOUNDS_MPS...)
    last_error = nothing
    for upper in bounds
        try
            f_hi = target_error(upper)
            f_hi == 0.0 && return continuous_action_from_delta_v(direction * upper)
            signbit(f_zero) == signbit(f_hi) && continue

            lo = 0.0
            hi = upper
            f_lo = f_zero
            while hi - lo > 2xtol
                mid = (lo + hi) / 2
                f_mid = target_error(mid)
                f_mid == 0.0 && return continuous_action_from_delta_v(direction * mid)
                if signbit(f_mid) == signbit(f_lo)
                    lo = mid
                    f_lo = f_mid
                else
                    hi = mid
                end
            end
            return continuous_action_from_delta_v(direction * ((lo + hi) / 2))
        catch err
            last_error = err
        end
    end

    detail = last_error === nothing ? "the target was not bracketed" : sprint(showerror, last_error)
    throw(ErrorException(
        "AADS could not bracket the $(policy.target_heat_rate_w_cm2) W/cm^2 heat-rate target " *
        "with maneuver magnitudes up to 6 m/s: $detail",
    ))
end

function policy_action_index(policy::AADSHeuristicPolicy, config, state,
                             observation::PaperObservation, rng::AbstractRNG)
    zero_action = continuous_action_from_delta_v(0.0)
    zero_heat = predicted_heat_rate(config, state, zero_action)
    if config.reward_config.heat_low_w_cm2 <= zero_heat <= config.reward_config.heat_high_w_cm2
        return zero_action
    end

    direction = zero_heat < config.reward_config.heat_low_w_cm2 ? -1.0 : 1.0
    return _aads_bisect_action(policy, config, state, direction, zero_heat)
end
