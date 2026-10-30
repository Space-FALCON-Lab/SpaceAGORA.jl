module PeriapsisCorridorController

export select_corridor_maneuver

const _DEFAULT_LOWER_M = 95_000.0
const _DEFAULT_UPPER_M = 100_000.0

function _result(dv, altitude, reason, feasible, evaluations, target, requested_target; error=nothing)
    return (; delta_v_mps=Float64(dv), predicted_periapsis_altitude_m=Float64(altitude),
        reason, feasible, evaluations, target_periapsis_altitude_m=Float64(target),
        requested_target_periapsis_altitude_m=Float64(requested_target), error)
end

"""
    select_corridor_maneuver(predictor; lower_m=95_000, upper_m=100_000,
        margin_m=250, target_m=97_500, policy=:near_passive,
        max_dv_mps=20, deadband_m=250)

Choose a signed tangential apoapsis-reset impulse using an injected predictor.
`predictor(dv::Float64)` must return the next physical periapsis GEODETIC
altitude (height above the planet ellipsoid) in meters after that hypothetical
reset; `lower_m`, `upper_m`, `target_m` and `margin_m` are in the same
altitude, so the policy itself is agnostic to how it is computed.
`monotone_tolerance_m` (default 0) widens the nondecreasing-response check by the
predictor's own noise: finite-burn execution makes the predicted periapsis
jitter by a few tenths of a meter at burns of ~1e-3 m/s. The policy checks only the sampled local
response on `[-max_dv_mps, max_dv_mps]`; it does not claim a global horizon
optimum. A nondecreasing sampled response is required, consistent with a local
prograde-increases-periapsis model.
"""
function select_corridor_maneuver(predictor;
        lower_m::Real=_DEFAULT_LOWER_M,
        upper_m::Real=_DEFAULT_UPPER_M,
        margin_m::Real=250.0,
        target_m::Real=(_DEFAULT_LOWER_M + _DEFAULT_UPPER_M) / 2,
        policy::Symbol=:near_passive,
        max_dv_mps::Real=20.0,
        deadband_m::Real=250.0,
        monotone_tolerance_m::Real=0.0)
    lo, hi, margin = Float64(lower_m), Float64(upper_m), Float64(margin_m)
    target, max_dv, deadband = Float64(target_m), Float64(max_dv_mps), Float64(deadband_m)
    all(isfinite, (lo, hi, margin, target, max_dv, deadband)) ||
        throw(ArgumentError("corridor policy parameters must be finite"))
    lo < hi || throw(ArgumentError("lower_m must be below upper_m"))
    margin >= 0.0 && 2.0 * margin < hi - lo ||
        throw(ArgumentError("margin_m must leave a nonempty tightened corridor"))
    safe_lo, safe_hi = lo + margin, hi - margin
    safe_lo <= target <= safe_hi ||
        throw(ArgumentError("target_m must lie inside the tightened physical corridor"))
    max_dv >= 0.0 || throw(ArgumentError("max_dv_mps must be nonnegative"))
    monotone_tolerance_m >= 0.0 || throw(ArgumentError("monotone_tolerance_m must be nonnegative"))
    deadband >= 0.0 || throw(ArgumentError("deadband_m must be nonnegative"))
    policy in (:near_passive, :frequent) ||
        throw(ArgumentError("policy must be :near_passive or :frequent"))

    values = Dict{Float64, Float64}()
    evaluation_count = 0
    function evaluate(dv::Float64)
        haskey(values, dv) && return values[dv]
        evaluation_count += 1
        altitude = try
            Float64(predictor(dv))
        catch err
            throw(ErrorException("predictor threw $(typeof(err)): $(sprint(showerror, err))"))
        end
        isfinite(altitude) || throw(ErrorException("predictor returned a non-finite altitude"))
        values[dv] = altitude
        return altitude
    end
    function explored_side_is_monotone()
        samples = sort!(collect(values); by=first)
        for idx in 2:length(samples)
            previous, current = samples[idx - 1][2], samples[idx][2]
            tolerance = 1e-8 * max(1.0, abs(previous), abs(current)) + monotone_tolerance_m
            current + tolerance >= previous || return false
        end
        return true
    end
    failed_result(reason, error="") = _result(0.0, NaN, reason, false,
        evaluation_count, target, target; error)

    y0 = try
        evaluate(0.0)
    catch err
        return failed_result(:predictor_error, sprint(showerror, err))
    end

    if policy === :near_passive && safe_lo <= y0 <= safe_hi
        return _result(0.0, y0, :inside_safety_band, true, evaluation_count, target, target)
    elseif policy === :frequent && safe_lo <= y0 <= safe_hi && abs(y0 - target) <= deadband
        return _result(0.0, y0, :inside_deadband, true, evaluation_count, target, target)
    end

    direction = target > y0 ? 1.0 : -1.0
    step = max_dv / 1024.0
    previous_dv = 0.0
    end_dv, end_altitude = 0.0, y0
    bracket = nothing
    try
        while step > 0.0
            dv = direction * min(step, max_dv)
            altitude = evaluate(dv)
            explored_side_is_monotone() ||
                return failed_result(:nonmonotonic_predictor,
                    "sampled response on the selected burn side is not monotonic")
            end_dv, end_altitude = dv, altitude
            if direction * (altitude - target) >= 0.0
                bracket = (previous_dv, dv)
                break
            end
            previous_dv = dv
            step >= max_dv && break
            step = min(2.0 * step, max_dv)
        end
    catch err
        return failed_result(:predictor_error, sprint(showerror, err))
    end

    if bracket === nothing
        inside_safe = safe_lo <= end_altitude <= safe_hi
        reason = inside_safe ? :target_limited_to_reachable : :actuator_limit_unsafe
        return _result(end_dv, end_altitude, reason, inside_safe,
            evaluation_count, end_altitude, target)
    end

    dv_lo, dv_hi = minmax(bracket...)
    altitude_lo, altitude_hi = evaluate(dv_lo), evaluate(dv_hi)
    try
        for _ in 1:80
            mid = (dv_lo + dv_hi) / 2.0
            altitude_mid = evaluate(mid)
            explored_side_is_monotone() ||
                return failed_result(:nonmonotonic_predictor,
                    "sampled response on the selected burn side is not monotonic")
            if altitude_mid < target
                dv_lo, altitude_lo = mid, altitude_mid
            else
                dv_hi, altitude_hi = mid, altitude_mid
            end
            abs(altitude_mid - target) <= 0.01 && break
            dv_hi - dv_lo <= 1e-10 && break
        end
    catch err
        return failed_result(:predictor_error, sprint(showerror, err))
    end

    dv = abs(altitude_lo - target) <= abs(altitude_hi - target) ? dv_lo : dv_hi
    altitude = evaluate(dv)
    return _result(dv, altitude, :target_reached, safe_lo <= altitude <= safe_hi,
        evaluation_count, target, target)
end

end # module PeriapsisCorridorController
