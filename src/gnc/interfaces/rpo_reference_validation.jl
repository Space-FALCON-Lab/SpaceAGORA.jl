"""
    validate_rpo_result(request, result; clearance_at=nothing, time_s=request.time_s)

Validate a candidate without mutating either input or installing a plan.
`clearance_at(point_rtn_m, geometry)` is a trusted geometry query supplied by the
caller, NOT the planner. It returns finite signed clearance after the declared
station/chaser approximation. Missing queries fail closed. Unexpected query
exceptions propagate, rather than being mislabeled as infeasibility.

Clearance is sampled on the reference polyline, independently of any diagnostic
control polygon. Speed is sampled velocity norm; acceleration is a finite
difference of velocity. These are reference checks, not physical tracking proof.
"""
function validate_rpo_result(request::RPOPlanningRequest, result::RPOPlanningResult;
                             clearance_at=nothing, time_s=request.time_s)
    result.request_id == request.request_id || return _reject(:request_mismatch)
    result.status === :candidate || return _reject(result.status)
    result.termination in (:completed, :iteration_limit, :time_budget) ||
        return _reject(:unusable_termination)
    result.termination === :time_budget && !request.validation.allow_time_budget_candidate &&
        return _reject(:time_budget_not_allowed)
    ref = result.reference::RPOReference
    ref.chaser_id == request.chaser_id && ref.target_id == request.target_id ||
        return _reject(:spacecraft_mismatch)
    ref.frame === request.frame || return _reject(:frame_mismatch)
    ref.geometry_revision == request.geometry_revision || return _reject(:geometry_mismatch)
    ref.terminal_policy === :repeat_last_sample || return _reject(:unsupported_terminal_policy)
    times, positions, velocities = ref.t_ref_s, ref.r_ref_rtn_m, ref.v_ref_rtn_mps
    n = length(times)
    n >= 2 && size(positions) == (3, n) && size(velocities) == (3, n) ||
        return _reject(:reference_shape)
    all(isfinite, times) && all(isfinite, positions) && all(isfinite, velocities) &&
        isfinite(ref.origin_time_s) && isfinite(ref.valid_until_s) || return _reject(:nonfinite_reference)
    settings = request.validation
    tol = settings.time_atol_s
    abs(ref.origin_time_s - request.time_s) <= tol || return _reject(:origin_mismatch)
    ref.valid_until_s > ref.origin_time_s && ref.valid_until_s <= request.valid_until_s + tol ||
        return _reject(:invalid_validity)
    times[1] == 0.0 || return _reject(:time_origin)
    for i in 2:n
        expected = (i - 1) * request.reference_dt_s
        isfinite(expected) && times[i] > times[i - 1] && abs(times[i] - expected) <= tol ||
            return _reject(:nonuniform_time_grid)
    end
    ref.origin_time_s + times[end] <= ref.valid_until_s + tol || return _reject(:reference_exceeds_validity)
    rpo_reference_is_current(request, ref, time_s) || return _reject(:expired_or_not_yet_valid)
    start_error = norm(positions[:, 1] .- request.x_rtn[1:3])
    goal_error = norm(positions[:, end] .- request.goal_rtn_m)
    max(start_error, goal_error) <= settings.endpoint_atol_m ||
        return _reject(:endpoint_mismatch; start_error_m=start_error, goal_error_m=goal_error)
    speed = maximum(norm, eachcol(velocities))
    acceleration = maximum(norm((velocities[:, i] - velocities[:, i-1]) / request.reference_dt_s) for i in 2:n)
    isfinite(speed) && isfinite(acceleration) || return _reject(:nonfinite_reference_metrics)
    constraints = request.constraints
    constraints.max_speed_mps !== nothing && speed > constraints.max_speed_mps &&
        return _reject(:speed_limit; max_speed_mps=speed)
    constraints.max_acceleration_mps2 !== nothing && acceleration > constraints.max_acceleration_mps2 &&
        return _reject(:acceleration_limit; max_acceleration_mps2=acceleration)
    clearance_at === nothing && return _reject(:clearance_not_checked)
    # Bound the complete validation work before calling the query. Never accept
    # a prefix of a reference when its remaining checks exceed the work cap.
    segments = Int[]
    sample_count = 1
    for i in 2:n
        count = norm(positions[:, i] - positions[:, i-1]) / settings.clearance_sample_ds_m
        isfinite(count) && count <= settings.max_clearance_samples - sample_count ||
            return _reject(:clearance_work_limit)
        steps = max(1, ceil(Int, count))
        sample_count += steps
        sample_count <= settings.max_clearance_samples || return _reject(:clearance_work_limit)
        push!(segments, steps)
    end
    first_point = ntuple(k -> positions[k, 1], 3)
    minimum_clearance = Float64(clearance_at(first_point, request.geometry))
    isfinite(minimum_clearance) || return _reject(:nonfinite_clearance)
    minimum_clearance >= constraints.clearance_m || return _reject(:clearance_limit; min_clearance_m=minimum_clearance)
    for i in 2:n
        for j in 1:segments[i-1]
            alpha = j / segments[i-1]
            p = ntuple(k -> (1 - alpha) * positions[k, i-1] + alpha * positions[k, i], 3)
            c = Float64(clearance_at(p, request.geometry))
            isfinite(c) || return _reject(:nonfinite_clearance)
            minimum_clearance = min(minimum_clearance, c)
            c >= constraints.clearance_m || return _reject(:clearance_limit; min_clearance_m=c)
        end
    end
    return RPOValidationResult(true, :validated,
        (start_error_m=start_error, goal_error_m=goal_error, max_speed_mps=speed,
         max_acceleration_mps2=acceleration, min_clearance_m=minimum_clearance,
         clearance_samples=sample_count, clearance_sample_ds_m=settings.clearance_sample_ds_m,
         geometry_revision=request.geometry_revision, clearance_method=:sampled_reference_polyline))
end
