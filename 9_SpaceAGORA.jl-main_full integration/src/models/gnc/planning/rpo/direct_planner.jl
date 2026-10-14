"""Internal direct-path baseline. Loads with the neutral contract and standard libraries only."""
module DirectRPOPlanning
using LinearAlgebra: norm
using ..RPOPlannerInterfaces
const P = RPOPlannerInterfaces
import ..RPOPlannerInterfaces: planner_capabilities, plan_rpo!

"""
Direct segment with a cubic rest-to-rest time law. The assembly supplies a trusted
analytic `segment_clearance_at(start, goal, geometry)` query for its declared
geometry. Missing queries reject; this baseline does not find obstacle detours.
`max_reference_samples` bounds allocation before building the timed reference.
"""
struct DirectRPOPlanner{F} <: P.AbstractRPOPlanner
    segment_clearance_at::F
    headroom::P.RPOPlanningHeadroom
    max_reference_samples::Int
    function DirectRPOPlanner(; segment_clearance_at=nothing,
                              headroom=P.RPOPlanningHeadroom(), max_reference_samples::Integer=100_000)
        2 <= max_reference_samples < typemax(Int) || throw(ArgumentError("Invalid reference sample budget."))
        new{typeof(segment_clearance_at)}(segment_clearance_at, headroom, Int(max_reference_samples))
    end
end
planner_capabilities(::DirectRPOPlanner) = P.RPOPlannerCapabilities(state_sources=(:truth,), frames=(:target_rtn,))

function plan_rpo!(::Nothing, planner::DirectRPOPlanner, request::P.RPOPlanningRequest, rng::P.AbstractRNG)
    reject(status, reason; diagnostics=NamedTuple()) = P.RPOPlanningResult(
        request_id=request.request_id, status=status, termination=reason, diagnostics=diagnostics)
    request.reason === :retime && return reject(:unsupported, :retiming_not_implemented)
    planner.segment_clearance_at === nothing && return reject(:unsupported, :segment_clearance_not_supplied)
    request.constraints.max_speed_mps === nothing && return reject(:unsupported, :speed_limit_required)
    a, b = collect(request.x_rtn[1:3]), collect(request.goal_rtn_m)
    d = b - a; distance = norm(d)
    isfinite(distance) || return reject(:unsupported, :nonfinite_displacement)
    budget = P.rpo_planning_budget(request, planner.headroom;
        position_scale_m=max(norm(a), norm(b)), velocity_scale_mps=request.constraints.max_speed_mps)
    budget.supported || return reject(:unsupported, budget.reason; diagnostics=(planning_budget=budget,))
    clearance = Float64(planner.segment_clearance_at(Tuple(a), Tuple(b), request.geometry))
    isfinite(clearance) || return reject(:failed, :nonfinite_segment_clearance)
    clearance >= request.constraints.clearance_m || return reject(:infeasible, :direct_segment_blocked;
        diagnostics=(segment_clearance_m=clearance,))
    # For s(q)=3q^2-2q^3, max s'=3/2 and max |s''|=6 on [0,1].
    duration = max(request.reference_dt_s, 1.5 * (distance / budget.max_speed_mps))
    if budget.max_acceleration_mps2 !== nothing
        duration = max(duration, sqrt(6 * (distance / budget.max_acceleration_mps2)))
    end
    intervals = duration / request.reference_dt_s
    isfinite(intervals) && intervals <= planner.max_reference_samples - 1 ||
        return reject(:failed, :reference_work_limit)
    count = max(1, ceil(Int, intervals))
    duration = count * request.reference_dt_s
    isfinite(duration) && request.time_s + duration <= request.valid_until_s ||
        return reject(:infeasible, :insufficient_reference_lifetime)
    times = [i * request.reference_dt_s for i in 0:count]
    positions, velocities = zeros(3, count + 1), zeros(3, count + 1)
    for i in eachindex(times)
        q = times[i] / duration
        positions[:, i] .= a .+ (q*q*(3 - 2q)) .* d
        velocities[:, i] .= (6q*(1 - q) / duration) .* d
    end
    positions[:, 1] .= a; positions[:, end] .= b
    velocities[:, 1] .= 0; velocities[:, end] .= 0
    actual_budget = P.rpo_planning_budget(request, planner.headroom;
        position_scale_m=maximum(norm, eachcol(positions)), velocity_scale_mps=maximum(norm, eachcol(velocities)))
    actual_budget.supported || return reject(:unsupported, actual_budget.reason; diagnostics=(planning_budget=actual_budget,))
    ref = P.RPOReference(t_ref_s=times, r_ref_rtn_m=positions, v_ref_rtn_mps=velocities,
        origin_time_s=request.time_s, valid_until_s=request.valid_until_s,
        chaser_id=request.chaser_id, target_id=request.target_id, frame=request.frame,
        geometry_revision=request.geometry_revision)
    return P.RPOPlanningResult(request_id=request.request_id, status=:candidate, termination=:completed, reference=ref,
        diagnostics=(planner=:direct, planning_budget=actual_budget, duration_s=duration,
            analytic_max_speed_mps=1.5*(distance/duration),
            analytic_max_acceleration_mps2=6*(distance/duration)/duration,
            segment_clearance_m=clearance, path_kind=:polyline, path=hcat(a,b)))
end
end
