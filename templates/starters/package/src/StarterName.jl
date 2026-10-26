"""
    StarterName

An RPO reference planner for SpaceAGORA. `StraightLinePlanner` is an example to
replace with the package's own method. Every `AbstractRPOPlanner` belongs to
SpaceAGORA's RPO reference family: the family builds the requests, validates each
proposed reference and installs only a validated copy, so a planner never
validates or installs its own output.
"""
module StarterName

using LinearAlgebra: norm
using Random: AbstractRNG
using SpaceAGORA

export StraightLinePlanner

"""
    StraightLinePlanner(; speed_fraction=0.5, acceleration_fraction=0.5)

A reference along the straight line from the chaser's position to the goal, in
the target's RTN frame: it accelerates, cruises and brakes within the given
fractions of the request's speed and acceleration limits. It ignores obstacles,
so the family's validation refuses a line that passes too close to the station.
It draws nothing from the random stream.
"""
Base.@kwdef struct StraightLinePlanner <: AbstractRPOPlanner
    speed_fraction::Float64 = 0.5
    acceleration_fraction::Float64 = 0.5
end

# What the planner supports: truth states in the target's RTN frame, no retiming.
SpaceAGORA.planner_capabilities(::StraightLinePlanner) =
    RPOPlannerCapabilities(state_sources=(:truth,), frames=(:target_rtn,))

# Deterministic: the conformance checks then require the stream to stay untouched.
SpaceAGORA.consumes_rng(::StraightLinePlanner) = false

# Distance travelled and speed at time `t` of a profile that accelerates at `a` to
# `v`, cruises and brakes to rest at `distance` (a triangle when it never cruises).
function line_profile(t, v, a, ramp, duration, distance)
    t <= 0 && return (0.0, 0.0)
    t >= duration && return (distance, 0.0)
    t < ramp && return (a * t^2 / 2, a * t)
    t <= duration - ramp && return (a * ramp^2 / 2 + v * (t - ramp), v)
    left = duration - t
    return (distance - a * left^2 / 2, a * left)
end

# `rng` is the only random stream a planner may draw from; this one draws nothing.
# Type it `AbstractRNG`, as SpaceAGORA's own method for every planner does.
function SpaceAGORA.plan_rpo!(state, planner::StraightLinePlanner, request::RPOPlanningRequest, rng::AbstractRNG)
    refuse(termination) = RPOPlanningResult(request_id=request.request_id, status=:infeasible, termination=termination)
    start = collect(Float64, request.x_rtn[1:3])
    goal = collect(Float64, request.goal_rtn_m)
    distance = norm(goal - start)
    limits = request.constraints
    v = planner.speed_fraction * something(limits.max_speed_mps, Inf)
    a = planner.acceleration_fraction * something(limits.max_acceleration_mps2, Inf)
    isfinite(v) && isfinite(a) && v > 0 && a > 0 || return refuse(:needs_speed_and_acceleration_limits)
    if distance == 0
        ramp, duration = 0.0, 0.0
    elseif distance >= v^2 / a
        ramp = v / a
        duration = distance / v + ramp
    else
        v = sqrt(distance * a)
        ramp = v / a
        duration = 2ramp
    end
    dt = request.reference_dt_s
    steps = max(1, ceil(Int, duration / dt))
    request.time_s + steps * dt <= request.valid_until_s || return refuse(:longer_than_the_request_is_valid)
    direction = distance == 0 ? zeros(3) : (goal - start) / distance
    times = [(i - 1) * dt for i in 1:steps+1]
    samples = [line_profile(t, v, a, ramp, duration, distance) for t in times]
    positions = reduce(hcat, [start + s * direction for (s, _) in samples])
    velocities = reduce(hcat, [speed * direction for (_, speed) in samples])
    positions[:, end] = goal
    reference = RPOReference(t_ref_s=times, r_ref_rtn_m=positions, v_ref_rtn_mps=velocities,
        origin_time_s=request.time_s, valid_until_s=request.valid_until_s, chaser_id=request.chaser_id,
        target_id=request.target_id, geometry_revision=request.geometry_revision)
    return RPOPlanningResult(request_id=request.request_id, status=:candidate, termination=:completed,
        reference=reference, diagnostics=(duration_s=duration,))
end

end # module StarterName
