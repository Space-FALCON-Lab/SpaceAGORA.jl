# Worked example C of "Writing a GNC Algorithm" (docs/src/user/gnc_planning.md):
# a planner of the RPO reference family, run in the RPO pilot. Each region between
# a `# snippet:` line and `# end snippet` appears on that page as it is here; the
# HPC and extensibility docs gate compares them.
using SpaceAGORA
using SpaceAGORA.Conformance: test_algorithm, rpo_conformance_requests, rpo_trusted_queries,
    rpo_conformance_simulation
using LinearAlgebra: norm
using Random: AbstractRNG
using Test

# snippet: planner
# A reference along the straight line to the goal, in the target's RTN frame: it
# accelerates, cruises and brakes within fractions of the request's limits. It
# ignores obstacles; the family's validation refuses a line that comes too close.
Base.@kwdef struct LinePlanner <: AbstractRPOPlanner
    speed_fraction::Float64 = 0.5            # of the request's speed limit
    acceleration_fraction::Float64 = 0.5     # of the request's acceleration limit
end

SpaceAGORA.planner_capabilities(::LinePlanner) =
    RPOPlannerCapabilities(state_sources=(:truth,), frames=(:target_rtn,))
SpaceAGORA.consumes_rng(::LinePlanner) = false     # it draws nothing from `rng`

function SpaceAGORA.plan_rpo!(state, planner::LinePlanner, request::RPOPlanningRequest, rng::AbstractRNG)
    refuse(why) = RPOPlanningResult(request_id=request.request_id, status=:infeasible, termination=why)
    start, goal = collect(Float64, request.x_rtn[1:3]), collect(Float64, request.goal_rtn_m)
    distance = norm(goal - start)
    v = planner.speed_fraction * something(request.constraints.max_speed_mps, Inf)
    a = planner.acceleration_fraction * something(request.constraints.max_acceleration_mps2, Inf)
    isfinite(v) && isfinite(a) && v > 0 && a > 0 || return refuse(:needs_speed_and_acceleration_limits)
    v = min(v, sqrt(distance * a))           # too short to reach cruise speed: a triangle
    ramp_s = v / a
    duration_s = distance == 0 ? 0.0 : distance / v + ramp_s
    dt = request.reference_dt_s
    steps = max(1, ceil(Int, duration_s / dt))
    request.time_s + steps * dt <= request.valid_until_s || return refuse(:outlives_the_request)
    # Distance along the line and speed at time t from the start.
    function along(t)
        t >= duration_s && return (distance, 0.0)
        t < ramp_s && return (a * t^2 / 2, a * t)
        t <= duration_s - ramp_s && return (v * (t - ramp_s / 2), v)
        left = duration_s - t
        return (distance - a * left^2 / 2, a * left)
    end
    direction = distance == 0 ? zeros(3) : (goal - start) / distance
    times = [k * dt for k in 0:steps]
    profile = along.(times)
    reference = RPOReference(t_ref_s=times,
        r_ref_rtn_m=reduce(hcat, [start + s * direction for (s, _) in profile]),
        v_ref_rtn_mps=reduce(hcat, [speed * direction for (_, speed) in profile]),
        origin_time_s=request.time_s, valid_until_s=request.valid_until_s,
        chaser_id=request.chaser_id, target_id=request.target_id,
        geometry_revision=request.geometry_revision)
    return RPOPlanningResult(request_id=request.request_id, status=:candidate,
        termination=:completed, reference=reference, diagnostics=(duration_s=duration_s,))
end
# end snippet

# snippet: run
# The RPO pilot: a chaser 3 m from the station plans to 5 m, and an MPC controller
# tracks the installed reference.
configuration = make_rpo_configuration(planner=LinePlanner(), mission_time_s=2.0)
solution = run_simulation(configuration; return_solution=true)
report = rpo_run_report(solution)
installed = only(filter(record -> record.event === :installed, only(report.planners).records))
println("installed request ", installed.request_id, ": ", installed.validation.reason,
        ", reference of ", last(installed.reference.t_ref_s), " s")
@test installed.validation.accepted
# end snippet

# snippet: conformance
# The family's representative requests, its trusted validation queries and a short
# simulation through the RPO binding.
@testset "LinePlanner conformance" begin
    test_algorithm(LinePlanner(), rpo_conformance_requests(); trusted=rpo_trusted_queries(),
                   simulation=rpo_conformance_simulation())
end
# end snippet
