# Disturbances on the RPO planner contract.
#
# A chaser flies a direct reference to a goal 3 m further along the radial axis of a
# passive target (`DirectRPOPlanner`). The planner lifecycle then handles:
# - an obstacle: a sphere appears on the remaining path at 10 s. The replan is
#   refused, because the direct segment is blocked, so the run installs a hold at
#   the chaser's position instead of stopping (`on_rejection=:install_fallback`).
#   Once the chaser settles on the hold, the planner is asked once more and refuses
#   again. When the sphere disappears at 30 s, the new geometry revision replaces
#   the hold with a new reference;
# - a displacement: at 55 s the chaser is moved 0.3 m normal to the orbit plane. The
#   tracking error exceeds 0.2 m, so the same guidance tick replans from there.
# `rpo_run_report` records every request, refusal, hold and installed reference.
# Runs in the examples environment (`--project=examples`); neither HYPR nor SPICE is needed.
include(joinpath(@__DIR__, "common.jl"))
using LinearAlgebra
using StaticArrays
using DiffEqCallbacks: PresetTimeCallback
using DiffEqBase: derivative_discontinuity!

const RPO_DISTURBANCE_FRAMES = SimulationModel.FrameTransforms
const RPO_DISTURBANCE_CLEARANCE = SimulationModel.GuidanceHooks.rpo_clearance_distance_to_station

"""
    build_rpo_disturbance_demo(; kwargs...)

Return the run configuration (`args`), the priority callback factories that apply
the displacement, and the scenario parameters. Every keyword has a default:
`start_rtn_m`, `goal_rtn_m`, `obstacle_center_rtn_m`, `obstacle_radius_m`,
`obstacle_appear_s`, `obstacle_disappear_s`, `displacement_time_s`,
`displacement_rtn_m`, `tracking_error_limit_m` and `mission_time_s`.
"""
function build_rpo_disturbance_demo(;
    start_rtn_m=(3.0, 0.0, 0.0),
    goal_rtn_m=(6.0, 0.0, 0.0),
    obstacle_center_rtn_m=(5.0, 0.0, 0.0),
    obstacle_radius_m=0.1,
    obstacle_appear_s=10.0,
    obstacle_disappear_s=30.0,
    displacement_time_s=55.0,
    displacement_rtn_m=(0.0, 0.0, 0.3),
    tracking_error_limit_m=0.2,
    mission_time_s=110.0,
)
    obstacles = RPOReplanningSpheres([(center_rtn=obstacle_center_rtn_m, radius_m=obstacle_radius_m,
        appear_time_s=obstacle_appear_s, disappear_time_s=obstacle_disappear_s, label=:debris)])
    args = make_rpo_configuration(
        planner=DirectRPOPlanner(),
        start_rtn_m=start_rtn_m,
        goal_rtn_m=goal_rtn_m,
        mission_time_s=mission_time_s,
        geometry_provider=obstacles,
        on_rejection=:install_fallback,
        tracking_error_limit_m=tracking_error_limit_m,
    )
    shift = SVector{3, Float64}(displacement_rtn_m)
    # The displacement keeps the chaser's velocity in the target's rotating frame.
    function displace!(integrator)
        chaser, target = integrator.u.sc[1], integrator.u.sc[2]
        rt, vt = SVector{3}(target.pos), SVector{3}(target.vel)
        x = RPO_DISTURBANCE_FRAMES.inertial_to_rtn_relative_state(SVector{3}(chaser.pos), SVector{3}(chaser.vel), rt, vt)
        r, v = RPO_DISTURBANCE_FRAMES.rtn_to_inertial_relative_state(SVector{3}(x[1:3]) + shift, SVector{3}(x[4:6]), rt, vt)
        chaser.pos .= r
        chaser.vel .= v
        derivative_discontinuity!(integrator, true)
        return nothing
    end
    # A priority callback runs before guidance, so the guidance tick at the
    # displacement time already sees the displaced chaser.
    factories = (_ -> PresetTimeCallback([Float64(displacement_time_s)], displace!; save_positions=(false, false)),)
    return (
        args=args,
        priority_callback_factories=factories,
        scenario=(start_rtn_m=start_rtn_m, goal_rtn_m=goal_rtn_m, obstacle_center_rtn_m=obstacle_center_rtn_m,
            obstacle_radius_m=obstacle_radius_m, obstacle_appear_s=obstacle_appear_s,
            obstacle_disappear_s=obstacle_disappear_s, displacement_time_s=displacement_time_s,
            displacement_rtn_m=displacement_rtn_m, tracking_error_limit_m=tracking_error_limit_m,
            mission_time_s=mission_time_s),
    )
end

"""
    rpo_disturbance_timeline(report)

The planner lifecycle's records without the per-tick guidance and control entries:
requests, triggers, refusals, holds (`:fallback_installed`) and installed references.
"""
function rpo_disturbance_timeline(report)
    records = only(report.planners).records
    return [x for x in records if x.event in (:request, :trigger, :rejected, :fallback_installed, :installed, :failure)]
end

"""
    rpo_disturbance_summary(solution, demo)

The end time, final distance to the goal, propellant used, and the smallest clearance
over the saved states to the station and to the obstacle while it is present.
"""
function rpo_disturbance_summary(solution, demo)
    binding = only(solution.prob.p.args.guidance_model.guidance_effectors)
    relative(u) = RPO_DISTURBANCE_FRAMES.inertial_to_rtn_relative_state(SVector{3}(u.sc[1].pos), SVector{3}(u.sc[1].vel),
        SVector{3}(u.sc[2].pos), SVector{3}(u.sc[2].vel))
    clearance = minimum(RPO_DISTURBANCE_CLEARANCE(relative(u)[1:3],
        rpo_geometry_snapshot(binding.geometry_provider, binding.geometry, t)) for (t, u) in zip(solution.t, solution.u))
    final = relative(last(solution.u))
    return (
        end_time_s=last(solution.t),
        retcode=solution.retcode,
        final_goal_error_m=norm(final[1:3] - collect(demo.scenario.goal_rtn_m)),
        final_speed_mps=norm(final[4:6]),
        propellant_used_kg=solution(0.0; continuity=:right).sc[1].mass - last(solution.u).sc[1].mass,
        min_saved_clearance_m=clearance,
    )
end

"""
    run_rpo_disturbance_demo(; kwargs...)

Build and run the scenario (`build_rpo_disturbance_demo` keywords). Returns the
solution, its `rpo_run_report`, the timeline and the summary.
"""
function run_rpo_disturbance_demo(; kwargs...)
    demo = build_rpo_disturbance_demo(; kwargs...)
    solution = run_simulation(demo.args; return_solution=true,
        priority_callback_factories=demo.priority_callback_factories)
    report = rpo_run_report(solution)
    return (demo=demo, solution=solution, report=report, timeline=rpo_disturbance_timeline(report),
        summary=rpo_disturbance_summary(solution, demo))
end

function _rpo_disturbance_detail(x)
    x.event === :request && return string(x.request.reason, " on ", x.request.geometry_revision)
    x.event === :trigger && return string(x.trigger, haskey(x.details, :action) ? string(" -> ", x.details.action) : "")
    x.event === :rejected && return string(x.termination)
    x.event === :installed && return string(x.reason, ", reference ", round(last(x.reference.t_ref_s); digits=1), " s")
    return ""
end

if abspath(PROGRAM_FILE) == @__FILE__
    out = run_rpo_disturbance_demo()
    for x in out.timeline
        println(rpad(string(round(x.time_s; digits=1), " s"), 9), rpad(string(x.event), 20), _rpo_disturbance_detail(x))
    end
    s = out.summary
    println("end ", s.end_time_s, " s (", s.retcode, "); final distance to goal ", s.final_goal_error_m,
        " m; propellant ", 1e3 * s.propellant_used_kg, " g; smallest saved clearance ", s.min_saved_clearance_m, " m")
end
