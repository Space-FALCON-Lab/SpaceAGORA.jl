# RPO, the first family: its representative requests, the trusted queries of its
# validation, its conformance link and a short simulation through the RPO binding.

const _RPO = RPOPlannerInterfaces

"""
    rpo_conformance_requests(; goal_rtn_m=(5.0, 0.0, 0.0)) -> Vector{RPOPlanningRequest}

Representative requests of the RPO reference family for
[`check_algorithm`](@ref): an initial request at 0 s from 3 m along the target's
radial axis, a forced request at 1 s from another state, and a retiming request at
2 s, with the pilot's chaser and station geometry (the station 50 m off the path)
and constraints. `goal_rtn_m=nothing` makes each request's goal its own start
position, for a planner valid only for an unchanged goal.
"""
function rpo_conformance_requests(; goal_rtn_m::Union{Nothing, Tuple{Real, Real, Real}}=(5.0, 0.0, 0.0))
    geometry = SimulationModel.RPOReferenceGeometry(
        SimulationModel.RPOStationGeometry(reshape([0.0, 0.0, 50.0], 3, 1); keepout_radius_m=0.25);
        chaser=SimulationModel.RPOCubeSatGeometry(dims_m=(0.1, 0.1, 0.3)))
    constraints = _RPO.RPOPlanningConstraints(clearance_m=0.2, max_speed_mps=0.5, max_acceleration_mps2=0.1)
    starts = ((3.0, 0.0, 0.0, 0.0, 0.0, 0.0), (4.0, 0.5, 0.0, 0.0, 0.0, 0.0), (4.0, 0.5, 0.0, 0.0, 0.0, 0.0))
    reasons = (:initial, :forced, :retime)
    return [_RPO.RPOPlanningRequest(request_id=k, reason=reasons[k], chaser_id=101, target_id=201,
        epoch="conformance", time_s=k - 1.0, x_rtn=starts[k], target_state_ii=(7.0e6, 0.0, 0.0, 0.0, 7500.0, 0.0),
        goal_rtn_m=goal_rtn_m === nothing ? starts[k][1:3] : Float64.(goal_rtn_m), geometry=geometry,
        geometry_revision="conformance-v1", constraints=constraints, reference_dt_s=0.1, preview_horizon_steps=2,
        valid_until_s=k + 99.0) for k in 1:3]
end

"""
    rpo_trusted_queries() -> NamedTuple

The queries the RPO family's validation needs, as the RPO binding supplies them:
`clearance_at`, the distance from a point to the station in a request's geometry.
"""
rpo_trusted_queries() =
    (clearance_at=(point, geometry) -> SimulationModel.rpo_clearance_distance_to_station(point, geometry),)

conformance_link(::_RPO.RPOReferenceFamily) = SimulationModel.RPOPlanBuffer()

"""
    rpo_conformance_simulation(; mission_time_s=0.2, kwargs...) -> Function

A `simulation` for [`check_algorithm`](@ref): a function of a planner that runs
`make_rpo_configuration(; planner, mission_time_s, kwargs...)` twice in sequence
and, with two or more threads, twice concurrently, and returns the planner
binding's run records of the first run. It throws when a run fails or the runs
differ, so `simulation_run` fails for a planner whose runs do not repeat.
"""
rpo_conformance_simulation(; mission_time_s::Real=0.2, kwargs...) =
    planner -> _rpo_simulation_records(planner; mission_time_s, kwargs...)

function _rpo_simulation_records(planner; kwargs...)
    configuration = RPOPlannerLifecycle.make_rpo_configuration(; planner, kwargs...)
    base = _solve(configuration)
    _failed(base) && throw(base.error)
    again = _solve(configuration)
    gap = _run_difference(base, again)
    gap === nothing || error("two runs of the RPO configuration differ: $(gap)")
    concurrency = _concurrency_check(configuration, base, again, FIRST_RELEASE)
    concurrency.status === :fail && error(concurrency.detail)
    return only(RPOPlannerLifecycle.rpo_run_report(base.solution).planners).records
end
