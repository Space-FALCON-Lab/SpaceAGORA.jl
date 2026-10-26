# A chaser on the RPO planner path that touches down mid-run: what its controller
# and planner do once the spacecraft is inactive. The right-hand side freezes an
# inactive spacecraft, and both models stop ticking for it (`skip_inactive`), so
# these runs pin where the command log and the planner's records end, that a
# replan event after the touchdown is not delivered, and the trajectory.
module InactiveChaserTests
using Test, LinearAlgebra, StaticArrays, SpaceAGORA
const SM = SpaceAGORA.SimulationModel
const ZERO3 = SVector(0.0, 0.0, 0.0)

# A control model that only declares a touchdown surface for one spacecraft: a
# sphere of radius `radius_m` about the planet's centre.
struct TouchdownSurface <: SM.AbstractControlEffectorModel
    spacecraft_index::Int
    radius_m::Float64
    touchdown_s::Vector{Float64}
end
SpaceAGORA.spacecraft_scope(s::TouchdownSurface) = BoundSpacecraft((s.spacecraft_index,))
SpaceAGORA.touchdown_spec(s::TouchdownSurface, i::Int) = i == s.spacecraft_index ?
    (terrain=NoTerrainModel(), reference_radius_m=s.radius_m, height_m=0.0,
     on_touchdown=(t, r_p, v_p, k) -> push!(s.touchdown_s, t)) : nothing
SpaceAGORA.update_control!(::TouchdownSurface, u, p, t, i) = nothing
SpaceAGORA.control_wrench(::TouchdownSurface, x::StateSample, p, t) = (ZERO3, ZERO3)

# The planner pilot with the chaser 3 m above the target, sent 2 m inward, over a
# surface 4 mm below its start: it touches down at about 1.31 s of the 2 s run.
function touchdown_args(; kwargs...)
    args = make_rpo_configuration(; planner=DirectRPOPlanner(), start_rtn_m=(3.0, 0.0, 0.0),
        goal_rtn_m=(1.0, 0.0, 0.0), kwargs...)
    start = args.dynamics_model.spacecraft[1].initial_condition.pos
    control = args.control_model
    return SM.SimConfig._with_configuration(args; control_model=SM.ControlModel(
        control_effectors=(control.control_effectors..., TouchdownSurface(1, norm(start) - 0.004, Float64[])),
        control_rates=[control.control_rates..., control.control_rates[1]]))
end
fly(; kwargs...) = run_simulation(touchdown_args(; kwargs...); return_solution=true)

# A replan requested after the touchdown.
const REPLAN = (planning_events=(RPOPlanningEvent(1.6; reason=:replan, token=1),),)

run_model(sol, T) = only(m for m in (sol.prob.p.args.guidance_model.guidance_effectors...,
    sol.prob.p.args.control_model.control_effectors...) if m isa T)
touchdown_time(sol) = only(run_model(sol, TouchdownSurface).touchdown_s)
record_times(records, event) = [r.time_s for r in records if r.event === event]

@testset "A chaser deactivated mid-run on the planner path" begin
    sol = fly(; REPLAN..., on_rejection=:install_fallback)
    t_down = touchdown_time(sol)
    @test 1.3 < t_down < 1.4
    @test sol.prob.p.is_active == [false, true] && sol.t[end] == 2.0
    # The right-hand side freezes the chaser from the touchdown on; the target flies on.
    frozen = sol.u[end].sc[1]
    after = findall(>(t_down), sol.t)
    @test length(after) > 1
    @test all(k -> sol.u[k].sc[1].pos == frozen.pos && sol.u[k].sc[1].vel == frozen.vel, after)
    @test sol.u[after[1]].sc[2].pos != sol.u[end].sc[2].pos

    report = rpo_run_report(sol)
    log, records = report.commands[1].log, report.planners[1].records
    # The controller and the planner stop at the touchdown (`skip_inactive`): the
    # command log ends with the last tick before it, and the held command is the
    # last one logged.
    @test length(log.t_s) == 13 && all(<(t_down), log.t_s)
    @test record_times(records, :control) == record_times(records, :guidance) == log.t_s
    @test held_command(run_model(sol, SM.RPOMPCControlModel), 1).thruster_forces_n == log.thruster_forces_n[end]
    # Nothing is recorded after the touchdown: the event is never delivered, and no
    # request is built from the frozen state.
    @test !any(r -> hasproperty(r, :time_s) && r.time_s > t_down, records)
    @test report.planners[1].request_id == 1
    @test !(1 in run_model(sol, SpaceAGORA.RPOPlannerLifecycle.PlannerGuidance).runtime.processed_tokens)

    # The trajectory is the one without the replan.
    plain = fly()
    @test touchdown_time(plain) == t_down && plain.t == sol.t && plain.u == sol.u
    # With the default `on_rejection=:fail` the run no longer stops: the same flight,
    # entry for entry and record for record.
    failing = fly(; REPLAN...)
    failing_report = rpo_run_report(failing)
    failing_log = failing_report.commands[1].log
    @test failing.t == sol.t && failing.u == sol.u
    @test all(name -> getfield(failing_log, name) == getfield(log, name), fieldnames(typeof(log)))
    @test isequal([(r.event, get(r, :time_s, NaN)) for r in failing_report.planners[1].records],
                  [(r.event, get(r, :time_s, NaN)) for r in records])
end
end
