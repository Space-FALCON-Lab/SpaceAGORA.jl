# The RPO examples on the planner contract: the main example plans inside the run,
# including retained-path replay and deterministic disturbance behavior.
module RPOExamplesPlannerContractTests
using Test
using LinearAlgebra
using Random
using StaticArrays
using JSON
using SpaceAGORA
const SM = SpaceAGORA.SimulationModel
const P = SpaceAGORA.RPOPlannerInterfaces

# The example scripts carry their own includes and constants; load them into a
# dedicated module so nothing leaks into Main or other test files. The
# batch example includes the main example.
module RPOExamplesOnContract
using SpaceAGORA
end
const RPO_EXAMPLES_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "examples"))
const RPO_EXAMPLES_SPICE = joinpath(RPO_EXAMPLES_DIR, "..", "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
module RPODisturbancesExample
using SpaceAGORA
end

# The disturbance example on the planner path: the obstacle's refused replans hold
# the chaser until the obstacle clears, and the displacement replans at once.
@testset "Disturbance example: obstacle, hold and displacement" begin
    Base.include(RPODisturbancesExample, joinpath(RPO_EXAMPLES_DIR, "Earth_RPO_Planner_Disturbances.jl"))
    out = Base.invokelatest(RPODisturbancesExample.run_rpo_disturbance_demo)
    timeline = out.timeline
    @test [x.event for x in timeline if x.event !== :request] == [:installed,
        :trigger, :rejected, :fallback_installed,   # the sphere blocks the path: hold
        :trigger, :rejected, :fallback_installed,   # settled on the hold, still blocked
        :trigger, :installed,                       # the sphere is gone: replanned
        :installed]                                 # displaced: tracking-error replan
    requests = [x.request for x in timeline if x.event === :request]
    @test [r.request_id for r in requests] == 1:5
    @test [r.reason for r in requests] == [:initial, :replan, :replan, :replan, :replan]
    @test [r.geometry_revision for r in requests] ==
        ["rpo-pilot-v1", "rpo-pilot-v1+spheres:1", "rpo-pilot-v1+spheres:1", "rpo-pilot-v1", "rpo-pilot-v1"]
    times = [x.time_s for x in timeline if x.event === :request]
    # Appearance at 10 s; the sphere stays active through 30 s, so its removal is seen
    # at the next tick; the displacement's priority callback precedes the 55 s tick.
    @test isapprox(times[[1, 2, 4, 5]], [0.0, 10.0, 30.1, 55.0]; atol=1e-9)
    # Settling depends on the braking flight; recorded at 24.2 s, pinned to one tick.
    @test isapprox(times[3], 24.2; atol=0.1 + 1e-9)
    marks = [x for x in timeline if x.event === :trigger]
    @test [x.trigger for x in marks] == [:geometry_revision, :hold_settled, :geometry_revision]
    @test marks[1].details.action === :replan && marks[1].details.min_clearance_m < 0.2
    @test marks[3].details.action === :replan && isnan(marks[3].details.min_clearance_m)  # a hold is replaced unchecked
    @test all(x -> x.termination === :direct_segment_blocked, filter(x -> x.event === :rejected, timeline))
    holds = [x.product for x in timeline if x.event === :fallback_installed]
    @test all(h -> h.r_ref_rtn_m[:, 1] == h.r_ref_rtn_m[:, end] && iszero(h.v_ref_rtn_mps), holds)
    installed = filter(x -> x.event === :installed, timeline)
    @test [x.request_id for x in installed] == [1, 4, 5]
    @test all(x -> x.validation.accepted && x.validation.metrics.min_clearance_m >= 0.2, installed)
    # The last reference starts at the displaced chaser, 0.3 m normal to the plane.
    @test isapprox(last(installed).reference.r_ref_rtn_m[3, 1], 0.3; atol=1e-3)
    @test last(installed).reference.r_ref_rtn_m[:, end] == [6.0, 0.0, 0.0]
    s = out.summary
    @test s.end_time_s == 110.0 && string(s.retcode) == "Success"
    @test s.final_goal_error_m < 1e-3 && s.final_speed_mps < 1e-3
    @test 0.004 < s.propellant_used_kg < 0.008   # recorded 5.84 g
    @test s.min_saved_clearance_m >= 0.2
end

@testset "Examples build the planner contract" begin
    source(name) = read(joinpath(RPO_EXAMPLES_DIR, name), String)
    for name in filter(endswith(".jl"), readdir(RPO_EXAMPLES_DIR))
        @test !occursin("RPOGuidanceModel", source(name))
    end
    @test occursin("PlannerGuidance", source("Earth_RPO_CubeSat_MPC.jl"))
    @test occursin("RetainedPathRPOPlanner", source("Earth_RPO_CubeSat_MPC.jl"))
    @test occursin("make_rpo_configuration", source("Earth_RPO_Planner_Disturbances.jl"))
end

if !isdir(RPO_EXAMPLES_SPICE)
    @info "SPICE kernels absent; skipping the RPO example planner-contract runs" RPO_EXAMPLES_SPICE
else
    Base.include(RPOExamplesOnContract, joinpath(RPO_EXAMPLES_DIR, "Earth_RPO_CubeSat_MPC_Batch.jl"))
    const X = RPOExamplesOnContract

    @testset "The main example plans, validates, flies and stops on the planner path" begin
        small = (pso_n_particles=8, pso_n_iters=2, n_station_points=256, verbose=false, results_directory=mktempdir(),
            mission_time=10.0, data_rate_s=1.0, start_rtn=SVector(-8.0, -4.0, 2.0), goal_rtn=SVector(-5.5, -2.5, 1.0),
            record_control_commands=true)
        demo = Base.invokelatest(X.build_rpo_cubesat_mpc_demo; small...)
        @test demo.guidance isa SpaceAGORA.RPOPlannerLifecycle.PlannerGuidance
        @test demo.planner isa HYPRRPOPlanner
        # The guidance plans with it under the cap check (see below).
        capped = demo.guidance.planner
        @test capped isa X.RPOMissionCapPlanner && capped.planner === demo.planner && capped.end_by_s == 1800.0
        @test planner_capabilities(capped) == planner_capabilities(demo.planner)
        @test consumes_rng(capped) && !consumes_rng(X.RPOMissionCapPlanner(RetainedPathRPOPlanner(zeros(3, 2), RPOPSOConfig()), 1.0))
        @test demo.guidance.constraints.clearance_m == 0.1
        @test demo.guidance.constraints.max_speed_mps === nothing
        @test demo.guidance.constraints.max_acceleration_mps2 === nothing   # the default retimer is not acceleration-limited
        @test demo.args.mission_configuration.mission_time == demo.mission.max_time_s == 1800.0
        @test !haskey(demo, :initial_plan)
        sol = Base.invokelatest(X.run_rpo_case, demo)
        plan = Base.invokelatest(X.rpo_installed_plan, sol)
        report = only(rpo_run_report(sol).planners)
        @test report.request_id == 1 && plan.validation.accepted && plan.planner === :hypr
        @test report.stream_words == gnc_stream_words("RPO/v1", 741, 101, 201)
        @test plan.r_ref_rtn[:, 1] ≈ [-8.0, -4.0, 2.0] atol=1e-9
        @test plan.r_ref_rtn[:, end] == [-5.5, -2.5, 1.0]
        # The run ends post_reference_hold_s after the installed reference, as the
        # builder's simulated time did before.
        @test last(sol.t) ≈ max(10.0, plan.t_ref_s[end] + 20.0) atol=0.06
        @test !isempty(only(rpo_run_report(sol).commands).log.t_s)
        # Its path and configuration replay through RetainedPathRPOPlanner to the same
        # reference and the same flight.
        replay = Base.invokelatest(X.build_rpo_cubesat_mpc_demo; small..., pso_config=plan.config,
            precomputed_plan=(path=plan.path_rtn, config=plan.config))
        @test replay.planner isa RetainedPathRPOPlanner
        sol2 = Base.invokelatest(X.run_rpo_case, replay)
        plan2 = Base.invokelatest(X.rpo_installed_plan, sol2)
        @test plan2.planner === :retained_path
        @test (plan2.t_ref_s, plan2.r_ref_rtn, plan2.v_ref_rtn) == (plan.t_ref_s, plan.r_ref_rtn, plan.v_ref_rtn)
        @test sol2.t == sol.t && sol2.u == sol.u
        # Planning-time hooks of the old builder are refused, not ignored.
        @test_throws ArgumentError Base.invokelatest(X.build_rpo_cubesat_mpc_demo; small..., pso_iteration_callback=identity)
        @test_throws ArgumentError Base.invokelatest(X.build_rpo_cubesat_mpc_demo; small..., max_mission_time_s=5.0)
        # The reference must end by max_mission_time_s. The plan's lifetime adds the
        # controller's preview, so the cap is checked separately: replayed, the same
        # reference is refused under a cap just before its end and accepted at its end.
        ends_s = plan.t_ref_s[end]
        replay_capped(cap) = Base.invokelatest(X.build_rpo_cubesat_mpc_demo; small..., pso_config=plan.config,
            precomputed_plan=(path=plan.path_rtn, config=plan.config), max_mission_time_s=cap)
        short = replay_capped(ends_s - 0.5)
        @test ends_s - 0.5 >= 10.0 && ends_s - 0.5 + (short.mpc.horizon + 1) * short.mpc.dt_s > ends_s  # the lifetime alone admits it
        refused = try run_simulation(short.args; return_solution=true); nothing catch e; e end
        @test refused isa RPOPlanningError && refused.reason === :infeasible && refused.request_id == 1
        failure = last(refused.report.records)
        @test failure.event === :failure && failure.details.termination === :reference_ends_after_mission_cap
        @test !any(x -> x.event === :installed, refused.report.records)
        first_step = X.DiscreteCallback((u, t, integrator) -> t > 0, X.terminate!; save_positions=(false, false))
        at_cap = run_simulation(replay_capped(ends_s).args; return_solution=true, extra_callbacks=(first_step,))
        @test Base.invokelatest(X.rpo_installed_plan, at_cap).t_ref_s == plan.t_ref_s
    end
end
end # module RPOExamplesPlannerContractTests
