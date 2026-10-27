using Test
using SpaceAGORA

# The RHS plan record carries the sweep time and a per-context sweep count, and
# a cache hit records neither.
const PP = SpaceAGORA.SimulationModel.ParallelPolicy

@testset "rhs plan telemetry: sweep time and count" begin
    PP.with_policy_context() do
        s0 = PP.policy_telemetry_snapshot()
        @test s0.rhs_plan_sweeps_total == 0
        @test s0.rhs_plan_sweep_ns == 0.0

        PP.record_rhs_plan_selection!(:cache, :flat, 4, :static)
        s1 = PP.policy_telemetry_snapshot()
        @test s1.rhs_plan_source == "cache"
        @test s1.rhs_plan_sweeps_total == 0
        @test s1.rhs_plan_sweep_ns == 0.0

        PP.record_rhs_plan_selection!(:sweep, :flat, 8, :dynamic; sweep_ns=2.5e9)
        s2 = PP.policy_telemetry_snapshot()
        @test s2.rhs_plan_source == "sweep"
        @test s2.rhs_plan_mode == "flat"
        @test s2.rhs_plan_allotment == 8
        @test s2.rhs_plan_scheduler == "dynamic"
        @test s2.rhs_plan_sweep_ns == 2.5e9
        @test s2.rhs_plan_sweeps_total == 1
    end
    # A fresh context starts clean: nothing leaks from the scope above.
    PP.with_policy_context() do
        @test PP.policy_telemetry_snapshot().rhs_plan_sweeps_total == 0
    end
end
