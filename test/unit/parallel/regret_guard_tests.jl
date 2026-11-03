using Test
using SpaceAGORA
const PP = SpaceAGORA.SimulationModel.ParallelPolicy

@testset "regret guard decision (pure)" begin
    d = PP.regret_guard_decision
    @test d(100.0, [130.0, 130.0, 130.0], 0.10, 3) == (true, :fire)
    @test d(100.0, [130.0, 130.0], 0.10, 3) == (false, :too_few_solves)
    # a cold first solve (JIT) with warm solves inside the margin does not fire
    @test d(100.0, [300.0, 105.0, 108.0], 0.10, 3) == (false, :within_margin)
    @test d(100.0, [111.0, 112.0, 115.0], 0.10, 3) == (true, :fire)
    @test d(100.0, [110.0, 112.0, 115.0], 0.10, 3) == (false, :within_margin)  # boundary: not strictly greater
    @test d(0.0, [1.0], 0.1, 1) == (false, :disabled_input)
    @test d(NaN, [1.0], 0.1, 1) == (false, :disabled_input)
    @test d(100.0, [NaN, Inf, 200.0], 0.1, 1) == (true, :fire)  # non-finite samples are dropped
end

@testset "regret guard env + observe" begin
    PP.reset_regret_guard!()
    withenv("SPACEAGORA_POLICY_REGRET_GUARD" => nothing) do
        @test PP.regret_guard_observe!("s", 1e9, 1.0) == false   # off: inert
        @test !PP.regret_guard_frozen()
    end
    withenv("SPACEAGORA_POLICY_REGRET_GUARD" => "1", "SPACEAGORA_POLICY_REGRET_MARGIN" => nothing) do
        @test_throws ArgumentError PP.regret_guard_observe!("s", 1.0, 1.0)
    end
    withenv("SPACEAGORA_POLICY_REGRET_GUARD" => "1", "SPACEAGORA_POLICY_REGRET_MARGIN" => "0.1",
            "SPACEAGORA_POLICY_REGRET_MIN_SOLVES" => "2") do
        PP.reset_policy_telemetry!()
        @test PP.regret_guard_observe!("s", 200.0, 100.0) == false
        @test PP.regret_guard_observe!("s", 200.0, 100.0) == true
        @test PP.regret_guard_frozen()
        snap = PP.regret_guard_snapshot()
        @test snap.fired && snap.reason == "fire" && snap.margin == 0.1 && snap.count == 2
        ts = PP.policy_telemetry_snapshot()    # the ppb CSV reads these
        @test ts.regret_guard_enabled && ts.regret_guard_fired && ts.regret_guard_frozen && ts.regret_guard_count == 2
        PP.reset_policy_telemetry!()           # freeze survives telemetry reset
        @test PP.policy_telemetry_snapshot().regret_guard_frozen && !PP.policy_telemetry_snapshot().regret_guard_fired
        @test PP.regret_guard_observe!("s", 1.0, 100.0) == true
    end
    PP.reset_regret_guard!()
end

@testset "regret guard hook arithmetic (pure)" begin
    SEng = SpaceAGORA.SimulationEngine
    b = SEng._rhs_calib_regret_baseline_ns
    @test b(10.0, 2.0, 4.0) == 5.0           # prev * ref / elapsed
    @test b(10.0, 4.0, 4.0) == 10.0
    @test b(0.0, 2.0, 4.0) == 0.0            # no earlier solve
    @test b(10.0, 0.0, 4.0) == 0.0           # no reference candidate
    @test b(10.0, 2.0, 0.0) == 0.0
    @test b(NaN, 2.0, 4.0) == 0.0
    @test b(10.0, Inf, 4.0) == 0.0
    base_row = Dict("mode" => "flat_constellation_effector_queue", "allotment" => 1, "scheduler" => "static")
    @test SEng._rhs_calib_is_baseline_row(base_row)
    @test !SEng._rhs_calib_is_baseline_row(merge(base_row, Dict("scheduler" => "dynamic")))
    @test SEng._rhs_calib_regret_ref_ns(Dict("candidates" => [Dict("mode" => "satellite_batch", "ns" => 1.0), merge(base_row, Dict("ns" => 7.0))])) == 7.0
    @test SEng._rhs_calib_regret_ref_ns(Dict()) == 0.0
end
