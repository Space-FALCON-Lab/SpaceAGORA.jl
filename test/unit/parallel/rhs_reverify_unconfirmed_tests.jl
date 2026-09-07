using Test
using SpaceAGORA
const SEng = SpaceAGORA.SimulationEngine

function _with_entry(f, entry)
    mktempdir() do dir
        withenv("SPACEAGORA_RHS_CALIBRATION_PATH" => joinpath(dir, "store.toml"),
                "SPACEAGORA_RHS_CALIBRATE" => "auto",
                "SPACEAGORA_RHS_CALIBRATE_REVERIFY_SHARE" => "0.05",
                "SPACEAGORA_RHS_CALIBRATE_REVERIFY_SHARE_UNCONFIRMED" => "0.20") do
            lock(SEng._rhs_calib_lock) do
                SEng._rhs_calib_loaded[] = true; empty!(SEng._rhs_calib_cache); SEng._rhs_calib_cache["sig"] = entry
            end
            try f() finally
                lock(SEng._rhs_calib_lock) do; empty!(SEng._rhs_calib_cache); SEng._rhs_calib_loaded[] = false; end
            end
        end
    end
end
_plan(; votes, honoured) = Dict{String, Any}("mode" => "flat_constellation_effector_queue", "allotment" => 4,
    "scheduler" => "static", "elapsed_mean_ns" => 1.0e6, "solve_ns" => 2.0e9, "heuristic_votes" => 0,
    "sweep_ns" => 0.5e9, "honoured_ns" => honoured, "plan_votes" => votes)

@testset "a plan pinned by one sweep is re-verified at 5x its cost; a confirmed one at 20x" begin
    _with_entry(_plan(votes = 1, honoured = 2.0e9)) do; @test !SEng._rhs_calib_reverify_due("sig"); end
    _with_entry(_plan(votes = 1, honoured = 2.5e9)) do; @test SEng._rhs_calib_reverify_due("sig"); end
    _with_entry(_plan(votes = 2, honoured = 2.5e9)) do; @test !SEng._rhs_calib_reverify_due("sig"); end
    _with_entry(_plan(votes = 2, honoured = 10.0e9)) do; @test SEng._rhs_calib_reverify_due("sig"); end
end

@testset "plan votes count consecutive sweeps that pinned the same plan" begin
    _with_entry(_plan(votes = 1, honoured = 0.0)) do
        SEng._rhs_calib_store!("sig", SEng._make_calib_flat_plan(4, :static), 1.0e6; sweep_ns = 0.5e9)
        @test lock(SEng._rhs_calib_lock) do; SEng._rhs_calib_cache["sig"]["plan_votes"]; end == 2
        SEng._rhs_calib_store!("sig", SEng._make_calib_flat_plan(8, :static), 1.0e6; sweep_ns = 0.5e9)
        @test lock(SEng._rhs_calib_lock) do; SEng._rhs_calib_cache["sig"]["plan_votes"]; end == 1
        SEng._rhs_calib_store!("sig", SEng._make_calib_satellite_batch_plan(), 1.0e6; sweep_ns = 0.5e9)
        @test lock(SEng._rhs_calib_lock) do; SEng._rhs_calib_cache["sig"]["plan_votes"]; end == 1
    end
end
