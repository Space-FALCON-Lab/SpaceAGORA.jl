using Test
using SpaceAGORA
const SEng = SpaceAGORA.SimulationEngine

# Drive _rhs_calib_cached_verdict through the in-memory store directly: the
# question is the gate's arithmetic, not the sweep.
function _with_entry(f, entry)
    mktempdir() do dir
        withenv("SPACEAGORA_RHS_CALIBRATION_PATH" => joinpath(dir, "store.toml"),
                "SPACEAGORA_RHS_CALIBRATE" => "auto",
                "SPACEAGORA_RHS_CALIBRATE_REVERIFY_SHARE" => "0.05") do
            lock(SEng._rhs_calib_lock) do
                SEng._rhs_calib_loaded[] = true
                empty!(SEng._rhs_calib_cache)
                SEng._rhs_calib_cache["sig"] = entry
            end
            try
                f()
            finally
                lock(SEng._rhs_calib_lock) do
                    empty!(SEng._rhs_calib_cache)
                    SEng._rhs_calib_loaded[] = false
                end
            end
        end
    end
end

# A CONFIRMED pin (two agreeing sweeps): the amortised budget is what this file
# tests; an unconfirmed pin's own rule is in rhs_reverify_unconfirmed_tests.
_plan_entry(; solve_ns, sweep_ns, honoured_ns) = Dict{String, Any}(
    "mode" => "flat_constellation_effector_queue", "allotment" => 4, "scheduler" => "static",
    "elapsed_mean_ns" => 1.0e6, "solve_ns" => solve_ns, "heuristic_votes" => 0,
    "sweep_ns" => sweep_ns, "honoured_ns" => honoured_ns, "plan_votes" => 2)
_heur_entry(; votes, solve_ns, sweep_ns, honoured_ns) = Dict{String, Any}(
    "mode" => SEng._CALIB_HEURISTIC_MODE, "allotment" => 1, "scheduler" => "auto",
    "elapsed_mean_ns" => 1.0e6, "solve_ns" => solve_ns, "heuristic_votes" => votes,
    "sweep_ns" => sweep_ns, "honoured_ns" => honoured_ns)

@testset "a cached plan on a long solve is honoured until the sweep has earned its share" begin
    long = 2.0e9
    # Fresh verdict, nothing run on it yet: honoured.
    _with_entry(_plan_entry(solve_ns = long, sweep_ns = 0.5e9, honoured_ns = 0.0)) do
        v = SEng._rhs_calib_cached_verdict("sig", true)
        @test v !== nothing && v !== :heuristic && v.allotment == 4
        @test !SEng._rhs_calib_reverify_due("sig")
    end
    # 20x the sweep's cost run on it (5 % share): one re-sweep is due.
    _with_entry(_plan_entry(solve_ns = long, sweep_ns = 0.5e9, honoured_ns = 10.0e9)) do
        @test SEng._rhs_calib_reverify_due("sig")
        @test SEng._rhs_calib_cached_verdict("sig", true) === nothing
    end
    # Sweep cost never measured (older store / width trial): re-verified, as shipped.
    _with_entry(_plan_entry(solve_ns = long, sweep_ns = 0.0, honoured_ns = 0.0)) do
        @test SEng._rhs_calib_cached_verdict("sig", true) === nothing
    end
    # Not V2: a long solve always re-sweeps, as shipped.
    _with_entry(_plan_entry(solve_ns = long, sweep_ns = 0.5e9, honoured_ns = 0.0)) do
        @test SEng._rhs_calib_cached_verdict("sig", false) === nothing
    end
    # Short solve: the cached plan is honoured either way.
    _with_entry(_plan_entry(solve_ns = 0.2e9, sweep_ns = 0.5e9, honoured_ns = 10.0e9)) do
        v = SEng._rhs_calib_cached_verdict("sig", false)
        @test v !== nothing && v !== :heuristic
    end
end

@testset "an unconfirmed heuristic verdict follows the same budget; a reproduced one is honoured outright" begin
    long = 2.0e9
    _with_entry(_heur_entry(votes = 1, solve_ns = long, sweep_ns = 0.5e9, honoured_ns = 0.0)) do
        @test SEng._rhs_calib_cached_verdict("sig", true) === :heuristic
    end
    _with_entry(_heur_entry(votes = 1, solve_ns = long, sweep_ns = 0.5e9, honoured_ns = 10.0e9)) do
        @test SEng._rhs_calib_cached_verdict("sig", true) === nothing
    end
    _with_entry(_heur_entry(votes = 3, solve_ns = long, sweep_ns = 0.5e9, honoured_ns = 10.0e9)) do
        @test SEng._rhs_calib_cached_verdict("sig", true) === :heuristic
    end
end

@testset "the store keeps the solve length across a re-pin and restarts the honoured clock" begin
    _with_entry(_plan_entry(solve_ns = 2.0e9, sweep_ns = 0.5e9, honoured_ns = 7.0e9)) do
        SEng._rhs_calib_store!("sig", SEng._make_calib_flat_plan(2, :static), 1.0e6; sweep_ns = 0.7e9)
        e = lock(SEng._rhs_calib_lock) do; SEng._rhs_calib_cache["sig"]; end
        @test e["solve_ns"] == 2.0e9
        @test e["sweep_ns"] == 0.7e9
        @test e["honoured_ns"] == 0.0
        SEng._rhs_calib_store_heuristic!("sig", 1.0e6; sweep_ns = 0.9e9)
        e = lock(SEng._rhs_calib_lock) do; SEng._rhs_calib_cache["sig"]; end
        @test e["solve_ns"] == 2.0e9 && e["sweep_ns"] == 0.9e9 && e["honoured_ns"] == 0.0
    end
end
