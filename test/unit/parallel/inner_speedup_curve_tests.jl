using Test
using SpaceAGORA
using TOML

# The RHS calibration store's per-candidate timings (schema 2), its code token,
# the inner-speedup curve assembled from those timings, and the predictive
# planner's plans whose samples run on more than one thread. No solves: the
# store is driven through its own functions and a temporary path.

const SEng = SpaceAGORA.SimulationEngine
const SCamp = SpaceAGORA.SimulationCampaigns

const TOKEN = SEng._RHS_CALIB_CODE_TOKEN

_sig(; budget = 8, outer = 0, eff = "abcd", code = true) =
    "v6|machine=testbox|budget=$(budget)|sats=5_8|effs=2|harm=0|eff=$(eff)|dens=ExponentialAtmosphereModel|outer=$(outer)" *
    (code ? "|code=$(TOKEN)" : "")

_stem(; eff = "abcd") = join(sort!([
    "machine=testbox", "sats=5_8", "effs=2", "harm=0", "eff=$(eff)",
    "dens=ExponentialAtmosphereModel", "code=$(TOKEN)"]), "|")

_flat(a) = SEng._make_calib_flat_plan(a, :static)

_record(plan, ns; reps = 4, round = 2, default = false) =
    SEng._rhs_sweep_timing_record(plan, Float64(ns), reps, round, default)

# Run `f` against a fresh store at `path`: the in-process mirror is dropped
# before and after, so nothing leaks into or out of the real store.
function _with_store(f, path)
    withenv("SPACEAGORA_RHS_CALIBRATION_PATH" => path) do
        lock(SEng._rhs_calib_lock) do
            empty!(SEng._rhs_calib_cache)
            SEng._rhs_calib_loaded[] = false
            SEng._rhs_calib_loaded_path[] = ""
        end
        try
            f()
        finally
            lock(SEng._rhs_calib_lock) do
                empty!(SEng._rhs_calib_cache)
                SEng._rhs_calib_loaded[] = false
                SEng._rhs_calib_loaded_path[] = ""
            end
        end
    end
end

function _reload!()
    lock(SEng._rhs_calib_lock) do
        empty!(SEng._rhs_calib_cache)
        SEng._rhs_calib_loaded[] = false
        SEng._rhs_calib_loaded_path[] = ""
    end
    SEng._rhs_calib_load!()
end

@testset "the store keeps per-candidate timings and still reads old rows" begin
    dir = mktempdir()
    path = joinpath(dir, "rhs_calibration_test.toml")
    # A schema-1 file as older code wrote it: no timings, no code token.
    open(path, "w") do io
        TOML.print(io, Dict{String, Any}("schema_version" => 1, "calibrations" => [Dict{String, Any}(
            "signature" => _sig(code = false), "mode" => "flat_constellation_effector_queue",
            "allotment" => 4, "scheduler" => "static", "elapsed_mean_ns" => 2.0e6,
            "solve_ns" => 5.0e9, "heuristic_votes" => 0, "sweep_ns" => 1.0e8,
            "honoured_ns" => 0.0, "plan_votes" => 1)]))
    end
    _with_store(path) do
        old = SEng._rhs_calib_lookup(_sig(code = false))
        @test old !== nothing && old !== :heuristic
        @test old.allotment == 4

        timings = [_record(_flat(1), 4.0e6), _record(_flat(2), 2.2e6),
                   _record(_flat(4), 1.3e6, default = true), _record(_flat(4), 1.25e6; reps = 10, round = 0)]
        SEng._rhs_calib_store!(_sig(), _flat(4), 1.25e6; sweep_ns = 1.0e8, timings = timings)
        SEng._rhs_calib_save!()

        parsed = TOML.parsefile(path)
        @test parsed["schema_version"] == 2
        rows = Dict(r["signature"] => r for r in parsed["calibrations"])
        @test haskey(rows, _sig(code = false))
        @test !haskey(rows[_sig(code = false)], "candidates")
        @test rows[_sig(code = false)]["allotment"] == 4
        new = rows[_sig()]
        @test new["mode"] == "flat_constellation_effector_queue"
        @test new["allotment"] == 4 && new["elapsed_mean_ns"] == 1.25e6
        @test length(new["candidates"]) == 4
        c1 = new["candidates"][1]
        @test c1["mode"] == "flat_constellation_effector_queue" && c1["allotment"] == 1
        @test c1["width"] == 1 && c1["ns"] == 4.0e6 && c1["reps"] == 4 && c1["round"] == 2
        @test count(c -> c["default"], new["candidates"]) == 1
        @test new["candidates"][4]["round"] == 0

        _reload!()
        entry = SEng._rhs_calib_cache[_sig()]
        @test length(entry["candidates"]) == 4
        @test [c["ns"] for c in entry["candidates"]] == [4.0e6, 2.2e6, 1.3e6, 1.25e6]
        @test !haskey(SEng._rhs_calib_cache[_sig(code = false)], "candidates")

        # A verdict formed without a sweep keeps the shape's timings.
        SEng._rhs_calib_store_heuristic!(_sig(), 1.0e6)
        @test length(SEng._rhs_calib_cache[_sig()]["candidates"]) == 4
    end
end

@testset "the code token makes a store written by older code cold" begin
    @test TOKEN isa String && !isempty(TOKEN)
    dir = mktempdir()
    path = joinpath(dir, "rhs_calibration_old.toml")
    _with_store(path) do
        # The same shape, as a store written before the token existed.
        SEng._rhs_calib_store!(_sig(code = false), _flat(8), 1.0e6;
                               timings = [_record(_flat(1), 8.0e6), _record(_flat(8), 1.0e6)])
        @test SEng._rhs_calib_lookup(_sig()) === nothing
        @test SEng.rhs_inner_speedup_curve(_stem()) === nothing
        # And the same row under the current token is found.
        SEng._rhs_calib_store!(_sig(), _flat(8), 1.0e6;
                               timings = [_record(_flat(1), 8.0e6), _record(_flat(8), 1.0e6)])
        @test SEng._rhs_calib_lookup(_sig()) !== nothing
        @test SEng.rhs_inner_speedup_curve(_stem()) !== nothing
    end
    # The planner's corrections are keyed by the same token.
    @test SCamp.campaign_corrections_code_token() == TOKEN
end

@testset "the inner-speedup curve is assembled across a shape's rows" begin
    dir = mktempdir()
    _with_store(joinpath(dir, "rhs_calibration_curve.toml")) do
        SEng._rhs_calib_store!(_sig(budget = 8, outer = 0), _flat(4), 1.0e6;
                               timings = [_record(_flat(1), 8.0e6), _record(_flat(2), 4.4e6),
                                          _record(_flat(4), 2.5e6), _record(_flat(8), 3.0e6)])
        # Another budget of the same shape, under an outer split: its width-2
        # reading is the faster one and wins at width 2.
        SEng._rhs_calib_store!(_sig(budget = 2, outer = 1), _flat(2), 1.0e6;
                               timings = [_record(_flat(1), 8.2e6), _record(_flat(2), 4.0e6)])
        # A different shape never contributes.
        SEng._rhs_calib_store!(_sig(eff = "ffff"), _flat(8), 1.0e6;
                               timings = [_record(_flat(1), 8.0e6), _record(_flat(8), 0.5e6)])
        curve = SEng.rhs_inner_speedup_curve(_stem())
        @test length(curve) == 8
        @test curve[1] == 1.0
        @test curve[2] ≈ 8.0 / 4.0
        @test curve[3] ≈ 8.0 / 4.0          # nothing measured at 3: the best at <= 3
        @test curve[4] ≈ 8.0 / 2.5
        @test curve[8] ≈ 8.0 / 2.5          # width 8 is slower; b = 8 can run width 4
        # No serial reading: no curve.
        SEng._rhs_calib_store!(_sig(eff = "noserial"), _flat(8), 1.0e6;
                               timings = [_record(_flat(2), 4.0e6), _record(_flat(8), 1.0e6)])
        @test SEng.rhs_inner_speedup_curve(_stem(eff = "noserial")) === nothing
        @test SEng.rhs_inner_speedup_curve(_stem(eff = "absent")) === nothing
    end
    # The curve type never falls and starts at one.
    c = SCamp.InnerSpeedupCurve([0.9, 1.8, 1.5, 3.0])
    @test c.speedup == [1.0, 1.8, 1.8, 3.0]
    @test SCamp.inner_speedup(c, 64) == 3.0
    @test_throws ArgumentError SCamp.InnerSpeedupCurve(Float64[])
end

_cfg(; threads) = SCamp.PredictivePlannerConfig(
    margin = 0.15, guard_factor = 1.5, local_slots_max = threads - 1, remote_overhead = 0.0,
    route_switch = false, heap_model = :locals, local_thrash = 3.0)

_planning(; n, threads, pool = 0, curve = nothing) = SCamp.predictive_plan(
    n_samples = n, threads = threads, process_workers = pool, threads_candidate = true,
    local_slots_cap = threads - 1, constants = nothing, config = _cfg(threads = threads),
    inner_curve = curve)

@testset "a strong curve spends the idle threads; a flat one does not" begin
    strong = SCamp.InnerSpeedupCurve(Float64.(1:8))          # linear to eight threads
    flat = SCamp.InnerSpeedupCurve(ones(32))
    for pool in (0, 4)
        p = _planning(n = 4, threads = 32, pool = pool, curve = strong)
        @test p.chosen.route === :threads
        @test p.chosen.workers == 4
        @test p.chosen.inner_thread_budget == 8
        @test p.reason === :predicted_gain
        @test SCamp._predictive_declared_budget(p.chosen) == 8
        @test SCamp.predictive_plan_key(p.chosen) == "threads@w4+l0+b8"

        q = _planning(n = 4, threads = 32, pool = pool, curve = flat)
        @test q.chosen.static_equivalent
        @test SCamp._predictive_declared_budget(q.chosen) == 1
    end
    # No curve: the v1 plan space, priced exactly as before.
    v1 = _planning(n = 32, threads = 16, pool = 16)
    @test all(p -> max(1, p.inner_thread_budget) == 1 || p.route === :none, v1.plans)
    @test all(p -> p.inner_thread_budget <= 1, v1.plans)
    # Every multi-thread plan fits the thread budget: W * b <= T for threads,
    # L * b <= T - 1 for mixed local slots.
    mixed = _planning(n = 32, threads = 16, pool = 8, curve = strong)
    for p in mixed.plans
        b = max(1, p.inner_thread_budget)
        b > 1 || continue
        @test !p.static_equivalent
        p.route === :threads && @test p.workers * b <= 16
        p.route === :process && @test p.local_slots * b <= 15
    end
    @test any(p -> p.route === :process && p.inner_thread_budget > 1, mixed.plans)
end

@testset "the leash walks the budget ladder one rung at a time" begin
    strong = SCamp.InnerSpeedupCurve(Float64.(1:8))
    # Four samples on 32 threads: the shape's ladder is budgets 1, 8, 16, 32
    # (threads@4, threads@4+b8, threads@2+b16, none+b32).
    p = _planning(n = 4, threads = 32, curve = strong)
    @test p.chosen.inner_thread_budget == 8
    plan, why = SCamp.predictive_leash(p, "threads@w4+l0")
    @test why === :within_leash && plan === p.chosen
    # From the serial plan on the whole pool, two rungs away: one rung.
    plan, why = SCamp.predictive_leash(p, "none@w1+l0+b32")
    @test why === :leashed
    @test plan.route === :threads && plan.workers == 2 && plan.inner_thread_budget == 16
    @test SCamp._predictive_parse_plan_key("threads@w8+l0+b4") == (:threads, 8, 0, 4)
    @test SCamp._predictive_parse_plan_key("process@w8+l3") == (:process, 8, 3, 1)
end
