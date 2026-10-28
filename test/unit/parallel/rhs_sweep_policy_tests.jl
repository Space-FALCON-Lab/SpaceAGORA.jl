using Test
using SpaceAGORA
const SEng = SpaceAGORA.SimulationEngine

@testset "skip-by-prediction decision is pure and off by default" begin
    skip = SEng._rhs_calib_skip_by_prediction
    @test skip(0.0, 1e9, 1e8) === nothing           # ratio off
    @test skip(5.0, 0.0, 1e8) === nothing           # no solve prediction
    @test skip(5.0, 1e9, 0.0) === nothing           # no sweep prediction
    @test skip(5.0, 4e8, 1e8) === true              # 4e8 < 5e8
    @test skip(5.0, 5e8, 1e8) === false             # boundary sweeps
    @test skip(5.0, 1e9, 1e8) === false
end

@testset "predicted solve time uses solve_ns, then evals x recorded default call" begin
    pred = SEng._rhs_calib_predicted_solve_ns
    cands = Any[Dict("default" => false, "ns" => 9.0), Dict("default" => true, "ns" => 2.0)]
    @test pred(Dict("solve_ns" => 3e9), 100.0) == 3e9
    @test pred(Dict("candidates" => cands), 100.0) == 200.0
    @test pred(Dict("candidates" => cands), 0.0) == 0.0
    @test pred(Dict{String, Any}(), 100.0) == 0.0
end

@testset "skip ratio env gates the solve threshold" begin
    entry = Dict{String, Any}("mode" => "heuristic", "solve_ns" => 4e8, "sweep_ns" => 1e8)
    mktempdir() do dir
        withenv("SPACEAGORA_RHS_CALIBRATION_PATH" => joinpath(dir, "s.toml"),
                "SPACEAGORA_RHS_CALIBRATE_MIN_SOLVE_S" => "1.0") do
            lock(SEng._rhs_calib_lock) do
                SEng._rhs_calib_loaded[] = true; SEng._rhs_calib_loaded_path[] = SEng._rhs_calib_path()
                empty!(SEng._rhs_calib_cache); SEng._rhs_calib_cache["sig"] = entry
            end
            try
                # 0.4 s solve < MIN_SOLVE_S: short either way with the ratio off.
                @test !SEng._rhs_calib_solve_exceeds_threshold("sig")
                # Ratio 2: 4e8 >= 2 * 1e8 -> long, which the fixed rule would call short.
                withenv("SPACEAGORA_RHS_CALIBRATE_SKIP_RATIO" => "2") do
                    @test SEng._rhs_calib_solve_exceeds_threshold("sig")
                end
                withenv("SPACEAGORA_RHS_CALIBRATE_SKIP_RATIO" => "5") do
                    @test !SEng._rhs_calib_solve_exceeds_threshold("sig")
                end
            finally
                lock(SEng._rhs_calib_lock) do
                    empty!(SEng._rhs_calib_cache); SEng._rhs_calib_loaded[] = false; SEng._rhs_calib_loaded_path[] = ""
                end
            end
        end
    end
end

@testset "seeded subset: nearest to the heuristic plus serial reference and heuristic" begin
    F = SEng._make_calib_flat_plan
    cands = Any[SEng._make_calib_satellite_batch_plan()]
    for s in (:static, :dynamic), a in (2, 4, 8, 16, 32)
        push!(cands, F(a, s))
    end
    push!(cands, F(1, :static))
    h = F(16, :static)
    sub = SEng._rhs_seed_subset(cands, h, 3)
    key(c) = (c.mode, c.allotment, c.scheduler)
    @test key(h) in key.(sub)                          # heuristic kept
    @test key(F(1, :static)) in key.(sub)              # serial reference kept
    @test key(F(16, :dynamic)) in key.(sub)            # same width, other scheduler
    @test key(F(8, :static)) in key.(sub) || key(F(32, :static)) in key.(sub)
    @test !(key(SEng._make_calib_satellite_batch_plan()) in key.(sub))
    @test length(sub) < length(cands)
    @test length(SEng._rhs_seed_subset(cands, h, 100)) == length(cands)
    @test key.(sub) == key.(filter(c -> key(c) in key.(sub), cands))   # original order kept
end

@testset "seeded result falls back to the full sweep when untrustworthy" begin
    rec(mode, a, s, ns; dflt = false, round = 0) = Dict{String, Any}(
        "mode" => mode, "allotment" => a, "scheduler" => s, "ns" => ns,
        "default" => dflt, "round" => round)
    flat = "flat_constellation_effector_queue"
    full = SEng._rhs_seed_needs_full_sweep
    ok = [rec(flat, 1, "static", 100.0), rec(flat, 16, "static", 50.0; dflt = true),
          rec(flat, 8, "static", 52.0)]
    @test !full(ok, 0.10, 0.10)
    # best does not clear the serial reference by the override margin
    @test full([rec(flat, 1, "static", 100.0), rec(flat, 16, "static", 95.0; dflt = true)], 0.10, 0.10)
    # a seeded plan beats the heuristic by more than the disagreement margin
    @test full([rec(flat, 1, "static", 100.0), rec(flat, 16, "static", 80.0; dflt = true),
                rec(flat, 8, "static", 50.0)], 0.10, 0.10)
    # ... but not when the margin is wide enough
    @test !full([rec(flat, 1, "static", 100.0), rec(flat, 16, "static", 80.0; dflt = true),
                 rec(flat, 8, "static", 50.0)], 0.10, 0.50)
    # a round-0 reading supersedes a coarser earlier one for the same plan
    @test !full([rec(flat, 1, "static", 100.0), rec(flat, 16, "static", 30.0; dflt = true, round = 1),
                 rec(flat, 16, "static", 50.0; dflt = true, round = 0)], 0.10, 0.10)
    @test !full(Dict{String, Any}[], 0.10, 0.10)
end

@testset "sweep_seeded round-trips through the store and is absent by default" begin
    mktempdir() do dir
        path = joinpath(dir, "s.toml")
        withenv("SPACEAGORA_RHS_CALIBRATION_PATH" => path) do
            lock(SEng._rhs_calib_lock) do
                SEng._rhs_calib_loaded[] = false; SEng._rhs_calib_loaded_path[] = ""; empty!(SEng._rhs_calib_cache)
            end
            try
                SEng._rhs_calib_store_heuristic!("a", 1.0; sweep_ns = 1.0, seeded = true)
                SEng._rhs_calib_store_heuristic!("b", 1.0; sweep_ns = 1.0)
                SEng._rhs_calib_save!()
                txt = read(path, String)
                @test count("sweep_seeded", txt) == 1
                lock(SEng._rhs_calib_lock) do
                    SEng._rhs_calib_loaded[] = false; SEng._rhs_calib_loaded_path[] = ""
                end
                SEng._rhs_calib_load!()
                @test get(SEng._rhs_calib_cache["a"], "sweep_seeded", false) === true
                @test !haskey(SEng._rhs_calib_cache["b"], "sweep_seeded")
            finally
                lock(SEng._rhs_calib_lock) do
                    empty!(SEng._rhs_calib_cache); SEng._rhs_calib_loaded[] = false; SEng._rhs_calib_loaded_path[] = ""
                end
            end
        end
    end
end

@testset "accept-winner waives only the beats-heuristic fallback, default off" begin
    fb = SEng._rhs_seed_should_fallback
    @test !fb(:none, false) && !fb(:none, true)
    @test fb(:beats_heuristic, false)
    @test !fb(:beats_heuristic, true)
    @test fb(:no_gain_over_reference, false) && fb(:no_gain_over_reference, true)
    withenv("SPACEAGORA_RHS_CALIBRATE_SEED_ACCEPT_WINNER" => nothing) do
        @test !SEng._rhs_calibrate_seed_accept_winner()
    end
    withenv("SPACEAGORA_RHS_CALIBRATE_SEED_ACCEPT_WINNER" => "1") do
        @test SEng._rhs_calibrate_seed_accept_winner()
    end
end

@testset "sweep_ns floor keeps max(new, previous), default off" begin
    f = SEng._rhs_stored_sweep_ns
    prev = Dict{String, Any}("sweep_ns" => 2.0e9)
    @test f(prev, 3.0e8, false) == 3.0e8
    @test f(prev, 3.0e8, true) == 2.0e9
    @test f(prev, 3.0e9, true) == 3.0e9
    @test f(nothing, 3.0e8, true) == 3.0e8
    @test f(Dict{String, Any}(), 3.0e8, true) == 3.0e8
    withenv("SPACEAGORA_RHS_CALIBRATE_SWEEP_NS_FLOOR" => nothing) do
        @test !SEng._rhs_calibrate_sweep_ns_floor()
    end
    withenv("SPACEAGORA_RHS_CALIBRATE_SWEEP_NS_FLOOR" => "1") do
        @test SEng._rhs_calibrate_sweep_ns_floor()
    end
    # store path: second (cheap) sweep keeps the first sweep's sweep_ns when on
    sig = "sweepns_floor_test"
    plan = (mode = :flat, allotment = 2, scheduler = :static)
    for (flag, want) in (("0", 3.0e8), ("1", 2.0e9))
        withenv("SPACEAGORA_RHS_CALIBRATE_SWEEP_NS_FLOOR" => flag) do
            delete!(SEng._rhs_calib_cache, sig)
            SEng._rhs_calib_store!(sig, plan, 1.0; sweep_ns = 2.0e9)
            SEng._rhs_calib_store!(sig, plan, 1.0; sweep_ns = 3.0e8)
            @test SEng._rhs_calib_cache[sig]["sweep_ns"] == want
            delete!(SEng._rhs_calib_cache, sig)
        end
    end
end

@testset "pin margin: off unchanged, inside rejected, outside pinned" begin
    acc = SEng._rhs_pin_accepts
    ov = 0.10
    # off (0.0): effective margin is the override margin, as before
    @test max(ov, 0.0) == ov
    @test acc(0.85, 1.0, max(ov, 0.0))          # 15% better, pinned
    @test !acc(0.95, 1.0, max(ov, 0.0))         # 5% better, rejected
    # pin margin 0.30: 15% winner now inside the margin, 40% winner outside
    @test !acc(0.85, 1.0, max(ov, 0.30))
    @test acc(0.60, 1.0, max(ov, 0.30))
    @test !acc(Inf, 1.0, 0.1) && acc(0.5, Inf, 0.1)
    withenv("SPACEAGORA_RHS_CALIBRATE_PIN_MARGIN" => nothing) do
        @test SEng._rhs_calibrate_pin_margin() == 0.0
    end
    withenv("SPACEAGORA_RHS_CALIBRATE_PIN_MARGIN" => "0.25") do
        @test SEng._rhs_calibrate_pin_margin() == 0.25
    end
    for bad in ("-0.1", "1", "1.5", "abc", "NaN")
        withenv("SPACEAGORA_RHS_CALIBRATE_PIN_MARGIN" => bad) do
            @test_throws ArgumentError SEng._rhs_calibrate_pin_margin()
        end
    end
end
