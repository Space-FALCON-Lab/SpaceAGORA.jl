using Test
using SpaceAGORA

# Pre-solve measured inner-width search (rhs_calibration.jl, "Pre-solve
# inner-width search"). Timing is mocked: no solve runs here.

const SEw = SpaceAGORA.SimulationEngine
const PPw = SpaceAGORA.SimulationModel.ParallelPolicy

# A mock `w -> (t_budget, t_w)` from a table of per-width times; counts calls.
function _mock_timing(times::Dict{Int, Float64}, budget::Int, calls::Vector{Int})
    return w -> (push!(calls, w); (times[budget], times[w]))
end

function _with_clean_store(f)
    mktempdir() do dir
        withenv("SPACEAGORA_RHS_CALIBRATION_PATH" => joinpath(dir, "store.toml")) do
            lock(SEw._rhs_calib_lock) do
                empty!(SEw._rhs_calib_cache); SEw._rhs_calib_loaded[] = false; SEw._rhs_calib_loaded_path[] = ""
            end
            try f() finally
                lock(SEw._rhs_calib_lock) do
                    empty!(SEw._rhs_calib_cache); SEw._rhs_calib_loaded[] = false; SEw._rhs_calib_loaded_path[] = ""
                end
            end
        end
    end
end

@testset "width ladder: B, 3B/4, B/2, B/4, ... rounded to integers >= 1" begin
    @test SEw._width_search_rungs(32) == [32, 24, 16, 8, 4, 2, 1]
    @test SEw._width_search_rungs(12) == [12, 9, 6, 3, 2, 1]
    @test SEw._width_search_rungs(2) == [2, 1]
    @test SEw._width_search_rungs(3) == [3, 2, 1]
    @test SEw._width_search_rungs(1) == [1]
    @test all(issorted(SEw._width_search_rungs(b); rev = true) && last(SEw._width_search_rungs(b)) >= 1 for b in 1:64)
end

@testset "search: descends while winning, stops at the first slower rung" begin
    B = 32
    calls = Int[]
    # Flat past 16 (GRAM-lookahead shaped): 24 and 16 win, 8 loses -> 16.
    t = Dict(32 => 10.95, 24 => 10.6, 16 => 10.3, 8 => 11.5, 4 => 20.0, 2 => 30.0, 1 => 60.0)
    chosen, ratios = SEw._search_inner_width(B, _mock_timing(t, B, calls))
    @test chosen == 16
    @test calls == [24, 16, 8]                 # 4 never measured: 8 lost
    @test last(ratios[end]) > 1.0

    # Scaling workload (gravity shaped): nothing narrower wins; stop after B/2.
    empty!(calls)
    t = Dict(32 => 1.04, 24 => 1.20, 16 => 1.37, 8 => 2.5, 4 => 5.0, 2 => 9.0, 1 => 18.0)
    chosen, _ = SEw._search_inner_width(B, _mock_timing(t, B, calls))
    @test chosen == B
    @test calls == [24, 16]                    # 3B/4 and B/2 always, then stop

    # Keeps winning all the way down.
    empty!(calls)
    t = Dict(32 => 10.0, 24 => 9.0, 16 => 8.0, 8 => 7.0, 4 => 6.0, 2 => 5.0, 1 => 4.0)
    chosen, _ = SEw._search_inner_width(B, _mock_timing(t, B, calls))
    @test chosen == 1
    @test calls == [24, 16, 8, 4, 2, 1]

    # 3B/4 slower than B but B/2 faster: the base rungs are always measured.
    empty!(calls)
    t = Dict(32 => 10.0, 24 => 10.4, 16 => 9.0, 8 => 9.5, 4 => 20.0, 2 => 30.0, 1 => 60.0)
    chosen, _ = SEw._search_inner_width(B, _mock_timing(t, B, calls))
    @test chosen == 16
    @test calls == [24, 16, 8]
end

@testset "search: a narrower width needs the 3% margin over the budget" begin
    @test SEw._WIDTH_SEARCH_MARGIN == 0.03
    B = 8
    for (ratio, want) in ((0.98, 8), (0.9705, 8), (0.97, 6), (0.90, 6), (1.0, 8))
        t = Dict(8 => 100.0, 6 => 100.0 * ratio, 4 => 200.0, 2 => 400.0, 1 => 800.0)
        chosen, _ = SEw._search_inner_width(B, _mock_timing(t, B, Int[]))
        @test chosen == want
    end
    # A failed or non-finite reading is a loss, never a win.
    chosen, _ = SEw._search_inner_width(8, w -> (100.0, NaN))
    @test chosen == 8
end

@testset "verdict persists in the RHS store and reloads from disk" begin
    _with_clean_store() do
        sig = SEw._WIDTH_SEARCH_PREFIX * "sig-a"
        @test SEw._inner_width_lookup(sig) === nothing
        SEw._inner_width_store!(sig, 16)
        @test SEw._inner_width_lookup(sig) == 16
        SEw._rhs_calib_save!()
        lock(SEw._rhs_calib_lock) do
            empty!(SEw._rhs_calib_cache); SEw._rhs_calib_loaded[] = false; SEw._rhs_calib_loaded_path[] = ""
        end
        @test SEw._inner_width_lookup(sig) == 16     # a later process reads it
        # A plan row under the same stem is never read as a width.
        SEw._rhs_calib_store!("plain", SEw._make_calib_flat_plan(4, :static), 1.0e6)
        @test SEw._inner_width_lookup("plain") === nothing
    end
end

# A stand-in for the solve's params: enough for the signature and the gates.
function _stub(; v2 = true, adaptive = true)
    penv = withenv("SPACEAGORA_PARALLEL_POLICY_V2" => v2 ? "1" : "0",
                   "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => adaptive ? "1" : "0") do
        PPw.snapshot_policy_decision_env()
    end
    sb = SpaceAGORA.SimulationModel.ConfigTypes.SharedBuffers(n_sats = 4)
    sb.policy_env_config[] = penv
    p = (shared_buffers = sb, is_active = fill(true, 4))
    args = (dynamics_model = (dynamic_effectors = (1,),),
            environment_model = (density_model = :stub,),
            simulation_settings = (verbose = false,))
    return p, args, sb
end

@testset "calibrate: search, persist, reuse, apply, switches" begin
    if Threads.nthreads() < 2
        @test_skip Threads.nthreads() >= 2
    else
    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => string(min(8, Threads.nthreads())),
            "SPACEAGORA_RHS_CALIBRATE" => "auto", "SPACEAGORA_RHS_WIDTH_SEARCH" => nothing,
            "SPACEAGORA_RHS_WIDTH_SEARCH_TRACE" => "0") do
        B = PPw.effective_inner_thread_budget()
        narrow = max(1, B ÷ 2)
        # Every width below B wins by 10%, but B/2's rung is the last that wins.
        times = Dict(w => (w > narrow ? 100.0 - 1.0 * (B - w) : w == narrow ? 80.0 : 1000.0)
                     for w in 1:B)
        times[B] = 100.0
        _with_clean_store() do
            p, args, sb = _stub()
            calls = Int[]
            SEw._calibrate_inner_width!(p, [0.0], args; timing = _mock_timing(times, B, calls))
            @test !isempty(calls)
            @test sb.inner_width_chosen[] == narrow
            @test sb.policy_env_config[].inner_thread_budget == narrow
            @test SEw._solve_inner_width(p) == narrow
            # A later solve (new buffers) reuses the stored verdict without timing.
            p2, args2, sb2 = _stub()
            calls2 = Int[]
            SEw._calibrate_inner_width!(p2, [0.0], args2; timing = _mock_timing(times, B, calls2))
            @test isempty(calls2)
            @test sb2.policy_env_config[].inner_thread_budget == narrow
            @test sb2.inner_width_chosen[] == narrow
            # The off switch leaves the solve at the full budget, store untouched.
            withenv("SPACEAGORA_RHS_WIDTH_SEARCH" => "off") do
                p3, args3, sb3 = _stub()
                SEw._calibrate_inner_width!(p3, [0.0], args3; timing = _mock_timing(times, B, Int[]))
                @test sb3.inner_width_chosen[] == 0
                @test sb3.policy_env_config[].inner_thread_budget == B
            end
        end
        # A pinned flat plan is clamped to the chosen width.
        _with_clean_store() do
            p, args, sb = _stub()
            sb.rhs_plan_override[] = SEw._make_calib_flat_plan(B, :static)
            SEw._calibrate_inner_width!(p, [0.0], args; timing = _mock_timing(times, B, Int[]))
            @test sb.rhs_plan_override[].allotment == narrow
        end
        # Static modes and calibration off are untouched.
        _with_clean_store() do
            for (kw, env) in ((:static, ()), (:calib_off, ("SPACEAGORA_RHS_CALIBRATE" => "off",)))
                p, args, sb = kw === :static ? _stub(v2 = false, adaptive = false) : _stub()
                withenv(env...) do
                    SEw._calibrate_inner_width!(p, [0.0], args; timing = _mock_timing(times, B, Int[]))
                end
                @test sb.inner_width_chosen[] == 0
                @test sb.policy_env_config[].inner_thread_budget == B
            end
        end
    end
    end
end

# A function, so the closure the spinners run is not newer than their world.
function _spin_dispatch_check(n)
    seen = zeros(Int, 40); ids = Set{Int}(); lk = ReentrantLock()
    PPw.threaded_foreach_worker_spin(:width_pool_test, 40, n; pool_width = 3) do w, i
        seen[i] += 1; lock(lk) do; push!(ids, w); end
    end
    return seen, ids
end

@testset "spin-barrier pool is sized to the chosen width" begin
    n = Threads.nthreads()
    if n < 4
        @test_skip n >= 4
    else
        scope = PPw._active_policy_scope_id()
        try
            @test PPw._spin_pool_worker_count(0) == n - 1          # shipped size
            @test PPw._spin_pool_worker_count(2) == 1
            @test PPw._spin_pool_worker_count(n + 10) == n - 1     # never past the pool
            full = PPw._spin_barrier_pool_for(:width_pool_test, 0)
            @test full.workers == n - 1
            narrow = PPw._spin_barrier_pool_for(:width_pool_test, 3)   # replaces the wider pool
            @test narrow.workers == 2
            @test full.stop[]                                          # old spinners stopped
            @test PPw._spin_barrier_pool_for(:width_pool_test, 3) === narrow
            @test PPw._spin_barrier_pool_for(:width_pool_test, 0) === narrow   # width 0 never resizes
            # Dispatch under a width runs every item exactly once, on <= width workers.
            seen, ids = _spin_dispatch_check(n)
            @test all(==(1), seen)
            @test length(ids) <= 3
        finally
            PPw._destroy_persistent_foreach_scope!(scope)
        end
    end
end
