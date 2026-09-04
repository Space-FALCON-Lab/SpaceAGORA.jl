using Test
using TOML

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

using SpaceAGORA
const SimulationModel = SpaceAGORA.SimulationModel
using .SimulationModel

const PP_PROBES = SimulationModel.ParallelPolicy

function _pp_probes_reset_hint_state!()
    lock(PP_PROBES._persistent_hint_lock) do
        state = PP_PROBES._persistent_hint_state[]
        state.loaded = false
        state.dirty = false
        state.path = ""
        empty!(state.history)
    end
    return nothing
end

function _pp_probes_prime_hint_state!(path::String)
    lock(PP_PROBES._persistent_hint_lock) do
        state = PP_PROBES._persistent_hint_state[]
        state.loaded = true
        state.dirty = false
        state.path = path
        empty!(state.history)
    end
    return nothing
end

@testset "Parallel Policy Coverage Debt Probes" begin
    policy = PP_PROBES
    @test Threads.nthreads() > 1

    @testset "policy contexts and scopes" begin
        # Outside any scope the active context is the global one.
        @test policy._active_policy_context() === policy._global_policy_context[]
        global_scope_id = policy._active_policy_scope_id()
        @test global_scope_id == policy._policy_scope_id(policy._global_policy_context[])

        inner_scope_id = Ref{UInt}(0)
        result = policy.with_policy_context() do
            ctx = policy._active_policy_context()
            @test ctx !== policy._global_policy_context[]
            inner_scope_id[] = policy._active_policy_scope_id()
            @test inner_scope_id[] == policy._policy_scope_id(ctx)
            return :scope_return_probe
        end
        @test result === :scope_return_probe
        @test inner_scope_id[] != global_scope_id
        # After the scope closes the TLS entry is gone again.
        @test policy._active_policy_scope_id() == global_scope_id
    end

    @testset "threaded_foreach / worker / reduce (static)" begin
        withenv(
            "SPACEAGORA_INNER_THREAD_BUDGET" => "2",
            "SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER" => "static",
            "SPACEAGORA_PARALLEL_POLICY_INNER_DYNAMIC_CHUNK" => nothing
        ) do
            # Empty input: body never runs.
            touched = Ref(false)
            @test policy.threaded_foreach(0, 4, idx -> (touched[] = true)) === nothing
            @test !touched[]

            # Serial fallback with allotment == 1 preserves ordering.
            serial_order = Int[]
            policy.threaded_foreach(4, 1, idx -> push!(serial_order, idx))
            @test serial_order == [1, 2, 3, 4]

            # Threaded static path, f-first argument order, num_items < allotment.
            acc_static = Base.Threads.Atomic{Int}(0)
            policy.threaded_foreach(idx -> Base.Threads.atomic_add!(acc_static, idx), 5, 8)
            @test acc_static[] == sum(1:5)

            # Errors inside worker bodies propagate out of the @sync.
            err_foreach = try
                policy.threaded_foreach(4, 2, idx -> error("foreach_static_probe"))
                nothing
            catch e
                e
            end
            @test err_foreach !== nothing

            # Worker-count helper edge cases.
            @test policy._thread_worker_count(0, 4) == 1
            @test policy._thread_worker_count(8, 1) == 1
            @test policy._thread_worker_count(8, 8) == 2
            @test policy._thread_worker_count(1, 8) == 1
            @test policy.thread_worker_count(3, 2) == 2

            # threaded_foreach_worker: empty, serial, threaded static partition.
            @test policy.threaded_foreach_worker(0, 2, (w, i) -> error("never")) === nothing
            worker_serial = Int[]
            policy.threaded_foreach_worker(3, 1, (w, i) -> push!(worker_serial, w))
            @test worker_serial == [1, 1, 1]
            worker_ids = zeros(Int, 5)
            policy.threaded_foreach_worker((w, i) -> (worker_ids[i] = w), 5, 2)
            @test all(i -> worker_ids[i] == mod1(i, 2), 1:5)

            # threaded_reduce: empty returns init, serial, threaded static.
            reduce_empty = policy.threaded_reduce(
                0, 2,
                () -> Ref(7),
                (a, i) -> (a[] += i; nothing),
                (d, s) -> (d[] += s[]; nothing)
            )
            @test reduce_empty[] == 7
            reduce_serial = policy.threaded_reduce(
                5, 1,
                () -> Ref(0),
                (a, i) -> (a[] += i; nothing),
                (d, s) -> (d[] += s[]; nothing)
            )
            @test reduce_serial[] == sum(1:5)
            reduce_static = policy.threaded_reduce(
                6, 2,
                () -> Ref(0),
                (a, i) -> (a[] += i; nothing),
                (d, s) -> (d[] += s[]; nothing)
            )
            @test reduce_static[] == sum(1:6)
        end

        # Budget 0 means "use all threads".
        withenv("SPACEAGORA_INNER_THREAD_BUDGET" => "0") do
            @test policy._thread_worker_count(8, 8) == Threads.nthreads()
        end
    end

    @testset "threaded_foreach / worker / reduce (dynamic)" begin
        withenv(
            "SPACEAGORA_INNER_THREAD_BUDGET" => "2",
            "SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER" => "dynamic",
            "SPACEAGORA_PARALLEL_POLICY_INNER_DYNAMIC_CHUNK" => "2"
        ) do
            acc_dyn = Base.Threads.Atomic{Int}(0)
            policy.threaded_foreach(7, 2, idx -> Base.Threads.atomic_add!(acc_dyn, idx))
            @test acc_dyn[] == sum(1:7)

            acc_dyn_w = Base.Threads.Atomic{Int}(0)
            seen_workers = zeros(Int, 7)
            policy.threaded_foreach_worker(7, 2, (w, i) -> begin
                seen_workers[i] = w
                Base.Threads.atomic_add!(acc_dyn_w, i)
            end)
            @test acc_dyn_w[] == sum(1:7)
            @test all(w -> w in (1, 2), seen_workers)

            reduce_dyn = policy.threaded_reduce(
                7, 2,
                () -> Ref(0),
                (a, i) -> (a[] += i; nothing),
                (d, s) -> (d[] += s[]; nothing)
            )
            @test reduce_dyn[] == sum(1:7)
        end
    end

    @testset "persistent channel pools" begin
        withenv(
            "SPACEAGORA_INNER_THREAD_BUDGET" => "2",
            "SPACEAGORA_CALLBACK_PERSISTENT_WORKERS" => "1",
            "SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER" => "static",
            "SPACEAGORA_PARALLEL_POLICY_INNER_DYNAMIC_CHUNK" => nothing
        ) do
            scope_used = policy.with_policy_context() do
                scope_id = policy._active_policy_scope_id()

                acc = Base.Threads.Atomic{Int}(0)
                policy.threaded_foreach_persistent(:probe_chan, 6, 2) do idx
                    Base.Threads.atomic_add!(acc, idx)
                end
                # Second dispatch reuses the same pool (get! hit path).
                policy.threaded_foreach_persistent(
                    :probe_chan, 6, 2, idx -> Base.Threads.atomic_add!(acc, idx)
                )
                @test acc[] == 2 * sum(1:6)

                # Worker-id persistent variant with static partition check.
                ids = zeros(Int, 6)
                policy.threaded_foreach_worker_persistent(:probe_work, 6, 2) do w, i
                    ids[i] = w
                end
                @test all(i -> ids[i] == mod1(i, 2), 1:6)
                policy.threaded_foreach_worker_persistent(
                    :probe_work, 6, 2, (w, i) -> (ids[i] = -w)
                )
                @test all(i -> ids[i] == -mod1(i, 2), 1:6)

                # Error propagation through persistent dispatch (both loop kinds).
                err_chan = try
                    policy.threaded_foreach_persistent(
                        :probe_chan, 4, 2, idx -> error("persistent_chan_probe")
                    )
                    nothing
                catch e
                    e
                end
                @test err_chan isa Base.CapturedException
                err_work = try
                    policy.threaded_foreach_worker_persistent(
                        :probe_work, 4, 2, (w, i) -> error("persistent_work_probe")
                    )
                    nothing
                catch e
                    e
                end
                @test err_work isa Base.CapturedException

                # Dynamic scheduler routed through both persistent worker loops.
                withenv(
                    "SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER" => "dynamic",
                    "SPACEAGORA_PARALLEL_POLICY_INNER_DYNAMIC_CHUNK" => "2"
                ) do
                    acc_dyn = Base.Threads.Atomic{Int}(0)
                    policy.threaded_foreach_persistent(
                        :probe_chan, 7, 2, idx -> Base.Threads.atomic_add!(acc_dyn, idx)
                    )
                    @test acc_dyn[] == sum(1:7)
                    acc_dyn_w = Base.Threads.Atomic{Int}(0)
                    policy.threaded_foreach_worker_persistent(
                        :probe_work, 7, 2, (w, i) -> Base.Threads.atomic_add!(acc_dyn_w, i)
                    )
                    @test acc_dyn_w[] == sum(1:7)
                end

                # Serial and empty fallbacks never touch the pool machinery.
                serial_hits = Ref(0)
                policy.threaded_foreach_persistent(:probe_chan, 3, 1, idx -> (serial_hits[] += 1))
                @test serial_hits[] == 3
                serial_wids = Int[]
                policy.threaded_foreach_worker_persistent(
                    :probe_work, 3, 1, (w, i) -> push!(serial_wids, w)
                )
                @test serial_wids == [1, 1, 1]
                @test policy.threaded_foreach_persistent(:probe_chan, 0, 2, idx -> error("never")) === nothing
                @test policy.threaded_foreach_worker_persistent(:probe_work, 0, 2, (w, i) -> error("never")) === nothing

                lock(policy._persistent_foreach_lock) do
                    @test haskey(policy._persistent_foreach_pools, (scope_id, :probe_chan))
                    @test haskey(policy._persistent_foreach_worker_pools, (scope_id, :probe_work))
                end
                return scope_id
            end

            # Scope teardown removed and shut down both pools.
            lock(policy._persistent_foreach_lock) do
                @test !haskey(policy._persistent_foreach_pools, (scope_used, :probe_chan))
                @test !haskey(policy._persistent_foreach_worker_pools, (scope_used, :probe_work))
            end
        end

        # Persistent workers disabled: falls back to plain threaded loops.
        withenv(
            "SPACEAGORA_INNER_THREAD_BUDGET" => "2",
            "SPACEAGORA_CALLBACK_PERSISTENT_WORKERS" => "0"
        ) do
            acc_off = Base.Threads.Atomic{Int}(0)
            policy.threaded_foreach_persistent(
                :probe_disabled, 5, 2, idx -> Base.Threads.atomic_add!(acc_off, idx)
            )
            @test acc_off[] == sum(1:5)
            acc_off_w = Base.Threads.Atomic{Int}(0)
            policy.threaded_foreach_worker_persistent(
                :probe_disabled, 5, 2, (w, i) -> Base.Threads.atomic_add!(acc_off_w, i)
            )
            @test acc_off_w[] == sum(1:5)
            global_id = policy._active_policy_scope_id()
            lock(policy._persistent_foreach_lock) do
                @test !haskey(policy._persistent_foreach_pools, (global_id, :probe_disabled))
                @test !haskey(policy._persistent_foreach_worker_pools, (global_id, :probe_disabled))
            end
        end

        # Direct pool constructors clamp to at least two workers.
        withenv(
            "SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER" => "static",
            "SPACEAGORA_PARALLEL_POLICY_INNER_DYNAMIC_CHUNK" => nothing
        ) do
            pool_min = policy._create_persistent_foreach_pool(1)
            @test pool_min.workers == 2
            acc_direct = Base.Threads.Atomic{Int}(0)
            policy._threaded_foreach_persistent!(
                pool_min, 3, 2, idx -> Base.Threads.atomic_add!(acc_direct, idx)
            )
            @test acc_direct[] == sum(1:3)
            policy._shutdown_persistent_foreach_pool!(pool_min)

            wpool_min = policy._create_persistent_foreach_worker_pool(1)
            @test wpool_min.workers == 2
            ids_direct = zeros(Int, 4)
            policy._threaded_foreach_persistent!(
                wpool_min, 4, 2, (w, i) -> (ids_direct[i] = w)
            )
            @test all(i -> ids_direct[i] == mod1(i, 2), 1:4)
            policy._shutdown_persistent_foreach_pool!(wpool_min)
        end
        sleep(0.1)  # allow worker tasks to observe :stop before the process exits
    end

    @testset "spin barrier pools" begin
        # Serial and empty paths never allocate a pool.
        spin_serial = Int[]
        policy.threaded_foreach_worker_spin(:spin_probe_serial, 3, 1, (w, i) -> push!(spin_serial, w))
        @test spin_serial == [1, 1, 1]
        @test policy.threaded_foreach_worker_spin(:spin_probe_serial, 0, 2, (w, i) -> error("never")) === nothing

        withenv(
            "SPACEAGORA_INNER_THREAD_BUDGET" => "2",
            "SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER" => "static",
            "SPACEAGORA_PARALLEL_POLICY_INNER_DYNAMIC_CHUNK" => nothing
        ) do
            spin_scope = policy.with_policy_context() do
                scope_id = policy._active_policy_scope_id()

                acc = Base.Threads.Atomic{Int}(0)
                seen = zeros(Int, 5)
                policy.threaded_foreach_worker_spin(:spin_probe, 5, 2) do w, i
                    seen[i] = w
                    Base.Threads.atomic_add!(acc, i)
                end
                @test acc[] == sum(1:5)
                @test all(i -> seen[i] == mod1(i, 2), 1:5)

                # Second dispatch reuses the cached pool (get! hit path).
                acc2 = Base.Threads.Atomic{Int}(0)
                policy.threaded_foreach_worker_spin(
                    :spin_probe, 4, 2, (w, i) -> Base.Threads.atomic_add!(acc2, i)
                )
                @test acc2[] == sum(1:4)

                # Pool-worker error surfaces as the original exception.
                @test_throws ErrorException policy.threaded_foreach_worker_spin(
                    :spin_probe, 4, 2,
                    (w, i) -> (w == 1 && error("spin_worker_probe"); nothing)
                )
                # Coordinator (last worker slot) error path.
                @test_throws ErrorException policy.threaded_foreach_worker_spin(
                    :spin_probe, 4, 2,
                    (w, i) -> (w == 2 && error("spin_coordinator_probe"); nothing)
                )
                # Both worker slots fail: exactly one error propagates and the
                # pool-worker error (lowest worker index) wins over the
                # coordinator's, per the per-worker error-slot ordering.
                both_err = try
                    policy.threaded_foreach_worker_spin(
                        :spin_probe, 4, 2,
                        (w, i) -> error("spin_both_probe_w$(w)")
                    )
                    nothing
                catch err
                    err
                end
                @test both_err isa ErrorException
                @test both_err.msg == "spin_both_probe_w1"
                # A prior round's error is not resurrected by a clean round.
                acc3 = Base.Threads.Atomic{Int}(0)
                policy.threaded_foreach_worker_spin(
                    :spin_probe, 4, 2, (w, i) -> Base.Threads.atomic_add!(acc3, i)
                )
                @test acc3[] == sum(1:4)

                # Dynamic scheduler dispatch through the spin barrier. A
                # rendezvous on each worker's first chunk guarantees both the
                # coordinator and the pool worker execute chunks (so the pool
                # worker's dynamic loop body is deterministically exercised).
                withenv(
                    "SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER" => "dynamic",
                    "SPACEAGORA_PARALLEL_POLICY_INNER_DYNAMIC_CHUNK" => "1"
                ) do
                    entered = (Base.Threads.Atomic{Int}(0), Base.Threads.Atomic{Int}(0))
                    acc_dyn = Base.Threads.Atomic{Int}(0)
                    seen_dyn = zeros(Int, 8)
                    policy.threaded_foreach_worker_spin(:spin_probe, 8, 2) do w, i
                        if entered[w][] == 0
                            entered[w][] = 1
                            while entered[1][] + entered[2][] < 2
                                GC.safepoint()
                            end
                        end
                        seen_dyn[i] = w
                        Base.Threads.atomic_add!(acc_dyn, i)
                    end
                    @test acc_dyn[] == sum(1:8)
                    @test any(==(1), seen_dyn)
                    @test any(==(2), seen_dyn)
                end

                lock(policy._spin_barrier_lock) do
                    @test haskey(policy._spin_barrier_pools, (scope_id, :spin_probe))
                end
                return scope_id
            end

            lock(policy._spin_barrier_lock) do
                @test !haskey(policy._spin_barrier_pools, (spin_scope, :spin_probe))
            end
        end

        # Direct construction clamps workers to [1, nthreads-1]: spin workers
        # never yield, so oversubscribing past nthreads-1 would deadlock the
        # coordinator. Explicit shutdown afterwards.
        spool = policy._create_spin_barrier_pool(0)
        @test spool.workers == 1
        policy._shutdown_spin_barrier_pool!(spool)
        @test spool.stop[] == true
        spool_big = policy._create_spin_barrier_pool(Base.Threads.nthreads() + 7)
        @test spool_big.workers == max(1, Base.Threads.nthreads() - 1)
        policy._shutdown_spin_barrier_pool!(spool_big)
        @test spool_big.stop[] == true
        sleep(0.1)  # allow the spinning workers to observe stop and exit
    end

    @testset "scope teardown covers every pool dictionary" begin
        withenv(
            "SPACEAGORA_INNER_THREAD_BUDGET" => "2",
            "SPACEAGORA_CALLBACK_PERSISTENT_WORKERS" => "1",
            "SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER" => "static"
        ) do
            all_scope = policy.with_policy_context() do
                scope_id = policy._active_policy_scope_id()
                acc = Base.Threads.Atomic{Int}(0)
                policy.threaded_foreach_persistent(
                    :probe_all, 4, 2, idx -> Base.Threads.atomic_add!(acc, idx)
                )
                policy.threaded_foreach_worker_persistent(
                    :probe_all, 4, 2, (w, i) -> Base.Threads.atomic_add!(acc, i)
                )
                policy.threaded_foreach_worker_spin(
                    :probe_all, 4, 2, (w, i) -> Base.Threads.atomic_add!(acc, i)
                )
                @test acc[] == 3 * sum(1:4)
                lock(policy._persistent_foreach_lock) do
                    @test haskey(policy._persistent_foreach_pools, (scope_id, :probe_all))
                    @test haskey(policy._persistent_foreach_worker_pools, (scope_id, :probe_all))
                end
                lock(policy._spin_barrier_lock) do
                    @test haskey(policy._spin_barrier_pools, (scope_id, :probe_all))
                end
                return scope_id
            end
            lock(policy._persistent_foreach_lock) do
                @test !haskey(policy._persistent_foreach_pools, (all_scope, :probe_all))
                @test !haskey(policy._persistent_foreach_worker_pools, (all_scope, :probe_all))
            end
            lock(policy._spin_barrier_lock) do
                @test !haskey(policy._spin_barrier_pools, (all_scope, :probe_all))
            end
        end
        sleep(0.1)
    end

    @testset "thread_execution env toggles" begin
        withenv("SPACEAGORA_CALLBACK_PERSISTENT_WORKERS" => nothing) do
            @test policy.callback_persistent_workers_enabled() == true
        end
        withenv("SPACEAGORA_CALLBACK_PERSISTENT_WORKERS" => "0") do
            @test policy.callback_persistent_workers_enabled() == false
        end
        withenv("SPACEAGORA_HARMONICS_BATCH_SPIN_BARRIER" => nothing) do
            @test policy.harmonics_batch_spin_barrier_enabled() == false
        end
        withenv("SPACEAGORA_HARMONICS_BATCH_SPIN_BARRIER" => "1") do
            @test policy.harmonics_batch_spin_barrier_enabled() == true
        end
    end

    @testset "persistent hints: helpers" begin
        # Bucketing thresholds.
        @test policy._hint_bucket(0) == "1"
        @test policy._hint_bucket(1) == "1"
        @test policy._hint_bucket(2) == "2"
        @test policy._hint_bucket(4) == "3_4"
        @test policy._hint_bucket(8) == "5_8"
        @test policy._hint_bucket(16) == "9_16"
        @test policy._hint_bucket(17) == "17p"

        # Workload signature sanitizes profile/machine tokens.
        withenv(
            "SPACEAGORA_PARALLEL_PROFILE" => "My Profile!",
            "SPACEAGORA_PERF_MACHINE_LABEL" => "Node#7"
        ) do
            sig = policy._hint_workload_signature(:density_callback, 12, 1, 2, true, false, true)
            @test sig == "profile=my_profile_|machine=node_7|src=density_callback|items=9_16|thr=1|budget=2|outer=1|heavy_only=0|heavy=1"
        end

        # Candidate allotments per cap.
        @test policy._hint_candidate_allotments(1, 1) == Int64[1]
        @test policy._hint_candidate_allotments(8, 1) == Int64[1]
        @test policy._hint_candidate_allotments(2, 8) == Int64[1, 2]
        @test policy._hint_candidate_allotments(5, 8) == Int64[1, 2, 3, 4, 5]
        @test policy._hint_candidate_allotments(16, 16) == Int64[1, 2, 4, 8, 16]
        @test policy._hint_candidate_allotments(12, 12) == Int64[1, 2, 4, 6, 8, 12]

        # Mean-and-width math.
        stats = policy.AdaptiveChoiceStats(
            samples=4, successes=4, failures=0, elapsed_sum_ns=400.0, elapsed_sq_sum_ns=40_000.0
        )
        mean_ns, width = policy._hint_mean_and_width(stats, Int64(8), 2.0)
        @test mean_ns == 100.0
        @test isapprox(width, 2.0 * sqrt(log(8.0) / 4); rtol=1e-12)

        bucket = Dict{Int64, policy.AdaptiveChoiceStats}(Int64(2) => stats)
        @test policy._hint_samples(bucket, Int64(2)) == 4
        @test policy._hint_samples(bucket, Int64(3)) == 0

        # Payload clamping and rejection.
        neg_stats = policy.AdaptiveChoiceStats(
            samples=-2, successes=-1, failures=-3, elapsed_sum_ns=-5.0, elapsed_sq_sum_ns=-1.0
        )
        payload = policy._hint_stats_payload(neg_stats)
        @test payload["samples"] == 0
        @test payload["successes"] == 0
        @test payload["failures"] == 0
        @test payload["elapsed_sum_ns"] == 0.0
        @test payload["elapsed_sq_sum_ns"] == 0.0
        @test policy._hint_payload_stats(nothing) === nothing
        @test policy._hint_payload_stats("bad") === nothing
        @test policy._hint_payload_stats(Dict("samples" => 0)) === nothing
        @test policy._hint_payload_stats(Dict("samples" => "oops")) === nothing
        clamped = policy._hint_payload_stats(Dict(
            "samples" => 2, "successes" => 5, "failures" => 9,
            "elapsed_sum_ns" => 10.0, "elapsed_sq_sum_ns" => 100.0
        ))
        @test clamped.samples == 2
        @test clamped.successes == 2
        @test clamped.failures == 0

        # Signature key extraction.
        @test policy._hint_signature_value("profile=r5|machine=m1", "profile") == "r5"
        @test policy._hint_signature_value("profile=|machine=m1", "profile") == ""
        @test policy._hint_signature_value("profile=r5|machine=m1", "absent") === nothing

        # Entry counting across buckets.
        counted = policy._PersistentHintState()
        counted.history["a"] = Dict(
            Int64(1) => policy.AdaptiveChoiceStats(),
            Int64(2) => policy.AdaptiveChoiceStats()
        )
        counted.history["b"] = Dict(Int64(1) => policy.AdaptiveChoiceStats())
        @test policy._hint_entry_count(counted) == 3
    end

    @testset "persistent hints: load / save / reset" begin
        mktempdir() do tmp
            # Missing file: loads to an empty enabled state.
            missing_path = joinpath(tmp, "missing_policy.toml")
            withenv(
                "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "1",
                "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => missing_path,
                "SPACEAGORA_PARALLEL_POLICY_STATE_RESET" => nothing
            ) do
                _pp_probes_reset_hint_state!()
                lock(policy._persistent_hint_lock) do
                    policy._load_persistent_hint_state_locked!()
                    state = policy._persistent_hint_state[]
                    @test state.loaded
                    @test state.path == normpath(missing_path)
                    @test isempty(state.history)
                    # Already-loaded early return leaves state untouched.
                    state.history["marker"] = Dict(Int64(1) => policy.AdaptiveChoiceStats(samples=1))
                    policy._load_persistent_hint_state_locked!()
                    @test haskey(state.history, "marker")
                end
            end

            # Broken TOML swallowed cleanly.
            broken_path = joinpath(tmp, "broken_policy.toml")
            write(broken_path, "history = [")
            withenv(
                "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "1",
                "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => broken_path
            ) do
                _pp_probes_reset_hint_state!()
                lock(policy._persistent_hint_lock) do
                    policy._load_persistent_hint_state_locked!()
                    state = policy._persistent_hint_state[]
                    @test state.loaded
                    @test isempty(state.history)
                end
            end

            # history key that is not an array.
            notvec_path = joinpath(tmp, "notvec_policy.toml")
            write(notvec_path, "schema_version = 1\nhistory = \"nope\"\n")
            withenv(
                "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "1",
                "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => notvec_path
            ) do
                _pp_probes_reset_hint_state!()
                lock(policy._persistent_hint_lock) do
                    policy._load_persistent_hint_state_locked!()
                    @test isempty(policy._persistent_hint_state[].history)
                end
            end

            # Well-formed file with malformed rows interleaved; duplicate rows blend.
            good_path = joinpath(tmp, "good_policy.toml")
            write(good_path, """
            schema_version = 1
            history = [
                "not_a_dict",
                { signature = "", allotment = 1, stats = { samples = 1 } },
                { signature = "sig_load", allotment = 0, stats = { samples = 1 } },
                { signature = "sig_load", allotment = "bad", stats = { samples = 1 } },
                { signature = "sig_load", allotment = 2, stats = "bad" },
                { signature = "sig_load", allotment = 2, stats = { samples = 0 } },
                { signature = "sig_load", allotment = 2, stats = { samples = 3, successes = 2, failures = 1, elapsed_sum_ns = 90.0, elapsed_sq_sum_ns = 2700.0 } },
                { signature = "sig_load", allotment = 2, stats = { samples = 2, successes = 5, failures = 9, elapsed_sum_ns = 10.0, elapsed_sq_sum_ns = 100.0 } },
                { signature = "sig_load", allotment = 4, stats = { samples = 2, successes = "x", failures = "y", elapsed_sum_ns = "z", elapsed_sq_sum_ns = "w" } },
            ]
            """)
            withenv(
                "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "1",
                "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => good_path
            ) do
                _pp_probes_reset_hint_state!()
                lock(policy._persistent_hint_lock) do
                    policy._load_persistent_hint_state_locked!()
                    state = policy._persistent_hint_state[]
                    @test collect(keys(state.history)) == ["sig_load"]
                    blended = state.history["sig_load"][Int64(2)]
                    @test blended.samples == 5
                    @test blended.successes == 4
                    @test blended.failures == 1
                    @test blended.elapsed_sum_ns == 100.0
                    @test blended.elapsed_sq_sum_ns == 2800.0
                    catch_row = state.history["sig_load"][Int64(4)]
                    @test catch_row.samples == 2
                    @test catch_row.successes == 0
                    @test catch_row.failures == 0
                    @test catch_row.elapsed_sum_ns == 0.0
                    @test policy._hint_entry_count(state) == 2
                end
            end

            # State reset request skips loading the existing file.
            withenv(
                "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "1",
                "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => good_path,
                "SPACEAGORA_PARALLEL_POLICY_STATE_RESET" => "1"
            ) do
                _pp_probes_reset_hint_state!()
                lock(policy._persistent_hint_lock) do
                    policy._load_persistent_hint_state_locked!()
                    state = policy._persistent_hint_state[]
                    @test state.loaded
                    @test !state.dirty
                    @test isempty(state.history)
                end
                @test policy.persistent_hints_state_reset_requested()
            end

            # Hints disabled: file present but not loaded.
            withenv(
                "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "0",
                "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => good_path
            ) do
                _pp_probes_reset_hint_state!()
                lock(policy._persistent_hint_lock) do
                    policy._load_persistent_hint_state_locked!()
                    state = policy._persistent_hint_state[]
                    @test state.loaded
                    @test isempty(state.history)
                end
            end

            # Save: writes only non-empty buckets with samples, creates directories,
            # clears dirty, and round-trips through the loader.
            save_path = joinpath(tmp, "nested", "deeper", "saved_policy.toml")
            withenv(
                "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "1",
                "SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => "1",
                "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => save_path
            ) do
                _pp_probes_prime_hint_state!(normpath(save_path))
                lock(policy._persistent_hint_lock) do
                    state = policy._persistent_hint_state[]
                    state.dirty = true
                    state.history["sig_save"] = Dict(
                        Int64(2) => policy.AdaptiveChoiceStats(
                            samples=3, successes=2, failures=1,
                            elapsed_sum_ns=90.5, elapsed_sq_sum_ns=2700.25
                        ),
                        Int64(1) => policy.AdaptiveChoiceStats()  # zero samples: skipped
                    )
                    state.history["sig_empty_bucket"] = Dict{Int64, policy.AdaptiveChoiceStats}()
                    policy._save_persistent_hint_state_locked!()
                    @test !state.dirty
                end
                @test isfile(save_path)
                @test !isfile(save_path * ".tmp")
                parsed = TOML.parsefile(save_path)
                @test parsed["schema_version"] == 1
                rows = parsed["history"]
                @test length(rows) == 1
                @test rows[1]["signature"] == "sig_save"
                @test rows[1]["allotment"] == 2
                @test rows[1]["stats"]["samples"] == 3
                @test rows[1]["stats"]["elapsed_sum_ns"] == 90.5

                # Round-trip: loader reproduces the saved stats.
                _pp_probes_reset_hint_state!()
                lock(policy._persistent_hint_lock) do
                    policy._load_persistent_hint_state_locked!()
                    reloaded = policy._persistent_hint_state[].history["sig_save"][Int64(2)]
                    @test reloaded.samples == 3
                    @test reloaded.successes == 2
                    @test reloaded.failures == 1
                    @test reloaded.elapsed_sum_ns == 90.5
                    @test reloaded.elapsed_sq_sum_ns == 2700.25
                end

                # Not dirty: save is a no-op.
                rm(save_path)
                lock(policy._persistent_hint_lock) do
                    policy._save_persistent_hint_state_locked!()
                end
                @test !isfile(save_path)
            end

            # Persist disabled: dirty state stays dirty and no file appears.
            withenv(
                "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "1",
                "SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => "0",
                "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => save_path
            ) do
                _pp_probes_prime_hint_state!(normpath(save_path))
                lock(policy._persistent_hint_lock) do
                    state = policy._persistent_hint_state[]
                    state.dirty = true
                    state.history["sig_nopersist"] = Dict(
                        Int64(1) => policy.AdaptiveChoiceStats(samples=1, successes=1)
                    )
                    policy._save_persistent_hint_state_locked!()
                    @test state.dirty
                end
                @test !isfile(save_path)
            end

            # ensure-loaded registers the atexit hook exactly once.
            withenv(
                "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "1",
                "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => missing_path
            ) do
                previous_registered = policy._persistent_hint_atexit_registered[]
                policy._persistent_hint_atexit_registered[] = false
                _pp_probes_reset_hint_state!()
                policy._ensure_persistent_hint_state_loaded!()
                @test policy._persistent_hint_atexit_registered[]
                policy._persistent_hint_atexit_registered[] =
                    previous_registered || policy._persistent_hint_atexit_registered[]
            end
        end

        # Public reset clears everything.
        lock(policy._persistent_hint_lock) do
            state = policy._persistent_hint_state[]
            state.loaded = true
            state.dirty = true
            state.path = "reset_probe.toml"
            state.history["reset_sig"] = Dict(Int64(1) => policy.AdaptiveChoiceStats(samples=1))
        end
        policy.reset_persistent_hint_state!()
        lock(policy._persistent_hint_lock) do
            state = policy._persistent_hint_state[]
            @test !state.loaded
            @test !state.dirty
            @test isempty(state.path)
            @test isempty(state.history)
        end
    end

    @testset "persistent hints: choose / record" begin
        mktempdir() do tmp
            hint_path = joinpath(tmp, "choose_probe.toml")

            # Disabled hints always pick allotment 1 without exploring.
            withenv(
                "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "0",
                "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => hint_path
            ) do
                _pp_probes_prime_hint_state!(hint_path)
                disabled = policy._hint_choose_allotment("sig_any", Int64[1, 2])
                @test disabled.allotment == 1
                @test disabled.confidence == 0.0
                @test disabled.regret_ns == 0.0
                @test !disabled.exploring
                # Record is a no-op while disabled.
                policy._hint_record_observation!("sig_disabled", Int64(2), Int64(15); success=true)
                lock(policy._persistent_hint_lock) do
                    @test isempty(policy._persistent_hint_state[].history)
                end
            end

            withenv(
                "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "1",
                "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => hint_path,
                "SPACEAGORA_PARALLEL_POLICY_HINT_MIN_SAMPLES" => "2",
                "SPACEAGORA_PARALLEL_POLICY_HINT_EXPLORATION" => "50.0"
            ) do
                _pp_probes_prime_hint_state!(hint_path)

                # Empty candidate vector degrades to [1]; unknown signature explores.
                unknown = policy._hint_choose_allotment("sig_unknown", Int64[])
                @test unknown.allotment == 1
                @test unknown.exploring
                @test unknown.confidence == Inf
                @test unknown.samples == 0
                # Negative candidates clamp to 1 as well.
                clamped = policy._hint_choose_allotment("sig_unknown", Int64[-3])
                @test clamped.allotment == 1

                lock(policy._persistent_hint_lock) do
                    state = policy._persistent_hint_state[]
                    state.history["sig_explore"] = Dict(
                        Int64(1) => policy.AdaptiveChoiceStats(
                            samples=3, successes=3, elapsed_sum_ns=300.0, elapsed_sq_sum_ns=30_000.0
                        ),
                        Int64(2) => policy.AdaptiveChoiceStats(
                            samples=1, successes=1, elapsed_sum_ns=10.0, elapsed_sq_sum_ns=100.0
                        )
                    )
                    state.history["sig_explore_inf"] = Dict(
                        Int64(1) => policy.AdaptiveChoiceStats(
                            samples=3, successes=3, elapsed_sum_ns=300.0, elapsed_sq_sum_ns=30_000.0
                        )
                    )
                    state.history["sig_best"] = Dict(
                        Int64(1) => policy.AdaptiveChoiceStats(
                            samples=100, successes=100,
                            elapsed_sum_ns=5000.0, elapsed_sq_sum_ns=250_000.0
                        ),
                        Int64(2) => policy.AdaptiveChoiceStats(
                            samples=2, successes=2,
                            elapsed_sum_ns=120.0, elapsed_sq_sum_ns=7200.0
                        )
                    )
                end

                # Under-sampled candidate is explored with a finite confidence width.
                explored = policy._hint_choose_allotment("sig_explore", Int64[1, 2])
                @test explored.allotment == 2
                @test explored.exploring
                @test isfinite(explored.confidence)
                @test explored.confidence > 0.0
                @test explored.samples == 1
                @test explored.regret_ns == 0.0

                # Candidate with zero samples explores with infinite confidence.
                explored_inf = policy._hint_choose_allotment("sig_explore_inf", Int64[1, 2])
                @test explored_inf.allotment == 2
                @test explored_inf.exploring
                @test explored_inf.confidence == Inf
                @test explored_inf.samples == 0

                # Fully sampled bucket: UCB-style score picks the wide-but-fast
                # candidate; regret is the gap to the best observed mean (50 vs 60).
                chosen = policy._hint_choose_allotment("sig_best", Int64[1, 2])
                @test chosen.allotment == 2
                @test !chosen.exploring
                @test chosen.samples == 2
                @test isapprox(chosen.regret_ns, 10.0; rtol=1e-9)
                @test chosen.confidence > 0.0

                # Record: allotment <= 0 ignored; success/failure and clamped
                # negative elapsed accumulate as expected; state goes dirty.
                policy._hint_record_observation!("sig_rec", Int64(0), Int64(9); success=true)
                lock(policy._persistent_hint_lock) do
                    @test !haskey(policy._persistent_hint_state[].history, "sig_rec")
                end
                policy._hint_record_observation!("sig_rec", Int64(3), Int64(10); success=true)
                policy._hint_record_observation!("sig_rec", Int64(3), Int64(-5); success=false)
                lock(policy._persistent_hint_lock) do
                    state = policy._persistent_hint_state[]
                    rec = state.history["sig_rec"][Int64(3)]
                    @test rec.samples == 2
                    @test rec.successes == 1
                    @test rec.failures == 1
                    @test rec.elapsed_sum_ns == 10.0
                    @test rec.elapsed_sq_sum_ns == 100.0
                    @test state.dirty
                end
            end
        end
        _pp_probes_reset_hint_state!()
    end

    @testset "persistent hints: layer stats snapshot" begin
        mktempdir() do tmp
            snap_path = joinpath(tmp, "snapshot_probe.toml")
            withenv(
                "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "1",
                "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => snap_path,
                "SPACEAGORA_PARALLEL_POLICY_HINT_EXPLORATION" => nothing,
                "SPACEAGORA_PARALLEL_POLICY_HINT_MIN_SAMPLES" => nothing
            ) do
                _pp_probes_prime_hint_state!(snap_path)
                @test isempty(policy.hint_layer_stats_snapshot())

                lock(policy._persistent_hint_lock) do
                    state = policy._persistent_hint_state[]
                    state.history["profile=p1|machine=m1|src=density_callback|items=1"] = Dict(
                        Int64(1) => policy.AdaptiveChoiceStats(
                            samples=2, successes=2, elapsed_sum_ns=100.0, elapsed_sq_sum_ns=5000.0
                        ),
                        Int64(2) => policy.AdaptiveChoiceStats(
                            samples=2, successes=1, failures=1,
                            elapsed_sum_ns=60.0, elapsed_sq_sum_ns=1800.0
                        )
                    )
                    state.history["profile=p1|machine=m1|src=density_callback|items=2"] = Dict(
                        Int64(1) => policy.AdaptiveChoiceStats(
                            samples=1, successes=1, elapsed_sum_ns=30.0, elapsed_sq_sum_ns=900.0
                        )
                    )
                    state.history["profile=p1|machine=m1|src=control_callback|items=1"] = Dict(
                        Int64(1) => policy.AdaptiveChoiceStats(
                            samples=4, successes=4, elapsed_sum_ns=200.0, elapsed_sq_sum_ns=10_000.0
                        )
                    )
                    # Missing profile/machine/src tokens map to unknown/other.
                    state.history["foo=bar"] = Dict(
                        Int64(1) => policy.AdaptiveChoiceStats(
                            samples=1, successes=1, elapsed_sum_ns=10.0, elapsed_sq_sum_ns=100.0
                        )
                    )
                    # Zero-sample bucket is skipped entirely.
                    state.history["profile=p1|machine=m1|src=thermal_callback|items=1"] = Dict(
                        Int64(1) => policy.AdaptiveChoiceStats()
                    )
                end

                rows = policy.hint_layer_stats_snapshot()
                @test length(rows) == 3
                @test [(r.profile, r.machine, r.layer) for r in rows] == [
                    ("p1", "m1", "control_callback"),
                    ("p1", "m1", "density_callback"),
                    ("unknown", "unknown", "other")
                ]
                density = rows[2]
                @test density.signature_count == 2
                @test density.choice_count == 3
                @test density.samples_total == 5
                @test density.successes_total == 4
                @test density.failures_total == 1
                @test density.elapsed_mean_ns == 38.0
                @test isapprox(density.elapsed_std_ns, sqrt(96.0); rtol=1e-9)
                @test density.confidence_mean >= 0.0
                @test density.regret_mean_ns >= 0.0
                @test density.exploration_c == 1.5
                @test density.min_samples == 2
                @test density.state_path == snap_path
                control = rows[1]
                @test control.samples_total == 4
                @test control.elapsed_mean_ns == 50.0

                filtered = policy.hint_layer_stats_snapshot(profile="p1", machine="m1")
                @test length(filtered) == 2
                @test all(r -> r.profile == "p1" && r.machine == "m1", filtered)
                @test isempty(policy.hint_layer_stats_snapshot(profile="p1", machine="elsewhere"))
                @test isempty(policy.hint_layer_stats_snapshot(profile="ghost"))
            end
        end
        _pp_probes_reset_hint_state!()
    end
end
