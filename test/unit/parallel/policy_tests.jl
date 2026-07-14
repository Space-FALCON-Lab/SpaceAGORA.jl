@testset "Coverage Threaded Probe Driver" begin
    if Base.JLOptions().code_coverage == 0
        @test true
    else
        probe_script = joinpath(REPO_ROOT, "test", "coverage_threaded_probes.jl")
        cmd = `$(Base.julia_cmd()) --startup-file=no --depwarn=error --project=$(REPO_ROOT) --code-coverage=user --threads=2 $(probe_script)`
        cmd = addenv(
            cmd,
            "SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0",
            "SPACEAGORA_WARN_NORMALIZE" => "0"
        )

        output = IOBuffer()
        proc = run(pipeline(ignorestatus(cmd), stdout=output, stderr=output))
        text = String(take!(output))
        if !success(proc)
            println(text)
        end

        @test success(proc)
        @test occursin("coverage_threaded_probes_ok", text)
    end
end

@testset "Multibody Parallel Policy Gates" begin
    use_threads = SimulationModel.DynamicEffectors._multibody_use_threads
    has_worker_threads = Threads.nthreads() > 1

    withenv(
        "SPACEAGORA_MULTIBODY_PARALLEL" => "auto",
        "SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "1",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1"
    ) do
        @test use_threads(64) == false
    end

    withenv(
        "SPACEAGORA_MULTIBODY_PARALLEL" => "auto",
        "SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "2",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0"
    ) do
        @test use_threads(64) == has_worker_threads
    end

    withenv(
        "SPACEAGORA_MULTIBODY_PARALLEL" => "on",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1"
    ) do
        @test use_threads(64) == has_worker_threads
    end
end

@testset "Parallel Policy Adaptive Controller" begin
    policy = SimulationModel.ParallelPolicy
    policy.reset_policy_telemetry!()

    withenv(
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
        "SPACEAGORA_PARALLEL_POLICY_WINDOW" => "1",
        "SPACEAGORA_PARALLEL_POLICY_TRIM_QUANTA" => "1",
        "SPACEAGORA_PARALLEL_POLICY_DELTA" => "0.8",
        "SPACEAGORA_PARALLEL_POLICY_RHO" => "1.5",
        "SPACEAGORA_INNER_THREAD_BUDGET" => "1"
    ) do
        @test policy.use_threads_policy(
            8;
            mode=:auto,
            threshold=1,
            source=:density_callback
        ) == false

        policy.record_policy_observation!(
            :density_callback;
            mode=:auto,
            num_items=1,
            use_threads=false,
            elapsed_ns=10
        )
        snap1 = policy.policy_telemetry_snapshot()
        @test snap1.last_classification == "efficient_satisfied"
        @test snap1.adaptation_updates_total >= 1
        @test snap1.last_desire >= 2

        policy.record_policy_observation!(
            :density_callback;
            mode=:auto,
            num_items=1,
            use_threads=false,
            elapsed_ns=11
        )
        snap2 = policy.policy_telemetry_snapshot()
        @test snap2.last_classification == "efficient_deprived"

        policy.record_policy_observation!(
            :density_callback;
            mode=:auto,
            num_items=0,
            use_threads=false,
            elapsed_ns=12
        )
        snap3 = policy.policy_telemetry_snapshot()
        @test snap3.last_classification == "inefficient"
        @test snap3.last_utilization <= 0.1
        @test snap3.serial_elapsed_ns_total >= 33
        @test snap3.quantum_length == 1
        @test snap3.trim_quanta_budget == 1
        @test snap3.quantums_total >= 3
        @test snap3.quantums_inefficient >= 1
        @test snap3.quantums_efficient_satisfied >= 1
        @test snap3.quantums_efficient_deprived >= 1
        @test snap3.accounted_fraction_proxy >= 0.0
        @test snap3.trimmed_accounted_fraction_proxy >= 0.0
    end

    if Threads.nthreads() > 1
        policy.reset_policy_telemetry!()
        withenv(
            "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
            "SPACEAGORA_PARALLEL_POLICY_WINDOW" => "1",
            "SPACEAGORA_PARALLEL_POLICY_DELTA" => "0.8",
            "SPACEAGORA_PARALLEL_POLICY_RHO" => "1.5",
            "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads())
        ) do
            use_history = Bool[]
            for _ in 1:8
                decision = policy.thread_policy_decision(
                    6;
                    mode=:auto,
                    threshold=1,
                    source=:other_source
                )
                push!(use_history, decision.use_threads)
                policy.record_policy_observation!(
                    :other_source;
                    mode=:auto,
                    num_items=6,
                    use_threads=decision.use_threads,
                    elapsed_ns=1
                )
            end
            @test any(use_history)
            @test use_history[end] == true
        end

        policy.reset_policy_telemetry!()
        withenv(
            "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
            "SPACEAGORA_PARALLEL_POLICY_BOOTSTRAP_THREADS" => "1",
            "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads())
        ) do
            decision = policy.thread_policy_decision(
                8;
                mode=:auto,
                threshold=2,
                source=:control_callback
            )
            @test decision.use_threads == true
            @test decision.allotment >= 2
        end

        policy.reset_policy_telemetry!()
        withenv(
            "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
            "SPACEAGORA_PARALLEL_POLICY_BOOTSTRAP_THREADS" => "0",
            "SPACEAGORA_PARALLEL_POLICY_CONTROL_TAIL_GUARD" => "0",
            "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads())
        ) do
            decision = policy.thread_policy_decision(
                8;
                mode=:auto,
                threshold=2,
                source=:control_callback
            )
            @test decision.use_threads == false
            @test decision.allotment == 1
        end

        policy.reset_policy_telemetry!()
        withenv(
            "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
            "SPACEAGORA_PARALLEL_POLICY_BOOTSTRAP_THREADS" => "0",
            "SPACEAGORA_PARALLEL_POLICY_CONTROL_TAIL_GUARD" => "1",
            "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads())
        ) do
            decision = policy.thread_policy_decision(
                8;
                mode=:auto,
                threshold=2,
                source=:control_callback
            )
            @test decision.use_threads == true
            @test decision.allotment == min(Threads.nthreads(), 8)
        end
    end

    withenv(
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "0",
        "SPACEAGORA_INNER_THREAD_BUDGET" => "1"
    ) do
        @test policy.use_threads_policy(4; mode=:on, threshold=1, source=:other_source) == false
        policy.record_policy_observation!(
            :other_source;
            mode=:on,
            num_items=4,
            use_threads=true,
            elapsed_ns=UInt(5)
        )
        snap = policy.policy_telemetry_snapshot()
        @test snap.threaded_elapsed_ns_total >= 5
        @test snap.other_decisions >= 1
    end

    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => "oops") do
        @test_throws ArgumentError policy.effective_inner_thread_budget()
    end
end

@testset "Parallel Policy threaded_reduce" begin
    policy = SimulationModel.ParallelPolicy
    budget = max(1, Threads.nthreads())
    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => string(budget)) do
        reduced = policy.threaded_reduce(
            16,
            budget,
            () -> MVector{2, Int}(0, 0),
            (local_acc, idx) -> begin
                local_acc[1] += idx
                local_acc[2] += 1
                return nothing
            end,
            (dest, src) -> begin
                dest[1] += src[1]
                dest[2] += src[2]
                return nothing
            end
        )
        @test reduced[1] == sum(1:16)
        @test reduced[2] == 16

        reduced_empty = policy.threaded_reduce(
            0,
            budget,
            () -> Ref(0),
            (local_acc, idx) -> begin
                local_acc[] += idx
                return nothing
            end,
            (dest, src) -> begin
                dest[] += src[]
                return nothing
            end
        )
        @test reduced_empty[] == 0
    end
end

@testset "Parallel Policy threaded_foreach_persistent" begin
    policy = SimulationModel.ParallelPolicy
    budget = max(1, Threads.nthreads())

    withenv(
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(budget),
        "SPACEAGORA_CALLBACK_PERSISTENT_WORKERS" => "1"
    ) do
        acc = Base.Threads.Atomic{Int}(0)
        policy.with_policy_context() do
            policy.threaded_foreach_persistent(:density_callback, 16, budget) do idx
                Base.Threads.atomic_add!(acc, idx)
            end
            policy.threaded_foreach_persistent(:density_callback, 16, budget) do idx
                Base.Threads.atomic_add!(acc, idx)
            end
        end
        @test acc[] == 2 * sum(1:16)
    end

    withenv(
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(budget),
        "SPACEAGORA_CALLBACK_PERSISTENT_WORKERS" => "0"
    ) do
        acc = Base.Threads.Atomic{Int}(0)
        policy.threaded_foreach_persistent(:density_callback, 8, budget) do idx
            Base.Threads.atomic_add!(acc, idx)
        end
        @test acc[] == sum(1:8)
    end

    withenv(
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(budget),
        "SPACEAGORA_CALLBACK_PERSISTENT_WORKERS" => "1"
    ) do
        err = try
            policy.with_policy_context() do
                policy.threaded_foreach_persistent(:thermal_callback, 8, budget) do idx
                    if idx == 3
                        error("threaded_foreach_persistent_probe")
                    end
                end
            end
            nothing
        catch e
            e
        end
        @test err !== nothing
    end
end

