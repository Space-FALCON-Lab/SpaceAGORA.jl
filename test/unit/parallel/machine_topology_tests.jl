using Test
using SpaceAGORA

const MT = SpaceAGORA.ParallelProfiles

@testset "machine topology" begin
    @testset "cpu list parsing" begin
        @test MT._parse_cpu_list("0-3") == Set([0, 1, 2, 3])
        @test MT._parse_cpu_list("0,2,4") == Set([0, 2, 4])
        @test MT._parse_cpu_list("0-1,8,12-13") == Set([0, 1, 8, 12, 13])
        @test MT._parse_cpu_list("7") == Set([7])
        @test isempty(MT._parse_cpu_list(""))
        # Malformed entries are dropped, never thrown on: an unreadable mask
        # must degrade to "unknown", which callers treat as no constraint.
        @test MT._parse_cpu_list("0-,x,3") == Set([3])
    end

    @testset "physical cores" begin
        n = MT.physical_core_count()
        @test n >= 1
        # SMT can only make the logical count larger, never smaller.
        @test n <= Sys.CPU_THREADS
    end

    @testset "cgroup quota" begin
        q = MT.cgroup_cpu_quota()
        # Unlimited/unreadable is -1.0; anything else must be a positive
        # fractional core count. Zero and negatives other than the sentinel
        # would silently produce a zero budget downstream.
        @test q == -1.0 || q > 0.0
    end

    @testset "snapshot is internally consistent" begin
        t = MT.refresh_machine_topology!()
        @test t.cpu_threads >= 1
        @test t.physical_cores >= 1
        @test t.usable_cores >= 1
        @test t.usable_cores <= t.physical_cores
        @test t.smt_ratio >= 1.0
        @test t.affinity_cores == -1 || 1 <= t.affinity_cores <= t.physical_cores
        @test t.source in (:physical, :affinity, :quota, :override)
        @test MT.usable_core_budget() == t.usable_cores
    end

    @testset "memory" begin
        t = MT.refresh_machine_topology!()
        @test t.total_memory > 0
        @test 0 < t.usable_memory <= t.total_memory
        @test t.cgroup_memory == -1 || t.cgroup_memory > 0
        @test t.memory_source in (:total, :cgroup)
        @test MT.process_rss_bytes() > 0
        @test MT.available_memory_bytes() > 0
        @test MT.memory_worker_cap() >= 0
        # A worker never costs less than the floor unless overridden.
        @test MT.worker_memory_estimate_bytes() >= MT._WORKER_MEMORY_FLOOR_BYTES
        @test MT.worker_memory_estimate_bytes(extra = 7) == MT.worker_memory_estimate_bytes() + 7
        withenv("SPACEAGORA_MEMORY_BUDGET_GB" => "1024", "SPACEAGORA_PERF_WORKER_MEMORY_GB" => "0.001") do
            @test MT.memory_budget_bytes() == 1024 * 2^30
            @test MT.memory_worker_cap() >= 2
        end
        # A budget smaller than this process's own footprint fits no worker.
        withenv("SPACEAGORA_MEMORY_BUDGET_GB" => "0.001") do
            @test MT.memory_worker_cap() == 0
        end
        for bad in ("0", "-1", "abc")
            withenv("SPACEAGORA_MEMORY_BUDGET_GB" => bad) do
                @test_throws ArgumentError MT.memory_budget_bytes()
            end
        end
        @test MT.native_gram_worker_extra_bytes(256) == 256 * MT._GRAM_SAT_MEMORY_BYTES
        withenv("SPACEAGORA_GRAM_SAT_MEMORY_MB" => "0") do
            @test MT.native_gram_worker_extra_bytes(256) == 0
        end
        @test MT.native_gram_worker_extra_bytes(0) == 0
    end

    @testset "snapshot is cached" begin
        a = MT.machine_topology()
        b = MT.machine_topology()
        @test a === b
    end

    @testset "overrides" begin
        try
            withenv("SPACEAGORA_PHYSICAL_CORES" => "8", "SPACEAGORA_CORE_BUDGET" => nothing) do
                t = MT.refresh_machine_topology!()
                @test t.physical_cores == 8
                @test t.usable_cores <= 8
            end
            withenv("SPACEAGORA_CORE_BUDGET" => "5", "SPACEAGORA_PHYSICAL_CORES" => nothing) do
                t = MT.refresh_machine_topology!()
                @test t.usable_cores == 5
                @test t.source === :override
            end
            # A budget override wins over every measured term, including one
            # that would have bound tighter.
            withenv("SPACEAGORA_CORE_BUDGET" => "3", "SPACEAGORA_PHYSICAL_CORES" => "2") do
                t = MT.refresh_machine_topology!()
                @test t.usable_cores == 3
            end
            for bad in ("0", "-2", "abc")
                withenv("SPACEAGORA_CORE_BUDGET" => bad) do
                    @test_throws ArgumentError MT.refresh_machine_topology!()
                end
            end
        finally
            # Leave the process holding a snapshot of the real machine: later
            # tests in the same session read this cache.
            MT.refresh_machine_topology!()
        end
    end

    @testset "the budget binds" begin
        # process_max_workers now sizes from usable cores rather than from
        # Sys.CPU_THREADS. On an SMT or affinity-constrained machine those
        # differ, and the old default asked for more parallelism than the
        # machine could deliver.
        t = MT.machine_topology()
        @test SpaceAGORA.ParallelProfiles.OuterRouteTuning().process_max_workers == t.usable_cores
        @test SpaceAGORA.ParallelProfiles.OuterRouteTuning().process_max_workers <= Sys.CPU_THREADS
        # Machine class is derived from the same budget, so a container or a
        # taskset cannot classify a machine by cores it may not touch.
        @test SpaceAGORA.ParallelProfiles._machine_parallel_class() in (:small, :medium, :large)
    end
end
