using Test
using Distributed

const _WP_REPO_ROOT = isdefined(Main, :REPO_ROOT) ? Main.REPO_ROOT : normpath(joinpath(@__DIR__, "..", ".."))

if !isdefined(@__MODULE__, :ParallelProcess)
    include(joinpath(_WP_REPO_ROOT, "src", "parallel", "process", "parallel_process.jl"))
    using .ParallelProcess
end

@testset "Process Pool Coverage Debt Probes" begin
    @testset "pure helpers" begin
        cmd = ParallelProcess._process_worker_exeflags(_WP_REPO_ROOT)
        @test cmd isa Cmd
        @test "--threads=1" in cmd.exec
        @test "--startup-file=no" in cmd.exec
        @test "--project=$(_WP_REPO_ROOT)" in cmd.exec

        fresh_pool = ProcessPool(_WP_REPO_ROOT)
        @test fresh_pool.workers == Int[]
        @test fresh_pool.project_path == _WP_REPO_ROOT

        # shutdown on an already-empty pool is a documented no-op, not an error.
        @test shutdown_process_pool!(fresh_pool) === nothing
        @test isempty(fresh_pool.workers)
    end

    @testset "campaign_process_pool singleton" begin
        p1 = campaign_process_pool()
        p2 = campaign_process_pool()
        @test p1 === p2
        @test p1 isa ProcessPool
    end

    @testset "ensure_process_workers! real bootstrap" begin
        pool = ProcessPool(_WP_REPO_ROOT)
        try
            ids1 = ensure_process_workers!(pool, 1)
            @test length(ids1) == 1
            @test ids1[1] in Distributed.procs()
            # Bootstrapped with `using SpaceAGORA` and --threads=1, regardless of
            # whether the best-effort GRAMSuite/SPICE furnishing steps succeeded
            # in this (GRAM-less) test environment.
            @test remotecall_fetch(() -> isdefined(Main, :SpaceAGORA), ids1[1])
            @test remotecall_fetch(() -> Threads.nthreads(), ids1[1]) == 1

            # Requesting <= the current worker count spawns nothing new.
            ids_same = ensure_process_workers!(pool, 1)
            @test ids_same == ids1

            # Growing the pool only spawns the shortfall; existing workers are reused.
            ids2 = ensure_process_workers!(pool, 2; warmup_fn=() -> 1 + 1)
            @test length(ids2) == 2
            @test ids2[1] == ids1[1]
            @test ids2[2] in Distributed.procs()
            @test remotecall_fetch(() -> isdefined(Main, :SpaceAGORA), ids2[2])

            # A warmup_fn that throws is best-effort: the new worker still joins the pool.
            ids3 = ensure_process_workers!(pool, 3; warmup_fn=() -> error("probe-injected warmup failure"))
            @test length(ids3) == 3
            @test ids3[1:2] == ids2
        finally
            shutdown_process_pool!(pool)
            @test isempty(pool.workers)
        end
    end
end
