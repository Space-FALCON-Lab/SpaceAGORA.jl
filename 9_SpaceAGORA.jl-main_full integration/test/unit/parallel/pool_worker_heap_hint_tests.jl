using Test
using SpaceAGORA

# Default GC heap-size-hint for process-pool workers and the
# parallelization_performance harness's per-point worker. See
# docs/architecture/heap_contention.md for the diagnostic this is derived
# from and the peak-RSS verification. Pure-function tests only: no worker is
# actually spawned here (that is process_pool_probes.jl's job).

const PWH = SpaceAGORA.ParallelProcess

@testset "_default_pool_worker_heap_hint_bytes: clamped share-of-total-memory" begin
    # Int, not UInt64 -- _default_pool_worker_heap_hint_bytes takes
    # pool_size::Int, and Sys.total_memory() returns UInt64. Note this reads
    # whatever memory limit is visible to this process: under this study's
    # own systemd-run --user --scope -p MemoryMax=16G wrapper (required for
    # every Julia invocation here, see docs/architecture/heap_contention.md),
    # Julia's cgroup-aware Sys.total_memory() reports the 16 GiB scope limit,
    # not the host's physical total -- which is itself a demonstration that
    # the formula tracks whatever memory budget the process actually has.
    total = Int(Sys.total_memory())
    # Budget is `share * total`, not `total`: the point of the share
    # constant is that pool_size workers plus the coordinator, each hinted
    # at their equal split of `budget`, sum to at most `budget` -- strictly
    # less than the whole machine at share < 1 -- rather than to `total` the
    # way `total ÷ (pool_size + 1)` alone would.
    budget = total * PWH._POOL_WORKER_HEAP_HINT_SHARE_DEFAULT
    @test budget < total

    # Result is always within [floor, ceil], for every pool size tried.
    for pool_size in (0, 1, 2, 8, 32, 10_000)
        v = PWH._default_pool_worker_heap_hint_bytes(pool_size)
        @test PWH._POOL_WORKER_HEAP_HINT_FLOOR_BYTES <= v <= PWH._POOL_WORKER_HEAP_HINT_CEIL_BYTES
    end

    # A pool_size large enough drives the naive share below the floor;
    # clamped there. (2 GiB is small enough that this is reachable under any
    # plausible total, including the 16 GiB test-scope cap above.)
    huge_pool = cld(budget, PWH._POOL_WORKER_HEAP_HINT_FLOOR_BYTES) + 1000
    @test PWH._default_pool_worker_heap_hint_bytes(round(Int, huge_pool)) == PWH._POOL_WORKER_HEAP_HINT_FLOOR_BYTES

    # A pool_size small enough drives the naive share above the ceiling;
    # clamped there. Only reachable when this process's visible budget
    # (share * total) exceeds the 64 GiB ceiling -- not the case under this
    # study's own 16 GiB test cap (budget = 8 GiB there), so this only
    # exercises on a larger host/uncapped run.
    if budget > PWH._POOL_WORKER_HEAP_HINT_CEIL_BYTES
        @test PWH._default_pool_worker_heap_hint_bytes(0) == PWH._POOL_WORKER_HEAP_HINT_CEIL_BYTES
    end

    # A pool_size chosen to land the raw share inside [floor, ceil] agrees
    # with the formula exactly (no clamping engaged), and the aggregate
    # across that pool plus the coordinator is at most `budget`.
    mid_point = (PWH._POOL_WORKER_HEAP_HINT_FLOOR_BYTES + PWH._POOL_WORKER_HEAP_HINT_CEIL_BYTES) / 2
    mid_pool = max(1, round(Int, budget / mid_point) - 1)
    expected = floor(Int, budget / (mid_pool + 1))
    if PWH._POOL_WORKER_HEAP_HINT_FLOOR_BYTES <= expected <= PWH._POOL_WORKER_HEAP_HINT_CEIL_BYTES
        @test PWH._default_pool_worker_heap_hint_bytes(mid_pool) == expected
        @test expected * (mid_pool + 1) <= budget
    end

    # Monotone: a larger pool never gets a larger per-worker share.
    @test PWH._default_pool_worker_heap_hint_bytes(8) <= PWH._default_pool_worker_heap_hint_bytes(1)
end

@testset "_pool_worker_heap_hint_share: default / env override" begin
    withenv("SPACEAGORA_POOL_WORKER_HEAP_HINT_SHARE" => nothing) do
        @test PWH._pool_worker_heap_hint_share() == PWH._POOL_WORKER_HEAP_HINT_SHARE_DEFAULT
    end
    withenv("SPACEAGORA_POOL_WORKER_HEAP_HINT_SHARE" => "0.25") do
        @test PWH._pool_worker_heap_hint_share() == 0.25
        # Halving the share halves the unclamped per-worker hint.
        total = Int(Sys.total_memory())
        pool_size = max(1, cld(total, PWH._POOL_WORKER_HEAP_HINT_CEIL_BYTES) + 1)  # avoid the ceiling clamp
        quarter_share = PWH._default_pool_worker_heap_hint_bytes(pool_size)
        withenv("SPACEAGORA_POOL_WORKER_HEAP_HINT_SHARE" => "0.5") do
            half_share = PWH._default_pool_worker_heap_hint_bytes(pool_size)
            @test quarter_share <= half_share
        end
    end
    withenv("SPACEAGORA_POOL_WORKER_HEAP_HINT_SHARE" => "0") do
        @test PWH._pool_worker_heap_hint_share() == PWH._POOL_WORKER_HEAP_HINT_SHARE_DEFAULT  # <= 0 falls back
    end
    withenv("SPACEAGORA_POOL_WORKER_HEAP_HINT_SHARE" => "not-a-number") do
        @test PWH._pool_worker_heap_hint_share() == PWH._POOL_WORKER_HEAP_HINT_SHARE_DEFAULT  # unparseable falls back
    end
end

@testset "_format_heap_size_hint" begin
    @test PWH._format_heap_size_hint(2 * 1024^3) == "2.0G"
    @test occursin("G", PWH._format_heap_size_hint(1))
end

@testset "_pool_worker_heap_size_hint: default / off / explicit override" begin
    withenv("SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT" => nothing) do
        hint = PWH._pool_worker_heap_size_hint(4)
        @test hint isa String
        @test endswith(hint, "G")
    end
    withenv("SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT" => "off") do
        @test PWH._pool_worker_heap_size_hint(4) === nothing
    end
    withenv("SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT" => "OFF") do  # case-insensitive
        @test PWH._pool_worker_heap_size_hint(4) === nothing
    end
    withenv("SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT" => "3.5G") do
        @test PWH._pool_worker_heap_size_hint(4) == "3.5G"  # passed through verbatim
    end
end

@testset "_process_worker_exeflags: hint on by default, off removes it, explicit passes through" begin
    project = Base.active_project()

    withenv("SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT" => nothing) do
        cmd = PWH._process_worker_exeflags(project, 4)
        @test "--threads=1" in cmd.exec
        @test "--startup-file=no" in cmd.exec
        @test "--project=$(project)" in cmd.exec
        hint_flags = filter(a -> startswith(a, "--heap-size-hint="), cmd.exec)
        @test length(hint_flags) == 1
    end

    withenv("SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT" => "off") do
        cmd = PWH._process_worker_exeflags(project, 4)
        @test !any(a -> startswith(a, "--heap-size-hint="), cmd.exec)
    end

    withenv("SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT" => "7G") do
        cmd = PWH._process_worker_exeflags(project, 4)
        @test "--heap-size-hint=7G" in cmd.exec
    end

    # Default (no pool_size argument) still carries a hint -- backward
    # compatible with the 1-arg call sites in process_pool_probes.jl.
    withenv("SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT" => nothing) do
        cmd = PWH._process_worker_exeflags(project)
        @test any(a -> startswith(a, "--heap-size-hint="), cmd.exec)
    end
end

# --- Harness (parallelization_performance) per-point worker: the same
# derived default, its own override knob, and `off`. execution.jl's function
# signatures reference PPCCaseSpec (defined in cases.jl, not cli.jl) as type
# annotations, which Julia resolves at `include` time, so all three of
# cli.jl/modes.jl/cases.jl are needed even though ppc_worker_cmd and
# _ppc_worker_gc_flags themselves only touch PPCConfig -- same include chain
# benchmarks/studies/heap_contention/dump_state.jl and alloc_profile.jl use.
const _PWH_PPC_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "benchmarks", "studies", "parallelization_performance"))
include(joinpath(_PWH_PPC_DIR, "cli.jl"))
include(joinpath(_PWH_PPC_DIR, "modes.jl"))
include(joinpath(_PWH_PPC_DIR, "cases.jl"))
include(joinpath(_PWH_PPC_DIR, "execution.jl"))

@testset "harness per-point worker: heap-size-hint default / off / explicit" begin
    cfg = PPCConfig(profile="full", process_workers=4)

    withenv("SPACEAGORA_PPC_WORKER_HEAP_SIZE_HINT" => nothing) do
        flags = _ppc_worker_gc_flags(cfg)
        hint_flags = filter(a -> startswith(a, "--heap-size-hint="), flags)
        @test length(hint_flags) == 1
        expected_bytes = SpaceAGORA.ParallelProcess._default_pool_worker_heap_hint_bytes(max(cfg.process_workers, 1))
        expected = "--heap-size-hint=$(SpaceAGORA.ParallelProcess._format_heap_size_hint(expected_bytes))"
        @test only(hint_flags) == expected

        cmd = ppc_worker_cmd(cfg; case="independent_1sat_1hr", mode="serial", threads=1,
                              repeat=1, seed=1, mc_samples=1, outfile="/tmp/does_not_matter.csv", parity=false)
        @test any(a -> startswith(a, "--heap-size-hint="), cmd.exec)
    end

    withenv("SPACEAGORA_PPC_WORKER_HEAP_SIZE_HINT" => "off") do
        flags = _ppc_worker_gc_flags(cfg)
        @test !any(a -> startswith(a, "--heap-size-hint="), flags)
    end

    withenv("SPACEAGORA_PPC_WORKER_HEAP_SIZE_HINT" => "5G") do
        flags = _ppc_worker_gc_flags(cfg)
        @test "--heap-size-hint=5G" in flags
    end

    withenv("SPACEAGORA_PPC_WORKER_GCTHREADS" => "2", "SPACEAGORA_PPC_WORKER_HEAP_SIZE_HINT" => "off") do
        flags = _ppc_worker_gc_flags(cfg)
        @test "--gcthreads=2" in flags
        @test !any(a -> startswith(a, "--heap-size-hint="), flags)
    end
end

# _ppc_pool_worker_exeflags backs ppc_ensure_process_workers!, the harness's
# OWN Distributed-pool spawn path for outer_process/GRAM-live MC batches --
# independent of SpaceAGORA.ParallelProcess.ensure_process_workers! and of
# _ppc_worker_gc_flags above (which governs the per-point controller
# subprocess's own heap, not the sub-workers it spawns). It reuses
# SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT, the same knob
# _process_worker_exeflags reads, not a third env var.
@testset "harness process-pool sub-workers: heap-size-hint default / off / explicit" begin
    withenv("SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT" => nothing) do
        cmd = _ppc_pool_worker_exeflags(4)
        @test "--threads=1" in cmd.exec
        @test "--startup-file=no" in cmd.exec
        @test "--project=$(PPC_REPO_ROOT)" in cmd.exec
        hint_flags = filter(a -> startswith(a, "--heap-size-hint="), cmd.exec)
        @test length(hint_flags) == 1
        expected_bytes = SpaceAGORA.ParallelProcess._default_pool_worker_heap_hint_bytes(4)
        @test only(hint_flags) == "--heap-size-hint=$(SpaceAGORA.ParallelProcess._format_heap_size_hint(expected_bytes))"
    end

    withenv("SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT" => "off") do
        cmd = _ppc_pool_worker_exeflags(4)
        @test !any(a -> startswith(a, "--heap-size-hint="), cmd.exec)
    end

    withenv("SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT" => "9G") do
        cmd = _ppc_pool_worker_exeflags(4)
        @test "--heap-size-hint=9G" in cmd.exec
    end
end
