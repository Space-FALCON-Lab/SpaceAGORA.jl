using Test
using Distributed
using SpaceAGORA

# SPACEAGORA_PROCESS_WORKER_PRELOAD names packages a new pool worker loads right
# after SpaceAGORA (the paper harness uses it for its precompile workload). A
# name that does not resolve, or is not precompiled for the worker's
# environment, is skipped with a warning rather than failing the bootstrap.

const PWPL = SpaceAGORA.ParallelProcess

@testset "process worker preload" begin
    withenv("SPACEAGORA_PROCESS_WORKER_PRELOAD" => nothing) do
        @test isempty(PWPL._process_worker_preload_names())
    end
    withenv("SPACEAGORA_PROCESS_WORKER_PRELOAD" => " Statistics, ,NotAPackageForPreloadTests ") do
        @test PWPL._process_worker_preload_names() == ["Statistics", "NotAPackageForPreloadTests"]
    end

    prior_workers = sort(workers())
    w = only(PWPL._spawn_process_workers(1, Base.active_project()))
    try
        stats_id = Base.identify_package("Statistics")
        loaded_on_worker() = remotecall_fetch(Core.eval, w, Main,
            :(Base.root_module_exists($(stats_id))))
        @test !loaded_on_worker()
        withenv("SPACEAGORA_PROCESS_WORKER_PRELOAD" => "Statistics,NotAPackageForPreloadTests") do
            @test_logs (:warn, r"NotAPackageForPreloadTests") match_mode=:any PWPL._preload_worker_packages!(w)
        end
        @test loaded_on_worker()
    finally
        rmprocs(w; waitfor=60)
    end
    @test sort(workers()) == prior_workers
end
