using Test
using Distributed
using SpaceAGORA

# Regression for the PR #139 coverage failure (job 105631685037): Distributed
# serialises the coordinator's LOAD_PATH into each worker's JULIA_LOAD_PATH, so
# a coordinator that had prepended the vendored GRAMSuite environment (the
# native probes and ensure_gramsuite_loaded! do this while GRAMSuite is only a
# weak dependency) handed every pool worker a stack in which that environment
# shadowed the pool's own project, and the worker's `using SpaceAGORA` failed
# while precompiling a package extension it resolved through the wrong
# manifest. Pool workers must resolve packages from the pool project first;
# other coordinator entries stay behind it as fallbacks. Native-free: the
# spawned worker is bare (no SpaceAGORA bootstrap) and is removed afterwards.

const PWLP = SpaceAGORA.ParallelProcess
const PWLP_PROJECT = Base.active_project()
const PWLP_REPO = dirname(dirname(pathof(SpaceAGORA)))
const PWLP_VENDORED = joinpath(PWLP_REPO, "data", "GRAMSuite.jl")

@testset "process workers start on a project-first load path" begin
    saved = copy(LOAD_PATH)
    prior_workers = sort(workers())
    extra = isdir(PWLP_VENDORED) ? PWLP_VENDORED : mktempdir()
    try
        pushfirst!(LOAD_PATH, extra)
        pathsep = Sys.iswindows() ? ";" : ":"
        entries = split(PWLP._process_worker_load_path(), pathsep)
        @test entries[1] == "@"
        @test count(==("@"), entries) == 1
        @test entries[2:end] == filter(!=("@"), LOAD_PATH)
        @test extra in entries[2:end]

        # The actual spawn used by ensure_process_workers!, on a bare worker.
        w = only(PWLP._spawn_process_workers(1, PWLP_PROJECT))
        try
            seen = remotecall_fetch(Core.eval, w, Main, quote
                stack = Base.load_path()
                sa = Base.identify_package_env("StaticArrays")
                (first=first(stack), extra_kept=$(extra) in dirname.(stack) || $(extra) in stack,
                 project=Base.active_project(), sa_env=sa === nothing ? nothing : sa[2])
            end)
            @test seen.project == PWLP_PROJECT
            @test seen.first == PWLP_PROJECT
            @test seen.extra_kept
            # StaticArrays is a direct dependency of the project and of the
            # vendored GRAMSuite project: the worker must find it in the project.
            @test seen.sa_env == PWLP_PROJECT
        finally
            rmprocs(w; waitfor=60)
        end
    finally
        empty!(LOAD_PATH)
        append!(LOAD_PATH, saved)
    end
    @test sort(workers()) == prior_workers
end
