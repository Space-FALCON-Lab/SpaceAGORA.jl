using Test
using Distributed
using SpaceAGORA

# A vendored environment prepended by the coordinator must not shadow the pool
# project on a fresh worker. This pins the reproduced mixed-resolution mechanism;
# the exact PR139 CI extension failure still needs current-head CI confirmation.
# Native-free: the worker is bare and is removed afterwards.

const PWLP = SpaceAGORA.ParallelProcess
const PWLP_PROJECT = Base.active_project()
const PWLP_REPO = dirname(dirname(pathof(SpaceAGORA)))
const PWLP_VENDORED = joinpath(PWLP_REPO, "data", "GRAMSuite.jl")

@testset "process workers start on a project-first load path" begin
    saved = copy(LOAD_PATH)
    prior_workers = sort(workers())
    temporary_extra = !isdir(PWLP_VENDORED)
    extra = temporary_extra ? mktempdir() : PWLP_VENDORED
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
        temporary_extra && rm(extra; recursive=true, force=true)
    end
    @test sort(workers()) == prior_workers
end
