using Test
using Distributed
using Base.ScopedValues: with
using SpaceAGORA

# PROCESS_WARMUP overrides the warm-up the campaign runner gives new process
# workers, and ensure_process_workers! bootstraps and warms new workers
# concurrently. A blank line once detached ensure_process_workers!'s docstring,
# so both docstrings are checked as attached.

const PW = SpaceAGORA.ParallelProcess

@testset "process warm-up" begin
    # Read the module's own doc table: the top-level @doc forwarding in
    # SpaceAGORA.jl would mask a detached docstring from Docs.hasdoc.
    doctext(s) = join((join(d.text) for d in values(get(Base.Docs.meta(PW), Base.Docs.Binding(PW, s), Base.Docs.MultiDoc()).docs)), "\n")
    @test occursin("warmup_fn=nothing", doctext(:ensure_process_workers!))
    @test occursin("Scoped override", doctext(:PROCESS_WARMUP))

    @test PW.PROCESS_WARMUP[] === nothing
    with(PW.PROCESS_WARMUP => false) do
        @test PW.PROCESS_WARMUP[] === false
    end
    @test PW.PROCESS_WARMUP[] === nothing

    prior_workers = sort(workers())
    pool = PW.ProcessPool(Base.active_project())
    try
        ids = PW.ensure_process_workers!(pool, 2; warmup_fn=() -> myid())
        @test length(ids) == 2 && !any(in(prior_workers), ids)
        @test all(w -> remotecall_fetch(Core.eval, w, Main, :(isdefined(Main, :SpaceAGORA))), ids)
        # Existing workers are kept and not warmed again.
        @test PW.ensure_process_workers!(pool, 2; warmup_fn=() -> error("rewarmed")) == ids
    finally
        PW.shutdown_process_pool!(pool)
    end
    @test sort(workers()) == prior_workers
end
