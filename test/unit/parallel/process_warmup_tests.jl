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
    prior_processes = procs()
    pool = PW.ProcessPool(Base.active_project())
    events = RemoteChannel(() -> Channel{Tuple{Symbol,Int}}(8))
    release = RemoteChannel(() -> Channel{Nothing}(2))
    startup_task = nothing
    released = false
    release_workers!() = if !released
        released = true
        put!(release, nothing)
        put!(release, nothing)
    end
    try
        warmup = let events=events, release=release
            () -> begin
                put!(events, (:start, myid()))
                timedwait(() -> isready(release), 55.0; pollint=0.02) === :ok ||
                    error("Process warm-up test barrier timed out")
                take!(release)
                put!(events, (:finish, myid()))
            end
        end
        startup_task = @async PW.ensure_process_workers!(pool, 2; warmup_fn=warmup)
        starts = Tuple{Symbol,Int}[]
        try
            # Bootstrap may include cold package loading. Once one warm-up has
            # started, both must reach the barrier before either is released.
            for timeout in (600.0, 30.0)
                status = timedwait(() -> isready(events) || istaskdone(startup_task),
                                   timeout; pollint=0.02)
                status === :ok && isready(events) || break
                push!(starts, take!(events))
            end
            @test length(starts) == 2
            @test all(e -> e[1] === :start, starts)
            @test length(unique(e[2] for e in starts)) == 2
        finally
            # The second-arrival timeout is shorter than the worker timeout,
            # so a serial implementation fails without stranding its workers.
            release_workers!()
        end
        timedwait(() -> istaskdone(startup_task), 120.0; pollint=0.02) === :ok ||
            error("Process startup did not finish after releasing warm-up")
        ids = fetch(startup_task)
        finishes = Tuple{Symbol,Int}[]
        while isready(events)
            push!(finishes, take!(events))
        end
        @test length(ids) == 2 && !any(in(prior_workers), ids)
        @test sort([e[2] for e in starts]) == sort(ids)
        @test sort(finishes) == sort([(:finish, w) for w in ids])
        @test all(w -> remotecall_fetch(Core.eval, w, Main, :(isdefined(Main, :SpaceAGORA))), ids)
        # A thrown warm-up error alone is not evidence: production catches it.
        # Record any rewarm before throwing, then require no such event.
        unexpected_rewarm = let events=events
            () -> begin
                put!(events, (:unexpected_rewarm, myid()))
                error("rewarmed")
            end
        end
        @test PW.ensure_process_workers!(pool, 2; warmup_fn=unexpected_rewarm) == ids
        @test !isready(events)
    finally
        try
            release_workers!()
            if startup_task !== nothing && !istaskdone(startup_task)
                # ensure holds pool.lock until it finishes. Remove this test's
                # processes first to unblock remote calls before taking that
                # lock in shutdown. Also covers workers not yet registered.
                spawned = setdiff(procs(), prior_processes)
                isempty(spawned) || rmprocs(spawned; waitfor=30)
                timedwait(() -> istaskdone(startup_task), 30.0; pollint=0.02) === :ok ||
                    error("Process startup task did not stop during cleanup")
            end
            PW.shutdown_process_pool!(pool)
            spawned = setdiff(procs(), prior_processes)
            isempty(spawned) || rmprocs(spawned; waitfor=30)
        finally
            close(events)
            close(release)
        end
    end
    @test sort(workers()) == prior_workers
end
