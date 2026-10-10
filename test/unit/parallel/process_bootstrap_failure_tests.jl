using Test
import Distributed

# Exercise the actual pool/ensure implementation without starting OS workers or
# modifying Distributed's methods. Unqualified rmprocs resolves to this fixture's
# own function; all other production code is included unchanged.
module ProcessBootstrapFailureFixture
import Distributed
using Distributed: CachingPool

include(joinpath(@__DIR__, "..", "..", "..", "src", "parallel", "process", "worker_pool.jl"))

const STATE = Ref{Any}()
const FAILURE_MESSAGE = "injected mandatory bootstrap failure"

function reset!(; mode=:blocked_pair, cleanup_error=false)
    STATE[] = (
        mode=Ref(mode), cleanup_error=cleanup_error, next_id=Ref(101),
        spawned=Vector{Int}[], spawn_requests=Tuple{Int,Int}[],
        started=Int[], finished=Int[], removed=Vector{Int}[],
        finished_at_cleanup=Vector{Int}[], lock_held_at_cleanup=Bool[],
        sibling_started=Channel{Nothing}(1), failure_started=Channel{Nothing}(1),
        release=Channel{Nothing}(1), released=Ref(false), pool=Ref{Any}(nothing),
    )
end

function _spawn_process_workers(n::Int, project_path::AbstractString, pool_size::Int)
    state = STATE[]
    ids = collect(state.next_id[]:(state.next_id[] + n - 1))
    state.next_id[] += n
    push!(state.spawned, copy(ids))
    push!(state.spawn_requests, (n, pool_size))
    return ids
end

function _bootstrap_process_worker!(worker::Int, project_path::String)::Nothing
    state = STATE[]
    push!(state.started, worker)
    if state.mode[] === :blocked_pair
        if worker == 101
            timedwait(() -> isready(state.sibling_started), 5.0; pollint=0.005) === :ok ||
                error("Bootstrap fixture: sibling did not start")
            put!(state.failure_started, nothing)
            error(FAILURE_MESSAGE)
        elseif worker == 102
            put!(state.sibling_started, nothing)
            timedwait(() -> isready(state.release), 10.0; pollint=0.005) === :ok ||
                error("Bootstrap fixture: sibling was not released")
            take!(state.release)
        end
    elseif state.mode[] === :immediate_failure
        error(FAILURE_MESSAGE)
    end
    push!(state.finished, worker)
    return nothing
end

function rmprocs(ids::AbstractVector{<:Integer}; kwargs...)
    state = STATE[]
    push!(state.removed, collect(Int, ids))
    push!(state.finished_at_cleanup, copy(state.finished))
    push!(state.lock_held_at_cleanup, islocked(state.pool[].lock))
    state.cleanup_error && error("injected removal failure")
    return nothing
end

function release_sibling!()
    state = STATE[]
    if !state.released[]
        state.released[] = true
        put!(state.release, nothing)
    end
end

function capture_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end
end

@testset "process bootstrap failure ownership" begin
    fixture = ProcessBootstrapFailureFixture
    prior_processes = Distributed.procs()
    fixture.reset!()
    state = fixture.STATE[]
    pool = fixture.ProcessPool(Base.active_project())
    push!(pool.workers, 41) # An existing member must survive a failed addition.
    state.pool[] = pool
    startup = @async fixture.ensure_process_workers!(pool, 3)
    try
        failure_seen = timedwait(() -> isready(state.failure_started), 5.0; pollint=0.005)
        @test failure_seen === :ok
        # The first bootstrap has thrown, while the sibling cannot finish until
        # released below. The old asyncmap returned exceptionally at this point.
        @test timedwait(() -> istaskdone(startup), 0.5; pollint=0.005) === :timed_out
        got_lock = trylock(pool.lock)
        got_lock && unlock(pool.lock)
        @test !got_lock
        @test sort(state.started) == [101, 102]
        @test isempty(state.finished)
        @test isempty(state.removed)
        @test pool.workers == [41]
    finally
        fixture.release_sibling!()
    end

    try
        @test timedwait(() -> istaskdone(startup), 5.0; pollint=0.005) === :ok
        err = fixture.capture_error(() -> fetch(startup))
        @test err !== nothing
        @test occursin(fixture.FAILURE_MESSAGE, sprint(showerror, err))
        @test state.finished == [102]
        @test state.removed == [[101, 102]]
        @test state.finished_at_cleanup == [[102]]
        @test state.lock_held_at_cleanup == [true]
        @test pool.workers == [41]
        got_lock = trylock(pool.lock)
        got_lock && unlock(pool.lock)
        @test got_lock

        # A retry must provision only the shortfall and register its own new
        # members, without reusing failed IDs or removing the existing member.
        state.mode[] = :success
        @test fixture.ensure_process_workers!(pool, 3) == [41, 103, 104]
        @test state.spawned == [[101, 102], [103, 104]]
        @test state.spawn_requests == [(2, 3), (2, 3)]
        @test sort(state.finished) == [102, 103, 104]
        @test state.removed == [[101, 102]]
        @test fixture.ensure_process_workers!(pool, 3) == [41, 103, 104]
        @test length(state.spawned) == 2
        @test sort(state.started) == [101, 102, 103, 104]
    finally
        fixture.release_sibling!()
        # Every injected wait is bounded. Ensure no gated task survives a test
        # failure, including when this file is checked against the old source.
        @test timedwait(() -> istaskdone(startup), 10.0; pollint=0.005) === :ok
        @test timedwait(() -> 102 in state.finished, 2.0; pollint=0.005) === :ok
    end
    @test Distributed.procs() == prior_processes
end

@testset "bootstrap cleanup error preserves original failure" begin
    fixture = ProcessBootstrapFailureFixture
    prior_processes = Distributed.procs()
    fixture.reset!(mode=:immediate_failure, cleanup_error=true)
    state = fixture.STATE[]
    pool = fixture.ProcessPool(Base.active_project())
    push!(pool.workers, 41)
    state.pool[] = pool
    err = @test_logs (:warn, "Could not remove new process workers after bootstrap failed.") begin
        fixture.capture_error(() -> fixture.ensure_process_workers!(pool, 2))
    end
    @test err !== nothing
    @test occursin(fixture.FAILURE_MESSAGE, sprint(showerror, err))
    @test !occursin("injected removal failure", sprint(showerror, err))
    @test state.removed == [[101]]
    @test state.lock_held_at_cleanup == [true]
    @test pool.workers == [41]
    got_lock = trylock(pool.lock)
    got_lock && unlock(pool.lock)
    @test got_lock
    @test Distributed.procs() == prior_processes
end
