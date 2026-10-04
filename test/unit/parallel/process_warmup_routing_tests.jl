# Exercise the production routing and pool methods with local worker transport.
# The existing process_warmup_tests.jl covers actual Distributed workers; this
# fixture makes same-sample overlap deterministic without another cold pool.
module ProcessWarmupRoutingTests
using Test
using SpaceAGORA
using Base.ScopedValues: with

const Campaigns = SpaceAGORA.SimulationCampaigns
const ProcessPool = SpaceAGORA.ParallelProcess.ProcessPool
const MonteCarloSpec = Campaigns.MonteCarloSpec
const MonteCarloResult = Campaigns.MonteCarloResult
const MonteCarloSampleResult = Campaigns.MonteCarloSampleResult
const DispatchAdmission = Campaigns.DispatchAdmission
const PROCESS_WARMUP = SpaceAGORA.PROCESS_WARMUP
const POOL = Ref(ProcessPool(Base.active_project()))
const BOOTSTRAPPED = Int[]
const WARMED = Int[]

# Named injection seams: worker allocation/bootstrap, remote transport, and
# the later real dispatch. No methods in SpaceAGORA or Distributed are replaced.
campaign_process_pool() = POOL[]
_spawn_process_workers(n::Int, project::AbstractString, size::Int) = collect(2:(n + 1))
_bootstrap_process_worker!(worker::Int, project::String) = push!(BOOTSTRAPPED, worker)
function remotecall_fetch(f, worker::Int)
    push!(WARMED, worker)
    return task_local_storage(() -> f(), :warmup_routing_worker, worker)
end
_maybe_probe_pool_dispatch(args...) = nothing
_dispatch_trace_enabled() = false
_collect_on_workers(workers) = nothing
function _run_monte_carlo_process(f, seeds, spec, workers; kwargs...)
    return MonteCarloSampleResult[
        MonteCarloSampleResult(index=i, seed=seed, success=true, elapsed_s=0.0,
                              value=f(seed)) for (i, seed) in enumerate(seeds)
    ]
end

# Parse the complete production function, including its keyword arguments and
# body. In particular, warm-up selection and the concurrency flag are not
# restated in the fixture.
function load_method(path, name)
    source = read(path, String)
    marker = findfirst("function $(name)(", source)
    marker === nothing && error("Missing production method: $(name)")
    method, _ = Meta.parse(source, first(marker))
    method isa Expr && method.head === :function || error("Expected function: $(name)")
    Core.eval(@__MODULE__, method)
end
const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
load_method(joinpath(REPO, "src", "parallel", "process", "worker_pool.jl"), :ensure_process_workers!)
load_method(joinpath(REPO, "src", "simulation", "campaigns", "adaptive_routing.jl"), :_run_campaign_with_route_env)

function reset_pool!()
    POOL[] = ProcessPool(Base.active_project())
    empty!(BOOTSTRAPPED)
    empty!(WARMED)
end
function run_route(f)
    spec = MonteCarloSpec(seeds=[11, 22], threads=2)
    return _run_campaign_with_route_env(f, spec, (route=:process, local_slots=0))
end

# A callable object, not just a Function, and a bounded rendezvous. Serializing
# this explicitly selected warm-up must fail the overlap assertion, not hang.
struct ConcurrentWarmup
    starts::Vector{Int}
    finishes::Vector{Int}
    timed_out::Vector{Int}
end
function (warm::ConcurrentWarmup)()
    worker = task_local_storage(:warmup_routing_worker)
    push!(warm.starts, worker)
    if timedwait(() -> length(warm.starts) == 2, 2.0; pollint=0.01) !== :ok
        push!(warm.timed_out, worker)
    end
    push!(warm.finishes, worker)
    return nothing
end

@testset "implicit warm-up preserves per-sample output exclusion" begin
    reset_pool!()
    previous = PROCESS_WARMUP[]
    mktempdir() do output
        calls = Tuple{Int, Int}[]
        collisions = Int[]
        writer = seed -> begin
            worker = get(task_local_storage(), :warmup_routing_worker, 0)
            push!(calls, (worker, seed))
            lockdir = joinpath(output, "member_$(seed).lock")
            # Record before throwing: production deliberately catches warm-up
            # errors, so successful route completion alone cannot test safety.
            if isdir(lockdir)
                push!(collisions, seed)
                error("Concurrent writer for member $(seed)")
            end
            mkdir(lockdir)
            try
                # Give an incorrectly concurrent second warm-up a bounded
                # opportunity to collide. The correct serial path waits once.
                if worker == 2
                    timedwait(() -> length(calls) >= 2, 1.0; pollint=0.01)
                end
                write(joinpath(output, "member_$(seed).txt"), string(seed))
            finally
                rm(lockdir)
            end
            return seed
        end
        result = with(PROCESS_WARMUP => nothing) do
            run_route(writer)
        end
        @test isempty(collisions)
        @test calls == [(2, 11), (3, 11), (0, 11), (0, 22)]
        @test sort(BOOTSTRAPPED) == [2, 3]
        @test WARMED == [2, 3]
        @test result.route === :process
        @test getproperty.(result.samples, :value) == [11, 22]
        @test all(seed -> read(joinpath(output, "member_$(seed).txt"), String) == string(seed), (11, 22))
        # Reusing the same pool must not repeat either warm-up or bootstrap.
        empty!(calls)
        with(PROCESS_WARMUP => nothing) do
            run_route(writer)
        end
        @test calls == [(0, 11), (0, 22)]
        @test WARMED == [2, 3]
        @test sort(BOOTSTRAPPED) == [2, 3]
    end
    @test PROCESS_WARMUP[] === previous
end

@testset "explicit callable overlaps and false skips warm-up" begin
    previous = PROCESS_WARMUP[]
    reset_pool!()
    custom = ConcurrentWarmup(Int[], Int[], Int[])
    calls = Int[]
    sample = seed -> (push!(calls, seed); seed)
    result = with(PROCESS_WARMUP => custom) do
        run_route(sample)
    end
    @test sort(custom.starts) == [2, 3]
    @test sort(custom.finishes) == [2, 3]
    @test isempty(custom.timed_out)
    @test sort(WARMED) == [2, 3]
    @test calls == [11, 22]
    @test getproperty.(result.samples, :value) == [11, 22]
    @test PROCESS_WARMUP[] === previous

    reset_pool!()
    empty!(calls)
    with(PROCESS_WARMUP => custom) do
        with(PROCESS_WARMUP => false) do
            run_route(sample)
        end
        @test PROCESS_WARMUP[] === custom
    end
    @test isempty(WARMED)
    @test sort(BOOTSTRAPPED) == [2, 3]
    @test calls == [11, 22]
    @test length(custom.starts) == 2
    @test PROCESS_WARMUP[] === previous
end
end
