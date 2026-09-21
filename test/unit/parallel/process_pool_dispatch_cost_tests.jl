using Test
using SpaceAGORA

# The cross-campaign closure cache ships OFF (see `_pool_dispatch_cache_enabled`);
# these tests exercise it ON, and restore the caller's setting afterwards.
const _PPDC_CACHE_ENV_BEFORE = get(ENV, "SPACEAGORA_POOL_DISPATCH_CACHE", nothing)
ENV["SPACEAGORA_POOL_DISPATCH_CACHE"] = "1"

# The process dispatch's fixed cost.
#
# `CachingPool` keys its worker-side cache on the IDENTITY of the function it
# is handed, so a dispatch that builds a fresh pool and a fresh wrapper closure
# is a guaranteed cache miss on every worker: the campaign's closure is
# serialized to each of them, a channel is created there to hold it, and the
# `clear!` afterwards tears both down. That is paid per dispatch, not per
# sample, which is why it is invisible on a long campaign and most of a short
# one. These tests are about the cache that stops it being paid twice for the
# same campaign function, and about the bound on what it retains.
#
# Everything here runs against the `worker_runner` seam, in one process, with
# no `Distributed` workers spawned: the identity rules are what this file is
# for and they are the same whichever side of the socket the sample runs on.

const SCamp = SpaceAGORA.SimulationCampaigns

_dispatch_pool() = SpaceAGORA.ProcessPool(Base.active_project())
_in_process_worker = (f, index, seed) -> SCamp._run_monte_carlo_sample(f, index, seed)

# Two separately built closures over the same captured value are `===` in
# Julia -- a closure is an immutable struct, and egality on those is by field.
# So is `objectid`, which is what the `IdDict` behind a `CachingPool` keys on:
# the cache's notion of "the same function" and this cache's `!==` test are the
# same notion, which is why the tests below assert on the `CachingPool` object
# (mutable, so identity is identity) wherever the question is whether a cache
# was REUSED rather than rebuilt identically.
@testset "the dispatch runner is cached per campaign function" begin
    pool = _dispatch_pool()
    f = x -> x * 2
    g = x -> x * 3

    first_use = SCamp._acquire_dispatch_runner(pool, f, [2, 3], true)
    @test first_use.shared
    @test first_use.pool === pool.dispatch_pool
    SCamp._release_dispatch_runner(pool, first_use)

    # Same f, same workers: the same callable through the same pool, so the
    # workers' copies of it are still the ones they hold.
    second_use = SCamp._acquire_dispatch_runner(pool, f, [2, 3], true)
    @test second_use.run === first_use.run
    @test second_use.pool === first_use.pool
    @test second_use.shared
    SCamp._release_dispatch_runner(pool, second_use)

    # A different campaign function drops the cached one rather than keeping
    # both: what the cache retains is one campaign's closure, never a set.
    other = SCamp._acquire_dispatch_runner(pool, g, [2, 3], true)
    @test other.run !== first_use.run
    @test other.pool !== first_use.pool
    @test pool.dispatch_f === g
    SCamp._release_dispatch_runner(pool, other)

    # A different worker set is a different cache too: a `CachingPool` hands
    # out the workers it was built with.
    widened = SCamp._acquire_dispatch_runner(pool, g, [2, 3, 4], true)
    @test widened.pool !== other.pool
    @test pool.dispatch_workers == [2, 3, 4]
    SCamp._release_dispatch_runner(pool, widened)

    # A dispatch whose worker class is stood in for needs no `CachingPool`,
    # and still caches the callable both classes run.
    seamed = SCamp._acquire_dispatch_runner(pool, g, [2, 3, 4], false)
    @test seamed.pool === nothing
    @test seamed.shared
    SCamp._release_dispatch_runner(pool, seamed)
end

@testset "a dispatch beside another one gets its own cache, not the shared one" begin
    pool = _dispatch_pool()
    f = x -> x * 2
    held = SCamp._acquire_dispatch_runner(pool, f, [2, 3], true)
    @test held.shared
    # One CachingPool cannot serve two dispatches: its feeders take workers
    # from one channel. The second dispatch pays the old cost instead of
    # racing the first one for its workers.
    beside = SCamp._acquire_dispatch_runner(pool, f, [2, 3], true)
    @test !beside.shared
    @test beside.pool !== held.pool
    # Releasing the private one leaves the shared cache alone.
    SCamp._release_dispatch_runner(pool, beside)
    @test pool.dispatch_pool === held.pool
    SCamp._release_dispatch_runner(pool, held)
    again = SCamp._acquire_dispatch_runner(pool, f, [2, 3], true)
    @test again.pool === held.pool
    SCamp._release_dispatch_runner(pool, again)
end

@testset "SPACEAGORA_POOL_DISPATCH_CACHE=0 restores the per-dispatch closure" begin
    pool = _dispatch_pool()
    f = x -> x * 2
    withenv("SPACEAGORA_POOL_DISPATCH_CACHE" => "0") do
        a = SCamp._acquire_dispatch_runner(pool, f, [2, 3], true)
        SCamp._release_dispatch_runner(pool, a)
        b = SCamp._acquire_dispatch_runner(pool, f, [2, 3], true)
        SCamp._release_dispatch_runner(pool, b)
        @test !a.shared && !b.shared
        @test a.pool !== b.pool
        @test pool.dispatch_runner === nothing      # nothing was retained
        @test pool.dispatch_pool === nothing
    end
end

@testset "two dispatches of one campaign function share one cached closure" begin
    pool = _dispatch_pool()
    seeds = collect(101:120)
    spec = SCamp.MonteCarloSpec(seeds = seeds, threads = 1)
    f = x -> x * 2
    out1 = SCamp._run_monte_carlo_mixed(f, seeds, spec, [2, 3], 2;
        worker_runner = _in_process_worker, cache_pool = pool)
    runner = pool.dispatch_runner
    @test runner !== nothing
    @test pool.dispatch_f === f
    @test !pool.dispatch_busy                       # released at the end of the dispatch

    out2 = SCamp._run_monte_carlo_mixed(f, seeds, spec, [2, 3], 2;
        worker_runner = _in_process_worker, cache_pool = pool)
    @test pool.dispatch_runner === runner           # the second dispatch reused it

    # The campaign is unchanged by the reuse: same samples, same order, same
    # values, whichever consumer ran them.
    for out in (out1, out2)
        @test [s.index for s in out] == collect(1:length(seeds))
        @test [s.seed for s in out] == seeds
        @test [s.value for s in out] == seeds .* 2
        @test all(s -> s.success, out)
    end

    # A campaign with a different function drops the previous one.
    g = x -> x * 3
    out3 = SCamp._run_monte_carlo_mixed(g, seeds, spec, [2, 3], 2;
        worker_runner = _in_process_worker, cache_pool = pool)
    @test pool.dispatch_runner !== runner
    @test pool.dispatch_f === g
    @test [s.value for s in out3] == seeds .* 3
end

@testset "the cached closure is dropped with the pool that holds it" begin
    pool = _dispatch_pool()
    f = x -> x * 2
    entry = SCamp._acquire_dispatch_runner(pool, f, [2, 3], true)
    SCamp._release_dispatch_runner(pool, entry)
    @test pool.dispatch_runner !== nothing
    @test pool.dispatch_pool !== nothing
    # Shutting the pool down is the other bound on retention (the first being
    # the next campaign with a different f); a pool with no workers still has a
    # closure to let go of.
    SpaceAGORA.shutdown_process_pool!(pool)
    @test pool.dispatch_runner === nothing
    @test pool.dispatch_f === nothing
    @test pool.dispatch_pool === nothing
    @test isempty(pool.dispatch_workers)

    # And the explicit drop, which is what clears the workers' copies.
    reacquired = SCamp._acquire_dispatch_runner(pool, f, [2, 3], true)
    SCamp._release_dispatch_runner(pool, reacquired)
    SCamp._drop_dispatch_cache!(pool)
    @test pool.dispatch_runner === nothing
    @test pool.dispatch_pool === nothing
end

@testset "ordinal_sink names the consumer, not just its class" begin
    seeds = collect(1:24)
    spec = SCamp.MonteCarloSpec(seeds = seeds, threads = 1)
    classes = fill(:unset, length(seeds))
    ordinals = zeros(Int, length(seeds))
    out = SCamp._run_monte_carlo_mixed(x -> x, seeds, spec, [2, 3], 2;
        class_sink = classes, ordinal_sink = ordinals,
        worker_runner = _in_process_worker)
    @test length(out) == length(seeds)
    @test all(c -> c in (:worker, :local), classes)
    # Two consumers per class, so every ordinal is one of them and every
    # sample has one.
    @test all(o -> 1 <= o <= 2, ordinals)
    single = SCamp._run_monte_carlo_mixed(x -> x, seeds, spec, [2], 1;
        ordinal_sink = zeros(Int, length(seeds)),
        worker_runner = _in_process_worker)
    @test length(single) == length(seeds)
    @test_throws ArgumentError SCamp._run_monte_carlo_mixed(
        x -> x, seeds, spec, [2], 1; ordinal_sink = zeros(Int, 3),
        worker_runner = _in_process_worker)
end

if _PPDC_CACHE_ENV_BEFORE === nothing
    delete!(ENV, "SPACEAGORA_POOL_DISPATCH_CACHE")
else
    ENV["SPACEAGORA_POOL_DISPATCH_CACHE"] = _PPDC_CACHE_ENV_BEFORE
end
