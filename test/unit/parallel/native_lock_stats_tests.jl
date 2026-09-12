using Test
using SpaceAGORA

const RS = SpaceAGORA.RuntimeServices

@testset "native lock occupancy" begin
    @testset "site indexing" begin
        for (i, site) in enumerate(RS.NATIVE_LOCK_SITES)
            @test RS.native_lock_site_index(site) == i
        end
        # An unknown site is charged to :other rather than rejected: a new call
        # site must not be able to fail a run by naming itself.
        @test RS.native_lock_site_index(:not_a_site) ==
              RS.native_lock_site_index(:other)
    end

    @testset "mutual exclusion is unchanged" begin
        # Different sites are views of ONE lock. If they ever stop excluding
        # each other, the CSPICE call-trace corruption that merged GRAM_LOCK
        # and SPICE_LOCK comes straight back.
        a = RS.tracked_lock(:gram_density)
        b = RS.tracked_lock(:spice_frame)
        lock(a)
        try
            @test islocked(b)
            @test islocked(RS.SPICE_LOCK)
        finally
            unlock(a)
        end
        @test !islocked(b)
    end

    @testset "records one acquisition and a hold" begin
        RS.reset_native_lock_stats!()
        RS.with_native_lock(:spice_frame) do
            # Enough work that the hold is unambiguously above timer noise.
            s = 0.0
            for i in 1:100_000
                s += sqrt(i)
            end
            s
        end
        snap = RS.native_lock_stats_snapshot()
        @test snap.acquisitions == 1
        @test snap.hold_ns > 0
        @test snap.wait_ns == 0
        @test snap.sites.spice_frame.acquisitions == 1
        @test snap.sites.spice_frame.hold_ns == snap.hold_ns
        @test snap.sites.gram_density.acquisitions == 0
    end

    @testset "reentrant acquisition counts once" begin
        # A ReentrantLock lets the owning task re-acquire. Charging each nesting
        # level would inflate occupancy by the call depth, and the depth is a
        # property of the call graph rather than of the contention.
        RS.reset_native_lock_stats!()
        RS.with_native_lock(:gram_density) do
            RS.with_native_lock(:spice_body) do
                RS.with_native_lock(:spice_body) do
                    nothing
                end
            end
        end
        snap = RS.native_lock_stats_snapshot()
        @test snap.acquisitions == 1
        # Attributed to the OUTERMOST site: that is the span other tasks are
        # actually excluded for.
        @test snap.sites.gram_density.acquisitions == 1
        @test snap.sites.spice_body.acquisitions == 0
    end

    @testset "wait time is recorded under contention" begin
        if Threads.nthreads() < 2
            @test_skip "needs at least two threads"
        else
            RS.reset_native_lock_stats!()
            holder = Threads.@spawn RS.with_native_lock(:gram_density) do
                sleep(0.25)
            end
            sleep(0.05)
            waiter = Threads.@spawn RS.with_native_lock(:gram_density) do
                nothing
            end
            wait(holder)
            wait(waiter)
            snap = RS.native_lock_stats_snapshot()
            @test snap.acquisitions == 2
            @test snap.contended == 1
            @test snap.wait_ns > 0
            # The waiter blocked for most of the holder's span, so the ratio has
            # to be a substantial fraction of one. This is the quantity that
            # reports over-width, so a zero here would be a silent failure of
            # the whole mechanism.
            @test snap.wait_hold_ratio > 0.3
        end
    end

    @testset "duty cycle" begin
        RS.reset_native_lock_stats!()
        RS.with_native_lock(:gram_density) do
            sleep(0.05)
        end
        snap = RS.native_lock_stats_snapshot()
        rho = RS.native_lock_duty_cycle(snap.hold_ns, 1)
        @test rho ≈ 1.0
        # Same hold spread over four workers' worth of time is a quarter of it.
        @test RS.native_lock_duty_cycle(snap.hold_ns, 4) ≈ 0.25 atol = 1e-9
        # A degenerate window reports "no constraint observed", never Inf/NaN.
        @test RS.native_lock_duty_cycle(0, 1) == 0.0
        @test RS.native_lock_duty_cycle(-1, 8) == 0.0
        # Occupancy can never exceed the window it is measured against.
        @test RS.native_lock_duty_cycle(1, 1) <= 1.0
    end

    @testset "disabling stops recording but still locks" begin
        previous = RS.set_native_lock_stats!(false)
        try
            RS.reset_native_lock_stats!()
            entered = false
            RS.with_native_lock(:gram_density) do
                entered = true
                @test islocked(RS.SPICE_LOCK)
            end
            @test entered
            @test !islocked(RS.SPICE_LOCK)
            snap = RS.native_lock_stats_snapshot()
            @test snap.acquisitions == 0
            @test !snap.enabled
        finally
            RS.set_native_lock_stats!(previous)
            RS.reset_native_lock_stats!()
        end
    end

    @testset "reset clears every site" begin
        RS.reset_native_lock_stats!()
        RS.with_native_lock(:gram_cache) do; nothing; end
        RS.with_native_lock(:spice_ephemeris) do; nothing; end
        @test RS.native_lock_stats_snapshot().acquisitions == 2
        RS.reset_native_lock_stats!()
        snap = RS.native_lock_stats_snapshot()
        @test snap.acquisitions == 0
        @test snap.hold_ns == 0
        @test snap.wait_ns == 0
        @test all(getproperty(snap.sites, s).acquisitions == 0 for s in RS.NATIVE_LOCK_SITES)
    end
end
