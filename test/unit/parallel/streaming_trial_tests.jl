using Test
using SpaceAGORA

const PT = SpaceAGORA.SimulationModel.ParallelCost

# Drive a trial from a cost function, the way a solve would: the trial says
# which arm, the caller reports what it cost.
function _drive!(trial, cost)
    n = 0
    while PT.trial_active(trial) && n < 10_000
        a = PT.next_arm(trial)
        PT.observe!(trial, cost(a, n))
        n += 1
    end
    return n
end

@testset "streaming paired trial" begin
    @testset "schedule" begin
        t = PT.StreamingPairedTrial(3; rounds = 4, warmup_rounds = 0)
        seen = Int[]
        while PT.trial_active(t)
            push!(seen, PT.next_arm(t))
            PT.observe!(t, 1.0)
        end
        # Every arm once per round.
        for chunk in Iterators.partition(seen, 3)
            @test sort(collect(chunk)) == [1, 2, 3]
        end
        # And the order rotates, so no arm is stuck in first position.
        @test seen[1:3] != seen[4:6]
    end

    @testset "warm-up rounds are discarded" begin
        # A switch of execution width re-warms caches; the first observation
        # after one measures the switch, not the arm.
        t = PT.StreamingPairedTrial(2; rounds = 5, warmup_rounds = 2)
        n = 0
        while PT.trial_active(t)
            PT.observe!(t, n < 4 ? 1000.0 : 10.0)   # first two rounds are slow
            n += 1
        end
        @test PT.trial_rounds(t) == 5
        @test all(all(==(10.0), o) for o in t.obs)
    end

    @testset "detects a real difference" begin
        t = PT.StreamingPairedTrial(2; rounds = 15, warmup_rounds = 1)
        _drive!(t, (arm, i) -> arm == 2 ? 50.0 + (i % 3) : 100.0 + (i % 3))
        v = PT.trial_verdict(t)
        @test v.arm == 2
        @test v.significant
        @test v.median_ratio ≈ 2.0 atol = 0.1
    end

    @testset "declines to call a tie" begin
        # Arms that differ only by noise must leave the incumbent in place. This
        # is the property that makes the trial safe to run by default.
        t = PT.StreamingPairedTrial(2; rounds = 15, warmup_rounds = 1)
        _drive!(t, (arm, i) -> 100.0 + (arm == 1 ? (i % 5) : ((i + 2) % 5)))
        v = PT.trial_verdict(t)
        @test !v.significant
    end

    @testset "a lucky median does not win without the rounds" begin
        # One wild round can drag a median below the incumbent's while the arm
        # loses nearly every round. Significance is judged on round wins, so
        # that arm must not be selected.
        t = PT.StreamingPairedTrial(2; rounds = 15, warmup_rounds = 0)
        n = 0
        while PT.trial_active(t)
            arm = PT.next_arm(t)
            PT.observe!(t, arm == 1 ? 100.0 : (n == 3 ? 1.0 : 101.0))
            n += 1
        end
        v = PT.trial_verdict(t)
        @test !v.significant
    end

    @testset "drift cancels within a round" begin
        # A cost that climbs steadily with time -- thermal drift, a co-tenant
        # arriving -- must not create a winner. Blocked sampling would report
        # the first-measured arm as faster by the whole drift.
        t = PT.StreamingPairedTrial(2; rounds = 15, warmup_rounds = 1)
        _drive!(t, (arm, i) -> 100.0 + 5.0 * i)
        v = PT.trial_verdict(t)
        @test !v.significant
    end

    @testset "speedups feed the USL fit" begin
        widths = [1, 2, 4, 8]
        alpha, beta = 0.08, 0.004
        t = PT.StreamingPairedTrial(length(widths); rounds = 6, warmup_rounds = 1)
        _drive!(t, (arm, _i) -> 1000.0 / PT.usl_speedup(alpha, beta, widths[arm]))
        sp = PT.trial_speedups(t)
        @test sp[1] ≈ 1.0 atol = 1e-9
        fa, fb = PT._fit_usl(widths, sp)
        @test fa ≈ alpha atol = 1e-3
        @test fb ≈ beta atol = 1e-4
    end

    @testset "no observations is not a verdict" begin
        t = PT.StreamingPairedTrial(2; rounds = 3, warmup_rounds = 0)
        v = PT.trial_verdict(t)
        @test v.arm == 1 && !v.significant && v.rounds == 0
        @test_throws ArgumentError PT.StreamingPairedTrial(1)
    end
end
