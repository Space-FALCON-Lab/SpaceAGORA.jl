using Test
using SpaceAGORA

const PC = SpaceAGORA.SimulationModel.ParallelCost

@testset "Sign test matches the binomial tail" begin
    # Hand-checked two-sided Binomial(n, 1/2) tail probabilities.
    @test PC._sign_test_two_sided(5, 5) ≈ 2 * (1 / 32) atol = 1e-12
    @test PC._sign_test_two_sided(4, 5) ≈ 2 * (6 / 32) atol = 1e-12
    @test PC._sign_test_two_sided(3, 5) ≈ 1.0 atol = 1e-12   # clamped
    @test PC._sign_test_two_sided(10, 10) ≈ 2 * (1 / 1024) atol = 1e-12
    # Symmetry: 2-of-5 is as extreme as 3-of-5.
    @test PC._sign_test_two_sided(2, 5) == PC._sign_test_two_sided(3, 5)
    @test PC._sign_test_two_sided(0, 0) == 1.0
    # 12 of 15 clears 0.05; 10 of 15 does not.
    @test PC._sign_test_two_sided(12, 15) <= 0.05
    @test PC._sign_test_two_sided(10, 15) > 0.05
end

# Deterministic work of a controllable size, reading from a runtime-filled
# buffer so the kernel cannot be constant-folded.
#
# This is not incidental. An earlier version of this file spun on a literal
# (`_spin(512)`), which Julia folded to a constant before timing, and the
# measurement came back at a fraction of a nanosecond -- reporting the cost of
# returning a constant rather than of doing the work. A 4x-work comparison then
# read as 137x. The `timed_min` sink defeats dead-code elimination but not
# constant folding, so a pure kernel must be given opaque inputs.
const _SPIN_BUF = [1.0 + i * 1e-9 for i in 1:8192]

function _spin(n::Int)::Float64
    acc = 0.0
    @inbounds @simd for i in 1:n
        acc = muladd(_SPIN_BUF[i], 1.0000001, acc)
    end
    return acc
end

@testset "Noise floor separates cheap from elided" begin
    floor_ns = PC.timing_noise_floor()
    @test floor_ns > 0.0
    # Real work must clear the harness floor by a wide margin, or the
    # measurement is not measuring anything.
    @test PC.timed_min(() -> _spin(4096)) > 10 * floor_ns

    # And the failure mode is detectable: a folded literal kernel sits at the
    # floor. This is the regression guard for the bug described above.
    folded() = 0.0
    @test PC.timed_min(folded) < 5 * floor_ns
end

@testset "timed_min recovers relative cost" begin
    t1 = PC.timed_min(() -> _spin(512))
    t4 = PC.timed_min(() -> _spin(2048))
    @test t1 > 0.0
    @test t4 > 0.0
    # Ratios are far more stable than absolute durations, so the model's
    # consumers use them and so does this test. Wide bounds: the point is that
    # 4x the work reads as several times the cost, not that it reads as exactly
    # 4.00x -- SIMD width and loop overhead make the true factor machine-specific.
    @test 2.0 < t4 / t1 < 8.0
    # Per-lane cost should be near-constant across sizes once the kernel is
    # opaque -- this is what a correctly-scaling estimator looks like.
    @test 0.6 < (t4 / 2048) / (t1 / 512) < 1.7
end

@testset "timed_min is not fooled by one-sided interference" begin
    # Interference is additive and one-sided: a disturbed sample is longer,
    # never shorter. So the minimum over windows should track the clean cost
    # while the mean is dragged upward. Deterministic, not random: every third
    # call does an order of magnitude more work.
    clean = PC.timed_min(() -> _spin(1024); k = 9)

    counter = Ref(0)
    disturbed = PC.timed_min(k = 9) do
        counter[] += 1
        counter[] % 3 == 0 ? _spin(10_240) : _spin(1024)
    end

    mean_inflation = disturbed / clean
    # The disturbed workload averages ~4x the clean one (1/3 of calls at 10x).
    # The minimum cannot fully escape it here, because interference recurs
    # inside every window -- which is the honest result and the reason the
    # estimator is paired with interleaving rather than relied on alone.
    # What must hold is that it never reads *below* clean, i.e. the estimator
    # is not optimistic.
    @test mean_inflation > 0.9
end

@testset "paired_compare identifies a large true difference" begin
    cmp = PC.paired_compare(() -> _spin(256), () -> _spin(4096); pairs = 15)
    @test cmp.wins_a > cmp.wins_b
    @test cmp.significant
    @test cmp.median_ratio < 1.0          # A is cheaper, so ta/tb < 1
end

@testset "paired_compare does not manufacture a winner from noise" begin
    # Same work on both sides. A clean interleave should split the pairs; a
    # systematic order bias would show as a sweep for whichever side is
    # measured first, which is exactly the failure this check is here to catch.
    cmp = PC.paired_compare(() -> _spin(1024), () -> _spin(1024); pairs = 15)
    @test cmp.wins_a > 0
    @test cmp.wins_b > 0
    @test 0.5 < cmp.median_ratio < 2.0
end

@testset "Reference kernels are positive and reproducible" begin
    for f in (PC.reference_kernel_ns, PC.reference_memory_kernel_ns)
        a = f()
        b = f()
        @test a > 0.0
        @test isfinite(a)
        # Reproducibility is the property the staleness canary depends on: if
        # two back-to-back readings disagree, the canary cannot distinguish
        # "machine changed" from "measurement noise".
        @test 0.5 < a / b < 2.0
    end
end

@testset "The two references measure different resources" begin
    # If they tracked the same bottleneck there would be no point carrying both.
    # A cache-line touch that misses L1 costs materially more than one lane of
    # L1-resident FMA; if this ever stops holding, the memory reference is not
    # doing its job and the stride-bound terms are effectively unnormalised.
    fma = PC.reference_kernel_ns()
    mem = PC.reference_memory_kernel_ns()
    @test mem > fma
end
