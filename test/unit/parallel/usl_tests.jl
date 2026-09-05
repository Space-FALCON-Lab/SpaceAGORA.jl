using Test
using SpaceAGORA

const PU = SpaceAGORA.SimulationModel.ParallelCost

@testset "universal scalability law" begin
    @testset "shape" begin
        # beta = 0 is Amdahl: saturating, never turning over.
        @test PU.usl_speedup(0.0, 0.0, 8) ≈ 8.0
        @test PU.usl_speedup(0.1, 0.0, 12) < 12.0
        @test PU.usl_peak_workers(0.1, 0.0) == Inf
        # beta > 0 turns the curve back down, which is the behaviour a
        # saturating model cannot express and the measured RHS curves require.
        peak = PU.usl_peak_workers(0.05, 0.01)
        @test 9.0 < peak < 10.0
        @test PU.usl_speedup(0.05, 0.01, 30) < PU.usl_speedup(0.05, 0.01, Int(round(peak)))
        # Never below one: declining to split is always available.
        @test PU.usl_speedup(0.9, 0.9, 64) >= 1.0
        @test PU.usl_speedup(0.5, 0.5, 1) == 1.0
    end

    @testset "fit inverts the model" begin
        for (a, b) in ((0.05, 0.01), (0.2, 0.0), (0.0, 0.005))
            ks = [1, 2, 4, 8, 12]
            sp = [PU.usl_speedup(a, b, k) for k in ks]
            fa, fb = PU._fit_usl(ks, sp)
            @test fa ≈ a atol = 1e-6
            @test fb ≈ b atol = 1e-6
        end
        # Degenerate ladders must not throw or produce a negative beta, which
        # would predict speedup rising without bound.
        @test PU._fit_usl(Int[], Float64[]) == (0.0, 0.0)
        @test PU._fit_usl([1], [1.0]) == (0.0, 0.0)
        a, b = PU._fit_usl([2], [1.5])
        @test b == 0.0 && a >= 0.0
        # A ladder that scales better than linearly (noise) clamps rather than
        # inverting the sign of either parameter.
        a2, b2 = PU._fit_usl([2, 4, 8], [2.2, 4.5, 9.0])
        @test a2 >= 0.0 && b2 >= 0.0
    end

    @testset "beta does not depend on the width it is evaluated at" begin
        # The parameters describe the workload. A beta that moved with the
        # candidate width made k* = sqrt((1-alpha)/beta) report a different
        # optimum depending on where it was evaluated -- an optimum that moves
        # when you look at it from a different place is not an optimum.
        mc = PU.MachineConstants(
            simd_lane = PU.RateCurve([3.0], [1.0]),
            coeff_touch = PU.RateCurve([3.0], [1.0]),
            parallel_speedup = PU.RateCurve([0.0], [1.0]),
            ns_per_scalar_item = 1.0, ns_per_queue_node = 1.0,
            dispatch_pool_ns_base = 0.0, dispatch_pool_ns_per_worker = 0.0,
            dispatch_batch_ns_base = 0.0, dispatch_batch_ns_per_worker = 0.0,
            ns_per_atomic = 1.0, reference_fma_ns = 1.0, reference_mem_ns = 1.0,
            usl_alpha_base = 0.05, usl_beta_arith = 0.001,
            usl_beta_alloc = 0.01, usl_beta_bw = 1.0,
            llc_bytes = 32.0 * 1024^2,
        )
        ci = PU.ContentionInputs(
            alloc_bytes_per_pass = 1024.0 * 1024, pass_ns = 1.0e6,
            workspace_bytes_per_sat = 22_000.0, n_active_sats = 256, valid = true)
        a4, b4 = PU.workload_usl_parameters(mc, ci; n_active_sats = 256)
        @test b4 > mc.usl_beta_arith
        @test PU.usl_peak_workers(a4, b4) > 0.0

        # An invalid probe must fall back to the machine baseline, not to zero
        # contention: "unmeasured" and "scales perfectly" are opposite claims.
        a0, b0 = PU.workload_usl_parameters(mc, PU.ContentionInputs(); n_active_sats = 256)
        @test a0 == mc.usl_alpha_base
        @test b0 == mc.usl_beta_arith
        @test PU.workload_usl_parameters(mc, nothing; n_active_sats = 256) == (a0, b0)
    end

    @testset "lock occupancy enters alpha, not beta" begin
        mc = PU.MachineConstants(
            simd_lane = PU.RateCurve([3.0], [1.0]),
            coeff_touch = PU.RateCurve([3.0], [1.0]),
            parallel_speedup = PU.RateCurve([0.0], [1.0]),
            ns_per_scalar_item = 1.0, ns_per_queue_node = 1.0,
            dispatch_pool_ns_base = 0.0, dispatch_pool_ns_per_worker = 0.0,
            dispatch_batch_ns_base = 0.0, dispatch_batch_ns_per_worker = 0.0,
            ns_per_atomic = 1.0, reference_fma_ns = 1.0, reference_mem_ns = 1.0,
            usl_alpha_base = 0.02, usl_beta_arith = 0.0,
            usl_beta_alloc = 0.0, usl_beta_bw = 0.0, llc_bytes = 32.0 * 1024^2,
        )
        # rho = 0.5 of the pass held in the native lock.
        locked = PU.ContentionInputs(pass_ns = 1.0e6, lock_hold_ns = 5.0e5, valid = true)
        a, b = PU.workload_usl_parameters(mc, locked; n_active_sats = 64)
        @test b == 0.0                      # a lock is not a contention term
        @test a ≈ 0.52 atol = 1e-9          # it is a serial fraction
        # ... and a serial fraction caps speedup at 1/alpha however wide.
        @test PU.usl_speedup(a, b, 64) < 1.0 / a + 1.0
        @test PU.usl_speedup(a, b, 1024) < 1.0 / a + 1.0
    end
end
