using Test
using SpaceAGORA

const RSC = SpaceAGORA.RuntimeServices
const SEC = SpaceAGORA.SimulationEngine

@testset "native-lock width cap" begin
    @testset "ceiling from occupancy" begin
        RSC.reset_native_lock_stats!()
        # No lock activity means no constraint -- not a ceiling of one.
        @test RSC.lock_width_ceiling() == typemax(Int)
        o = RSC.native_lock_occupancy()
        @test o.rho == 0.0
        @test o.window_s >= 0.0

        RSC.reset_native_lock_stats!()
        RSC.with_native_lock(:gram_density) do
            sleep(0.05)
        end
        sleep(0.05)
        occ = RSC.native_lock_occupancy()
        # Roughly half the window was spent holding, so the ceiling is small.
        @test 0.2 < occ.rho < 0.8
        c = RSC.lock_width_ceiling()
        @test 1 <= c <= 6
        # The ceiling is 1/rho: a lock caps speedup there however wide the split.
        @test c == max(1, ceil(Int, 1 / occ.rho))
    end

    @testset "floor_rho refuses to build a ceiling out of noise" begin
        RSC.reset_native_lock_stats!()
        RSC.with_native_lock(:spice_frame) do; nothing; end
        sleep(0.05)
        # A handful of nanoseconds of hold over a 50 ms window implies a ceiling
        # far above any real core budget; reporting "no constraint" is both
        # cheaper and more honest than a large number derived from one
        # acquisition.
        @test RSC.lock_width_ceiling() == typemax(Int)
        # ... but a low floor lets the same window produce one.
        @test RSC.lock_width_ceiling(floor_rho = 1e-12) < typemax(Int)
    end

    @testset "the clamp only removes width, never adds it" begin
        params = (; shared_buffers = (; rhs_width_ceiling = Ref{Int}(3)))
        wide = SEC._make_calib_flat_plan(12, :static)
        @test SEC._clamp_plan_to_lock_ceiling(wide, params).allotment == 3
        # Already under the ceiling: untouched, and the same object's other
        # fields survive.
        narrow = SEC._make_calib_flat_plan(2, :static)
        clamped = SEC._clamp_plan_to_lock_ceiling(narrow, params)
        @test clamped.allotment == 2
        @test clamped.scheduler === narrow.scheduler
        @test clamped.mode === narrow.mode
        # Unconstrained ceiling changes nothing.
        open_params = (; shared_buffers = (; rhs_width_ceiling = Ref{Int}(typemax(Int))))
        @test SEC._clamp_plan_to_lock_ceiling(wide, open_params).allotment == 12
        # No params, no clamp -- callers outside a solve are unaffected.
        @test SEC._clamp_plan_to_lock_ceiling(wide, nothing).allotment == 12
    end

    @testset "satellite_batch is exempt" begin
        # It takes its width from Polyester's own pool and honours neither
        # `allotment` nor the inner thread budget, so clamping the field would
        # change what the plan reports without changing what it runs.
        params = (; shared_buffers = (; rhs_width_ceiling = Ref{Int}(1)))
        batch = SEC._make_calib_satellite_batch_plan()
        @test SEC._clamp_plan_to_lock_ceiling(batch, params) === batch
    end
end
