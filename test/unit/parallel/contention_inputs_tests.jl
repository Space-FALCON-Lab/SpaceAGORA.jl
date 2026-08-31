using Test
using SpaceAGORA

const PCC = SpaceAGORA.SimulationModel.ParallelCost

@testset "contention inputs" begin
    @testset "an unmeasured workload is not a well-behaved one" begin
        # The default must read as "no evidence", never as "allocates nothing,
        # holds no lock, scales perfectly" -- those are the same numbers and
        # opposite conclusions.
        ci = PCC.ContentionInputs()
        @test !ci.valid
        @test PCC.alloc_bytes_per_second(ci) == 0.0
        @test PCC.lock_duty_cycle(ci) == 0.0
        @test PCC.lock_wait_hold_ratio(ci) == 0.0
    end

    @testset "allocation rate" begin
        ci = PCC.ContentionInputs(
            alloc_bytes_per_pass = 2_000.0,
            pass_ns = 1_000_000.0,
            valid = true,
        )
        # 2 kB per millisecond = 2 MB/s.
        @test PCC.alloc_bytes_per_second(ci) ≈ 2.0e6
        # A zero-length window yields no rate rather than Inf.
        @test PCC.alloc_bytes_per_second(
            PCC.ContentionInputs(alloc_bytes_per_pass = 2_000.0, pass_ns = 0.0, valid = true)
        ) == 0.0
        # A genuinely allocation-free workload is distinguishable from an
        # unmeasured one only by `valid`, so that distinction has to hold.
        free = PCC.ContentionInputs(alloc_bytes_per_pass = 0.0, pass_ns = 1.0e6, valid = true)
        @test free.valid
        @test PCC.alloc_bytes_per_second(free) == 0.0
    end

    @testset "lock duty cycle and the width signal" begin
        ci = PCC.ContentionInputs(
            pass_ns = 1_000_000.0,
            lock_hold_ns = 250_000.0,
            lock_wait_ns = 750_000.0,
            lock_acquisitions = 12,
            valid = true,
        )
        @test PCC.lock_duty_cycle(ci) ≈ 0.25
        # Spread over four workers' worth of time the same hold is a quarter as
        # much of the budget. 1/rho is the speedup ceiling, so this is the term
        # that decides how wide the inner split may usefully go.
        @test PCC.lock_duty_cycle(ci, 4) ≈ 0.0625
        # Three times as much waiting as holding: already past what the lock
        # admits.
        @test PCC.lock_wait_hold_ratio(ci) ≈ 3.0
        # Occupancy can never exceed its own window.
        @test PCC.lock_duty_cycle(
            PCC.ContentionInputs(pass_ns = 1.0, lock_hold_ns = 1.0e9, valid = true)
        ) == 1.0
    end

    @testset "footprint is carried, not re-measured" begin
        # Re-measuring what WorkCounts already counts is how the two come to
        # disagree; this field exists to be filled from the counts.
        ci = PCC.ContentionInputs(
            workspace_bytes_per_sat = 4096.0, n_active_sats = 1024, valid = true,
        )
        @test ci.workspace_bytes_per_sat * ci.n_active_sats == 4_194_304.0
    end
end
