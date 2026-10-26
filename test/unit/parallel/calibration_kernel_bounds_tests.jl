using Test
using SpaceAGORA

# The calibration kernels read their tables under `@inbounds`. In a normal
# session a read past the table is not an error: it returns whatever memory
# follows the table, and now and then it faults. These tests keep every
# calibration read inside its buffer.

const CKB_PC = SpaceAGORA.SimulationModel.ParallelCost
const CKB_COST_DIR = joinpath(pkgdir(SpaceAGORA), "src", "orchestration", "parallel", "cost")

# The child loads no package. Under `--check-bounds=yes` Julia cannot use the
# package caches compiled for the default flags, and compiling SpaceAGORA's
# dependencies again would take far longer than the test. So it includes the
# calibration sources themselves (timing, kernels, ladders) into one module,
# with a serial stand-in for the persistent worker pool, and runs the ladders
# with the calibration's own arguments.
const CKB_CHILD_SCRIPT = raw"""
const COST_DIR = ARGS[1]

module CalibrationLadders
    module ParallelPolicy
        # Serial stand-in for the persistent worker pool. Every ladder gives
        # each of its workers one item, so item `i` runs as worker `i`.
        function threaded_foreach_worker_persistent(f, ::Symbol, num_items::Int, allotment::Int; kwargs...)
            for i in 1:num_items
                f(mod1(i, max(1, allotment)), i)
            end
            return nothing
        end
    end

    # Named in method signatures only; no ladder builds one.
    struct MachineConstants end

    include(joinpath(Main.COST_DIR, "robust_timing.jl"))
    include(joinpath(Main.COST_DIR, "calibration_kernels.jl"))
    include(joinpath(Main.COST_DIR, "machine_calibration.jl"))
end

const CL = CalibrationLadders
CL._calibrate_usl(1, false)
println("CALIBRATION_LADDERS_IN_BOUNDS threads=", Threads.nthreads(),
        " check_bounds=", Base.JLOptions().check_bounds)
"""

@testset "calibration ladders stay inside their buffers under --check-bounds=yes" begin
    dir = mktempdir()
    script = joinpath(dir, "calibration_ladders.jl")
    write(script, CKB_CHILD_SCRIPT)
    # Four threads, so the ladders also run their two- and four-worker rungs.
    cmd = `$(Base.julia_cmd()) --startup-file=no --check-bounds=yes --threads=4 $(script) $(CKB_COST_DIR)`
    output = IOBuffer()
    ok = success(pipeline(setenv(ignorestatus(cmd), ENV; dir = dir); stdout = output, stderr = output))
    text = String(take!(output))
    ok || println(text)
    @test ok
    # check_bounds 1 is `--check-bounds=yes`.
    @test occursin("CALIBRATION_LADDERS_IN_BOUNDS threads=4 check_bounds=1", text)
end

@testset "the batch kernel walks its table again when it has more touches than columns" begin
    tab = CKB_PC.calib_table(6, 5)
    row, spt = 3, 2
    for n_touch in (1, 5, 6, 13)
        # The same touches on a table that repeats the columns, read without
        # wrapping, give the same bits.
        wide = repeat(tab, 1, cld(n_touch, size(tab, 2)))
        ws_wrapped = zeros(Float64, 8, 4)
        ws_wide = zeros(Float64, 8, 4)
        a = CKB_PC.calib_batch_kernel!(ws_wrapped, tab, row, n_touch, spt)
        b = CKB_PC.calib_batch_kernel!(ws_wide, wide, row, n_touch, spt)
        @test a === b
        @test ws_wrapped == ws_wide
    end
end

@testset "the batch kernel refuses a row outside the table and an empty workspace" begin
    tab = CKB_PC.calib_table(6, 5)
    ws = zeros(Float64, 8, 4)
    @test_throws BoundsError CKB_PC.calib_batch_kernel!(ws, tab, 7, 4, 2)
    @test_throws BoundsError CKB_PC.calib_batch_kernel!(ws, tab, 0, 4, 2)
    @test_throws BoundsError CKB_PC.calib_batch_kernel!(zeros(Float64, 0, 4), tab, 3, 4, 2)
end
