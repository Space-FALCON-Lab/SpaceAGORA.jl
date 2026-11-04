using Test
using SpaceAGORA

# Automatic machine calibration (`ensure_machine_constants!`) must not keep
# constants measured with fewer threads than the session using them, must not
# let several processes measure at once in one working directory, and every
# persisted policy file must be written through a unique temporary name.

const MC_PC = SpaceAGORA.SimulationModel.ParallelCost
const MC_RS = SpaceAGORA.RuntimeServices
const MC_SC = SpaceAGORA.SimulationCampaigns

function mc_constants(; threads::Int)
    return MC_PC.MachineConstants(
        simd_lane = MC_PC.RateCurve([4.0, 8.0], [0.24, 0.026]),
        coeff_touch = MC_PC.RateCurve([12.0, 16.0], [0.10, 0.17]),
        parallel_speedup = MC_PC.RateCurve([0.0, 1.0], [1.0, 2.0]),
        ns_per_scalar_item = 1.0,
        ns_per_queue_node = 0.05,
        dispatch_pool_ns_base = 3000.0,
        dispatch_pool_ns_per_worker = 1000.0,
        dispatch_batch_ns_base = 400.0,
        dispatch_batch_ns_per_worker = 25.0,
        ns_per_atomic = 20.0,
        reference_fma_ns = 20.0,
        reference_mem_ns = 1.0,
        fingerprint = "mc-test",
        schema_version = MC_PC.CALIBRATION_SCHEMA_VERSION,
        measured_threads = threads,
    )
end

mc_leftover_temps(dir) = filter(f -> endswith(f, ".tmp") || endswith(f, ".lock"), readdir(dir))

@testset "the measuring thread count is recorded and round-trips" begin
    mktempdir() do dir
        path = joinpath(dir, "constants.toml")
        MC_PC.save_machine_constants(mc_constants(threads = 7), path)
        @test MC_PC.load_machine_constants(path).measured_threads == 7
        # A file written before the count was recorded reads as 0 (unknown).
        text = replace(read(path, String), r"measured_threads = 7\n" => "")
        write(path, text)
        @test MC_PC.load_machine_constants(path).measured_threads == 0
        @test isempty(mc_leftover_temps(dir))
    end
end

@testset "constants measured with fewer threads are re-measured, never the reverse" begin
    mktempdir() do dir
        calls = Ref(0)
        calibrate(threads) = () -> (calls[] += 1; mc_constants(threads = threads))

        narrow = joinpath(dir, "narrow.toml")
        MC_PC.save_machine_constants(mc_constants(threads = 1), narrow)
        @test MC_PC.ensure_machine_constants!(path = narrow, calibrate = calibrate(8), threads = 8) === :calibrated
        @test calls[] == 1
        @test MC_PC.load_machine_constants(narrow).measured_threads == 8

        wide = joinpath(dir, "wide.toml")
        MC_PC.save_machine_constants(mc_constants(threads = 16), wide)
        @test MC_PC.ensure_machine_constants!(path = wide, calibrate = calibrate(4), threads = 4) === :present
        @test calls[] == 1
        @test MC_PC.load_machine_constants(wide).measured_threads == 16

        legacy = joinpath(dir, "legacy.toml")
        MC_PC.save_machine_constants(mc_constants(threads = 0), legacy)
        @test MC_PC.ensure_machine_constants!(path = legacy, calibrate = calibrate(64), threads = 64) === :present
        @test calls[] == 1

        # A campaign table survives the re-measurement.
        table = joinpath(dir, "table.toml")
        MC_PC.save_machine_constants(mc_constants(threads = 2), table)
        write(table, read(table, String) * "\n[campaign]\npool_cold_start_s = 1.5\n")
        @test MC_PC.ensure_machine_constants!(path = table, calibrate = calibrate(4), threads = 4) === :calibrated
        @test occursin("pool_cold_start_s", read(table, String))
        @test isempty(mc_leftover_temps(dir))
    end
end

@testset "a calibration another process holds is waited for, not repeated" begin
    mktempdir() do dir
        calls = Ref(0)
        calibrate = () -> (calls[] += 1; mc_constants(threads = 4))

        # The holder finishes while this process waits: its file is used.
        path = joinpath(dir, "shared.toml")
        write(path * ".lock", "otherhost 1\n")
        holder = @async begin
            sleep(0.5)
            MC_PC.save_machine_constants(mc_constants(threads = 4), path)
            rm(path * ".lock")
        end
        @test MC_PC.ensure_machine_constants!(path = path, calibrate = calibrate, threads = 4,
                                              lock_wait_s = 30.0) === :present
        wait(holder)
        @test calls[] == 0

        # The holder never finishes within the wait: nothing is measured here.
        stuck = joinpath(dir, "stuck.toml")
        write(stuck * ".lock", "otherhost 2\n")
        @test MC_PC.ensure_machine_constants!(path = stuck, calibrate = calibrate, threads = 4,
                                              lock_wait_s = 0.3) === :deferred
        @test calls[] == 0
        @test !isfile(stuck)
        @test isfile(stuck * ".lock")   # someone else's lock is left alone

        # A lock left by a process that died is taken over.
        dead = joinpath(dir, "dead.toml")
        write(dead * ".lock", "otherhost 3\n")
        @test MC_PC.ensure_machine_constants!(path = dead, calibrate = calibrate, threads = 4,
                                              lock_wait_s = 0.0, lock_stale_s = -1.0) === :calibrated
        @test calls[] == 1
        @test !isfile(dead * ".lock")
    end
end

@testset "concurrent tasks settle one path with one measurement" begin
    mktempdir() do dir
        path = joinpath(dir, "together.toml")
        calls = Threads.Atomic{Int}(0)
        calibrate = () -> (Threads.atomic_add!(calls, 1); sleep(0.2); mc_constants(threads = 2))
        outcomes = fetch.([Threads.@spawn MC_PC.ensure_machine_constants!(
                               path = path, calibrate = calibrate, threads = 2) for _ in 1:4])
        @test calls[] == 1
        @test count(==(:calibrated), outcomes) == 1
        @test count(==(:checked), outcomes) == 3
        @test isempty(mc_leftover_temps(dir))
    end
end

@testset "atomic writes use a unique temporary name and leave nothing behind" begin
    mktempdir() do dir
        path = joinpath(dir, "state.toml")
        seen = String[]
        seen_lock = ReentrantLock()
        tasks = [Threads.@spawn MC_RS.write_file_atomically(path) do io
                     lock(() -> push!(seen, io.name), seen_lock)
                     yield()
                     print(io, "value = ", i, "\n")
                 end for i in 1:16]
        foreach(wait, tasks)
        @test length(unique(seen)) == 16
        @test occursin(r"^value = \d+\n$", read(path, String))
        @test isempty(mc_leftover_temps(dir))

        # A failing writer removes its temporary file and leaves the target.
        @test_throws ErrorException MC_RS.write_file_atomically(io -> error("writer failed"), path)
        @test occursin(r"^value = \d+\n$", read(path, String))
        @test isempty(mc_leftover_temps(dir))
    end
end

@testset "campaign corrections are written the same way" begin
    mktempdir() do dir
        path = joinpath(dir, "corrections.toml")
        c = MC_SC.CampaignCorrections("mc-test", "t")
        c.campaigns = 3
        foreach(wait, [Threads.@spawn MC_SC.save_campaign_corrections(c, path) for _ in 1:8])
        @test MC_SC.load_campaign_corrections(path; fingerprint = "mc-test", code_token = "t").campaigns == 3
        @test isempty(mc_leftover_temps(dir))
    end
end
