using Test
using SpaceAGORA

const PC = SpaceAGORA.SimulationModel.ParallelCost

@testset "RateCurve interpolates in log2 and clamps" begin
    c = PC.RateCurve([log2(1024.0), log2(4096.0), log2(16384.0)], [0.10, 0.20, 0.30])

    # Knots reproduce exactly.
    @test PC.rate_at(c, 1024) ≈ 0.10
    @test PC.rate_at(c, 4096) ≈ 0.20
    @test PC.rate_at(c, 16384) ≈ 0.30

    # Midpoint in log2, not in linear size: 2048 is halfway between 1024 and
    # 4096 on a log2 axis, so it must read 0.15 rather than the linear-interp
    # value. Cache behaviour is a function of log size, which is why the curve
    # is stored that way.
    @test PC.rate_at(c, 2048) ≈ 0.15

    # Clamped outside the measured range. Extrapolating a curve whose shape is
    # only known where it was sampled would invent cache behaviour.
    @test PC.rate_at(c, 1) ≈ 0.10
    @test PC.rate_at(c, 10^9) ≈ 0.30

    # Degenerate single-knot curve is constant.
    @test PC.rate_at(PC.RateCurve([0.0], [0.42]), 12345) ≈ 0.42
end

@testset "RateCurve rejects malformed knots" begin
    @test_throws ArgumentError PC.RateCurve([1.0, 2.0], [0.1])
    @test_throws ArgumentError PC.RateCurve(Float64[], Float64[])
    @test_throws ArgumentError PC.RateCurve([2.0, 1.0], [0.1, 0.2])
end

@testset "Machine fingerprint" begin
    a = PC.machine_fingerprint()
    @test a == PC.machine_fingerprint()          # stable within a session
    @test length(a) == 16                        # 8 bytes hex

    # An explicit machine label must produce a distinct identity, so two
    # machines that happen to share a CPU model do not silently share constants.
    b = withenv("SPACEAGORA_PERF_MACHINE_LABEL" => "some-other-box") do
        PC.machine_fingerprint()
    end
    @test a != b
end

function _synthetic_constants(; ref_fma = 20.0)
    return PC.MachineConstants(
        simd_lane = PC.RateCurve([4.0, 8.0, 12.0], [0.24, 0.026, 0.044]),
        coeff_touch = PC.RateCurve([12.0, 16.0, 20.0], [0.10, 0.17, 0.30]),
        ns_per_scalar_item = 1.38,
        ns_per_queue_node = 0.055,
        dispatch_ns_base = 2396.0,
        dispatch_ns_per_worker = 447.0,
        ns_per_atomic = 19.9,
        reference_fma_ns = ref_fma,
        reference_mem_ns = 0.9,
        fingerprint = "testfingerprint",
        schema_version = PC.CALIBRATION_SCHEMA_VERSION,
    )
end

@testset "Constants round-trip through TOML" begin
    mktempdir() do dir
        path = joinpath(dir, "constants.toml")
        mc = _synthetic_constants()
        PC.save_machine_constants(mc, path)
        @test isfile(path)

        back = PC.load_machine_constants(path)
        @test back !== nothing
        @test back.fingerprint == mc.fingerprint
        @test back.ns_per_atomic ≈ mc.ns_per_atomic
        @test back.dispatch_ns_base ≈ mc.dispatch_ns_base
        @test back.simd_lane.ns ≈ mc.simd_lane.ns
        @test back.coeff_touch.log2_size ≈ mc.coeff_touch.log2_size
        # The curves must survive, not just the scalars: they are what
        # discriminates the routing candidates.
        @test PC.rate_at(back.coeff_touch, 2^14) ≈ PC.rate_at(mc.coeff_touch, 2^14)
    end
end

@testset "Stale or damaged constants are declined, not reinterpreted" begin
    mktempdir() do dir
        @test PC.load_machine_constants(joinpath(dir, "absent.toml")) === nothing

        garbage = joinpath(dir, "garbage.toml")
        write(garbage, "this is not toml {{{")
        @test PC.load_machine_constants(garbage) === nothing

        # A file from a different schema version must be discarded rather than
        # read: term semantics changed by definition when the version moved, so
        # its numbers do not mean what the current predictor would take them to.
        stale = joinpath(dir, "stale.toml")
        PC.save_machine_constants(_synthetic_constants(), stale)
        body = replace(read(stale, String),
                       "schema_version = $(PC.CALIBRATION_SCHEMA_VERSION)" => "schema_version = 999")
        write(stale, body)
        @test PC.load_machine_constants(stale) === nothing
    end
end

@testset "Staleness canary" begin
    # Constants recorded against this machine's actual reference must pass.
    live = PC.reference_kernel_ns()
    ok = PC.constants_are_current(_synthetic_constants(ref_fma = live))
    @test ok.ok
    @test abs(ok.drift) < 0.25

    # Constants from a machine four times faster must fail: that is the
    # categorical error the canary exists to catch -- a constants file that does
    # not describe this box, or a cgroup quota making its throughput
    # unreachable.
    bad = PC.constants_are_current(_synthetic_constants(ref_fma = live / 4))
    @test !bad.ok
    @test bad.drift > 0.25

    # A missing reference cannot be validated, so it is treated as stale.
    @test !PC.constants_are_current(_synthetic_constants(ref_fma = 0.0)).ok
end
