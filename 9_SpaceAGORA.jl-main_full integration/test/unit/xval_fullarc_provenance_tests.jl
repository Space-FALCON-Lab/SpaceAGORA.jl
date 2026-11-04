using Test, JSON
module FullarcPreflightFixture
    include(joinpath(@__DIR__, "..", "..", "scripts", "xval_provenance.jl"))
end
const XVAL_PREFLIGHT = FullarcPreflightFixture.XvalProvenance
@testset "Full-arc preflight preserves published rejection cases" begin
    keys = ["XVAL_BASE", "XVAL_J0_GM", "XVAL_GMAT_DIR", "XVAL_EARTH_FIELD_FILE",
            "XVAL_MOON_FIELD_FILE", "XVAL_MOON_FIELD_GM", "XVAL_PLANETARY_KERNEL",
            "XVAL_FRAME_TABLE_DIR", "SPACEAGORA_SPICE_PCK_OVERRIDES",
            "SPACEAGORA_SPICE_PLANETARY_KERNEL_RELPATH", "XVAL_FUTURE_MODEL_OVERRIDE",
            "SPACEAGORA_TELEMETRY_J2_SOURCE_DEFAULT", "SPACEAGORA_TELEMETRY_HARMONICS_NORMALIZED_DEFAULT",
            "SPACEAGORA_TELEMETRY_SOLVER_MAXITERS", "SPACEAGORA_GMAT_PARITY_SOLVER",
            "SPACEAGORA_SOLVER_SAVE_EVERYSTEP"]
    for key in keys
        env = Dict(key => "diagnostic-input")
        @test_throws ArgumentError XVAL_PREFLIGHT.check_primary_controls("gmat", "committed"; env)
        @test XVAL_PREFLIGHT.check_primary_controls("gmat", "frame_diagnostic"; env) === nothing
        @test XVAL_PREFLIGHT.controls(env)[key] == "diagnostic-input"
    end
    @test XVAL_PREFLIGHT.check_primary_controls("stk", "committed";
        env=Dict("XVAL_SCENARIOS" => "earth_j0_tbfalse", "XVAL_FRAME_TABLE_DIR" => "")) === nothing
    @test_throws ArgumentError XVAL_PREFLIGHT.check_primary_controls("gmat", "committed";
        env=Dict("XVAL_FRAME_TABLE_DIR" => " "))
    mktempdir() do root
        run = joinpath(root, "gmat_committed")
        withenv("XVAL_GMAT_DIR" => "diagnostic") do
            @test_throws ArgumentError XVAL_PREFLIGHT.start_record(root, run, "gmat", "committed", [])
            @test !ispath(run)
        end
    end
end
