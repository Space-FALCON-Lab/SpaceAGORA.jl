using Test
using TOML

module FullarcProvenanceFixture
    include(joinpath(@__DIR__, "..", "..", "scripts", "xval_fullarc_provenance.jl"))
end
const XVAL_PROVENANCE_TEST = FullarcProvenanceFixture.XvalFullarcProvenance

@testset "committed full-arc runs reject sensitivity inputs" begin
    keys = ["XVAL_BASE", "XVAL_J0_GM", "XVAL_GMAT_DIR", "XVAL_EARTH_FIELD_FILE",
            "XVAL_MOON_FIELD_FILE", "XVAL_MOON_FIELD_GM", "XVAL_PLANETARY_KERNEL",
            "XVAL_FRAME_TABLE_DIR", "SPACEAGORA_SPICE_PCK_OVERRIDES",
            "SPACEAGORA_SPICE_PLANETARY_KERNEL_RELPATH", "XVAL_FUTURE_MODEL_OVERRIDE",
            "SPACEAGORA_TELEMETRY_J2_SOURCE_DEFAULT", "SPACEAGORA_TELEMETRY_HARMONICS_NORMALIZED_DEFAULT",
            "SPACEAGORA_TELEMETRY_SOLVER_MAXITERS", "SPACEAGORA_GMAT_PARITY_SOLVER",
            "SPACEAGORA_SOLVER_SAVE_EVERYSTEP"]
    for key in keys
        env = Dict(key => "diagnostic-input")
        @test_throws ArgumentError XVAL_PROVENANCE_TEST.run_info("gmat", "committed", ["earth_j0_tbfalse"]; env)
        info = XVAL_PROVENANCE_TEST.run_info("gmat", "frame_diagnostic", ["earth_j0_tbfalse"]; env)
        @test !info["primary_eligible"]
        @test info["input_overrides"][key] == "diagnostic-input"
    end
    info = XVAL_PROVENANCE_TEST.run_info("stk", "committed", ["earth_j0_tbfalse"];
        env=Dict("XVAL_SCENARIOS" => "earth_j0_tbfalse", "XVAL_FRAME_TABLE_DIR" => ""))
    @test_throws ArgumentError XVAL_PROVENANCE_TEST.run_info("gmat", "committed", String[];
        env=Dict("XVAL_FRAME_TABLE_DIR" => " "))
    @test isempty(info["input_overrides"])
    @test info["scenarios"] == ["earth_j0_tbfalse"]
end

@testset "completion binds artifacts and a rerun invalidates it first" begin
    mktempdir() do dir
        info = XVAL_PROVENANCE_TEST.run_info("gmat", "committed", ["earth_j0_tbfalse"]; env=Dict())
        path = joinpath(dir, "run_info.toml")
        XVAL_PROVENANCE_TEST.write_info(path, info)
        @test TOML.parsefile(path)["status"] == "running"
        mkpath(joinpath(dir, "earth_j0_tbfalse"))
        for relative in ("results.csv", "earth_j0_tbfalse/manifest.toml", "earth_j0_tbfalse/series.arrow")
            write(joinpath(dir, relative), "synthetic fixture")
        end
        XVAL_PROVENANCE_TEST.finish_run!(dir, info)
        complete = TOML.parsefile(path)
        @test complete["status"] == "complete"
        @test length(complete["artifacts_sha256"]) == 3
        @test complete["artifacts_sha256"]["results.csv"] == XVAL_PROVENANCE_TEST.file_digest(joinpath(dir, "results.csv"))
        fresh = XVAL_PROVENANCE_TEST.run_info("gmat", "committed", ["moon_j0_tbfalse"]; env=Dict())
        XVAL_PROVENANCE_TEST.write_info(path, fresh)
        @test TOML.parsefile(path)["status"] == "running"
        @test !haskey(TOML.parsefile(path), "artifacts_sha256")
    end
end
