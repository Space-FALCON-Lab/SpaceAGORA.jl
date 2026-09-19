using Test, Libdl

# A dedicated process prevents native-wrapper state from other suites masking
# accidental initialization. Load ordinary packages and their extension only.
const GRID_TEST_REPO = normpath(joinpath(@__DIR__, "..", ".."))
const GRID_TEST_VENDORED = joinpath(GRID_TEST_REPO, "data", "GRAMSuite.jl")
if Base.find_package("GRAMSuite") === nothing && isfile(joinpath(GRID_TEST_VENDORED, "Project.toml"))
    pushfirst!(LOAD_PATH, GRID_TEST_VENDORED)
end
ENV["SPACEAGORA_GRAM_STATIC_GRID"] = "off"
ENV["SPACEAGORA_GRAM_OFFLINE_SURROGATE"] = "off"

using SpaceAGORA, GRAMSuite

@testset "Normal grid-adapter package import" begin
    @test realpath(dirname(dirname(pathof(SpaceAGORA)))) == realpath(GRID_TEST_REPO)
    if haskey(ENV, "EXPECTED_GRAMSUITE_ROOT")
        @test realpath(dirname(dirname(pathof(GRAMSuite)))) == realpath(ENV["EXPECTED_GRAMSUITE_ROOT"])
    end
    @test Base.get_extension(SpaceAGORA, :SpaceAGORAGRAMSuiteExt) !== nothing
    @test :GRAMGridAtmosphereModel in names(SpaceAGORA)
    @test SpaceAGORA.GRAMGridAtmosphereModel === SpaceAGORA.SimulationModel.EnvironmentModels.GRAMGridAtmosphereModel
    for symbol in (:with_density_model_epoch, :MeshAeroSurrogate, :NoTerrainModel)
        @test symbol in names(SpaceAGORA)
    end
    @test Threads.nthreads() >= 4
    @test GRAMSuite._GRAM_WRAPPER[] === nothing
    @test !any(path -> occursin("libgram", lowercase(path)), Libdl.dllist())
end

module GridInterfaceTests
using SpaceAGORA, GRAMSuite
include("interface_tests.jl")
end

module GridSchemaTests
using SpaceAGORA, GRAMSuite
include("schema_tests.jl")
end

# CI uses small synthetic payloads. The 654-point retained Odyssey fixture is
# opt-in and must identify all three externally supplied inputs by SHA256.
const GRID_REFERENCE_KEYS = ["SPACEAGORA_TEST_GRID_" * name * suffix
    for name in ("FILE", "POINTS", "REFERENCE") for suffix in ("", "_SHA256")]
const GRID_REFERENCE_REQUESTED = any(key -> haskey(ENV, key), GRID_REFERENCE_KEYS)
if GRID_REFERENCE_REQUESTED
    all(key -> haskey(ENV, key) && !isempty(ENV[key]), GRID_REFERENCE_KEYS) ||
        error("Optional Odyssey regression requires paths and SHA256 for FILE, POINTS and REFERENCE")
    @eval module GridOdysseyReferenceTests
        using SpaceAGORA, GRAMSuite
        include("odyssey_reference_tests.jl")
    end
else
    @info "Optional retained Odyssey payload regression was not requested"
end

@testset "Grid suite never initialized native GRAM" begin
    @test GRAMSuite._GRAM_WRAPPER[] === nothing
    @test GRAMSuite._GRAM_WRAPPER_FILE[] == ""
    @test !any(path -> occursin("libgram", lowercase(path)), Libdl.dllist())
end
