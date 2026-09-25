using Test, TOML, Libdl, SpaceAGORA, GRAMSuite

@testset "Public preset CLI and data" begin
    @test GRAMSuite._GRAM_WRAPPER[] === nothing
    @test !any(p -> occursin("libgram", lowercase(p)), Libdl.dllist())
    listing = sprint(io -> @test SpaceAGORA.run_cli(["assets", "list"]; io) == 0)
    @test occursin("odyssey_p20_frozen_v1@1.0.0", listing)
    @test_throws ArgumentError SpaceAGORA.run_cli(["assets", "fetch", "--preset=odyssey_p20_frozen_v1"])
    @test_throws ArgumentError SpaceAGORA.run_cli(["assets", "check", "--preset=x", "--version=latest"])
    @test_throws ArgumentError SpaceAGORA.run_cli(["assets", "fetch", "--preset=x", "--preset=y", "--version=1.0.0"])
    model = surrogate_preset_model("odyssey_p20_frozen_v1"; version="1.0.0", offline=true)
    metadata = atmosphere_provenance(model)
    @test metadata["backend"] == "gram_grid_surrogate"
    assets = odyssey_surrogate_assets(; offline=true)
    @test length(assets.kernels) == 4
    @test isfile(assets.gravity)
    validation = sprint(io -> @test SpaceAGORA.run_cli(["assets", "check", "--preset=odyssey_p20_frozen_v1", "--version=1.0.0"]; io) == 0)
    @test occursin("3643c9116b75c511", validation)
    # The ordinary results writer includes the resolved atmosphere record.
    output = mktempdir()
    args = (; environment_model=(; density_model=model),
        mission_configuration=(; mission_time=1.0, orientation_sim=false),
        dynamics_model=(; spacecraft=(nothing,)),
        simulation_settings=(; save_csv=false, results_directory=output))
    table = SpaceAGORA.SimulationModel.DataFrames.DataFrame(time=[0.0, 1.0])
    SpaceAGORA.SimulationModel.IOOutputs._write_results_bundle!(table, [0.0, 1.0], args, 1)
    saved = TOML.parsefile(joinpath(output, "simulation_results.manifest.toml"))
    @test saved["atmosphere"] == metadata
    @test saved["atmosphere"]["preset_id"] == "odyssey_p20_frozen_v1"
    @test saved["atmosphere"]["preset_version"] == "1.0.0"
    @test saved["atmosphere"]["source_sha256"] == "3643c9116b75c511d20edee2866e2b5ba06baaeba9e693c4f3ad1055a7820e8c"
    @test GRAMSuite._GRAM_WRAPPER[] === nothing
    @test !any(p -> occursin("libgram", lowercase(p)), Libdl.dllist())
end

@testset "Active public Odyssey outputs" begin
    directory = joinpath(@__DIR__, "..", "output", "odyssey-surrogate")
    comparison = TOML.parsefile(joinpath(directory, "comparison.toml"))
    @test comparison["effect_assertions_passed"]
    @test !comparison["native_gram_loaded"]
    @test comparison["comparison_time_s"] == 600.0
    @test comparison["position_difference_m"] > 0
    @test comparison["panel_heat_load_difference_norm_J_cm2"] > 0
    @test comparison["baseline_cap_deg"] == 90.0 && comparison["variant_cap_deg"] == 30.0
    for cap in (90, 30)
        summary = TOML.parsefile(joinpath(directory, "cap_$(cap)_deg", "summary.toml"))
        @test summary["panel_cap_active"]
        @test !summary["thermal_feedback_tested"]
        @test startswith(summary["state_frame"], "Mars-centred J2000 inertial")
        @test summary["atmosphere"]["preset_version"] == "1.0.0"
        @test summary["atmosphere"]["source_sha256"] == "3643c9116b75c511d20edee2866e2b5ba06baaeba9e693c4f3ad1055a7820e8c"
        @test summary["scenario_assets"]["version"] == "1.0.0"
    end
end

@testset "Documented CLI entrypoint activates the adapter" begin
    script = joinpath(@__DIR__, "..", "src", "cli", "main.jl")
    project = dirname(Base.active_project())
    command = `$(Base.julia_cmd()) --startup-file=no --project=$project $script assets check --preset=odyssey_p20_frozen_v1 --version=1.0.0`
    text = read(command, String)
    @test occursin("Validated odyssey_p20_frozen_v1@1.0.0", text)
    @test occursin("3643c9116b75c511", text)
end

# Top-level include, so the testset below runs in a world that sees the example's methods.
include(joinpath(@__DIR__, "..", "examples", "odyssey_surrogate.jl"))

@testset "Odyssey example options fail before any run" begin
    E = OdysseySurrogateExample
    @test E.parse_arguments(["--offline", "--help"]) === nothing
    defaults = E.parse_arguments(String[])
    @test (defaults.baseline_cap_deg, defaults.variant_cap_deg, defaults.offline) == (90.0, 30.0, false)
    @test defaults.output_dir == joinpath(pwd(), "odyssey_surrogate_results")
    @test E.parse_arguments(["--cap=60", "--baseline-cap=85", "--output=new", "--offline"]) ==
        (; output_dir="new", offline=true, baseline_cap_deg=85.0, variant_cap_deg=60.0)
    for argv in (["--cap=0"], ["--cap=90.5"], ["--cap=NaN"], ["--cap=Inf"], ["--cap=sixty"], ["--cap="],
            ["--baseline-cap=-1"], ["--cap=60", "--cap=45"], ["--output="], ["--caps=60"])
        @test_throws ArgumentError E.parse_arguments(argv)
    end
    @test E.cap_directory(90.0) == "cap_90_deg" && E.cap_directory(22.5) == "cap_22.5_deg"
    mktempdir() do directory
        unused = joinpath(directory, "unused")
        @test_throws ArgumentError E.compare_panel_caps(; output_dir=unused, variant_cap_deg=95)
        @test_throws ArgumentError E.compare_panel_caps(; output_dir=unused, baseline_cap_deg=60, variant_cap_deg=60.0)
        @test_throws ArgumentError E.run_case(; panel_cap_deg=0, output_dir=unused)
        @test !ispath(unused)
        existing = try E.compare_panel_caps(; output_dir=directory) catch err err end
        @test existing isa ArgumentError && occursin("--output=DIR", existing.msg)
    end
end
