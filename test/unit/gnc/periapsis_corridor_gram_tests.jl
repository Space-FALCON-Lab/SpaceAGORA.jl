module PeriapsisCorridorGramTests

using Test
using CSV
using DataFrames
using StaticArrays
using TOML

# Keep the example definitions in a module so the script's PROGRAM_FILE guard
# cannot launch a mission while this test file is included by runtests.jl.
include(joinpath(@__DIR__, "../../../examples/AGORA_Earth_Periapsis_Corridor.jl"))

function _config_with(; backend=nothing)
    values = deepcopy(load_corridor_config().values)
    backend === nothing || (values["atmosphere"]["backend"] = backend)
    path = joinpath(mktempdir(), "periapsis_corridor_config.toml")
    open(io -> TOML.print(io, values), path, "w")
    return load_corridor_config(path)
end

@testset "periapsis corridor atmosphere backend selection (no GRAM needed)" begin
    @testset "default is the US76 fit and nothing else changes" begin
        config = load_corridor_config()
        @test corridor_atmosphere_backend(config) == "us76_piecewise"
        values = deepcopy(config.values)
        delete!(values["atmosphere"], "backend")
        path = joinpath(mktempdir(), "c.toml")
        open(io -> TOML.print(io, values), path, "w")
        @test corridor_atmosphere_backend(load_corridor_config(path)) == "us76_piecewise"
        @test config.values["mission"]["horizon_s"] == 5_184_000.0
        @test config.values["mission"]["corridor_min_altitude_m"] == 95_000.0
        @test config.values["mission"]["corridor_max_altitude_m"] == 100_000.0
        cli = _corridor_cli(String[])
        @test cli.atmosphere === nothing
        @test cli.output_dir == joinpath(REPO_ROOT, "output", "periapsis_corridor")
    end

    @testset "backend validation and CLI flag" begin
        @test corridor_atmosphere_backend(_config_with(backend="gram")) == "gram"
        @test_throws ArgumentError _config_with(backend="exponential")
        values = deepcopy(load_corridor_config().values)
        values["atmosphere"]["backend"] = "gram"
        delete!(values, "atmosphere_gram")
        path = joinpath(mktempdir(), "c.toml")
        open(io -> TOML.print(io, values), path, "w")
        @test_throws ArgumentError load_corridor_config(path)

        gram_cli = _corridor_cli(["--atmosphere=gram"])
        @test gram_cli.atmosphere == "gram"
        @test gram_cli.output_dir == joinpath(REPO_ROOT, "output", "periapsis_corridor_gram")
        @test _corridor_cli(["--atmosphere=gram", "--output=/tmp/x"]).output_dir == abspath("/tmp/x")
        @test _corridor_cli(["--atmosphere=us76_piecewise"]).output_dir ==
            joinpath(REPO_ROOT, "output", "periapsis_corridor")
        @test_throws ArgumentError _corridor_cli(["--atmosphere=nope"])
    end

    @testset "drag reads the sampled density only for the :sampled source" begin
        atmosphere = fitted_atmosphere_parameters(load_corridor_config())
        radius = 6_378_137.0
        pos = SVector(radius + 97_500.0, 0.0, 0.0)
        state = SimulationModel.StateSample(pos, SVector(0.0, 7000.0, 0.0), 820.0)
        planet_frame = SimulationModel.PlanetFrameSample(
            SMatrix{3, 3, Float64, 9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
            pos, SVector(0.0, 7000.0, 0.0), 97_500.0, 0.0, 0.0)
        sampled_rho = 2.0e-7
        env = SimulationModel.EnvironmentSample(nothing; planet_frame=planet_frame,
            atmosphere=SimulationModel.AtmosphereSample(sampled_rho, 190.0, SVector(0.0, 0.0, 0.0)))
        fit_drag = CorridorBallisticDrag(1000.0, 820.0, atmosphere)
        sampled_drag = CorridorBallisticDrag(1000.0, 820.0, atmosphere, :sampled)
        @test fit_drag.density_source === :fit
        f_fit, _ = SimulationModel.wrench(fit_drag, state, env, 0.0)
        f_sampled, _ = SimulationModel.wrench(sampled_drag, state, env, 0.0)
        # Same ballistic-coefficient treatment: force = -0.5 rho m v^2 / beta, along -v.
        @test f_sampled[2] ≈ -0.5 * sampled_rho * 820.0 * 7000.0^2 / 1000.0 rtol=1e-12
        @test f_fit[2] ≈ -0.5 * corridor_density(atmosphere, 97_500.0) * 820.0 * 7000.0^2 / 1000.0 rtol=1e-12
        @test f_sampled[2] / f_fit[2] ≈ sampled_rho / corridor_density(atmosphere, 97_500.0) rtol=1e-12
        no_atmosphere = SimulationModel.EnvironmentSample(nothing; planet_frame=planet_frame)
        @test_throws ErrorException SimulationModel.wrench(sampled_drag, state, no_atmosphere, 0.0)
    end

    @testset "missing GRAM fails clearly and never falls back to the fit" begin
        config = _config_with(backend="gram")
        original = _GRAM_LOADER[]
        try
            _GRAM_LOADER[] = () -> error("simulated: GRAMSuite is not importable")
            output_root = mktempdir()
            err = try
                run_corridor_case(config, "favorable", output_root; horizon_s=1000.0)
                nothing
            catch caught
                caught
            end
            @test err isa ErrorException
            @test occursin("atmosphere backend \"gram\" is unavailable", err.msg)
            @test occursin("NOT substituted", err.msg)
            @test occursin("simulated: GRAMSuite is not importable", err.msg)
            @test !isdir(joinpath(output_root, "favorable"))
            @test_throws ErrorException _corridor_gram_density_model(config, _CORRIDOR_EPOCH)
            @test_throws ArgumentError run_corridor_case(config, "vacuum", mktempdir(); horizon_s=1000.0)
        finally
            _GRAM_LOADER[] = original
        end
    end

    @testset "GRAM requires direct sampling and outputs cannot be mixed" begin
        for name in ("SPACEAGORA_GRAM_TRACK_CACHE", "SPACEAGORA_VACUUM_GRAM_CACHE",
                     "SPACEAGORA_DENSITY_FREEZE_PER_STEP")
            withenv(name => "on") do
                @test_throws ErrorException _corridor_gram_require_direct_sampling()
            end
        end
        withenv("SPACEAGORA_GRAM_TRACK_CACHE" => nothing, "SPACEAGORA_GRAM_SEGMENT_CACHE" => nothing,
                "SPACEAGORA_VACUUM_GRAM_CACHE" => nothing, "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => nothing) do
            @test _corridor_gram_require_direct_sampling() === nothing
        end

        case_dir = mktempdir()
        write_manifest(joinpath(case_dir, "run_manifest.toml"), _config_with();
            case_name="favorable", gravity_path=joinpath(REPO_ROOT, "data", "planets", "gravity", "EarthGGM05C.csv"))
        @test _assert_case_dir_backend(case_dir, "us76_piecewise") === nothing
        @test_throws ErrorException _assert_case_dir_backend(case_dir, "gram")
        @test _assert_case_dir_backend(joinpath(case_dir, "absent"), "gram") === nothing
    end

    @testset "manifest records backend, density source and provenance" begin
        gravity = joinpath(REPO_ROOT, "data", "planets", "gravity", "EarthGGM05C.csv")
        us76_path = joinpath(mktempdir(), "m.toml")
        write_manifest(us76_path, _config_with(); case_name="favorable", gravity_path=gravity)
        us76 = TOML.parsefile(us76_path)
        @test us76["atmosphere_backend"] == "us76_piecewise"
        @test us76["atmosphere_fit_used_for_drag"] == true
        gram_path = joinpath(mktempdir(), "m.toml")
        write_manifest(gram_path, _config_with(backend="gram"); case_name="favorable", gravity_path=gravity,
            atmosphere_provenance=Dict{String, Any}("backend" => "gram", "gram_core_epoch_utc" => "2020-1-1 0:0:0.0 UTC"))
        gram = TOML.parsefile(gram_path)
        @test gram["atmosphere_backend"] == "gram"
        @test gram["atmosphere_fit_used_for_drag"] == false
        @test occursin("GRAM", gram["density_source"])
        @test gram["atmosphere_provenance"]["gram_core_epoch_utc"] == "2020-1-1 0:0:0.0 UTC"
    end

    @testset "predictor epoch shifts by elapsed time" begin
        shifted = _shift_initial_time(_CORRIDOR_EPOCH, 86_400.5)
        @test (shifted.year, shifted.month, shifted.day, shifted.hour, shifted.minute) == (2020, 1, 2, 0, 0)
        @test shifted.second ≈ 0.5 atol=1e-9
        # 60 days from 2020-01-01 (leap year: 31 + 29 days) is 2020-03-01.
        march = _shift_initial_time(_CORRIDOR_EPOCH, 5_184_000.0)
        @test (march.year, march.month, march.day) == (2020, 3, 1)
        @test _shift_initial_time(_CORRIDOR_EPOCH, 0.0) === _CORRIDOR_EPOCH
    end
end

# GRAM-gated: needs the GRAMSuite extension, the native libGRAM and SPICE kernels.
const _GRAM_PROBE = try
    (model=_corridor_gram_density_model(_config_with(backend="gram"), _CORRIDOR_EPOCH), reason="")
catch err
    (model=nothing, reason=sprint(showerror, err))
end

if _GRAM_PROBE.model === nothing
    @info "SKIPPED GRAM-gated periapsis corridor tests: GRAM unavailable" reason = _GRAM_PROBE.reason
    @testset "GRAM-backed periapsis corridor (skipped: GRAM unavailable)" begin
        @test_skip false
    end
else
    @testset "GRAM-backed periapsis corridor" begin
        config = _config_with(backend="gram")
        model = _GRAM_PROBE.model
        @testset "GRAM is built at the configured run epoch, not 2000-01-01" begin
            @test _utc_string(model.core.initial_time) == _utc_string(_CORRIDOR_EPOCH)
            @test _utc_string(model.core.initial_time) == "2020-1-1 0:0:0.0 UTC"
            provenance = _corridor_gram_provenance(config, model, _CORRIDOR_EPOCH)
            @test provenance["gram_epoch_matches_run_epoch"] == true
            @test provenance["wind_applied_to_drag"] == false
            @test provenance["gram_density_ceiling_m"] == 2.0e6
        end

        @testset "mean density in the corridor band is finite, positive and decreasing" begin
            rho(h) = Base.invokelatest(SimulationModel.EnvironmentModels._gram_point_density,
                model, h, 0.0, 0.0, 0.0, false)[1]
            @test 0.0 < rho(100_000.0) < rho(97_500.0) < rho(95_000.0) < 1e-3
        end

        @testset "short GRAM run uses sampled density and records provenance" begin
            output_root = mktempdir()
            # Past the first periapsis (about 34,953 s) of the favorable case.
            summary = run_corridor_case(config, "favorable", output_root; horizon_s=40_000.0)
            @test summary.atmosphere_backend == "gram"
            @test summary.gram_density_evaluations > 0
            @test summary.mission_result == "reduced_horizon_diagnostic"
            passes = CSV.read(joinpath(output_root, "favorable", "passes.csv"), DataFrame)
            @test nrow(passes) >= 1
            @test all(passes.atmosphere_backend .== "gram")
            @test all(passes.density_kg_m3 .> 0.0)
            manifest = TOML.parsefile(joinpath(output_root, "favorable", "run_manifest.toml"))
            @test manifest["atmosphere_backend"] == "gram"
            @test manifest["atmosphere_provenance"]["gram_core_epoch_utc"] == "2020-1-1 0:0:0.0 UTC"
            @test manifest["atmosphere_provenance"]["gram_epoch_used_by_engine_utc"] == "2020-1-1 0:0:0.0 UTC"
        end
    end
end

end # module PeriapsisCorridorGramTests
