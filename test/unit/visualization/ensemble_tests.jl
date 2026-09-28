using Test
using SpaceAGORA
using StaticArrays
using Arrow
using DataFrames
using JSON
using Base64

import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft

const SM = SpaceAGORA.SimulationModel
const SV = SM.SceneVisualization

# Mars, two-body + J2, no atmosphere; the seed only moves the starting true
# anomaly and the mission length, so samples end at different times.
function _sample_args(seed::Integer; results_directory::String=mktempdir())
    planet = make_no_gram_planet(:mars)
    spacecraft = make_three_body_spacecraft(
        bus_dims=(2.0, 2.0, 2.5),
        panel_dims=(0.01, 2.5, 1.0),
        bus_mass=500.0,
        panel_mass_each=10.0,
        panel_offset_y=2.3,
        ic=SM.InitialCondition(ra=4_500.0e3, rp=3_800.0e3, i=30.0, ω=0.0, Ω=45.0, ν=10.0 * seed),
        prop_mass=0.0,
        id=1
    )
    return make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=300.0 + 100.0 * seed,
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=(SM.InverseSquaredJ2GravityModel(),),
        density_model=SM.NoAtmosphereModel(),
        ephemerides_model=SM.SimpleEphemeridesModel(),
        orientation_sim=false,
        keplerian=true,
        EI_km=120.0,
        verbose=false,
        results=true,
        results_directory=results_directory
    )
end

_decode(b64, T) = collect(reinterpret(T, base64decode(b64)))

function _payload_of(html_path)
    html = read(html_path, String)
    start = findfirst("window.SPACEAGORA_VIEWER = ", html)
    stop = findnext(";\n</script>", html, last(start))
    return JSON.parse(html[last(start)+1:first(stop)-1])
end

@testset "Ensembles" begin
    @testset "time axis and resampling helpers" begin
        @test SV.ensemble_time_axis(Float64[], Float64[]) == [0.0]
        t = SV.ensemble_time_axis([100.0, 250.0], [10.0, 5.0]; max_frames=1000)
        @test t[1] == 0.0 && t[end] == 250.0
        @test all(diff(t) .≈ 5.0)
        coarse = SV.ensemble_time_axis([1000.0], [1.0]; max_frames=11)
        @test length(coarse) == 11 && coarse[end] == 1000.0

        out = zeros(5)
        SV._resample_linear!(out, [0.0, 10.0, 20.0], [0.0, 100.0, 0.0], [-1.0, 0.0, 5.0, 20.0, 25.0])
        @test isnan(out[1]) && out[2] == 0.0 && out[3] == 50.0 && out[4] == 0.0 && isnan(out[5])
        q = zeros(4, 3)
        src = [0.0 0.0; 0.0 0.0; 0.0 sin(pi / 4); 1.0 cos(pi / 4)]
        SV._resample_quaternions!(q, [0.0, 1.0], src, [0.0, 0.5, 2.0])
        @test q[:, 1] ≈ [0.0, 0.0, 0.0, 1.0]
        @test q[3, 2] ≈ sin(pi / 8) && q[4, 2] ≈ cos(pi / 8)
        @test all(isnan, q[:, 3])
    end

    @testset "sample directories and manifest round trip" begin
        dir = mktempdir()
        @test SV.sample_results_directory(dir, 7) == joinpath(dir, "sample_0007")
        args = with_results_directory(_sample_args(1), joinpath(dir, "sample_0001"))
        @test args.simulation_settings.results_directory == joinpath(dir, "sample_0001")
        @test args.simulation_settings.save_visualization_scene
        @test args.simulation_settings.results
        samples = [
            EnsembleSample(1, "11", true, 250.5, "first", joinpath(dir, "sample_0001")),
            EnsembleSample(2, "12", false, NaN, "second", joinpath(dir, "sample_0002")),
        ]
        path = write_ensemble_manifest(dir, samples; scalar_name="periapsis (km)")
        @test isfile(path)
        raw = JSON.parsefile(path)
        @test raw["samples"][1]["directory"] == "sample_0001"
        @test raw["samples"][2]["scalar"] === nothing
        back, name = SV.read_ensemble_manifest(dir)
        @test name == "periapsis (km)"
        @test back[1].directory == joinpath(dir, "sample_0001")
        @test back[1].scalar == 250.5 && isnan(back[2].scalar) && !back[2].success
        @test SV.default_sample_scalar(DataFrame(time=[0.0])) |> isnan
        @test SV.default_sample_scalar(DataFrame(time=[0.0, 1.0], sc1_periapsis_altitude=[1.0e5, 2.5e5])) == 250.0
    end

    @testset "Monte Carlo campaign writes samples, manifest and page" begin
        campaign = mktempdir()
        result, page = run_monte_carlo_visualization(_sample_args, [1, 2], campaign; threads=1, max_frames=200)
        @test length(result.successful) == 2
        @test isfile(page) && page == joinpath(campaign, "ensemble_viewer.html")
        @test isfile(joinpath(campaign, "sample_0001", "simulation_results_scene.json"))
        @test isfile(joinpath(campaign, "sample_0002", "simulation_results.feather"))
        samples, name = SV.read_ensemble_manifest(campaign)
        @test length(samples) == 2 && name == "final spherical periapsis altitude (km)"
        @test all(s -> isfinite(s.scalar) && s.success, samples)
        @test samples[1].seed == "1" && samples[2].seed == "2"
        @test all(s -> isfinite(s.value), result.successful)

        payload = _payload_of(page)
        @test payload["textures"]["mars"]["resolution"] == "4k"
        detailed = export_ensemble_visualization(campaign; out=joinpath(campaign, "8k.html"), texture_resolution="8k")
        @test _payload_of(detailed)["textures"]["mars"]["resolution"] == "8k"
        @test payload["ensemble"]["count"] == 2
        @test payload["ensemble"]["spacecraft_per_sample"] == 1
        @test payload["ensemble"]["samples"][2]["label"] == "sample 2 (seed 2)"
        frames = payload["frames"]
        @test frames["sats"] == 2
        t = _decode(frames["t_s"], Float64)
        @test t[1] == 0.0 && t[end] == 500.0      # the longer sample (seed 2) sets the axis
        @test frames["pos_dtype"] == "f64"
        pos = _decode(frames["pos_km"], Float64)
        N = frames["count"]
        @test length(pos) == N * 2 * 3
        # Sample 1 ends at 400 s: present before, absent after; sample 2 present throughout.
        k_before = findlast(<=(390.0), t)
        k_after = findfirst(>(400.0), t)
        @test isfinite(pos[((k_before - 1) * 2 + 0) * 3 + 1])
        @test isnan(pos[((k_after - 1) * 2 + 0) * 3 + 1])
        @test isfinite(pos[((k_after - 1) * 2 + 1) * 3 + 1])
        @test frames["link_pose"]["counts"] == [2, 2]
        @test length(payload["scene"]["spacecraft"]) == 2
        @test startswith(payload["scene"]["spacecraft"][2]["name"], "sample 2")
        @test occursin("ensemble", read(page, String))

        # Without the manifest the sample directories are discovered by name.
        rm(joinpath(campaign, "ensemble_manifest.json"))
        found, fname = SV.discover_ensemble_samples(campaign)
        @test length(found) == 2 && fname == ""
        @test all(s -> isnan(s.scalar), found)
        @test found[2].label == "sample_0002"
        again = export_ensemble_visualization(campaign; out=joinpath(campaign, "again.html"), textures=false, scalar_name="none")
        @test isfile(again)
        @test _payload_of(again)["ensemble"]["scalar_name"] == "none"

        # Constellation-ensemble member directories use the same exporter.
        members = mktempdir()
        cp(joinpath(campaign, "sample_0001"), joinpath(members, "sat_1_id_1"))
        cp(joinpath(campaign, "sample_0002"), joinpath(members, "sat_2_id_2"))
        member_page = export_ensemble_visualization(members; textures=false)
        @test _payload_of(member_page)["ensemble"]["samples"][1]["label"] == "sat_1_id_1"

        @test_throws ArgumentError export_ensemble_visualization(mktempdir())

        cli_out = IOBuffer()
        @test run_cli(["visualize", "--run=$(members)", "--ensemble", "--no-textures", "--out=$(joinpath(members, "cli.html"))"]; io=cli_out) == 0
        @test isfile(joinpath(members, "cli.html"))
    end
    @testset "checkpoint files and resumed states stay with each sample" begin
        mktempdir() do root
            campaign = joinpath(root, "campaign")
            checkpoint_root = joinpath(root, "shared_checkpoints")
            solver = SM.SolverConfig(solver_mode=:tsit5)
            function checkpoint_args(seed, duration; resume=false, checkpoint=true)
                base = _sample_args(seed; results_directory=joinpath(root, "builder_output"))
                return SM.SimConfig._with_configuration(base;
                    mission_configuration=SM.MissionConfiguration(
                        mission_time=duration, orientation_sim=false, keplerian=true,
                        data_rate=0.5,
                    ),
                    simulation_settings=SM.SimulationSettings(
                        results=true, verbose=false, generate_plots=false, save_csv=false,
                        results_directory=base.simulation_settings.results_directory,
                        checkpoint_enabled=checkpoint, checkpoint_interval_s=1.0,
                        checkpoint_directory=checkpoint_root, resume_from_checkpoint=resume,
                    ),
                    solver_config=solver,
                )
            end
            first_args = [checkpoint_args(seed, 2.0) for seed in 1:2]
            result, page = run_monte_carlo_visualization(
                seed -> first_args[seed], [1, 2], campaign;
                threads=1, fail_fast=true, export_page=false,
            )
            @test length(result.successful) == 2
            @test page === nothing
            sample_dirs = [SV.sample_results_directory(campaign, seed) for seed in 1:2]
            checkpoint_dirs = [joinpath(checkpoint_root, basename(dir)) for dir in sample_dirs]
            @test checkpoint_dirs[1] != checkpoint_dirs[2]
            @test !isfile(joinpath(checkpoint_root, "simulation_checkpoint.bin"))
            @test all(dir -> isfile(joinpath(dir, "simulation_checkpoint.bin")), checkpoint_dirs)
            @test all(dir -> isfile(joinpath(dir, "simulation_checkpoint.manifest.toml")), checkpoint_dirs)
            first_tables = [DataFrame(Arrow.Table(joinpath(dir, "simulation_results.feather"))) for dir in sample_dirs]
            position(df, row) = Float64[df[row, Symbol("sc1_pos_", axis)] for axis in 1:3]
            first_positions = [position(df, nrow(df)) for df in first_tables]
            @test all(df -> df.time[end] ≈ 2.0, first_tables)
            # The seeds start far apart: reusing the other member's checkpoint
            # cannot pass the continuation comparison accidentally.
            @test sum(abs2, first_positions[1] - first_positions[2]) > 1e8

            # Resume-only mode also needs isolated checkpoint paths, even though
            # checkpoint_enabled is false. Run each seed from 2 s through 4 s.
            resume_args = [checkpoint_args(seed, 4.0; resume=true, checkpoint=false) for seed in 1:2]
            resumed, _ = run_monte_carlo_visualization(
                seed -> resume_args[seed], [1, 2], campaign;
                threads=1, fail_fast=true, export_page=false,
            )
            @test length(resumed.successful) == 2
            for seed in 1:2
                df = DataFrame(Arrow.Table(joinpath(sample_dirs[seed], "simulation_results.feather")))
                @test nrow(df) >= 2
                @test issorted(df.time)
                @test df.time[1] ≈ 2.0
                @test df.time[end] ≈ 4.0
                @test isapprox(position(df, 1), first_positions[seed]; rtol=0.0, atol=1e-6)
            end
            # Campaign wrapping must leave builder-owned settings and the typed
            # solver untouched in both the first run and the continuation.
            for args in (first_args..., resume_args...)
                @test args.simulation_settings.checkpoint_directory == checkpoint_root
                @test args.simulation_settings.results_directory == joinpath(root, "builder_output")
                @test !args.simulation_settings.save_visualization_scene
                @test args.solver_config === solver
            end
        end
    end

end
