using Test
include("test16_options.jl")

@testset "OracleOptions controls the animation case" begin
    for feather_only in (false, true)
        mktempdir() do directory
            opts = OracleOptions(helpers=1, orbits=0.002,
                schedule=feather_only ? :none : :gve_sma,
                dt_max_s=0.25, mass_kg=123.0, laser_power_w=5000.0,
                magnification=50.0, use_J2=false, use_los=false,
                feather_only=feather_only, result_plots=false,
                animate=true, show_earth=false, duration_seconds=1.0,
                animation_fps=2.0, output_dir=directory)
            result = run_animation_case(opts)
            expected_duration = opts.orbits * 2pi * sqrt((R_EARTH + opts.target_altitude_km*1e3)^3 / MU)
            @test string(result.sol.retcode) == "Success"
            @test result.sol.t[end] ≈ expected_duration
            @test result.sol.t == output_times(expected_duration)
            @test result.sol.destats.naccept >= floor(Int, expected_duration / opts.dt_max_s)
            @test result.params[:gve_schedule] == string(opts.schedule)
            @test !result.params[:use_J2]
            @test !result.params[:use_los]
            @test result.params[:masses] == fill(opts.mass_kg, 2)
            @test result.params[:cavity][(1, 2)][:Pin] == opts.laser_power_w
            @test result.params[:cavity][(1, 2)][:B] == opts.magnification
            @test occursin("_$(opts.schedule)_J2false_", result.paths.scenario)
            @test isfile(result.feather_path)
            @test isfile(joinpath(result.paths.feather, "laser_on.feather"))
            @test isdir(result.paths.csv) == !feather_only
            video = joinpath(directory, "videos", result.paths.scenario, "animation.mp4")
            @test isfile(video) == !feather_only
            feather_only || @test filesize(video) > 0
        end
    end
    for opts in (OracleOptions(planet=:mars), OracleOptions(paper_grid=true),
                 OracleOptions(beta=0.5), OracleOptions(eta=0.5),
                 OracleOptions(timeseries_points=200), OracleOptions(schedule=:invalid),
                 OracleOptions(orbits=0.0), OracleOptions(dt_max_s=0.0))
        @test_throws ArgumentError run_animation_case(opts)
    end
end