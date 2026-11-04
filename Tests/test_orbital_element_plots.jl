import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "9_SpaceAGORA.jl-main_full integration")))
using Test, Arrow, TOML, LinearAlgebra
include(joinpath(@__DIR__, "..", "Simulation", "Code", "plot_orbital_elements.jl"))
using .OrbitalElementPlots
const OEP = OrbitalElementPlots

@testset "Combined helper-count and target-altitude campaigns" begin
    @test OEP.COUNTS == (1, 50, 100, 150, 200, 250, 300, 400, 500, 600, 800)
    @test OEP.ALTITUDES == (1150, 1100, 1050, 1000, 950, 900, 850)
    @test OEP.CASE_COUNT == 462
end

@testset "Completed grid included, other helper altitudes and pending 700-helper cases excluded" begin
    mktempdir() do root
        for helper_altitude in (900, 1000, 1100), target_altitude in (900, 1050, 1100)
            directory = joinpath(root, "h$(helper_altitude)km_t$(target_altitude)km")
            mkpath(directory)
            case = OEP.StudyConfig.StudyCase(target_altitude, 1, 50, :gve_sma;
                helper_altitude_km=helper_altitude)
            OEP.StudyConfig.atomic_toml(joinpath(directory, "run_info.toml"),
                Dict("settings" => OEP.StudyConfig.case_settings(case)))
        end
        for helpers in (600, 700, 800)
            directory = joinpath(root, "N$helpers")
            mkpath(directory)
            case = OEP.StudyConfig.StudyCase(1050, 1, helpers, :gve_sma)
            OEP.StudyConfig.atomic_toml(joinpath(directory, "run_info.toml"),
                Dict("settings" => OEP.StudyConfig.case_settings(case)))
        end
        @test_throws r"Expected the complete 462-case study, found 5" OEP.load_grid(root, [80, 100])
        duplicate = joinpath(root, "duplicate")
        mkpath(duplicate)
        cp(joinpath(root, "h1000km_t1050km", "run_info.toml"), joinpath(duplicate, "run_info.toml"))
        @test_throws r"Ambiguous duplicate simulation" OEP.load_grid(root, [80, 100])
    end
end

@testset "Continuous-access population" begin
    radius = 7.3781366e6
    result = OEP.continuous_access(radius, radius, 200e3, 0.0)
    @test result.alpha_min_rad ≈ 2asin(200e3 / (2radius))
    @test result.minimum_helpers == 116
    @test 2pi / result.minimum_helpers <= 2result.alpha_min_rad
    @test 2pi / (result.minimum_helpers - 1) > 2result.alpha_min_rad
    @test OEP.continuous_access(radius, radius + 150e3, 200e3, deg2rad(1)).minimum_helpers == 972
    @test OEP.continuous_access(radius, radius - 150e3, 200e3, deg2rad(1)).minimum_helpers == 648
    @test ismissing(OEP.continuous_access(radius, radius + 200e3, 200e3, 0).minimum_helpers)
    @test isinf(OEP.continuous_access(radius, radius + 201e3, 200e3, 0).continuous_helpers)
    @test ismissing(OEP.continuous_access(radius, radius, 200e3, deg2rad(2)).minimum_helpers)
    @test_throws ArgumentError OEP.continuous_access(radius, radius, 200e3, pi / 2)
    @test_throws ArgumentError OEP.continuous_access(radius, radius, 200e3, -0.1)
    @test_throws ArgumentError OEP.continuous_access(radius, radius, 2radius, 0)
    @test_throws ArgumentError OEP.continuous_access(NaN, radius, 200e3, 0)
end

@testset "Orbital elements and angle branch" begin
    radius, mu = 7.3781366e6, OEP.StudyConfig.MU
    speed = sqrt(mu / radius)
    oe = OEP.orbital_elements([radius, 0, 0], [0, speed, 0], mu)
    @test oe.a ≈ radius
    @test oe.e < 1e-14
    @test oe.i == 0
    @test isnan(oe.raan)
    inclination = deg2rad(1)
    tilted = OEP.orbital_elements([radius, 0, 0],
        [0, speed * cos(inclination), speed * sin(inclination)], mu)
    @test tilted.i ≈ inclination
    @test tilted.raan ≈ 0
    @test rad2deg(OEP.signed_angle(deg2rad(359))) ≈ -1
    @test rad2deg(OEP.signed_angle(deg2rad(-359))) ≈ 1
    @test_throws ArgumentError OEP.orbital_elements(zeros(3), zeros(3), mu)
end

@testset "Exact saved windows and options" begin
    for altitude in OEP.ALTITUDES
        case = OEP.StudyConfig.StudyCase(altitude, 0.5, 50, :gve_sma)
        times = OEP.StudyConfig.trajectory_times(case)
        @test OEP.exact_sample_index(times, 80 * last(times) / 100) == 801
        @test OEP.exact_sample_index(times, last(times)) == 1001
        @test_throws ErrorException OEP.exact_sample_index(times, times[801] + 1)
    end
    @test OEP.parse_options(String[]).windows == [80, 100]
    @test OEP.parse_options(["--orbits", "100"]).windows == [100]
    @test_throws ErrorException OEP.parse_options(["--orbits", "101"])
    @test_throws ErrorException OEP.parse_options(["--orbits", "80,80"])
    @test_throws ErrorException OEP.parse_options(["--orbits"])
    @test_throws ErrorException OEP.parse_options(["--unknown", "x"])
end

@testset "Saved target endpoints (h=$altitude, N=$helpers, not helper sc1)" for
        altitude in (900, 1050, 1100), helpers in (50, 400, 500, 600, 800)
    case = OEP.StudyConfig.StudyCase(altitude, 1, helpers, :gve_sma)
    directory = OEP.StudyConfig.case_directory(OEP.StudyConfig.DEFAULT_OUTPUT, case)
    info = TOML.parsefile(joinpath(directory, "run_info.toml"))
    results = OEP.read_case(directory, info, [80, 100])
    table = Arrow.Table(joinpath(directory, "trajectory.feather"))
    function elements_at(index)
        target = helpers + 1
        position = [getproperty(table, Symbol("sc$(target)_pos_$axis"))[index] for axis in 1:3]
        velocity = [getproperty(table, Symbol("sc$(target)_vel_$axis"))[index] for axis in 1:3]
        OEP.orbital_elements(position, velocity, info["settings"]["mu_m3_s2"])
    end
    initial = elements_at(1)
    for (result, row) in zip(results, (801, 1001))
        final = elements_at(row)
        @test result.row == row
        @test result.time_s == table.time[row]
        @test result.da_m ≈ final.a - initial.a atol=1e-8
        @test result.de ≈ final.e - initial.e atol=1e-14
        @test result.di_deg ≈ rad2deg(final.i - initial.i) atol=1e-12
        @test result.draan_deg ≈ rad2deg(OEP.signed_angle(final.raan - initial.raan)) atol=1e-9
        @test result.raan_defined
    end
    @test_throws ErrorException OEP.read_case(directory, info, [101])
    bad_info = deepcopy(info)
    bad_info["status"] = "running"
    @test_throws ErrorException OEP.read_case(directory, bad_info, [80])
    bad_info = deepcopy(info)
    bad_info["checksums"]["trajectory.feather"] = "invalid"
    @test_throws ErrorException OEP.read_case(directory, bad_info, [80])
end
