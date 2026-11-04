import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "9_SpaceAGORA.jl-main_full integration")))
using Test, Arrow, TOML, Plots
include(joinpath(@__DIR__, "..", "Simulation", "Code", "plot_orbital_elements.jl"))
const OE = OrbitalElementPlots

@testset "Exact saved cumulative RTN endpoints" begin
    for (altitude, inclination, helpers, policy) in
            ((850, 1, 800, :gve_sma), (1100, 0.5, 600, :gve_inc), (1000, 0, 50, :gve_inc))
        case = OE.StudyConfig.StudyCase(altitude, inclination, helpers, policy)
        directory = OE.StudyConfig.case_directory(OE.StudyConfig.DEFAULT_OUTPUT, case)
        info = TOML.parsefile(joinpath(directory, "run_info.toml"))
        results = OE.read_case(directory, info, [80, 100]; include_rtn=true)
        original = OE.read_case(directory, info, [80, 100])
        table = Arrow.Table(joinpath(directory, "trajectory.feather"))
        for (result, old, index) in zip(results, original, (801, 1001))
            @test result.row == index
            @test result.time_s == table.time[index]
            @test result.dv_R_mps == table.target_dv_R_mps[index]
            @test result.dv_T_mps == table.target_dv_T_mps[index]
            @test result.dv_N_mps == table.target_dv_N_mps[index]
            @test all(name -> isequal(getproperty(result, name), getproperty(old, name)), keys(old))
        end
    end
end

@testset "Three-row RTN rendering" begin
    rows = [(scheduler=policy, target_altitude_km=altitude, delta_i_deg=inclination,
        helpers=helpers, orbits=window, dv_R_mps=(altitude - 1000) * helpers * 1e-4,
        dv_T_mps=helpers * 0.01, dv_N_mps=inclination * helpers * -1e-3)
        for window in (80, 100) for policy in OE.POLICIES for altitude in OE.ALTITUDES
        for inclination in OE.INCLINATIONS for helpers in OE.COUNTS]
    thresholds = [merge((target_altitude_km=altitude, delta_i_deg=inclination,),
        OE.continuous_access(OE.EARTH_RADIUS_M + 1000e3,
            OE.EARTH_RADIUS_M + altitude * 1000, OE.RANGE_M, deg2rad(inclination)))
        for altitude in OE.ALTITUDES for inclination in OE.INCLINATIONS]
    fields = (:dv_R_mps, :dv_T_mps, :dv_N_mps)
    colors = ("#652884", "#4E5689", "#018B38", "#B46DA9", "#41B9C1", "#57AF37")
    gr()
    mktempdir() do output
        for window in (80, 100)
            figure = OE.make_plot(rows, thresholds, "combined", window, output; quantity=:rtn)
            @test length(figure.subplots) == 22
            for (index, panel) in enumerate(figure.subplots[1:21])
                r, c = divrem(index - 1, 7) .+ 1
                curves = filter(s -> s[:seriestype] == :path, panel.series_list)
                @test length(curves) == 6
                @test panel[:xaxis][:grid] && panel[:yaxis][:grid]
                @test panel[:xaxis][:guide] == (r == 3 ? "N_helpers" : "")
                for (j, policy) in enumerate(("gve_sma", "gve_inc")),
                        (k, inclination) in enumerate(OE.INCLINATIONS)
                    series = curves[(j - 1) * 3 + k]
                    expected = sort(filter(row -> row.orbits == window &&
                        row.scheduler == policy && row.delta_i_deg == inclination &&
                        row.target_altitude_km == OE.ALTITUDES[c], rows); by=row -> row.helpers)
                    @test collect(series[:x]) == collect(OE.COUNTS)
                    @test collect(series[:y]) == [getproperty(row, fields[r]) for row in expected]
                    @test series[:linecolor] == Plots.plot_color(colors[(j - 1) * 3 + k])
                    @test series[:linestyle] == (policy == "gve_inc" ? :dot : :solid)
                end
            end
            for extension in ("png", "pdf")
                @test filesize(joinpath(output,
                    "dv_RTN_from_feather_combined_e0.0000_nu1.0000_$(window)orbits_h_cols.$extension")) > 10_000
            end
        end
        @test_throws ArgumentError OE.make_plot(rows, thresholds, "combined", 80, output; quantity=:unknown)
    end
end
