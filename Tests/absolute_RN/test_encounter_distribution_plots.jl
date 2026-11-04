import Pkg
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(ROOT, "9_SpaceAGORA.jl-main_full integration"))
using Test, CSV, DataFrames, Plots
include(joinpath(ROOT, "Simulation_absolute_RN", "Code", "plot_encounter_distributions.jl"))
const EP = EncounterDistributionPlots

@testset "Exact bin inclusion and empirical CDF" begin
    @test EP.histogram_counts([0., 0.5, 1., 2.], [0., 1., 2.]) == [2, 2]
    @test EP.histogram_counts([1., 10., 100.], [1., 10., 100.]) == [1, 2]
    @test_throws ErrorException EP.histogram_counts([3.], [0., 2.])
    @test_throws ErrorException EP.histogram_counts([NaN], [0., 2.])
    @test_throws ErrorException EP.histogram_counts([1.], [0., 0., 2.])
    @test_throws ErrorException EP.bin_edges([0.]; logarithmic=true)
    @test all(diff(EP.bin_edges([2., 2.])) .> 0)
    @test all(diff(EP.bin_edges([1., 1.]; logarithmic=true)) .> 0)
    cdf = EP.empirical_cdf([3., 1., 1., 2.])
    @test cdf.x == [1., 2., 3.]
    @test cdf.probability == [0.5, 0.75, 1.]
    @test_throws ErrorException EP.empirical_cdf(Float64[])
    edges = EP.bin_edges([0., 0., 10.]; logarithmic=true, include_zero=true)
    @test first(edges) == 0
    @test sum(EP.histogram_counts([0., 0., 10.], edges)) == 3
    @test first(EP.empirical_cdf([0., 0., 10.]).probability) == 2/3
end

@testset "Selected grid and plot data (per encounter=$per_encounter)" for per_encounter in (false, true)
    data = EP.load_data(EP.EV.OUTPUT)
    @test nrow(data.cases) == 162
    @test nrow(data.encounters) == 74130
    @test nrow(data.laser) == 72005
    @test count(iszero, data.cases.laser_interval_count) == 27
    specs = per_encounter ? EP.PER_ENCOUNTER_SPECS : EP.SPECS
    prepared = EP.plot_data(data; per_encounter)
    @test nrow(prepared.histograms) == (per_encounter ? 12 : 21) * EP.BIN_COUNT
    @test nrow(prepared.quantiles) == (per_encounter ? 12 : 21)
    @test length(prepared.cdfs) == (per_encounter ? 12 : 24)
    for group in groupby(prepared.histograms, [:quantity, :target_altitude_km, :metric, :scheduler])
        @test sum(group.count) == first(group.n)
        @test sum(group.probability) ≈ 1
        @test all(group.lower[2:end] .== group.upper[1:end-1])
    end
    for curve in prepared.cdfs
        @test issorted(curve.x)
        @test all(diff(curve.probability) .> 0)
        @test last(curve.probability) == 1
        spec = only(filter(s -> s.name == curve.quantity, specs))
        table = spec.kind == :encounter ? data.encounters : data.laser
        selected = filter(r -> r.target_altitude_km == curve.target_altitude_km &&
            (curve.scheduler == "pooled" || r.scheduler == curve.scheduler) &&
            (curve.population == "all_observed" || !(r.start_clipped || r.end_clipped)), table)
        values = selected[!, Symbol(curve.metric)]
        for index in (1, cld(length(curve.x), 2), length(curve.x))
            @test curve.probability[index] ≈ count(<=(curve.x[index]), values) / length(values)
        end
    end
    for row in eachrow(prepared.quantiles)
        table = !per_encounter && row.quantity == "laser_on_duration" ? data.laser : data.encounters
        selected = filter(r -> r.target_altitude_km == row.target_altitude_km &&
            (row.scheduler == "pooled" || r.scheduler == row.scheduler), table)
        values = selected[!, Symbol(row.metric)]
        @test row.n == length(values)
        @test row.median ≈ EP.quantile(values, 0.5)
        @test row.p95 ≈ EP.quantile(values, 0.95)
        if per_encounter
            @test row.n == sum(filter(r -> r.target_altitude_km == row.target_altitude_km,
                data.cases).encounter_count)
        end
    end
    laser_stats = filter(r -> r.quantity == "laser_on_duration", prepared.quantiles)
    if per_encounter
        @test sum(laser_stats.zero_events) == sum(data.cases.zero_laser_encounter_count)
    else
        @test Set(laser_stats.scheduler) == Set(["gve_inc", "gve_sma"])
        @test sum(laser_stats.n) == nrow(data.laser)
    end
    mktempdir() do output
        for spec in specs
            figure = EP.make_plot(spec, prepared, data.cases; per_encounter)
            @test length(figure.subplots) == 8
            @test figure[:plot_title] == spec.title
            @test all(panel -> xlims(panel) == xlims(figure.subplots[1]), figure.subplots[1:6])
            @test all(panel -> ylims(panel) == ylims(figure.subplots[1]), figure.subplots[1:3])
            @test all(panel -> ylims(panel) == (0, 1.03), figure.subplots[4:6])
            @test all(panel -> panel[:xaxis][:scale] == (spec.log && !spec.zero ? :log10 : :identity),
                figure.subplots[1:6])
            if per_encounter || spec.kind == :laser
                expected_series = per_encounter ? 1 : 2
                @test all(panel -> length(panel.series_list) == expected_series, figure.subplots[1:6])
            end
            if spec.kind == :laser
                @test occursin("schedulers separated", figure.subplots[7][:annotations][1][3].str)
            end
            for extension in ("png", "pdf")
                path = joinpath(output, "$(spec.name).$extension")
                savefig(figure, path)
                @test filesize(path) > 1000
                if extension == "png"
                    @test read(path, 8) == UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]
                else
                    @test String(read(path, 4)) == "%PDF"
                end
            end
        end
    end
end
