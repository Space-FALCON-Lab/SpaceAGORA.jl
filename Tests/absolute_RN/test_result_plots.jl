import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..",
    "9_SpaceAGORA.jl-main_full integration")))
using Test, Arrow, CSV, TOML
include(joinpath(@__DIR__, "..", "..", "Simulation_absolute_RN", "Code", "plot_all_results.jl"))
using Plots
const AP = AbsoluteRNPlots
const OE = AP.OE
const CONFIG = OE.StudyConfig

@testset "Campaign and CLI" begin
    @test OE.CASE_COUNT == 270
    @test OE.ALTITUDES == (1100, 1050, 1000, 950, 900)
    @test OE.COUNTS == (1, 100, 200, 300, 400, 500, 600, 700, 800)
    @test OE.parse_options([]).windows == [80, 100]
    @test_throws ErrorException OE.parse_options(["--orbits", "80,80"])
    @test_throws ErrorException OE.parse_options(["--orbits", "101"])
    @test_throws ErrorException OE.exact_sample_index([0.0, 1.0], 0.5)
    @test_throws ArgumentError OE.quantity_spec(:unknown)
    @test_throws ErrorException OE.load_grid(joinpath(@__DIR__, "nonexistent_plot_input"), [80])
    mktempdir() do input
        @test_throws ErrorException OE.load_grid(input, [80])
        settings = CONFIG.case_settings(CONFIG.StudyCase(1000, 0, 1, :gve_sma))
        for name in ("first", "duplicate")
            directory = joinpath(input, name)
            mkpath(directory)
            open(joinpath(directory, "run_info.toml"), "w") do io
                TOML.print(io, Dict("settings" => settings))
            end
        end
        @test_throws r"Ambiguous duplicate simulation" OE.load_grid(input, [80])
    end
end

@testset "Branch semantics and validation" begin
    zero = NamedTuple{OE.RTN_FIELDS}(ntuple(_ -> 0.0, length(OE.RTN_FIELDS)))
    values = merge(zero, (dv_R_mps=-1.0, dv_R_positive_mps=2.0, dv_R_negative_mps=-3.0,
        dv_T_mps=3.0, dv_T_positive_mps=4.0, dv_T_negative_mps=-1.0,
        dv_N_mps=-2.0, dv_N_positive_mps=5.0, dv_N_negative_mps=-7.0))
    batch = NamedTuple{CONFIG.TARGET_DIAGNOSTIC_COLUMNS}(
        Tuple([getproperty(values, name)] for name in OE.RTN_FIELDS))
    result = OE.rtn_values(batch, 1)
    @test result.dv_R_absolute_mps == 5.0
    @test result.dv_T_absolute_mps == 5.0
    @test result.dv_N_absolute_mps == 12.0
    @test result.dv_N_negative_mps == -7.0
    @test result.dv_N_absolute_mps != abs(result.dv_N_mps)
    @test OE.check_rtn(result, zero) == (0.0, 0.0, 0.0)
    @test_throws ErrorException OE.rtn_values(merge(batch, (target_dv_R_mps=[NaN],)), 1)
    @test_throws ErrorException OE.check_rtn(merge(values, (dv_R_positive_mps=-1.0,)), zero)
    @test_throws ErrorException OE.check_rtn(merge(values, (dv_N_negative_mps=1.0,)), zero)
    @test_throws ErrorException OE.check_rtn(merge(values, (dv_N_mps=0.0,)), zero)
    @test_throws ErrorException OE.check_rtn(zero, values)
    @test OE.check_rtn(zero, merge(zero, (dv_R_positive_mps=5e-9,))) == (5e-9, 0.0, 0.0)
end

@testset "Exact Feather endpoints and target indexing" begin
    for (altitude, inclination, helpers, policy) in
            ((900, 1, 800, :gve_sma), (1100, 0.5, 700, :gve_inc), (1000, 0, 1, :gve_inc))
        case = CONFIG.StudyCase(altitude, inclination, helpers, policy)
        directory = CONFIG.case_directory(CONFIG.DEFAULT_OUTPUT, case)
        info = TOML.parsefile(joinpath(directory, "run_info.toml"))
        @test info["settings"]["target_index"] == helpers + 1
        old_schema = deepcopy(info)
        old_schema["settings"]["schema_version"] = 2
        @test_throws ErrorException OE.read_case(directory, old_schema, [80]; include_rtn=true)
        bad_hash = deepcopy(info)
        bad_hash["checksums"]["trajectory.feather"] = "invalid"
        @test_throws ErrorException OE.read_case(directory, bad_hash, [80]; include_rtn=true)
        results = OE.read_case(directory, info, [80, 100]; include_rtn=true)
        oe_only = OE.read_case(directory, info, [80, 100])
        table = Arrow.Table(joinpath(directory, "trajectory.feather"))
        for (result, original, index) in zip(results, oe_only, (801, 1001))
            @test result.row == index
            @test result.time_s == table.time[index]
            @test all(name -> isequal(getproperty(result, name), getproperty(original, name)), keys(original))
            for (field, column) in zip(OE.RTN_FIELDS, CONFIG.TARGET_DIAGNOSTIC_COLUMNS)
                @test getproperty(result, field) == getproperty(table, column)[index]
            end
            for component in ("R", "T", "N")
                positive = getproperty(result, Symbol("dv_$(component)_positive_mps"))
                negative = getproperty(result, Symbol("dv_$(component)_negative_mps"))
                @test getproperty(result, Symbol("dv_$(component)_mps")) ≈ positive + negative atol=1e-8 rtol=0
                @test getproperty(result, Symbol("dv_$(component)_absolute_mps")) == positive - negative
            end
        end
    end
end

thresholds = [merge((target_altitude_km=altitude, delta_i_deg=inclination,),
    OE.continuous_access(CONFIG.EARTH_RADIUS_M + 1000e3,
        CONFIG.EARTH_RADIUS_M + altitude * 1000, CONFIG.RANGE_M, deg2rad(inclination)))
    for altitude in OE.ALTITUDES for inclination in OE.INCLINATIONS]

@testset "Continuous-access thresholds" begin
    @test length(thresholds) == 15
    for threshold in thresholds
        @test threshold.minimum_helpers isa Int
        @test threshold.minimum_helpers == ceil(Int, pi / threshold.alpha_min_rad)
        @test 2pi / threshold.minimum_helpers <= 2threshold.alpha_min_rad
        @test 2pi / (threshold.minimum_helpers - 1) > 2threshold.alpha_min_rad
    end
    inaccessible = OE.continuous_access(7000e3, 7100e3, 100e3, 0.0)
    @test ismissing(inaccessible.minimum_helpers)
    @test isinf(inaccessible.continuous_helpers)
end

function synthetic_row(policy, altitude, inclination, helpers, window)
    scale = helpers * window / 100 * (1 + (altitude - 900) / 200) *
        (policy == "gve_sma" ? 1.0 : 0.5)
    (scheduler=policy, target_altitude_km=altitude, delta_i_deg=inclination,
        helpers=helpers, orbits=window, time_s=Float64(window), row=10window + 1,
        source_trajectory="synthetic", da_m=scale, de=scale * 1e-5,
        di_deg=inclination * scale * 1e-3,
        draan_deg=inclination == 0 ? NaN : -scale * 0.001, raan_defined=inclination != 0,
        dv_R_mps=-scale, dv_T_mps=scale, dv_N_mps=-inclination * scale,
        dv_R_positive_mps=scale, dv_R_negative_mps=-2scale,
        dv_T_positive_mps=2scale, dv_T_negative_mps=-scale,
        dv_N_positive_mps=inclination * scale, dv_N_negative_mps=-2inclination * scale,
        dv_R_absolute_mps=3scale, dv_T_absolute_mps=3scale, dv_N_absolute_mps=3inclination * scale)
end

rows = [synthetic_row(policy, altitude, inclination, helpers, window)
    for window in (80, 100) for policy in OE.POLICIES for altitude in OE.ALTITUDES
    for inclination in OE.INCLINATIONS for helpers in OE.COUNTS]

@testset "Orbital-element bounds and undefined RAAN" begin
    for (field, values) in (
            (:da_m, [-715.0, 18420.0]), (:de, [1e-5, 0.00099]),
            (:di_deg, [-0.0026, 0.022]), (:draan_deg, [-45.2, NaN, -32.1]))
        lo, hi = OE.ylimits(values, field)
        finite = filter(isfinite, values)
        @test lo < hi
        @test lo <= minimum(finite) <= maximum(finite) <= hi
    end
    @test OE.ylimits([-45.2, NaN, -32.1], :draan_deg) ==
        OE.ylimits([-45.2, -32.1], :draan_deg)
    @test last(OE.ylimits([-45.2, NaN, -32.1], :draan_deg)) < 0
    for field in OE.quantity_spec(:oe).fields
        lo, hi = OE.ylimits([0.0, 0.0], field)
        @test lo < 0 < hi
    end
    @test_throws ErrorException OE.ylimits(Float64[], :da_m)
    @test_throws ErrorException OE.ylimits([NaN], :draan_deg)
    @test_throws ErrorException OE.ylimits([NaN, 1.0], :de)
    @test_throws ErrorException OE.ylimits([Inf, 1.0], :draan_deg)
end

@testset "Delta-v bounds enclose data without changing branch signs" begin
    @test OE.dv_ylimits([0.0, 7.24], :dv_R_positive_mps) == (0.0, 8.0)
    @test OE.dv_ylimits([-6.0, 8.0], :dv_R_mps) == (-6.0, 8.0)
    for (field, values) in (
            (:dv_R_mps, [-6.290155758434805, 7.231988120548606]),
            (:dv_T_mps, [-0.3551522586733772, 9.240851720180482]),
            (:dv_N_mps, [-0.6989699463815524, 0.17279314407506108]),
            (:dv_N_absolute_mps, [0.0, 5.26840155809547]),
            (:dv_T_negative_mps, [-1.93637628376923, 0.0]),
            (:dv_R_positive_mps, [-5e-9, 2.0]),
            (:dv_N_negative_mps, [-2.0, 5e-9]))
        lo, hi = OE.dv_ylimits(values, field)
        @test lo < hi
        @test lo <= minimum(values) <= maximum(values) <= hi
        @test lo <= 0 <= hi
    end
    for component in ("R", "T", "N")
        positive = OE.dv_ylimits([0.0], Symbol("dv_$(component)_positive_mps"))
        negative = OE.dv_ylimits([0.0], Symbol("dv_$(component)_negative_mps"))
        absolute = OE.dv_ylimits([0.0], Symbol("dv_$(component)_absolute_mps"))
        net = OE.dv_ylimits([0.0], Symbol("dv_$(component)_mps"))
        @test positive[1] == 0 < positive[2]
        @test negative[1] < 0 == negative[2]
        @test absolute[1] == 0 < absolute[2]
        @test net[1] < 0 < net[2]
    end
    @test_throws ErrorException OE.dv_ylimits(Float64[], :dv_R_mps)
    @test_throws ErrorException OE.dv_ylimits([NaN, 1.0], :dv_R_mps)
    @test_throws ErrorException OE.dv_ylimits([0.0, Inf], :dv_R_mps)
end

@testset "Grid layouts, labels, colors, curves and thresholds" begin
    gr()
    mktempdir() do output
        seen_limits = Dict{Symbol,Tuple{Float64,Float64}}()
        for quantity in (:oe, :rtn, :absolute, :split_n, :branches), window in (80, 100)
            spec = OE.quantity_spec(quantity)
            figure = OE.make_plot(rows, thresholds, "combined", window, output; quantity)
            @test length(figure.subplots) == 5length(spec.fields) + 1
            for (index, panel) in enumerate(figure.subplots[1:end-1])
                r, c = divrem(index - 1, 5) .+ 1
                field = spec.fields[r]
                curves = filter(s -> s[:seriestype] == :path, panel.series_list)
                @test length(curves) == (field == :draan_deg ? 4 : 6)
                @test panel[:xaxis][:grid] && panel[:yaxis][:grid]
                @test panel[:yaxis][:guide] == (c == 1 ? spec.labels[r] : "")
                @test panel[:xaxis][:guide] == (r == length(spec.fields) ? "N_helpers" : "")
                if quantity != :oe
                    limits = ylims(panel)
                    values = [getproperty(row, field) for row in rows]
                    @test limits == OE.dv_ylimits(values, field)
                    @test limits[1] <= minimum(values) <= maximum(values) <= limits[2]
                    @test limits == get!(seen_limits, field, limits)
                    endswith(string(field), "_negative_mps") && @test limits[2] == 0
                    if endswith(string(field), "_positive_mps") || endswith(string(field), "_absolute_mps")
                        @test limits[1] == 0
                    end
                else
                    values = [getproperty(row, field) for row in rows]
                    finite = filter(isfinite, values)
                    limits = ylims(panel)
                    @test limits == OE.ylimits(values, field)
                    @test limits[1] <= minimum(finite) <= maximum(finite) <= limits[2]
                    @test limits == get!(seen_limits, field, limits)
                end
                series_index = 0
                for policy in ("gve_sma", "gve_inc"), (k, inclination) in enumerate(OE.INCLINATIONS)
                    field == :draan_deg && inclination == 0 && continue
                    series_index += 1
                    series = curves[series_index]
                    expected = sort(filter(row -> row.orbits == window &&
                        row.scheduler == policy && row.delta_i_deg == inclination &&
                        row.target_altitude_km == OE.ALTITUDES[c], rows); by=row -> row.helpers)
                    @test collect(series[:x]) == collect(OE.COUNTS)
                    @test collect(series[:y]) == [getproperty(row, field) for row in expected]
                    @test series[:linecolor] == Plots.plot_color(getproperty(OE.COMBINED_COLORS, Symbol(policy))[k])
                    @test series[:linestyle] == (policy == "gve_inc" ? :dot : :solid)
                end
                lines = filter(s -> s[:seriestype] == :straightline, panel.series_list)
                @test length(lines) == 3
                for (line, inclination) in zip(lines, OE.INCLINATIONS)
                    threshold = only(filter(t -> t.target_altitude_km == OE.ALTITUDES[c] &&
                        t.delta_i_deg == inclination, thresholds))
                    @test all(==(threshold.minimum_helpers), line[:x])
                    @test line[:linestyle] == :dash
                end
            end
            for extension in ("png", "pdf")
                @test filesize(joinpath(output,
                    "$(spec.prefix)_from_feather_combined_e0.0000_nu1.0000_$(window)orbits_h_cols.$extension")) > 10_000
            end
        end
        for policy in OE.POLICIES, window in (80, 100)
            single = OE.make_plot(rows, thresholds, policy, window, output)
            @test length(single.subplots) == 21
            @test length(filter(s -> s[:seriestype] == :path, single.subplots[1].series_list)) == 3
            for (index, panel) in enumerate(single.subplots[1:end-1])
                field = OE.quantity_spec(:oe).fields[div(index - 1, 5) + 1]
                @test ylims(panel) == seen_limits[field]
            end
        end
        @test_throws ErrorException OE.make_plot(rows[2:end], thresholds, "combined", 80, output)
        for window in (80, 100), xfield in (:helpers, :target_altitude_km)
            figure = AP.PooledOrbitalElementPlots.make_plot(rows, window, xfield)
            @test length(figure.subplots) == 4
            for (index, panel) in enumerate(figure.subplots)
                field = OE.quantity_spec(:oe).fields[index]
                @test ylims(panel) == seen_limits[field]
                points = filter(s -> s[:seriestype] == :scatter, panel.series_list)
                @test sum(length(s[:x]) for s in points) == (index == 4 ? 180 : 270)
                @test panel[:xaxis][:grid] && panel[:yaxis][:grid]
                if index == 4
                    @test panel[:title] == "Unwrapped RAAN; equatorial cases omitted (undefined)"
                    @test isempty(panel[:annotations])
                end
            end
        end
    end
end

@testset "Endpoint CSV shape and signed values" begin
    mktempdir() do output
        for window in (80, 100), family in AP.RTN_FAMILIES
            fields = OE.quantity_spec(family.quantity).fields
            path = joinpath(output, "endpoints.csv")
            AP.write_endpoints(path, rows, window, fields)
            records = CSV.File(path)
            @test length(records) == 270
            @test propertynames(records) == collect((AP.IDENTIFIERS..., fields...))
            @test all(row -> row.row == window * 10 + 1, records)
            @test all(row -> row.orbits == window, records)
            for (record, original) in zip(records, filter(row -> row.orbits == window, rows))
                @test all(field -> getproperty(record, field) == getproperty(original, field), fields)
            end
        end
        @test_throws ErrorException AP.write_endpoints(joinpath(output, "incomplete.csv"),
            rows[2:end], 80, OE.quantity_spec(:rtn).fields)
    end
end

@testset "Selected-altitude folders, ascending columns and CSVs" begin
    @test AP.SELECTED_ALTITUDES == (1000, 1050, 1100)
    @test [f.quantity for f in AP.SELECTED_FAMILIES] == [:oe, :absolute, :split_n, :branches]
    @test [first(basename(f.folder)) for f in AP.SELECTED_FAMILIES] == ['4', '5', '6', '7']
    mktempdir() do output
        AP.render_selected_altitudes(rows, thresholds, [80, 100], output)
        for family in AP.SELECTED_FAMILIES
            directory = joinpath(output, family.folder)
            spec = OE.quantity_spec(family.quantity)
            threshold_csv = CSV.File(joinpath(directory, "continuous_access_thresholds.csv"))
            @test length(threshold_csv) == 9
            @test Set(t.target_altitude_km for t in threshold_csv) == Set(AP.SELECTED_ALTITUDES)
            for window in (80, 100)
                records = CSV.File(joinpath(directory, "$(family.csv)_$(window)orbits.csv"))
                expected = filter(row -> row.orbits == window &&
                    row.target_altitude_km in AP.SELECTED_ALTITUDES, rows)
                @test length(records) == length(expected) == 162
                @test Set(r.target_altitude_km for r in records) == Set(AP.SELECTED_ALTITUDES)
                for (record, original) in zip(records, expected)
                    @test record.row == window * 10 + 1
                    @test all(field -> isequal(getproperty(record, field), getproperty(original, field)),
                        (AP.IDENTIFIERS..., spec.fields...))
                end
                for extension in ("png", "pdf")
                    @test filesize(joinpath(directory,
                        "$(spec.prefix)_from_feather_combined_e0.0000_nu1.0000_$(window)orbits_h_cols.$extension")) > 10_000
                end
            end
            figure = OE.make_plot(rows, thresholds, "combined", 80, directory;
                quantity=family.quantity, altitudes=AP.SELECTED_ALTITUDES)
            @test length(figure.subplots) == length(spec.fields) * 3 + 1
            for (r, field) in enumerate(spec.fields), (c, altitude) in enumerate(AP.SELECTED_ALTITUDES)
                panel = figure.subplots[(r - 1) * 3 + c]
                @test panel[:title] == (r == 1 ? "h_t = $altitude km" : "")
                limit_function = family.quantity == :oe ? OE.ylimits : OE.dv_ylimits
                @test ylims(panel) == limit_function([getproperty(row, field) for row in rows], field)
                curves = filter(s -> s[:seriestype] == :path, panel.series_list)
                @test length(curves) == (field == :draan_deg ? 4 : 6)
                index = 0
                for policy in ("gve_sma", "gve_inc"), (k, inclination) in enumerate(OE.INCLINATIONS)
                    field == :draan_deg && inclination == 0 && continue
                    index += 1
                    curve = curves[index]
                    expected = sort(filter(row -> row.orbits == 80 &&
                        row.target_altitude_km == altitude && row.scheduler == policy &&
                        row.delta_i_deg == inclination, rows); by=row -> row.helpers)
                    @test collect(curve[:y]) == [getproperty(row, field) for row in expected]
                    @test collect(curve[:x]) == collect(OE.COUNTS)
                    @test curve[:linecolor] == Plots.plot_color(getproperty(OE.COMBINED_COLORS, Symbol(policy))[k])
                    @test curve[:linestyle] == (policy == "gve_inc" ? :dot : :solid)
                end
                lines = filter(s -> s[:seriestype] == :straightline, panel.series_list)
                @test length(lines) == 3
                for (line, inclination) in zip(lines, OE.INCLINATIONS)
                    expected = only(filter(t -> t.target_altitude_km == altitude &&
                        t.delta_i_deg == inclination, thresholds))
                    @test all(==(expected.minimum_helpers), line[:x])
                end
            end
        end
        for altitudes in ((), (1000, 1000), (850,))
            @test_throws ArgumentError OE.make_plot(rows, thresholds, "combined", 80, output; altitudes)
        end
        @test_throws ErrorException AP.write_endpoints(joinpath(output, "incomplete.csv"),
            rows[2:end], 80, OE.quantity_spec(:rtn).fields; altitudes=AP.SELECTED_ALTITUDES)
        @test_throws ErrorException AP.render_selected_altitudes(rows, thresholds[2:end], [80], output)
    end
end
