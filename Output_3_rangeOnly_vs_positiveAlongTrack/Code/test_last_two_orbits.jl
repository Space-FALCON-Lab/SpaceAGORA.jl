include("ComparisonResults.jl")
using .ComparisonResults
using Arrow, CSV, DataFrames, LinearAlgebra, StaticArrays, Test, TOML

const CR = ComparisonResults
const SC = CR.SC
const ROOTS = (CR.CampaignPaths.OUTPUT, CR.CampaignPaths.NO_LASER_OUTPUT,
    CR.CampaignPaths.DIFFERENCE_OUTPUT)

function fixture(case)
    time_s = SC.trajectory_times(case)
    x = collect(range(0.0, case.orbits; length=length(time_s)))
    scale = case.policy === :no_laser ? 0.0 : case.helpers / 100
    elements = [(a=7e6 + 10cospi(2orbit) + scale * orbit,
        e=0.001 + 1e-6scale * orbit, i=deg2rad(case.inclination_deg) + 1e-7scale * orbit)
        for orbit in x]
    angles = [(raan=deg2rad(170 + 0.2orbit + 0.001scale * orbit),) for orbit in x]
    raan = case.inclination_deg == 0 ? fill(NaN, length(x)) : CR.unwrap_raan(angles)
    positions = repeat([7e6, 0.0, 0.0], 1, length(x))
    velocities = reduce(hcat, [SVector(scale * orbit, 7500 - scale * orbit, 2scale * orbit) for orbit in x])
    rtn = reduce(hcat, [[(occursin("negative", string(field)) ? -1 :
        occursin("positive", string(field)) ? 2 :
        occursin("absolute", string(field)) ? 3 : 1) * scale * orbit
        for field in CR.RTN_FIELDS] for orbit in x])
    info = Dict("run_id" => SC.case_signature(case), "settings" => SC.case_settings(case),
        "checksums" => Dict("trajectory.feather" => "synthetic fixture only"),
        "initial_states" => [], "source_hashes" => Dict("fixture" => "same source"),
        "firing_intervals" => 0)
    (; case, time_s, elements, raan, positions, velocities, rtn, info,
        source_trajectory="synthetic fixture only")
end

function unit_tests()
    @testset "Trapezoidal time-weighted mean" begin
        @test CR.time_average([0.0, 1.0, 4.0], [3.0, 3.0, 3.0]) == 3.0
        @test CR.time_average([0.0, 1.0, 4.0], [0.0, 1.0, 4.0]) == 2.0
        @test CR.time_average([0.0, 1.0, 4.0], [0.0, 1.0, 16.0]) == 6.5
        @test isnan(CR.time_average([0.0, 1.0], [NaN, NaN]))
        @test_throws ErrorException CR.time_average([0.0], [1.0])
        @test_throws ErrorException CR.time_average([0.0, 1.0], [1.0])
        @test_throws ErrorException CR.time_average([0.0, 0.0], [1.0, 1.0])
        @test_throws ErrorException CR.time_average([0.0, -1.0], [1.0, 1.0])
        @test_throws ErrorException CR.time_average([0.0, Inf], [1.0, 1.0])
        @test_throws ErrorException CR.time_average([0.0, 1.0], [1.0, NaN])
        @test_throws ErrorException CR.time_average([0.0, 1.0], [1.0, Inf])
    end
    @testset "Window boundaries, changes from initial state and paired differences" begin
        laser = fixture(SC.StudyCase(1050, 0.5, 100, :naive_next_entering))
        baseline = fixture(SC.StudyCase(1050, 0.5, 100, :no_laser))
        on = CR.orbit_averages(laser, [80, 100])
        off = CR.orbit_averages(baseline, [80, 100])
        difference = CR.orbit_averages(laser, [80, 100]; baseline)
        for (window, row, zero_row, diff_row) in zip((80, 100), on, off, difference)
            @test row.sample_count == 21
            @test row.window_start_row == 10(window - 2) + 1
            @test row.window_end_row == row.row == 10window + 1
            @test row.window_start_orbits == window - 2
            @test row.window_end_orbits == window
            @test row.window_start_time_s == laser.time_s[row.window_start_row]
            @test row.window_end_time_s == row.time_s == laser.time_s[row.row]
            @test row.da_m ≈ window - 1 - 10 atol=1e-8
            @test zero_row.da_m ≈ -10 atol=1e-8
            @test all(field -> getproperty(zero_row, field) == 0, CR.RTN_FIELDS)
            @test row.dv_T_mps ≈ window - 1
            @test row.dv_N_negative_mps ≈ -(window - 1)
            for field in (CR.OE_FIELDS..., CR.RTN_FIELDS...)
                @test getproperty(diff_row, field) ≈
                    getproperty(row, field) - getproperty(zero_row, field) atol=1e-9
            end
            @test diff_row.velocity_difference_R_mps ≈ window - 1
            @test diff_row.velocity_difference_T_mps ≈ -(window - 1)
            @test diff_row.velocity_difference_N_mps ≈ 2(window - 1)
            @test diff_row.velocity_difference_norm_mps ≈ sqrt(6) * (window - 1)
        end
        @test only(CR.orbit_averages(laser, [2])).window_start_row == 1
        @test_throws ErrorException CR.orbit_averages(laser, [])
        @test_throws ErrorException CR.orbit_averages(laser, [1])
        @test_throws ErrorException CR.orbit_averages(laser, [101])
        @test_throws ErrorException CR.orbit_averages(laser, [80, 80])
        @test_throws ErrorException CR.orbit_averages(laser, [80.01])
        equatorial = fixture(SC.StudyCase(1000, 0, 1, :no_laser))
        @test all(row -> isnan(row.draan_deg) && !row.raan_defined,
            CR.orbit_averages(equatorial, [80, 100]))
        @test CR.plot_spec(:oe, 100).stem == CR.PP.SPECS.oe.stem
        @test CR.plot_spec(:dv, 80).csv == "cumulative_dv_RT_N_branches_80orbits.csv"
        @test CR.plot_spec(:velocity, 80; difference=true).stem ==
            "velocity_difference_RTN_laser_minus_no_laser_80orbits_h_cols"
        @test startswith(CR.plot_spec(:oe, 80; averaged=true).stem, "mean_last2orbits_")
    end
end

function independently_averaged(row, fields; difference=false)
    first_row, last_row = row.window_start_row, row.window_end_row
    table = Arrow.Table(difference ? row.source_difference_trajectory : row.source_trajectory)
    time = collect(table[difference ? :time_s : :time])[first_row:last_row]
    widths = diff(time)
    weights = vcat(widths[1], widths[1:end-1] + widths[2:end], widths[end]) /
        (2(last(time) - first(time)))
    expected = Dict{Symbol,Float64}()
    target = Int(row.helpers) + 1
    function oe(index)
        r = SVector{3}(ntuple(k -> table[Symbol("sc$(target)_pos_$k")][index], 3))
        v = SVector{3}(ntuple(k -> table[Symbol("sc$(target)_vel_$k")][index], 3))
        CR.OE.orbital_elements(r, v, SC.MU)
    end
    initial = difference ? nothing : oe(1)
    for field in fields
        values = if difference
            collect(table[field])[first_row:last_row]
        elseif field in CR.OE_FIELDS
            if field == :draan_deg
                # Unwrap the full history before selecting the averaging interval.
                angle = CR.unwrap_raan([oe(k) for k in 1:last_row])
                rad2deg.(angle[first_row:last_row] .- angle[1])
            else
                name = field == :da_m ? :a : field == :de ? :e : :i
                scale = field == :di_deg ? 180 / pi : 1.0
                [scale * (getproperty(oe(k), name) - getproperty(initial, name))
                    for k in first_row:last_row]
            end
        elseif occursin("absolute", string(field))
            positive = Symbol("target_" * replace(string(field), "absolute" => "positive"))
            negative = Symbol("target_" * replace(string(field), "absolute" => "negative"))
            [table[positive][k] - table[negative][k] for k in first_row:last_row]
        else
            collect(table[Symbol("target_$field")])[first_row:last_row]
        end
        expected[field] = dot(weights, values)
    end
    expected
end

function saved_tests()
    @testset "All saved averages and exact plotted curves" begin
        for (root_index, root) in enumerate(ROOTS)
            difference = root_index == 3
            rows = CSV.read(joinpath(root, "last_two_orbits_average_summary.csv"), DataFrame)
            @test nrow(rows) == (root_index == 2 ? 162 : 324)
            sources = CSV.read(joinpath(root, "production_endpoint_summary.csv"), DataFrame)
            keys(row) = (row.scheduler, row.target_altitude_km, row.delta_i_deg, row.helpers, row.orbits)
            lookup = Dict(keys(row) => row for row in eachrow(sources))
            fields = (CR.OE_FIELDS..., CR.RTN_FIELDS..., (difference ? CR.VELOCITY_FIELDS : ())...)
            for row in eachrow(rows)
                @test row.aggregation == "time_weighted_mean" && row.sample_count == 21
                @test row.window_start_orbits == row.orbits - 2
                original = lookup[keys(row)]
                reference = merge(NamedTuple(original), (window_start_row=row.window_start_row,
                    window_end_row=row.window_end_row))
                expected = independently_averaged(reference, fields; difference)
                for field in fields
                    actual = row[field]
                    @test isnan(expected[field]) ? isnan(actual) :
                        isapprox(actual, expected[field]; atol=1e-9, rtol=1e-12)
                end
                root_index == 2 && @test all(field -> row[field] == 0, CR.RTN_FIELDS)
            end
            for window in (80, 100)
                directory = joinpath(root, "images", "last_two_orbits_average", "$(window)orbits")
                manifest = TOML.parsefile(joinpath(directory, "figure_manifest.toml"))
                @test manifest["aggregation"] == "time_weighted_mean"
                @test manifest["sample_count"] == [21]
                for section in ("figure_checksums", "csv_checksums"), (file, digest) in manifest[section]
                    @test SC.filehash(joinpath(directory, file)) == digest
                end
                for kind in (difference ? (:oe, :dv, :velocity) : (:oe, :dv))
                    spec = CR.plot_spec(kind, window; difference, averaged=true)
                    data = CSV.read(joinpath(directory, spec.csv), DataFrame)
                    selected = filter(row -> row.orbits == window, rows)
                    @test nrow(data) == nrow(selected)
                    for column in propertynames(data)
                        @test isequal(data[!, column], selected[!, column])
                    end
                    mktempdir() do output
                        figure = CR.PP.make_plot(kind == :oe ? :oe : :dv; directory, output,
                            policies=Tuple(unique(data.scheduler)), policy_labels=CR.POLICY_LABELS,
                            colors=CR.COLORS, styles=CR.STYLES, specification=spec,
                            orbit_window=window, endpoint_row=10window + 1, averaging_orbits=2,
                            source_limits=CR.plot_limits(eachrow(selected), spec))
                        @test occursin("$(window-2)-$window", figure[:plot_title])
                        for (r, field) in enumerate(spec.fields), (c, altitude) in enumerate(CR.PP.ALTITUDES)
                            panel = figure.subplots[(r-1)*3+c]
                            curves = filter(series -> series[:primary], panel.series_list)
                            expected_count = length(unique(data.scheduler)) *
                                (field == :draan_deg ? 2 : 3)
                            @test length(curves) >= expected_count
                            ordinal = 0
                            for policy in unique(data.scheduler), inclination in CR.PP.INCLINATIONS
                                field == :draan_deg && inclination == 0 && continue
                                ordinal += 1
                                values = sort(filter(row -> row.scheduler == policy &&
                                    row.target_altitude_km == altitude && row.delta_i_deg == inclination, data),
                                    :helpers)
                                @test collect(curves[ordinal][:x]) == values.helpers
                                @test collect(curves[ordinal][:y]) == values[!, field] .* spec.scales[r]
                            end
                        end
                    end
                end
            end
            GC.gc()
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 && only(ARGS) in ("unit", "saved") ||
        error("Usage: test_last_two_orbits.jl unit|saved")
    only(ARGS) == "unit" ? unit_tests() : saved_tests()
end
