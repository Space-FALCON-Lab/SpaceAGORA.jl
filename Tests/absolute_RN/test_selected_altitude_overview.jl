import Pkg
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(ROOT, "9_SpaceAGORA.jl-main_full integration"))
using Test, CSV, DataFrames, Markdown
include(joinpath(ROOT, "Simulation_absolute_RN", "Code", "write_selected_altitude_overview.jl"))
const SO = SelectedAltitudeOverview

@testset "All-scenario selected-altitude summary" begin
    data = SO.EP.load_data(SO.EP.EV.OUTPUT)
    rows = SO.summarize(data)
    combined = SO.summarize_all_altitudes(data)
    @test nrow(data.cases) == 162
    @test nrow(rows) == 15
    @test nrow(combined) == 5
    @test allunique(zip(rows.target_altitude_km, rows.metric))
    for row in eachrow(rows)
        values = SO.metric_values(data, row.target_altitude_km, Symbol(row.metric))
        @test row.events == length(values)
        @test row.cases == (row.scheduler == "pooled" ? 54 : 27)
        @test row.min == minimum(values)
        @test row.average ≈ sum(values) / length(values)
        @test row.max == maximum(values)
    end
    for row in eachrow(combined)
        values = vcat([SO.metric_values(data, altitude, Symbol(row.metric))
            for altitude in SO.ALTITUDES]...)
        @test row.events == length(values)
        @test row.min == minimum(values)
        @test row.average ≈ sum(values) / length(values)
        @test row.max == maximum(values)
        @test row.cases == (row.scheduler == "pooled" ? 162 : 81)
    end
    @test only(filter(r -> r.metric == "entry_speed", combined)).events == 74130
    @test only(filter(r -> r.metric == "laser_sma", combined)).events == 39328
    @test only(filter(r -> r.metric == "laser_inc", combined)).events == 32677
    @test SO.presentation_number("entry_speed", 245.524406) == "245.52"
    @test SO.presentation_number("entry_range_rate", -105.00506) == "-105.01"
    @test SO.presentation_number("encounter_duration", 4261.160985) == "4261.16"
    @test SO.presentation_number("laser_sma", 1035.79396) == "1035.8"
    @test SO.presentation_number("laser_inc", 1.936267) == "1.9363"
    @test SO.presentation_number("laser_inc", missing) == "N/A"
    @test SO.full_precision(245.52440607433067) == string(245.52440607433067)
    @test SO.full_precision(missing) == "N/A"
    for altitude in SO.ALTITUDES
        subset = filter(r -> r.target_altitude_km == altitude, rows)
        @test only(filter(r -> r.metric == "entry_speed", subset)).events ==
            sum(filter(r -> r.target_altitude_km == altitude, data.cases).encounter_count)
        @test only(filter(r -> r.metric == "entry_range_rate", subset)).events ==
            only(filter(r -> r.metric == "encounter_duration", subset)).events
        for scheduler in ("gve_sma", "gve_inc")
            laser = only(filter(r -> r.metric == "laser_$(scheduler == "gve_sma" ? "sma" : "inc")", subset))
            source = filter(r -> r.target_altitude_km == altitude && r.scheduler == scheduler, data.laser)
            @test laser.events == nrow(source)
            @test laser.clipped_events == count(source.start_clipped .| source.end_clipped)
        end
    end
    mktempdir() do output
        CSV.write(joinpath(output, "selected_altitude_statistics.csv"), rows)
        CSV.write(joinpath(output, "selected_altitude_combined_statistics.csv"), combined)
        report = read(SO.write_report(output, SO.OUTPUT), String)
        @test length(collect(eachmatch(r"(?m)^## Target altitude:", report))) == 3
        @test length(collect(eachmatch(r"(?m)^## Full-precision combined summary:", report))) == 1
        @test length(collect(eachmatch(r"(?m)^\| Statistic \|", report))) == 5
        @test length(collect(eachmatch(r"(?m)^\| (Max|Avg|Min) \|", report))) == 15
        @test length(collect(eachmatch(r"(?m)^\|:---\|", report))) == 5
        @test occursin("event-weighted", report)
        @test occursin("### Presentation-rounded combined table", report)
        @test occursin("5 significant figures** for encounter speed, signed range rate (assumed)", report)
        @test occursin("Boundary flags", report)
        full_table = first(split(report, "### Presentation-rounded combined table"))
        for (statistic, field) in (("Max", :max), ("Avg", :average), ("Min", :min))
            line = only(filter(line -> startswith(line, "| $statistic |"),
                split(full_table, '\n')))
            cells = strip.(split(strip(line), '|'))
            @test length(cells) == length(SO.DISPLAY) + 3
            for (index, (metric, _, _)) in enumerate(SO.DISPLAY)
                row = only(filter(r -> r.metric == metric, combined))
                @test cells[index + 2] == SO.full_precision(getproperty(row, field))
            end
        end
        html = sprint(show, MIME"text/html"(), Markdown.parse(report))
        @test length(collect(eachmatch(r"<table>", html))) == 5
        @test !occursin("<pre>", html)
        selected = filter(r -> r.target_altitude_km == 1000, rows)
        for (metric, label, _) in SO.DISPLAY
            row = only(filter(r -> r.metric == metric, selected))
            @test occursin(label, report)
            @test occursin(SO.number(row.average), report)
        end
        rows.average[1] = 987654.0
        CSV.write(joinpath(output, "selected_altitude_statistics.csv"), rows)
        combined.average[1] = 876543.0
        CSV.write(joinpath(output, "selected_altitude_combined_statistics.csv"), combined)
        regenerated = SO.write_report(output, SO.OUTPUT)
        @test occursin("987654", read(regenerated, String))
        @test occursin("876543", read(regenerated, String))
    end
end
