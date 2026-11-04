import Pkg
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(ROOT, "9_SpaceAGORA.jl-main_full integration"))
using Test, CSV, DataFrames, Markdown
include(joinpath(ROOT, "Simulation_absolute_RN", "Code", "write_encounter_parameter_pages.jl"))
const PP = EncounterParameterPages

@testset "Four parameter pages from all selected events" begin
    data = PP.EP.load_data(PP.EP.EV.OUTPUT)
    rows = PP.summarize(data)
    @test nrow(rows) == 72
    @test all(==(9), rows.cases)
    @test allunique(zip(rows.page, rows.scheduler, rows.target_altitude_km, rows.target_inclination_deg))
    for row in eachrow(rows)
        spec = PP.PAGES[row.page]
        table = spec.kind == :laser ? data.laser : data.encounters
        values = [getproperty(r, spec.field) for r in eachrow(table) if
            r.scheduler == row.scheduler && r.target_altitude_km == row.target_altitude_km &&
            r.target_inclination_deg == row.target_inclination_deg]
        @test row.events == length(values)
        @test isequal(row.min, isempty(values) ? missing : minimum(values))
        @test isequal(row.max, isempty(values) ? missing : maximum(values))
        @test ismissing(row.avg) ? isempty(values) : row.avg ≈ sum(values) / length(values)
    end
    for page in 1:4
        selected = filter(r -> r.page == page, rows)
        @test sum(selected.events) == (page == 4 ? 72005 : 74130)
    end
    empty = filter(r -> r.events == 0, rows)
    @test nrow(empty) == 3
    @test all(r -> r.page == 4 && r.scheduler == "gve_inc" && r.target_inclination_deg == 0, eachrow(empty))
    @test sum(empty.cases_without_events) == 27
    @test PP.number(missing) == "N/A"
    @test PP.number(-12.34) == "-12.34"
    mktempdir() do output
        CSV.write(joinpath(output, "parameter_pages_min_max.csv"), rows)
        text = read(PP.write_report(output), String)
        @test length(collect(eachmatch(r"# Page \d\.", text))) == 4
        @test length(collect(eachmatch(r"(?m)^---$", text))) == 3
        @test !occursin(r"<[^>]+>", text)
        @test length(collect(eachmatch(r"(?m)^!\[", text))) == 4
        @test occursin("arithmetic event-weighted mean", text)
        @test length(collect(eachmatch(r"## Conclusion", text))) == 4
        rendered = sprint(show, MIME"text/html"(), Markdown.parse(text))
        @test length(collect(eachmatch(r"<img\b", rendered))) == 4
        @test !occursin("<pre>", rendered)
        @test !occursin("&lt;table", rendered)
        pages = split(text, r"# Page \d\.")[2:end]
        for (page, section) in enumerate(pages)
            @test occursin("parameter_page_$(page)_table.png", section)
            for extension in ("png", "pdf")
                @test filesize(joinpath(output, "parameter_page_$(page)_table.$extension")) > 1000
            end
            selected = filter(r -> r.page == page, rows)
            figure = PP.combined_table(selected, PP.PAGES[page])
            annotations = figure.subplots[1][:annotations]
            for (s, scheduler) in enumerate(PP.SCHEDULERS),
                    (i, inclination) in enumerate(PP.INCLINATIONS),
                    (h, altitude) in enumerate(PP.EP.ALTITUDES)
                row = only(eachrow(filter(r -> r.scheduler == scheduler &&
                    r.target_inclination_deg == inclination && r.target_altitude_km == altitude, selected)))
                for (k, field) in enumerate((:min, :avg, :max))
                    x, y = 2.5 + 3(h - 1) + k - 1, 5.5 - 3(s - 1) - (i - 1)
                    annotation = only(filter(a -> a[1] == x && a[2] == y, annotations))
                    @test annotation[3].str == PP.number(getproperty(row, field))
                end
            end
        end
        rows.max[1] = 987654.0
        CSV.write(joinpath(output, "parameter_pages_min_max.csv"), rows)
        @test occursin("987654", read(PP.write_report(output), String))
    end
end
