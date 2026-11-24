#!/usr/bin/env julia

if abspath(PROGRAM_FILE) == @__FILE__
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..",
        "9_SpaceAGORA.jl-main_full integration")))
end

module EncounterParameterPages

include("plot_encounter_distributions.jl")
using .EncounterDistributionPlots
using CSV, DataFrames, Printf, Statistics, Plots

const EP = EncounterDistributionPlots
const OUTPUT = EP.EV.OUTPUT
const INCLINATIONS = (0.0, 0.5, 1.0)
const SCHEDULERS = ("gve_sma", "gve_inc")
const PAGES = (
    (title="Encounter speed at entry", field=:entry_relative_speed_mps, kind=:encounter,
        unit="m/s", definition="One relative-speed value per geometric encounter, calculated from the saved entry velocities. This is not the minimum or maximum speed during the encounter."),
    (title="Signed range rate at entry", field=:entry_range_rate_mps, kind=:encounter,
        unit="m/s", definition="One signed range-rate value per geometric encounter at its saved entry state: negative means approaching, positive means receding. This is not the within-encounter range-rate envelope."),
    (title="Encounter duration", field=:duration_s, kind=:encounter,
        unit="s", definition="Observed time from entry to exit of the 200 km geometric encounter range, with one duration per target-helper encounter."),
    (title="Uninterrupted laser-on duration", field=:duration_s, kind=:laser,
        unit="s", definition="Duration of each uninterrupted firing interval, including held out-of-range firing. This matches the scheduler-separated laser plot, not summed laser-on time per geometric encounter."),
)

function summarize(data)
    rows = NamedTuple[]
    for (page, spec) in enumerate(PAGES), scheduler in SCHEDULERS,
            inclination in INCLINATIONS, altitude in EP.ALTITUDES
        cases = filter(r -> r.scheduler == scheduler && r.target_inclination_deg == inclination &&
            r.target_altitude_km == altitude, data.cases)
        nrow(cases) == 9 || error("Expected all nine helper counts in each table cell")
        table = spec.kind == :encounter ? data.encounters : data.laser
        ids = Set(cases.run_id)
        events = filter(r -> r.run_id in ids, table)
        values = events[!, spec.field]
        all(isfinite, values) || error("Nonfinite $(spec.field)")
        counts = spec.kind == :encounter ? cases.encounter_count : cases.laser_interval_count
        length(values) == sum(counts) || error("Cell count mismatch")
        push!(rows, (page=page, parameter=spec.title, metric=string(spec.field), unit=spec.unit,
            scheduler=scheduler, target_altitude_km=altitude, target_inclination_deg=inclination,
            cases=nrow(cases), cases_without_events=count(iszero, counts), events=length(values),
            min=isempty(values) ? missing : minimum(values),
            avg=isempty(values) ? missing : mean(values),
            max=isempty(values) ? missing : maximum(values),
            clipped_events=count(events.start_clipped .| events.end_clipped),
            full_run_events=count(events.start_clipped .& events.end_clipped),
            positive_entry_events=page == 2 ? count(>(0), values) : 0,
            positive_entry_start_clipped=page == 2 ?
                count((values .> 0) .& events.start_clipped) : 0))
    end
    DataFrame(rows)
end

number(x::Real) = @sprintf("%.6g", x)
number(::Missing) = "N/A"

function combined_table(rows, spec)
    figure = plot(; xlims=(-0.12, 11.12), ylims=(-0.12, 8.2), size=(1700, 950), dpi=180,
        framestyle=:none, ticks=nothing, grid=false, legend=false,
        title="$(spec.title) ($(spec.unit))", titlefontsize=19,
        background_color=:white, foreground_color=:black,
        margin=4Plots.mm, fontfamily="Computer Modern")
    line(x, y; width=1.2) = plot!(figure, x, y; color=:black, linewidth=width, label="")
    for y in (0, 3, 6)
        line([0, 11], [y, y]; width=2)
    end
    for y in (1, 2, 4, 5)
        line([1, 11], [y, y])
    end
    for y in (7, 8)
        line([1, 11], [y, y]; width=2)
    end
    line([0, 0], [0, 6]; width=2)
    for x in (1, 2, 5, 8, 11)
        line([x, x], [0, 8]; width=2)
    end
    for x in (3, 4, 6, 7, 9, 10)
        line([x, x], [0, 7])
    end
    line([1, 2], [7, 6])
    annotate!(figure, 1.5, 7.5, text("h_t", 16))
    annotate!(figure, 1.65, 6.78, text("Statistic", 10))
    annotate!(figure, 1.32, 6.22, text("i_t (deg)", 10))
    for (h, altitude) in enumerate(EP.ALTITUDES)
        annotate!(figure, 3.5 + 3(h - 1), 7.5, text("$altitude km", 18))
        for (k, label) in enumerate(("Min", "Avg", "Max"))
            annotate!(figure, 2.5 + 3(h - 1) + k - 1, 6.5, text(label, 15))
        end
    end
    for (s, scheduler) in enumerate(SCHEDULERS)
        label = scheduler == "gve_sma" ? "SMA scheduler" : "INC scheduler"
        annotate!(figure, 0.5, 4.5 - 3(s - 1), text(label, 15; rotation=90))
        for (i, inclination) in enumerate(INCLINATIONS)
            y = 5.5 - 3(s - 1) - (i - 1)
            annotate!(figure, 1.5, y, text(@sprintf("%.1f", inclination), 15))
            for (h, altitude) in enumerate(EP.ALTITUDES)
                row = only(eachrow(filter(r -> r.scheduler == scheduler &&
                    r.target_inclination_deg == inclination && r.target_altitude_km == altitude, rows)))
                for (k, field) in enumerate((:min, :avg, :max))
                    annotate!(figure, 2.5 + 3(h - 1) + k - 1, y,
                        text(number(getproperty(row, field)), 13))
                end
            end
        end
    end
    figure
end

function conclusions(io, rows, page)
    populated = filter(r -> r.events > 0, rows)
    isempty(populated) && error("No observations on page $page")
    low, high = minimum(populated.min), maximum(populated.max)
    top = populated[argmax(populated.max), :]
    unit = PAGES[page].unit
    println(io, "## Conclusion\n")
    println(io, "- Across all cells, the observed values span **$(number(low)) to $(number(high)) $unit**. ",
        "One cell attaining the upper bound is **$(top.scheduler), h_t = $(top.target_altitude_km) km, i = $(top.target_inclination_deg) deg**.")
    if page == 1
        groups = groupby(rows, [:scheduler, :target_altitude_km])
        increasing = count(g -> all(diff(sort(g, :target_inclination_deg).max) .> 0), groups)
        println(io, "- The maximum observed entry speed increases from i = 0 to 0.5 to 1.0 deg in **$increasing of $(length(groups)) scheduler/altitude groups**. ",
            "This describes the upper bounds, not an increase in every encounter's speed.")
        println(io, "- These bounds describe **speed at encounter entry**, not closest-approach speed, encounter-average speed, or the maximum speed reached later. ",
            "For encounters already in range at simulation start, entry refers to the initial saved state.")
    elseif page == 2
        positive = sum(rows.positive_entry_events)
        clipped = sum(rows.positive_entry_start_clipped)
        println(io, "- There are **$positive positive entry range-rate observations**, of which **$clipped** are start-clipped. ",
            "A positive value therefore must not automatically be interpreted as a normal inward range crossing; an encounter may already be in range at t = 0.")
    elseif page == 3
        full_window = filter(r -> r.events > 0 && r.full_run_events == r.events, rows)
        for r in eachrow(full_window)
            println(io, "- For **$(r.scheduler), h_t = $(r.target_altitude_km) km, i = $(r.target_inclination_deg) deg**, ",
                "all **$(r.events)** saved encounters span the complete run. Its equal Min/Max values describe the observation window, not measured complete encounter lifetimes.")
        end
        println(io, "- **$(sum(rows.clipped_events)) of $(sum(rows.events)) encounters are boundary-clipped**, including **$(sum(rows.full_run_events)) spanning the full run**. ",
            "Long observed durations can be limited by the simulation window; they are not necessarily complete encounter lifetimes.")
    else
        same_altitude = filter(r -> r.target_altitude_km == 1000 && r.events > 0, rows)
        sma_max = maximum(filter(r -> r.scheduler == "gve_sma", same_altitude).max)
        inc_max = maximum(filter(r -> r.scheduler == "gve_inc", same_altitude).max)
        println(io, "- At **h_t = 1000 km**, the longest observed uninterrupted firing is **$(number(sma_max)) s for SMA**, ",
            "versus **$(number(inc_max)) s for INC**. SMA reaches much longer intervals here, but these maxima do not quantify how frequently they occur.")
        empty = filter(r -> r.events == 0, rows)
        if !isempty(empty)
            regimes = join(["$(r.scheduler), h_t = $(r.target_altitude_km) km, i = $(r.target_inclination_deg) deg"
                for r in eachrow(empty)], "; ")
            println(io, "- **No firing intervals exist in $(nrow(empty)) cells** ($regimes), representing **$(sum(empty.cases)) cases**. ",
                "Their Min/Max entries are N/A, not zero-duration firing intervals.")
        end
        println(io, "- Longer uninterrupted firing is not, by itself, evidence of better orbital control. ",
            "Compare with total duty fraction and the scheduler's intended orbital-change objective.")
    end
    overlap, compared = 0, 0
    for altitude in EP.ALTITUDES, inclination in INCLINATIONS
        pair = filter(r -> r.target_altitude_km == altitude &&
            r.target_inclination_deg == inclination && r.events > 0, rows)
        nrow(pair) == 2 || continue
        compared += 1
        maximum(pair.min) <= minimum(pair.max) && (overlap += 1)
    end
    println(io, "- SMA and INC min-max ranges overlap in **$overlap of $compared** matched altitude/inclination cells with observations from both. ",
        "Min/max alone describe extremes, not typical behavior, statistical uncertainty, or scheduler superiority.\n")
end

function write_report(output)
    rows = CSV.read(joinpath(output, "parameter_pages_min_max.csv"), DataFrame)
    nrow(rows) == 72 || error("Expected 72 page/scheduler/altitude/inclination rows")
    path = joinpath(output, "parameter_pages.md")
    open(path, "w") do io
        for (page, spec) in enumerate(PAGES)
            page > 1 && println(io, "\n---\n")
            println(io, "# Page $page. $(spec.title)\n")
            if page == 1
                println(io, "**Ordering assumption:** the request repeated encounter duration on pages 1 and 3; page 1 is interpreted as encounter speed to cover all four quantities.\n")
            end
            println(io, "**Definition:** $(spec.definition)\n")
            println(io, "**Scope:** target-helper pairs; 100 initial target periods; initial target h_t and i. ",
                "Each Min/Avg/Max group pools all events over helper counts 1, 100, 200, 300, 400, 500, 600, 700 and 800 ",
                "at that scheduler/h_t/i setting (nine cases). Boundary-clipped events are included. ",
                "Avg is the arithmetic event-weighted mean, not an equal-case mean or a within-encounter time average. ",
                "Values use $(spec.unit), rounded to six significant digits; N/A means no observations.\n")
            selected = filter(r -> r.page == page, rows)
            stem = "parameter_page_$(page)_table"
            figure = combined_table(selected, spec)
            for extension in ("png", "pdf")
                savefig(figure, joinpath(output, "$stem.$extension"))
            end
            println(io, "![SMA and INC $(spec.title): Min, Avg, Max by altitude and inclination]($stem.png)\n")
            println(io, "[Open full-resolution table]($stem.png) | [PDF table]($stem.pdf)\n")
            conclusions(io, selected, page)
            println(io, "*Data: [full-precision table values and counts](parameter_pages_min_max.csv); ",
                "[source CSV checksums](parameter_pages_sources.csv). Tables and conclusions are generated from the exported table CSV.*\n")
        end
    end
    path
end

function main(args=ARGS)
    input, output = EP.EV.OUTPUT, OUTPUT
    iseven(length(args)) || error("Expected --input PATH and/or --output PATH")
    for k in 1:2:length(args)
        args[k] in ("--input", "--output") || error("Unknown option: $(args[k])")
        args[k] == "--input" ? (input = abspath(args[k + 1])) : (output = abspath(args[k + 1]))
    end
    data = EP.load_data(input)
    rows = summarize(data)
    mkpath(output)
    CSV.write(joinpath(output, "parameter_pages_min_max.csv"), rows)
    CSV.write(joinpath(output, "parameter_pages_sources.csv"),
        [(source_file=joinpath(input, name), sha256=EP.EV.StudyConfig.filehash(joinpath(input, name)))
            for name in ("case_summary.csv", "target_helper_encounters.csv", "target_helper_laser_intervals.csv")])
    println("Saved ", write_report(output))
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    EncounterParameterPages.main()
end
