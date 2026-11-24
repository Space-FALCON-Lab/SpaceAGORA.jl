#!/usr/bin/env julia

if abspath(PROGRAM_FILE) == @__FILE__
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..",
        "9_SpaceAGORA.jl-main_full integration")))
end

module SelectedAltitudeOverview

include("plot_encounter_distributions.jl")
using .EncounterDistributionPlots
using CSV, DataFrames, Statistics, Printf

const EP = EncounterDistributionPlots
const OUTPUT = EP.EV.OUTPUT
const ALTITUDES = EP.ALTITUDES

function metric_values(data, altitude, metric)
    if metric == :entry_speed
        data.encounters[data.encounters.target_altitude_km .== altitude, :entry_relative_speed_mps]
    elseif metric == :entry_range_rate
        data.encounters[data.encounters.target_altitude_km .== altitude, :entry_range_rate_mps]
    elseif metric == :encounter_duration
        data.encounters[data.encounters.target_altitude_km .== altitude, :duration_s]
    elseif metric == :laser_sma
        data.laser[(data.laser.target_altitude_km .== altitude) .&
            (data.laser.scheduler .== "gve_sma"), :duration_s]
    elseif metric == :laser_inc
        data.laser[(data.laser.target_altitude_km .== altitude) .&
            (data.laser.scheduler .== "gve_inc"), :duration_s]
    else
        throw(ArgumentError("Unknown metric $metric"))
    end
end

function summarize(data)
    rows = NamedTuple[]
    for altitude in ALTITUDES, metric in (:entry_speed, :entry_range_rate,
            :encounter_duration, :laser_sma, :laser_inc)
        values = metric_values(data, altitude, metric)
        all(isfinite, values) || error("Nonfinite values for $metric at $altitude km")
        events = metric in (:laser_sma, :laser_inc) ?
            filter(r -> r.target_altitude_km == altitude &&
                r.scheduler == (metric == :laser_sma ? "gve_sma" : "gve_inc"), data.laser) :
            filter(r -> r.target_altitude_km == altitude, data.encounters)
        length(values) == nrow(events) || error("Event count mismatch for $metric at $altitude km")
        cells = metric in (:laser_sma, :laser_inc) ?
            filter(r -> r.target_altitude_km == altitude &&
                r.scheduler == (metric == :laser_sma ? "gve_sma" : "gve_inc"), data.cases) :
            filter(r -> r.target_altitude_km == altitude, data.cases)
        clipped_cases = metric in (:laser_sma, :laser_inc) ?
            sum(cells.clipped_laser_interval_count) : sum(cells.clipped_encounter_count)
        push!(rows, (target_altitude_km=altitude, metric=string(metric),
            label=metric == :entry_speed ? "Encounter speed at entry" :
                metric == :entry_range_rate ? "Signed range rate at entry" :
                metric == :encounter_duration ? "Encounter duration" :
                metric == :laser_sma ? "Laser-on duration (gve_sma)" : "Laser-on duration (gve_inc)",
            unit=metric in (:entry_speed, :entry_range_rate) ? "m/s" : "s",
            scheduler=metric == :laser_sma ? "gve_sma" :
                metric == :laser_inc ? "gve_inc" : "pooled",
            cases=nrow(cells), events=length(values),
            clipped_events=count(events.start_clipped .| events.end_clipped),
            full_run_events=count(events.start_clipped .& events.end_clipped),
            min=isempty(values) ? missing : minimum(values),
            average=isempty(values) ? missing : mean(values),
            max=isempty(values) ? missing : maximum(values)))
    end
    DataFrame(rows)
end

function summarize_all_altitudes(data)
    rows = NamedTuple[]
    for metric in (:entry_speed, :entry_range_rate, :encounter_duration, :laser_sma, :laser_inc)
        if metric in (:laser_sma, :laser_inc)
            scheduler = metric == :laser_sma ? "gve_sma" : "gve_inc"
            events = filter(r -> r.scheduler == scheduler, data.laser)
            cells = filter(r -> r.scheduler == scheduler, data.cases)
        else
            scheduler = "pooled"
            events = data.encounters
            cells = data.cases
        end
        field = metric == :entry_speed ? :entry_relative_speed_mps :
            metric == :entry_range_rate ? :entry_range_rate_mps : :duration_s
        values = metric in (:laser_sma, :laser_inc) ? events.duration_s : events[!, field]
        all(isfinite, values) || error("Nonfinite values for $metric across selected altitudes")
        push!(rows, (metric=string(metric), label=metric == :entry_speed ?
                "Encounter speed at entry" :
                metric == :entry_range_rate ? "Signed range rate at entry" :
                metric == :encounter_duration ? "Encounter duration" :
                metric == :laser_sma ? "Laser-on duration (gve_sma)" :
                "Laser-on duration (gve_inc)",
            unit=metric in (:entry_speed, :entry_range_rate) ? "m/s" : "s",
            scheduler, cases=nrow(cells), events=nrow(events),
            clipped_events=count(events.start_clipped .| events.end_clipped),
            full_run_events=count(events.start_clipped .& events.end_clipped),
            min=isempty(values) ? missing : minimum(values),
            average=isempty(values) ? missing : mean(values),
            max=isempty(values) ? missing : maximum(values)))
    end
    DataFrame(rows)
end

number(x::Real) = @sprintf("%.6g", x)
number(::Missing) = "N/A"
full_precision(x::Real) = string(x)
full_precision(::Missing) = "N/A"

function presentation_number(metric, value)
    value isa Missing && return "N/A"
    digits = metric == "encounter_duration" ? 6 : 5
    @sprintf("%.*g", digits, value)
end

const DISPLAY = (
    ("entry_speed", "Encounter speed at entry", "m/s"),
    ("entry_range_rate", "Signed range rate at entry", "m/s"),
    ("encounter_duration", "Encounter duration", "s"),
    ("laser_sma", "Laser-on time (gve_sma)", "s"),
    ("laser_inc", "Laser-on time (gve_inc)", "s"),
)

function write_report(output, source_directory)
    rows = CSV.read(joinpath(output, "selected_altitude_statistics.csv"), DataFrame)
    combined = CSV.read(joinpath(output, "selected_altitude_combined_statistics.csv"), DataFrame)
    nrow(rows) == length(ALTITUDES) * length(DISPLAY) ||
        error("Incomplete selected-altitude summary")
    nrow(combined) == length(DISPLAY) || error("Incomplete combined-altitude summary")
    path = joinpath(output, "selected_altitude_summary.md")
    open(path, "w") do io
        println(io, "# All-scenario encounter and laser-on summary\n")
        println(io, "## Scope and definitions\n")
        println(io, "- Covers **all 162 verified 100-orbit cases** at initial target altitudes 1000, 1050 and 1100 km: all three inclinations, all nine helper counts, and both schedulers.")
        println(io, "- For each altitude, speed, signed range rate and encounter duration pool every target-helper encounter across both schedulers. Laser-on duration is the duration of each uninterrupted firing interval, tabulated separately for `gve_sma` and `gve_inc`.")
        println(io, "- **Avg is event-weighted**: each encounter/firing interval is one observation. Scenarios producing more events contribute more; these are descriptive aggregates, not equal-weight averages of scenarios and not independent-replicate estimates.")
        println(io, "- Speed and signed range rate are measured **at encounter entry** from saved start-boundary states. If already in range at simulation start, that initial state is used. Range rate is signed: negative means approaching; positive means receding.")
        println(io, "- Encounter duration is observed time inside 200 km. Firing duration is uninterrupted selected-link firing and may include held firing outside that range. Boundary-clipped events remain included; maxima that equal a whole run's duration are observation-window bounds, not complete encounter lifetimes.")
        println(io, "- Min, Avg and Max are across the observations in each column. `N/A` indicates no observations, not zero. Units are shown in each column header. The first combined table uses full round-trip Float64 values; the duplicate below is presentation-rounded. Per-altitude tables use six significant figures. CSVs retain full precision.\n")
        println(io, "## Full-precision combined summary: all three target altitudes\n")
        println(io, "This table pools every event at **1000, 1050, and 1100 km together**. Entry speed, entry range rate and encounter duration include both schedulers. Laser-on intervals remain separated by scheduler.\n")
        println(io, "| Statistic | ",
            join(["$label ($unit)" for (_, label, unit) in DISPLAY], " | "), " |")
        println(io, "|:---|", repeat(" ---: |", length(DISPLAY)))
        for (statistic, field) in (("Max", :max), ("Avg", :average), ("Min", :min))
            print(io, "| $statistic |")
            for (metric, _, _) in DISPLAY
                row = only(eachrow(filter(r -> r.metric == metric, combined)))
                print(io, " ", full_precision(getproperty(row, field)), " |")
            end
            println(io)
        end
        println(io)
        println(io, "### Presentation-rounded combined table\n")
        println(io, "This repeats the combined statistics above, rounded to the requested precision: **5 significant figures** for encounter speed, signed range rate (assumed), and laser-on durations; **6 significant figures** for encounter duration.\n")
        println(io, "| Statistic | ",
            join(["$label ($unit)" for (_, label, unit) in DISPLAY], " | "), " |")
        println(io, "|:---|", repeat(" ---: |", length(DISPLAY)))
        for (statistic, field) in (("Max", :max), ("Avg", :average), ("Min", :min))
            print(io, "| $statistic |")
            for (metric, _, _) in DISPLAY
                row = only(eachrow(filter(r -> r.metric == metric, combined)))
                print(io, " ", presentation_number(metric, getproperty(row, field)), " |")
            end
            println(io)
        end
        println(io)
        println(io, "**Combined coverage:** ", join(["$(r.events) $(r.label) events" for r in eachrow(combined)], "; "), ".")
        println(io, "**Combined boundary flags:** ", join(["$(r.clipped_events) clipped $(r.label) events ($(r.full_run_events) span the full run)"
            for r in eachrow(combined) if r.clipped_events > 0], "; "), ".\n")
        for altitude in ALTITUDES
            selected = filter(r -> r.target_altitude_km == altitude, rows)
            println(io, "## Target altitude: $altitude km\n")
            println(io, "| Statistic | ",
                join(["$label ($unit)" for (_, label, unit) in DISPLAY], " | "), " |")
            println(io, "|:---|", repeat(" ---: |", length(DISPLAY)))
            for (statistic, field) in (("Max", :max), ("Avg", :average), ("Min", :min))
                print(io, "| $statistic |")
                for (metric, _, _) in DISPLAY
                    row = only(eachrow(filter(r -> r.metric == metric, selected)))
                    print(io, " ", number(getproperty(row, field)), " |")
                end
                println(io)
            end
            println(io)
            println(io, "**Coverage:** ", join(["$(r.events) $(r.label) events" for r in eachrow(selected)], "; "), ".")
            println(io, "**Boundary flags:** ", join(["$(r.clipped_events) clipped $(r.label) events ($(r.full_run_events) span the full run)"
                for r in eachrow(selected) if r.clipped_events > 0], "; "), ".\n")
        end
        println(io, "## Interpretation\n")
        println(io, "This is a compact description of the simulated orbital/encounter regimes, not a scheduler-effect estimate. The encounter columns mix scheduler-conditioned trajectories and weight them by event count. The scheduler split in laser-on time describes the resulting firing intervals, not control success. For claims about scheduler performance, compare matched cases and report per-case orbital-change/control outcomes.\n")
        println(io, "Sources were the verified evaluator CSVs in `$(relpath(source_directory, output))`. ",
            "[Full-precision statistics and counts](selected_altitude_statistics.csv); ",
            "[combined three-altitude statistics](selected_altitude_combined_statistics.csv); ",
            "[source file checksums](selected_altitude_sources.csv).\n")
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
    combined = summarize_all_altitudes(data)
    mkpath(output)
    CSV.write(joinpath(output, "selected_altitude_statistics.csv"), rows)
    CSV.write(joinpath(output, "selected_altitude_combined_statistics.csv"), combined)
    CSV.write(joinpath(output, "selected_altitude_sources.csv"),
        [(source_file=joinpath(input, name), sha256=EP.EV.StudyConfig.filehash(joinpath(input, name)))
            for name in ("case_summary.csv", "target_helper_encounters.csv", "target_helper_laser_intervals.csv")])
    println("Saved ", write_report(output, input))
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    SelectedAltitudeOverview.main()
end
