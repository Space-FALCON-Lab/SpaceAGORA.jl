#!/usr/bin/env julia

if abspath(PROGRAM_FILE) == @__FILE__
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..",
        "9_SpaceAGORA.jl-main_full integration")))
end

module EncounterDistributionPlots

include("evaluate_encounters.jl")
using .EncounterEvaluation
using CSV, DataFrames, Statistics, Printf
get!(ENV, "GKSwstype", "100")
using Plots

const EV = EncounterEvaluation
const ALTITUDES = (1000, 1050, 1100)
const OUTPUT = joinpath(EV.StudyConfig.DEFAULT_OUTPUT, "images",
    "encounter_distributions_1000_1050_1100km")
const COLORS = ("#0072B2", "#D55E00")
const BIN_COUNT = 40
const SPECS = (
    (name="encounter_speed", title="Encounter speed: sampled per-encounter extrema",
        kind=:encounter, fields=(:min_relative_speed_mps, :max_relative_speed_mps),
        labels=("Per-encounter minimum", "Per-encounter maximum"),
        xlabel="Relative speed (m/s)", log=false, zero=false),
    (name="range_rate", title="Signed range rate: sampled per-encounter extrema",
        kind=:encounter, fields=(:min_range_rate_mps, :max_range_rate_mps),
        labels=("Per-encounter minimum", "Per-encounter maximum"),
        xlabel="Signed range rate (m/s)", log=false, zero=false),
    (name="encounter_duration", title="Geometric encounter duration",
        kind=:encounter, fields=(:duration_s,), labels=("All observed encounters",),
        xlabel="Encounter duration (s; log scale)", log=true, zero=false),
    (name="laser_on_duration", title="Uninterrupted laser-on duration",
        kind=:laser, fields=(:duration_s,), labels=("All observed firing intervals",),
        xlabel="Laser-on duration (s; log scale)", log=true, zero=false),
)
const PER_ENCOUNTER_SPECS = (
    (name="encounter_speed", title="Relative speed at encounter entry",
        kind=:encounter, fields=(:entry_relative_speed_mps,), labels=("Every encounter: entry speed",),
        xlabel="Entry relative speed (m/s)", log=false, zero=false),
    (name="range_rate", title="Signed range rate at encounter entry",
        kind=:encounter, fields=(:entry_range_rate_mps,), labels=("Every encounter: entry range rate",),
        xlabel="Entry signed range rate (m/s)", log=false, zero=false),
    SPECS[3],
    (name="laser_on_duration", title="Total laser-on time within each encounter",
        kind=:encounter, fields=(:laser_on_within_encounter_s,), labels=("Every encounter, including zeros",),
        xlabel="Laser-on per encounter (s; log10(1 + t) spacing)", log=true, zero=true),
)

axis_values(spec, values) = spec.zero ? log10.(1 .+ values) : values

function load_data(input)
    cases = CSV.read(joinpath(input, "case_summary.csv"), DataFrame)
    encounters = CSV.read(joinpath(input, "target_helper_encounters.csv"), DataFrame)
    laser = CSV.read(joinpath(input, "target_helper_laser_intervals.csv"), DataFrame)
    cases = filter(r -> r.target_altitude_km in ALTITUDES, cases)
    expected = EV.StudyConfig.grid_cases()
    expected = filter(c -> c.altitude_km in ALTITUDES, expected)
    expected_ids = Set(EV.StudyConfig.case_signature(c) for c in expected)
    nrow(cases) == length(expected_ids) && Set(cases.run_id) == expected_ids ||
        error("Expected the complete 162-case selected-altitude grid")
    for c in expected
        row = only(eachrow(filter(r -> r.run_id == EV.StudyConfig.case_signature(c), cases)))
        (row.target_altitude_km, row.target_inclination_deg, row.helper_count, row.scheduler) ==
            (c.altitude_km, c.inclination_deg, c.helpers, string(c.policy)) ||
            error("Case metadata mismatch")
        row.duration_s == EV.StudyConfig.mission_duration(c) || error("Expected 100-orbit cases")
    end
    metadata = Dict(row.run_id => row for row in eachrow(cases))
    tables = DataFrame[]
    for (table, id, count_field, clipped_field) in (
            (encounters, :encounter_id, :encounter_count, :clipped_encounter_count),
            (laser, :laser_id, :laser_interval_count, :clipped_laser_interval_count))
        selected = filter(r -> r.target_altitude_km in ALTITUDES, table)
        allunique(zip(selected.run_id, selected[!, id])) || error("Duplicate event IDs")
        for row in eachrow(selected)
            haskey(metadata, row.run_id) || error("Event references unknown case")
            case = metadata[row.run_id]
            for field in (:target_altitude_km, :target_inclination_deg, :helper_count, :scheduler)
                getproperty(row, field) == getproperty(case, field) || error("Event metadata mismatch")
            end
            1 <= row.sc_a < row.sc_b == case.helper_count + 1 || error("Not a target-helper event")
        end
        grouped = Dict(first(g.run_id) => g for g in groupby(selected, :run_id))
        for case in eachrow(cases)
            group = get(grouped, case.run_id, selected[1:0, :])
            nrow(group) == getproperty(case, count_field) || error("Event count mismatch")
            count(group.start_clipped .| group.end_clipped) == getproperty(case, clipped_field) ||
                error("Clipped count mismatch")
            EV.validate_intervals(group, id, case.duration_s, case.helper_count + 1)
        end
        push!(tables, selected)
    end
    for spec in SPECS
        table = spec.kind == :encounter ? tables[1] : tables[2]
        for field in spec.fields
            all(isfinite, table[!, field]) || error("Nonfinite $field")
            spec.log && !all(>(0), table[!, field]) && error("Log durations must be positive")
        end
        if length(spec.fields) == 2
            all(table[!, spec.fields[1]] .<= table[!, spec.fields[2]]) ||
                error("Minimum exceeds maximum")
        end
    end
    (encounters=tables[1], laser=tables[2], cases=cases)
end

function bin_edges(values; logarithmic=false, include_zero=false, bins=BIN_COUNT)
    isempty(values) && error("Cannot bin an empty distribution")
    all(isfinite, values) || error("Cannot bin nonfinite values")
    bins > 0 || error("Bin count must be positive")
    low, high = extrema(values)
    if include_zero
        low >= 0 || error("Zero-inclusive duration bins require nonnegative values")
        return 10.0 .^ range(0, max(1, ceil(log10(1 + high))); length=bins + 1) .- 1
    end
    if logarithmic
        low > 0 || error("Logarithmic bins require positive values")
        lower, upper = floor(log10(low)), ceil(log10(high))
        upper == lower && (upper += 1)
        return 10.0 .^ range(lower, upper; length=bins + 1)
    end
    span = max(high - low, 1.0)
    collect(range(low - 0.02span, high + 0.02span; length=bins + 1))
end

function histogram_counts(values, edges)
    length(edges) >= 2 && all(isfinite, edges) && all(diff(edges) .> 0) ||
        error("Invalid histogram edges")
    counts = zeros(Int, length(edges) - 1)
    for value in values
        isfinite(value) && first(edges) <= value <= last(edges) ||
            error("Histogram value outside shared bins")
        index = min(searchsortedlast(edges, value), length(counts))
        counts[index] += 1
    end
    counts
end

function empirical_cdf(values)
    isempty(values) && error("Cannot form an empty ECDF")
    all(isfinite, values) || error("Cannot form a nonfinite ECDF")
    ordered = sort(values)
    last_indices = findall(vcat(diff(ordered) .!= 0, true))
    (x=ordered[last_indices], probability=last_indices ./ length(ordered))
end

function plot_data(data; per_encounter=false)
    histograms, quantiles, cdfs = NamedTuple[], NamedTuple[], NamedTuple[]
    for spec in (per_encounter ? PER_ENCOUNTER_SPECS : SPECS)
        table = spec.kind == :encounter ? data.encounters : data.laser
        values = vcat([table[!, field] for field in spec.fields]...)
        edges = bin_edges(values; logarithmic=spec.log, include_zero=spec.zero)
        schedulers = spec.kind == :laser ? ("gve_inc", "gve_sma") : ("pooled",)
        for altitude in ALTITUDES, scheduler in schedulers
            subset = filter(r -> r.target_altitude_km == altitude &&
                (scheduler == "pooled" || r.scheduler == scheduler), table)
            for (field, label) in zip(spec.fields, spec.labels)
                values = subset[!, field]
                n = length(values)
                n > 0 || error("No events for $(spec.name), $altitude km")
                counts = histogram_counts(values, edges)
                metadata = (quantity=spec.name, target_altitude_km=altitude, metric=string(field),
                    scheduler=scheduler, label=scheduler == "pooled" ? label : scheduler,
                    n=n, zero_events=count(iszero, values),
                    clipped_events=count(subset.start_clipped .| subset.end_clipped))
                for k in eachindex(counts)
                    push!(histograms, (; metadata..., bin_index=k, lower=edges[k],
                        upper=edges[k + 1], count=counts[k], probability=counts[k] / n))
                end
                push!(quantiles, (; metadata..., EV.describe(values)...))
                cdf = empirical_cdf(values)
                push!(cdfs, (; metadata..., population="all_observed", cdf...))
                if spec.log && !per_encounter && spec.kind != :laser
                    complete = subset[.!(subset.start_clipped .| subset.end_clipped), field]
                    isempty(complete) && error("No unclipped events for $(spec.name), $altitude km")
                    push!(cdfs, (; metadata..., population="unclipped_only",
                        empirical_cdf(complete)...))
                end
            end
        end
    end
    (histograms=DataFrame(histograms), quantiles=DataFrame(quantiles), cdfs=cdfs)
end

function make_plot(spec, prepared, cases; per_encounter=false)
    hist = filter(r -> r.quantity == spec.name, prepared.histograms)
    limits = (minimum(hist.lower), maximum(hist.upper))
    max_percent = 100maximum(hist.probability)
    panels = Plots.Plot[]
    for row in 1:2, altitude in ALTITUDES
        altitude_hist = filter(r -> r.target_altitude_km == altitude, hist)
        samples = [first(eachrow(g)) for g in groupby(altitude_hist, :scheduler)]
        n = sum(s.n for s in samples)
        clipped = sum(s.clipped_events for s in samples)
        title = row == 1 ? "Target altitude: $altitude km\nn = $n; clipped = $clipped" : ""
        panel = plot(; title, titlefontsize=11, xlabel=spec.xlabel,
            ylabel=row == 1 ? "Events per bin (%)" : "Cumulative fraction",
            xscale=spec.log && !spec.zero ? :log10 : :identity,
            xlims=Tuple(axis_values(spec, collect(limits))),
            ylims=row == 1 ? (0, max_percent * 1.12) : (0, 1.03),
            legend=:bottomright, legendfontsize=8, guidefontsize=10, tickfontsize=9,
            gridalpha=0.2, framestyle=:box, margin=4Plots.mm)
        if spec.zero
            ticks = vcat(0.0, 10.0 .^ (0:floor(Int, log10(limits[2]))))
            plot!(panel; xticks=(axis_values(spec, ticks), [@sprintf("%g", x) for x in ticks]))
        end
        for (index, selected) in enumerate(groupby(altitude_hist, [:metric, :scheduler]))
            sample = first(eachrow(selected))
            edges = vcat(selected.lower, last(selected.upper))
            color = COLORS[index]
            label = sample.scheduler == "pooled" ? sample.label : "$(sample.label) (n=$(sample.n))"
            if row == 1
                plot!(panel, axis_values(spec, edges), vcat(100selected.probability, 0); seriestype=:steppost,
                    fillrange=0, fillalpha=0.2, color, linewidth=1.5, label)
            else
                curves = filter(r -> r.quantity == spec.name && r.target_altitude_km == altitude &&
                    r.metric == sample.metric && r.scheduler == sample.scheduler, prepared.cdfs)
                for curve in curves
                    complete = curve.population == "unclipped_only"
                    plot!(panel, axis_values(spec, vcat(limits[1], curve.x, limits[2])),
                        vcat(0, curve.probability, 1); seriestype=:steppost,
                        color=complete ? COLORS[2] : color, linestyle=complete ? :dash : :solid,
                        linewidth=1.7, label=complete ? "Unclipped only" : label)
                end
            end
        end
        if spec.name == "range_rate" && !per_encounter
            vline!(panel, [0]; color=:gray, linestyle=:dot, linewidth=0.8, label="")
        end
        push!(panels, panel)
    end
    footer = plot(; framestyle=:none, ticks=nothing, legend=false, grid=false,
        xlims=(0, 1), ylims=(0, 1), margin=0Plots.mm)
    pooling = spec.kind == :laser ? "pooled inclinations and helper counts; schedulers separated" :
        "pooled inclinations, helper counts and schedulers"
    annotate!(footer, 0.5, 0.78, text(
        "162 cases | 100 initial target periods | target-helper pairs | $pooling", 10))
    note = if per_encounter
        spec.zero ? "One value per encounter, including zero. Shared log10(1 + t) bins retain zeros and full tails; percentages are NOT density per second." :
            "One value per encounter; 40 shared bins, normalized as event percentages. Every observed encounter included, including boundary-clipped events."
    elseif spec.kind == :laser
        "Schedulers normalized separately using shared log-width bins. Blue: gve_inc; orange: gve_sma. Percentages are NOT density per second."
    else
        spec.log ?
        "40 equal log-width bins; heights are event percentages, NOT density per second. Full tails retained; clipped events included." :
        "40 shared linear bins; each minimum/maximum distribution normalized separately. Extrema sampled at 10 s plus encounter boundaries."
    end
    annotate!(footer, 0.5, 0.48, text(note, 10))
    last_note = if per_encounter
        spec.zero ? "Laser-on is summed WITHIN each geometric encounter, not individual firing intervals. Encounters with no firing contribute zero." :
            spec.name == "encounter_duration" ? "Duration is the observed time within 200 km; clipped intervals are not extrapolated to complete lifetimes." :
            "Entry uses saved start boundary states; encounters already in range at t=0 use that initial state. Not a time average or an extremum."
    elseif spec.name == "laser_on_duration"
        zero = count(iszero, cases.laser_interval_count)
        "$zero / $(nrow(cases)) cases never fire: all gve_inc, 9 at each altitude. No artificial zero-duration intervals. Held out-of-range firing included."
    elseif spec.name == "range_rate"
        "Negative range rate: approaching; positive: receding. These are distributions of encounter extrema, NOT instantaneous range rates."
    elseif spec.name == "encounter_duration"
        "Dashed ECDF excludes boundary-clipped events as a sensitivity comparison; it is NOT a censoring-corrected lifetime estimate."
    else
        "Relative speed = norm(v_target - v_helper). These are distributions of encounter extrema, NOT instantaneous speeds."
    end
    annotate!(footer, 0.5, 0.18, text(last_note, 10))
    layout = @layout [grid(2, 3); footer{0.13h}]
    plot(panels..., footer; layout, size=(1680, 1040), dpi=180,
        plot_title=spec.title, plot_titlefontsize=16, fontfamily="Computer Modern",
        left_margin=7Plots.mm, right_margin=4Plots.mm, bottom_margin=6Plots.mm)
end

function main(args=ARGS)
    input, output, per_encounter = EV.OUTPUT, nothing, false
    k = 1
    while k <= length(args)
        if args[k] == "--per-encounter"
            per_encounter = true
        else
            args[k] in ("--input", "--output") || error("Unknown option: $(args[k])")
            k < length(args) || error("Missing path for $(args[k])")
            args[k] == "--input" ? (input = abspath(args[k + 1])) : (output = abspath(args[k + 1]))
            k += 1
        end
        k += 1
    end
    output === nothing && (output = per_encounter ? joinpath(OUTPUT, "per_encounter") : OUTPUT)
    data = load_data(input)
    prepared = plot_data(data; per_encounter)
    mkpath(output)
    CSV.write(joinpath(output, "histogram_bins.csv"), prepared.histograms)
    CSV.write(joinpath(output, "distribution_quantiles.csv"), prepared.quantiles)
    CSV.write(joinpath(output, "selected_cases.csv"), data.cases)
    provenance = [(source_file=joinpath(input, name),
        sha256=EV.StudyConfig.filehash(joinpath(input, name)))
        for name in ("case_summary.csv", "target_helper_encounters.csv", "target_helper_laser_intervals.csv")]
    CSV.write(joinpath(output, "input_checksums.csv"), provenance)
    gr()
    for spec in (per_encounter ? PER_ENCOUNTER_SPECS : SPECS)
        figure = make_plot(spec, prepared, data.cases; per_encounter)
        for extension in ("png", "pdf")
            path = joinpath(output, "$(spec.name)_histogram_ecdf.$extension")
            savefig(figure, path)
            println("Saved $path")
        end
    end
    println("Validated $(nrow(data.cases)) cases, $(nrow(data.encounters)) encounters, ",
        "$(nrow(data.laser)) firing intervals. Histograms include every selected event.")
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    EncounterDistributionPlots.main()
end
