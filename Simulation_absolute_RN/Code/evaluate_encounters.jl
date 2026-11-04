#!/usr/bin/env julia

if abspath(PROGRAM_FILE) == @__FILE__
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..",
        "9_SpaceAGORA.jl-main_full integration")))
end

module EncounterEvaluation

using Arrow, CSV, DataFrames, Statistics, Printf, TOML
include("StudyConfig.jl")
using .StudyConfig
include("EncounterProducts.jl")
using .EncounterProducts: boundary_state, encounter_quantities

const OUTPUT = joinpath(DEFAULT_OUTPUT, "encounter_evaluation")
const DIMENSIONS = (
    (:target_inclination_deg, "Target inclination i (deg)"),
    (:target_altitude_km, "Initial target altitude h (km)"),
    (:helper_count, "Helper number"),
    (:scheduler, "Scheduler"),
)
const METRICS = (
    (:encounter, :min_relative_speed_mps, "Encounter speed: sampled minimum", "m/s"),
    (:encounter, :max_relative_speed_mps, "Encounter speed: sampled maximum", "m/s"),
    (:encounter, :min_range_rate_mps, "Signed range rate: sampled minimum", "m/s"),
    (:encounter, :max_range_rate_mps, "Signed range rate: sampled maximum", "m/s"),
    (:encounter, :duration_s, "Encounter duration", "s"),
    (:laser, :duration_s, "Uninterrupted laser-on duration", "s"),
    (:encounter, :laser_on_within_encounter_s, "Total laser-on within an encounter", "s"),
)
const INTERVAL_FIELDS = (:sc_a, :sc_b, :start_time_s, :end_time_s, :duration_s,
    :start_clipped, :end_clipped)

function checked_table(directory, name, info)
    path = joinpath(directory, name)
    filehash(path) == info["checksums"][name] || error("Checksum mismatch: $path")
    table = DataFrame(Arrow.Table(path))
    all(==(info["run_id"]), table.run_id) || error("Run ID mismatch: $path")
    table
end

function validate_intervals(table, id, duration, target)
    allunique(table[!, id]) || error("Duplicate $id")
    for row in eachrow(table)
        1 <= row.sc_a < row.sc_b <= target || error("Invalid spacecraft pair")
        all(isfinite, (row.start_time_s, row.end_time_s, row.duration_s)) ||
            error("Nonfinite interval")
        0 <= row.start_time_s < row.end_time_s <= duration || error("Invalid interval bounds")
        isapprox(row.duration_s, row.end_time_s - row.start_time_s; atol=1e-8, rtol=1e-12) ||
            error("Inconsistent interval duration")
        row.start_clipped == (row.start_time_s == 0) || error("Invalid start clipping flag")
        row.end_clipped == (row.end_time_s == duration) || error("Invalid end clipping flag")
    end
    for group in groupby(table, [:sc_a, :sc_b])
        ordered = sort(group, :start_time_s)
        all(ordered.end_time_s[1:end-1] .<= ordered.start_time_s[2:end]) ||
            error("Overlapping intervals for one pair")
    end
end

function attach_firing!(encounters, laser, audit)
    encounters.laser_on_within_encounter_s = zeros(nrow(encounters))
    laser.within_range_s = zeros(nrow(laser))
    laser.outside_range_s = zeros(nrow(laser))
    eindex = Dict(id => k for (k, id) in enumerate(encounters.encounter_id))
    lindex = Dict(id => k for (k, id) in enumerate(laser.laser_id))
    for segment in eachrow(audit)
        haskey(lindex, segment.laser_id) || error("Unknown laser ID in range audit")
        k = lindex[segment.laser_id]
        firing = laser[k, :]
        (segment.sc_a, segment.sc_b) == (firing.sc_a, firing.sc_b) ||
            error("Range audit pair mismatch")
        firing.start_time_s <= segment.start_time_s < segment.end_time_s <= firing.end_time_s ||
            error("Range audit outside firing interval")
        isfinite(segment.duration_s) && isapprox(segment.duration_s,
            segment.end_time_s - segment.start_time_s; atol=1e-8, rtol=1e-12) ||
            error("Invalid range audit duration")
        if segment.within_range
            ismissing(segment.encounter_id) && error("Missing in-range encounter ID")
            haskey(eindex, segment.encounter_id) || error("Unknown in-range encounter ID")
            j = eindex[segment.encounter_id]
            encounter = encounters[j, :]
            (segment.sc_a, segment.sc_b) == (encounter.sc_a, encounter.sc_b) ||
                error("Range audit encounter pair mismatch")
            encounter.start_time_s <= segment.start_time_s < segment.end_time_s <= encounter.end_time_s ||
                error("Range audit outside encounter")
            encounters.laser_on_within_encounter_s[j] += segment.duration_s
            laser.within_range_s[k] += segment.duration_s
        else
            ismissing(segment.encounter_id) || error("Outside-range segment has an encounter ID")
            laser.outside_range_s[k] += segment.duration_s
        end
    end
    for group in groupby(audit, :laser_id)
        ordered = sort(group, :start_time_s)
        firing = laser[lindex[first(group.laser_id)], :]
        first(ordered.start_time_s) == firing.start_time_s &&
            last(ordered.end_time_s) == firing.end_time_s &&
            all(ordered.end_time_s[1:end-1] .== ordered.start_time_s[2:end]) ||
            error("Range audit has gaps or overlaps")
    end
    all(isapprox.(laser.within_range_s + laser.outside_range_s, laser.duration_s;
        atol=1e-7, rtol=1e-10)) || error("Range audit does not cover all firing")
    all(encounters.laser_on_within_encounter_s .<= encounters.duration_s .+ 1e-7) ||
        error("Encounter firing exceeds encounter duration")
    encounters, laser
end

function add_metadata!(table, metadata)
    for (field, value) in pairs(metadata)
        table[!, field] = fill(value, nrow(table))
    end
    table
end

function read_case(root, c)
    directory = case_directory(root, c)
    info_path = joinpath(directory, "run_info.toml")
    verification_path = joinpath(directory, "verification.toml")
    info = TOML.parsefile(info_path)
    verification = TOML.parsefile(verification_path)
    info["status"] == "verified" && verification["passed"] ||
        error("Case is not verified: $directory")
    for (field, value) in case_settings(c)
        info["settings"][field] == value || error("Unexpected $field in $directory")
    end
    info["run_id"] == case_signature(c) || error("Unexpected case signature: $directory")
    geometry = checked_table(directory, "encounters.feather", info)
    laser = checked_table(directory, "laser_firing.feather", info)
    audit = checked_table(directory, "firing_range_audit.feather", info)
    nrow(geometry) == info["encounters"] == verification["encounters"] ||
        error("Encounter count mismatch")
    nrow(laser) == info["firing_intervals"] == verification["firing_intervals"] ||
        error("Firing count mismatch")
    target = c.helpers + 1
    duration = mission_duration(c)
    validate_intervals(geometry, :encounter_id, duration, target)
    validate_intervals(laser, :laser_id, duration, target)
    all(==(target), laser.sc_b) && all(==(target), audit.sc_b) ||
        error("Firing is not target-helper only")
    # One terminal on the target means firing intervals cannot overlap across helpers.
    ordered_laser = sort(laser, :start_time_s)
    all(ordered_laser.end_time_s[1:end-1] .<= ordered_laser.start_time_s[2:end]) ||
        error("Concurrent target firing")
    target_geometry = filter(row -> row.sc_b == target, geometry)
    target_geometry.entry_relative_speed_mps = zeros(nrow(target_geometry))
    target_geometry.entry_range_rate_mps = zeros(nrow(target_geometry))
    for row in eachrow(target_geometry)
        speed, rate, _ = encounter_quantities(
            boundary_state(row, "start", "a"), boundary_state(row, "start", "b"))
        isfinite(speed) && isfinite(rate) || error("Nonfinite encounter entry metric")
        row.entry_relative_speed_mps = speed
        row.entry_range_rate_mps = rate
        row.sample_count >= 2 || error("Encounter lacks boundary samples")
        for quantity in ("relative_speed_mps", "range_rate_mps")
            low = getproperty(row, Symbol("min_$quantity"))
            high = getproperty(row, Symbol("max_$quantity"))
            isfinite(low) && isfinite(high) && low <= high || error("Invalid sampled extrema")
            for kind in ("min", "max")
                time = getproperty(row, Symbol("$(kind)_$(quantity)_time_s"))
                row.start_time_s <= time <= row.end_time_s || error("Extremum time outside encounter")
            end
        end
        row.min_relative_speed_mps >= 0 || error("Negative relative speed")
        max(abs(row.min_range_rate_mps), abs(row.max_range_rate_mps)) <=
            row.max_relative_speed_mps + 1e-8 || error("Range rate exceeds relative speed")
    end
    attach_firing!(target_geometry, laser, audit)
    isapprox(sum(laser.duration_s), info["actual_firing_s"]; atol=1e-6, rtol=1e-10) ||
        error("Total firing mismatch")
    isapprox(sum(laser.outside_range_s), info["outside_range_firing_s"]; atol=1e-6, rtol=1e-10) ||
        error("Outside-range firing mismatch")
    metadata = (run_id=info["run_id"], target_inclination_deg=c.inclination_deg,
        target_altitude_km=c.altitude_km, helper_count=c.helpers, scheduler=string(c.policy))
    encounter_fields = (:encounter_id, INTERVAL_FIELDS..., :sample_count,
        :min_relative_speed_mps, :max_relative_speed_mps, :min_range_rate_mps,
        :max_range_rate_mps, :entry_relative_speed_mps, :entry_range_rate_mps,
        :laser_on_within_encounter_s)
    encounters = add_metadata!(select(target_geometry, collect(encounter_fields)), metadata)
    firing = add_metadata!(select(laser, [:laser_id, :start_encounter_id,
        INTERVAL_FIELDS..., :within_range_s, :outside_range_s]), metadata)
    case_row = (; metadata..., duration_s=duration, encounter_count=nrow(encounters),
        laser_interval_count=nrow(firing),
        clipped_encounter_count=count(encounters.start_clipped .| encounters.end_clipped),
        clipped_laser_interval_count=count(firing.start_clipped .| firing.end_clipped),
        zero_laser_encounter_count=count(iszero, encounters.laser_on_within_encounter_s),
        total_encounter_pair_s=sum(encounters.duration_s),
        total_laser_on_s=sum(firing.duration_s), within_range_laser_on_s=sum(firing.within_range_s),
        outside_range_laser_on_s=sum(firing.outside_range_s),
        laser_duty_fraction=sum(firing.duration_s) / duration,
        source_directory=directory, run_info_sha256=filehash(info_path),
        verification_sha256=filehash(verification_path),
        encounters_sha256=info["checksums"]["encounters.feather"],
        laser_firing_sha256=info["checksums"]["laser_firing.feather"],
        firing_range_audit_sha256=info["checksums"]["firing_range_audit.feather"])
    encounters, firing, case_row
end

function describe(values)
    all(isfinite, values) || error("Nonfinite distribution input")
    n = length(values)
    fields = (:mean, :std_population, :min, :p05, :p25, :median, :p75, :p95, :max)
    n == 0 && return (; n, NamedTuple{fields}(ntuple(_ -> missing, length(fields)))...)
    q = quantile(values, [0, 0.05, 0.25, 0.5, 0.75, 0.95, 1])
    (; n, mean=mean(values), std_population=std(values; corrected=false),
        min=q[1], p05=q[2], p25=q[3], median=q[4], p75=q[5], p95=q[6], max=q[7])
end

function distributions(encounters, laser, cases)
    rows = NamedTuple[]
    groups = [(dimension=:all, value="all", ids=Set(cases.run_id))]
    for (dimension, _) in DIMENSIONS, value in sort(unique(cases[!, dimension]))
        selected = filter(row -> getproperty(row, dimension) == value, cases)
        push!(groups, (dimension=dimension, value=string(value), ids=Set(selected.run_id)))
    end
    for group in groups, population in ("all_observed", "unclipped_only"),
            (kind, field, label, unit) in METRICS
        table = kind == :encounter ? encounters : laser
        selected = filter(row -> row.run_id in group.ids &&
            (population == "all_observed" || !(row.start_clipped || row.end_clipped)), table)
        push!(rows, (; dimension=string(group.dimension), group=group.value, population,
            event_kind=string(kind), metric=string(field), label, unit,
            cases=length(group.ids), contributing_cases=length(unique(selected.run_id)),
            clipped_events=count(selected.start_clipped .| selected.end_clipped),
            zero_events=count(iszero, selected[!, field]), describe(selected[!, field])...))
    end
    DataFrame(rows)
end

fmt(value::Real) = @sprintf("%.6g", value)
fmt(::Missing) = "NA"

function distribution_table(io, rows; grouped=false)
    prefix = grouped ? "Group | " : ""
    println(io, "| $(prefix)Metric (unit) | n | Mean | SD | Min | P05 | P25 | Median | P75 | P95 | Max |")
    println(io, "|", grouped ? "---|" : "", "---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for row in eachrow(rows)
        values = join(fmt.([row.mean, row.std_population, row.min, row.p05, row.p25,
            row.median, row.p75, row.p95, row.max]), " | ")
        println(io, "| ", grouped ? "$(row.group) | " : "", row.label,
            " (", row.unit, ") | ", row.n, " | ", values, " |")
    end
end

function write_report(output)
    # Deliberately reload CSVs: every numeric report value comes from exported data.
    stats = CSV.read(joinpath(output, "distributions.csv"), DataFrame; types=Dict(:group => String))
    cases = CSV.read(joinpath(output, "case_summary.csv"), DataFrame)
    observed = filter(row -> row.population == "all_observed", stats)
    open(joinpath(output, "summary.md"), "w") do io
        println(io, "# Target-helper encounter and laser-on distributions\n")
        println(io, "## Are these quantities collected?\n")
        println(io, "**Yes, all four are collected**, with the following definitions:\n")
        println(io, "| Quantity | Saved product | Definition |\n|---|---|---|")
        println(io, "| Encounter speed | `encounters.feather` | Relative ECI speed `norm(v_b - v_a)`, m/s; sampled minimum and maximum for each encounter. |")
        println(io, "| Range rate | `encounters.feather` | Signed `dot(r_b-r_a, v_b-v_a) / norm(r_b-r_a)`, m/s; negative approaching, positive receding; sampled minimum and maximum. |")
        println(io, "| Encounter duration | `encounters.feather` | End minus start time inside the 200 km geometric range, seconds. |")
        println(io, "| Laser-on duration | `laser_firing.feather` | Duration of each uninterrupted selected-pair firing interval, seconds. Includes any firing outside 200 km while selection is held between updates. |")
        println(io, "| Laser-on within an encounter | `firing_range_audit.feather` | Derived by summing saved within-range segment durations for each encounter; includes zero for encounters never fired upon. |\n")
        println(io, "## Scope and interpretation\n")
        println(io, "- All **$(nrow(cases)) verified cases**, each spanning **100 initial target periods**; target-helper pairs only. Helpers start at 1000 km and 0 deg. Here i and h are the **initial target** inclination and altitude, not their evolving values.")
        println(io, "- Exact grid: i = 0, 0.5, 1 deg; h = 900, 950, 1000, 1050, 1100 km; helpers = 1, 100, 200, 300, 400, 500, 600, 700, 800; schedulers = `gve_sma`, `gve_inc`.")
        println(io, "- Each table pools events over the other three dimensions. **Event-weighted, not equal-case-weighted or time-weighted**: cases producing more events contribute more. These are marginal descriptive comparisons, not isolated causal scheduler effects.")
        println(io, "- Speed/range-rate extrema use **10-second samples plus exact boundary states**. They are sampled estimates, not guaranteed continuous-time extrema. Full speed/range-rate time histories for every encounter are not retained; only the longest-encounter tracks are saved. Do not interpret the extrema distributions as instantaneous-speed distributions.")
        println(io, "- Main tables include boundary-clipped intervals as observed, without extrapolating full lifetimes. Both start- and end-clipping flags are exported; `unclipped_only` distributions are also in the CSV.")
        println(io, "- Each geometric encounter contributes one duration, two speed extrema and two range-rate extrema. Each firing interval contributes one uninterrupted laser duration; it need not match an encounter. Accumulated in-encounter laser-on may combine multiple firing intervals. No-firing cases have no fictional zero-length firing interval; their case total is zero.")
        println(io, "- SD is population SD; percentiles use Julia's default linearly interpolated quantile (Hyndman-Fan type 7). NA means no events, not zero. Numbers shown to six significant digits; CSVs retain full precision.")
        println(io, "- All input settings, verified status, run IDs, and SHA-256 checksums of the three consumed Feather products were checked. This is an evaluation of saved results, not a new propagation or a rerun of trajectory physics verification.\n")
        println(io, "## Campaign totals\n")
        println(io, "- Target-helper encounters: **$(sum(cases.encounter_count))**, of which **$(sum(cases.clipped_encounter_count))** are boundary-clipped and **$(sum(cases.zero_laser_encounter_count))** have no within-encounter laser-on time.")
        println(io, "- Uninterrupted laser-on intervals: **$(sum(cases.laser_interval_count))**, of which **$(sum(cases.clipped_laser_interval_count))** are boundary-clipped.")
        println(io, "- Cases without encounters: **$(count(iszero, cases.encounter_count))**; cases without firing: **$(count(iszero, cases.laser_interval_count))**.")
        println(io, "- Summed firing across independent runs: **$(fmt(sum(cases.total_laser_on_s))) s**, including **$(fmt(sum(cases.outside_range_laser_on_s))) s** outside range. These are run-summed totals, not one mission's duration.")
        println(io, "- Geometric pair-time can overlap among helpers and must not be interpreted as target access time or duty fraction. The case CSV's laser duty fraction is total actual target firing divided by that case's duration.\n")
        distribution_table(io, filter(row -> row.dimension == "all", observed))
        for (number, (dimension, title)) in enumerate(DIMENSIONS)
            println(io, "\n## $number) Distribution by $title\n")
            selected = filter(row -> row.dimension == string(dimension), observed)
            println(io, "### Coverage and case totals\n")
            println(io, "| Group | Cases | Encounters | Clipped encounters | Firing intervals | Clipped firing | Zero-firing cases | Median case total laser-on (s) |")
            println(io, "|---|---:|---:|---:|---:|---:|---:|---:|")
            for value in sort(unique(cases[!, dimension]))
                subset = filter(row -> getproperty(row, dimension) == value, cases)
                println(io, "| $value | $(nrow(subset)) | $(sum(subset.encounter_count)) | ",
                    "$(sum(subset.clipped_encounter_count)) | $(sum(subset.laser_interval_count)) | ",
                    "$(sum(subset.clipped_laser_interval_count)) | $(count(iszero, subset.laser_interval_count)) | ",
                    "$(fmt(median(subset.total_laser_on_s))) |")
            end
            for (heading, indices) in (("Encounter speed", 1:2), ("Signed range rate", 3:4),
                    ("Encounter duration", 5:5), ("Laser-on duration", 6:7))
                println(io, "\n### $heading\n")
                labels = [METRICS[k][3] for k in indices]
                distribution_table(io, filter(row -> row.label in labels, selected); grouped=true)
            end
        end
        println(io, "\n## Reproducibility and exported data\n")
        println(io, "- [Encounter events](target_helper_encounters.csv): one row per target-helper encounter, with clipping flags, sampled extrema, duration, and accumulated within-encounter laser-on.")
        println(io, "- [Laser-on events](target_helper_laser_intervals.csv): one row per uninterrupted firing interval, with within/outside-range durations.")
        println(io, "- [Case summary and provenance](case_summary.csv): all 270 cases, including zero-event cases, totals, source directories, and source checksums.")
        println(io, "- [Distributions](distributions.csv): all numerical distribution tables, plus unclipped-only sensitivity results. `contributing_cases`, `clipped_events` and `zero_events` accompany each distribution.")
        println(io, "- [Evaluation code](../../Code/evaluate_encounters.jl); [tests](../../../Tests/absolute_RN/test_encounter_evaluation.jl).\n")
        println(io, "From the project root:\n\n```bash\njulia --startup-file=no --threads=1 --gcthreads=1 Simulation_absolute_RN/Code/evaluate_encounters.jl\n```")
        println(io, "\nThe pipeline writes event CSVs and case summaries, reloads them to compute distributions, then reloads the distribution/case CSVs to generate this Markdown. No simulation files are modified. Use `--report-only` to regenerate Markdown from existing CSVs; `--input PATH` and `--output PATH` override the campaign and report directories.")
    end
end

function main(args=ARGS)
    input, output, report_only = DEFAULT_OUTPUT, OUTPUT, false
    k = 1
    while k <= length(args)
        option = args[k]
        if option == "--report-only"
            report_only = true
        elseif option in ("--input", "--output")
            k < length(args) || error("Missing value for $option")
            k += 1
            option == "--input" ? (input = abspath(args[k])) : (output = abspath(args[k]))
        else
            error("Unknown option: $option")
        end
        k += 1
    end
    if !report_only
        encounters, firing, cases = DataFrame(), DataFrame(), NamedTuple[]
        for (index, c) in enumerate(grid_cases())
            e, l, row = read_case(input, c)
            append!(encounters, e; cols=:union)
            append!(firing, l; cols=:union)
            push!(cases, row)
            index % 30 == 0 && println("Validated $index / 270 cases")
        end
        length(cases) == 270 && allunique(row.run_id for row in cases) ||
            error("Expected 270 unique cases")
        mkpath(output)
        CSV.write(joinpath(output, "target_helper_encounters.csv"), encounters)
        CSV.write(joinpath(output, "target_helper_laser_intervals.csv"), firing)
        CSV.write(joinpath(output, "case_summary.csv"), cases)
        encounters = CSV.read(joinpath(output, "target_helper_encounters.csv"), DataFrame)
        firing = CSV.read(joinpath(output, "target_helper_laser_intervals.csv"), DataFrame)
        case_table = CSV.read(joinpath(output, "case_summary.csv"), DataFrame)
        CSV.write(joinpath(output, "distributions.csv"), distributions(encounters, firing, case_table))
    end
    write_report(output)
    println("Saved ", joinpath(output, "summary.md"))
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    EncounterEvaluation.main()
end
