using Statistics
using TOML

length(ARGS) in (2, 4) ||
    error("Usage: compare_sma_runs.jl helper_count orbit_count [gve_inc target_inclination_deg]")
const N_HELPERS = parse(Int, ARGS[1])
const ORBITS = parse(Int, ARGS[2])
N_HELPERS > 0 && ORBITS > 0 || error("Helper and orbit counts must be positive.")
const INC_CASE = length(ARGS) == 4
const TARGET_INC_DEG = INC_CASE ? parse(Float64, ARGS[4]) : 0.0
!INC_CASE || ARGS[3] == "gve_inc" ||
    error("The optional comparison policy must be gve_inc.")
isfinite(TARGET_INC_DEG) && 0 <= TARGET_INC_DEG <= 180 ||
    error("Target inclination must be finite and between 0 and 180 degrees.")
const RESULT_DIR = joinpath(@__DIR__, "results", "$(N_HELPERS)_helpers_$(ORBITS)_orbits" *
    (INC_CASE ? "_gve_inc_target_inc_$(TARGET_INC_DEG)_deg" : ""))
const CODES = ("prototype", "4_v1", "9_full")

function load_run(code)
    summary = TOML.parsefile(joinpath(RESULT_DIR, code, "summary.toml"))
    @assert summary["retcode"] == "Success"
    @assert summary["helper_count"] == N_HELPERS
    @assert summary["requested_target_orbits"] == ORBITS
    if INC_CASE
        @assert summary["schedule"] == "gve_inc"
        @assert isapprox(summary["initial_inclination_deg"], TARGET_INC_DEG; atol=1e-10)
        @assert summary["initial_target_altitude_m"] == 1050e3
        @assert summary["initial_helper_altitude_m"] == 1000e3
        @assert summary["initial_helper_inclination_deg"] == 0.0
        @assert isfinite(summary["final_inclination_deg"])
        @assert isapprox(summary["final_inclination_deg"] - summary["initial_inclination_deg"],
            summary["delta_inclination_deg"]; atol=1e-14)
    end
    @assert isapprox(summary["duration_s"], ORBITS * summary["target_period_s"])
    @assert isapprox(summary["final_sma_m"] - summary["initial_sma_m"],
        summary["delta_sma_m"]; atol=1e-8)
    @assert all(isfinite(summary[key]) for key in
        ("initial_sma_m", "final_sma_m", "delta_sma_m", "final_eccentricity", "runtime_s"))
    lines = readlines(joinpath(RESULT_DIR, code, "schedule.csv"))
    @assert first(lines) == "time_s,target_orbits,active_helper"
    rows = map(lines[2:end]) do line
        fields = split(line, ',')
        @assert length(fields) == 3
        time, orbits = parse.(Float64, fields[1:2])
        helper = parse(Int, fields[3])
        @assert isfinite(time) && isfinite(orbits)
        @assert isapprox(orbits, time / summary["target_period_s"]; atol=1e-12)
        @assert helper == 0 || (code == "prototype" ?
            1 <= helper <= N_HELPERS : 2 <= helper <= N_HELPERS + 1)
        slot = helper == 0 ? 0 : code == "prototype" ? helper : helper - 1
        return (time=time, slot=slot)
    end
    times = getproperty.(rows, :time)
    slots = getproperty.(rows, :slot)
    @assert first(times) == 0.0
    @assert last(times) == summary["duration_s"]
    @assert all(>(0), diff(times))
    @assert length(rows) == summary[code == "prototype" ? "selection_samples" : "scheduler_updates"]
    @assert count(>(0), slots) == summary["active_updates"]
    expected_final = summary["final_helper"]
    @assert last(slots) == (expected_final == 0 ? 0 :
        code == "prototype" ? expected_final : expected_final - 1)

    events = NamedTuple{(:slot, :start_s, :end_s, :censored), Tuple{Int, Float64, Float64, Bool}}[]
    start = 1
    for i in 2:length(rows)
        if slots[i] != slots[i - 1]
            slots[i - 1] > 0 && push!(events,
                (slot=slots[i - 1], start_s=times[start], end_s=times[i], censored=false))
            start = i
        end
    end
    last(slots) > 0 && push!(events,
        (slot=last(slots), start_s=times[start], end_s=last(times), censored=true))
    # Label only a selection at initialization that disappears within one second.
    startup(e) = e.start_s == 0.0 && e.end_s < 1.0
    contacts = filter(e -> !startup(e), events)
    if code == "4_v1"
        @assert summary["link_activations"] == length(events)
        @assert summary["physical_link_activations"] == length(contacts)
        @assert summary["startup_link_activations"] == count(startup, events)
    end
    active_s = sum(e.end_s - e.start_s for e in events)
    metrics = Dict(
        "physical_helper_sequence" => getproperty.(contacts, :slot),
        "recorded_active_segments" => length(events),
        "nonstartup_segments" => length(contacts),
        "startup_segments" => count(startup, events),
        "estimated_active_time_s" => active_s,
        "estimated_duty_fraction" => active_s / last(times),
        "median_sample_gap_s" => median(diff(times)),
        "maximum_sample_gap_s" => maximum(diff(times)),
        "final_physical_helper_slot" => last(slots),
    )
    return (; summary, times, slots, events, contacts, metrics)
end

function compare_pair(left, right)
    @assert left.summary["duration_s"] == right.summary["duration_s"]
    @assert left.summary["initial_sma_m"] == right.summary["initial_sma_m"]
    boundaries = sort!(unique(vcat(left.times, right.times)))
    disagreement_s = 0.0
    for i in 1:(length(boundaries) - 1)
        t = boundaries[i]
        if left.slots[searchsortedlast(left.times, t)] != right.slots[searchsortedlast(right.times, t)]
            disagreement_s += boundaries[i + 1] - t
        end
    end
    matching_sequence = getproperty.(left.contacts, :slot) == getproperty.(right.contacts, :slot)
    metrics = Dict{String, Any}(
        "same_nonstartup_helper_sequence" => matching_sequence,
        "absolute_delta_sma_difference_m" =>
            abs(left.summary["delta_sma_m"] - right.summary["delta_sma_m"]),
        "held_schedule_disagreement_s" => disagreement_s,
        "held_schedule_agreement_fraction" => 1 - disagreement_s / last(boundaries),
    )
    if matching_sequence && !isempty(left.contacts)
        metrics["max_absolute_start_time_difference_s"] = maximum(
            abs(a.start_s - b.start_s) for (a, b) in zip(left.contacts, right.contacts))
        metrics["max_absolute_end_time_difference_s"] = maximum(
            abs(a.end_s - b.end_s) for (a, b) in zip(left.contacts, right.contacts))
    end
    return metrics
end

runs = Dict(code => load_run(code) for code in CODES)
comparison = Dict{String, Any}(code => runs[code].metrics for code in CODES)
gains = [runs[code].summary["delta_sma_m"] for code in CODES]
comparison["sma_gain_spread_m"] = maximum(gains) - minimum(gains)
if !INC_CASE
    comparison["sma_gain_spread_percent_of_mean"] =
        100 * (maximum(gains) - minimum(gains)) / mean(gains)
end
if INC_CASE
    inc_gains = [runs[code].summary["delta_inclination_deg"] for code in CODES]
    comparison["inclination_gain_spread_deg"] = maximum(inc_gains) - minimum(inc_gains)
    comparison["inclination_gain_spread_percent_of_mean"] =
        100 * (maximum(inc_gains) - minimum(inc_gains)) / mean(inc_gains)
end
for i in 1:(length(CODES) - 1), j in (i + 1):length(CODES)
    left, right = CODES[i], CODES[j]
    comparison[left * "_vs_" * right] = compare_pair(runs[left], runs[right])
end

open(joinpath(RESULT_DIR, "contacts.csv"), "w") do io
    println(io, "code,event,physical_helper_slot,initial_phase_deg,start_s,end_s,duration_s,censored_at_end")
    for code in CODES
        for (i, e) in enumerate(runs[code].events)
            println(io, join((code, i, e.slot, 360 * (e.slot - 1) / N_HELPERS,
                e.start_s, e.end_s, e.end_s - e.start_s, e.censored), ','))
        end
    end
end
open(joinpath(RESULT_DIR, "comparison.toml"), "w") do io
    TOML.print(io, comparison; sorted=true)
end
println("Validated all $(length(CODES)) summaries and schedules.")
TOML.print(stdout, comparison; sorted=true)
