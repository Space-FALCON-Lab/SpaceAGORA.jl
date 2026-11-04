#!/usr/bin/env julia

if abspath(PROGRAM_FILE) == @__FILE__
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..",
        "9_SpaceAGORA.jl-main_full integration")))
end

module SMA500ClosestApproach

include("StudyRunner.jl")
using .OracleStudy
using Arrow, CSV, DataFrames, LinearAlgebra, Printf, Statistics, TOML
using .OracleStudy.EncounterProducts: boundary_state, encounter_quantities
using .OracleStudy.EncounterTimingRecorder: bisect_crossing

const SC = OracleStudy.StudyConfig
const OUTPUT = joinpath(SC.DEFAULT_OUTPUT, "encounter_evaluation", "sma500_closest_approach")
const REPORT = joinpath(SC.DEFAULT_OUTPUT, "images", "for_paper",
    "parameter_pages_sma_500_helpers.md")
const ALTITUDES = (1000, 1050, 1100)
const INCLINATIONS = (0.0, 0.5, 1.0)

function quantities(state, a, b)
    r = state.sc[b].pos - state.sc[a].pos
    v = state.sc[b].vel - state.sc[a].vel
    distance = norm(r)
    distance > 0 || error("Coincident spacecraft at closest approach")
    result = (speed=norm(v), rate=dot(r, v) / distance, distance=distance)
    all(isfinite, result) || error("Nonfinite closest-approach state")
    result
end

function consider!(rows, index, time, q)
    if q.distance < rows.closest_separation_m[index] ||
            (q.distance == rows.closest_separation_m[index] && time < rows.closest_time_s[index])
        rows.closest_time_s[index] = time
        rows.closest_separation_m[index] = q.distance
        rows.closest_relative_speed_mps[index] = q.speed
        rows.closest_range_rate_mps[index] = q.rate
    end
end

function scan_interval!(rows, index, left, right, at)
    left <= right || error("Reversed closest-approach interval")
    ql, qr = at(left), at(right)
    consider!(rows, index, left, ql)
    consider!(rows, index, right, qr)
    if ql.rate < 0 < qr.rate
        time = bisect_crossing(t -> at(t).rate, left, right)
        consider!(rows, index, time, at(time))
    end
end

function closest_observer(geometry)
    rows = select(geometry, :run_id, :encounter_id, :sc_a, :sc_b,
        :start_time_s, :end_time_s, :start_clipped, :end_clipped)
    rows.closest_time_s = fill(NaN, nrow(rows))
    rows.closest_separation_m = fill(Inf, nrow(rows))
    rows.closest_relative_speed_mps = fill(NaN, nrow(rows))
    rows.closest_range_rate_mps = fill(NaN, nrow(rows))
    order = sortperm(rows.start_time_s)
    active = Set{Int}()
    cursor = 1
    last_stop = 0.0
    boundary_checks = 0
    max_position_error_m = 0.0
    max_velocity_error_mps = 0.0

    function observe(integrator, start, pairs)
        stop = Float64(integrator.t)
        start == last_stop && 0 < stop - start <= SC.STEP_S + 1e-8 ||
            error("Replay has a gap or invalid step")
        while cursor <= length(order) && rows.start_time_s[order[cursor]] <= stop
            push!(active, order[cursor])
            cursor += 1
        end
        filter!(j -> rows.end_time_s[j] >= start, active)
        state(t) = t == start ? integrator.uprev : t == stop ? integrator.u : integrator(t)
        for j in active
            event = geometry[j, :]
            left, right = max(start, event.start_time_s), min(stop, event.end_time_s)
            left <= right || error("Active encounter does not intersect step")
            for boundary in ("start", "end")
                time = getproperty(event, Symbol("$(boundary)_time_s"))
                if start < time <= stop || (start == 0 && time == 0)
                    saved = state(time)
                    for (endpoint, index) in (("a", event.sc_a), ("b", event.sc_b))
                        original = boundary_state(event, boundary, endpoint)
                        max_position_error_m = max(max_position_error_m,
                            norm(saved.sc[index].pos - original[1:3]))
                        max_velocity_error_mps = max(max_velocity_error_mps,
                            norm(saved.sc[index].vel - original[4:6]))
                    end
                    boundary_checks += 1
                end
            end
            at(t) = quantities(state(t), event.sc_a, event.sc_b)
            scan_interval!(rows, j, left, right, at)
        end
        last_stop = stop
        nothing
    end
    checks() = (; boundary_checks, max_position_error_m, max_velocity_error_mps, last_stop)
    (; observe, rows, checks)
end

case_file(c) = joinpath(OUTPUT, @sprintf("h%d_i%.1f.csv", c.altitude_km, c.inclination_deg))

function extract_case(c)
    source = case_directory(SC.DEFAULT_OUTPUT, c)
    original = TOML.parsefile(joinpath(source, "run_info.toml"))
    original["status"] == "verified" || error("Original case is not verified")
    original["run_id"] == SC.case_signature(c) || error("Original run ID mismatch")
    for name in ("encounters.feather", "trajectory.feather")
        SC.filehash(joinpath(source, name)) == original["checksums"][name] ||
            error("Original checksum mismatch: $name")
    end
    geometry = DataFrame(Arrow.Table(joinpath(source, "encounters.feather")))
    geometry = filter(r -> r.sc_b == c.helpers + 1, geometry)
    isempty(geometry) && error("Missing target-helper encounters")
    observer = closest_observer(geometry)
    replay_root = joinpath(OUTPUT, "replay")
    replay = run_case(c; output_root=replay_root, record_events=false, observer=observer.observe)
    info = TOML.parsefile(joinpath(replay, "run_info.toml"))
    info["schedule_sha256"] == original["schedule_sha256"] ||
        error("Replay scheduler history differs from original")
    info["checksums"]["trajectory.feather"] == original["checksums"]["trajectory.feather"] ||
        error("Replay trajectory is not bit-for-bit identical to original")
    check = observer.checks()
    check.last_stop == SC.mission_duration(c) || error("Incomplete replay")
    check.boundary_checks == 2nrow(geometry) || error("Missing encounter boundary checks")
    check.max_position_error_m < 1e-5 && check.max_velocity_error_mps < 1e-8 ||
        error("Replay boundary states differ from original: $check")
    rows = observer.rows
    all(isfinite, rows.closest_relative_speed_mps) || error("Missing closest-approach speeds")
    all(rows.start_time_s .<= rows.closest_time_s .<= rows.end_time_s) ||
        error("Closest approach outside encounter")
    all(rows.closest_separation_m .<= geometry.min_separation_m .+ 1e-5) ||
        error("Refined closest approach exceeds original sampled minimum")
    interior = (rows.closest_time_s .> rows.start_time_s) .&
        (rows.closest_time_s .< rows.end_time_s)
    all(abs.(rows.closest_range_rate_mps[interior]) .< 1e-5) ||
        error("Interior closest approach has nonzero range rate")
    rows.target_altitude_km = fill(c.altitude_km, nrow(rows))
    rows.target_inclination_deg = fill(c.inclination_deg, nrow(rows))
    rows.helper_count = fill(c.helpers, nrow(rows))
    rows.scheduler = fill(string(c.policy), nrow(rows))
    mkpath(OUTPUT)
    CSV.write(case_file(c), rows)
    SC.atomic_toml(case_file(c) * ".verification.toml", Dict(
        "passed" => true, "run_id" => original["run_id"],
        "source_directory" => source, "replay_directory" => replay,
        "events" => nrow(rows), "boundary_checks" => check.boundary_checks,
        "max_boundary_position_error_m" => check.max_position_error_m,
        "max_boundary_velocity_error_mps" => check.max_velocity_error_mps,
        "trajectory_sha256" => original["checksums"]["trajectory.feather"],
        "schedule_sha256" => original["schedule_sha256"],
        "source_encounters_sha256" => original["checksums"]["encounters.feather"],
        "result_sha256" => SC.filehash(case_file(c))))
    println("VERIFIED closest approaches: h=$(c.altitude_km) i=$(c.inclination_deg), $(nrow(rows)) events")
    flush(stdout)
    nothing
end

function write_report()
    cases = CSV.read(joinpath(SC.DEFAULT_OUTPUT, "encounter_evaluation", "case_summary.csv"), DataFrame)
    summary = DataFrame(target_altitude_km=Int[], target_inclination_deg=Float64[],
        events=Int[], min=Float64[], avg=Float64[], max=Float64[])
    records = DataFrame[]
    for inclination in INCLINATIONS, altitude in ALTITUDES
        c = StudyCase(altitude, inclination, 500, :gve_sma)
        verification = TOML.parsefile(case_file(c) * ".verification.toml")
        verification["passed"] && verification["result_sha256"] == SC.filehash(case_file(c)) ||
            error("Unverified closest-approach result")
        rows = CSV.read(case_file(c), DataFrame)
        expected = only(eachrow(filter(r -> r.run_id == SC.case_signature(c), cases)))
        nrow(rows) == expected.encounter_count || error("Encounter count mismatch")
        all(==(SC.case_signature(c)), rows.run_id) || error("Wrong case in closest-approach CSV")
        allunique(rows.encounter_id) || error("Duplicate closest-approach records")
        values = rows.closest_relative_speed_mps
        all(isfinite, values) || error("Invalid speed")
        push!(summary, (altitude, inclination, length(values), minimum(values),
            mean(values), maximum(values)))
        push!(records, rows)
    end
    CSV.write(joinpath(OUTPUT, "encounters.csv"), vcat(records...))
    CSV.write(joinpath(OUTPUT, "summary.csv"), summary)
    io = IOBuffer()
    println(io, "## 1. Encounter speed at closest approach (m/s)\n")
    println(io, "One relative-speed magnitude per geometric encounter, evaluated when helper-target separation is smallest within the observed encounter (`closest_relative_speed_mps`). This is not minimum relative speed or entry speed.\n")
    println(io, "Closest approach is found from the replayed solver's dense states: compare encounter boundaries and accepted-step endpoints, and refine approaching-to-receding range-rate zero crossings within each accepted step (at most 10 s). The smallest separation over the encounter determines the reported speed. Boundary-clipped encounters use the closest point within the observed interval. Replay trajectories and scheduler histories match the original runs bit-for-bit; encounter boundary states are also checked.\n")
    println(io, "| Inclination (deg) | 1000 km Min | 1000 km Avg | 1000 km Max | 1050 km Min | 1050 km Avg | 1050 km Max | 1100 km Min | 1100 km Avg | 1100 km Max |")
    println(io, "|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for inclination in INCLINATIONS
        cells = String[]
        for altitude in ALTITUDES
            row = only(eachrow(filter(r -> r.target_altitude_km == altitude &&
                r.target_inclination_deg == inclination, summary)))
            append!(cells, [@sprintf("%.6g", row[k]) for k in (:min, :avg, :max)])
        end
        println(io, "| ", @sprintf("%.1f", inclination), " | ", join(cells, " | "), " |")
    end
    println(io, "\n**Closest-approach data:** [Per-encounter times, separations and speeds](../../encounter_evaluation/sma500_closest_approach/encounters.csv); [full-precision summary and counts](../../encounter_evaluation/sma500_closest_approach/summary.csv). Signed range rate in parameter 2 remains measured at entry.\n")
    text = read(REPORT, String)
    pattern = r"(?s)## 1\. Encounter speed.*?(?=## 2\.)"
    length(collect(eachmatch(pattern, text))) == 1 || error("Expected one encounter-speed section")
    text = replace(text, pattern => String(take!(io)))
    text = replace(text, "for parameters 1-3;" => "for parameters 2-3;")
    write(REPORT, text)
    println("UPDATED $REPORT")
end

function main(args)
    if args == ["--report-only"]
        write_report()
    elseif length(args) == 3 && args[1] == "--case"
        c = StudyCase(parse(Float64, args[2]), parse(Float64, args[3]), 500, :gve_sma)
        c.altitude_km in ALTITUDES || error("Unsupported report altitude")
        extract_case(c)
    elseif isempty(args)
        for inclination in INCLINATIONS, altitude in ALTITUDES
            extract_case(StudyCase(altitude, inclination, 500, :gve_sma))
            GC.gc()
        end
        write_report()
    else
        error("Usage: write_sma500_closest_approach.jl [--case ALTITUDE_KM INCLINATION_DEG | --report-only]")
    end
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    SMA500ClosestApproach.main(ARGS)
end
