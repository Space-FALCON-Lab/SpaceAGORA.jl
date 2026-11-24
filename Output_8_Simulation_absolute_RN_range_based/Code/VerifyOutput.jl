module StudyVerification

include("StudyConfig.jl")
using .StudyConfig
import Pkg
Pkg.activate(LIBRARY_ROOT; io=devnull)
using Arrow, DataFrames, Dates, LinearAlgebra, StaticArrays, TOML
include("EncounterProducts.jl")
using .EncounterProducts

export verify_output

const SIGNED_T_MONOTONICITY_TOLERANCE_MPS = 1e-7
const SIGNED_T_SUM_TOLERANCE_MPS = 1e-8
const BRANCH_TOLERANCE_MPS = Dict("R" => 1e-8, "T" => SIGNED_T_MONOTONICITY_TOLERANCE_MPS,
    "N" => 1e-8)

require(condition, message) = condition || error(message)
vector_at(table, sc, field, row) = SVector{3,Float64}(ntuple(
    k -> table[Symbol("sc$(sc)_$(field)_$k")][row], 3))

function check_signed_totals(positive, negative, net, last_positive, last_negative, component="T")
    # RK stages and dense interpolation can slightly reverse a branch at a zero crossing.
    tolerance = BRANCH_TOLERANCE_MPS[component]
    require(positive >= -tolerance && negative <= tolerance,
        "$component branches have the wrong signs")
    residual = max(0.0, last_positive - positive, negative - last_negative)
    require(residual <= tolerance, "$component branches have the wrong monotonicity")
    require(abs(net - positive - negative) < SIGNED_T_SUM_TOLERANCE_MPS,
        "$component positive and negative branches do not sum to the signed net")
    residual
end

function check_intervals(table, id, duration, count)
    require(length(unique(table[!, id])) == nrow(table), "Duplicate $id")
    for event in eachrow(table)
        require(1 <= event.sc_a < event.sc_b <= count, "Invalid spacecraft pair")
        require(0 <= event.start_time_s < event.end_time_s <= duration, "Invalid event times")
        require(event.duration_s == event.end_time_s - event.start_time_s, "Invalid event duration")
        require(!event.start_clipped || event.start_time_s == 0, "Invalid start clipping")
        require(!event.end_clipped || event.end_time_s == duration, "Invalid end clipping")
    end
    for group in groupby(table, [:sc_a, :sc_b])
        ordered = sort(group, :start_time_s)
        require(all(ordered.end_time_s[1:end-1] .<= ordered.start_time_s[2:end]),
            "Overlapping same-pair intervals")
    end
end

function check_audit(geometry, laser, segments, run_id)
    require(all(segments.run_id .== run_id), "Range audit run-ID mismatch")
    require(all(segments.laser_id .∈ Ref(Set(laser.laser_id))), "Unknown firing ID in audit")
    require(all(segments.within_range .== .!ismissing.(segments.encounter_id)), "Invalid range flags")
    require(all(segments.end_time_s .> segments.start_time_s), "Nonpositive audit segment")
    require(all(segments.duration_s .== segments.end_time_s .- segments.start_time_s), "Invalid segment duration")
    by_id = Dict(row.encounter_id => row for row in eachrow(geometry))
    for firing in eachrow(laser)
        parts = sort(filter(s -> s.laser_id == firing.laser_id, segments), :start_time_s)
        require(!isempty(parts), "Firing interval has no audit segments")
        require(first(parts.start_time_s) == firing.start_time_s &&
            last(parts.end_time_s) == firing.end_time_s, "Audit coverage does not match firing")
        require(parts.end_time_s[1:end-1] == parts.start_time_s[2:end], "Audit gap or overlap")
        require(isapprox(sum(parts.duration_s), firing.duration_s; atol=1e-8, rtol=1e-12),
            "Inside/outside time does not sum to actual firing")
        for part in eachrow(parts)
            require((part.sc_a, part.sc_b) == (firing.sc_a, firing.sc_b), "Audit pair mismatch")
            if part.within_range
                require(haskey(by_id, part.encounter_id), "Unknown encounter ID")
                event = by_id[part.encounter_id]
                require((part.sc_a, part.sc_b) == (event.sc_a, event.sc_b) &&
                    event.start_time_s <= part.start_time_s < part.end_time_s <= event.end_time_s,
                    "Within-range segment is not inside its encounter")
            else
                require(!any(event -> event.sc_a == part.sc_a && event.sc_b == part.sc_b &&
                    event.start_time_s < part.end_time_s && event.end_time_s > part.start_time_s,
                    eachrow(geometry)), "Outside-range segment overlaps an encounter")
            end
        end
    end
end

function check_encounter_products(geometry, tracks, target, duration)
    selected = longest_ids(geometry, target)
    require(Set(tracks.encounter_id) == Set(filter(!iszero, collect(selected))), "Missing longest encounter track")
    for event in eachrow(geometry)
        expected_count = 2 + max(0, ceil(Int, event.end_time_s / ENCOUNTER_SAMPLE_S) -
            floor(Int, event.start_time_s / ENCOUNTER_SAMPLE_S) - 1)
        require(event.sample_count == expected_count, "Encounter extrema sample count mismatch")
        require(0 <= event.min_relative_speed_mps <= event.max_relative_speed_mps,
            "Invalid encounter speed extrema")
        require(event.min_range_rate_mps <= event.max_range_rate_mps, "Invalid range-rate extrema")
        require(0 < event.min_separation_m <= event.max_separation_m <= RANGE_M + 1e-3,
            "Invalid encounter range extrema")
        for quantity in ("relative_speed_mps", "range_rate_mps", "separation_m"), kind in ("min", "max")
            value = getproperty(event, Symbol("$(kind)_$quantity"))
            time = getproperty(event, Symbol("$(kind)_$(quantity)_time_s"))
            require(isfinite(value) && event.start_time_s <= time <= event.end_time_s,
                "Invalid sampled extremum")
            require(time in (event.start_time_s, event.end_time_s) ||
                rem(time, ENCOUNTER_SAMPLE_S) == 0, "Extremum time is not a boundary or dense sample")
        end
    end
    for group in groupby(tracks, :encounter_id)
        id = first(group.encounter_id)
        event = only(eachrow(filter(e -> e.encounter_id == id, geometry)))
        expected_times = vcat(event.start_time_s,
            collect((floor(Int, event.start_time_s / ENCOUNTER_SAMPLE_S) + 1) * ENCOUNTER_SAMPLE_S:
                ENCOUNTER_SAMPLE_S:(ceil(Int, event.end_time_s / ENCOUNTER_SAMPLE_S) - 1) * ENCOUNTER_SAMPLE_S),
            event.end_time_s)
        require(group.time == expected_times, "Longest track cadence or endpoints mismatch")
        require(all(group.sc_a .== event.sc_a) && all(group.sc_b .== event.sc_b), "Longest track pair mismatch")
        require(all(group.is_all_pair_longest .== (id == selected[1])) &&
            all(group.is_target_helper_longest .== (id == selected[2])), "Longest track scope mismatch")
        values = zeros(3, nrow(group))
        for (i, row) in enumerate(eachrow(group))
            a = SVector{6,Float64}(ntuple(k -> getproperty(row,
                Symbol("sc_a_$(k <= 3 ? "pos" : "vel")_$(mod1(k, 3))")), 6))
            b = SVector{6,Float64}(ntuple(k -> getproperty(row,
                Symbol("sc_b_$(k <= 3 ? "pos" : "vel")_$(mod1(k, 3))")), 6))
            values[:, i] .= encounter_quantities(a, b)
            if i in (1, nrow(group))
                boundary = i == 1 ? "start" : "end"
                require(a == boundary_state(event, boundary, "a") &&
                    b == boundary_state(event, boundary, "b"), "Longest track boundary state mismatch")
            end
        end
        for (i, quantity) in enumerate(("relative_speed_mps", "range_rate_mps", "separation_m"))
            require(minimum(values[i, :]) == getproperty(event, Symbol("min_$quantity")) &&
                maximum(values[i, :]) == getproperty(event, Symbol("max_$quantity")),
                "Saved longest track does not reproduce its sampled extrema")
        end
    end
end

function verify_output(directory; require_events=true)
    info_path = joinpath(directory, "run_info.toml")
    info = TOML.parsefile(info_path)
    require(info["status"] in ("complete", "verified"), "Run is not complete")
    require(!require_events || info["record_events"], "Production event recording was disabled")
    settings = info["settings"]
    require(settings["schema_version"] == SCHEMA_VERSION, "Unexpected output schema version")
    case = StudyCase(settings["target_altitude_km"], settings["target_inclination_deg"],
        settings["helper_count"], Symbol(settings["scheduler"]);
        orbits=settings["initial_target_periods"], power_W=settings["laser_power_W"],
        helper_altitude_km=settings["helper_altitude_km"])
    require(settings == case_settings(case) && info["run_id"] == case_signature(case),
        "Run settings or diagnostic contract do not match this campaign")
    count, duration = settings["helper_count"] + 1, settings["duration_s"]
    no_authority = settings["laser_power_W"] == 0 ||
        (settings["scheduler"] == "gve_inc" && settings["target_inclination_deg"] == 0)
    checksums = info["checksums"]
    for name in ("trajectory.feather", "encounters.feather", "laser_firing.feather",
        "firing_range_audit.feather", "longest_encounter.feather")
        require(haskey(checksums, name), "Missing output checksum: $name")
    end
    for (name, checksum) in checksums
        require(isfile(joinpath(directory, name)), "Missing output: $name")
        require(filehash(joinpath(directory, name)) == checksum, "Checksum mismatch: $name")
    end
    function read_product(name)
        table = Arrow.Table(joinpath(directory, name))
        metadata = Arrow.getmetadata(table)
        require(metadata !== nothing, "Missing product metadata: $name")
        require(get(metadata, "schema_version", "") == string(SCHEMA_VERSION) &&
            get(metadata, "run_id", "") == info["run_id"], "Product metadata mismatch: $name")
        DataFrame(table)
    end
    geometry = read_product("encounters.feather")
    laser = read_product("laser_firing.feather")
    segments = read_product("firing_range_audit.feather")
    tracks = read_product("longest_encounter.feather")
    check_intervals(geometry, :encounter_id, duration, count)
    check_intervals(laser, :laser_id, duration, count)
    require(all(geometry.run_id .== info["run_id"]) && all(laser.run_id .== info["run_id"]),
        "Event run-ID mismatch")
    require(all(laser.sc_b .== count), "Firing includes helper-helper links")
    require(all(laser.power_W .== settings["laser_power_W"]) &&
        all(laser.force_N .== settings["force_N"]), "Firing power/force mismatch")
    require(isempty(laser) || settings["force_N"] > 0, "Zero-force selection recorded as firing")
    check_audit(geometry, laser, segments, info["run_id"])
    check_encounter_products(geometry, tracks, count, duration)
    require(all(tracks.run_id .== info["run_id"]), "Longest track run-ID mismatch")
    require(!isfile(joinpath(directory, "encounter_states.bin.partial")), "Temporary encounter spool was not removed")
    require(nrow(geometry) == info["encounters"] && nrow(laser) == info["firing_intervals"],
        "Event count mismatch")
    require(sum(laser.duration_s) == info["actual_firing_s"], "Firing total mismatch")
    require(sum(segments.duration_s[.!segments.within_range]) == info["outside_range_firing_s"],
        "Outside-range total mismatch")

    max_boundary_residual = 0.0
    for event in eachrow(geometry), boundary in ("start", "end")
        getproperty(event, Symbol("$(boundary)_clipped")) && continue
        a = SVector{3,Float64}(ntuple(k -> getproperty(event,
            Symbol("$(boundary)_sc_a_r_$("xyz"[k])")), 3))
        b = SVector{3,Float64}(ntuple(k -> getproperty(event,
            Symbol("$(boundary)_sc_b_r_$("xyz"[k])")), 3))
        max_boundary_residual = max(max_boundary_residual, abs(norm(b - a) - RANGE_M))
    end
    require(max_boundary_residual < 1e-3, "Encounter boundary state is not on the range sphere")

    rows = 0
    expected_times = collect(range(0.0, duration; length=TRAJECTORY_ROWS))
    previous_time = -Inf
    initial_energy, initial_Hz, initial_Hnorm = 0.0, 0.0, 0.0
    initial_sma, initial_inc = 0.0, 0.0
    final_sma, final_inc = 0.0, 0.0
    max_energy_residual, max_axial_drift, max_impulse_balance = 0.0, 0.0, 0.0
    max_monotonicity = zeros(3)
    max_sum_residual = zeros(3)
    max_sign_residual = zeros(3)
    last_positive, last_negative = zeros(3), zeros(3)
    final_rtn = zeros(3)
    final_work = 0.0
    for table in Arrow.Stream(joinpath(directory, "trajectory.feather"))
        metadata = Arrow.getmetadata(table)
        require(metadata !== nothing, "Missing trajectory metadata")
        require(get(metadata, "run_id", "") == info["run_id"], "Trajectory run-ID mismatch")
        require(get(metadata, "schema_version", "") == string(SCHEMA_VERSION) &&
            get(metadata, "target_id", "") == string(count), "Trajectory schema/target mismatch")
        require(get(metadata, "units", "") == "time:s; pos:m; vel,laser_dv,RTN:m/s; work:J" &&
            get(metadata, "frame", "") == "J2000 ECI", "Trajectory units/frame mismatch")
        require(get(metadata, "rtn_branch_definition", "") == RTN_BRANCH_DEFINITION,
            "Trajectory R/T/N branch definition mismatch")
        require(collect(propertynames(table)) == trajectory_columns(count), "Unexpected trajectory schema")
        positive_columns = [table[Symbol("target_dv_$(component)_positive_mps")] for component in ("R", "T", "N")]
        negative_columns = [table[Symbol("target_dv_$(component)_negative_mps")] for component in ("R", "T", "N")]
        for name in propertynames(table)
            require(all(isfinite, table[name]), "Nonfinite trajectory column $name")
            require(!occursin("delta_sma", string(name)), "Unwanted direct-SMA column")
        end
        for row in eachindex(table.time)
            t = table.time[row]
            rows += 1
            require(t > previous_time && 0 <= t <= duration, "Invalid trajectory time")
            require(rows <= TRAJECTORY_ROWS, "Too many trajectory rows")
            expected = expected_times[rows]
            require(t == expected, "Trajectory time is not on the requested output grid")
            energy, work, Hz = 0.0, 0.0, 0.0
            H = zero(SVector{3,Float64})
            impulse = zero(SVector{3,Float64})
            for sc in 1:count
                r, v = vector_at(table, sc, "pos", row), vector_at(table, sc, "vel", row)
                radius = norm(r)
                require(radius > EARTH_RADIUS_M + RANGE_M, "Outside verified LOS domain")
                potential = -MU / radius + MU * J2 * EARTH_RADIUS_M^2 / (2radius^3) *
                    (3(r[3] / radius)^2 - 1)
                energy += MASS_KG * (dot(v, v) / 2 + potential)
                H += MASS_KG * cross(r, v)
                impulse += MASS_KG * vector_at(table, sc, "laser_dv", row)
                work += table[Symbol("sc$(sc)_laser_work_J")][row]
                if rows == 1
                    saved = info["initial_states"][sc]
                    require(saved["id"] == sc && saved["index"] == sc, "Initial ID mismatch")
                    require(r == saved["position_m"] && v == saved["velocity_mps"], "Initial state mismatch")
                end
                if rows == 1 || no_authority
                    require(iszero(vector_at(table, sc, "laser_dv", row)) &&
                        table[Symbol("sc$(sc)_laser_work_J")][row] == 0,
                        "Nonzero initial or no-authority spacecraft diagnostic")
                end
                if sc == count
                    final_sma = inv(2 / radius - dot(v, v) / MU)
                    h = cross(r, v)
                    final_inc = rad2deg(atan(hypot(h[1], h[2]), h[3]))
                end
            end
            Hz = H[3]
            if rows == 1
                initial_energy, initial_Hz, initial_Hnorm = energy, Hz, norm(H)
                initial_sma, initial_inc = final_sma, final_inc
            end
            residual = abs(energy - initial_energy - work)
            require(residual <= max(20.0, Float64(count)) + 1e-4abs(work),
                "J2 energy change does not match optical work")
            max_energy_residual = max(max_energy_residual, residual)
            max_axial_drift = max(max_axial_drift, abs(Hz - initial_Hz) / initial_Hnorm)
            max_impulse_balance = max(max_impulse_balance, norm(impulse))
            final_rtn .= (table.target_dv_R_mps[row], table.target_dv_T_mps[row], table.target_dv_N_mps[row])
            for (index, component) in enumerate(("R", "T", "N"))
                positive, negative = positive_columns[index][row], negative_columns[index][row]
                max_monotonicity[index] = max(max_monotonicity[index],
                    check_signed_totals(positive, negative, final_rtn[index],
                        last_positive[index], last_negative[index], component))
                max_sum_residual[index] = max(max_sum_residual[index], abs(final_rtn[index] - positive - negative))
                max_sign_residual[index] = max(max_sign_residual[index], -positive, negative)
                if rows == 1 || no_authority
                    require(final_rtn[index] == positive == negative == 0,
                        "Nonzero initial or no-authority $component diagnostic")
                end
                last_positive[index], last_negative[index] = positive, negative
            end
            no_authority && require(isempty(laser), "Unexpected no-authority firing")
            final_work = work
            previous_time = t
        end
    end
    require(previous_time == duration && rows == info["trajectory_rows"], "Incomplete trajectory")
    require(rows == TRAJECTORY_ROWS, "Wrong row count")
    require(max_axial_drift < 1e-9, "J2 axial angular momentum is not balanced")
    require(max_impulse_balance < 1e-8, "Mass-weighted laser impulse is not balanced")
    require(info["maximum_step_s"] <= STEP_S + 2eps(duration), "Step cap was exceeded")
    for (index, component) in enumerate(("R", "T", "N"))
        if max_monotonicity[index] > BRANCH_TOLERANCE_MPS[component] / 10
            @warn "Small $component branch reversal within the recorded tolerance; values were not clamped" residual_mps=max_monotonicity[index] tolerance_mps=BRANCH_TOLERANCE_MPS[component]
        end
    end
    report = Dict{String,Any}("passed" => true, "verified_utc" => string(now(UTC)),
        "rows" => rows, "encounters" => nrow(geometry), "firing_intervals" => nrow(laser),
        "range_segments" => nrow(segments), "max_J2_energy_work_residual_J" => max_energy_residual,
        "energy_work_absolute_tolerance_J" => max(20.0, Float64(count)),
        "energy_work_relative_tolerance" => 1e-4,
        "max_relative_axial_angular_momentum_drift" => max_axial_drift,
        "max_pair_impulse_residual_kg_mps" => max_impulse_balance,
        "max_boundary_residual_m" => max_boundary_residual,
        "longest_encounter_track_rows" => nrow(tracks),
        "max_signed_T_monotonicity_residual_mps" => max_monotonicity[2],
        "signed_T_tolerance_mps" => SIGNED_T_MONOTONICITY_TOLERANCE_MPS,
        "signed_T_sum_tolerance_mps" => SIGNED_T_SUM_TOLERANCE_MPS,
        "signed_T_monotonicity_warning" => max_monotonicity[2] > 1e-8,
        "RTN_branches" => Dict(component => Dict(
            "max_monotonicity_residual_mps" => max_monotonicity[index],
            "max_sum_to_net_residual_mps" => max_sum_residual[index],
            "max_sign_residual_mps" => max_sign_residual[index],
            "branch_tolerance_mps" => BRANCH_TOLERANCE_MPS[component],
            "sum_tolerance_mps" => SIGNED_T_SUM_TOLERANCE_MPS,
            "final_positive_mps" => last_positive[index],
            "final_negative_mps" => last_negative[index],
            "final_absolute_mps" => last_positive[index] - last_negative[index])
            for (index, component) in enumerate(("R", "T", "N"))),
        "total_delta_sma_m" => final_sma - initial_sma,
        "total_delta_inclination_deg" => final_inc - initial_inc,
        "final_target_RTN_mps" => final_rtn, "final_total_optical_work_J" => final_work,
        "verification_source_signature" => source_signature(source_hashes()))
    atomic_toml(joinpath(directory, "verification.toml"), report)
    info["status"] = "verified"
    atomic_toml(info_path, info)
    println("VERIFIED $directory rows=$rows encounters=$(nrow(geometry)) firing=$(nrow(laser))")
    report
end

end
