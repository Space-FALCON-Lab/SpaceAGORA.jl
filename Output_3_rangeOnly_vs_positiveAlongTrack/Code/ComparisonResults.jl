module ComparisonResults

include("CampaignPaths.jl")
import Pkg
Pkg.activate(CampaignPaths.LIBRARY_ROOT; io=devnull)
include(joinpath(CampaignPaths.BASE_CODE, "plot_orbital_elements.jl"))
include(joinpath(CampaignPaths.PAPER_CODE, "plot_paper_results.jl"))
include(joinpath(CampaignPaths.BASE_CODE, "evaluate_encounters.jl"))
using Arrow, CSV, DataFrames, Dates, LinearAlgebra, SHA, StaticArrays, TOML
using Plots

const OE = OrbitalElementPlots
const SC = OE.StudyConfig
const PP = PaperPlots
const OE_FIELDS = (:da_m, :de, :di_deg, :draan_deg)
const RTN_FIELDS = (OE.RTN_FIELDS..., :dv_R_absolute_mps, :dv_T_absolute_mps, :dv_N_absolute_mps)
const VELOCITY_FIELDS = (:velocity_difference_R_mps, :velocity_difference_T_mps,
    :velocity_difference_N_mps, :velocity_difference_norm_mps)
const PAIR_EXCEPTIONS = Set(("scheduler", "laser_power_W", "force_N"))
const POLICY_LABELS = (naive_next_entering="Range-only",
    positive_along_track="Along-track+", no_laser="No laser")
const COLORS = (naive_next_entering=PP.LINE_COLORS.gve_sma,
    positive_along_track=PP.LINE_COLORS.gve_inc, no_laser=PP.LINE_COLORS.gve_sma)
const STYLES = (naive_next_entering=PP.LINE_STYLES.gve_sma,
    positive_along_track=PP.LINE_STYLES.gve_inc, no_laser=PP.LINE_STYLES.gve_sma)
const IDENTIFIERS = (:scheduler, :target_altitude_km, :delta_i_deg, :helpers,
    :orbits, :time_s, :row, :run_id)
const AVERAGE_FIELDS = (:aggregation, :averaging_orbits, :window_start_orbits,
    :window_end_orbits, :window_start_row, :window_end_row, :window_start_time_s,
    :window_end_time_s, :sample_count)

function unwrap_raan(elements)
    angles = [element.raan for element in elements]
    if isnan(first(angles))
        all(isnan, angles) || error("RAAN became defined after an undefined initial node")
        return angles
    end
    all(isfinite, angles) || error("RAAN became undefined during a non-equatorial run")
    result = copy(angles)
    for index in 2:length(result)
        result[index] = result[index - 1] + OE.signed_angle(angles[index] - angles[index - 1])
    end
    result
end

function read_target(directory, case)
    info = TOML.parsefile(joinpath(directory, "run_info.toml"))
    verification = TOML.parsefile(joinpath(directory, "verification.toml"))
    info["status"] == "verified" && verification["passed"] ||
        error("Unverified case: $directory")
    info["settings"] == SC.case_settings(case) && info["run_id"] == SC.case_signature(case) ||
        error("Case settings or run ID mismatch: $directory")
    info["record_events"] || error("Encounter recording was disabled: $directory")
    info["trajectory_rows"] == SC.TRAJECTORY_ROWS ||
        error("Unexpected trajectory row count: $directory")
    path = joinpath(directory, "trajectory.feather")
    SC.filehash(path) == info["checksums"]["trajectory.feather"] ||
        error("Trajectory checksum mismatch: $path")
    times = SC.trajectory_times(case)
    positions, velocities = zeros(3, length(times)), zeros(3, length(times))
    rtn = zeros(length(RTN_FIELDS), length(times))
    target = case.helpers + 1
    previous = NamedTuple{RTN_FIELDS}(ntuple(_ -> 0.0, length(RTN_FIELDS)))
    max_reversal = zeros(3)
    row = 0
    for batch in Arrow.Stream(path)
        collect(propertynames(batch)) == SC.trajectory_columns(target) ||
            error("Trajectory schema mismatch: $path")
        metadata = Arrow.getmetadata(batch)
        metadata !== nothing && get(metadata, "run_id", "") == info["run_id"] &&
            get(metadata, "schema_version", "") == string(SC.SCHEMA_VERSION) &&
            get(metadata, "target_id", "") == string(target) ||
            error("Trajectory metadata mismatch: $path")
        for index in eachindex(batch.time)
            row += 1
            row <= length(times) && batch.time[index] == times[row] ||
                error("Time-grid mismatch at row $row: $path")
            for axis in 1:3
                positions[axis, row] = batch[Symbol("sc$(target)_pos_$axis")][index]
                velocities[axis, row] = batch[Symbol("sc$(target)_vel_$axis")][index]
            end
            values = OE.rtn_values(batch, index)
            max_reversal .= max.(max_reversal, OE.check_rtn(values, previous; context=path))
            rtn[:, row] .= Tuple(values)
            previous = values
        end
    end
    row == length(times) || error("Incomplete trajectory: $path")
    all(isfinite, positions) && all(isfinite, velocities) || error("Nonfinite target state: $path")
    initial = info["initial_states"][target]
    positions[:, 1] == initial["position_m"] && velocities[:, 1] == initial["velocity_mps"] ||
        error("Initial target state mismatch: $path")
    all(iszero, rtn[:, 1]) || error("Nonzero initial laser diagnostics: $path")
    if case.policy === :no_laser
        all(iszero, rtn) && info["firing_intervals"] == 0 && info["actual_firing_s"] == 0 ||
            error("No-laser case contains laser diagnostics or firing: $path")
    end
    any(max_reversal .> (1e-9, 1e-8, 1e-9)) &&
        @warn "Small RTN branch reversals within existing tolerances; no clamping" path max_reversal
    elements = [OE.orbital_elements(SVector{3}(positions[:, index]),
        SVector{3}(velocities[:, index]), SC.MU) for index in eachindex(times)]
    (; directory, case, info, time_s=times, positions, velocities, rtn, elements,
        raan=unwrap_raan(elements), source_trajectory=path)
end

function identifiers(track, window, index)
    case = track.case
    (; scheduler=string(case.policy), target_altitude_km=case.altitude_km,
        delta_i_deg=case.inclination_deg, helpers=case.helpers,
        orbits=window, time_s=track.time_s[index], row=index, run_id=track.info["run_id"],
        source_trajectory=track.source_trajectory,
        source_trajectory_sha256=track.info["checksums"]["trajectory.feather"])
end

function endpoint_indices(track, windows)
    !isempty(windows) && allunique(windows) &&
        all(window -> 0 < window <= track.case.orbits, windows) ||
        error("Orbit windows must be distinct and within the run duration")
    period = SC.mission_duration(track.case) / track.case.orbits
    [OE.exact_sample_index(track.time_s, window * period) for window in windows]
end

function absolute_endpoints(track, windows)
    initial = first(track.elements)
    [merge(identifiers(track, window, index),
        (da_m=track.elements[index].a - initial.a,
            de=track.elements[index].e - initial.e,
            di_deg=rad2deg(track.elements[index].i - initial.i),
            draan_deg=rad2deg(track.raan[index] - first(track.raan)),
            raan_defined=isfinite(track.raan[index])),
        NamedTuple{RTN_FIELDS}(Tuple(track.rtn[:, index])))
        for (window, index) in zip(windows, endpoint_indices(track, windows))]
end

function validate_pair(laser, baseline)
    laser.case.policy in (:naive_next_entering, :positive_along_track) ||
        error("Expected a laser-enabled comparison policy")
    baseline.case.policy === :no_laser && baseline.case.power_W == 0 &&
        laser.case.power_W == 10_000 || error("Invalid on/off laser powers or baseline policy")
    on_settings, off_settings = laser.info["settings"], baseline.info["settings"]
    filter(pair -> first(pair) ∉ PAIR_EXCEPTIONS, on_settings) ==
        filter(pair -> first(pair) ∉ PAIR_EXCEPTIONS, off_settings) ||
        error("Laser and baseline physical/numerical settings do not match")
    laser.info["initial_states"] == baseline.info["initial_states"] ||
        error("Laser and baseline initial spacecraft states do not match")
    laser.info["source_hashes"] == baseline.info["source_hashes"] ||
        error("Laser and baseline were generated from different source code")
    laser.time_s == baseline.time_s ||
        error("Laser and baseline times do not match; interpolation is not permitted")
    all(iszero, baseline.rtn) && baseline.info["firing_intervals"] == 0 ||
        error("The no-laser baseline has nonzero laser diagnostics")
    nothing
end

function velocity_difference_in_baseline_rtn(on_velocity, off_position, off_velocity)
    radial = normalize(off_position)
    normal = normalize(cross(off_position, off_velocity))
    transverse = cross(normal, radial)
    difference = on_velocity - off_velocity
    result = SVector(dot(difference, radial), dot(difference, transverse),
        dot(difference, normal), norm(difference))
    all(isfinite, result) || error("Cannot form a finite baseline RTN velocity difference")
    result
end

function difference_series(laser, baseline)
    validate_pair(laser, baseline)
    count = length(laser.time_s)
    elements = zeros(4, count)
    velocity = zeros(4, count)
    for index in 1:count
        on, off = laser.elements[index], baseline.elements[index]
        elements[:, index] .= (on.a - off.a, on.e - off.e, rad2deg(on.i - off.i),
            rad2deg(laser.raan[index] - baseline.raan[index]))
        velocity[:, index] .= velocity_difference_in_baseline_rtn(
            SVector{3}(laser.velocities[:, index]),
            SVector{3}(baseline.positions[:, index]), SVector{3}(baseline.velocities[:, index]))
    end
    rtn = laser.rtn - baseline.rtn
    all(iszero, velocity[:, 1]) && all(iszero, elements[1:3, 1]) ||
        error("A paired difference is nonzero at initialization")
    (; elements, rtn, velocity)
end

function save_difference(laser, baseline, series, root)
    directory = SC.case_directory(root, laser.case)
    mkpath(directory)
    trajectory = joinpath(directory, "laser_minus_no_laser.feather")
    run_id = bytes2hex(sha256(join((laser.info["run_id"], baseline.info["run_id"],
        laser.info["checksums"]["trajectory.feather"],
        baseline.info["checksums"]["trajectory.feather"]), "\n")))
    table = DataFrame(merge(
        (time_s=laser.time_s, reference_orbits=laser.time_s .* (laser.case.orbits / last(laser.time_s)),
            row=collect(eachindex(laser.time_s)), raan_defined=isfinite.(series.elements[4, :])),
        NamedTuple{OE_FIELDS}(Tuple(series.elements[k, :] for k in 1:4)),
        NamedTuple{RTN_FIELDS}(Tuple(series.rtn[k, :] for k in eachindex(RTN_FIELDS))),
        NamedTuple{VELOCITY_FIELDS}(Tuple(series.velocity[k, :] for k in 1:4))))
    Arrow.write(trajectory * ".partial", table; file=true, compress=:zstd,
        metadata=Dict("run_id" => run_id, "comparison" => "laser_on_minus_no_laser",
            "velocity_difference_frame" => "instantaneous no-laser target RTN",
            "units" => "time:s; da:m; de:1; di,draan:deg; RTN,velocity:m/s"))
    mv(trajectory * ".partial", trajectory; force=true)
    info = Dict("status" => "complete", "comparison_id" => run_id,
        "rows" => nrow(table), "difference_order" => "laser_on_minus_no_laser",
        "laser_run_id" => laser.info["run_id"], "no_laser_run_id" => baseline.info["run_id"],
        "laser_trajectory" => laser.source_trajectory, "no_laser_trajectory" => baseline.source_trajectory,
        "laser_trajectory_sha256" => laser.info["checksums"]["trajectory.feather"],
        "no_laser_trajectory_sha256" => baseline.info["checksums"]["trajectory.feather"],
        "difference_trajectory_sha256" => SC.filehash(trajectory),
        "identical_initial_states" => true, "exact_time_alignment" => true,
        "laser_RTN_baseline_is_zero" => true,
        "velocity_difference_definition" => "project v_on(t)-v_off(t) into the same RTN basis built from r_off(t),v_off(t); includes induced gravitational trajectory changes",
        "OE_definition" => "osculating OE_on(t)-OE_off(t); RAAN independently unwrapped over all saved rows; equatorial RAAN is undefined",
        "created_utc" => string(now(UTC)))
    SC.atomic_toml(joinpath(directory, "comparison_info.toml"), info)
    trajectory
end

function difference_endpoints(laser, baseline, series, trajectory, windows)
    [merge(identifiers(laser, window, index),
        (source_no_laser_trajectory=baseline.source_trajectory,
            source_no_laser_trajectory_sha256=baseline.info["checksums"]["trajectory.feather"],
            source_no_laser_run_id=baseline.info["run_id"], source_difference_trajectory=trajectory,
            raan_defined=isfinite(series.elements[4, index])),
        NamedTuple{OE_FIELDS}(Tuple(series.elements[:, index])),
        NamedTuple{RTN_FIELDS}(Tuple(series.rtn[:, index])),
        NamedTuple{VELOCITY_FIELDS}(Tuple(series.velocity[:, index])))
        for (window, index) in zip(windows, endpoint_indices(laser, windows))]
end

function time_average(times, values)
    length(times) == length(values) >= 2 ||
        error("A time average requires at least two matching time/value samples")
    all(isfinite, times) && all(>(0), diff(times)) ||
        error("Averaging times must be finite and strictly increasing")
    # Equatorial RAAN is undefined over the entire window, not zero.
    all(isnan, values) && return NaN
    all(isfinite, values) || error("Cannot average partially undefined or nonfinite values")
    sum((times[k + 1] - times[k]) * (values[k] + values[k + 1]) / 2
        for k in 1:length(times)-1) / (last(times) - first(times))
end

function orbit_averages(track, windows; baseline=nothing)
    all(window -> window >= 2, windows) ||
        error("The last-two-orbit window must not start before initialization")
    endpoints = endpoint_indices(track, windows)
    period = SC.mission_duration(track.case) / track.case.orbits
    series = if baseline === nothing
        initial = first(track.elements)
        elements = reduce(hcat, [SVector(element.a - initial.a, element.e - initial.e,
            rad2deg(element.i - initial.i), rad2deg(track.raan[k] - first(track.raan)))
            for (k, element) in enumerate(track.elements)])
        (; elements, rtn=track.rtn)
    else
        difference_series(track, baseline)
    end
    rows = NamedTuple[]
    for (window, stop) in zip(windows, endpoints)
        start = OE.exact_sample_index(track.time_s, (window - 2) * period)
        indices = start:stop
        times = @view track.time_s[indices]
        average_fields(fields, matrix) = NamedTuple{fields}(Tuple(
            time_average(times, @view matrix[k, indices]) for k in eachindex(fields)))
        oe = average_fields(OE_FIELDS, series.elements)
        row = merge(identifiers(track, window, stop),
            (aggregation="time_weighted_mean", averaging_orbits=2,
                window_start_orbits=window - 2, window_end_orbits=window,
                window_start_row=start, window_end_row=stop,
                window_start_time_s=first(times), window_end_time_s=last(times),
                sample_count=length(indices), raan_defined=isfinite(oe.draan_deg)),
            oe, average_fields(RTN_FIELDS, series.rtn))
        if baseline !== nothing
            row = merge(row,
                (source_no_laser_trajectory=baseline.source_trajectory,
                    source_no_laser_trajectory_sha256=baseline.info["checksums"]["trajectory.feather"],
                    source_no_laser_run_id=baseline.info["run_id"]),
                average_fields(VELOCITY_FIELDS, series.velocity))
        end
        push!(rows, row)
    end
    rows
end

function access_thresholds()
    [merge((target_altitude_km=altitude, delta_i_deg=inclination,
        helper_radius_m=SC.EARTH_RADIUS_M + 1000e3,
        target_radius_m=SC.EARTH_RADIUS_M + altitude * 1000, range_m=SC.RANGE_M),
        OE.continuous_access(SC.EARTH_RADIUS_M + 1000e3,
            SC.EARTH_RADIUS_M + altitude * 1000, SC.RANGE_M, deg2rad(inclination)))
        for altitude in PP.ALTITUDES for inclination in PP.INCLINATIONS]
end

function plot_spec(kind, window; difference=false, averaged=false)
    if averaged
        spec = plot_spec(kind, window; difference)
        return merge(spec, (csv="mean_last2orbits_" * spec.csv,
            stem="mean_last2orbits_" * spec.stem))
    end
    if kind == :velocity
        return (csv="velocity_difference_RTN_$(window)orbits.csv",
            stem="velocity_difference_RTN_laser_minus_no_laser_$(window)orbits_h_cols",
            fields=VELOCITY_FIELDS, scales=(1.0, 1.0, 1.0, 1.0),
            labels=("\$\\delta v_R\\; (\\mathrm{m/s})\$", "\$\\delta v_T\\; (\\mathrm{m/s})\$",
                "\$\\delta v_N\\; (\\mathrm{m/s})\$", "\$\\Vert\\delta\\mathbf{v}\\Vert\\; (\\mathrm{m/s})\$"))
    end
    spec = getproperty(PP.SPECS, kind)
    prefix = difference ? "laser_minus_no_laser_" : ""
    labels = if difference && kind == :oe
        ("\$\\delta a_{\\mathrm{on-off}}\\; (\\mathrm{km})\$",
            "\$\\delta e_{\\mathrm{on-off}}\\quad (\\times10^{-4})\$",
            "\$\\delta i_{\\mathrm{on-off}}\\quad (\\times10^{-3}\\,\\mathrm{deg})\$",
            "\$\\delta\\Omega_{\\mathrm{on-off}}\\; (\\mathrm{deg})\$")
    elseif difference
        ("\$\\delta(\\Delta v_R)\\; (\\mathrm{m/s})\$", "\$\\delta(\\Delta v_T)\\; (\\mathrm{m/s})\$",
            "\$\\delta(\\Delta v_{N+})\\; (\\mathrm{m/s})\$",
            "\$\\delta(\\Delta v_{N-})\\; (\\mathrm{m/s})\$")
    else
        spec.labels
    end
    merge(spec, (csv=prefix * replace(spec.csv, "100orbits" => "$(window)orbits"),
        stem=prefix * replace(spec.stem, "100orbits" => "$(window)orbits"), labels))
end

function plot_limits(rows, spec)
    NamedTuple{spec.fields}(Tuple(
        field in OE_FIELDS ? OE.ylimits([getproperty(row, field) for row in rows], field) :
            OE.dv_ylimits([getproperty(row, field) for row in rows], field)
        for field in spec.fields))
end

function write_figures(rows, root, windows; difference=false, averaged=false)
    mkpath(root)
    summary = averaged ? "last_two_orbits_average_summary.csv" : "production_endpoint_summary.csv"
    CSV.write(joinpath(root, summary), rows)
    outputs = String[]
    for window in windows
        selected = filter(row -> row.orbits == window, rows)
        isempty(selected) && error("Missing $window-orbit results")
        policies = Tuple(unique(row.scheduler for row in selected))
        directory = averaged ?
            joinpath(root, "images", "last_two_orbits_average", "$(window)orbits") :
            joinpath(root, "images", window == 100 ? "for_paper" : "$(window)orbits")
        mkpath(directory)
        CSV.write(joinpath(directory, "continuous_access_thresholds.csv"), access_thresholds())
        for kind in (difference ? (:oe, :dv, :velocity) : (:oe, :dv))
            spec = plot_spec(kind, window; difference, averaged)
            source_fields = Tuple(field for field in propertynames(first(selected))
                if startswith(string(field), "source_"))
            fields = (IDENTIFIERS..., (averaged ? AVERAGE_FIELDS : ())..., source_fields..., spec.fields...,
                (kind == :oe ? (:raan_defined,) : ())...)
            CSV.write(joinpath(directory, spec.csv),
                [NamedTuple{fields}(Tuple(getproperty(row, field) for field in fields)) for row in selected])
            PP.make_plot(kind == :oe ? :oe : :dv; directory, policies,
                policy_labels=POLICY_LABELS, colors=COLORS, styles=STYLES, specification=spec,
                orbit_window=window, endpoint_row=round(Int, 10window) + 1,
                averaging_orbits=averaged ? 2 : 0,
                source_limits=plot_limits(selected, spec))
            append!(outputs, [joinpath(directory, "$(spec.stem).$extension") for extension in ("png", "pdf")])
        end
        paths = filter(path -> dirname(path) == directory, outputs)
        all(path -> isfile(path) && filesize(path) > 10_000, paths) ||
            error("A paper figure is missing or empty")
        manifest = Dict(
            "status" => "complete", "orbits" => window, "cases" => length(selected),
            "schedulers" => collect(policies), "difference" => difference,
            "OE_definition" => difference ? "OE_laser(t)-OE_no_laser(t)" : "OE(t)-OE(0)",
            "RAAN_definition" => "independently unwrapped at every saved row; undefined equatorial RAAN omitted",
            "laser_RTN_definition" => difference ?
                "laser-only RTN integrals minus zero no-laser integrals; not a state-velocity difference" :
                "laser-only integrals in instantaneous target RTN; signed R/T and positive/negative N",
            "state_velocity_difference_definition" =>
                "v_laser(t)-v_no_laser(t), projected into the no-laser target RTN basis; norm is the ECI difference magnitude",
            "figure_checksums" => Dict(basename(path) => SC.filehash(path) for path in paths))
        if averaged
            merge!(manifest, Dict(
                "aggregation" => "time_weighted_mean", "averaging_orbits" => 2,
                "window_start_orbits" => window - 2, "window_end_orbits" => window,
                "method" => "trapezoidal integration of saved samples divided by elapsed window time",
                "period_definition" => "initial Keplerian target period; not counted revolutions",
                "operation_order" => "compute instantaneous changes/differences and RTN projections first, then average",
                "sample_count" => unique(row.sample_count for row in selected),
                "csv_checksums" => Dict(file => SC.filehash(joinpath(directory, file))
                    for file in readdir(directory) if endswith(file, ".csv"))))
        end
        SC.atomic_toml(joinpath(directory, "figure_manifest.toml"), manifest)
    end
    outputs
end

function render_campaign(campaign, root; windows=[80, 100])
    campaign in ("range_along", "no_laser") || error("Unknown comparison campaign: $campaign")
    rows, summaries = NamedTuple[], NamedTuple[]
    for (index, case) in enumerate(SC.grid_cases(campaign))
        track = read_target(SC.case_directory(root, case), case)
        append!(rows, absolute_endpoints(track, windows))
        push!(summaries, (scheduler=string(case.policy), target_altitude_km=case.altitude_km,
            target_inclination_deg=case.inclination_deg, helper_count=case.helpers,
            run_id=track.info["run_id"], source_directory=track.directory,
            trajectory_sha256=track.info["checksums"]["trajectory.feather"],
            encounters=track.info["encounters"], firing_intervals=track.info["firing_intervals"],
            actual_firing_s=track.info["actual_firing_s"], output_bytes=track.info["output_bytes"]))
        println("Read $campaign case $index")
    end
    CSV.write(joinpath(root, "production_case_summary.csv"), summaries)
    EncounterEvaluation.main(["--campaign", campaign, "--input", root,
        "--output", joinpath(root, "encounter_evaluation")])
    outputs = write_figures(rows, root, windows)
    println("Completed $campaign: $(length(summaries)) verified cases, $(length(outputs)) PNG/PDF files.")
    outputs
end

function render_differences(laser_root, baseline_root, output; windows=[80, 100])
    rows = NamedTuple[]
    for (index, baseline_case) in enumerate(SC.grid_cases("no_laser"))
        baseline = read_target(SC.case_directory(baseline_root, baseline_case), baseline_case)
        for policy in SC.campaign_settings("range_along").policies
            case = SC.StudyCase(baseline_case.altitude_km, baseline_case.inclination_deg,
                baseline_case.helpers, policy)
            laser = read_target(SC.case_directory(laser_root, case), case)
            series = difference_series(laser, baseline)
            trajectory = save_difference(laser, baseline, series, output)
            append!(rows, difference_endpoints(laser, baseline, series, trajectory, windows))
        end
        println("Paired baseline $index / 81 with both laser schedulers")
    end
    outputs = write_figures(rows, output, windows; difference=true)
    println("Completed 162 exact-time laser-minus-no-laser comparisons, $(length(outputs)) PNG/PDF files.")
    outputs
end

function render_last_two_orbits(laser_root, baseline_root, difference_root; windows=[80, 100])
    laser_rows, baseline_rows, difference_rows = NamedTuple[], NamedTuple[], NamedTuple[]
    for (index, baseline_case) in enumerate(SC.grid_cases("no_laser"))
        baseline = read_target(SC.case_directory(baseline_root, baseline_case), baseline_case)
        append!(baseline_rows, orbit_averages(baseline, windows))
        for policy in SC.campaign_settings("range_along").policies
            case = SC.StudyCase(baseline_case.altitude_km, baseline_case.inclination_deg,
                baseline_case.helpers, policy)
            laser = read_target(SC.case_directory(laser_root, case), case)
            append!(laser_rows, orbit_averages(laser, windows))
            append!(difference_rows, orbit_averages(laser, windows; baseline))
        end
        println("Averaged geometry $index / 81 for all three outputs")
        flush(stdout)
    end
    outputs = vcat(write_figures(laser_rows, laser_root, windows; averaged=true),
        write_figures(baseline_rows, baseline_root, windows; averaged=true),
        write_figures(difference_rows, difference_root, windows; difference=true, averaged=true))
    println("Completed last-two-orbit averages: $(length(outputs)) PNG/PDF files.")
    outputs
end

end
