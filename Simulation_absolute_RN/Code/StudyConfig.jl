module StudyConfig

using Dates, Printf, SHA, TOML

export StudyCase, grid_cases, campaign_settings, case_directory, case_settings, case_signature,
    mission_duration, atomic_toml, filehash, source_hashes, source_signature,
    free_disk_bytes, require_disk_space, peak_rss_bytes, DEFAULT_OUTPUT, LIBRARY_ROOT,
    MU, EARTH_RADIUS_M, J2, MASS_KG, RANGE_M, B, LIGHT_SPEED, STEP_S,
    ENCOUNTER_SAMPLE_S, TRAJECTORY_ROWS, trajectory_times, trajectory_columns, SCHEMA_VERSION,
    TARGET_DIAGNOSTIC_COLUMNS, RTN_BRANCH_DEFINITION

const LIBRARY_ROOT = normpath(joinpath(@__DIR__, "..", "..",
    "9_SpaceAGORA.jl-main_full integration"))
const DEFAULT_OUTPUT = normpath(joinpath(@__DIR__, "..", "Output"))
const MU = 3.98600436233e14
const EARTH_RADIUS_M = 6378136.6
const J2 = 0.00108263
const MASS_KG = 227.0
const RANGE_M = 200_000.0
const B = 100.0
const LIGHT_SPEED = 299_792_458.0
const STEP_S = 10.0
const ENCOUNTER_SAMPLE_S = 10.0
const TRAJECTORY_ROWS = 1001
const SCHEMA_VERSION = 3
const TARGET_ALTITUDES_KM = (900, 950, 1000, 1050, 1100)
const HELPER_COUNTS = (1, 100, 200, 300, 400, 500, 600, 700, 800)
const TARGET_DIAGNOSTIC_COLUMNS = (
    :target_dv_R_mps, :target_dv_T_mps, :target_dv_N_mps,
    :target_dv_T_positive_mps, :target_dv_T_negative_mps,
    :target_dv_R_positive_mps, :target_dv_R_negative_mps,
    :target_dv_N_positive_mps, :target_dv_N_negative_mps)
const RTN_BRANCH_DEFINITION =
    "laser-only instantaneous RTN; positive=integral(max(a,0)); negative=integral(min(a,0)); net=positive+negative; absolute=positive-negative"

struct StudyCase
    altitude_km::Float64
    inclination_deg::Float64
    helpers::Int
    policy::Symbol
    orbits::Float64
    power_W::Float64
    helper_altitude_km::Float64
    function StudyCase(altitude, inclination, helpers, policy;
        orbits=100.0, power_W=10_000.0, helper_altitude_km=1000.0)
        altitude in TARGET_ALTITUDES_KM ||
            throw(ArgumentError("Target altitude is outside the study grid"))
        helper_altitude_km == 1000 ||
            throw(ArgumentError("This campaign requires helpers initially at 1000 km"))
        inclination in (0, 0.5, 1) ||
            throw(ArgumentError("Target inclination is outside the study grid"))
        helpers in HELPER_COUNTS ||
            throw(ArgumentError("Helper count is outside the study grid"))
        policy in (:gve_sma, :gve_inc) || throw(ArgumentError("Unsupported scheduler"))
        isfinite(orbits) && orbits > 0 || throw(ArgumentError("Orbits must be positive and finite"))
        isfinite(power_W) && power_W >= 0 || throw(ArgumentError("Power must be nonnegative and finite"))
        new(altitude, inclination, helpers, policy, orbits, power_W, helper_altitude_km)
    end
end

function campaign_settings(campaign::AbstractString)
    campaign == "absolute_rn" || throw(ArgumentError("Unknown campaign '$campaign'; use absolute_rn"))
    (target_altitudes_km=TARGET_ALTITUDES_KM, helper_altitudes_km=(1000,),
        helper_counts=HELPER_COUNTS, batch_directory="_batch_absolute_RN",
        progress_file="simulation_progress_absolute_RN.log")
end

function grid_cases(campaign::AbstractString="absolute_rn")
    selection = campaign_settings(campaign)
    [StudyCase(h, i, n, p; helper_altitude_km=helper_h)
        for helper_h in selection.helper_altitudes_km for h in selection.target_altitudes_km
        for i in (0, 0.5, 1) for n in selection.helper_counts for p in (:gve_sma, :gve_inc)]
end

mission_duration(c::StudyCase) = c.orbits * 2pi * sqrt((EARTH_RADIUS_M + c.altitude_km * 1000)^3 / MU)
trajectory_times(c::StudyCase) = collect(range(0.0, mission_duration(c); length=TRAJECTORY_ROWS))

function trajectory_columns(spacecraft_count)
    columns = Symbol[:time]
    for sc in 1:spacecraft_count
        append!(columns, [Symbol("sc$(sc)_$(field)_$axis")
            for field in ("pos", "vel", "laser_dv") for axis in 1:3])
        push!(columns, Symbol("sc$(sc)_laser_work_J"))
    end
    append!(columns, TARGET_DIAGNOSTIC_COLUMNS)
    columns
end

function case_directory(root, c::StudyCase)
    joinpath(abspath(root), "paper_plot_mode",
        @sprintf("h%.0fkm_t%.0fkm", c.helper_altitude_km, c.altitude_km),
        @sprintf("ih0.0deg_it%.1fdeg", c.inclination_deg), "N$(c.helpers)",
        "T$(round(Int, mission_duration(c)))s", "$(c.policy)_e0.0000_nu1.0000")
end

function case_settings(c::StudyCase)
    Dict{String,Any}(
        "schema_version" => SCHEMA_VERSION, "target_altitude_km" => c.altitude_km,
        "target_inclination_deg" => c.inclination_deg, "helper_altitude_km" => c.helper_altitude_km,
        "helper_inclination_deg" => 0.0, "helper_count" => c.helpers,
        "target_index" => c.helpers + 1, "target_id" => c.helpers + 1,
        "target_true_anomaly_deg" => 1.0, "helper_phase_rule" => "360*(k-1)/N degrees",
        "eccentricity" => 0.0, "RAAN_deg" => 0.0, "argument_of_periapsis_deg" => 0.0,
        "initial_target_periods" => c.orbits, "duration_s" => mission_duration(c),
        "mass_kg" => MASS_KG, "laser_power_W" => c.power_W, "magnification_B" => B,
        "speed_of_light_mps" => LIGHT_SPEED, "force_N" => B * c.power_W / LIGHT_SPEED,
        "range_m" => RANGE_M, "terminals_per_spacecraft" => 1,
        "link_scope" => "target-helper", "encounter_scope" => "all spacecraft pairs",
        "scheduler" => string(c.policy), "scheduler_updates" => "initialization and accepted steps",
        "range_rule" => "hold selected links between accepted steps, including outside-range firing",
        "pointing_acquisition_losses" => false, "J2" => J2, "mu_m3_s2" => MU,
        "earth_equatorial_radius_m" => EARTH_RADIUS_M, "drag" => false,
        "epoch" => "2026-01-01T00:00:00", "frame" => "J2000 ECI",
        "solver" => "Tsit5", "reltol" => 1e-12, "abstol" => 1e-12,
        "dt_max_s" => STEP_S, "save_interval_s" => mission_duration(c) / (TRAJECTORY_ROWS - 1),
        "trajectory_rows" => TRAJECTORY_ROWS, "encounter_sample_interval_s" => ENCOUNTER_SAMPLE_S,
        "encounter_extrema" => "10 s samples plus exact boundary states; sampled estimates",
        "longest_tracks" => "longest all-pair encounter and longest target-helper encounter, deduplicated",
        "units" => "SI: m, s, kg, m/s, J; orbital setup angles in degrees",
        "diagnostics" => "solver-integrated laser-only RTN, positive/negative R/T/N parts, optical work",
        "RTN_branch_definition" => RTN_BRANCH_DEFINITION,
        "target_diagnostic_columns" => string.(collect(TARGET_DIAGNOSTIC_COLUMNS)),
        "diagnostic_initial_values" => 0.0,
        "compression" => "Arrow IPC file / Feather V2, zstd, chunked Float64")
end

function case_signature(c::StudyCase)
    io = IOBuffer()
    TOML.print(io, case_settings(c); sorted=true)
    bytes2hex(sha256(take!(io)))
end

filehash(path) = open(io -> bytes2hex(sha256(io)), path)

function source_hashes()
    result = Dict{String,String}()
    for root in (joinpath(LIBRARY_ROOT, "src"), joinpath(LIBRARY_ROOT, "ext"), @__DIR__)
        for (folder, _, files) in walkdir(root), name in sort(files)
            endswith(name, ".jl") || endswith(name, ".sh") || continue
            path = joinpath(folder, name)
            result[relpath(path, dirname(@__DIR__))] = filehash(path)
        end
    end
    for name in ("Project.toml", "Manifest.toml", "Artifacts.toml", "CondaPkg.toml", "examples/common.jl")
        path = joinpath(LIBRARY_ROOT, name)
        isfile(path) && (result[relpath(path, dirname(@__DIR__))] = filehash(path))
    end
    result
end

source_signature(hashes) = bytes2hex(sha256(join(
    ["$key=$(hashes[key])" for key in sort!(collect(keys(hashes)))], "\n")))

function atomic_toml(path, values)
    temporary = path * ".tmp"
    open(temporary, "w") do io
        TOML.print(io, values; sorted=true)
        flush(io)
    end
    mv(temporary, path; force=true)
    nothing
end

function free_disk_bytes(path)
    text = readchomp(`df --output=avail -B1 $path`)
    parse(Int, strip(last(split(text, '\n'))))
end

function require_disk_space(path; reserve=5 * 1024^3)
    available = free_disk_bytes(path)
    available >= reserve || error(
        "Insufficient disk space at $path: $available bytes free; reserve is $reserve. No files were deleted.")
    available
end

function peak_rss_bytes()
    line = only(filter(line -> startswith(line, "VmHWM:"), readlines("/proc/self/status")))
    parse(Int, split(line)[2]) * 1024
end

end
