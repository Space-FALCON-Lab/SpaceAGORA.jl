module StudyConfig

using Dates, Printf, SHA, TOML

export StudyCase, grid_cases, campaign_settings, case_directory, case_settings, case_signature,
    mission_duration, atomic_toml, filehash, source_hashes, source_signature,
    free_disk_bytes, require_disk_space, peak_rss_bytes, DEFAULT_OUTPUT, LIBRARY_ROOT,
    MU, EARTH_RADIUS_M, J2, MASS_KG, RANGE_M, B, LIGHT_SPEED, STEP_S,
    ENCOUNTER_SAMPLE_S, TRAJECTORY_ROWS, trajectory_times, SCHEMA_VERSION

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
const SCHEMA_VERSION = 2
const ORIGINAL_TARGET_ALTITUDES_KM = (850, 950, 1000, 1050, 1150)

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
        altitude in (ORIGINAL_TARGET_ALTITUDES_KM..., 900, 1100) ||
            throw(ArgumentError("Target altitude is outside the study grid"))
        helper_altitude_km in (900, 1000, 1100) ||
            throw(ArgumentError("Helper altitude is outside the study grid"))
        inclination in (0, 0.5, 1) ||
            throw(ArgumentError("Target inclination is outside the study grid"))
        helpers in (1, 50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800) ||
            throw(ArgumentError("Helper count is outside the study grid"))
        policy in (:gve_sma, :gve_inc) || throw(ArgumentError("Unsupported scheduler"))
        isfinite(orbits) && orbits > 0 || throw(ArgumentError("Orbits must be positive and finite"))
        isfinite(power_W) && power_W >= 0 || throw(ArgumentError("Power must be nonnegative and finite"))
        new(altitude, inclination, helpers, policy, orbits, power_W, helper_altitude_km)
    end
end

function campaign_settings(campaign::AbstractString)
    campaign == "main" && return (
        target_altitudes_km=ORIGINAL_TARGET_ALTITUDES_KM,
        helper_altitudes_km=(1000,),
        helper_counts=(1, 50, 100, 150, 200, 250, 300),
        batch_directory="_batch", progress_file="simulation_progress.log")
    campaign == "helpers400500" && return (
        target_altitudes_km=ORIGINAL_TARGET_ALTITUDES_KM,
        helper_altitudes_km=(1000,),
        helper_counts=(400, 500), batch_directory="_batch_helpers400500",
        progress_file="simulation_progress_helpers400500.log")
    campaign == "helpers9001100km" && return (
        target_altitudes_km=ORIGINAL_TARGET_ALTITUDES_KM,
        helper_altitudes_km=(900, 1100),
        helper_counts=(1, 50, 100, 150, 200, 250, 300, 400, 500),
        batch_directory="_batch_helpers9001100km",
        progress_file="simulation_progress_helpers9001100km.log")
    campaign == "targets9001100km" && return (
        target_altitudes_km=(900, 1100), helper_altitudes_km=(1000,),
        helper_counts=(1, 50, 100, 150, 200, 250, 300, 400, 500),
        batch_directory="_batch_targets9001100km",
        progress_file="simulation_progress_targets9001100km.log")
    campaign == "helpers600800" && return (
        target_altitudes_km=(850, 900, 950, 1000, 1050, 1100, 1150),
        helper_altitudes_km=(1000,), helper_counts=(600, 800),
        batch_directory="_batch_helpers600800",
        progress_file="simulation_progress_helpers600800.log")
    campaign == "helpers700" && return (
        target_altitudes_km=(850, 900, 950, 1000, 1050, 1100, 1150),
        helper_altitudes_km=(1000,), helper_counts=(700,),
        batch_directory="_batch_helpers700",
        progress_file="simulation_progress_helpers700.log")
    throw(ArgumentError("Unknown campaign '$campaign'; use main, helpers400500, helpers600800, helpers700, targets9001100km or helpers9001100km"))
end

function grid_cases(campaign::AbstractString="main")
    selection = campaign_settings(campaign)
    [StudyCase(h, i, n, p; helper_altitude_km=helper_h)
        for helper_h in selection.helper_altitudes_km for h in selection.target_altitudes_km
        for i in (0, 0.5, 1) for n in selection.helper_counts for p in (:gve_sma, :gve_inc)]
end

mission_duration(c::StudyCase) = c.orbits * 2pi * sqrt((EARTH_RADIUS_M + c.altitude_km * 1000)^3 / MU)
trajectory_times(c::StudyCase) = collect(range(0.0, mission_duration(c); length=TRAJECTORY_ROWS))

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
        "diagnostics" => "solver-integrated laser-only RTN, signed along-track parts, optical work",
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
    for root in (joinpath(LIBRARY_ROOT, "src"), @__DIR__)
        for (folder, _, files) in walkdir(root), name in sort(files)
            endswith(name, ".jl") || endswith(name, ".sh") || continue
            path = joinpath(folder, name)
            result[relpath(path, dirname(@__DIR__))] = filehash(path)
        end
    end
    for name in ("Project.toml", "Manifest.toml")
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
