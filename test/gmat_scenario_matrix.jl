using Test
using TOML
using Arrow
using CSV
using DataFrames
using Statistics
using LinearAlgebra
using Dates
using Plots
using StaticArrays

const _GMAT_REPO_ROOT = isdefined(Main, :REPO_ROOT) ? Main.REPO_ROOT : normpath(joinpath(@__DIR__, ".."))

if !isdefined(@__MODULE__, :SimulationModel)
    include(joinpath(_GMAT_REPO_ROOT, "src", "core", "simulation_model.jl"))
end

if !isdefined(@__MODULE__, :SimulationEngine)
    include(joinpath(_GMAT_REPO_ROOT, "src", "simulation", "engine", "simulation_engine.jl"))
end

if !isdefined(@__MODULE__, :TelemetryVerification)
    include(joinpath(_GMAT_REPO_ROOT, "src", "analysis", "verification", "telemetry_verification.jl"))
end

const TV = TelemetryVerification
const SM = SimulationModel

# Basilisk parity references and STK results are not tracked here; the
# references are synced from the lab-org verification-data repo by
# scripts/dev/fetch_private_telemetry.sh (see data/telemetry/PRIVATE_TELEMETRY.md).
# The parity testsets skip cleanly when they are absent.
const _BASILISK_REFERENCE_DIR = joinpath(
    _GMAT_REPO_ROOT,
    "data",
    "telemetry",
    "Basilisk_Examples_Full"
)
_basilisk_reference_available() = isdir(_BASILISK_REFERENCE_DIR)
const _STK_RESULTS_DIR = joinpath(
    _GMAT_REPO_ROOT,
    "data",
    "telemetry",
    "stk_results"
)
_stk_reference_available() = isdir(_STK_RESULTS_DIR)

const _GMAT_HARMONICS_EARTH_FILE = "data/Gravity_harmonics_data/EarthGGM05C.csv" # For internal GMAT parity, matches the file used in the GMAT scenarios
# Must match data/telemetry/gmat_matrix_parity_locked.py, which generated the
# locked Mars references with Mars-50c (GMM2B.csv stays available as an asset).
const _GMAT_HARMONICS_MARS_FILE = "data/Gravity_harmonics_data/Mars50c.csv"
const _GMAT_HARMONICS_VENUS_FILE = "data/Gravity_harmonics_data/MGNP180U.csv"
const _GMAT_HARMONICS_MOON_FILE = "data/Gravity_harmonics_data/LP165P.csv"
const _GMAT_PLANETARY_KERNEL_CANDIDATES = (
    "spk/planets/de430.bsp",
)
const _CYGNSS_48HR_TELEMETRY_FEATHER = joinpath(_GMAT_REPO_ROOT, "data", "telemetry", "CYGNSS", "cygnss_data_48hr.feather")
const _CYGNSS_96HR_TELEMETRY_FEATHER = joinpath(_GMAT_REPO_ROOT, "data", "telemetry", "CYGNSS", "cyg04_nasa_pvt_96hr.feather")
const _CYGNSS_CYG04_96HR_TELEMETRY_FEATHER = joinpath(_GMAT_REPO_ROOT, "data", "telemetry", "CYGNSS", "cyg04_nasa_pvt_96hr.feather")
# CYGNSS flight telemetry is access-restricted and gitignored (see
# data/telemetry/PRIVATE_TELEMETRY.md). The CYGNSS testsets below skip cleanly
# when it has not been synced into this checkout.
_cygnss_private_data_available() =
    isfile(_CYGNSS_48HR_TELEMETRY_FEATHER) && isfile(_CYGNSS_CYG04_96HR_TELEMETRY_FEATHER)
const _CYGNSS_GMAT_COMPARISON_PATH = let
    basilisk_path = joinpath(_BASILISK_REFERENCE_DIR, "Sim_CYGNSS_Comparison.feather")
    isfile(basilisk_path) ? basilisk_path : joinpath(_GMAT_REPO_ROOT, "data", "telemetry", "GMAT_Examples", "Sim_CYGNSS_Comparison.feather")
end
const _GMAT_MATRIX_EXCLUDED_FILES = Set([
    "Sim_CYGNSS_Comparison.feather",
    "Sim_Hubble_Comp.feather"
])
const _CYGNSS_96HR_J2000_CACHE = Ref{Any}(nothing)
const _CYGNSS_CYG04_96HR_INERTIAL_CACHE = Ref{Any}(nothing)

const TEST_MODE::Symbol = :quick

@inline function _parse_bool_env(name::String, default::Bool)::Bool
    token = lowercase(strip(get(ENV, name, default ? "1" : "0")))
    if token in ("1", "true", "yes", "on")
        return true
    elseif token in ("0", "false", "no", "off")
        return false
    end
    return default
end

@inline function _gmat_parity_solver_mode_enabled()::Bool
    return _parse_bool_env("SPACEAGORA_GMAT_PARITY_SOLVER", false)
end

function _gmat_planetary_kernel_relpath()::String
    for relpath in _GMAT_PLANETARY_KERNEL_CANDIDATES
        if isfile(joinpath(TV.SPICE_PATH, relpath))
            return relpath
        end
    end
    throw(ArgumentError("Unable to find a supported planetary SPICE kernel in $(TV.SPICE_PATH). Tried: $(join(_GMAT_PLANETARY_KERNEL_CANDIDATES, ", "))"))
end

function _telemetry_solver_env_overrides()::Dict{String, String}
    if _gmat_parity_solver_mode_enabled()
        return Dict(
            "SPACEAGORA_TELEMETRY_SOLVER_MODE" => "auto_stiff",
            "SPACEAGORA_TELEMETRY_DT_MAX_ORBIT" => "1.0",
            "SPACEAGORA_TELEMETRY_RELTOL_ORBIT" => "1e-13",
            "SPACEAGORA_TELEMETRY_ABSTOL_ORBIT" => "1e-13",
            "SPACEAGORA_TELEMETRY_RELTOL_ATM" => "1e-13",
            "SPACEAGORA_TELEMETRY_ABSTOL_ATM" => "1e-13"
        )
    end
    return Dict(
        "SPACEAGORA_TELEMETRY_SOLVER_MODE" => "dp8",
        "SPACEAGORA_TELEMETRY_DT_MAX_ORBIT" => "10.0",
        "SPACEAGORA_TELEMETRY_RELTOL_ORBIT" => "1e-12",
        "SPACEAGORA_TELEMETRY_ABSTOL_ORBIT" => "1e-12",
        "SPACEAGORA_TELEMETRY_RELTOL_ATM" => "1e-12",
        "SPACEAGORA_TELEMETRY_ABSTOL_ATM" => "1e-12"
    )
end

@inline function _gmat_example_file_to_scenario_name(file_name::String)::String
    stem = replace(file_name, r"\.feather$" => "")
    stem = replace(stem, r"^Sim_" => "")
    stem = replace(stem, r"_1M_" => "_")
    stem = replace(lowercase(stem), " " => "_")
    # GMAT example assets use "Luna" in filenames, while the internal scenario
    # naming throughout this matrix remains "moon_*".
    stem = replace(stem, r"^luna_" => "moon_")
    return stem
end

function _initial_time_dict(initial_time::SM.InitialTime)::Dict{String, Any}
    return Dict{String, Any}(
        "year" => Int(initial_time.year),
        "month" => Int(initial_time.month),
        "day" => Int(initial_time.day),
        "hour" => Int(initial_time.hour),
        "minute" => Int(initial_time.minute),
        "second" => Float64(initial_time.second)
    )
end

function _initial_time_from_unix(unix_seconds::Float64)::SM.InitialTime
    dt = unix2datetime(unix_seconds)
    return SM.InitialTime(
        year=Int32(year(dt)),
        month=Int16(month(dt)),
        day=Int16(day(dt)),
        hour=Int16(hour(dt)),
        minute=Int16(minute(dt)),
        second=Float32(second(dt) + millisecond(dt) / 1.0e3)
    )
end

function _load_cygnss_96hr_j2000_series()
    cached = _CYGNSS_96HR_J2000_CACHE[]
    cached !== nothing && return cached

    # Ensure the leapseconds and Earth orientation kernels are in the pool
    # before calling utc2et/sxform on the 96hr ECEF telemetry.
    TV._planet_from_name("earth")

    df = DataFrame(Arrow.Table(_CYGNSS_96HR_TELEMETRY_FEATHER))
    @test nrow(df) > 0

    t_unix = Float64.(df[!, "pvt_unix_seconds"])
    x_ecef_km = Float64.(df[!, "sc_pos_x_pvt_m"]) .* 1.0e-3
    y_ecef_km = Float64.(df[!, "sc_pos_y_pvt_m"]) .* 1.0e-3
    z_ecef_km = Float64.(df[!, "sc_pos_z_pvt_m"]) .* 1.0e-3
    vx_ecef_kmps = Float64.(df[!, "sc_vel_x_pvt_mps"]) .* 1.0e-3
    vy_ecef_kmps = Float64.(df[!, "sc_vel_y_pvt_mps"]) .* 1.0e-3
    vz_ecef_kmps = Float64.(df[!, "sc_vel_z_pvt_mps"]) .* 1.0e-3

    perm = sortperm(t_unix)
    t_unix = t_unix[perm]
    x_ecef_km = x_ecef_km[perm]
    y_ecef_km = y_ecef_km[perm]
    z_ecef_km = z_ecef_km[perm]
    vx_ecef_kmps = vx_ecef_kmps[perm]
    vy_ecef_kmps = vy_ecef_kmps[perm]
    vz_ecef_kmps = vz_ecef_kmps[perm]

    t_rel = t_unix .- t_unix[1]
    initial_time = _initial_time_from_unix(t_unix[1])
    et0 = TV._initial_time_et(initial_time)

    n = length(t_rel)
    x_j2000_km = Vector{Float64}(undef, n)
    y_j2000_km = Vector{Float64}(undef, n)
    z_j2000_km = Vector{Float64}(undef, n)
    vx_j2000_kmps = Vector{Float64}(undef, n)
    vy_j2000_kmps = Vector{Float64}(undef, n)
    vz_j2000_kmps = Vector{Float64}(undef, n)

    @inbounds for i in eachindex(t_rel)
        r_j2000_m, v_j2000_mps = TV._planet_fixed_to_j2000_state(
            "earth",
            et0 + t_rel[i],
            SVector{3, Float64}(x_ecef_km[i], y_ecef_km[i], z_ecef_km[i]) .* 1.0e3,
            SVector{3, Float64}(vx_ecef_kmps[i], vy_ecef_kmps[i], vz_ecef_kmps[i]) .* 1.0e3
        )
        x_j2000_km[i] = r_j2000_m[1] * 1.0e-3
        y_j2000_km[i] = r_j2000_m[2] * 1.0e-3
        z_j2000_km[i] = r_j2000_m[3] * 1.0e-3
        vx_j2000_kmps[i] = v_j2000_mps[1] * 1.0e-3
        vy_j2000_kmps[i] = v_j2000_mps[2] * 1.0e-3
        vz_j2000_kmps[i] = v_j2000_mps[3] * 1.0e-3
    end

    series = (
        initial_time=initial_time,
        t_rel=t_rel,
        x_km=x_j2000_km,
        y_km=y_j2000_km,
        z_km=z_j2000_km,
        vx_kmps=vx_j2000_kmps,
        vy_kmps=vy_j2000_kmps,
        vz_kmps=vz_j2000_kmps
    )
    _CYGNSS_96HR_J2000_CACHE[] = series
    return series
end

function _load_cygnss_cyg04_96hr_inertial_series()
    cached = _CYGNSS_CYG04_96HR_INERTIAL_CACHE[]
    cached !== nothing && return cached

    TV._planet_from_name("earth")

    df = DataFrame(Arrow.Table(_CYGNSS_CYG04_96HR_TELEMETRY_FEATHER))
    @test nrow(df) > 0

    t_unix = Float64.(df[!, "pvt_unix_seconds"])
    t_native = Float64.(df[!, "time"])
    x_km = Float64.(df[!, "pos_ii_1"]) .* 1.0e-3
    y_km = Float64.(df[!, "pos_ii_2"]) .* 1.0e-3
    z_km = Float64.(df[!, "pos_ii_3"]) .* 1.0e-3
    vx_kmps = Float64.(df[!, "vel_ii_1"]) .* 1.0e-3
    vy_kmps = Float64.(df[!, "vel_ii_2"]) .* 1.0e-3
    vz_kmps = Float64.(df[!, "vel_ii_3"]) .* 1.0e-3

    perm = sortperm(t_unix)
    t_unix = t_unix[perm]
    t_native = t_native[perm]
    x_km = x_km[perm]
    y_km = y_km[perm]
    z_km = z_km[perm]
    vx_kmps = vx_kmps[perm]
    vy_kmps = vy_kmps[perm]
    vz_kmps = vz_kmps[perm]

    series = (
        initial_time=_initial_time_from_unix(t_unix[1]),
        t_rel=t_native .- t_native[1],
        t_abs=t_native,
        x_km=x_km,
        y_km=y_km,
        z_km=z_km,
        vx_kmps=vx_kmps,
        vy_kmps=vy_kmps,
        vz_kmps=vz_kmps
    )
    _CYGNSS_CYG04_96HR_INERTIAL_CACHE[] = series
    return series
end

function _gmat_matrix_expected_scenario_names()::Set{String}
    files = filter(
        f -> startswith(f, "Sim_") && endswith(lowercase(f), ".feather") && !(f in _GMAT_MATRIX_EXCLUDED_FILES),
        readdir(_BASILISK_REFERENCE_DIR)
    )
    return Set(_gmat_example_file_to_scenario_name.(files))
end

function _selected_gmat_scenario_names()::Union{Nothing, Set{String}}
    raw = strip(get(ENV, "SPACEAGORA_GMAT_SCENARIOS", ""))
    isempty(raw) && return nothing
    tokens = Set(filter!(!isempty, strip.(split(raw, ","))))
    return tokens
end

@inline function _gmat_matrix_cache_key()::String
    selected = _selected_gmat_scenario_names()
    return selected === nothing ? "__all__" : join(sort!(collect(selected)), ",")
end

@inline function _filter_gmat_scenarios(scenarios::Vector{Dict{String, Any}})::Vector{Dict{String, Any}}
    selected = _selected_gmat_scenario_names()
    selected === nothing && return scenarios
    return filter(scenario -> String(scenario["name"]) in selected, scenarios)
end

@inline function _active_gmat_expected_scenario_names()::Set{String}
    selected = _selected_gmat_scenario_names()
    return selected === nothing ? _gmat_matrix_expected_scenario_names() : selected
end

@inline function _strict_position_rmse_limit_km(scenario_name::String, profile::Symbol)::Float64
    if startswith(scenario_name, "earth_j50_")
        return profile == :full ? 0.45 : 0.05
    elseif startswith(scenario_name, "earth_j2_")
        return profile == :full ? 1e-2 : 1e-3
    elseif startswith(scenario_name, "earth_j0_")
        return profile == :full ? 1e-2 : 1e-3
    elseif startswith(scenario_name, "mars_j50_")
        return profile == :full ? 0.5 : 0.05
    elseif startswith(scenario_name, "mars_j2_")
        return profile == :full ? 5e-2 : 5e-3
    elseif startswith(scenario_name, "mars_j0_")
        return profile == :full ? 5e-2 : 5e-3
    elseif startswith(scenario_name, "venus_j50_")
        return profile == :full ? 0.5 : 0.05
    elseif startswith(scenario_name, "venus_j2_")
        return profile == :full ? 5e-2 : 5e-3
    elseif startswith(scenario_name, "venus_j0_")
        return profile == :full ? 5e-2 : 5e-3
    elseif startswith(scenario_name, "moon_j50_")
        return profile == :full ? 0.5 : 0.05
    elseif startswith(scenario_name, "moon_j2_")
        return profile == :full ? 5e-2 : 5e-3
    elseif startswith(scenario_name, "moon_j0_")
        return profile == :full ? 5e-2 : 5e-3
    end
    return 1e-3
end

function _required_column(df::DataFrame, candidates::Vector{String})::String
    for cname in candidates
        if cname in names(df)
            return cname
        end
    end
    throw(ArgumentError("None of the candidate columns were found: $(join(candidates, ", "))"))
end

@inline function _find_column(df::DataFrame, candidates::Vector{String})::Union{Nothing, String}
    for cname in candidates
        if cname in names(df)
            return cname
        end
    end
    return nothing
end

@inline function _read_tabular(path::String)::DataFrame
    if endswith(lowercase(path), ".feather")
        return DataFrame(Arrow.Table(path))
    end
    return CSV.read(path, DataFrame)
end

function _build_time_aligned_reference(
    csv_path::String,
    planet_name::String,
    outdir::String,
    stem::String;
    sma_km::Float64=7000.0,
    ecc::Float64=1.0e-5,
    inc_deg::Float64=28.5,
    aop_deg::Float64=0.0,
    raan_deg::Float64=45.0,
    ta_deg::Float64=0.0,
    raan_offset_deg::Float64=0.0,
    ta_offset_deg::Float64=0.0
)
    df = _read_tabular(csv_path)

    time_col = _required_column(df, ["ElapsedSecs", "Sat.ElapsedSecs", "DefaultSC.ElapsedSecs"])
    x_col = _required_column(df, ["X", "Sat.X", "Sat.PlanetInertial.X", "DefaultSC.EarthMJ2000Eq.X"])
    y_col = _required_column(df, ["Y", "Sat.Y", "Sat.PlanetInertial.Y", "DefaultSC.EarthMJ2000Eq.Y"])
    z_col = _required_column(df, ["Z", "Sat.Z", "Sat.PlanetInertial.Z", "DefaultSC.EarthMJ2000Eq.Z"])
    vx_col = _required_column(df, ["VX", "Sat.VX", "Sat.PlanetInertial.VX", "DefaultSC.EarthMJ2000Eq.VX"])
    vy_col = _required_column(df, ["VY", "Sat.VY", "Sat.PlanetInertial.VY", "DefaultSC.EarthMJ2000Eq.VY"])
    vz_col = _required_column(df, ["VZ", "Sat.VZ", "Sat.PlanetInertial.VZ", "DefaultSC.EarthMJ2000Eq.VZ"])
    planet_fixed_x_col = _find_column(df, ["PlanetFixedX", "Sat.PlanetFixed.X"])
    planet_fixed_y_col = _find_column(df, ["PlanetFixedY", "Sat.PlanetFixed.Y"])
    planet_fixed_z_col = _find_column(df, ["PlanetFixedZ", "Sat.PlanetFixed.Z"])

    t = Float64.(df[!, time_col])
    x = Float64.(df[!, x_col])
    y = Float64.(df[!, y_col])
    z = Float64.(df[!, z_col])
    vx = Float64.(df[!, vx_col])
    vy = Float64.(df[!, vy_col])
    vz = Float64.(df[!, vz_col])
    planet_fixed_x = planet_fixed_x_col === nothing ? similar(x) : Float64.(df[!, planet_fixed_x_col])
    planet_fixed_y = planet_fixed_y_col === nothing ? similar(y) : Float64.(df[!, planet_fixed_y_col])
    planet_fixed_z = planet_fixed_z_col === nothing ? similar(z) : Float64.(df[!, planet_fixed_z_col])
    planet = TV._planet_from_name(planet_name)

    perm = sortperm(t)
    t = t[perm]
    x = x[perm]
    y = y[perm]
    z = z[perm]
    vx = vx[perm]
    vy = vy[perm]
    vz = vz[perm]
    planet_fixed_x = planet_fixed_x[perm]
    planet_fixed_y = planet_fixed_y[perm]
    planet_fixed_z = planet_fixed_z[perm]

    t_rel = t .- t[1]

    # GMAT exports in this folder can mix km and m positions across scenarios.
    # Detect oversized radii relative to the planet and convert m -> km when needed.
    rp_km = planet.Rp_e * 1e-3
    r_raw = sqrt.(x .^ 2 .+ y .^ 2 .+ z .^ 2)
    if median(r_raw) > 100.0 * rp_km
        x .*= 1.0e-3
        y .*= 1.0e-3
        z .*= 1.0e-3
        vx .*= 1.0e-3
        vy .*= 1.0e-3
        vz .*= 1.0e-3
        if planet_fixed_x_col !== nothing
            planet_fixed_x .*= 1.0e-3
            planet_fixed_y .*= 1.0e-3
            planet_fixed_z .*= 1.0e-3
        end
    end

    if planet_fixed_x_col === nothing || planet_fixed_y_col === nothing || planet_fixed_z_col === nothing
        ephemerides_model = SM.SpiceEphemeridesModel()
        et0 = SM.ephemerides_time_seconds(_default_matrix_initial_time(), ephemerides_model)
        for i in eachindex(t_rel)
            lpi = SM.planet_frame_lpi(planet, et0 + t_rel[i], ephemerides_model)
            r_pf = lpi * SVector{3, Float64}(x[i], y[i], z[i])
            planet_fixed_x[i] = r_pf[1]
            planet_fixed_y[i] = r_pf[2]
            planet_fixed_z[i] = r_pf[3]
        end
    end

    r_km = sqrt.(x .^ 2 .+ y .^ 2 .+ z .^ 2)
    alt_km = r_km .- planet.Rp_e * 1e-3

    raan_deg = _mod360(raan_deg + raan_offset_deg)
    ta_deg = _mod360(ta_deg + ta_offset_deg)

    telemetry_df = DataFrame(
        time_s=t_rel,
        altitude_km=alt_km,
        x_km=x,
        y_km=y,
        z_km=z,
        planet_fixed_x_km=planet_fixed_x,
        planet_fixed_y_km=planet_fixed_y,
        planet_fixed_z_km=planet_fixed_z,
        sma_km=fill(sma_km, length(t_rel)),
        ecc=fill(ecc, length(t_rel)),
        inc_deg=fill(inc_deg, length(t_rel)),
        aop_deg=fill(aop_deg, length(t_rel)),
        raan_deg=fill(raan_deg, length(t_rel)),
        ta_deg=fill(ta_deg, length(t_rel)),
        x_ic_km=fill(x[1], length(t_rel)),
        y_ic_km=fill(y[1], length(t_rel)),
        z_ic_km=fill(z[1], length(t_rel)),
        vx_ic_kmps=fill(vx[1], length(t_rel)),
        vy_ic_kmps=fill(vy[1], length(t_rel)),
        vz_ic_kmps=fill(vz[1], length(t_rel))
    )

    telemetry_path = joinpath(outdir, "$(stem)_time_aligned.arrow")
    Arrow.write(telemetry_path, telemetry_df)

    return (
        telemetry_path=telemetry_path,
        x_ic_km=x[1],
        y_ic_km=y[1],
        z_ic_km=z[1],
        vx_ic_kmps=vx[1],
        vy_ic_kmps=vy[1],
        vz_ic_kmps=vz[1],
        sma_km=sma_km,
        ecc=ecc,
        inc_deg=inc_deg,
        aop_deg=aop_deg,
        raan_deg=raan_deg,
        ta_deg=ta_deg
    )
end

function _build_cygnss_48hr_reference(outdir::String, stem::String)
    df = DataFrame(Arrow.Table(_CYGNSS_48HR_TELEMETRY_FEATHER))
    @test nrow(df) > 0

    time_col = _required_column(df, ["TIME OFFSET", "time"])
    x_col = _required_column(df, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_X (m)", "pos_ii_1"])
    y_col = _required_column(df, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_Y (m)", "pos_ii_2"])
    z_col = _required_column(df, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_Z (m)", "pos_ii_3"])
    vx_col = _required_column(df, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_X (m/s)", "vel_ii_1"])
    vy_col = _required_column(df, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_Y (m/s)", "vel_ii_2"])
    vz_col = _required_column(df, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_Z (m/s)", "vel_ii_3"])

    t = Float64.(df[!, time_col])
    x_km = Float64.(df[!, x_col]) .* 1.0e-3
    y_km = Float64.(df[!, y_col]) .* 1.0e-3
    z_km = Float64.(df[!, z_col]) .* 1.0e-3
    vx_kmps = Float64.(df[!, vx_col]) .* 1.0e-3
    vy_kmps = Float64.(df[!, vy_col]) .* 1.0e-3
    vz_kmps = Float64.(df[!, vz_col]) .* 1.0e-3

    perm = sortperm(t)
    t = t[perm]
    x_km = x_km[perm]
    y_km = y_km[perm]
    z_km = z_km[perm]
    vx_kmps = vx_kmps[perm]
    vy_kmps = vy_kmps[perm]
    vz_kmps = vz_kmps[perm]

    t_rel = t .- t[1]
    planet = TV._planet_from_name("earth")
    # Altitude in the same (oblate/geodetic) convention the scenario scoring
    # uses; the vacuum r - Rp_e form carries a latitude-dependent bias of up
    # to ~21 km against the oblate channel.
    alt_km = [
        TV._telemetry_altitude_km(
            SVector{3, Float64}(x_km[i], y_km[i], z_km[i]) .* 1.0e3,
            SVector{3, Float64}(vx_kmps[i], vy_kmps[i], vz_kmps[i]) .* 1.0e3,
            planet,
            :oblate
        ) for i in eachindex(t_rel)
    ]

    # Published CYGNSS IC recipe: average the vis-viva SMA over the first 5
    # samples and scale the first-sample velocity to that orbital energy. The
    # single-sample GPS state carries ~63 m of SMA noise (~3.5 cm/s velocity
    # noise at 1786 m SMA per m/s), which drifts ~18 km along-track over 48 h
    # and quadruples the entry-point RMSE (6.2 km vs the published 1.578 km).
    mu_kmc = planet.μ * 1.0e-9
    a_samples_km = [
        let r = sqrt(x_km[i]^2 + y_km[i]^2 + z_km[i]^2),
            v = sqrt(vx_kmps[i]^2 + vy_kmps[i]^2 + vz_kmps[i]^2)
            1.0 / (2.0 / r - v^2 / mu_kmc)
        end for i in 1:5
    ]
    a_target_km = sum(a_samples_km) / length(a_samples_km)
    r0_km = sqrt(x_km[1]^2 + y_km[1]^2 + z_km[1]^2)
    v_target_kmps = sqrt(mu_kmc * (2.0 / r0_km - 1.0 / a_target_km))
    v_scale = v_target_kmps / sqrt(vx_kmps[1]^2 + vy_kmps[1]^2 + vz_kmps[1]^2)
    vx_ic = vx_kmps[1] * v_scale
    vy_ic = vy_kmps[1] * v_scale
    vz_ic = vz_kmps[1] * v_scale

    oe0 = TV.rvtoorbitalelement(
        SVector{3, Float64}(x_km[1], y_km[1], z_km[1]) .* 1.0e3,
        SVector{3, Float64}(vx_ic, vy_ic, vz_ic) .* 1.0e3,
        planet
    )
    sma_km = oe0[1] * 1.0e-3
    ecc = oe0[2]
    inc_deg = rad2deg(oe0[3])
    raan_deg = rad2deg(oe0[4])
    aop_deg = rad2deg(oe0[5])
    ta_deg = rad2deg(oe0[6])

    telemetry_df = DataFrame(
        time_s=t_rel,
        altitude_km=alt_km,
        x_km=x_km,
        y_km=y_km,
        z_km=z_km,
        sma_km=fill(sma_km, length(t_rel)),
        ecc=fill(ecc, length(t_rel)),
        inc_deg=fill(inc_deg, length(t_rel)),
        aop_deg=fill(aop_deg, length(t_rel)),
        raan_deg=fill(raan_deg, length(t_rel)),
        ta_deg=fill(ta_deg, length(t_rel)),
        x_ic_km=fill(x_km[1], length(t_rel)),
        y_ic_km=fill(y_km[1], length(t_rel)),
        z_ic_km=fill(z_km[1], length(t_rel)),
        vx_ic_kmps=fill(vx_ic, length(t_rel)),
        vy_ic_kmps=fill(vy_ic, length(t_rel)),
        vz_ic_kmps=fill(vz_ic, length(t_rel))
    )

    telemetry_path = joinpath(outdir, "$(stem)_time_aligned.arrow")
    Arrow.write(telemetry_path, telemetry_df)

    return (
        telemetry_path=telemetry_path,
        t0_s=t_rel[1],
        tf_s=t_rel[end],
        x_ic_km=x_km[1],
        y_ic_km=y_km[1],
        z_ic_km=z_km[1],
        vx_ic_kmps=vx_ic,
        vy_ic_kmps=vy_ic,
        vz_ic_kmps=vz_ic,
        sma_km=sma_km,
        ecc=ecc,
        inc_deg=inc_deg,
        aop_deg=aop_deg,
        raan_deg=raan_deg,
        ta_deg=ta_deg
    )
end

function _build_cygnss_96hr_reference(outdir::String, stem::String)
    series = _load_cygnss_96hr_j2000_series()
    t_rel = series.t_rel
    x_km = series.x_km
    y_km = series.y_km
    z_km = series.z_km
    vx_kmps = series.vx_kmps
    vy_kmps = series.vy_kmps
    vz_kmps = series.vz_kmps
    planet = TV._planet_from_name("earth")
    # Same oblate/geodetic altitude convention as the scenario scoring channel.
    alt_km = [
        TV._telemetry_altitude_km(
            SVector{3, Float64}(x_km[i], y_km[i], z_km[i]) .* 1.0e3,
            SVector{3, Float64}(vx_kmps[i], vy_kmps[i], vz_kmps[i]) .* 1.0e3,
            planet,
            :oblate
        ) for i in eachindex(x_km)
    ]

    oe0 = TV.rvtoorbitalelement(
        SVector{3, Float64}(x_km[1], y_km[1], z_km[1]) .* 1.0e3,
        SVector{3, Float64}(vx_kmps[1], vy_kmps[1], vz_kmps[1]) .* 1.0e3,
        planet
    )
    sma_km  = oe0[1] * 1.0e-3
    ecc     = oe0[2]
    inc_deg  = rad2deg(oe0[3])
    raan_deg = rad2deg(oe0[4])
    aop_deg  = rad2deg(oe0[5])
    ta_deg   = rad2deg(oe0[6])

    telemetry_df = DataFrame(
        time_s=t_rel,
        altitude_km=alt_km,
        x_km=x_km,
        y_km=y_km,
        z_km=z_km,
        sma_km=fill(sma_km, length(t_rel)),
        ecc=fill(ecc, length(t_rel)),
        inc_deg=fill(inc_deg, length(t_rel)),
        aop_deg=fill(aop_deg, length(t_rel)),
        raan_deg=fill(raan_deg, length(t_rel)),
        ta_deg=fill(ta_deg, length(t_rel)),
        x_ic_km=fill(x_km[1], length(t_rel)),
        y_ic_km=fill(y_km[1], length(t_rel)),
        z_ic_km=fill(z_km[1], length(t_rel)),
        vx_ic_kmps=fill(vx_kmps[1], length(t_rel)),
        vy_ic_kmps=fill(vy_kmps[1], length(t_rel)),
        vz_ic_kmps=fill(vz_kmps[1], length(t_rel))
    )

    telemetry_path = joinpath(outdir, "$(stem)_time_aligned.arrow")
    Arrow.write(telemetry_path, telemetry_df)

    return (
        telemetry_path=telemetry_path,
        t0_s=t_rel[1],
        tf_s=t_rel[end],
        initial_time=series.initial_time,
        x_ic_km=x_km[1],
        y_ic_km=y_km[1],
        z_ic_km=z_km[1],
        vx_ic_kmps=vx_kmps[1],
        vy_ic_kmps=vy_kmps[1],
        vz_ic_kmps=vz_kmps[1],
        sma_km=sma_km,
        ecc=ecc,
        inc_deg=inc_deg,
        aop_deg=aop_deg,
        raan_deg=raan_deg,
        ta_deg=ta_deg
    )
end

function _build_cygnss_cyg04_96hr_inertial_reference(outdir::String, stem::String)
    series = _load_cygnss_cyg04_96hr_inertial_series()
    # t_abs: seconds since 2025-06-06T00:00:00 (same epoch as the 48hr IC).
    # Using t_abs as the time axis aligns the comparison with the propagation
    # started from the 48hr IC at t=0.
    t_abs = series.t_abs
    x_km = series.x_km
    y_km = series.y_km
    z_km = series.z_km

    planet = TV._planet_from_name("earth")
    # Same oblate/geodetic altitude convention as the scenario scoring channel.
    # The planet-fixed rotation of position does not depend on velocity, so a
    # zero velocity is passed for this position-only product.
    alt_km = [
        TV._telemetry_altitude_km(
            SVector{3, Float64}(x_km[i], y_km[i], z_km[i]) .* 1.0e3,
            SVector{3, Float64}(0.0, 0.0, 0.0),
            planet,
            :oblate
        ) for i in eachindex(x_km)
    ]

    # Use the 48hr IC: the cygnss_data_48hr.feather velocity is orbit-determined
    # and gives a much more accurate initial orbital energy than the raw GPS
    # Doppler velocity in the CYG04 file.  A ~0.6 m/s along-track error in the
    # CYG04 raw velocity shifts the SMA by ~1 km, producing ~120 km along-track
    # drift at 48 hours even though the reference positions agree to < 2 m.
    #
    # NOTE: the matching single-sample-noise fix applied to _build_cygnss_48hr_reference
    # above (average vis-viva SMA over a short window, rescale velocity magnitude to
    # match) was ALSO tried here and made this scenario noticeably worse (RMSE 3.70km ->
    # 15.44km against the CYG04 96hr reference), even though it substantially improved
    # the 48hr_pvt scenario against ITS OWN reference. The two scenarios compare against
    # different, independently-sourced reference trajectories with evidently different
    # noise/bias characteristics, so a correction tuned against one is not safe to reuse
    # against the other -- reverted here; kept only where it was actually validated.
    df48 = DataFrame(Arrow.Table(_CYGNSS_48HR_TELEMETRY_FEATHER))
    sort!(df48, _required_column(df48, ["TIME OFFSET", "time"]))
    x_ic_km   = Float64(df48[!, _required_column(df48, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_X (m)", "pos_ii_1"])][1]) * 1.0e-3
    y_ic_km   = Float64(df48[!, _required_column(df48, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_Y (m)", "pos_ii_2"])][1]) * 1.0e-3
    z_ic_km   = Float64(df48[!, _required_column(df48, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_Z (m)", "pos_ii_3"])][1]) * 1.0e-3
    vx_ic_kmps = Float64(df48[!, _required_column(df48, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_X (m/s)", "vel_ii_1"])][1]) * 1.0e-3
    vy_ic_kmps = Float64(df48[!, _required_column(df48, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_Y (m/s)", "vel_ii_2"])][1]) * 1.0e-3
    vz_ic_kmps = Float64(df48[!, _required_column(df48, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_Z (m/s)", "vel_ii_3"])][1]) * 1.0e-3

    oe0 = TV.rvtoorbitalelement(
        SVector{3, Float64}(x_ic_km, y_ic_km, z_ic_km) .* 1.0e3,
        SVector{3, Float64}(vx_ic_kmps, vy_ic_kmps, vz_ic_kmps) .* 1.0e3,
        planet
    )
    sma_km = oe0[1] * 1.0e-3
    ecc = oe0[2]
    inc_deg = rad2deg(oe0[3])
    raan_deg = rad2deg(oe0[4])
    aop_deg = rad2deg(oe0[5])
    ta_deg = rad2deg(oe0[6])

    telemetry_df = DataFrame(
        time_s=t_abs,
        altitude_km=alt_km,
        x_km=x_km,
        y_km=y_km,
        z_km=z_km,
        sma_km=fill(sma_km, length(t_abs)),
        ecc=fill(ecc, length(t_abs)),
        inc_deg=fill(inc_deg, length(t_abs)),
        aop_deg=fill(aop_deg, length(t_abs)),
        raan_deg=fill(raan_deg, length(t_abs)),
        ta_deg=fill(ta_deg, length(t_abs)),
        x_ic_km=fill(x_ic_km, length(t_abs)),
        y_ic_km=fill(y_ic_km, length(t_abs)),
        z_ic_km=fill(z_ic_km, length(t_abs)),
        vx_ic_kmps=fill(vx_ic_kmps, length(t_abs)),
        vy_ic_kmps=fill(vy_ic_kmps, length(t_abs)),
        vz_ic_kmps=fill(vz_ic_kmps, length(t_abs))
    )

    telemetry_path = joinpath(outdir, "$(stem)_time_aligned.arrow")
    Arrow.write(telemetry_path, telemetry_df)

    return (
        telemetry_path=telemetry_path,
        t0_s=t_abs[1],
        tf_s=t_abs[end],
        x_ic_km=x_ic_km,
        y_ic_km=y_ic_km,
        z_ic_km=z_ic_km,
        vx_ic_kmps=vx_ic_kmps,
        vy_ic_kmps=vy_ic_kmps,
        vz_ic_kmps=vz_ic_kmps,
        sma_km=sma_km,
        ecc=ecc,
        inc_deg=inc_deg,
        aop_deg=aop_deg,
        raan_deg=raan_deg,
        ta_deg=ta_deg
    )
end

@inline function _mod360(x::Float64)::Float64
    y = mod(x, 360.0)
    return y < 0.0 ? y + 360.0 : y
end

function _base_scenario_dict(name::String, telemetry_path::String)
    return Dict{String, Any}(
        "name" => name,
        "kind" => "time_aligned_state",
        "events" => Any["altitude_time", "state_x_time", "state_y_time", "state_z_time"],
        "telemetry" => telemetry_path,
        "telemetry_columns" => Dict{String, Any}(
            "time" => "time_s",
            "altitude" => "altitude_km",
            "x" => "x_km",
            "y" => "y_km",
            "z" => "z_km",
            "x_ic" => "x_ic_km",
            "y_ic" => "y_ic_km",
            "z_ic" => "z_ic_km",
            "vx_ic" => "vx_ic_kmps",
            "vy_ic" => "vy_ic_kmps",
            "vz_ic" => "vz_ic_kmps"
        ),
        "max_points_quick" => 10000,
        "max_points_full" => 100000,
        "min_eval_points" => 2,
        "units" => Dict{String, Any}(
            "x" => "s",
            "altitude_time" => "km",
            "state_x_time" => "km",
            "state_y_time" => "km",
            "state_z_time" => "km"
        ),
        "tolerances_quick" => Dict{String, Any}(
            "altitude_time" => Dict("max_abs_km" => 1.0e6, "max_nmae" => 1.0e6, "max_rmse_km" => 1.0e6),
            "state_x_time" => Dict("max_abs_km" => 1.0e6, "max_nmae" => 1.0e6, "max_rmse_km" => 1.0e6),
            "state_y_time" => Dict("max_abs_km" => 1.0e6, "max_nmae" => 1.0e6, "max_rmse_km" => 1.0e6),
            "state_z_time" => Dict("max_abs_km" => 1.0e6, "max_nmae" => 1.0e6, "max_rmse_km" => 1.0e6)
        ),
        "tolerances_full" => Dict{String, Any}(
            "altitude_time" => Dict("max_abs_km" => 1.0e6, "max_nmae" => 1.0e6, "max_rmse_km" => 1.0e6),
            "state_x_time" => Dict("max_abs_km" => 1.0e6, "max_nmae" => 1.0e6, "max_rmse_km" => 1.0e6),
            "state_y_time" => Dict("max_abs_km" => 1.0e6, "max_nmae" => 1.0e6, "max_rmse_km" => 1.0e6),
            "state_z_time" => Dict("max_abs_km" => 1.0e6, "max_nmae" => 1.0e6, "max_rmse_km" => 1.0e6)
        ),
        "initial_time" => Dict{String, Any}(
            "year" => 2026,
            "month" => 1,
            "day" => 1,
            "hour" => 12,
            "minute" => 0,
            "second" => 0.0
        ),
        "spacecraft" => Dict{String, Any}(
            "bus_dims_m" => Any[2.05e-1, 3.7e-1, 0.8e-1],
            "panel_dims_m" => Any[10e-3, 28.5e-3, 0.0001],
            "bus_mass_kg" => 29.0,
            "panel_mass_each_kg" => 0.0,
            "panel_offset_y_m" => 2.45,
            "prop_mass_kg" => 0.0,
            "id" => 1002
        ),
        "atmosphere_truth" => Dict{String, Any}(
            "assumption_id" => "earth_gmat_gram_deterministic_v1",
            "atmosphere_model" => "GRAM",
            "atmosphere_dataset" => "MERRA2 All Mean",
            "space_weather_model" => "EarthGRAM MERRA2 climatology (deterministic)",
            "solar_flux_model" => "EarthGRAM/MERRA2 epoch-fixed defaults",
            "gram_seed" => 1001,
            "gram_perturbation_scales" => Any[0.0, 0.0, 0.0, 0.0],
            "gram_offline_surrogate" => "auto",
            "gram_static_grid" => false,
            "gram_track_cache" => false,
            "gram_global_lock" => "on"
        ),
        "calibration" => Dict{String, Any}("enabled" => false),
        "drag_enabled" => false,
        "EI_km" => 300.0,
        "comparison_mode" => "time_aligned_state"
    )
end

function _scenario_rmse(summary::DataFrame, name::String)::Float64
    rows = summary[(summary.scenario .== name) .& in.(summary.event, Ref(["state_x_time", "state_y_time", "state_z_time"])), :]
    @test nrow(rows) == 3
    rmse = mean(Float64.(rows.rmse_km))
    @test isfinite(rmse)
    return rmse
end

@inline function _scenario_planet_name(scenario_name::String)::String
    stem = first(split(scenario_name, "_"; limit=2))
    return uppercasefirst(stem)
end

@inline function _scenario_basilisk_file_name(scenario_name::String)::String
    parts = split(scenario_name, "_")
    @test length(parts) == 3
    body = parts[1] == "moon" ? "Luna" : uppercasefirst(parts[1])
    jtag = uppercase(parts[2])
    tbtag = parts[3] == "tbtrue" ? "TBTrue" : "TBFalse"
    basename = "Sim_$(body)_$(jtag)_$(tbtag).feather"
    if isfile(joinpath(_BASILISK_REFERENCE_DIR, basename))
        return basename
    end
    full_basename = "Sim_$(body)_1M_$(jtag)_$(tbtag).feather"
    if isfile(joinpath(_BASILISK_REFERENCE_DIR, full_basename))
        return full_basename
    end
    return basename
end

@inline function _scenario_basilisk_path(scenario_name::String)::String
    return joinpath(_BASILISK_REFERENCE_DIR, _scenario_basilisk_file_name(scenario_name))
end

@inline function _scenario_stk_file_name(scenario_name::String)::String
    parts = split(scenario_name, "_")
    @test length(parts) == 3
    body = parts[1] == "moon" ? "Luna" : uppercasefirst(parts[1])
    jtag = uppercase(parts[2])
    tbtag = parts[3] == "tbtrue" ? "TB1" : "TB0"
    return "Sim_$(body)_$(jtag)_$(tbtag).csv"
end

@inline function _scenario_stk_path(scenario_name::String)::String
    return joinpath(_STK_RESULTS_DIR, _scenario_stk_file_name(scenario_name))
end

@inline function _default_matrix_initial_time()::SM.InitialTime
    return SM.InitialTime(
        year=2026,
        month=1,
        day=1,
        hour=12,
        minute=0,
        second=0.0
    )
end

function _matrix_initial_conditions(scenario_name::String)
    planet = first(split(scenario_name, "_"; limit=2))
    if planet == "earth"
        return (
            sma_km=7000.0,
            ecc=1.0e-5,
            inc_deg=28.5,
            aop_deg=0.0,
            raan_deg=45.0,
            ta_deg=0.0
        )
    elseif planet == "mars"
        return (
            sma_km=3999.999999970649,
            ecc=0.001999999999785415,
            inc_deg=45.00000000009245,
            aop_deg=360.0,
            raan_deg=360.0,
            ta_deg=90.00000010508394
        )
    elseif planet == "venus"
        return (
            sma_km=7000.0,
            ecc=0.001,
            inc_deg=30.0,
            aop_deg=0.0,
            raan_deg=0.0,
            ta_deg=180.0
        )
    elseif planet == "moon"
        return (
            sma_km=5000.0,
            ecc=0.6,
            inc_deg=42.0,
            aop_deg=0.0,
            raan_deg=0.0,
            ta_deg=0.0
        )
    end
    error("Unsupported matrix scenario planet in $scenario_name")
end

function _matrix_scenario_overrides(scenario_name::String)::Dict{String, Any}
    parts = split(scenario_name, "_")
    @test length(parts) == 3
    planet, gravity_tag, tb_tag = parts

    gravity_degree, gravity_order = if gravity_tag == "j0"
        (0, 0)
    elseif gravity_tag == "j2"
        planet == "earth" ? (2, 0) : (2, 2)
    elseif gravity_tag == "j50"
        (50, 50)
    else
        error("Unsupported gravity tag in $scenario_name")
    end

    nbody_bodies = if tb_tag == "tbtrue"
        if planet == "earth"
            Any["sun", "moon"]
        elseif planet in ("mars", "venus")
            Any["sun"]
        elseif planet == "moon"
            Any["sun", "earth"]
        else
            Any[]
        end
    else
        Any[]
    end

    harmonics_file = if planet == "earth"
        _GMAT_HARMONICS_EARTH_FILE
    elseif planet == "mars"
        _GMAT_HARMONICS_MARS_FILE
    elseif planet == "venus"
        _GMAT_HARMONICS_VENUS_FILE
    elseif planet == "moon"
        _GMAT_HARMONICS_MOON_FILE
    else
        error("Unsupported planet in $scenario_name")
    end

    overrides = Dict{String, Any}(
        "planet" => planet,
        "gravity_model" => "inverse_squared",
        "gravity_harmonics_degree" => gravity_degree,
        "gravity_harmonics_order" => gravity_order,
        "gravity_harmonics_file" => harmonics_file,
        "nbody_bodies" => nbody_bodies,
        "orbit_altitude_mode" => "oblate"
    )

    if scenario_name == "earth_j0_tbtrue"
        merge!(overrides, Dict{String, Any}(
            "srp_enabled" => false,
            "srp_area_m2" => 5.0,
            "srp_cr" => 1.35
        ))
    end

    if planet == "moon"
        merge!(overrides, Dict{String, Any}(
            "drag_enabled" => false,
            "include_wind" => false
        ))
    end

    return overrides
end

function _scenario_planet_fixed_position_rmse(errors::DataFrame, scenario_name::String)
    rows = errors[errors.scenario .== scenario_name, :]
    xrows = rows[rows.event .== "state_x_time", :]
    yrows = rows[rows.event .== "state_y_time", :]
    zrows = rows[rows.event .== "state_z_time", :]
    @test nrow(xrows) >= 2
    @test nrow(yrows) == nrow(xrows)
    @test nrow(zrows) == nrow(xrows)

    t = Float64.(xrows.telemetry_axis)
    @test t == Float64.(yrows.telemetry_axis)
    @test t == Float64.(zrows.telemetry_axis)

    n = min(nrow(xrows), nrow(yrows), nrow(zrows))
    @test n >= 2

    x_sim = Float64.(xrows.sim_interp_value_km[1:n])
    y_sim = Float64.(yrows.sim_interp_value_km[1:n])
    z_sim = Float64.(zrows.sim_interp_value_km[1:n])
    basilisk_path = joinpath(_BASILISK_REFERENCE_DIR, _scenario_basilisk_file_name(scenario_name))
    @test isfile(basilisk_path)
    telemetry_df = _read_tabular(basilisk_path)
    x_tel_pf_col = _required_column(telemetry_df, ["PlanetFixedX", "Sat.PlanetFixed.X"])
    y_tel_pf_col = _required_column(telemetry_df, ["PlanetFixedY", "Sat.PlanetFixed.Y"])
    z_tel_pf_col = _required_column(telemetry_df, ["PlanetFixedZ", "Sat.PlanetFixed.Z"])
    time_tel_col = _required_column(telemetry_df, ["ElapsedSecs", "Sat.ElapsedSecs", "DefaultSC.ElapsedSecs"])
    perm_tel = sortperm(Float64.(telemetry_df[!, time_tel_col]))
    x_tel_pf = Float64.(telemetry_df[perm_tel, x_tel_pf_col])[1:n]
    y_tel_pf = Float64.(telemetry_df[perm_tel, y_tel_pf_col])[1:n]
    z_tel_pf = Float64.(telemetry_df[perm_tel, z_tel_pf_col])[1:n]
    planet = TV._planet_from_name(_scenario_planet_name(scenario_name))
    rp_km = planet.Rp_e * 1e-3
    if median(sqrt.(x_tel_pf .^ 2 .+ y_tel_pf .^ 2 .+ z_tel_pf .^ 2)) > 100.0 * rp_km
        x_tel_pf .*= 1.0e-3
        y_tel_pf .*= 1.0e-3
        z_tel_pf .*= 1.0e-3
    end

    ephemerides_model = SM.SpiceEphemeridesModel()
    et0 = SM.ephemerides_time_seconds(_default_matrix_initial_time(), ephemerides_model)

    pos_sq = Vector{Float64}(undef, n)
    for i in 1:n
        lpi = SM.planet_frame_lpi(planet, et0 + t[i], ephemerides_model)
        r_sim_pf = lpi * SVector{3, Float64}(x_sim[i], y_sim[i], z_sim[i])
        r_tel_pf = SVector{3, Float64}(x_tel_pf[i], y_tel_pf[i], z_tel_pf[i])
        dr_pf = r_sim_pf - r_tel_pf
        pos_sq[i] = dr_pf[1]^2 + dr_pf[2]^2 + dr_pf[3]^2
    end

    first_step_error_km = sqrt(pos_sq[2])
    full_rmse_km = sqrt(mean(pos_sq))
    return (first_step_error_km=first_step_error_km, full_rmse_km=full_rmse_km, n_points=n)
end

const _GMAT_MATRIX_SUMMARY_CACHE = Ref{Union{Nothing, DataFrame}}(nothing)
const _GMAT_MATRIX_RESULT_CACHE = Ref{Union{Nothing, TV.VerificationResult}}(nothing)
const _GMAT_MATRIX_CACHE_KEY = Ref{String}("")
const _STK_MATRIX_SUMMARY_CACHE = Ref{Union{Nothing, DataFrame}}(nothing)
const _STK_MATRIX_RESULT_CACHE = Ref{Union{Nothing, TV.VerificationResult}}(nothing)
const _STK_MATRIX_CACHE_KEY = Ref{String}("")
const _CYGNSS_48HR_SUMMARY_CACHE = Ref{Union{Nothing, DataFrame}}(nothing)
const _CYGNSS_48HR_RESULT_CACHE = Ref{Union{Nothing, TV.VerificationResult}}(nothing)
const _CYGNSS_48HR_TIMESPAN_CACHE = Ref{Union{Nothing, NamedTuple{(:t0_s, :tf_s), Tuple{Float64, Float64}}}}(nothing)
const _CYGNSS_GMAT_SUMMARY_CACHE = Ref{Union{Nothing, DataFrame}}(nothing)
const _CYGNSS_GMAT_RESULT_CACHE = Ref{Union{Nothing, TV.VerificationResult}}(nothing)
const _CYGNSS_GMAT_TIMESPAN_CACHE = Ref{Union{Nothing, NamedTuple{(:t0_s, :tf_s), Tuple{Float64, Float64}}}}(nothing)
const _CYGNSS_96HR_SUMMARY_CACHE = Ref{Union{Nothing, DataFrame}}(nothing)
const _CYGNSS_96HR_RESULT_CACHE = Ref{Union{Nothing, TV.VerificationResult}}(nothing)
const _CYGNSS_96HR_TIMESPAN_CACHE = Ref{Union{Nothing, NamedTuple{(:t0_s, :tf_s), Tuple{Float64, Float64}}}}(nothing)
const _CYGNSS_CYG04_96HR_SUMMARY_CACHE = Ref{Union{Nothing, DataFrame}}(nothing)
const _CYGNSS_CYG04_96HR_RESULT_CACHE = Ref{Union{Nothing, TV.VerificationResult}}(nothing)
const _CYGNSS_CYG04_96HR_TIMESPAN_CACHE = Ref{Union{Nothing, NamedTuple{(:t0_s, :tf_s), Tuple{Float64, Float64}}}}(nothing)

function _run_reference_scenario_matrix_result_once(
    path_resolver::Function,
    result_cache::Base.RefValue{Union{Nothing, TV.VerificationResult}},
    summary_cache::Base.RefValue{Union{Nothing, DataFrame}},
    cache_key_ref::Base.RefValue{String}
)::TV.VerificationResult
    cache_key = _gmat_matrix_cache_key()
    if result_cache[] !== nothing && cache_key_ref[] == cache_key
        return result_cache[]
    end

    active_scenarios = sort!(collect(_active_gmat_expected_scenario_names()))
    for scenario_name in active_scenarios
        @test isfile(path_resolver(scenario_name))
    end
    @test isfile(joinpath(TV.SPICE_PATH, _gmat_planetary_kernel_relpath()))
    @test isfile(joinpath(_GMAT_REPO_ROOT, _GMAT_HARMONICS_EARTH_FILE))
    @test isfile(joinpath(_GMAT_REPO_ROOT, _GMAT_HARMONICS_MARS_FILE))
    @test isfile(joinpath(_GMAT_REPO_ROOT, _GMAT_HARMONICS_VENUS_FILE))
    @test isfile(joinpath(_GMAT_REPO_ROOT, _GMAT_HARMONICS_MOON_FILE))

    result = mktempdir() do tmp
        trajectories = Dict{String, Any}()
        initial_state_logs = String[]
        for scenario_name in active_scenarios
            csv_path = path_resolver(scenario_name)
            ic = _matrix_initial_conditions(scenario_name)
            traj = _build_time_aligned_reference(
                csv_path,
                _scenario_planet_name(scenario_name),
                tmp,
                scenario_name;
                sma_km=ic.sma_km,
                ecc=ic.ecc,
                inc_deg=ic.inc_deg,
                aop_deg=ic.aop_deg,
                raan_deg=ic.raan_deg,
                ta_deg=ic.ta_deg
            )
            trajectories[scenario_name] = traj

            if startswith(scenario_name, "earth_")
                push!(
                    initial_state_logs,
                    "$scenario_name=($(traj.x_ic_km),$(traj.y_ic_km),$(traj.z_ic_km),$(traj.vx_ic_kmps),$(traj.vy_ic_kmps),$(traj.vz_ic_kmps))"
                )
            end
        end

        if !isempty(initial_state_logs)
            println(
                "SpaceAGORA initial Cartesian states from reference [x_km,y_km,z_km,vx_kmps,vy_kmps,vz_kmps]: " *
                join(initial_state_logs, "; ")
            )
        end

        scenarios = Dict{String, Any}[]
        for scenario_name in active_scenarios
            scenario = _base_scenario_dict(scenario_name, trajectories[scenario_name].telemetry_path)
            merge!(scenario, _matrix_scenario_overrides(scenario_name))
            push!(scenarios, scenario)
        end

        manifest_path = joinpath(tmp, "gmat_scenario_matrix_manifest.toml")
        open(manifest_path, "w") do io
            TOML.print(io, Dict{String, Any}("scenarios" => Any[scenarios...]))
        end

        req = TV.VerificationRequest(
            profile=TEST_MODE,
            out_summary=joinpath(tmp, "summary.csv"),
            out_errors=joinpath(tmp, "errors.csv"),
            manifest_path=manifest_path,
            enforce=false,
            generate_plots=false
        )
        solver_env = _telemetry_solver_env_overrides()
        gmat_kernel_env = Pair{String, String}(
            "SPACEAGORA_SPICE_PLANETARY_KERNEL_RELPATH",
            _gmat_planetary_kernel_relpath()
        )
        return withenv(gmat_kernel_env, pairs(solver_env)...) do
            TV.run_verification(req)
        end
    end

    result_cache[] = result
    summary_cache[] = result.summary
    cache_key_ref[] = cache_key
    return result
end

function _run_gmat_scenario_matrix_result_once()::TV.VerificationResult
    return _run_reference_scenario_matrix_result_once(
        _scenario_basilisk_path,
        _GMAT_MATRIX_RESULT_CACHE,
        _GMAT_MATRIX_SUMMARY_CACHE,
        _GMAT_MATRIX_CACHE_KEY
    )
end

function _run_gmat_scenario_matrix_once()::DataFrame
    return _run_gmat_scenario_matrix_result_once().summary
end

function _run_stk_scenario_matrix_result_once()::TV.VerificationResult
    return _run_reference_scenario_matrix_result_once(
        _scenario_stk_path,
        _STK_MATRIX_RESULT_CACHE,
        _STK_MATRIX_SUMMARY_CACHE,
        _STK_MATRIX_CACHE_KEY
    )
end

function _run_stk_scenario_matrix_once()::DataFrame
    return _run_stk_scenario_matrix_result_once().summary
end

function _run_cygnss_48hr_result_once()::TV.VerificationResult
    if _CYGNSS_48HR_RESULT_CACHE[] !== nothing
        return _CYGNSS_48HR_RESULT_CACHE[]
    end

    @test isfile(_CYGNSS_48HR_TELEMETRY_FEATHER)

    result = mktempdir() do tmp
        traj = _build_cygnss_48hr_reference(tmp, "cygnss_48hr_pvt")
        _CYGNSS_48HR_TIMESPAN_CACHE[] = (t0_s=traj.t0_s, tf_s=traj.tf_s)

        scenario = _base_scenario_dict("cygnss_48hr_pvt", traj.telemetry_path)
        merge!(scenario, Dict{String, Any}(
            "planet" => "earth",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 50,
            "gravity_harmonics_order" => 50,
            "gravity_harmonics_file" => _GMAT_HARMONICS_EARTH_FILE,
            "nbody_bodies" => Any[],
            "orbit_altitude_mode" => "oblate",
            "drag_enabled" => false,
            "include_wind" => false,
            "max_points_quick" => 200000,
            "max_points_full" => 200000,
            "initial_time" => Dict{String, Any}(
                "year" => 2025,
                "month" => 6,
                "day" => 6,
                "hour" => 0,
                "minute" => 0,
                "second" => 0.0
            )
        ))

        manifest_path = joinpath(tmp, "cygnss_48hr_manifest.toml")
        open(manifest_path, "w") do io
            TOML.print(io, Dict{String, Any}("scenarios" => Any[scenario]))
        end

        req = TV.VerificationRequest(
            profile=TEST_MODE,
            out_summary=joinpath(tmp, "summary.csv"),
            out_errors=joinpath(tmp, "errors.csv"),
            manifest_path=manifest_path,
            enforce=false,
            generate_plots=false
        )

        solver_env = _telemetry_solver_env_overrides()
        return withenv(pairs(solver_env)...) do
            TV.run_verification(req)
        end
    end

    _CYGNSS_48HR_RESULT_CACHE[] = result
    _CYGNSS_48HR_SUMMARY_CACHE[] = result.summary
    return result
end

function _run_cygnss_48hr_once()::DataFrame
    return _run_cygnss_48hr_result_once().summary
end

function _run_cygnss_96hr_result_once()::TV.VerificationResult
    if _CYGNSS_96HR_RESULT_CACHE[] !== nothing
        return _CYGNSS_96HR_RESULT_CACHE[]
    end

    @test isfile(_CYGNSS_96HR_TELEMETRY_FEATHER)

    result = mktempdir() do tmp
        traj = _build_cygnss_96hr_reference(tmp, "cygnss_96hr_pvt")
        _CYGNSS_96HR_TIMESPAN_CACHE[] = (t0_s=traj.t0_s, tf_s=traj.tf_s)

        scenario = _base_scenario_dict("cygnss_96hr_pvt", traj.telemetry_path)
        merge!(scenario, Dict{String, Any}(
            "planet" => "earth",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 50,
            "gravity_harmonics_order" => 50,
            "gravity_harmonics_file" => _GMAT_HARMONICS_EARTH_FILE,
            "nbody_bodies" => Any[],
            "orbit_altitude_mode" => "oblate",
            "drag_enabled" => false,
            "include_wind" => false,
            "max_points_quick" => 500000,
            "max_points_full" => 500000,
            "initial_time" => _initial_time_dict(traj.initial_time)
        ))

        manifest_path = joinpath(tmp, "cygnss_96hr_manifest.toml")
        open(manifest_path, "w") do io
            TOML.print(io, Dict{String, Any}("scenarios" => Any[scenario]))
        end

        req = TV.VerificationRequest(
            profile=TEST_MODE,
            out_summary=joinpath(tmp, "summary.csv"),
            out_errors=joinpath(tmp, "errors.csv"),
            manifest_path=manifest_path,
            enforce=false,
            generate_plots=false
        )

        solver_env = _telemetry_solver_env_overrides()
        return withenv(pairs(solver_env)...) do
            TV.run_verification(req)
        end
    end

    _CYGNSS_96HR_RESULT_CACHE[] = result
    _CYGNSS_96HR_SUMMARY_CACHE[] = result.summary
    return result
end

function _run_cygnss_96hr_once()::DataFrame
    return _run_cygnss_96hr_result_once().summary
end

function _run_cygnss_cyg04_96hr_result_once()::TV.VerificationResult
    if _CYGNSS_CYG04_96HR_RESULT_CACHE[] !== nothing
        return _CYGNSS_CYG04_96HR_RESULT_CACHE[]
    end

    @test isfile(_CYGNSS_CYG04_96HR_TELEMETRY_FEATHER)

    result = mktempdir() do tmp
        traj = _build_cygnss_cyg04_96hr_inertial_reference(tmp, "cyg04_96hr_inertial")
        _CYGNSS_CYG04_96HR_TIMESPAN_CACHE[] = (t0_s=traj.t0_s, tf_s=traj.tf_s)

        scenario = _base_scenario_dict("cyg04_96hr_inertial", traj.telemetry_path)
        merge!(scenario, Dict{String, Any}(
            "planet" => "earth",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 50,
            "gravity_harmonics_order" => 50,
            "gravity_harmonics_file" => _GMAT_HARMONICS_EARTH_FILE,
            "nbody_bodies" => Any[],
            "orbit_altitude_mode" => "oblate",
            "drag_enabled" => false,
            "include_wind" => false,
            "max_points_quick" => 500000,
            "max_points_full" => 500000,
            "initial_time" => Dict{String, Any}(
                "year" => 2025,
                "month" => 6,
                "day" => 6,
                "hour" => 0,
                "minute" => 0,
                "second" => 0.0
            )
        ))

        manifest_path = joinpath(tmp, "cyg04_96hr_inertial_manifest.toml")
        open(manifest_path, "w") do io
            TOML.print(io, Dict{String, Any}("scenarios" => Any[scenario]))
        end

        req = TV.VerificationRequest(
            profile=TEST_MODE,
            out_summary=joinpath(tmp, "summary.csv"),
            out_errors=joinpath(tmp, "errors.csv"),
            manifest_path=manifest_path,
            enforce=false,
            generate_plots=false
        )

        solver_env = _telemetry_solver_env_overrides()
        return withenv(pairs(solver_env)...) do
            TV.run_verification(req)
        end
    end

    _CYGNSS_CYG04_96HR_RESULT_CACHE[] = result
    _CYGNSS_CYG04_96HR_SUMMARY_CACHE[] = result.summary
    return result
end

function _run_cygnss_cyg04_96hr_once()::DataFrame
    return _run_cygnss_cyg04_96hr_result_once().summary
end

function _run_cygnss_gmat_csv_result_once()::TV.VerificationResult
    if _CYGNSS_GMAT_RESULT_CACHE[] !== nothing
        return _CYGNSS_GMAT_RESULT_CACHE[]
    end

    @test isfile(_CYGNSS_GMAT_COMPARISON_PATH)

    result = mktempdir() do tmp
        traj = _build_time_aligned_reference(
            _CYGNSS_GMAT_COMPARISON_PATH,
            "earth",
            tmp,
            "cygnss_48hr_gmat_csv";
            sma_km=6818.860956945965,
            ecc=0.0004790007068221837,
            inc_deg=34.93573031327011,
            aop_deg=140.63502513151334,
            raan_deg=177.37117999662348,
            ta_deg=276.65207622217173
        )

        telem_df = _read_tabular(_CYGNSS_GMAT_COMPARISON_PATH)
        tcol = _required_column(telem_df, ["DefaultSC.ElapsedSecs", "Sat.ElapsedSecs"])
        t = sort(Float64.(telem_df[!, tcol]))
        _CYGNSS_GMAT_TIMESPAN_CACHE[] = (t0_s=t[1], tf_s=t[end])

        scenario = _base_scenario_dict("cygnss_48hr_gmat_csv", traj.telemetry_path)
        merge!(scenario, Dict{String, Any}(
            "planet" => "earth",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 50,
            "gravity_harmonics_order" => 50,
            "gravity_harmonics_file" => _GMAT_HARMONICS_EARTH_FILE,
            "nbody_bodies" => Any[],
            "orbit_altitude_mode" => "oblate",
            "drag_enabled" => false,
            "include_wind" => false,
            "max_points_quick" => 200000,
            "max_points_full" => 200000,
            "initial_time" => Dict{String, Any}(
                "year" => 2025,
                "month" => 6,
                "day" => 6,
                "hour" => 0,
                "minute" => 0,
                "second" => 0.0
            )
        ))

        manifest_path = joinpath(tmp, "cygnss_48hr_gmat_manifest.toml")
        open(manifest_path, "w") do io
            TOML.print(io, Dict{String, Any}("scenarios" => Any[scenario]))
        end

        req = TV.VerificationRequest(
            profile=TEST_MODE,
            out_summary=joinpath(tmp, "summary.csv"),
            out_errors=joinpath(tmp, "errors.csv"),
            manifest_path=manifest_path,
            enforce=false,
            generate_plots=false
        )

        solver_env = _telemetry_solver_env_overrides()
        return withenv(pairs(solver_env)...) do
            TV.run_verification(req)
        end
    end

    _CYGNSS_GMAT_RESULT_CACHE[] = result
    _CYGNSS_GMAT_SUMMARY_CACHE[] = result.summary
    return result
end

function _run_cygnss_gmat_csv_once()::DataFrame
    return _run_cygnss_gmat_csv_result_once().summary
end

function _scenario_position_rmse(errors::DataFrame, scenario_name::String)
    rows = errors[(errors.scenario .== scenario_name) .& in.(errors.event, Ref(["state_x_time", "state_y_time", "state_z_time"])), :]
    xrow = sort(rows[rows.event .== "state_x_time", :], :idx)
    yrow = sort(rows[rows.event .== "state_y_time", :], :idx)
    zrow = sort(rows[rows.event .== "state_z_time", :], :idx)

    n = min(nrow(xrow), nrow(yrow), nrow(zrow))
    @test n >= 2

    ex = Float64.(xrow.error_km[1:n])
    ey = Float64.(yrow.error_km[1:n])
    ez = Float64.(zrow.error_km[1:n])

    pos_sq = ex .^ 2 .+ ey .^ 2 .+ ez .^ 2
    first_step_error_km = sqrt(pos_sq[2])
    full_rmse_km = sqrt(mean(pos_sq))

    return (first_step_error_km=first_step_error_km, full_rmse_km=full_rmse_km, n_points=n)
end

@inline function _wrap_angle_deg(θ::Float64)::Float64
    θw = mod(θ + 180.0, 360.0)
    return θw < 0.0 ? θw + 180.0 : θw - 180.0
end

function _differentiate_position_series(pos::Vector{Float64}, time_s::Vector{Float64})::Vector{Float64}
    n = length(pos)
    n == length(time_s) || throw(ArgumentError("position/time length mismatch: $n vs $(length(time_s))"))
    n >= 2 || throw(ArgumentError("Need at least two samples to differentiate position series."))

    vel = Vector{Float64}(undef, n)
    if n == 2
        dt = time_s[2] - time_s[1]
        vel[1] = (pos[2] - pos[1]) / dt
        vel[2] = vel[1]
        return vel
    end

    vel[1] = (pos[2] - pos[1]) / (time_s[2] - time_s[1])
    for i in 2:(n - 1)
        vel[i] = (pos[i + 1] - pos[i - 1]) / (time_s[i + 1] - time_s[i - 1])
    end
    vel[n] = (pos[n] - pos[n - 1]) / (time_s[n] - time_s[n - 1])
    return vel
end

function _scenario_orbital_element_rmse(errors::DataFrame, scenario_name::String)
    rows = errors[(errors.scenario .== scenario_name) .& in.(errors.event, Ref(["state_x_time", "state_y_time", "state_z_time"])), :]
    xrow = sort(rows[rows.event .== "state_x_time", :], :idx)
    yrow = sort(rows[rows.event .== "state_y_time", :], :idx)
    zrow = sort(rows[rows.event .== "state_z_time", :], :idx)

    n = min(nrow(xrow), nrow(yrow), nrow(zrow))
    @test n >= 3

    t_s = Float64.(xrow.telemetry_axis[1:n])
    x_tel = Float64.(xrow.telemetry_value_km[1:n])
    y_tel = Float64.(yrow.telemetry_value_km[1:n])
    z_tel = Float64.(zrow.telemetry_value_km[1:n])
    x_sim = Float64.(xrow.sim_interp_value_km[1:n])
    y_sim = Float64.(yrow.sim_interp_value_km[1:n])
    z_sim = Float64.(zrow.sim_interp_value_km[1:n])

    vx_tel = _differentiate_position_series(x_tel, t_s)
    vy_tel = _differentiate_position_series(y_tel, t_s)
    vz_tel = _differentiate_position_series(z_tel, t_s)
    vx_sim = _differentiate_position_series(x_sim, t_s)
    vy_sim = _differentiate_position_series(y_sim, t_s)
    vz_sim = _differentiate_position_series(z_sim, t_s)

    planet = TV._planet_from_name(_scenario_planet_name(scenario_name))
    sma_sq = Vector{Float64}(undef, n)
    ecc_sq = Vector{Float64}(undef, n)
    inc_sq = Vector{Float64}(undef, n)
    aop_sq = Vector{Float64}(undef, n)
    raan_sq = Vector{Float64}(undef, n)
    ta_sq = Vector{Float64}(undef, n)

    for i in 1:n
        oe_tel = TV.rvtoorbitalelement(
            SVector{3, Float64}(x_tel[i], y_tel[i], z_tel[i]) .* 1.0e3,
            SVector{3, Float64}(vx_tel[i], vy_tel[i], vz_tel[i]) .* 1.0e3,
            planet
        )
        oe_sim = TV.rvtoorbitalelement(
            SVector{3, Float64}(x_sim[i], y_sim[i], z_sim[i]) .* 1.0e3,
            SVector{3, Float64}(vx_sim[i], vy_sim[i], vz_sim[i]) .* 1.0e3,
            planet
        )

        dsma_km = (oe_sim[1] - oe_tel[1]) * 1.0e-3
        decc = oe_sim[2] - oe_tel[2]
        dinc_deg = _wrap_angle_deg(rad2deg(oe_sim[3] - oe_tel[3]))
        daop_deg = _wrap_angle_deg(rad2deg(oe_sim[5] - oe_tel[5]))
        draan_deg = _wrap_angle_deg(rad2deg(oe_sim[4] - oe_tel[4]))
        dta_deg = _wrap_angle_deg(rad2deg(oe_sim[6] - oe_tel[6]))

        sma_sq[i] = dsma_km^2
        ecc_sq[i] = decc^2
        inc_sq[i] = dinc_deg^2
        aop_sq[i] = daop_deg^2
        raan_sq[i] = draan_deg^2
        ta_sq[i] = dta_deg^2
    end

    return (
        first_step=(
            sma_km=sqrt(sma_sq[2]),
            ecc=sqrt(ecc_sq[2]),
            inc_deg=sqrt(inc_sq[2]),
            aop_deg=sqrt(aop_sq[2]),
            raan_deg=sqrt(raan_sq[2]),
            ta_deg=sqrt(ta_sq[2])
        ),
        full_rmse=(
            sma_km=sqrt(mean(sma_sq)),
            ecc=sqrt(mean(ecc_sq)),
            inc_deg=sqrt(mean(inc_sq)),
            aop_deg=sqrt(mean(aop_sq)),
            raan_deg=sqrt(mean(raan_sq)),
            ta_deg=sqrt(mean(ta_sq))
        ),
        n_points=n
    )
end

function _extract_position_error_series(errors::DataFrame, scenario_name::String)
    rows = errors[errors.scenario .== scenario_name, :]
    xrow = sort(rows[rows.event .== "state_x_time", :], :idx)
    yrow = sort(rows[rows.event .== "state_y_time", :], :idx)
    zrow = sort(rows[rows.event .== "state_z_time", :], :idx)

    n = min(nrow(xrow), nrow(yrow), nrow(zrow))
    n >= 2 || error("Need at least two aligned position-error samples for plotting scenario $(scenario_name).")

    t_s = Float64.(xrow.telemetry_axis[1:n])
    ex = Float64.(xrow.error_km[1:n])
    ey = Float64.(yrow.error_km[1:n])
    ez = Float64.(zrow.error_km[1:n])

    total_pos_err = sqrt.(ex .^ 2 .+ ey .^ 2 .+ ez .^ 2)
    running_total_rmse = sqrt.(cumsum(total_pos_err .^ 2) ./ collect(1:n))
    total_rmse = sqrt(mean(total_pos_err .^ 2))

    return (
        t_s=t_s,
        ex=ex,
        ey=ey,
        ez=ez,
        total_pos_err=total_pos_err,
        running_total_rmse=running_total_rmse,
        total_rmse=total_rmse
    )
end

function _iqr_inlier_mask(values::Vector{Float64}; whisker_scale::Float64=1.5)
    finite_idx = findall(isfinite, values)
    isempty(finite_idx) && return falses(length(values))

    finite_vals = values[finite_idx]
    q1 = quantile(finite_vals, 0.25)
    q3 = quantile(finite_vals, 0.75)
    iqr = q3 - q1
    if !isfinite(iqr) || iqr <= 0.0
        mask = falses(length(values))
        mask[finite_idx] .= true
        return mask
    end

    lo = q1 - whisker_scale * iqr
    hi = q3 + whisker_scale * iqr
    mask = falses(length(values))
    @inbounds for i in eachindex(values)
        v = values[i]
        mask[i] = isfinite(v) && lo <= v <= hi
    end
    return mask
end

function _filter_position_series_for_plot(series)
    mask = _iqr_inlier_mask(series.total_pos_err)
    sum(mask) >= 2 || (mask .= true)
    return (
        t_s=series.t_s[mask],
        ex=series.ex[mask],
        ey=series.ey[mask],
        ez=series.ez[mask],
        total_pos_err=series.total_pos_err[mask],
        running_total_rmse=sqrt.(cumsum(series.total_pos_err[mask] .^ 2) ./ collect(1:sum(mask))),
        total_rmse=series.total_rmse,
        inlier_count=sum(mask),
        total_count=length(mask)
    )
end

function _plot_cygnss_rtn_from_reference(
    errors::DataFrame,
    scenario_name::String,
    outpath::String,
    ref_series;
    title::String
)::String
    rows = errors[errors.scenario .== scenario_name, :]
    xrow = sort(rows[rows.event .== "state_x_time", :], :idx)
    yrow = sort(rows[rows.event .== "state_y_time", :], :idx)
    zrow = sort(rows[rows.event .== "state_z_time", :], :idx)

    n = min(nrow(xrow), nrow(yrow), nrow(zrow))
    n >= 2 || error("Need at least two aligned CYGNSS error samples for RTN plotting.")

    t_s = Float64.(xrow.telemetry_axis[1:n])
    t_hr = t_s ./ 3600.0
    ex = Float64.(xrow.error_km[1:n])
    ey = Float64.(yrow.error_km[1:n])
    ez = Float64.(zrow.error_km[1:n])

    rx = _interp_series_linear(ref_series.t_rel, ref_series.x_km, t_s)
    ry = _interp_series_linear(ref_series.t_rel, ref_series.y_km, t_s)
    rz = _interp_series_linear(ref_series.t_rel, ref_series.z_km, t_s)
    vx = _interp_series_linear(ref_series.t_rel, ref_series.vx_kmps, t_s)
    vy = _interp_series_linear(ref_series.t_rel, ref_series.vy_kmps, t_s)
    vz = _interp_series_linear(ref_series.t_rel, ref_series.vz_kmps, t_s)

    er = zeros(n)
    et = zeros(n)
    en = zeros(n)
    for i in 1:n
        r = [rx[i], ry[i], rz[i]]
        v = [vx[i], vy[i], vz[i]]
        e = [ex[i], ey[i], ez[i]]
        rhat = r / norm(r)
        h = cross(r, v)
        nhat = h / norm(h)
        that = cross(nhat, rhat)
        er[i] = dot(e, rhat)
        et[i] = dot(e, that)
        en[i] = dot(e, nhat)
    end

    total_pos_err = sqrt.(er .^ 2 .+ et .^ 2 .+ en .^ 2)
    total_rmse = sqrt(mean(total_pos_err .^ 2))
    mask = _iqr_inlier_mask(total_pos_err)
    sum(mask) >= 2 || (mask .= true)
    t_hr = t_hr[mask]
    er = er[mask]
    et = et[mask]
    en = en[mask]
    total_pos_err = total_pos_err[mask]
    running_total_rmse = sqrt.(cumsum(total_pos_err .^ 2) ./ collect(1:length(total_pos_err)))
    p = plot(t_hr, er, label="R error", xlabel="Time (hr)", ylabel="Error (km)", title=title, lw=1.2, alpha=0.9, legend=:topright)
    plot!(p, t_hr, et, label="T error", lw=1.2, alpha=0.9)
    plot!(p, t_hr, en, label="N error", lw=1.2, alpha=0.9)
    plot!(p, t_hr, running_total_rmse, label="filtered running RMSE", lw=2.0, color=:black)
    hline!(p, [total_rmse], label="total RMSE", linestyle=:dash, color=:black)

    mkpath(dirname(outpath))
    savefig(p, outpath)
    return outpath
end

function _plot_cygnss_error_timeseries(errors::DataFrame, outpath::String)::String
    series = _filter_position_series_for_plot(_extract_position_error_series(errors, "cygnss_48hr_pvt"))
    t_hr = series.t_s ./ 3600.0

    p = plot(
        t_hr,
        series.ex,
        label="x error",
        xlabel="Time (hr)",
        ylabel="Error (km)",
        title="CYGNSS 48hr Position Error Time Series",
        lw=1.2,
        alpha=0.9,
        legend=:topright
    )
    plot!(p, t_hr, series.ey, label="y error", lw=1.2, alpha=0.9)
    plot!(p, t_hr, series.ez, label="z error", lw=1.2, alpha=0.9)
    plot!(p, t_hr, series.running_total_rmse, label="filtered running RMSE", lw=2.0, color=:black)
    hline!(p, [series.total_rmse], label="total RMSE", linestyle=:dash, color=:black)

    mkpath(dirname(outpath))
    savefig(p, outpath)
    return outpath
end

function _plot_position_error_timeseries(
    errors::DataFrame,
    scenario_name::String,
    outpath::String;
    title::String="Position Error Time Series",
    time_unit::Symbol=:hr
)::String
    series = _filter_position_series_for_plot(_extract_position_error_series(errors, scenario_name))
    t_s = series.t_s
    t_axis, x_label = if time_unit == :s
        (t_s, "Time (s)")
    else
        (t_s ./ 3600.0, "Time (hr)")
    end

    p = plot(
        t_axis,
        series.ex,
        label="x error",
        xlabel=x_label,
        ylabel="Error (km)",
        title=title,
        lw=1.2,
        alpha=0.9,
        legend=:topright
    )
    plot!(p, t_axis, series.ey, label="y error", lw=1.2, alpha=0.9)
    plot!(p, t_axis, series.ez, label="z error", lw=1.2, alpha=0.9)
    plot!(p, t_axis, series.running_total_rmse, label="total running RMSE", lw=2.0, color=:black)
    hline!(p, [series.total_rmse], label="total RMSE", linestyle=:dash, color=:black)

    mkpath(dirname(outpath))
    savefig(p, outpath)
    return outpath
end

function _plot_cygnss_reference_error_comparison(
    cygnss_errors::DataFrame,
    gmat_errors::DataFrame,
    outpath::String
)::String
    cygnss_series = _filter_position_series_for_plot(_extract_position_error_series(cygnss_errors, "cygnss_48hr_pvt"))
    gmat_series = _filter_position_series_for_plot(_extract_position_error_series(gmat_errors, "cygnss_48hr_gmat_csv"))

    t0_s = max(cygnss_series.t_s[1], gmat_series.t_s[1])
    tf_s = min(cygnss_series.t_s[end], gmat_series.t_s[end])
    tf_s > t0_s || error("CYGNSS and GMAT comparison series do not overlap in time.")

    t_common_s = sort(unique(vcat(cygnss_series.t_s, gmat_series.t_s)))
    t_common_s = t_common_s[(t_common_s .>= t0_s) .& (t_common_s .<= tf_s)]
    length(t_common_s) >= 2 || error("Need at least two overlapping samples for CYGNSS reference comparison plotting.")

    cygnss_total = _interp_series_linear(cygnss_series.t_s, cygnss_series.total_pos_err, t_common_s)
    gmat_total = _interp_series_linear(gmat_series.t_s, gmat_series.total_pos_err, t_common_s)
    cygnss_running = _interp_series_linear(cygnss_series.t_s, cygnss_series.running_total_rmse, t_common_s)
    gmat_running = _interp_series_linear(gmat_series.t_s, gmat_series.running_total_rmse, t_common_s)

    t_hr = t_common_s ./ 3600.0

    p = plot(
        t_hr,
        cygnss_total,
        label="SpaceAGORA vs CYGNSS total error",
        xlabel="Time (hr)",
        ylabel="Error (km)",
        title="CYGNSS 48hr Error Comparison",
        lw=2.0,
        color=:blue,
        legend=:topright
    )
    plot!(p, t_hr, gmat_total, label="SpaceAGORA vs GMAT total error", lw=2.0, color=:red)
    plot!(p, t_hr, cygnss_running, label="SpaceAGORA vs CYGNSS running RMSE", lw=1.8, linestyle=:dash, color=:blue)
    plot!(p, t_hr, gmat_running, label="SpaceAGORA vs GMAT running RMSE", lw=1.8, linestyle=:dash, color=:red)

    mkpath(dirname(outpath))
    savefig(p, outpath)
    return outpath
end

function _interp_series_linear(x::Vector{Float64}, y::Vector{Float64}, xq::Vector{Float64})::Vector{Float64}
    n = length(x)
    n == length(y) || error("Interpolation grid/value length mismatch.")
    n >= 2 || error("Need at least two interpolation grid points.")

    out = similar(xq)
    for i in eachindex(xq)
        xi = xq[i]
        if xi <= x[1]
            out[i] = y[1]
            continue
        elseif xi >= x[end]
            out[i] = y[end]
            continue
        end

        j = searchsortedlast(x, xi)
        x0 = x[j]
        x1 = x[j + 1]
        y0 = y[j]
        y1 = y[j + 1]
        α = (xi - x0) / (x1 - x0)
        out[i] = (1.0 - α) * y0 + α * y1
    end
    return out
end

function _plot_cygnss_error_timeseries_rtn(errors::DataFrame, outpath::String)::String
    return _plot_cygnss_rtn_from_reference(
        errors,
        "cygnss_48hr_pvt",
        outpath,
        _load_cygnss_cyg04_96hr_inertial_series();
        title="CYGNSS 48hr Position Error Time Series (RTN)"
    )
end

function _plot_cygnss_96hr_error_timeseries_rtn(errors::DataFrame, outpath::String; title::String="CYGNSS 96hr Position Error Time Series (RTN)")::String
    return _plot_cygnss_rtn_from_reference(
        errors,
        "cygnss_96hr_pvt",
        outpath,
        _load_cygnss_96hr_j2000_series();
        title=title
    )
end

function _plot_cyg04_96hr_error_timeseries_rtn(errors::DataFrame, outpath::String; title::String="CYG04 96hr Position Error Time Series (RTN)")::String
    return _plot_cygnss_rtn_from_reference(
        errors,
        "cyg04_96hr_inertial",
        outpath,
        _load_cygnss_cyg04_96hr_inertial_series();
        title=title
    )
end

@inline function _wrap_pm180(x_deg::Float64)::Float64
    y = mod(x_deg + 180.0, 360.0)
    return y - 180.0
end

function _unwrap_deg(series::Vector{Float64})::Vector{Float64}
    n = length(series)
    n == 0 && return Float64[]
    out = copy(series)
    offset = 0.0
    for i in 2:n
        if !isfinite(series[i]) || !isfinite(series[i - 1])
            out[i] = series[i]
            continue
        end
        jump = series[i] - series[i - 1]
        if jump > 180.0
            offset -= 360.0
        elseif jump < -180.0
            offset += 360.0
        end
        out[i] = series[i] + offset
    end
    return out
end

function _plot_cygnss_orbital_elements_from_reference(
    errors::DataFrame,
    scenario_name::String,
    ref_series,
    outpath_elements::String,
    outpath_errors::String;
    title_suffix::String=""
)::Tuple{String, String}
    rows = errors[errors.scenario .== scenario_name, :]
    xrow = sort(rows[rows.event .== "state_x_time", :], :idx)
    yrow = sort(rows[rows.event .== "state_y_time", :], :idx)
    zrow = sort(rows[rows.event .== "state_z_time", :], :idx)

    n = min(nrow(xrow), nrow(yrow), nrow(zrow))
    n >= 2 || error("Need at least two aligned CYGNSS samples for orbital element plotting.")

    t_s = Float64.(xrow.telemetry_axis[1:n])
    t_hr = t_s ./ 3600.0

    rx = _interp_series_linear(ref_series.t_rel, ref_series.x_km, t_s)
    ry = _interp_series_linear(ref_series.t_rel, ref_series.y_km, t_s)
    rz = _interp_series_linear(ref_series.t_rel, ref_series.z_km, t_s)
    vx = _interp_series_linear(ref_series.t_rel, ref_series.vx_kmps, t_s)
    vy = _interp_series_linear(ref_series.t_rel, ref_series.vy_kmps, t_s)
    vz = _interp_series_linear(ref_series.t_rel, ref_series.vz_kmps, t_s)

    planet = TV._planet_from_name("earth")

    sma_ref_sim = similar(t_s)
    ecc_ref_sim = similar(t_s)
    inc_ref_sim = similar(t_s)
    aop_ref_sim = similar(t_s)
    raan_ref_sim = similar(t_s)
    ta_ref_sim = similar(t_s)

    for i in 1:n
        r = SVector{3, Float64}(rx[i], ry[i], rz[i]) .* 1.0e3
        v = SVector{3, Float64}(vx[i], vy[i], vz[i]) .* 1.0e3
        try
            oe = TV.rvtoorbitalelement(r, v, planet)
            sma_ref_sim[i] = oe[1] * 1.0e-3
            ecc_ref_sim[i] = oe[2]
            inc_ref_sim[i] = rad2deg(oe[3])
            raan_ref_sim[i] = rad2deg(oe[4])
            aop_ref_sim[i] = rad2deg(oe[5])
            ta_ref_sim[i] = rad2deg(oe[6])
        catch
            sma_ref_sim[i] = NaN
            ecc_ref_sim[i] = NaN
            inc_ref_sim[i] = NaN
            aop_ref_sim[i] = NaN
            raan_ref_sim[i] = NaN
            ta_ref_sim[i] = NaN
        end
    end

    ex = Float64.(xrow.error_km[1:n])
    ey = Float64.(yrow.error_km[1:n])
    ez = Float64.(zrow.error_km[1:n])
    rx_sim = rx .+ ex
    ry_sim = ry .+ ey
    rz_sim = rz .+ ez

    sma_sim = similar(t_s)
    ecc_sim = similar(t_s)
    inc_sim = similar(t_s)
    aop_sim = similar(t_s)
    raan_sim = similar(t_s)
    ta_sim = similar(t_s)

    for i in 1:n
        r = SVector{3, Float64}(rx_sim[i], ry_sim[i], rz_sim[i]) .* 1.0e3
        v = SVector{3, Float64}(vx[i], vy[i], vz[i]) .* 1.0e3
        try
            oe = TV.rvtoorbitalelement(r, v, planet)
            sma_sim[i] = oe[1] * 1.0e-3
            ecc_sim[i] = oe[2]
            inc_sim[i] = rad2deg(oe[3])
            raan_sim[i] = rad2deg(oe[4])
            aop_sim[i] = rad2deg(oe[5])
            ta_sim[i] = rad2deg(oe[6])
        catch
            sma_sim[i] = NaN
            ecc_sim[i] = NaN
            inc_sim[i] = NaN
            aop_sim[i] = NaN
            raan_sim[i] = NaN
            ta_sim[i] = NaN
        end
    end

    finite_ref = count(isfinite, sma_ref_sim)
    finite_sim = count(isfinite, sma_sim)
    finite_ref >= max(2, Int(floor(0.95 * n))) || error("Too many non-finite reference orbital element samples: $finite_ref / $n")
    finite_sim >= max(2, Int(floor(0.95 * n))) || error("Too many non-finite simulated orbital element samples: $finite_sim / $n")

    raan_ref_wrap = mod.(raan_ref_sim, 360.0)
    aop_ref_wrap = mod.(aop_ref_sim, 360.0)
    ta_ref_wrap = mod.(ta_ref_sim, 360.0)
    raan_sim_wrap = mod.(raan_sim, 360.0)
    aop_sim_wrap = mod.(aop_sim, 360.0)
    ta_sim_wrap = mod.(ta_sim, 360.0)

    da = sma_sim .- sma_ref_sim
    de = ecc_sim .- ecc_ref_sim
    di = inc_sim .- inc_ref_sim
    dw = _wrap_pm180.(aop_sim_wrap .- aop_ref_wrap)
    dW = _wrap_pm180.(raan_sim_wrap .- raan_ref_wrap)
    df = _wrap_pm180.(ta_sim_wrap .- ta_ref_wrap)

    pos_err = sqrt.(ex .^ 2 .+ ey .^ 2 .+ ez .^ 2)
    mask = _iqr_inlier_mask(pos_err)
    sum(mask) >= 2 || (mask .= true)
    t_hr = t_hr[mask]
    sma_ref_sim = sma_ref_sim[mask]
    ecc_ref_sim = ecc_ref_sim[mask]
    inc_ref_sim = inc_ref_sim[mask]
    aop_ref_wrap = aop_ref_wrap[mask]
    raan_ref_wrap = raan_ref_wrap[mask]
    ta_ref_wrap = ta_ref_wrap[mask]
    sma_sim = sma_sim[mask]
    ecc_sim = ecc_sim[mask]
    inc_sim = inc_sim[mask]
    aop_sim_wrap = aop_sim_wrap[mask]
    raan_sim_wrap = raan_sim_wrap[mask]
    ta_sim_wrap = ta_sim_wrap[mask]
    da = da[mask]
    de = de[mask]
    di = di[mask]
    dw = dw[mask]
    dW = dW[mask]
    df = df[mask]

    layout_3x2 = @layout [a b; c d; e f]
    pe1 = plot(t_hr, sma_sim, label="simulated", xlabel="Time (hr)", ylabel="SMA (km)", title="Semi-Major Axis$title_suffix", legend=:best, lw=1.2)
    plot!(pe1, t_hr, sma_ref_sim, label="reference", lw=1.2, alpha=0.75)
    pe2 = plot(t_hr, ecc_sim, label="simulated", xlabel="Time (hr)", ylabel="Ecc", title="Eccentricity$title_suffix", legend=:best, lw=1.2)
    plot!(pe2, t_hr, ecc_ref_sim, label="reference", lw=1.2, alpha=0.75)
    pe3 = plot(t_hr, inc_sim, label="simulated", xlabel="Time (hr)", ylabel="Inc (deg)", title="Inclination$title_suffix", legend=:best, lw=1.2)
    plot!(pe3, t_hr, inc_ref_sim, label="reference", lw=1.2, alpha=0.75)
    pe4 = plot(t_hr, aop_sim_wrap, label="simulated", xlabel="Time (hr)", ylabel="AoP (deg)", title="Argument of Perigee$title_suffix", legend=:best, lw=1.2, ylim=(0.0, 360.0))
    plot!(pe4, t_hr, aop_ref_wrap, label="reference", lw=1.2, alpha=0.75)
    pe5 = plot(t_hr, raan_sim_wrap, label="simulated", xlabel="Time (hr)", ylabel="RAAN (deg)", title="RAAN$title_suffix", legend=:best, lw=1.2, ylim=(0.0, 360.0))
    plot!(pe5, t_hr, raan_ref_wrap, label="reference", lw=1.2, alpha=0.75)
    pe6 = plot(t_hr, ta_sim_wrap, label="simulated", xlabel="Time (hr)", ylabel="TA (deg)", title="True Anomaly$title_suffix", legend=:best, lw=1.2, ylim=(0.0, 360.0))
    plot!(pe6, t_hr, ta_ref_wrap, label="reference", lw=1.2, alpha=0.75)
    p_elements = plot(pe1, pe2, pe3, pe4, pe5, pe6, layout=layout_3x2, size=(1500, 1100))

    per1 = plot(t_hr, da, xlabel="Time (hr)", ylabel="Delta SMA (km)", title="SMA Error", lw=1.2, legend=false, color=:red)
    hline!(per1, [0.0], linestyle=:dash, color=:black, alpha=0.5)
    per2 = plot(t_hr, de, xlabel="Time (hr)", ylabel="Delta Ecc", title="Eccentricity Error", lw=1.2, legend=false, color=:red)
    hline!(per2, [0.0], linestyle=:dash, color=:black, alpha=0.5)
    per3 = plot(t_hr, di, xlabel="Time (hr)", ylabel="Delta Inc (deg)", title="Inclination Error", lw=1.2, legend=false, color=:red)
    hline!(per3, [0.0], linestyle=:dash, color=:black, alpha=0.5)
    per4 = plot(t_hr, dw, xlabel="Time (hr)", ylabel="Delta AoP (deg)", title="AoP Error", lw=1.2, legend=false, color=:red)
    hline!(per4, [0.0], linestyle=:dash, color=:black, alpha=0.5)
    per5 = plot(t_hr, dW, xlabel="Time (hr)", ylabel="Delta RAAN (deg)", title="RAAN Error", lw=1.2, legend=false, color=:red)
    hline!(per5, [0.0], linestyle=:dash, color=:black, alpha=0.5)
    per6 = plot(t_hr, df, xlabel="Time (hr)", ylabel="Delta TA (deg)", title="True Anomaly Error", lw=1.2, legend=false, color=:red)
    hline!(per6, [0.0], linestyle=:dash, color=:black, alpha=0.5)
    p_errors = plot(per1, per2, per3, per4, per5, per6, layout=layout_3x2, size=(1500, 1100))

    mkpath(dirname(outpath_elements))
    savefig(p_elements, outpath_elements)
    mkpath(dirname(outpath_errors))
    savefig(p_errors, outpath_errors)
    return (outpath_elements, outpath_errors)
end

function _plot_cygnss_orbital_elements(
    errors::DataFrame,
    outpath_elements::String,
    outpath_errors::String
)::Tuple{String, String}
    return _plot_cygnss_orbital_elements_from_reference(
        errors,
        "cygnss_48hr_pvt",
        _load_cygnss_cyg04_96hr_inertial_series(),
        outpath_elements,
        outpath_errors
    )
end

function _plot_cygnss_96hr_orbital_elements(
    errors::DataFrame,
    outpath_elements::String,
    outpath_errors::String
)::Tuple{String, String}
    return _plot_cygnss_orbital_elements_from_reference(
        errors,
        "cygnss_96hr_pvt",
        _load_cygnss_96hr_j2000_series(),
        outpath_elements,
        outpath_errors
    )
end

function _plot_cyg04_96hr_orbital_elements(
    errors::DataFrame,
    outpath_elements::String,
    outpath_errors::String
)::Tuple{String, String}
    return _plot_cygnss_orbital_elements_from_reference(
        errors,
        "cyg04_96hr_inertial",
        _load_cygnss_cyg04_96hr_inertial_series(),
        outpath_elements,
        outpath_errors;
        title_suffix=" (CYG04)"
    )
end

function _plot_gmat_matrix_error_timeseries(errors::DataFrame, outpath::String)::String
    planets = ["earth", "mars", "venus", "moon"]
    cases   = ["j0_tbfalse", "j2_tbfalse", "j50_tbfalse", "j0_tbtrue", "j2_tbtrue", "j50_tbtrue"]

    n_rows = length(cases)
    n_cols = length(planets)
    subplots = Vector{Plots.Plot}(undef, n_rows * n_cols)

    for (ci, case) in enumerate(cases)
        for (pi, planet) in enumerate(planets)
            scenario_name = "$(planet)_$(case)"
            idx = (ci - 1) * n_cols + pi

            rows = errors[errors.scenario .== scenario_name, :]
            xrow = sort(rows[rows.event .== "state_x_time", :], :idx)
            yrow = sort(rows[rows.event .== "state_y_time", :], :idx)
            zrow = sort(rows[rows.event .== "state_z_time", :], :idx)

            n = min(nrow(xrow), nrow(yrow), nrow(zrow))

            if n < 2
                subplots[idx] = plot(title=scenario_name, legend=false)
                continue
            end

            t_hr = Float64.(xrow.telemetry_axis[1:n]) ./ 3600.0
            ex   = Float64.(xrow.error_km[1:n])
            ey   = Float64.(yrow.error_km[1:n])
            ez   = Float64.(zrow.error_km[1:n])
            pos_err = sqrt.(ex .^ 2 .+ ey .^ 2 .+ ez .^ 2)
            running_rmse = sqrt.(cumsum(pos_err .^ 2) ./ collect(1:n))

            ylabel_str = ci == 1 ? "Error (km)" : ""
            xlabel_str = ci == n_rows ? "Time (hr)" : ""

            sp = plot(
                t_hr, pos_err;
                label="3D error",
                title=scenario_name,
                xlabel=xlabel_str,
                ylabel=ylabel_str,
                lw=1.2,
                alpha=0.85,
                legend=:topright
            )
            plot!(sp, t_hr, running_rmse; label="running RMSE", lw=1.8, color=:black)
            subplots[idx] = sp
        end
    end

    fig = plot(
        subplots...;
        layout=(n_rows, n_cols),
        size=(420 * n_cols, 260 * n_rows),
        titlefontsize=8,
        tickfontsize=6,
        labelfontsize=7,
        legendfontsize=6
    )
    mkpath(dirname(outpath))
    savefig(fig, outpath)
    return outpath
end

@inline function _rmse(a::AbstractVector, b::AbstractVector)::Float64
    @test length(a) == length(b)
    @test !any(ismissing, a)
    @test !any(ismissing, b)
    d = Float64.(collect(a)) .- Float64.(collect(b))
    return sqrt(mean(d .^ 2))
end

@inline function _max_quaternion_norm_deviation(df::DataFrame)::Float64
    q1 = Float64.(df.q1)
    q2 = Float64.(df.q2)
    q3 = Float64.(df.q3)
    q4 = Float64.(df.q4)
    qnorm = sqrt.(q1 .^ 2 .+ q2 .^ 2 .+ q3 .^ 2 .+ q4 .^ 2)
    return maximum(abs.(qnorm .- 1.0))
end

if !_parse_bool_env("SPACEAGORA_SKIP_GMAT_MATRIX", false)

if !_basilisk_reference_available()
    @testset "Basilisk parity matrix" begin
        @test_skip "Basilisk parity references not present under data/telemetry/Basilisk_Examples_Full/; run scripts/dev/fetch_private_telemetry.sh references to sync them, then re-run."
    end
else

@testset "GMAT Early vs Full Error" begin
    result = _run_gmat_scenario_matrix_result_once()
    scenario_names = unique(String.(result.summary.scenario))

    for scenario_name in scenario_names
        metrics = _scenario_position_rmse(result.errors, scenario_name)
        pf_metrics = _scenario_planet_fixed_position_rmse(result.errors, scenario_name)
        oe_metrics = _scenario_orbital_element_rmse(result.errors, scenario_name)
        @test isfinite(metrics.first_step_error_km)
        @test isfinite(metrics.full_rmse_km)
        @test isfinite(pf_metrics.first_step_error_km)
        @test isfinite(pf_metrics.full_rmse_km)
        @test isfinite(oe_metrics.first_step.sma_km)
        @test isfinite(oe_metrics.first_step.ecc)
        @test isfinite(oe_metrics.first_step.inc_deg)
        @test isfinite(oe_metrics.first_step.aop_deg)
        @test isfinite(oe_metrics.first_step.raan_deg)
        @test isfinite(oe_metrics.first_step.ta_deg)
        @test isfinite(oe_metrics.full_rmse.sma_km)
        @test isfinite(oe_metrics.full_rmse.ecc)
        @test isfinite(oe_metrics.full_rmse.inc_deg)
        @test isfinite(oe_metrics.full_rmse.aop_deg)
        @test isfinite(oe_metrics.full_rmse.raan_deg)
        @test isfinite(oe_metrics.full_rmse.ta_deg)
        @test metrics.n_points >= 2
        @test pf_metrics.n_points >= 2
        @test oe_metrics.n_points >= 3
        println("$scenario_name position error [km]: first_step=$(metrics.first_step_error_km), full_rmse=$(metrics.full_rmse_km), planet_fixed_first_step=$(pf_metrics.first_step_error_km), planet_fixed_full_rmse=$(pf_metrics.full_rmse_km)")
        println("$scenario_name orbital element error: first_step=(sma_km=$(oe_metrics.first_step.sma_km), ecc=$(oe_metrics.first_step.ecc), inc_deg=$(oe_metrics.first_step.inc_deg), aop_deg=$(oe_metrics.first_step.aop_deg), raan_deg=$(oe_metrics.first_step.raan_deg), ta_deg=$(oe_metrics.first_step.ta_deg)), full_rmse=(sma_km=$(oe_metrics.full_rmse.sma_km), ecc=$(oe_metrics.full_rmse.ecc), inc_deg=$(oe_metrics.full_rmse.inc_deg), aop_deg=$(oe_metrics.full_rmse.aop_deg), raan_deg=$(oe_metrics.full_rmse.raan_deg), ta_deg=$(oe_metrics.full_rmse.ta_deg))")
    end
end

try
    @testset "GMAT Strict Acceptance All Cases" begin
        summary = _run_gmat_scenario_matrix_once()
        profile = TEST_MODE
        scenario_names = unique(String.(summary.scenario))
        expected_scenarios = _active_gmat_expected_scenario_names()

        @test Set(scenario_names) == expected_scenarios

        for scenario_name in scenario_names
            rows = summary[summary.scenario .== scenario_name, :]
            @test nrow(rows) == 4

            xrow = rows[rows.event .== "state_x_time", :]
            yrow = rows[rows.event .== "state_y_time", :]
            zrow = rows[rows.event .== "state_z_time", :]
            @test nrow(xrow) == 1
            @test nrow(yrow) == 1
            @test nrow(zrow) == 1

            println("$scenario_name strict trajectory error [km]: rmse=(x=$(xrow.rmse_km[1]), y=$(yrow.rmse_km[1]), z=$(zrow.rmse_km[1])) max_abs=(x=$(xrow.max_abs_km[1]), y=$(yrow.max_abs_km[1]), z=$(zrow.max_abs_km[1]))")

            @test sqrt(xrow.rmse_km[1]^2 + yrow.rmse_km[1]^2 + zrow.rmse_km[1]^2) < _strict_position_rmse_limit_km(scenario_name, profile)
        end

        result = _run_gmat_scenario_matrix_result_once()
        matrix_plot_path = joinpath(_GMAT_REPO_ROOT, "output", "gmat_matrix", "gmat_matrix_error_timeseries.png")
        _plot_gmat_matrix_error_timeseries(result.errors, matrix_plot_path)
        @test isfile(matrix_plot_path)
        println("GMAT matrix error timeseries plot: $(matrix_plot_path)")
    end
catch err
    if err isa Test.TestSetException
        println("GMAT Strict Acceptance All Cases reported failures; continuing with remaining testsets.")
        Base.display_error(stderr, err, catch_backtrace())
    else
        rethrow(err)
    end
end

end # Basilisk parity references present

if !_stk_reference_available()
    @testset "STK parity matrix" begin
        @test_skip "STK reference results not present under data/telemetry/stk_results/; skipping."
    end
else

try
    @testset "STK Strict Acceptance All Cases" begin
        summary = _run_stk_scenario_matrix_once()
        profile = TEST_MODE
        scenario_names = unique(String.(summary.scenario))
        expected_scenarios = _active_gmat_expected_scenario_names()

        @test Set(scenario_names) == expected_scenarios

        for scenario_name in scenario_names
            rows = summary[summary.scenario .== scenario_name, :]
            @test nrow(rows) == 4

            xrow = rows[rows.event .== "state_x_time", :]
            yrow = rows[rows.event .== "state_y_time", :]
            zrow = rows[rows.event .== "state_z_time", :]
            @test nrow(xrow) == 1
            @test nrow(yrow) == 1
            @test nrow(zrow) == 1

            println("$scenario_name STK trajectory error [km]: rmse=(x=$(xrow.rmse_km[1]), y=$(yrow.rmse_km[1]), z=$(zrow.rmse_km[1])) max_abs=(x=$(xrow.max_abs_km[1]), y=$(yrow.max_abs_km[1]), z=$(zrow.max_abs_km[1]))")

            @test sqrt(xrow.rmse_km[1]^2 + yrow.rmse_km[1]^2 + zrow.rmse_km[1]^2) < _strict_position_rmse_limit_km(scenario_name, profile)
        end

        result = _run_stk_scenario_matrix_result_once()
        matrix_plot_path = joinpath(_GMAT_REPO_ROOT, "output", "stk_matrix", "stk_matrix_error_timeseries.png")
        _plot_gmat_matrix_error_timeseries(result.errors, matrix_plot_path)
        @test isfile(matrix_plot_path)
        println("STK matrix error timeseries plot: $(matrix_plot_path)")
    end
catch err
    if err isa Test.TestSetException
        println("STK Strict Acceptance All Cases reported failures; continuing with remaining testsets.")
        Base.display_error(stderr, err, catch_backtrace())
    else
        rethrow(err)
    end
end

end # STK reference results present

end # SPACEAGORA_SKIP_GMAT_MATRIX (GMAT Early vs Full Error)

if !_cygnss_private_data_available()
    @testset "CYGNSS scenarios" begin
        @test_skip "CYGNSS private telemetry not present under data/telemetry/CYGNSS/; run scripts/dev/fetch_private_telemetry.sh to sync it, then re-run."
    end
else

@testset "CYGNSS Telemetry Folder Data" begin
    cygnss_feather_path = _CYGNSS_48HR_TELEMETRY_FEATHER
    @test isfile(cygnss_feather_path)

    cygnss_feather = DataFrame(Arrow.Table(cygnss_feather_path))

    @test nrow(cygnss_feather) > 0

    time_col = _required_column(cygnss_feather, ["TIME OFFSET", "time"])
    x_col = _required_column(cygnss_feather, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_X (m)", "pos_ii_1"])
    y_col = _required_column(cygnss_feather, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_Y (m)", "pos_ii_2"])
    z_col = _required_column(cygnss_feather, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_Z (m)", "pos_ii_3"])
    vx_col = _required_column(cygnss_feather, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_X (m/s)", "vel_ii_1"])
    vy_col = _required_column(cygnss_feather, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_Y (m/s)", "vel_ii_2"])
    vz_col = _required_column(cygnss_feather, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_Z (m/s)", "vel_ii_3"])

    t = Float64.(cygnss_feather[!, time_col])
    x = Float64.(cygnss_feather[!, x_col])
    y = Float64.(cygnss_feather[!, y_col])
    z = Float64.(cygnss_feather[!, z_col])
    vx = Float64.(cygnss_feather[!, vx_col])
    vy = Float64.(cygnss_feather[!, vy_col])
    vz = Float64.(cygnss_feather[!, vz_col])

    @test all(isfinite, t)
    @test !any(ismissing, cygnss_feather[!, x_col])
    @test !any(ismissing, cygnss_feather[!, y_col])
    @test !any(ismissing, cygnss_feather[!, z_col])
    @test !any(ismissing, cygnss_feather[!, vx_col])
    @test !any(ismissing, cygnss_feather[!, vy_col])
    @test !any(ismissing, cygnss_feather[!, vz_col])

    @test issorted(t)

    r = sqrt.(x .^ 2 .+ y .^ 2 .+ z .^ 2)
    v = sqrt.(vx .^ 2 .+ vy .^ 2 .+ vz .^ 2)

    @test minimum(r) > 6.0e6
    @test maximum(r) < 8.0e6
    @test minimum(v) > 7.0e3
    @test maximum(v) < 8.0e3
end

@testset "CYGNSS Legacy 48hr Entry Point" begin
    result = _run_cygnss_48hr_result_once()
    summary = result.summary
    errors = result.errors
    @test nrow(summary) == 4
    @test all(summary.scenario .== "cygnss_48hr_pvt")
    @test all(in.(summary.event, Ref(["altitude_time", "state_x_time", "state_y_time", "state_z_time"])))

    span = _CYGNSS_48HR_TIMESPAN_CACHE[]
    @test span !== nothing
    sim_span_s = span.tf_s - span.t0_s
    @test abs(sim_span_s - 48.0 * 3600.0) <= 10.0

    # Ensure the verification run used the full telemetry horizon (not a truncated quick subset).
    telemetry_end_s = maximum(Float64.(summary.telemetry_axis_end))
    @test abs(telemetry_end_s - span.tf_s) <= 1.0

    @test all(summary.n_sim .>= summary.min_eval_points)
    @test all(isfinite.(Float64.(summary.rmse_km)))
    @test all(isfinite.(Float64.(summary.max_abs_km)))

    plot_path = joinpath(_GMAT_REPO_ROOT, "output", "cygnss", "cygnss_48hr_error_timeseries.png")
    _plot_position_error_timeseries(errors, "cygnss_48hr_pvt", plot_path; title="CYGNSS 96hr Position Error Time Series")
    @test isfile(plot_path)
    plot_path_rtn = joinpath(_GMAT_REPO_ROOT, "output", "cygnss", "cygnss_48hr_error_timeseries_rtn.png")
    _plot_cygnss_error_timeseries_rtn(errors, plot_path_rtn)
    @test isfile(plot_path_rtn)
    plot_path_oe = joinpath(_GMAT_REPO_ROOT, "output", "cygnss", "cygnss_48hr_orbital_elements.png")
    plot_path_oe_err = joinpath(_GMAT_REPO_ROOT, "output", "cygnss", "cygnss_48hr_orbital_element_errors.png")
    _plot_cygnss_orbital_elements(errors, plot_path_oe, plot_path_oe_err)
    @test isfile(plot_path_oe)
    @test isfile(plot_path_oe_err)
    gmat_comparison_plot_path = joinpath(_GMAT_REPO_ROOT, "output", "cygnss", "cygnss_48hr_reference_error_comparison.png")
    _plot_cygnss_reference_error_comparison(errors, _run_cygnss_gmat_csv_result_once().errors, gmat_comparison_plot_path)
    @test isfile(gmat_comparison_plot_path)

    pos_rmse = _scenario_rmse(summary, "cygnss_48hr_pvt")
    println("cygnss_48hr_pvt mean position-axis RMSE [km]: $(pos_rmse)")
    # Pin the published drag-free baseline (IEEE Aerospace 2026 Table 5 lineage,
    # reproduced during the July 2026 verification campaign to 3 mm). A drift
    # here means the reference builder's IC recipe or the force model changed.
    @test abs(pos_rmse - 1.5777551) < 0.005
    println("cygnss_48hr error plot: $(plot_path)")
    println("cygnss_48hr RTN error plot: $(plot_path_rtn)")
    println("cygnss_48hr orbital elements plot: $(plot_path_oe)")
    println("cygnss_48hr orbital element errors plot: $(plot_path_oe_err)")
    println("cygnss_48hr reference comparison plot: $(gmat_comparison_plot_path)")
    println("cygnss_48hr_pvt_position max error [km]: $(max(Float64.(summary[summary.event .== "state_x_time", :].max_abs_km[1]), Float64.(summary[summary.event .== "state_y_time", :].max_abs_km[1]), Float64.(summary[summary.event .== "state_z_time", :].max_abs_km[1])))")
    @test pos_rmse < 1.0e4
end

@testset "CYGNSS 96hr Orbit Simulation" begin
    result = _run_cygnss_cyg04_96hr_result_once()
    summary = result.summary
    errors = result.errors
    @test nrow(summary) == 4
    @test all(summary.scenario .== "cyg04_96hr_inertial")
    @test all(in.(summary.event, Ref(["altitude_time", "state_x_time", "state_y_time", "state_z_time"])))

    # Validate that the telemetry span is approximately 96 hours.
    span = _CYGNSS_CYG04_96HR_TIMESPAN_CACHE[]
    @test span !== nothing
    sim_span_s = span.tf_s - span.t0_s
    @test abs(sim_span_s - 96.0 * 3600.0) <= 10.0

    # Ensure the verification run used the full telemetry horizon.
    telemetry_end_s = maximum(Float64.(summary.telemetry_axis_end))
    @test abs(telemetry_end_s - span.tf_s) <= 1.0

    @test all(summary.n_sim .>= summary.min_eval_points)
    @test all(isfinite.(Float64.(summary.rmse_km)))
    @test all(isfinite.(Float64.(summary.max_abs_km)))

    plot_path = joinpath(_GMAT_REPO_ROOT, "output", "cygnss", "cygnss_96hr_error_timeseries.png")
    _plot_position_error_timeseries(errors, "cyg04_96hr_inertial", plot_path; title="CYGNSS 96hr Position Error Time Series")
    @test isfile(plot_path)

    plot_path_rtn = joinpath(_GMAT_REPO_ROOT, "output", "cygnss", "cygnss_96hr_error_timeseries_rtn.png")
    _plot_cyg04_96hr_error_timeseries_rtn(errors, plot_path_rtn)
    @test isfile(plot_path_rtn)

    plot_path_oe = joinpath(_GMAT_REPO_ROOT, "output", "cygnss", "cygnss_96hr_orbital_elements.png")
    plot_path_oe_err = joinpath(_GMAT_REPO_ROOT, "output", "cygnss", "cygnss_96hr_orbital_element_errors.png")
    _plot_cyg04_96hr_orbital_elements(errors, plot_path_oe, plot_path_oe_err)
    @test isfile(plot_path_oe)
    @test isfile(plot_path_oe_err)

    pos_rmse = _scenario_rmse(summary, "cyg04_96hr_inertial")
    println("cyg04_96hr_inertial mean position-axis RMSE [km]: $(pos_rmse)")
    println("cygnss_96hr error plot: $(plot_path)")
    println("cygnss_96hr RTN error plot: $(plot_path_rtn)")
    println("cygnss_96hr orbital elements plot: $(plot_path_oe)")
    println("cygnss_96hr orbital element errors plot: $(plot_path_oe_err)")
    println("cyg04_96hr_inertial_position max error [km]: $(max(Float64.(summary[summary.event .== "state_x_time", :].max_abs_km[1]), Float64.(summary[summary.event .== "state_y_time", :].max_abs_km[1]), Float64.(summary[summary.event .== "state_z_time", :].max_abs_km[1])))")
    @test pos_rmse < 1.0e4
end

@testset "CYGNSS GMAT CSV Comparison" begin
    result = _run_cygnss_gmat_csv_result_once()
    summary = result.summary
    errors = result.errors
    @test nrow(summary) == 4
    @test all(summary.scenario .== "cygnss_48hr_gmat_csv")
    @test all(in.(summary.event, Ref(["altitude_time", "state_x_time", "state_y_time", "state_z_time"])))

    span = _CYGNSS_GMAT_TIMESPAN_CACHE[]
    @test span !== nothing
    sim_span_s = span.tf_s - span.t0_s
    @test abs(sim_span_s - 48.0 * 3600.0) <= 1.0

    telemetry_end_s = maximum(Float64.(summary.telemetry_axis_end))
    @test abs(telemetry_end_s - span.tf_s) <= 1.0

    @test all(summary.n_sim .>= summary.min_eval_points)
    @test all(isfinite.(Float64.(summary.rmse_km)))
    @test all(isfinite.(Float64.(summary.max_abs_km)))

    plot_path = joinpath(_GMAT_REPO_ROOT, "output", "cygnss", "cygnss_48hr_gmat_position_error_timeseries.png")
    _plot_position_error_timeseries(
        errors,
        "cygnss_48hr_gmat_csv",
        plot_path;
        title="CYGNSS 48hr Position Error Time Series (GMAT vs SpaceAGORA)",
        time_unit=:hr
    )
    @test isfile(plot_path)
    comparison_plot_path = joinpath(_GMAT_REPO_ROOT, "output", "cygnss", "cygnss_48hr_reference_error_comparison.png")
    _plot_cygnss_reference_error_comparison(_run_cygnss_48hr_result_once().errors, errors, comparison_plot_path)
    @test isfile(comparison_plot_path)

    pos_rmse = _scenario_rmse(summary, "cygnss_48hr_gmat_csv")
    println("cygnss_48hr_gmat_csv mean position-axis RMSE [km]: $(pos_rmse)")
    println("cygnss_48hr_gmat_csv error plot: $(plot_path)")
    println("cygnss_48hr reference comparison plot: $(comparison_plot_path)")
    @test pos_rmse < 1.0e4
end

function _plot_cygnss_drag_force_timeseries(
    t_s::Vector{Float64},
    Fx::Vector{Float64},
    Fy::Vector{Float64},
    Fz::Vector{Float64},
    F_mag::Vector{Float64},
    inlier_mask::BitVector,
    outpath::String
)::String
    t_hr    = t_s ./ 3600.0
    Fx_mN   = Fx .* 1000.0
    Fy_mN   = Fy .* 1000.0
    Fz_mN   = Fz .* 1000.0
    Fmag_mN = F_mag .* 1000.0

    outlier_mask = .!inlier_mask

    # Panel 1: per-axis force components (inliers only for clarity)
    p_axes = plot(
        t_hr[inlier_mask], Fx_mN[inlier_mask];
        label="Fx", xlabel="Time (hr)", ylabel="Force (mN)",
        title="Inferred Drag Force — Components (IQR inliers)",
        lw=1.0, alpha=0.8, legend=:topright
    )
    plot!(p_axes, t_hr[inlier_mask], Fy_mN[inlier_mask]; label="Fy", lw=1.0, alpha=0.8)
    plot!(p_axes, t_hr[inlier_mask], Fz_mN[inlier_mask]; label="Fz", lw=1.0, alpha=0.8)

    # Panel 2: force magnitude with inlier/outlier split and statistics
    mean_mN = mean(Fmag_mN[inlier_mask])
    q1_mN   = quantile(Fmag_mN, 0.25)
    q3_mN   = quantile(Fmag_mN, 0.75)
    iqr_mN  = q3_mN - q1_mN
    lo_mN   = q1_mN - 1.5 * iqr_mN
    hi_mN   = q3_mN + 1.5 * iqr_mN

    p_mag = plot(
        t_hr[inlier_mask], Fmag_mN[inlier_mask];
        label="|F| inliers ($(sum(inlier_mask)))",
        xlabel="Time (hr)", ylabel="|Force| (mN)",
        title="Inferred Drag Force Magnitude",
        lw=0.8, alpha=0.7, color=:steelblue, legend=:topright
    )
    if any(outlier_mask)
        scatter!(
            p_mag, t_hr[outlier_mask], Fmag_mN[outlier_mask];
            label="outliers ($(sum(outlier_mask)))",
            markersize=3, alpha=0.5, color=:red, markerstrokewidth=0
        )
    end
    hline!(p_mag, [mean_mN]; label="mean $(round(mean_mN, digits=2)) mN",
           linestyle=:dash, color=:black, lw=1.5)
    hline!(p_mag, [lo_mN, hi_mN]; label="IQR ±1.5 bounds",
           linestyle=:dot, color=:orange, lw=1.2)

    fig = plot(p_axes, p_mag; layout=(2, 1), size=(1000, 700))
    mkpath(dirname(outpath))
    savefig(fig, outpath)
    return outpath
end

function _plot_cygnss_drag_tangential_timeseries(
    t_s::Vector{Float64},
    F_tang::Vector{Float64},
    inlier_mask::BitVector,
    outpath::String
)::String
    t_hr         = t_s ./ 3600.0
    Ftang_mN     = F_tang .* 1000.0
    outlier_mask = .!inlier_mask

    mean_mN = mean(Ftang_mN[inlier_mask])
    q1_mN   = quantile(Ftang_mN, 0.25)
    q3_mN   = quantile(Ftang_mN, 0.75)
    iqr_mN  = q3_mN - q1_mN
    lo_mN   = q1_mN - 1.5 * iqr_mN
    hi_mN   = q3_mN + 1.5 * iqr_mN

    p = plot(
        t_hr[inlier_mask], Ftang_mN[inlier_mask];
        label="F_tang inliers ($(sum(inlier_mask)))",
        xlabel="Time (hr)", ylabel="Tangential Force (mN)",
        title="Inferred Tangential Drag Force (along-track)\nnegative = opposing motion",
        lw=0.8, alpha=0.7, color=:steelblue, legend=:topright
    )
    if any(outlier_mask)
        scatter!(
            p, t_hr[outlier_mask], Ftang_mN[outlier_mask];
            label="outliers ($(sum(outlier_mask)))",
            markersize=3, alpha=0.5, color=:red, markerstrokewidth=0
        )
    end
    hline!(p, [mean_mN]; label="mean $(round(mean_mN, digits=2)) mN",
           linestyle=:dash, color=:black, lw=1.5)
    hline!(p, [lo_mN, hi_mN]; label="IQR ±1.5 bounds",
           linestyle=:dot, color=:orange, lw=1.2)
    hline!(p, [0.0]; label="zero", linestyle=:solid, color=:gray, lw=0.8, alpha=0.5)

    fig = plot(p; size=(1000, 450))
    mkpath(dirname(outpath))
    savefig(fig, outpath)
    return outpath
end

@testset "CYGNSS Drag Force Estimation" begin
    @test isfile(_CYGNSS_48HR_TELEMETRY_FEATHER)

    mktempdir() do tmp
        # --- 1. Build CYGNSS telemetry reference for ICs ---
        traj = _build_cygnss_48hr_reference(tmp, "cygnss_drag_est")

        # --- 2. Load raw GPS telemetry for per-sample velocity ---
        raw_df = DataFrame(Arrow.Table(_CYGNSS_48HR_TELEMETRY_FEATHER))
        t_col  = _required_column(raw_df, ["TIME OFFSET", "time"])
        x_col  = _required_column(raw_df, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_X (m)", "pos_ii_1"])
        y_col  = _required_column(raw_df, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_Y (m)", "pos_ii_2"])
        z_col  = _required_column(raw_df, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_Z (m)", "pos_ii_3"])
        vx_col = _required_column(raw_df, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_X (m/s)", "vel_ii_1"])
        vy_col = _required_column(raw_df, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_Y (m/s)", "vel_ii_2"])
        vz_col = _required_column(raw_df, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_Z (m/s)", "vel_ii_3"])

        perm   = sortperm(Float64.(raw_df[!, t_col]))
        t_tel  = Float64.(raw_df[perm, t_col])
        t_tel  = t_tel .- t_tel[1]    # relative seconds, same as traj
        vx_tel = Float64.(raw_df[perm, vx_col])    # m/s
        vy_tel = Float64.(raw_df[perm, vy_col])
        vz_tel = Float64.(raw_df[perm, vz_col])

        # --- 3. Build no-drag scenario: J50 + n-body Sun+Moon + SRP ---
        scenario = _base_scenario_dict("cygnss_drag_est", traj.telemetry_path)
        merge!(scenario, Dict{String, Any}(
            "planet"                   => "earth",
            "gravity_model"            => "inverse_squared",
            "gravity_harmonics_degree" => 50,
            "gravity_harmonics_order"  => 50,
            "gravity_harmonics_file"   => _GMAT_HARMONICS_EARTH_FILE,
            "nbody_bodies"             => Any["sun", "moon"],
            "srp_enabled"              => true,
            "srp_cr"                   => 1.3,
            "orbit_altitude_mode"      => "oblate",
            "drag_enabled"             => false,
            "include_wind"             => false,
            "max_points_quick"         => 200000,
            "max_points_full"          => 200000,
            "initial_time"             => Dict{String, Any}(
                "year" => 2025, "month" => 6, "day" => 6,
                "hour" => 0, "minute" => 0, "second" => 0.0
            )
        ))

        # --- 4. Run the no-drag simulation ---
        manifest_path = joinpath(tmp, "cygnss_drag_est_manifest.toml")
        open(manifest_path, "w") do io
            TOML.print(io, Dict{String, Any}("scenarios" => Any[scenario]))
        end

        scenarios_cfg = TV._load_scenarios_from_manifest(manifest_path)
        cfg           = only(scenarios_cfg)
        telemetry_cfg = TV._load_time_aligned_telemetry(cfg, 200000)
        ic            = TV._initial_condition_from_time_aligned_telemetry(cfg, telemetry_cfg)
        mission_s     = max(telemetry_cfg.time_s[end] - telemetry_cfg.time_s[1], 1.0)

        args = TV._make_time_aligned_args(cfg, mission_s, ic)
        args = TV._with_study_settings(args; quick=(TEST_MODE == :quick))

        gmat_kernel_env = Pair{String, String}(
            "SPACEAGORA_SPICE_PLANETARY_KERNEL_RELPATH",
            _gmat_planetary_kernel_relpath()
        )
        sim_run = withenv(gmat_kernel_env, pairs(_telemetry_solver_env_overrides())...) do
            TV._run_simulation_dataframe(args, cfg.name, cfg.atmosphere_truth, TEST_MODE)
        end
        sim_df = sim_run.results_df

        # --- 5. Extract simulation velocity time series ---
        t_sim  = Float64.(sim_df[!, _required_column(sim_df, ["time"])])
        vx_sim = Float64.(sim_df[!, _required_column(sim_df, ["sc1_vel_1", "sc1_velocity_1"])])
        vy_sim = Float64.(sim_df[!, _required_column(sim_df, ["sc1_vel_2", "sc1_velocity_2"])])
        vz_sim = Float64.(sim_df[!, _required_column(sim_df, ["sc1_vel_3", "sc1_velocity_3"])])

        @test length(t_sim) >= 10
        @test issorted(t_sim)

        # --- 6. Differentiate to get accelerations ---
        # Gravitational acceleration from no-drag sim (m/s²): dv_sim/dt
        ax_grav = TV._differentiate_series(vx_sim, t_sim)
        ay_grav = TV._differentiate_series(vy_sim, t_sim)
        az_grav = TV._differentiate_series(vz_sim, t_sim)

        # Total telemetry acceleration (m/s²): dv_tel/dt
        ax_tel_acc = TV._differentiate_series(vx_tel, t_tel)
        ay_tel_acc = TV._differentiate_series(vy_tel, t_tel)
        az_tel_acc = TV._differentiate_series(vz_tel, t_tel)

        # --- 7. Interpolate gravitational acceleration to telemetry sample times ---
        function _lininterp1(xq::Float64, x::Vector{Float64}, y::Vector{Float64})::Float64
            i = searchsortedfirst(x, xq)
            i <= 1 && return y[1]
            i > length(x) && return y[end]
            frac = (xq - x[i - 1]) / (x[i] - x[i - 1])
            return y[i - 1] * (1.0 - frac) + y[i] * frac
        end

        ax_grav_interp = [_lininterp1(t, t_sim, ax_grav) for t in t_tel]
        ay_grav_interp = [_lininterp1(t, t_sim, ay_grav) for t in t_tel]
        az_grav_interp = [_lininterp1(t, t_sim, az_grav) for t in t_tel]

        # --- 8. Compute inferred drag force ---
        spacecraft_mass_kg = 29.0    # CYGNSS bus_mass_kg from _base_scenario_dict
        Fx = spacecraft_mass_kg .* (ax_tel_acc .- ax_grav_interp)    # N
        Fy = spacecraft_mass_kg .* (ay_tel_acc .- ay_grav_interp)
        Fz = spacecraft_mass_kg .* (az_tel_acc .- az_grav_interp)
        F_mag = sqrt.(Fx .^ 2 .+ Fy .^ 2 .+ Fz .^ 2)    # N

        # --- 8b. Tangential (along-track) component: F · v̂ ---
        # Negative values mean the force opposes motion (drag direction)
        v_mag_tel = sqrt.(vx_tel .^ 2 .+ vy_tel .^ 2 .+ vz_tel .^ 2)
        vx_hat    = vx_tel ./ v_mag_tel
        vy_hat    = vy_tel ./ v_mag_tel
        vz_hat    = vz_tel ./ v_mag_tel
        F_tang    = Fx .* vx_hat .+ Fy .* vy_hat .+ Fz .* vz_hat    # N, signed

        # --- 9. IQR outlier filtering ---
        q1      = quantile(F_mag, 0.25)
        q3      = quantile(F_mag, 0.75)
        iqr_val = q3 - q1
        inlier_mask = (F_mag .>= q1 - 1.5 * iqr_val) .& (F_mag .<= q3 + 1.5 * iqr_val)
        F_filtered  = F_mag[inlier_mask]
        n_outliers  = count(.!inlier_mask)

        # --- 10. Assertions and reporting ---
        @test length(F_filtered) >= 2
        @test all(isfinite.(F_filtered))
        @test all(F_filtered .>= 0.0)

        mean_F_mN   = mean(F_filtered) * 1000.0
        median_F_mN = median(F_filtered) * 1000.0
        rms_F_mN    = sqrt(mean(F_filtered .^ 2)) * 1000.0
        max_F_mN    = maximum(F_filtered) * 1000.0

        println("CYGNSS inferred drag force (IQR-filtered, no-drag sim: J50+n-body+SRP):")
        println("  Mean:      $(round(mean_F_mN,   digits=3)) mN")
        println("  Median:    $(round(median_F_mN, digits=3)) mN")
        println("  RMS:       $(round(rms_F_mN,    digits=3)) mN")
        println("  Max:       $(round(max_F_mN,    digits=3)) mN")
        println("  Inliers:   $(length(F_filtered)) / $(length(F_mag))  ($(n_outliers) outliers removed)")

        # --- 10b. Tangential IQR filtering and reporting ---
        q1_tang     = quantile(F_tang, 0.25)
        q3_tang     = quantile(F_tang, 0.75)
        iqr_tang    = q3_tang - q1_tang
        inlier_tang = (F_tang .>= q1_tang - 1.5 * iqr_tang) .& (F_tang .<= q3_tang + 1.5 * iqr_tang)
        F_tang_filt = F_tang[inlier_tang]
        n_out_tang  = count(.!inlier_tang)

        @test length(F_tang_filt) >= 2
        @test all(isfinite.(F_tang_filt))

        mean_tang_mN   = mean(F_tang_filt) * 1000.0
        median_tang_mN = median(F_tang_filt) * 1000.0
        rms_tang_mN    = sqrt(mean(F_tang_filt .^ 2)) * 1000.0
        maxabs_tang_mN = maximum(abs.(F_tang_filt)) * 1000.0

        println("CYGNSS inferred tangential force (IQR-filtered, negative = drag):")
        println("  Mean:      $(round(mean_tang_mN,   digits=3)) mN")
        println("  Median:    $(round(median_tang_mN, digits=3)) mN")
        println("  RMS:       $(round(rms_tang_mN,    digits=3)) mN")
        println("  |Max|:     $(round(maxabs_tang_mN, digits=3)) mN")
        println("  Inliers:   $(length(F_tang_filt)) / $(length(F_tang))  ($(n_out_tang) outliers removed)")

        # --- 11. Plot drag force timeseries ---
        drag_plot_path = joinpath(_GMAT_REPO_ROOT, "output", "cygnss", "cygnss_drag_force_timeseries.png")
        _plot_cygnss_drag_force_timeseries(t_tel, Fx, Fy, Fz, F_mag, inlier_mask, drag_plot_path)
        @test isfile(drag_plot_path)
        println("cygnss drag force timeseries plot: $(drag_plot_path)")

        # --- 12. Plot tangential component timeseries ---
        tang_plot_path = joinpath(_GMAT_REPO_ROOT, "output", "cygnss", "cygnss_drag_tangential_timeseries.png")
        _plot_cygnss_drag_tangential_timeseries(t_tel, F_tang, inlier_tang, tang_plot_path)
        @test isfile(tang_plot_path)
        println("cygnss drag tangential timeseries plot: $(tang_plot_path)")
    end
end

end # !_cygnss_private_data_available() guard around the CYGNSS testsets

# ---------------------------------------------------------------------------
# SpaceAGORA Examples Export
# ---------------------------------------------------------------------------

function _export_spaceagora_examples(result::TV.VerificationResult, outdir::String)
    mkpath(outdir)

    for scenario_name in sort!(collect(_active_gmat_expected_scenario_names()))
        errors = result.errors

        # Filter to x/y/z position events for this scenario, sorted by idx
        rows = errors[(errors.scenario .== scenario_name) .& in.(errors.event, Ref(["state_x_time", "state_y_time", "state_z_time"])), :]
        xrows = sort(rows[rows.event .== "state_x_time", :], :idx)
        yrows = sort(rows[rows.event .== "state_y_time", :], :idx)
        zrows = sort(rows[rows.event .== "state_z_time", :], :idx)

        n = min(nrow(xrows), nrow(yrows), nrow(zrows))
        n >= 2 || continue

        t_s  = Float64.(xrows.telemetry_axis[1:n])
        x_km = Float64.(xrows.sim_interp_value_km[1:n])
        y_km = Float64.(yrows.sim_interp_value_km[1:n])
        z_km = Float64.(zrows.sim_interp_value_km[1:n])

        vx_km_s = _differentiate_position_series(x_km, t_s)
        vy_km_s = _differentiate_position_series(y_km, t_s)
        vz_km_s = _differentiate_position_series(z_km, t_s)

        planet = TV._planet_from_name(_scenario_planet_name(scenario_name))

        sma_km = Vector{Float64}(undef, n)
        ecc    = Vector{Float64}(undef, n)
        inc_deg = Vector{Float64}(undef, n)

        for i in 1:n
            r = SVector{3, Float64}(x_km[i], y_km[i], z_km[i]) .* 1.0e3
            v = SVector{3, Float64}(vx_km_s[i], vy_km_s[i], vz_km_s[i]) .* 1.0e3
            try
                oe = TV.rvtoorbitalelement(r, v, planet)
                sma_km[i]  = oe[1] * 1.0e-3
                ecc[i]     = oe[2]
                inc_deg[i] = rad2deg(oe[3])
            catch
                sma_km[i]  = NaN
                ecc[i]     = NaN
                inc_deg[i] = NaN
            end
        end

        df = DataFrame(
            ElapsedSecs = t_s,
            SMA         = sma_km,
            ECC         = ecc,
            INC         = inc_deg,
            X           = x_km,
            Y           = y_km,
            Z           = z_km,
            VX          = vx_km_s,
            VY          = vy_km_s,
            VZ          = vz_km_s,
        )

        fname = _scenario_basilisk_file_name(scenario_name)
        Arrow.write(joinpath(outdir, fname), df)
    end
end

if !_parse_bool_env("SPACEAGORA_SKIP_GMAT_MATRIX", false)

@testset "SpaceAGORA Examples Export" begin
    result = _run_gmat_scenario_matrix_result_once()
    outdir = joinpath(_GMAT_REPO_ROOT, "data", "telemetry", "SpaceAGORA_Examples")
    _export_spaceagora_examples(result, outdir)
    for scenario_name in sort!(collect(_active_gmat_expected_scenario_names()))
        fname = _scenario_basilisk_file_name(scenario_name)
        fpath = joinpath(outdir, fname)
        @test isfile(fpath)
        println("SpaceAGORA exported: $fpath")
    end
end

end # SPACEAGORA_SKIP_GMAT_MATRIX (SpaceAGORA Examples Export)
