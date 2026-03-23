using Test
using TOML
using Arrow
using CSV
using DataFrames
using Statistics
using LinearAlgebra
using Plots
using StaticArrays

const _GMAT_REPO_ROOT = isdefined(Main, :REPO_ROOT) ? Main.REPO_ROOT : normpath(joinpath(@__DIR__, ".."))

if !isdefined(@__MODULE__, :TelemetryVerification)
    include(joinpath(_GMAT_REPO_ROOT, "src", "analysis", "verification", "TelemetryVerification.jl"))
end

if !isdefined(@__MODULE__, :SimulationModel)
    include(joinpath(_GMAT_REPO_ROOT, "src", "simulation_model", "SimulationModel.jl"))
end

const TV = TelemetryVerification
const SM = SimulationModel

const _GMAT_EXAMPLES_DIR = joinpath(
    _GMAT_REPO_ROOT,
    "data",
    "telemetry",
    "GMAT_Examples"
)

const _GMAT_HARMONICS_EARTH_FILE = "data/Gravity_harmonics_data/egm96.csv"
const _GMAT_HARMONICS_MARS_FILE = "data/Gravity_harmonics_data/Mars50c.csv"
const _GMAT_HARMONICS_VENUS_FILE = "data/Gravity_harmonics_data/MGNP180U.csv"
const _CYGNSS_48HR_TELEMETRY_FEATHER = joinpath(_GMAT_REPO_ROOT, "data", "telemetry", "CYGNSS", "cygnss_data_48hr.feather")
const _CYGNSS_GMAT_COMPARISON_PATH = joinpath(_GMAT_EXAMPLES_DIR, "Sim_CYGNSS_Comparison.feather")

# ============================================================================
# TEST MODE CONFIGURATION
# ============================================================================
# Change this to :full for comprehensive testing with more evaluation points
# Change to :quick for faster testing with fewer evaluation points
const TEST_MODE::Symbol = :quick
# ============================================================================

@inline function _verification_profile_from_env()::Symbol
    return TEST_MODE
end

@inline function _strict_position_rmse_limit_km(scenario_name::String, profile::Symbol)::Float64
    if startswith(scenario_name, "earth_j50_")
        return profile == :full ? 0.45 : 0.05
    elseif startswith(scenario_name, "earth_j2_")
        return profile == :full ? 1e-2 : 1e-3
    elseif startswith(scenario_name, "earth_j0_")
        return profile == :full ? 1e-2 : 1e-3
    elseif startswith(scenario_name, "earth_j2_tbfalse_manual_copy")
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

    time_col = _required_column(df, ["Sat.ElapsedSecs", "DefaultSC.ElapsedSecs"])
    x_col = _required_column(df, ["Sat.X", "Sat.PlanetInertial.X", "DefaultSC.EarthMJ2000Eq.X"])
    y_col = _required_column(df, ["Sat.Y", "Sat.PlanetInertial.Y", "DefaultSC.EarthMJ2000Eq.Y"])
    z_col = _required_column(df, ["Sat.Z", "Sat.PlanetInertial.Z", "DefaultSC.EarthMJ2000Eq.Z"])

    t = Float64.(df[!, time_col])
    x = Float64.(df[!, x_col])
    y = Float64.(df[!, y_col])
    z = Float64.(df[!, z_col])
    planet = TV._planet_from_name(planet_name)

    perm = sortperm(t)
    t = t[perm]
    x = x[perm]
    y = y[perm]
    z = z[perm]

    t_rel = t .- t[1]

    # # GMAT exports in this folder can mix km and m positions across scenarios.
    # # Detect oversized radii relative to the planet and convert m -> km when needed.
    # rp_km = planet.Rp_e * 1e-3
    # r_raw = sqrt.(x .^ 2 .+ y .^ 2 .+ z .^ 2)
    # if median(r_raw) > 100.0 * rp_km
    #     x .*= 1.0e-3
    #     y .*= 1.0e-3
    #     z .*= 1.0e-3
    # end

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
        sma_km=fill(sma_km, length(t_rel)),
        ecc=fill(ecc, length(t_rel)),
        inc_deg=fill(inc_deg, length(t_rel)),
        aop_deg=fill(aop_deg, length(t_rel)),
        raan_deg=fill(raan_deg, length(t_rel)),
        ta_deg=fill(ta_deg, length(t_rel))
    )

    telemetry_path = joinpath(outdir, "$(stem)_time_aligned.arrow")
    Arrow.write(telemetry_path, telemetry_df)

    return (
        telemetry_path=telemetry_path,
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

    t = Float64.(df[!, time_col])
    x_km = Float64.(df[!, x_col]) .* 1.0e-3
    y_km = Float64.(df[!, y_col]) .* 1.0e-3
    z_km = Float64.(df[!, z_col]) .* 1.0e-3

    perm = sortperm(t)
    t = t[perm]
    x_km = x_km[perm]
    y_km = y_km[perm]
    z_km = z_km[perm]

    t_rel = t .- t[1]
    planet = TV._planet_from_name("earth")
    r_km = sqrt.(x_km .^ 2 .+ y_km .^ 2 .+ z_km .^ 2)
    alt_km = r_km .- planet.Rp_e * 1.0e-3

    # User-provided CYGNSS initial orbital elements.
    sma_km = 6818.860956945965
    ecc = 0.0004790007068221837
    # ecc = 1e-5
    inc_deg = 34.93573031327011
    aop_deg = 140.63502513151334
    raan_deg = 177.37117999662348
    ta_deg = 276.65207622217173

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
        ta_deg=fill(ta_deg, length(t_rel))
    )

    telemetry_path = joinpath(outdir, "$(stem)_time_aligned.arrow")
    Arrow.write(telemetry_path, telemetry_df)

    return (
        telemetry_path=telemetry_path,
        t0_s=t_rel[1],
        tf_s=t_rel[end],
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
            "sma" => "sma_km",
            "ecc" => "ecc",
            "inc" => "inc_deg",
            "aop" => "aop_deg",
            "raan" => "raan_deg",
            "ta" => "ta_deg"
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
            # "bus_mass_kg" => 620.0,
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

const _GMAT_MATRIX_SUMMARY_CACHE = Ref{Union{Nothing, DataFrame}}(nothing)
const _GMAT_MATRIX_RESULT_CACHE = Ref{Union{Nothing, TV.VerificationResult}}(nothing)
const _CYGNSS_48HR_SUMMARY_CACHE = Ref{Union{Nothing, DataFrame}}(nothing)
const _CYGNSS_48HR_RESULT_CACHE = Ref{Union{Nothing, TV.VerificationResult}}(nothing)
const _CYGNSS_48HR_TIMESPAN_CACHE = Ref{Union{Nothing, NamedTuple{(:t0_s, :tf_s), Tuple{Float64, Float64}}}}(nothing)
const _CYGNSS_GMAT_SUMMARY_CACHE = Ref{Union{Nothing, DataFrame}}(nothing)
const _CYGNSS_GMAT_RESULT_CACHE = Ref{Union{Nothing, TV.VerificationResult}}(nothing)
const _CYGNSS_GMAT_TIMESPAN_CACHE = Ref{Union{Nothing, NamedTuple{(:t0_s, :tf_s), Tuple{Float64, Float64}}}}(nothing)

function _run_gmat_scenario_matrix_result_once()::TV.VerificationResult
    if _GMAT_MATRIX_RESULT_CACHE[] !== nothing
        return _GMAT_MATRIX_RESULT_CACHE[]
    end

    earth_j2_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Earth_J2_TBFalse.feather")
    earth_j0_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Earth_J0_TBFalse.feather")
    earth_j50_tbfalse_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Earth_J50_TBFalse.feather")
    earth_j50_tbtrue_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Earth_J50_TBTrue.feather")
    earth_j2_tbfalse_manual_copy_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Earth_J2_TBFalse_manual copy.feather")
    earth_j2_tbtrue_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Earth_J2_TBTrue.feather")
    earth_j0_tbtrue_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Earth_J0_TBTrue.feather")
    mars_j2_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Mars_J2_TBFalse.feather")
    mars_j0_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Mars_J0_TBFalse.feather")
    mars_j2_tbtrue_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Mars_J2_TBTrue.feather")
    mars_j0_tbtrue_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Mars_J0_TBTrue.feather")
    mars_j50_tbfalse_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Mars_J50_TBFalse.feather")
    mars_j50_tbtrue_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Mars_J50_TBTrue.feather")
    venus_j2_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Venus_J2_TBFalse.feather")
    venus_j0_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Venus_J0_TBFalse.feather")
    venus_j2_tbtrue_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Venus_J2_TBTrue.feather")
    venus_j0_tbtrue_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Venus_J0_TBTrue.feather")
    venus_j50_tbfalse_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Venus_J50_TBFalse.feather")
    venus_j50_tbtrue_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Venus_J50_TBTrue.feather")
    @test isfile(earth_j2_csv)
    @test isfile(earth_j0_csv)
    @test isfile(earth_j50_tbfalse_csv)
    @test isfile(earth_j50_tbtrue_csv)
    @test isfile(earth_j2_tbfalse_manual_copy_csv)
    @test isfile(earth_j2_tbtrue_csv)
    @test isfile(earth_j0_tbtrue_csv)
    @test isfile(mars_j2_csv)
    @test isfile(mars_j0_csv)
    @test isfile(mars_j2_tbtrue_csv)
    @test isfile(mars_j0_tbtrue_csv)
    @test isfile(mars_j50_tbfalse_csv)
    @test isfile(mars_j50_tbtrue_csv)
    @test isfile(venus_j2_csv)
    @test isfile(venus_j0_csv)
    @test isfile(venus_j2_tbtrue_csv)
    @test isfile(venus_j0_tbtrue_csv)
    @test isfile(venus_j50_tbfalse_csv)
    @test isfile(venus_j50_tbtrue_csv)
    @test isfile(joinpath(_GMAT_REPO_ROOT, _GMAT_HARMONICS_EARTH_FILE))
    @test isfile(joinpath(_GMAT_REPO_ROOT, _GMAT_HARMONICS_MARS_FILE))
    @test isfile(joinpath(_GMAT_REPO_ROOT, _GMAT_HARMONICS_VENUS_FILE))

    result = mktempdir() do tmp
        traj_j2 = _build_time_aligned_reference(earth_j2_csv, "earth", tmp, "earth_j2_tbfalse")
        traj_j2_ic = _build_time_aligned_reference(
            earth_j2_csv,
            "earth",
            tmp,
            "earth_j2_tbfalse_ic_perturbed";
            raan_offset_deg=35.0,
            ta_offset_deg=110.0
        )
        traj_j0 = _build_time_aligned_reference(earth_j0_csv, "earth", tmp, "earth_j0_tbfalse")
        traj_j50_tbfalse = _build_time_aligned_reference(earth_j50_tbfalse_csv, "earth", tmp, "earth_j50_tbfalse")
        traj_j50_tbtrue = _build_time_aligned_reference(earth_j50_tbtrue_csv, "earth", tmp, "earth_j50_tbtrue")
        traj_j2_tbfalse_manual_copy = _build_time_aligned_reference(earth_j2_tbfalse_manual_copy_csv, "earth", tmp, "earth_j2_tbfalse_manual_copy")
        traj_j2_tbtrue = _build_time_aligned_reference(earth_j2_tbtrue_csv, "earth", tmp, "earth_j2_tbtrue")
        traj_j0_tbtrue = _build_time_aligned_reference(earth_j0_tbtrue_csv, "earth", tmp, "earth_j0_tbtrue")
        traj_mars_j2_tbfalse = _build_time_aligned_reference(
            mars_j2_csv,
            "mars",
            tmp,
            "mars_j2_tbfalse";
            sma_km=3700.0,
            ecc=0.002,
            inc_deg=45.0,
            aop_deg=0.0,
            raan_deg=0.0,
            ta_deg=90.0
        )
        traj_mars_j0_tbfalse = _build_time_aligned_reference(
            mars_j0_csv,
            "mars",
            tmp,
            "mars_j0_tbfalse";
            sma_km=3700.0,
            ecc=0.002,
            inc_deg=45.0,
            aop_deg=0.0,
            raan_deg=0.0,
            ta_deg=90.0
        )
        traj_mars_j2_tbtrue = _build_time_aligned_reference(
            mars_j2_tbtrue_csv,
            "mars",
            tmp,
            "mars_j2_tbtrue";
            sma_km=3700.0,
            ecc=0.002,
            inc_deg=45.0,
            aop_deg=0.0,
            raan_deg=0.0,
            ta_deg=90.0
        )
        traj_mars_j0_tbtrue = _build_time_aligned_reference(
            mars_j0_tbtrue_csv,
            "mars",
            tmp,
            "mars_j0_tbtrue";
            sma_km=3700.0,
            ecc=0.002,
            inc_deg=45.0,
            aop_deg=0.0,
            raan_deg=0.0,
            ta_deg=90.0
        )
        traj_mars_j50_tbfalse = _build_time_aligned_reference(
            mars_j50_tbfalse_csv,
            "mars",
            tmp,
            "mars_j50_tbfalse";
            sma_km=3700.0,
            ecc=0.002,
            inc_deg=45.0,
            aop_deg=0.0,
            raan_deg=0.0,
            ta_deg=90.0
        )
        traj_mars_j50_tbtrue = _build_time_aligned_reference(
            mars_j50_tbtrue_csv,
            "mars",
            tmp,
            "mars_j50_tbtrue";
            sma_km=3700.0,
            ecc=0.002,
            inc_deg=45.0,
            aop_deg=0.0,
            raan_deg=0.0,
            ta_deg=90.0
        )
        traj_venus_j2_tbfalse = _build_time_aligned_reference(
            venus_j2_csv,
            "venus",
            tmp,
            "venus_j2_tbfalse";
            sma_km=6000.0,
            ecc=0.001,
            inc_deg=30.0,
            aop_deg=0.0,
            raan_deg=0.0,
            ta_deg=180.0
        )
        traj_venus_j0_tbfalse = _build_time_aligned_reference(
            venus_j0_csv,
            "venus",
            tmp,
            "venus_j0_tbfalse";
            sma_km=6000.0,
            ecc=0.001,
            inc_deg=30.0,
            aop_deg=0.0,
            raan_deg=0.0,
            ta_deg=180.0
        )
        traj_venus_j2_tbtrue = _build_time_aligned_reference(
            venus_j2_tbtrue_csv,
            "venus",
            tmp,
            "venus_j2_tbtrue";
            sma_km=6000.0,
            ecc=0.001,
            inc_deg=30.0,
            aop_deg=0.0,
            raan_deg=0.0,
            ta_deg=180.0
        )
        traj_venus_j0_tbtrue = _build_time_aligned_reference(
            venus_j0_tbtrue_csv,
            "venus",
            tmp,
            "venus_j0_tbtrue";
            sma_km=6000.0,
            ecc=0.001,
            inc_deg=30.0,
            aop_deg=0.0,
            raan_deg=0.0,
            ta_deg=180.0
        )
        traj_venus_j50_tbfalse = _build_time_aligned_reference(
            venus_j50_tbfalse_csv,
            "venus",
            tmp,
            "venus_j50_tbfalse";
            sma_km=6000.0,
            ecc=0.001,
            inc_deg=30.0,
            aop_deg=0.0,
            raan_deg=0.0,
            ta_deg=180.0
        )
        traj_venus_j50_tbtrue = _build_time_aligned_reference(
            venus_j50_tbtrue_csv,
            "venus",
            tmp,
            "venus_j50_tbtrue";
            sma_km=6000.0,
            ecc=0.001,
            inc_deg=30.0,
            aop_deg=0.0,
            raan_deg=0.0,
            ta_deg=180.0
        )

        println(
            "SpaceAGORA initial orbital elements [a_km,e,i_deg,aop_deg,raan_deg,ta_deg]: " *
            "earth_j2_tbfalse=($(traj_j2.sma_km),$(traj_j2.ecc),$(traj_j2.inc_deg),$(traj_j2.aop_deg),$(traj_j2.raan_deg),$(traj_j2.ta_deg)); " *
            # "earth_j2_tbfalse_ic_perturbed=($(traj_j2_ic.sma_km),$(traj_j2_ic.ecc),$(traj_j2_ic.inc_deg),$(traj_j2_ic.aop_deg),$(traj_j2_ic.raan_deg),$(traj_j2_ic.ta_deg)); " *
            "earth_j0_tbfalse=($(traj_j0.sma_km),$(traj_j0.ecc),$(traj_j0.inc_deg),$(traj_j0.aop_deg),$(traj_j0.raan_deg),$(traj_j0.ta_deg)); " *
            "earth_j50_tbfalse=($(traj_j50_tbfalse.sma_km),$(traj_j50_tbfalse.ecc),$(traj_j50_tbfalse.inc_deg),$(traj_j50_tbfalse.aop_deg),$(traj_j50_tbfalse.raan_deg),$(traj_j50_tbfalse.ta_deg)); " *
            "earth_j2_tbtrue=($(traj_j2_tbtrue.sma_km),$(traj_j2_tbtrue.ecc),$(traj_j2_tbtrue.inc_deg),$(traj_j2_tbtrue.aop_deg),$(traj_j2_tbtrue.raan_deg),$(traj_j2_tbtrue.ta_deg)); " *
            "earth_j0_tbtrue=($(traj_j0_tbtrue.sma_km),$(traj_j0_tbtrue.ecc),$(traj_j0_tbtrue.inc_deg),$(traj_j0_tbtrue.aop_deg),$(traj_j0_tbtrue.raan_deg),$(traj_j0_tbtrue.ta_deg))"
        )

        baseline = _base_scenario_dict("earth_j0_tbfalse", traj_j0.telemetry_path)
        merge!(baseline, Dict{String, Any}(
            "planet" => "earth",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 0,
            "gravity_harmonics_order" => 0,
            "gravity_harmonics_file" => _GMAT_HARMONICS_EARTH_FILE,
            "nbody_bodies" => Any[],
            "orbit_altitude_mode" => "oblate"
        ))

        j2_tbfalse = _base_scenario_dict("earth_j2_tbfalse", traj_j2.telemetry_path)
        merge!(j2_tbfalse, Dict{String, Any}(
            "planet" => "earth",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 2,
            "gravity_harmonics_order" => 0,
            "gravity_harmonics_file" => _GMAT_HARMONICS_EARTH_FILE,
            "nbody_bodies" => Any[],
            "orbit_altitude_mode" => "oblate"
        ))

        # ic_perturbed = deepcopy(baseline)
        # ic_perturbed["name"] = "earth_j2_tbfalse_ic_perturbed"
        # ic_perturbed["telemetry"] = traj_j2_ic.telemetry_path

        # force_variant = _base_scenario_dict("earth_j0_tbfalse", traj_j0.telemetry_path)
        # merge!(force_variant, Dict{String, Any}(
        #     "planet" => "earth",
        #     "gravity_model" => "inverse_squared",
        #     "gravity_harmonics_degree" => 0,
        #     "gravity_harmonics_order" => 0,
        #     "gravity_harmonics_file" => "",
        #     "nbody_bodies" => Any[],
        #     "orbit_altitude_mode" => "oblate",
        #     "srp_enabled" => false, # intentionally disable SRP for this case to isolate the effect of third bodies vs SRP in the j2_tbtrue case
        #     "srp_area_m2" => 5.0,
        #     "srp_cr" => 1.35
        # ))

        j2_tbtrue = _base_scenario_dict("earth_j2_tbtrue", traj_j2_tbtrue.telemetry_path)
        merge!(j2_tbtrue, Dict{String, Any}(
            "planet" => "earth",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 2,
            "gravity_harmonics_order" => 0,
            "gravity_harmonics_file" => _GMAT_HARMONICS_EARTH_FILE,
            "nbody_bodies" => Any["sun", "moon"],
            "orbit_altitude_mode" => "oblate"
        ))

        j50_tbfalse = _base_scenario_dict("earth_j50_tbfalse", traj_j50_tbfalse.telemetry_path)
        merge!(j50_tbfalse, Dict{String, Any}(
            "planet" => "earth",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 50,
            "gravity_harmonics_order" => 50,
            "gravity_harmonics_file" => _GMAT_HARMONICS_EARTH_FILE,
            "nbody_bodies" => Any[],
            "orbit_altitude_mode" => "oblate"
        ))

        j50_tbtrue = _base_scenario_dict("earth_j50_tbtrue", traj_j50_tbtrue.telemetry_path)
        merge!(j50_tbtrue, Dict{String, Any}(
            "planet" => "earth",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 50,
            "gravity_harmonics_order" => 50,
            "gravity_harmonics_file" => _GMAT_HARMONICS_EARTH_FILE,
            "nbody_bodies" => Any["sun", "moon"],
            "orbit_altitude_mode" => "oblate"
        ))

        j2_tbfalse_manual_copy = _base_scenario_dict("earth_j2_tbfalse_manual_copy", traj_j2_tbfalse_manual_copy.telemetry_path)
        merge!(j2_tbfalse_manual_copy, Dict{String, Any}(
            "planet" => "earth",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 2,
            "gravity_harmonics_order" => 0,
            "gravity_harmonics_file" => _GMAT_HARMONICS_EARTH_FILE,
            "nbody_bodies" => Any[],
            "orbit_altitude_mode" => "oblate"
        ))

        j0_tbtrue = _base_scenario_dict("earth_j0_tbtrue", traj_j0_tbtrue.telemetry_path)
        merge!(j0_tbtrue, Dict{String, Any}(
            "planet" => "earth",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 0,
            "gravity_harmonics_order" => 0,
            "gravity_harmonics_file" => _GMAT_HARMONICS_EARTH_FILE,
            "nbody_bodies" => Any["sun", "moon"],
            "orbit_altitude_mode" => "oblate",
            "srp_enabled" => false, # intentionally disable SRP for this case to isolate the effect of third bodies vs SRP in the j2_tbtrue case
            "srp_area_m2" => 5.0,
            "srp_cr" => 1.35
        ))

        mars_j0_tbfalse = _base_scenario_dict("mars_j0_tbfalse", traj_mars_j0_tbfalse.telemetry_path)
        merge!(mars_j0_tbfalse, Dict{String, Any}(
            "planet" => "mars",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 0,
            "gravity_harmonics_order" => 0,
            "gravity_harmonics_file" => _GMAT_HARMONICS_MARS_FILE,
            "nbody_bodies" => Any[],
            "orbit_altitude_mode" => "oblate"
        ))

        mars_j2_tbfalse = _base_scenario_dict("mars_j2_tbfalse", traj_mars_j2_tbfalse.telemetry_path)
        merge!(mars_j2_tbfalse, Dict{String, Any}(
            "planet" => "mars",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 2,
            "gravity_harmonics_order" => 0,
            "gravity_harmonics_file" => _GMAT_HARMONICS_MARS_FILE,
            "nbody_bodies" => Any[],
            "orbit_altitude_mode" => "oblate"
        ))

        mars_j0_tbtrue = _base_scenario_dict("mars_j0_tbtrue", traj_mars_j0_tbtrue.telemetry_path)
        merge!(mars_j0_tbtrue, Dict{String, Any}(
            "planet" => "mars",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 0,
            "gravity_harmonics_order" => 0,
            "gravity_harmonics_file" => _GMAT_HARMONICS_MARS_FILE,
            "nbody_bodies" => Any["sun"],
            "orbit_altitude_mode" => "oblate"
        ))

        mars_j2_tbtrue = _base_scenario_dict("mars_j2_tbtrue", traj_mars_j2_tbtrue.telemetry_path)
        merge!(mars_j2_tbtrue, Dict{String, Any}(
            "planet" => "mars",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 2,
            "gravity_harmonics_order" => 0,
            "gravity_harmonics_file" => _GMAT_HARMONICS_MARS_FILE,
            "nbody_bodies" => Any["sun"],
            "orbit_altitude_mode" => "oblate"
        ))

        venus_j0_tbfalse = _base_scenario_dict("venus_j0_tbfalse", traj_venus_j0_tbfalse.telemetry_path)
        merge!(venus_j0_tbfalse, Dict{String, Any}(
            "planet" => "venus",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 0,
            "gravity_harmonics_order" => 0,
            "gravity_harmonics_file" => _GMAT_HARMONICS_VENUS_FILE,
            "nbody_bodies" => Any[],
            "orbit_altitude_mode" => "oblate"
        ))

        venus_j2_tbfalse = _base_scenario_dict("venus_j2_tbfalse", traj_venus_j2_tbfalse.telemetry_path)
        merge!(venus_j2_tbfalse, Dict{String, Any}(
            "planet" => "venus",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 2,
            "gravity_harmonics_order" => 0,
            "gravity_harmonics_file" => _GMAT_HARMONICS_VENUS_FILE,
            "nbody_bodies" => Any[],
            "orbit_altitude_mode" => "oblate"
        ))

        venus_j0_tbtrue = _base_scenario_dict("venus_j0_tbtrue", traj_venus_j0_tbtrue.telemetry_path)
        merge!(venus_j0_tbtrue, Dict{String, Any}(
            "planet" => "venus",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 0,
            "gravity_harmonics_order" => 0,
            "gravity_harmonics_file" => _GMAT_HARMONICS_VENUS_FILE,
            "nbody_bodies" => Any["sun"],
            "orbit_altitude_mode" => "oblate"
        ))

        venus_j2_tbtrue = _base_scenario_dict("venus_j2_tbtrue", traj_venus_j2_tbtrue.telemetry_path)
        merge!(venus_j2_tbtrue, Dict{String, Any}(
            "planet" => "venus",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 2,
            "gravity_harmonics_order" => 0,
            "gravity_harmonics_file" => _GMAT_HARMONICS_VENUS_FILE,
            "nbody_bodies" => Any["sun"],
            "orbit_altitude_mode" => "oblate"
        ))

        mars_j50_tbfalse = _base_scenario_dict("mars_j50_tbfalse", traj_mars_j50_tbfalse.telemetry_path)
        merge!(mars_j50_tbfalse, Dict{String, Any}(
            "planet" => "mars",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 50,
            "gravity_harmonics_order" => 50,
            "gravity_harmonics_file" => _GMAT_HARMONICS_MARS_FILE,
            "nbody_bodies" => Any[],
            "orbit_altitude_mode" => "oblate"
        ))

        mars_j50_tbtrue = _base_scenario_dict("mars_j50_tbtrue", traj_mars_j50_tbtrue.telemetry_path)
        merge!(mars_j50_tbtrue, Dict{String, Any}(
            "planet" => "mars",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 50,
            "gravity_harmonics_order" => 50,
            "gravity_harmonics_file" => _GMAT_HARMONICS_MARS_FILE,
            "nbody_bodies" => Any["sun"],
            "orbit_altitude_mode" => "oblate"
        ))

        venus_j50_tbfalse = _base_scenario_dict("venus_j50_tbfalse", traj_venus_j50_tbfalse.telemetry_path)
        merge!(venus_j50_tbfalse, Dict{String, Any}(
            "planet" => "venus",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 50,
            "gravity_harmonics_order" => 50,
            "gravity_harmonics_file" => _GMAT_HARMONICS_VENUS_FILE,
            "nbody_bodies" => Any[],
            "orbit_altitude_mode" => "oblate"
        ))

        venus_j50_tbtrue = _base_scenario_dict("venus_j50_tbtrue", traj_venus_j50_tbtrue.telemetry_path)
        merge!(venus_j50_tbtrue, Dict{String, Any}(
            "planet" => "venus",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 50,
            "gravity_harmonics_order" => 50,
            "gravity_harmonics_file" => _GMAT_HARMONICS_VENUS_FILE,
            "nbody_bodies" => Any["sun"],
            "orbit_altitude_mode" => "oblate"
        ))

        manifest_path = joinpath(tmp, "gmat_scenario_matrix_manifest.toml")
        open(manifest_path, "w") do io
            TOML.print(io, Dict{String, Any}("scenarios" => Any[
                baseline,
                j2_tbfalse,
                j2_tbtrue,
                j0_tbtrue,
                j50_tbfalse,
                j50_tbtrue,
                j2_tbfalse_manual_copy,
                mars_j0_tbfalse,
                mars_j2_tbfalse,
                mars_j0_tbtrue,
                mars_j2_tbtrue,
                mars_j50_tbfalse,
                mars_j50_tbtrue,
                venus_j0_tbfalse,
                venus_j2_tbfalse,
                venus_j0_tbtrue,
                venus_j2_tbtrue,
                venus_j50_tbfalse,
                venus_j50_tbtrue
            ]))
        end

        req = TV.VerificationRequest(
            profile=_verification_profile_from_env(),
            out_summary=joinpath(tmp, "summary.csv"),
            out_errors=joinpath(tmp, "errors.csv"),
            manifest_path=manifest_path,
            enforce=false,
            generate_plots=false
        )
        return withenv(
            "SPACEAGORA_TELEMETRY_SOLVER_MODE" => "tsit5",
            "SPACEAGORA_TELEMETRY_DT_MAX_ORBIT" => "10.0",
            "SPACEAGORA_TELEMETRY_RELTOL_ORBIT" => "1e-12",
            "SPACEAGORA_TELEMETRY_ABSTOL_ORBIT" => "1e-12",
            "SPACEAGORA_TELEMETRY_RELTOL_ATM" => "1e-12",
            "SPACEAGORA_TELEMETRY_ABSTOL_ATM" => "1e-12"
        ) do
            TV.run_verification(req)
        end
    end

    _GMAT_MATRIX_RESULT_CACHE[] = result
    _GMAT_MATRIX_SUMMARY_CACHE[] = result.summary
    return result
end

function _run_gmat_scenario_matrix_once()::DataFrame
    return _run_gmat_scenario_matrix_result_once().summary
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
            "drag_enabled" => true,
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
            profile=_verification_profile_from_env(),
            out_summary=joinpath(tmp, "summary.csv"),
            out_errors=joinpath(tmp, "errors.csv"),
            manifest_path=manifest_path,
            enforce=false,
            generate_plots=false
        )

        return withenv(
            "SPACEAGORA_TELEMETRY_SOLVER_MODE" => "tsit5",
            "SPACEAGORA_TELEMETRY_DT_MAX_ORBIT" => "10.0",
            "SPACEAGORA_TELEMETRY_RELTOL_ORBIT" => "1e-12",
            "SPACEAGORA_TELEMETRY_ABSTOL_ORBIT" => "1e-12",
            "SPACEAGORA_TELEMETRY_RELTOL_ATM" => "1e-12",
            "SPACEAGORA_TELEMETRY_ABSTOL_ATM" => "1e-12"
        ) do
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
            "drag_enabled" => true,
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
            profile=_verification_profile_from_env(),
            out_summary=joinpath(tmp, "summary.csv"),
            out_errors=joinpath(tmp, "errors.csv"),
            manifest_path=manifest_path,
            enforce=false,
            generate_plots=false
        )

        return withenv(
            "SPACEAGORA_TELEMETRY_SOLVER_MODE" => "tsit5",
            "SPACEAGORA_TELEMETRY_DT_MAX_ORBIT" => "10.0",
            "SPACEAGORA_TELEMETRY_RELTOL_ORBIT" => "1e-12",
            "SPACEAGORA_TELEMETRY_ABSTOL_ORBIT" => "1e-12",
            "SPACEAGORA_TELEMETRY_RELTOL_ATM" => "1e-12",
            "SPACEAGORA_TELEMETRY_ABSTOL_ATM" => "1e-12"
        ) do
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

function _scenario_position_rmse(errors::DataFrame, scenario_name::String; nfirst::Int=10)
    rows = errors[(errors.scenario .== scenario_name) .& in.(errors.event, Ref(["state_x_time", "state_y_time", "state_z_time"])), :]
    xrow = sort(rows[rows.event .== "state_x_time", :], :idx)
    yrow = sort(rows[rows.event .== "state_y_time", :], :idx)
    zrow = sort(rows[rows.event .== "state_z_time", :], :idx)

    n = min(nrow(xrow), nrow(yrow), nrow(zrow))
    @test n >= nfirst

    ex = Float64.(xrow.error_km[1:n])
    ey = Float64.(yrow.error_km[1:n])
    ez = Float64.(zrow.error_km[1:n])

    pos_sq = ex .^ 2 .+ ey .^ 2 .+ ez .^ 2
    first10_rmse_km = sqrt(mean(pos_sq[1:nfirst]))
    full_rmse_km = sqrt(mean(pos_sq))

    return (first10_rmse_km=first10_rmse_km, full_rmse_km=full_rmse_km, n_points=n)
end

function _plot_cygnss_error_timeseries(errors::DataFrame, outpath::String)::String
    rows = errors[errors.scenario .== "cygnss_48hr_pvt", :]
    xrow = sort(rows[rows.event .== "state_x_time", :], :idx)
    yrow = sort(rows[rows.event .== "state_y_time", :], :idx)
    zrow = sort(rows[rows.event .== "state_z_time", :], :idx)

    n = min(nrow(xrow), nrow(yrow), nrow(zrow))
    n >= 2 || error("Need at least two aligned CYGNSS error samples for plotting.")

    t_hr = Float64.(xrow.telemetry_axis[1:n]) ./ 3600.0
    ex = Float64.(xrow.error_km[1:n])
    ey = Float64.(yrow.error_km[1:n])
    ez = Float64.(zrow.error_km[1:n])

    total_pos_err = sqrt.(ex .^ 2 .+ ey .^ 2 .+ ez .^ 2)
    running_total_rmse = sqrt.(cumsum(total_pos_err .^ 2) ./ collect(1:n))
    total_rmse = sqrt(mean(total_pos_err .^ 2))

    p = plot(
        t_hr,
        ex,
        label="x error",
        xlabel="Time (hr)",
        ylabel="Error (km)",
        title="CYGNSS 48hr Position Error Time Series",
        lw=1.2,
        alpha=0.9,
        legend=:topright
    )
    plot!(p, t_hr, ey, label="y error", lw=1.2, alpha=0.9)
    plot!(p, t_hr, ez, label="z error", lw=1.2, alpha=0.9)
    plot!(p, t_hr, running_total_rmse, label="total running RMSE", lw=2.0, color=:black)
    hline!(p, [total_rmse], label="total RMSE", linestyle=:dash, color=:black)

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
    rows = errors[errors.scenario .== scenario_name, :]
    xrow = sort(rows[rows.event .== "state_x_time", :], :idx)
    yrow = sort(rows[rows.event .== "state_y_time", :], :idx)
    zrow = sort(rows[rows.event .== "state_z_time", :], :idx)

    n = min(nrow(xrow), nrow(yrow), nrow(zrow))
    n >= 2 || error("Need at least two aligned position-error samples for plotting scenario $(scenario_name).")

    t_s = Float64.(xrow.telemetry_axis[1:n])
    t_axis, x_label = if time_unit == :s
        (t_s, "Time (s)")
    else
        (t_s ./ 3600.0, "Time (hr)")
    end

    ex = Float64.(xrow.error_km[1:n])
    ey = Float64.(yrow.error_km[1:n])
    ez = Float64.(zrow.error_km[1:n])

    total_pos_err = sqrt.(ex .^ 2 .+ ey .^ 2 .+ ez .^ 2)
    running_total_rmse = sqrt.(cumsum(total_pos_err .^ 2) ./ collect(1:n))
    total_rmse = sqrt(mean(total_pos_err .^ 2))

    p = plot(
        t_axis,
        ex,
        label="x error",
        xlabel=x_label,
        ylabel="Error (km)",
        title=title,
        lw=1.2,
        alpha=0.9,
        legend=:topright
    )
    plot!(p, t_axis, ey, label="y error", lw=1.2, alpha=0.9)
    plot!(p, t_axis, ez, label="z error", lw=1.2, alpha=0.9)
    plot!(p, t_axis, running_total_rmse, label="total running RMSE", lw=2.0, color=:black)
    hline!(p, [total_rmse], label="total RMSE", linestyle=:dash, color=:black)

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
    rows = errors[errors.scenario .== "cygnss_48hr_pvt", :]
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

    telem = DataFrame(Arrow.Table(_CYGNSS_48HR_TELEMETRY_FEATHER))
    time_col = _required_column(telem, ["TIME OFFSET", "time"])
    x_col = _required_column(telem, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_X (m)", "pos_ii_1"])
    y_col = _required_column(telem, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_Y (m)", "pos_ii_2"])
    z_col = _required_column(telem, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_Z (m)", "pos_ii_3"])
    vx_col = _required_column(telem, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_X (m/s)", "vel_ii_1"])
    vy_col = _required_column(telem, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_Y (m/s)", "vel_ii_2"])
    vz_col = _required_column(telem, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_Z (m/s)", "vel_ii_3"])

    t_tel = Float64.(telem[!, time_col])
    rx_tel = Float64.(telem[!, x_col]) .* 1.0e-3
    ry_tel = Float64.(telem[!, y_col]) .* 1.0e-3
    rz_tel = Float64.(telem[!, z_col]) .* 1.0e-3
    vx_tel = Float64.(telem[!, vx_col]) .* 1.0e-3
    vy_tel = Float64.(telem[!, vy_col]) .* 1.0e-3
    vz_tel = Float64.(telem[!, vz_col]) .* 1.0e-3

    perm = sortperm(t_tel)
    t_tel = t_tel[perm]
    rx_tel = rx_tel[perm]
    ry_tel = ry_tel[perm]
    rz_tel = rz_tel[perm]
    vx_tel = vx_tel[perm]
    vy_tel = vy_tel[perm]
    vz_tel = vz_tel[perm]

    rx = _interp_series_linear(t_tel, rx_tel, t_s)
    ry = _interp_series_linear(t_tel, ry_tel, t_s)
    rz = _interp_series_linear(t_tel, rz_tel, t_s)
    vx = _interp_series_linear(t_tel, vx_tel, t_s)
    vy = _interp_series_linear(t_tel, vy_tel, t_s)
    vz = _interp_series_linear(t_tel, vz_tel, t_s)

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
    running_total_rmse = sqrt.(cumsum(total_pos_err .^ 2) ./ collect(1:n))
    total_rmse = sqrt(mean(total_pos_err .^ 2))

    p = plot(
        t_hr,
        er,
        label="R error",
        xlabel="Time (hr)",
        ylabel="Error (km)",
        title="CYGNSS 48hr Position Error Time Series (RTN)",
        lw=1.2,
        alpha=0.9,
        legend=:topright
    )
    plot!(p, t_hr, et, label="T error", lw=1.2, alpha=0.9)
    plot!(p, t_hr, en, label="N error", lw=1.2, alpha=0.9)
    plot!(p, t_hr, running_total_rmse, label="total running RMSE", lw=2.0, color=:black)
    hline!(p, [total_rmse], label="total RMSE", linestyle=:dash, color=:black)

    mkpath(dirname(outpath))
    savefig(p, outpath)
    return outpath
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

function _plot_cygnss_orbital_elements(
    errors::DataFrame,
    outpath_elements::String,
    outpath_errors::String
)::Tuple{String, String}
    rows = errors[errors.scenario .== "cygnss_48hr_pvt", :]
    xrow = sort(rows[rows.event .== "state_x_time", :], :idx)
    yrow = sort(rows[rows.event .== "state_y_time", :], :idx)
    zrow = sort(rows[rows.event .== "state_z_time", :], :idx)

    n = min(nrow(xrow), nrow(yrow), nrow(zrow))
    n >= 2 || error("Need at least two aligned CYGNSS samples for orbital element plotting.")

    t_s = Float64.(xrow.telemetry_axis[1:n])
    t_hr = t_s ./ 3600.0

    telem = DataFrame(Arrow.Table(_CYGNSS_48HR_TELEMETRY_FEATHER))
    time_col = _required_column(telem, ["TIME OFFSET", "time"])
    x_col = _required_column(telem, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_X (m)", "pos_ii_1"])
    y_col = _required_column(telem, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_Y (m)", "pos_ii_2"])
    z_col = _required_column(telem, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_Z (m)", "pos_ii_3"])
    vx_col = _required_column(telem, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_X (m/s)", "vel_ii_1"])
    vy_col = _required_column(telem, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_Y (m/s)", "vel_ii_2"])
    vz_col = _required_column(telem, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_Z (m/s)", "vel_ii_3"])

    t_tel = Float64.(telem[!, time_col])
    rx_tel = Float64.(telem[!, x_col]) .* 1.0e-3
    ry_tel = Float64.(telem[!, y_col]) .* 1.0e-3
    rz_tel = Float64.(telem[!, z_col]) .* 1.0e-3
    vx_tel = Float64.(telem[!, vx_col]) .* 1.0e-3
    vy_tel = Float64.(telem[!, vy_col]) .* 1.0e-3
    vz_tel = Float64.(telem[!, vz_col]) .* 1.0e-3

    perm = sortperm(t_tel)
    t_tel = t_tel[perm]
    rx_tel = rx_tel[perm]
    ry_tel = ry_tel[perm]
    rz_tel = rz_tel[perm]
    vx_tel = vx_tel[perm]
    vy_tel = vy_tel[perm]
    vz_tel = vz_tel[perm]

    # Interpolate position/velocity to error times
    rx = _interp_series_linear(t_tel, rx_tel, t_s)
    ry = _interp_series_linear(t_tel, ry_tel, t_s)
    rz = _interp_series_linear(t_tel, rz_tel, t_s)
    vx = _interp_series_linear(t_tel, vx_tel, t_s)
    vy = _interp_series_linear(t_tel, vy_tel, t_s)
    vz = _interp_series_linear(t_tel, vz_tel, t_s)

    planet = TV._planet_from_name("earth")

    # Compute orbital elements from reference state (telemetry r/v)
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
            oe = SM.ControlEffectors.rvtoorbitalelement(r, v, planet)
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

    # Simulated position with errors applied
    ex = Float64.(xrow.error_km[1:n])
    ey = Float64.(yrow.error_km[1:n])
    ez = Float64.(zrow.error_km[1:n])
    rx_sim = rx .+ ex
    ry_sim = ry .+ ey
    rz_sim = rz .+ ez

    # Compute orbital elements from simulated state
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
            oe = SM.ControlEffectors.rvtoorbitalelement(r, v, planet)
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

    # Wrap angular elements for display in canonical 0..360 degree range.
    raan_ref_wrap = mod.(raan_ref_sim, 360.0)
    aop_ref_wrap = mod.(aop_ref_sim, 360.0)
    ta_ref_wrap = mod.(ta_ref_sim, 360.0)
    raan_sim_wrap = mod.(raan_sim, 360.0)
    aop_sim_wrap = mod.(aop_sim, 360.0)
    ta_sim_wrap = mod.(ta_sim, 360.0)

    # Compute orbital element errors as shortest signed angular differences.
    da = sma_sim .- sma_ref_sim
    de = ecc_sim .- ecc_ref_sim
    di = inc_sim .- inc_ref_sim
    dw = _wrap_pm180.(aop_sim_wrap .- aop_ref_wrap)
    dW = _wrap_pm180.(raan_sim_wrap .- raan_ref_wrap)
    df = _wrap_pm180.(ta_sim_wrap .- ta_ref_wrap)

    # Elements plot: 3 x 2
    layout_3x2 = @layout [a b; c d; e f]
    pe1 = plot(t_hr, sma_sim, label="simulated", xlabel="Time (hr)", ylabel="SMA (km)", title="Semi-Major Axis", legend=:best, lw=1.2)
    plot!(pe1, t_hr, sma_ref_sim, label="reference", lw=1.2, alpha=0.75)
    pe2 = plot(t_hr, ecc_sim, label="simulated", xlabel="Time (hr)", ylabel="Ecc", title="Eccentricity", legend=:best, lw=1.2)
    plot!(pe2, t_hr, ecc_ref_sim, label="reference", lw=1.2, alpha=0.75)
    pe3 = plot(t_hr, inc_sim, label="simulated", xlabel="Time (hr)", ylabel="Inc (deg)", title="Inclination", legend=:best, lw=1.2)
    plot!(pe3, t_hr, inc_ref_sim, label="reference", lw=1.2, alpha=0.75)
    pe4 = plot(t_hr, aop_sim_wrap, label="simulated", xlabel="Time (hr)", ylabel="AoP (deg)", title="Argument of Perigee", legend=:best, lw=1.2, ylim=(0.0, 360.0))
    plot!(pe4, t_hr, aop_ref_wrap, label="reference", lw=1.2, alpha=0.75)
    pe5 = plot(t_hr, raan_sim_wrap, label="simulated", xlabel="Time (hr)", ylabel="RAAN (deg)", title="RAAN", legend=:best, lw=1.2, ylim=(0.0, 360.0))
    plot!(pe5, t_hr, raan_ref_wrap, label="reference", lw=1.2, alpha=0.75)
    pe6 = plot(t_hr, ta_sim_wrap, label="simulated", xlabel="Time (hr)", ylabel="TA (deg)", title="True Anomaly", legend=:best, lw=1.2, ylim=(0.0, 360.0))
    plot!(pe6, t_hr, ta_ref_wrap, label="reference", lw=1.2, alpha=0.75)
    p_elements = plot(pe1, pe2, pe3, pe4, pe5, pe6, layout=layout_3x2, size=(1500, 1100))

    # Errors plot: 3 x 2
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

@testset "GMAT Early vs Full Error" begin
    result = _run_gmat_scenario_matrix_result_once()
    scenario_names = unique(String.(result.summary.scenario))

    for scenario_name in scenario_names
        metrics = _scenario_position_rmse(result.errors, scenario_name; nfirst=10)
        @test isfinite(metrics.first10_rmse_km)
        @test isfinite(metrics.full_rmse_km)
        @test metrics.n_points >= 10
        println("$scenario_name position RMSE [km]: first10=$(metrics.first10_rmse_km), full=$(metrics.full_rmse_km)")
    end
end

# @testset "GMAT Scenario Matrix Regression" begin
#     summary = _run_gmat_scenario_matrix_once()
#     @test nrow(summary) == 20

#     for scenario_name in ("earth_j2_tbfalse", "earth_j2_tbfalse_ic_perturbed", "earth_j0_tbfalse", "earth_j2_tbtrue", "earth_j0_tbtrue")
#         rows = summary[summary.scenario .== scenario_name, :]
#         @test nrow(rows) == 4
#         @test all(in.(rows.event, Ref(["altitude_time", "state_x_time", "state_y_time", "state_z_time"])))
#         @test all(rows.n_sim .>= rows.min_eval_points)
#         @test all(isfinite.(Float64.(rows.rmse_km)))
#         @test all(isfinite.(Float64.(rows.max_abs_km)))
#     end

#     rmse_baseline = _scenario_rmse(summary, "earth_j2_tbfalse")
#     rmse_ic = _scenario_rmse(summary, "earth_j2_tbfalse_ic_perturbed")
#     rmse_force = _scenario_rmse(summary, "earth_j0_tbfalse")
#     rmse_j2_tbtrue = _scenario_rmse(summary, "earth_j2_tbtrue")
#     rmse_j0_tbtrue = _scenario_rmse(summary, "earth_j0_tbtrue")

#     @test rmse_baseline < rmse_ic
#     @test abs(rmse_force - rmse_baseline) > 1e-6
#     @test abs(rmse_j2_tbtrue - rmse_baseline) > 1e-7
#     @test abs(rmse_j0_tbtrue - rmse_force) > 1e-7
# end

@testset "GMAT Strict Acceptance All Cases" begin
    summary = _run_gmat_scenario_matrix_once()
    profile = _verification_profile_from_env()
    scenario_names = unique(String.(summary.scenario))

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
        # @test xrow.rmse_km[1] <= 0.1
        # @test yrow.rmse_km[1] <= 0.1
        # @test zrow.rmse_km[1] <= 0.1
        # @test xrow.max_abs_km[1] <= 0.1
        # @test yrow.max_abs_km[1] <= 0.1
        # @test zrow.max_abs_km[1] <= 0.1
    end
end

@testset "CYGNSS Telemetry Folder Data" begin
    cygnss_feather_path = joinpath(_GMAT_REPO_ROOT, "data", "telemetry", "CYGNSS", "cygnss_data_48hr.feather")
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

    @test !any(ismissing, cygnss_feather[!, time_col])
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

@testset "CYGNSS 48hr Orbit Simulation" begin
    result = _run_cygnss_48hr_result_once()
    summary = result.summary
    errors = result.errors
    @test nrow(summary) == 4
    @test all(summary.scenario .== "cygnss_48hr_pvt")
    @test all(in.(summary.event, Ref(["altitude_time", "state_x_time", "state_y_time", "state_z_time"])))

    # Validate that the telemetry span is approximately 48 hours.
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
    _plot_cygnss_error_timeseries(errors, plot_path)
    @test isfile(plot_path)
    plot_path_rtn = joinpath(_GMAT_REPO_ROOT, "output", "cygnss", "cygnss_48hr_error_timeseries_rtn.png")
    _plot_cygnss_error_timeseries_rtn(errors, plot_path_rtn)
    @test isfile(plot_path_rtn)
    plot_path_oe = joinpath(_GMAT_REPO_ROOT, "output", "cygnss", "cygnss_48hr_orbital_elements.png")
    plot_path_oe_err = joinpath(_GMAT_REPO_ROOT, "output", "cygnss", "cygnss_48hr_orbital_element_errors.png")
    _plot_cygnss_orbital_elements(errors, plot_path_oe, plot_path_oe_err)
    @test isfile(plot_path_oe)
    @test isfile(plot_path_oe_err)

    pos_rmse = _scenario_rmse(summary, "cygnss_48hr_pvt")
    println("cygnss_48hr_pvt mean position-axis RMSE [km]: $(pos_rmse)")
    println("cygnss_48hr error plot: $(plot_path)")
    println("cygnss_48hr RTN error plot: $(plot_path_rtn)")
    println("cygnss_48hr orbital elements plot: $(plot_path_oe)")
    println("cygnss_48hr orbital element errors plot: $(plot_path_oe_err)")
    println("cygnss_48hr_pvt_position max error [km]: $(max(Float64.(summary[summary.event .== "state_x_time", :].max_abs_km[1]), Float64.(summary[summary.event .== "state_y_time", :].max_abs_km[1]), Float64.(summary[summary.event .== "state_z_time", :].max_abs_km[1])))")
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

    pos_rmse = _scenario_rmse(summary, "cygnss_48hr_gmat_csv")
    println("cygnss_48hr_gmat_csv mean position-axis RMSE [km]: $(pos_rmse)")
    println("cygnss_48hr_gmat_csv error plot: $(plot_path)")
    @test pos_rmse < 1.0e4
end
