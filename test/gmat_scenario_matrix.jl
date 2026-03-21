using Test
using TOML
using Arrow
using CSV
using DataFrames
using Statistics

const _GMAT_REPO_ROOT = isdefined(Main, :REPO_ROOT) ? Main.REPO_ROOT : normpath(joinpath(@__DIR__, ".."))

if !isdefined(@__MODULE__, :TelemetryVerification)
    include(joinpath(_GMAT_REPO_ROOT, "src", "analysis", "verification", "TelemetryVerification.jl"))
end

const TV = TelemetryVerification

const _GMAT_EXAMPLES_DIR = joinpath(
    _GMAT_REPO_ROOT,
    "data",
    "telemetry",
    "GMAT_Examples"
)

const _GMAT_HARMONICS_FILE = "data/Gravity_harmonics_data/egm96.csv"

# ============================================================================
# TEST MODE CONFIGURATION
# ============================================================================
# Change this to :full for comprehensive testing with more evaluation points
# Change to :quick for faster testing with fewer evaluation points
const TEST_MODE::Symbol = :full
# ============================================================================

@inline function _verification_profile_from_env()::Symbol
    return TEST_MODE
end

@inline function _strict_position_rmse_limit_km(scenario_name::String, profile::Symbol)::Float64
    if scenario_name == "earth_harm50_50_tbfalse"
        return profile == :full ? 0.45 : 0.05
    elseif startswith(scenario_name, "earth_j2_")
        return profile == :full ? 0.20 : 0.03
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
    df = CSV.read(csv_path, DataFrame)

    time_col = _required_column(df, ["Sat.ElapsedSecs", "DefaultSC.ElapsedSecs"])
    x_col = _required_column(df, ["Sat.X", "DefaultSC.EarthMJ2000Eq.X"])
    y_col = _required_column(df, ["Sat.Y", "DefaultSC.EarthMJ2000Eq.Y"])
    z_col = _required_column(df, ["Sat.Z", "DefaultSC.EarthMJ2000Eq.Z"])

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
            "bus_dims_m" => Any[2.05, 3.7, 2.8],
            "panel_dims_m" => Any[0.01, 2.85, 1.0],
            "bus_mass_kg" => 620.0,
            "panel_mass_each_kg" => 10.0,
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

function _run_gmat_scenario_matrix_result_once()::TV.VerificationResult
    if _GMAT_MATRIX_RESULT_CACHE[] !== nothing
        return _GMAT_MATRIX_RESULT_CACHE[]
    end

    earth_j2_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Earth_J2_TBFalse.csv")
    earth_j0_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Earth_J0_TBFalse.csv")
    earth_harm50_50_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Earth_Harm50_50_TBfalse.csv")
    earth_j2_tbtrue_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Earth_J2_TBTrue.csv")
    earth_j0_tbtrue_csv = joinpath(_GMAT_EXAMPLES_DIR, "Sim_Earth_J0_TBTrue.csv")
    @test isfile(earth_j2_csv)
    @test isfile(earth_j0_csv)
    @test isfile(earth_harm50_50_csv)
    @test isfile(earth_j2_tbtrue_csv)
    @test isfile(earth_j0_tbtrue_csv)
    @test isfile(joinpath(_GMAT_REPO_ROOT, _GMAT_HARMONICS_FILE))

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
        traj_harm50_50 = _build_time_aligned_reference(earth_harm50_50_csv, "earth", tmp, "earth_harm50_50_tbfalse")
        traj_j2_tbtrue = _build_time_aligned_reference(earth_j2_tbtrue_csv, "earth", tmp, "earth_j2_tbtrue")
        traj_j0_tbtrue = _build_time_aligned_reference(earth_j0_tbtrue_csv, "earth", tmp, "earth_j0_tbtrue")

        println(
            "SpaceAGORA initial orbital elements [a_km,e,i_deg,aop_deg,raan_deg,ta_deg]: " *
            "earth_j2_tbfalse=($(traj_j2.sma_km),$(traj_j2.ecc),$(traj_j2.inc_deg),$(traj_j2.aop_deg),$(traj_j2.raan_deg),$(traj_j2.ta_deg)); " *
            # "earth_j2_tbfalse_ic_perturbed=($(traj_j2_ic.sma_km),$(traj_j2_ic.ecc),$(traj_j2_ic.inc_deg),$(traj_j2_ic.aop_deg),$(traj_j2_ic.raan_deg),$(traj_j2_ic.ta_deg)); " *
            "earth_j0_tbfalse=($(traj_j0.sma_km),$(traj_j0.ecc),$(traj_j0.inc_deg),$(traj_j0.aop_deg),$(traj_j0.raan_deg),$(traj_j0.ta_deg)); " *
            "earth_harm50_50_tbfalse=($(traj_harm50_50.sma_km),$(traj_harm50_50.ecc),$(traj_harm50_50.inc_deg),$(traj_harm50_50.aop_deg),$(traj_harm50_50.raan_deg),$(traj_harm50_50.ta_deg)); " *
            "earth_j2_tbtrue=($(traj_j2_tbtrue.sma_km),$(traj_j2_tbtrue.ecc),$(traj_j2_tbtrue.inc_deg),$(traj_j2_tbtrue.aop_deg),$(traj_j2_tbtrue.raan_deg),$(traj_j2_tbtrue.ta_deg)); " *
            "earth_j0_tbtrue=($(traj_j0_tbtrue.sma_km),$(traj_j0_tbtrue.ecc),$(traj_j0_tbtrue.inc_deg),$(traj_j0_tbtrue.aop_deg),$(traj_j0_tbtrue.raan_deg),$(traj_j0_tbtrue.ta_deg))"
        )

        baseline = _base_scenario_dict("earth_j0_tbfalse", traj_j0.telemetry_path)
        merge!(baseline, Dict{String, Any}(
            "planet" => "earth",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 0,
            "gravity_harmonics_order" => 0,
            "gravity_harmonics_file" => _GMAT_HARMONICS_FILE,
            "nbody_bodies" => Any[],
            "orbit_altitude_mode" => "oblate"
        ))

        j2_tbfalse = _base_scenario_dict("earth_j2_tbfalse", traj_j2.telemetry_path)
        merge!(j2_tbfalse, Dict{String, Any}(
            "planet" => "earth",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 2,
            "gravity_harmonics_order" => 0,
            "gravity_harmonics_file" => _GMAT_HARMONICS_FILE,
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
            "gravity_harmonics_file" => _GMAT_HARMONICS_FILE,
            "nbody_bodies" => Any["sun", "moon"],
            "orbit_altitude_mode" => "oblate"
        ))

        harm50_50_tbfalse = _base_scenario_dict("earth_harm50_50_tbfalse", traj_harm50_50.telemetry_path)
        merge!(harm50_50_tbfalse, Dict{String, Any}(
            "planet" => "earth",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 50,
            "gravity_harmonics_order" => 50,
            "gravity_harmonics_file" => _GMAT_HARMONICS_FILE,
            "nbody_bodies" => Any[],
            "orbit_altitude_mode" => "oblate"
        ))

        j0_tbtrue = _base_scenario_dict("earth_j0_tbtrue", traj_j0_tbtrue.telemetry_path)
        merge!(j0_tbtrue, Dict{String, Any}(
            "planet" => "earth",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 0,
            "gravity_harmonics_order" => 0,
            "gravity_harmonics_file" => _GMAT_HARMONICS_FILE,
            "nbody_bodies" => Any["sun", "moon"],
            "orbit_altitude_mode" => "oblate",
            "srp_enabled" => false, # intentionally disable SRP for this case to isolate the effect of third bodies vs SRP in the j2_tbtrue case
            "srp_area_m2" => 5.0,
            "srp_cr" => 1.35
        ))

        manifest_path = joinpath(tmp, "gmat_scenario_matrix_manifest.toml")
        open(manifest_path, "w") do io
            TOML.print(io, Dict{String, Any}("scenarios" => Any[baseline, j2_tbfalse, j2_tbtrue, j0_tbtrue, harm50_50_tbfalse]))
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
