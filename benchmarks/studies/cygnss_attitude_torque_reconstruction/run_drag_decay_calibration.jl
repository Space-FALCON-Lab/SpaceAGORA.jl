##
# Drag rung of the CYGNSS 48-hour verification: fixed-cd-scale sweep with the
# full physics stack (50x50 harmonics + Sun/Moon + SRP + free-molecular drag)
# for decay-rate calibration, which — unlike a position-RMSE fit — is
# insensitive to initial-condition error. Post-process the per-run error
# tables with decay_postprocess.py to extract SMA decay slopes, cancel the
# common estimator bias against a drag-free run, and solve for the
# decay-matched cd-scale.
#
# Density source: an along-track density table (time_s, rho columns; e.g.
# NOAA WAM-IPE sampled along the flight track — see
# build_wam_density_table.jl), time-interpolated into the drag model. The
# spacecraft geometry comes from a private-side TOML so mission-derived
# figures stay out of the public repo.
#
# Required private inputs under data/telemetry/CYGNSS/ (see
# data/telemetry/PRIVATE_TELEMETRY.md):
#   cygnss_data_48hr.feather        48 h PVT telemetry (J2000)
#   wam_density_table.csv           along-track density table
#   cygnss_science_geometry.toml    spacecraft dict incl. bus_ram_face
#
# Env overrides:
#   SPACEAGORA_CYGNSS_CD_VALUES     comma list of cd scales (default 0.6,0.8,1.0,1.2)
#
# Usage: julia --project=. benchmarks/studies/cygnss_attitude_torque_reconstruction/run_drag_decay_calibration.jl
##
include(joinpath(@__DIR__, "harness_48hr_common.jl"))
include(joinpath(@__DIR__, "gram_wiring.jl"))

using DelimitedFiles

const _DENSITY_TABLE_PATH = joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS", "wam_density_table.csv")
const _GEOMETRY_PATH = joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS", "cygnss_science_geometry.toml")
for (p, hint) in (
    (_DENSITY_TABLE_PATH, "build_wam_density_table.jl"),
    (_GEOMETRY_PATH, "scripts/dev/fetch_private_telemetry.sh"),
)
    isfile(p) || error("Missing required private input $(p) — see $(hint) and data/telemetry/PRIVATE_TELEMETRY.md")
end

const _DENSITY_TABLE = let
    raw, _ = readdlm(_DENSITY_TABLE_PATH, ',', header=true)
    (t = Float64.(raw[:, 1]), rho = Float64.(raw[:, 2]))
end

@inline function _table_rho(el_time::Float64)::Float64
    t, rho = _DENSITY_TABLE.t, _DENSITY_TABLE.rho
    el_time <= t[1] && return rho[1]
    el_time >= t[end] && return rho[end]
    i = searchsortedfirst(t, el_time)
    w = (el_time - t[i-1]) / (t[i] - t[i-1])
    return rho[i-1] * (1 - w) + rho[i] * w
end

# Density source override: the scenario requests a GRAM model (so the truth
# metadata stays honest about provenance), and this method serves the
# along-track measured table instead of the surrogate grid.
function SimulationModel.EnvironmentModels.getDensity(
    model::SimulationModel.EnvironmentModels.GRAMAtmosphereModelSurrogate,
    h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p::params
)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    EI = p.args.environment_model.EI * 1e3
    drag_state = h - EI <= 0.0
    if h > 2000.0e3
        return 0.0, p.args.environment_model.planet.T_ref, SVector{3, Float64}(0.0, 0.0, 0.0)
    elseif !drag_state && !p.args.mission_configuration.keplerian
        return SimulationModel.EnvironmentModels.density_polyfit(h, p)
    end
    return _table_rho(el_time), 900.0, SVector{3, Float64}(0.0, 0.0, 0.0)
end

const _SPACECRAFT = TOML.parsefile(_GEOMETRY_PATH)["spacecraft"]
const _CD_VALUES = [parse(Float64, strip(s)) for s in
    split(get(ENV, "SPACEAGORA_CYGNSS_CD_VALUES", "0.6,0.8,1.0,1.2"), ",")]

println("=== CYGNSS drag decay-calibration sweep: cd scales $(_CD_VALUES) ===")
results = NamedTuple[]
mktempdir() do tmp
    traj = _build_cygnss_48hr_reference(tmp, "cygnss_48hr_pvt")
    for cd in _CD_VALUES
        label = replace(string(cd), "." => "p")
        scenario = _base_scenario_dict("cygnss_decay_cd$(label)", traj.telemetry_path)
        merge!(scenario, Dict{String, Any}(
            "planet" => "earth",
            "gravity_model" => "inverse_squared",
            "gravity_harmonics_degree" => 50,
            "gravity_harmonics_order" => 50,
            "gravity_harmonics_file" => _GMAT_HARMONICS_EARTH_FILE,
            "nbody_bodies" => Any["sun", "moon"],
            "srp_enabled" => true, "srp_cr" => 1.3, "srp_area_m2" => 1.0,
            "orbit_altitude_mode" => "oblate",
            "drag_enabled" => true,
            "include_wind" => false,
            "EI_km" => 600.0,
            "max_points_quick" => 200000,
            "max_points_full" => 200000,
            "initial_time" => Dict{String, Any}(
                "year" => 2025, "month" => 6, "day" => 6,
                "hour" => 0, "minute" => 0, "second" => 0.0),
            "spacecraft" => deepcopy(_SPACECRAFT),
            "calibration" => Dict{String, Any}(
                "enabled" => true,
                "profiles" => Any["quick"],
                "objective" => "mean_rmse_km",
                "fit_cd_scale" => true,
                "cd_scale_steps" => 1,
                "cd_scale_min" => cd,
                "cd_scale_max" => cd,
                "cr_steps" => 1,
            ),
        ))
        manifest_path = joinpath(tmp, "manifest_decay_$(label).toml")
        open(manifest_path, "w") do io
            TOML.print(io, Dict{String, Any}("scenarios" => Any[scenario]))
        end
        req = TV.VerificationRequest(
            profile=TEST_MODE,
            out_summary=joinpath(tmp, "summary_decay_$(label).csv"),
            out_errors=joinpath(tmp, "errors_decay_$(label).csv"),
            manifest_path=manifest_path,
            enforce=false, generate_plots=false)
        r = withenv(pairs(merge(_telemetry_solver_env_overrides(),
                Dict("SPACEAGORA_GRAM_OFFLINE_SURROGATE" => "auto")))...) do
            TV.run_verification(req)
        end
        rows = r.summary[in.(r.summary.event, Ref(["state_x_time", "state_y_time", "state_z_time"])), :]
        pos = sum(rows.rmse_km) / 3
        println("DECAY_CD >>> cd=$(cd)  pos RMSE $(round(pos, digits=5)) km  retcode $(rows.solver_retcode[1])")
        cp(joinpath(tmp, "errors_decay_$(label).csv"),
           joinpath(STUDY_DIR, "data", "decay_errors_cd$(label).csv"), force=true)
        push!(results, (cd=cd, pos=pos))
    end
end
println("=== SWEEP SUMMARY (position RMSE is IC-contaminated; use decay_postprocess.py for the decay-matched scale) ===")
for r in results
    println("cd_scale ", r.cd, "  pos RMSE ", round(r.pos, digits=4), " km")
end
println("errors tables in $(joinpath(STUDY_DIR, "data")) — next: python3 decay_postprocess.py")
