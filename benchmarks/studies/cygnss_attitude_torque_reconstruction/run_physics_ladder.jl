##
# Conservative-force ladder for the 48-hour CYGNSS PVT arc: 50x50 harmonics
# baseline, +Sun/Moon third-body, +SRP, each scored against flight telemetry
# via the real TV.run_verification pipeline. Requires the private telemetry
# (data/telemetry/CYGNSS/cygnss_data_48hr.feather; see
# data/telemetry/PRIVATE_TELEMETRY.md).
#
# Usage: julia --project=. benchmarks/studies/cygnss_attitude_torque_reconstruction/run_physics_ladder.jl
##
include(joinpath(@__DIR__, "harness_48hr_common.jl"))

# ==============================================================================
# Physics ladder: baseline -> +Sun/Moon third-body -> +SRP, scored vs telemetry.
# ==============================================================================
println("=== CYGNSS 48hr physics ladder ===")

configs = [
    ("baseline_harmonics50", Dict{String, Any}()),
    ("plus_sun_moon", Dict{String, Any}("nbody_bodies" => Any["sun", "moon"])),
    ("plus_sunmoon_srp", Dict{String, Any}(
        "nbody_bodies" => Any["sun", "moon"],
        "srp_enabled" => true, "srp_cr" => 1.3, "srp_area_m2" => 1.0)),
]

ladder = NamedTuple[]
mktempdir() do tmp
    traj = _build_cygnss_48hr_reference(tmp, "cygnss_48hr_pvt")
    for (label, extra) in configs
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
                "year" => 2025, "month" => 6, "day" => 6,
                "hour" => 0, "minute" => 0, "second" => 0.0)))
        merge!(scenario, extra)
        manifest_path = joinpath(tmp, "manifest_$(label).toml")
        open(manifest_path, "w") do io
            TOML.print(io, Dict{String, Any}("scenarios" => Any[scenario]))
        end
        req = TV.VerificationRequest(
            profile=TEST_MODE,
            out_summary=joinpath(tmp, "summary_$(label).csv"),
            out_errors=joinpath(tmp, "errors_$(label).csv"),
            manifest_path=manifest_path,
            enforce=false, generate_plots=false)
        r = withenv(pairs(_telemetry_solver_env_overrides())...) do
            TV.run_verification(req)
        end
        pos = _scenario_rmse(r.summary, "cygnss_48hr_pvt")
        vel = _scenario_velocity_rmse(r.summary, "cygnss_48hr_pvt")
        rows = r.summary[in.(r.summary.event, Ref(["state_x_time", "state_y_time", "state_z_time"])), :]
        maxabs = maximum(rows.max_abs_km)
        push!(ladder, (label=label, pos_rmse_km=pos, vel_rmse_ms=vel*1e3, max_abs_km=maxabs))
        println("LADDER >>> $(rpad(label, 24)) pos RMSE $(round(pos, digits=6)) km | vel RMSE $(round(vel*1e3, digits=4)) m/s | max|err| $(round(maxabs, digits=3)) km")
        try
            cp(joinpath(tmp, "errors_$(label).csv"), joinpath(STUDY_DIR, "data", "ladder_errors_$(label).csv"), force=true)
        catch e
            println("(errors table copy skipped: $(e))")
        end
    end
end
println()
println("=== LADDER SUMMARY (48 h, quick profile, rebuilt telemetry) ===")
for r in ladder
    println(rpad(r.label, 26), " pos ", round(r.pos_rmse_km, digits=4), " km   vel ",
        round(r.vel_rmse_ms, digits=4), " m/s   max|err| ", round(r.max_abs_km, digits=3), " km")
end
