Base.@kwdef struct CYGNSSEstimatorValidationRequest
    truth_path::String = joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS", "cyg04_nasa_pvt_96hr.feather")
    measurement_path::String = joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS", "cygnss_data_48hr.feather")
    ic_path::String = joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS", "cygnss_data_48hr.feather")
    gravity_harmonics_file::String = joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv")
    out_summary::String = joinpath(DEFAULT_OUTPUT_DIR, "cygnss_estimator_validation_summary.csv")
    out_errors::String = joinpath(DEFAULT_OUTPUT_DIR, "cygnss_estimator_validation_errors.csv")
    profile::Symbol = :quick
    max_points::Int = 0
    position_noise_m::Float64 = 100.0
    velocity_noise_mps::Float64 = 10.0
    process_accel_noise_mps2::Float64 = 1.0e-5
    initial_position_sigma_m::Float64 = 1_000.0
    initial_velocity_sigma_mps::Float64 = 10.0
    innovation_gate_sigma::Float64 = 0.0
    position_innovation_gate_km::Float64 = 25.0
    velocity_innovation_gate_mps::Float64 = 100.0
    include_velocity_measurements::Bool = false
    run_posvel_case::Bool = true
    run_sensitivity::Bool = false
    truth_time_offset_s::Float64 = 0.0
    generate_plots::Bool = true
    plots_dir::String = joinpath(DEFAULT_OUTPUT_DIR, "cygnss_estimator_validation_plots")
end

Base.@kwdef struct CYGNSSEstimatorValidationResult
    summary::DataFrame
    errors::DataFrame
    summary_path::String
    errors_path::String
    plots_dir::String
end

Base.@kwdef struct CYGNSSTimingOffsetSweepRequest
    truth_path::String = joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS", "cyg04_nasa_pvt_96hr.feather")
    measurement_path::String = joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS", "cygnss_data_48hr.feather")
    out_summary::String = joinpath(DEFAULT_OUTPUT_DIR, "cygnss_timing_offset_sweep.csv")
    plots_dir::String = joinpath(DEFAULT_OUTPUT_DIR, "cygnss_estimator_validation_plots")
    min_offset_s::Float64 = -120.0
    max_offset_s::Float64 = 120.0
    step_s::Float64 = 1.0
    max_points::Int = 0
    generate_plots::Bool = true
end

Base.@kwdef struct CYGNSSTimingOffsetSweepResult
    summary::DataFrame
    summary_path::String
    plots_dir::String
    best_raw_gps_rmse_offset_s::Float64
    best_raw_gps_initial_offset_s::Float64
end

@inline function _cygnss_required_column(df::DataFrame, name::String, context::String)
    name in names(df) || throw(ArgumentError("Missing required CYGNSS $context column '$name'."))
    return _to_float_vector(df[!, name], "cygnss-$context:$name")
end

@inline function _cygnss_first_existing_column(df::DataFrame, candidates::Vector{String}, context::String)
    for name in candidates
        name in names(df) && return _to_float_vector(df[!, name], "cygnss-$context:$name")
    end
    throw(ArgumentError("Missing required CYGNSS $context column. Tried: $(join(candidates, ", "))"))
end

function _cygnss_keep_indices(n::Int, max_points::Int)
    n >= 2 || throw(ArgumentError("CYGNSS validation requires at least 2 samples, got $n."))
    if max_points <= 0 || max_points >= n
        return collect(1:n)
    end
    max_points >= 2 || throw(ArgumentError("max_points must be 0 or >= 2, got $max_points."))
    # Cover the full measurement-supported span even for smoke runs.
    return unique(round.(Int, range(1, n, length=max_points)))
end

function _load_cygnss_truth_series(req::CYGNSSEstimatorValidationRequest)
    isfile(req.truth_path) || throw(ArgumentError("Missing CYGNSS truth file: $(req.truth_path)"))
    df = DataFrame(Arrow.Table(req.truth_path))
    nrow(df) > 0 || throw(ArgumentError("CYGNSS truth file is empty: $(req.truth_path)"))

    t_sort = _cygnss_first_existing_column(df, ["time", "pvt_unix_seconds"], "truth-sort-time")
    perm = sortperm(t_sort)
    t_native = _cygnss_first_existing_column(df, ["time", "pvt_unix_seconds"], "truth-time")[perm]
    t_rel = t_native .- t_native[1]

    truth_r_km = hcat(
        _cygnss_required_column(df, "pos_ii_1", "truth-x")[perm],
        _cygnss_required_column(df, "pos_ii_2", "truth-y")[perm],
        _cygnss_required_column(df, "pos_ii_3", "truth-z")[perm]
    ) .* 1.0e-3
    truth_v_kmps = hcat(
        _cygnss_required_column(df, "vel_ii_1", "truth-vx")[perm],
        _cygnss_required_column(df, "vel_ii_2", "truth-vy")[perm],
        _cygnss_required_column(df, "vel_ii_3", "truth-vz")[perm]
    ) .* 1.0e-3

    return (time_s=t_rel, r_km=truth_r_km, v_kmps=truth_v_kmps)
end

function _load_cygnss_measurement_series(req::CYGNSSEstimatorValidationRequest)
    isfile(req.measurement_path) || throw(ArgumentError("Missing CYGNSS GPS measurement file: $(req.measurement_path)"))
    # Loads SPICE leapseconds/PCK/EOP kernels needed for ECEF <-> J2000 transforms.
    _planet_from_name("earth")
    df = DataFrame(Arrow.Table(req.measurement_path))
    nrow(df) > 0 || throw(ArgumentError("CYGNSS GPS measurement file is empty: $(req.measurement_path)"))

    t_sort = _cygnss_first_existing_column(df, ["time", "pvt_unix_seconds"], "measurement-sort-time")
    perm = sortperm(t_sort)
    keep = _cygnss_keep_indices(length(perm), req.max_points)
    idx = perm[keep]

    t_native = _cygnss_first_existing_column(df, ["time", "pvt_unix_seconds"], "measurement-time")[idx]
    t_rel = t_native .- t_native[1]

    if all(name -> name in names(df), ("pos_ii_1", "pos_ii_2", "pos_ii_3", "vel_ii_1", "vel_ii_2", "vel_ii_3"))
        meas_r_km = hcat(
            _cygnss_required_column(df, "pos_ii_1", "measurement-x")[idx],
            _cygnss_required_column(df, "pos_ii_2", "measurement-y")[idx],
            _cygnss_required_column(df, "pos_ii_3", "measurement-z")[idx]
        ) .* 1.0e-3
        meas_v_kmps = hcat(
            _cygnss_required_column(df, "vel_ii_1", "measurement-vx")[idx],
            _cygnss_required_column(df, "vel_ii_2", "measurement-vy")[idx],
            _cygnss_required_column(df, "vel_ii_3", "measurement-vz")[idx]
        ) .* 1.0e-3
    else
        pvt_t_unix = _cygnss_required_column(df, "pvt_unix_seconds", "pvt-unix")[idx]
        pvt_t_rel = pvt_t_unix .- pvt_t_unix[1]
        pvt_et = if "pvt_et_seconds" in names(df)
            _cygnss_required_column(df, "pvt_et_seconds", "pvt-et")[idx]
        else
            _initial_time_et(InitialTime(2025, 6, 6, 0, 0, 0.0)) .+ pvt_t_rel
        end

        meas_r_ecef_m = hcat(
            _cygnss_required_column(df, "sc_pos_x_pvt_m", "measurement-x")[idx],
            _cygnss_required_column(df, "sc_pos_y_pvt_m", "measurement-y")[idx],
            _cygnss_required_column(df, "sc_pos_z_pvt_m", "measurement-z")[idx]
        )
        meas_v_ecef_mps = hcat(
            _cygnss_required_column(df, "sc_vel_x_pvt_mps", "measurement-vx")[idx],
            _cygnss_required_column(df, "sc_vel_y_pvt_mps", "measurement-vy")[idx],
            _cygnss_required_column(df, "sc_vel_z_pvt_mps", "measurement-vz")[idx]
        )

        meas_r_km = Matrix{Float64}(undef, length(idx), 3)
        meas_v_kmps = Matrix{Float64}(undef, length(idx), 3)
        @inbounds for i in eachindex(idx)
            r_j2000_m, v_j2000_mps = _planet_fixed_to_j2000_state(
                "earth",
                pvt_et[i],
                SVector{3, Float64}(meas_r_ecef_m[i, 1], meas_r_ecef_m[i, 2], meas_r_ecef_m[i, 3]),
                SVector{3, Float64}(meas_v_ecef_mps[i, 1], meas_v_ecef_mps[i, 2], meas_v_ecef_mps[i, 3])
            )
            meas_r_km[i, :] .= r_j2000_m .* 1.0e-3
            meas_v_kmps[i, :] .= v_j2000_mps .* 1.0e-3
        end
    end

    return (time_s=t_rel, r_km=meas_r_km, v_kmps=meas_v_kmps, source_indices=idx)
end

function _interp_matrix_columns(source_time::Vector{Float64}, values::Matrix{Float64}, query_time::Vector{Float64})
    out = Matrix{Float64}(undef, length(query_time), size(values, 2))
    @inbounds for j in 1:size(values, 2)
        out[:, j] .= _interp_linear(source_time, collect(values[:, j]), query_time)
    end
    return out
end

function _load_cygnss_truth_and_measurements(req::CYGNSSEstimatorValidationRequest)
    truth = _load_cygnss_truth_series(req)
    measurements = _load_cygnss_measurement_series(req)
    truth_query_time = measurements.time_s .+ req.truth_time_offset_s
    truth_r_km = _interp_matrix_columns(truth.time_s, truth.r_km, truth_query_time)
    truth_v_kmps = _interp_matrix_columns(truth.time_s, truth.v_kmps, truth_query_time)

    return (
        time_s=measurements.time_s,
        truth_r_km=truth_r_km,
        truth_v_kmps=truth_v_kmps,
        measurement_time_s=measurements.time_s,
        measurement_r_km=measurements.r_km,
        measurement_v_kmps=measurements.v_kmps,
        source_indices=measurements.source_indices
    )
end

function _cygnss_48hr_initial_condition(req::CYGNSSEstimatorValidationRequest)
    isfile(req.ic_path) || throw(ArgumentError("Missing CYGNSS 48-hour IC file: $(req.ic_path)"))
    df = DataFrame(Arrow.Table(req.ic_path))
    sort_col = _cygnss_first_existing_column(df, ["TIME OFFSET", "time"], "ic-time")
    perm = sortperm(sort_col)
    i0 = perm[1]
    x_col = _cygnss_first_existing_column(df, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_X (m)", "pos_ii_1"], "ic-x")
    y_col = _cygnss_first_existing_column(df, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_Y (m)", "pos_ii_2"], "ic-y")
    z_col = _cygnss_first_existing_column(df, ["OBS4.ENG_PVT.DDMI_PVT_SCPOS_Z (m)", "pos_ii_3"], "ic-z")
    vx_col = _cygnss_first_existing_column(df, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_X (m/s)", "vel_ii_1"], "ic-vx")
    vy_col = _cygnss_first_existing_column(df, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_Y (m/s)", "vel_ii_2"], "ic-vy")
    vz_col = _cygnss_first_existing_column(df, ["OBS4.ENG_PVT.DDMI_PVT_SCVEL_Z (m/s)", "vel_ii_3"], "ic-vz")
    return (
        r_km=SVector{3, Float64}(x_col[i0], y_col[i0], z_col[i0]) .* 1.0e-3,
        v_kmps=SVector{3, Float64}(vx_col[i0], vy_col[i0], vz_col[i0]) .* 1.0e-3
    )
end

function _write_cygnss_truth_arrow(req::CYGNSSEstimatorValidationRequest, series, path::String)
    planet = _planet_from_name("earth")
    ic = _cygnss_48hr_initial_condition(req)
    oe0 = rvtoorbitalelement(ic.r_km .* 1.0e3, ic.v_kmps .* 1.0e3, planet)
    r_norm = sqrt.(series.truth_r_km[:, 1] .^ 2 .+ series.truth_r_km[:, 2] .^ 2 .+ series.truth_r_km[:, 3] .^ 2)
    alt_km = r_norm .- planet.Rp_e * 1.0e-3
    n = length(series.time_s)
    df = DataFrame(
        time_s=series.time_s,
        altitude_km=alt_km,
        x_km=series.truth_r_km[:, 1],
        y_km=series.truth_r_km[:, 2],
        z_km=series.truth_r_km[:, 3],
        sma_km=fill(oe0[1] * 1.0e-3, n),
        ecc=fill(oe0[2], n),
        inc_deg=fill(rad2deg(oe0[3]), n),
        raan_deg=fill(rad2deg(oe0[4]), n),
        aop_deg=fill(rad2deg(oe0[5]), n),
        ta_deg=fill(rad2deg(oe0[6]), n),
        x_ic_km=fill(ic.r_km[1], n),
        y_ic_km=fill(ic.r_km[2], n),
        z_ic_km=fill(ic.r_km[3], n),
        vx_ic_kmps=fill(ic.v_kmps[1], n),
        vy_ic_kmps=fill(ic.v_kmps[2], n),
        vz_ic_kmps=fill(ic.v_kmps[3], n)
    )
    Arrow.write(path, df)
    return path
end

function _cygnss_estimator_scenario(req::CYGNSSEstimatorValidationRequest, truth_arrow::String)
    return Dict{String, Any}(
        "name" => "cyg04_96hr_inertial_truth",
        "kind" => "time_aligned_state",
        "events" => Any["altitude_time", "state_x_time", "state_y_time", "state_z_time"],
        "telemetry" => truth_arrow,
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
        "max_points_quick" => 500000,
        "max_points_full" => 500000,
        "min_eval_points" => 2,
        "units" => Dict{String, Any}(
            "x" => "s",
            "altitude_time" => "km",
            "state_x_time" => "km",
            "state_y_time" => "km",
            "state_z_time" => "km"
        ),
        "tolerances_quick" => _cygnss_loose_tolerances(),
        "tolerances_full" => _cygnss_loose_tolerances(),
        "planet" => "earth",
        "gravity_model" => "inverse_squared",
        "gravity_harmonics_degree" => 50,
        "gravity_harmonics_order" => 50,
        "gravity_harmonics_file" => req.gravity_harmonics_file,
        "nbody_bodies" => Any[],
        "orbit_altitude_mode" => "oblate",
        "drag_enabled" => false,
        "include_wind" => false,
        "cartesian_ic_frame" => "inertial",
        "comparison_frame" => "inertial",
        "comparison_mode" => "time_aligned_state",
        "initial_time" => Dict{String, Any}(
            "year" => 2025,
            "month" => 6,
            "day" => 6,
            "hour" => 0,
            "minute" => 0,
            "second" => 0.0
        ),
        "EI_km" => 300.0,
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
            "assumption_id" => "cygnss_estimator_validation_v1",
            "atmosphere_model" => "none",
            "atmosphere_dataset" => "none",
            "space_weather_model" => "none",
            "solar_flux_model" => "none",
            "gram_seed" => 1001,
            "gram_perturbation_scales" => Any[0.0, 0.0, 0.0, 0.0],
            "gram_offline_surrogate" => "off",
            "gram_static_grid" => false,
            "gram_track_cache" => false,
            "gram_global_lock" => "on"
        ),
        "calibration" => Dict{String, Any}("enabled" => false)
    )
end

function _cygnss_loose_tolerances()
    row = Dict("max_abs_km" => 1.0e9, "max_nmae" => 1.0e9, "max_rmse_km" => 1.0e9)
    return Dict{String, Any}(
        "altitude_time" => row,
        "state_x_time" => row,
        "state_y_time" => row,
        "state_z_time" => row
    )
end

function _run_cygnss_open_loop(req::CYGNSSEstimatorValidationRequest, series)
    return mktempdir() do tmp
        truth_arrow = _write_cygnss_truth_arrow(req, series, joinpath(tmp, "cyg04_96hr_truth.arrow"))
        manifest_path = joinpath(tmp, "cygnss_estimator_manifest.toml")
        open(manifest_path, "w") do io
            TOML.print(io, Dict{String, Any}("scenarios" => Any[_cygnss_estimator_scenario(req, truth_arrow)]))
        end
        sim_summary = joinpath(tmp, "open_loop_summary.csv")
        sim_errors = joinpath(tmp, "open_loop_errors.csv")
        result = withenv(
            "SPACEAGORA_TELEMETRY_SOLVER_MODE" => get(ENV, "SPACEAGORA_TELEMETRY_SOLVER_MODE", "dp8"),
            "SPACEAGORA_TELEMETRY_DT_MAX_ORBIT" => get(ENV, "SPACEAGORA_TELEMETRY_DT_MAX_ORBIT", "10.0"),
            "SPACEAGORA_TELEMETRY_RELTOL_ORBIT" => get(ENV, "SPACEAGORA_TELEMETRY_RELTOL_ORBIT", "1e-12"),
            "SPACEAGORA_TELEMETRY_ABSTOL_ORBIT" => get(ENV, "SPACEAGORA_TELEMETRY_ABSTOL_ORBIT", "1e-12")
        ) do
            run_verification(VerificationRequest(
                profile=req.profile,
                out_summary=sim_summary,
                out_errors=sim_errors,
                manifest_path=manifest_path,
                enforce=false,
                generate_plots=false
            ))
        end
        return result
    end
end

function _extract_open_loop_state(errors::DataFrame, n::Int)
    x = _extract_event_interp(errors, "state_x_time", n)
    y = _extract_event_interp(errors, "state_y_time", n)
    z = _extract_event_interp(errors, "state_z_time", n)
    return hcat(x, y, z)
end

function _extract_event_interp(errors::DataFrame, event::String, n::Int)
    rows = errors[errors.event .== event, :]
    nrow(rows) == n || throw(ArgumentError("Expected $n rows for $event, got $(nrow(rows))."))
    sort!(rows, :idx)
    return Float64.(rows.sim_interp_value_km)
end

function _state_from_position_series(r_km::Matrix{Float64}, time_s::Vector{Float64})
    return hcat(
        r_km,
        _differentiate_series(collect(r_km[:, 1]), time_s),
        _differentiate_series(collect(r_km[:, 2]), time_s),
        _differentiate_series(collect(r_km[:, 3]), time_s)
    )
end

function _cygnss_process_q(dt::Float64, sigma_a_mps2::Float64)
    q = zeros(Float64, 6, 6)
    sigma = sigma_a_mps2 * 1.0e-3
    dt = max(dt, eps(Float64))
    qpp = sigma^2 * dt^3 / 3.0
    qpv = sigma^2 * dt^2 / 2.0
    qvv = sigma^2 * dt
    @inbounds for axis in 1:3
        vaxis = axis + 3
        q[axis, axis] = qpp
        q[axis, vaxis] = qpv
        q[vaxis, axis] = qpv
        q[vaxis, vaxis] = qvv
    end
    return q
end

function _cygnss_error_state_transition(dt::Float64)
    F = Matrix{Float64}(I, 6, 6)
    dt = max(dt, eps(Float64))
    @inbounds for axis in 1:3
        F[axis, axis + 3] = dt
    end
    return F
end

function _run_error_state_filter(
    time_s::Vector{Float64},
    sim_state::Matrix{Float64},
    measurement_state::Matrix{Float64};
    position_noise_m::Float64,
    velocity_noise_mps::Float64,
    process_accel_noise_mps2::Float64,
    initial_position_sigma_m::Float64,
    initial_velocity_sigma_mps::Float64,
    innovation_gate_sigma::Float64,
    position_innovation_gate_km::Float64,
    velocity_innovation_gate_mps::Float64,
    include_velocity_measurements::Bool
)
    n = length(time_s)
    n == size(sim_state, 1) == size(measurement_state, 1) || throw(ArgumentError("Filter input row count mismatch."))
    est = Matrix{Float64}(undef, n, 6)
    innovations = Matrix{Float64}(undef, n, include_velocity_measurements ? 6 : 3)
    nis = Vector{Float64}(undef, n)
    accepted = Vector{Bool}(undef, n)

    δx = zeros(Float64, 6)
    p_pos = (initial_position_sigma_m * 1.0e-3)^2
    p_vel = (initial_velocity_sigma_mps * 1.0e-3)^2
    P = Diagonal(vcat(fill(p_pos, 3), fill(p_vel, 3))) |> Matrix

    H = include_velocity_measurements ? Matrix{Float64}(I, 6, 6) : hcat(Matrix{Float64}(I, 3, 3), zeros(Float64, 3, 3))
    r_pos = (position_noise_m * 1.0e-3)^2
    r_vel = (velocity_noise_mps * 1.0e-3)^2
    R = include_velocity_measurements ? Diagonal(vcat(fill(r_pos, 3), fill(r_vel, 3))) |> Matrix : Diagonal(fill(r_pos, 3)) |> Matrix
    I6 = Matrix{Float64}(I, 6, 6)

    @inbounds for k in 1:n
        if k > 1
            dt = time_s[k] - time_s[k - 1]
            F = _cygnss_error_state_transition(dt)
            δx = F * δx
            P = F * P * F' + _cygnss_process_q(dt, process_accel_noise_mps2)
        end
        x_nom = collect(sim_state[k, :])
        x_pred = x_nom .+ δx
        z = include_velocity_measurements ? collect(measurement_state[k, :]) : collect(measurement_state[k, 1:3])
        y = z .- H * x_pred
        S = H * P * H' + R
        nis_k = dot(y, S \ y)
        position_innovation_norm_km = norm(y[1:3])
        gate_sigma = innovation_gate_sigma <= 0.0 ? Inf : innovation_gate_sigma
        nis_limit = isfinite(gate_sigma) ? gate_sigma^2 * length(y) : Inf
        pos_limit = position_innovation_gate_km <= 0.0 ? Inf : position_innovation_gate_km
        vel_limit = velocity_innovation_gate_mps <= 0.0 ? Inf : velocity_innovation_gate_mps * 1.0e-3
        velocity_innovation_norm_kmps = include_velocity_measurements ? norm(y[4:6]) : 0.0
        use_measurement = isfinite(nis_k) &&
            nis_k <= nis_limit &&
            position_innovation_norm_km <= pos_limit &&
            velocity_innovation_norm_kmps <= vel_limit
        if use_measurement
            K = (P * H') / S
            δx = δx + K * y
            P = (I6 - K * H) * P * (I6 - K * H)' + K * R * K'
        end
        est[k, :] .= x_nom .+ δx
        innovations[k, :] .= y
        nis[k] = nis_k
        accepted[k] = use_measurement
    end
    return (state=est, innovations=innovations, nis=nis, accepted=accepted)
end

@inline function _safe_quantile(v::Vector{Float64}, p::Float64)
    isempty(v) && return NaN
    return quantile(v, p)
end

function _rtn_errors_km(r_truth::Matrix{Float64}, v_truth::Matrix{Float64}, err_r::Matrix{Float64})
    n = size(r_truth, 1)
    out = Matrix{Float64}(undef, n, 3)
    @inbounds for i in 1:n
        r = SVector{3, Float64}(r_truth[i, 1], r_truth[i, 2], r_truth[i, 3])
        v = SVector{3, Float64}(v_truth[i, 1], v_truth[i, 2], v_truth[i, 3])
        e = SVector{3, Float64}(err_r[i, 1], err_r[i, 2], err_r[i, 3])
        rhat = r / norm(r)
        nhat = cross(r, v)
        nhat = nhat / norm(nhat)
        that = cross(nhat, rhat)
        out[i, 1] = dot(e, rhat)
        out[i, 2] = dot(e, that)
        out[i, 3] = dot(e, nhat)
    end
    return out
end

function _drift_rate_km_per_day(time_s::Vector{Float64}, pos_norm_km::Vector{Float64})
    n = length(time_s)
    n < 2 && return NaN
    t_day = time_s ./ 86400.0
    t_mean = mean(t_day)
    y_mean = mean(pos_norm_km)
    denom = sum((t_day .- t_mean) .^ 2)
    denom <= eps(Float64) && return NaN
    return sum((t_day .- t_mean) .* (pos_norm_km .- y_mean)) / denom
end

function _method_summary(method::String, time_s::Vector{Float64}, truth_state::Matrix{Float64}, state::Matrix{Float64}, truth_v_kmps::Matrix{Float64}; nis=nothing, accepted=nothing)
    err = state .- truth_state
    pos_err = err[:, 1:3]
    vel_err = err[:, 4:6]
    pos_norm = vec(sqrt.(sum(pos_err .^ 2; dims=2)))
    vel_norm = vec(sqrt.(sum(vel_err .^ 2; dims=2)))
    rtn = _rtn_errors_km(truth_state[:, 1:3], truth_v_kmps, pos_err)
    return (
        method=method,
        n_samples=length(time_s),
        duration_hr=(time_s[end] - time_s[1]) / 3600.0,
        pos_rmse_km=sqrt(mean(pos_norm .^ 2)),
        pos_median_km=median(pos_norm),
        pos_p95_km=_safe_quantile(pos_norm, 0.95),
        pos_max_km=maximum(pos_norm),
        pos_drift_rate_km_per_day=_drift_rate_km_per_day(time_s, pos_norm),
        x_rmse_km=sqrt(mean(pos_err[:, 1] .^ 2)),
        y_rmse_km=sqrt(mean(pos_err[:, 2] .^ 2)),
        z_rmse_km=sqrt(mean(pos_err[:, 3] .^ 2)),
        radial_rmse_km=sqrt(mean(rtn[:, 1] .^ 2)),
        transverse_rmse_km=sqrt(mean(rtn[:, 2] .^ 2)),
        normal_rmse_km=sqrt(mean(rtn[:, 3] .^ 2)),
        vel_rmse_kmps=sqrt(mean(vel_norm .^ 2)),
        nis_mean=nis === nothing ? NaN : mean(nis),
        nis_median=nis === nothing ? NaN : median(nis),
        measurement_updates_accepted=accepted === nothing ? missing : count(identity, accepted),
        measurement_updates_rejected=accepted === nothing ? missing : count(!, accepted),
        measurement_rejection_fraction=accepted === nothing ? NaN : count(!, accepted) / length(accepted),
        timestamp_utc=string(now(UTC))
    )
end

function _method_errors(method::String, time_s::Vector{Float64}, truth_state::Matrix{Float64}, state::Matrix{Float64}, truth_v_kmps::Matrix{Float64}; innovations=nothing, nis=nothing, accepted=nothing)
    err = state .- truth_state
    pos_err = err[:, 1:3]
    vel_err = err[:, 4:6]
    pos_norm = vec(sqrt.(sum(pos_err .^ 2; dims=2)))
    vel_norm = vec(sqrt.(sum(vel_err .^ 2; dims=2)))
    rtn = _rtn_errors_km(truth_state[:, 1:3], truth_v_kmps, pos_err)
    n = length(time_s)
    df = DataFrame(
        method=fill(method, n),
        idx=collect(1:n),
        time_s=time_s,
        truth_x_km=truth_state[:, 1],
        truth_y_km=truth_state[:, 2],
        truth_z_km=truth_state[:, 3],
        truth_vx_kmps=truth_state[:, 4],
        truth_vy_kmps=truth_state[:, 5],
        truth_vz_kmps=truth_state[:, 6],
        state_x_km=state[:, 1],
        state_y_km=state[:, 2],
        state_z_km=state[:, 3],
        state_vx_kmps=state[:, 4],
        state_vy_kmps=state[:, 5],
        state_vz_kmps=state[:, 6],
        error_x_km=pos_err[:, 1],
        error_y_km=pos_err[:, 2],
        error_z_km=pos_err[:, 3],
        error_vx_kmps=vel_err[:, 1],
        error_vy_kmps=vel_err[:, 2],
        error_vz_kmps=vel_err[:, 3],
        position_error_km=pos_norm,
        velocity_error_kmps=vel_norm,
        radial_error_km=rtn[:, 1],
        transverse_error_km=rtn[:, 2],
        normal_error_km=rtn[:, 3]
    )
    if innovations !== nothing
        df.innovation_x_km = innovations[:, 1]
        df.innovation_y_km = innovations[:, 2]
        df.innovation_z_km = innovations[:, 3]
        if size(innovations, 2) == 6
            df.innovation_vx_kmps = innovations[:, 4]
            df.innovation_vy_kmps = innovations[:, 5]
            df.innovation_vz_kmps = innovations[:, 6]
        end
    end
    if nis !== nothing
        df.nis = nis
    end
    if accepted !== nothing
        df.measurement_update_accepted = accepted
        df.measurement_update_rejected = .!accepted
    end
    return df
end

function _run_filter_case(req::CYGNSSEstimatorValidationRequest, name::String, time_s, sim_state, measurement_state; include_velocity::Bool, multiplier::Float64=1.0)
    return _run_error_state_filter(
        time_s,
        sim_state,
        measurement_state;
        position_noise_m=req.position_noise_m * multiplier,
        velocity_noise_mps=req.velocity_noise_mps * multiplier,
        process_accel_noise_mps2=req.process_accel_noise_mps2 * multiplier,
        initial_position_sigma_m=req.initial_position_sigma_m * multiplier,
        initial_velocity_sigma_mps=req.initial_velocity_sigma_mps * multiplier,
        innovation_gate_sigma=req.innovation_gate_sigma,
        position_innovation_gate_km=req.position_innovation_gate_km,
        velocity_innovation_gate_mps=req.velocity_innovation_gate_mps,
        include_velocity_measurements=include_velocity
    )
end

function _cygnss_plots_module()
    return Plots
end

function _save_cygnss_summary_metric_plot(summary::DataFrame, outdir::String)
    Plots = _cygnss_plots_module()
    methods = String.(summary.method)
    x = collect(1:length(methods))
    p = Plots.plot(
        title="CYGNSS 48hr Truth-Referenced Summary Metrics",
        xlabel="method",
        ylabel="km / km per day",
        size=(1700, 950),
        dpi=150,
        legend=:topright,
        xrotation=25
    )
    Plots.plot!(p, x, Float64.(summary.pos_rmse_km); label="position RMSE [km]", marker=:circle, linewidth=2)
    Plots.plot!(p, x, Float64.(summary.pos_p95_km); label="position p95 [km]", marker=:diamond, linewidth=2)
    Plots.plot!(p, x, Float64.(summary.pos_max_km); label="position max [km]", marker=:star5, linewidth=2)
    Plots.plot!(p, x, abs.(Float64.(summary.pos_drift_rate_km_per_day)); label="|drift rate| [km/day]", marker=:utriangle, linewidth=2)
    Plots.plot!(p; xticks=(x, methods), gridalpha=0.25)
    out = joinpath(outdir, "summary_position_metrics.png")
    Plots.savefig(p, out)
    return out
end

function _save_cygnss_rtn_metric_plot(summary::DataFrame, outdir::String)
    Plots = _cygnss_plots_module()
    methods = String.(summary.method)
    x = collect(1:length(methods))
    p = Plots.plot(
        title="CYGNSS 48hr RTN Position RMSE",
        xlabel="method",
        ylabel="RMSE [km]",
        size=(1700, 950),
        dpi=150,
        legend=:topright,
        xrotation=25
    )
    Plots.plot!(p, x, Float64.(summary.radial_rmse_km); label="radial", marker=:circle, linewidth=2)
    Plots.plot!(p, x, Float64.(summary.transverse_rmse_km); label="transverse", marker=:diamond, linewidth=2)
    Plots.plot!(p, x, Float64.(summary.normal_rmse_km); label="normal", marker=:star5, linewidth=2)
    Plots.plot!(p; xticks=(x, methods), gridalpha=0.25)
    out = joinpath(outdir, "summary_rtn_rmse.png")
    Plots.savefig(p, out)
    return out
end

function _save_cygnss_velocity_nis_metric_plot(summary::DataFrame, outdir::String)
    Plots = _cygnss_plots_module()
    methods = String.(summary.method)
    x = collect(1:length(methods))
    p1 = Plots.plot(
        x,
        Float64.(summary.vel_rmse_kmps) .* 1.0e3;
        title="Velocity Error RMSE",
        xlabel="method",
        ylabel="m/s",
        marker=:circle,
        linewidth=2,
        legend=false,
        xrotation=25,
        xticks=(x, methods),
        gridalpha=0.25
    )
    nis = Float64.(summary.nis_mean)
    finite_idx = findall(isfinite, nis)
    p2 = Plots.plot(
        title="Mean Normalized Innovation Squared",
        xlabel="method",
        ylabel="NIS",
        legend=false,
        xrotation=25,
        gridalpha=0.25
    )
    if !isempty(finite_idx)
        Plots.scatter!(p2, finite_idx, nis[finite_idx]; marker=:diamond, markersize=6)
    end
    Plots.plot!(p2; xticks=(x, methods))
    p = Plots.plot(p1, p2; layout=(2, 1), size=(1700, 1200), dpi=150)
    out = joinpath(outdir, "summary_velocity_and_nis.png")
    Plots.savefig(p, out)
    return out
end

function _save_cygnss_position_error_timeseries(errors::DataFrame, outdir::String)
    Plots = _cygnss_plots_module()
    p = Plots.plot(
        title="CYGNSS 48hr Position Error Norm vs 96hr Truth",
        xlabel="time [hr]",
        ylabel="position error [km]",
        size=(1700, 950),
        dpi=150,
        legend=:topright,
        gridalpha=0.25
    )
    for method in unique(String.(errors.method))
        sub = errors[errors.method .== method, :]
        Plots.plot!(p, Float64.(sub.time_s) ./ 3600.0, Float64.(sub.position_error_km); label=method, linewidth=1.8)
    end
    out = joinpath(outdir, "position_error_timeseries.png")
    Plots.savefig(p, out)
    return out
end

function _save_cygnss_velocity_error_timeseries(errors::DataFrame, outdir::String)
    Plots = _cygnss_plots_module()
    p = Plots.plot(
        title="CYGNSS 48hr Velocity Error Norm vs 96hr Truth",
        xlabel="time [hr]",
        ylabel="velocity error [m/s]",
        size=(1700, 950),
        dpi=150,
        legend=:topright,
        gridalpha=0.25
    )
    for method in unique(String.(errors.method))
        sub = errors[errors.method .== method, :]
        Plots.plot!(p, Float64.(sub.time_s) ./ 3600.0, Float64.(sub.velocity_error_kmps) .* 1.0e3; label=method, linewidth=1.8)
    end
    out = joinpath(outdir, "velocity_error_timeseries.png")
    Plots.savefig(p, out)
    return out
end

function _save_cygnss_cartesian_error_timeseries(errors::DataFrame, outdir::String)
    Plots = _cygnss_plots_module()
    panels = []
    for (col, label) in ((:error_x_km, "x"), (:error_y_km, "y"), (:error_z_km, "z"))
        p = Plots.plot(title="Cartesian $label error", xlabel="time [hr]", ylabel="error [km]", legend=:topright, gridalpha=0.25)
        for method in unique(String.(errors.method))
            sub = errors[errors.method .== method, :]
            Plots.plot!(p, Float64.(sub.time_s) ./ 3600.0, Float64.(sub[!, col]); label=method, linewidth=1.5)
        end
        push!(panels, p)
    end
    p_all = Plots.plot(panels...; layout=(3, 1), size=(1700, 1300), dpi=150)
    out = joinpath(outdir, "cartesian_error_timeseries.png")
    Plots.savefig(p_all, out)
    return out
end

function _save_cygnss_rtn_error_timeseries(errors::DataFrame, outdir::String)
    Plots = _cygnss_plots_module()
    panels = []
    for (col, label) in ((:radial_error_km, "radial"), (:transverse_error_km, "transverse"), (:normal_error_km, "normal"))
        p = Plots.plot(title="RTN $label error", xlabel="time [hr]", ylabel="error [km]", legend=:topright, gridalpha=0.25)
        for method in unique(String.(errors.method))
            sub = errors[errors.method .== method, :]
            Plots.plot!(p, Float64.(sub.time_s) ./ 3600.0, Float64.(sub[!, col]); label=method, linewidth=1.5)
        end
        push!(panels, p)
    end
    p_all = Plots.plot(panels...; layout=(3, 1), size=(1700, 1300), dpi=150)
    out = joinpath(outdir, "rtn_error_timeseries.png")
    Plots.savefig(p_all, out)
    return out
end

function _save_cygnss_innovation_plots(errors::DataFrame, outdir::String)
    Plots = _cygnss_plots_module()
    estimator_rows = errors[startswith.(String.(errors.method), "estimator"), :]
    if nrow(estimator_rows) == 0 || !hasproperty(estimator_rows, :innovation_x_km)
        return String[]
    end

    p = Plots.plot(
        title="Estimator Position Innovations",
        xlabel="time [hr]",
        ylabel="innovation [km]",
        size=(1700, 950),
        dpi=150,
        legend=:topright,
        gridalpha=0.25
    )
    for method in unique(String.(estimator_rows.method))
        sub = estimator_rows[estimator_rows.method .== method, :]
        Plots.plot!(p, Float64.(sub.time_s) ./ 3600.0, Float64.(sub.innovation_x_km); label="$method x", linewidth=1.2)
        Plots.plot!(p, Float64.(sub.time_s) ./ 3600.0, Float64.(sub.innovation_y_km); label="$method y", linewidth=1.2)
        Plots.plot!(p, Float64.(sub.time_s) ./ 3600.0, Float64.(sub.innovation_z_km); label="$method z", linewidth=1.2)
    end
    innovation_out = joinpath(outdir, "position_innovations.png")
    Plots.savefig(p, innovation_out)

    p_nis = Plots.plot(
        title="Estimator Normalized Innovation Squared",
        xlabel="time [hr]",
        ylabel="NIS",
        size=(1700, 950),
        dpi=150,
        legend=:topright,
        gridalpha=0.25
    )
    for method in unique(String.(estimator_rows.method))
        sub = estimator_rows[(estimator_rows.method .== method) .& isfinite.(Float64.(estimator_rows.nis)), :]
        nrow(sub) == 0 && continue
        Plots.plot!(p_nis, Float64.(sub.time_s) ./ 3600.0, Float64.(sub.nis); label=method, linewidth=1.8)
    end
    nis_out = joinpath(outdir, "nis_timeseries.png")
    Plots.savefig(p_nis, nis_out)
    return [innovation_out, nis_out]
end

function _generate_cygnss_estimator_plots(summary::DataFrame, errors::DataFrame, outdir::String)
    mkpath(outdir)
    paths = String[]
    push!(paths, _save_cygnss_summary_metric_plot(summary, outdir))
    push!(paths, _save_cygnss_rtn_metric_plot(summary, outdir))
    push!(paths, _save_cygnss_velocity_nis_metric_plot(summary, outdir))
    push!(paths, _save_cygnss_position_error_timeseries(errors, outdir))
    push!(paths, _save_cygnss_velocity_error_timeseries(errors, outdir))
    push!(paths, _save_cygnss_cartesian_error_timeseries(errors, outdir))
    push!(paths, _save_cygnss_rtn_error_timeseries(errors, outdir))
    append!(paths, _save_cygnss_innovation_plots(errors, outdir))
    return paths
end

function _cygnss_offset_values(min_offset_s::Float64, max_offset_s::Float64, step_s::Float64)
    step_s > 0.0 || throw(ArgumentError("Offset sweep step_s must be > 0, got $step_s."))
    min_offset_s <= max_offset_s || throw(ArgumentError("Offset sweep min_offset_s must be <= max_offset_s."))
    vals = collect(min_offset_s:step_s:max_offset_s)
    if isempty(vals) || vals[end] < max_offset_s - 1.0e-9
        push!(vals, max_offset_s)
    end
    return Float64.(vals)
end

function _offset_metrics(offset_s::Float64, truth, measurements)
    truth_time = measurements.time_s .+ offset_s
    truth_r = _interp_matrix_columns(truth.time_s, truth.r_km, truth_time)
    truth_v = _interp_matrix_columns(truth.time_s, truth.v_kmps, truth_time)
    err_r = measurements.r_km .- truth_r
    err_v = measurements.v_kmps .- truth_v
    pos_norm = vec(sqrt.(sum(err_r .^ 2; dims=2)))
    vel_norm = vec(sqrt.(sum(err_v .^ 2; dims=2)))
    rtn = _rtn_errors_km(truth_r, truth_v, err_r)
    return (
        offset_s=offset_s,
        initial_position_error_km=pos_norm[1],
        initial_radial_error_km=rtn[1, 1],
        initial_transverse_error_km=rtn[1, 2],
        initial_normal_error_km=rtn[1, 3],
        position_rmse_km=sqrt(mean(pos_norm .^ 2)),
        position_median_km=median(pos_norm),
        position_p95_km=quantile(pos_norm, 0.95),
        position_max_km=maximum(pos_norm),
        velocity_rmse_mps=sqrt(mean(vel_norm .^ 2)) * 1.0e3,
        radial_rmse_km=sqrt(mean(rtn[:, 1] .^ 2)),
        transverse_rmse_km=sqrt(mean(rtn[:, 2] .^ 2)),
        normal_rmse_km=sqrt(mean(rtn[:, 3] .^ 2))
    )
end

function _save_cygnss_offset_sweep_plot(summary::DataFrame, outdir::String)
    Plots = _cygnss_plots_module()
    mkpath(outdir)
    p1 = Plots.plot(
        Float64.(summary.offset_s),
        Float64.(summary.position_rmse_km);
        title="CYGNSS Timing Offset Sweep",
        xlabel="truth time offset [s]",
        ylabel="raw GPS position RMSE [km]",
        label="position RMSE",
        marker=:circle,
        linewidth=2,
        gridalpha=0.25
    )
    best_rmse_idx = argmin(Float64.(summary.position_rmse_km))
    Plots.vline!(p1, [Float64(summary.offset_s[best_rmse_idx])]; label="best RMSE", linestyle=:dash)

    p2 = Plots.plot(
        Float64.(summary.offset_s),
        Float64.(summary.initial_position_error_km);
        xlabel="truth time offset [s]",
        ylabel="initial raw GPS error [km]",
        label="initial error",
        marker=:diamond,
        linewidth=2,
        gridalpha=0.25
    )
    best_initial_idx = argmin(Float64.(summary.initial_position_error_km))
    Plots.vline!(p2, [Float64(summary.offset_s[best_initial_idx])]; label="best initial", linestyle=:dash)

    p3 = Plots.plot(
        Float64.(summary.offset_s),
        Float64.(summary.radial_rmse_km);
        xlabel="truth time offset [s]",
        ylabel="RTN RMSE [km]",
        label="radial",
        linewidth=2,
        gridalpha=0.25
    )
    Plots.plot!(p3, Float64.(summary.offset_s), Float64.(summary.transverse_rmse_km); label="transverse", linewidth=2)
    Plots.plot!(p3, Float64.(summary.offset_s), Float64.(summary.normal_rmse_km); label="normal", linewidth=2)

    p = Plots.plot(p1, p2, p3; layout=(3, 1), size=(1700, 1500), dpi=150)
    out = joinpath(outdir, "timing_offset_sweep.png")
    Plots.savefig(p, out)
    return out
end

function run_cygnss_timing_offset_sweep(req::CYGNSSTimingOffsetSweepRequest=CYGNSSTimingOffsetSweepRequest())::CYGNSSTimingOffsetSweepResult
    truth = _load_cygnss_truth_series(CYGNSSEstimatorValidationRequest(truth_path=req.truth_path))
    measurements = _load_cygnss_measurement_series(CYGNSSEstimatorValidationRequest(measurement_path=req.measurement_path, max_points=req.max_points))
    rows = [_offset_metrics(offset_s, truth, measurements) for offset_s in _cygnss_offset_values(req.min_offset_s, req.max_offset_s, req.step_s)]
    summary = DataFrame(rows)
    mkpath(dirname(req.out_summary))
    CSV.write(req.out_summary, summary)
    plot_path = req.generate_plots ? _save_cygnss_offset_sweep_plot(summary, req.plots_dir) : ""

    best_rmse_idx = argmin(Float64.(summary.position_rmse_km))
    best_initial_idx = argmin(Float64.(summary.initial_position_error_km))
    println("CYGNSS timing offset sweep")
    println("truth=$(req.truth_path)")
    println("measurements=$(req.measurement_path)")
    println("offset_range_s=$(req.min_offset_s):$(req.step_s):$(req.max_offset_s)")
    println("summary=$(req.out_summary)")
    !isempty(plot_path) && println("plot=$(plot_path)")
    println("best raw GPS RMSE offset [s] = $(summary.offset_s[best_rmse_idx])")
    println("best raw GPS initial-error offset [s] = $(summary.offset_s[best_initial_idx])")

    return CYGNSSTimingOffsetSweepResult(
        summary=summary,
        summary_path=req.out_summary,
        plots_dir=req.generate_plots ? req.plots_dir : "",
        best_raw_gps_rmse_offset_s=Float64(summary.offset_s[best_rmse_idx]),
        best_raw_gps_initial_offset_s=Float64(summary.offset_s[best_initial_idx])
    )
end

function run_cygnss_estimator_validation(req::CYGNSSEstimatorValidationRequest=CYGNSSEstimatorValidationRequest())::CYGNSSEstimatorValidationResult
    series = _load_cygnss_truth_and_measurements(req)
    open_loop = _run_cygnss_open_loop(req, series)
    n = length(series.time_s)
    sim_r_km = _extract_open_loop_state(open_loop.errors, n)
    sim_state = _state_from_position_series(sim_r_km, series.time_s)
    truth_state = hcat(series.truth_r_km, series.truth_v_kmps)
    measurement_state = hcat(series.measurement_r_km, series.measurement_v_kmps)

    summary_rows = NamedTuple[]
    error_tables = DataFrame[]
    push!(summary_rows, _method_summary("sim_open_loop", series.time_s, truth_state, sim_state, series.truth_v_kmps))
    push!(error_tables, _method_errors("sim_open_loop", series.time_s, truth_state, sim_state, series.truth_v_kmps))
    push!(summary_rows, _method_summary("raw_gps", series.time_s, truth_state, measurement_state, series.truth_v_kmps))
    push!(error_tables, _method_errors("raw_gps", series.time_s, truth_state, measurement_state, series.truth_v_kmps))

    pos_filter = _run_filter_case(req, "estimator_position", series.time_s, sim_state, measurement_state; include_velocity=false)
    push!(summary_rows, _method_summary("estimator_position", series.time_s, truth_state, pos_filter.state, series.truth_v_kmps; nis=pos_filter.nis, accepted=pos_filter.accepted))
    push!(error_tables, _method_errors("estimator_position", series.time_s, truth_state, pos_filter.state, series.truth_v_kmps; innovations=pos_filter.innovations, nis=pos_filter.nis, accepted=pos_filter.accepted))

    if req.run_posvel_case || req.include_velocity_measurements
        posvel_filter = _run_filter_case(req, "estimator_posvel", series.time_s, sim_state, measurement_state; include_velocity=true)
        push!(summary_rows, _method_summary("estimator_posvel", series.time_s, truth_state, posvel_filter.state, series.truth_v_kmps; nis=posvel_filter.nis, accepted=posvel_filter.accepted))
        push!(error_tables, _method_errors("estimator_posvel", series.time_s, truth_state, posvel_filter.state, series.truth_v_kmps; innovations=posvel_filter.innovations, nis=posvel_filter.nis, accepted=posvel_filter.accepted))
    end

    if req.run_sensitivity
        for multiplier in (0.1, 1.0, 10.0)
            filt = _run_filter_case(req, "estimator_position_sweep_x$(multiplier)", series.time_s, sim_state, measurement_state; include_velocity=false, multiplier=multiplier)
            method = "estimator_position_sweep_x$(multiplier)"
            push!(summary_rows, _method_summary(method, series.time_s, truth_state, filt.state, series.truth_v_kmps; nis=filt.nis, accepted=filt.accepted))
            push!(error_tables, _method_errors(method, series.time_s, truth_state, filt.state, series.truth_v_kmps; innovations=filt.innovations, nis=filt.nis, accepted=filt.accepted))
        end
    end

    summary = DataFrame(summary_rows)
    errors = vcat(error_tables...; cols=:union)
    mkpath(dirname(req.out_summary))
    mkpath(dirname(req.out_errors))
    CSV.write(req.out_summary, summary)
    CSV.write(req.out_errors, errors)
    plot_paths = req.generate_plots ? _generate_cygnss_estimator_plots(summary, errors, req.plots_dir) : String[]

    println("CYGNSS estimator validation")
    println("truth=$(req.truth_path)")
    println("measurements=$(req.measurement_path)")
    println("summary=$(req.out_summary)")
    println("errors=$(req.out_errors)")
    println("plots=$(req.generate_plots ? req.plots_dir : "disabled")")
    show(summary, allrows=true, allcols=true)
    if !isempty(plot_paths)
        println("\nPlot files:")
        foreach(path -> println("  $path"), plot_paths)
    end
    println()

    return CYGNSSEstimatorValidationResult(
        summary=summary,
        errors=errors,
        summary_path=req.out_summary,
        errors_path=req.out_errors,
        plots_dir=req.generate_plots ? req.plots_dir : ""
    )
end

function _cygnss_estimator_parse_cli(args::Vector{String})
    vals = Dict{String, String}()
    flags = Set{String}()
    for arg in args
        if startswith(arg, "--") && occursin("=", arg)
            key, val = split(arg[3:end], "=", limit=2)
            vals[key] = val
        elseif startswith(arg, "--")
            push!(flags, arg[3:end])
        else
            throw(ArgumentError("Unsupported CYGNSS estimator validation argument '$arg'. Use --key=value flags."))
        end
    end
    bool_value(key, default) = haskey(vals, key) ? _safe_parse_bool(vals[key], default) : (key in flags ? true : default)
    sym_value(key, default) = Symbol(lowercase(get(vals, key, String(default))))
    return CYGNSSEstimatorValidationRequest(
        truth_path=abspath(get(vals, "truth", joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS", "cyg04_nasa_pvt_96hr.feather"))),
        measurement_path=abspath(get(vals, "measurements", joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS", "cygnss_data_48hr.feather"))),
        ic_path=abspath(get(vals, "ic", joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS", "cygnss_data_48hr.feather"))),
        gravity_harmonics_file=abspath(get(vals, "gravity-harmonics", joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv"))),
        out_summary=abspath(get(vals, "out-summary", joinpath(DEFAULT_OUTPUT_DIR, "cygnss_estimator_validation_summary.csv"))),
        out_errors=abspath(get(vals, "out-errors", joinpath(DEFAULT_OUTPUT_DIR, "cygnss_estimator_validation_errors.csv"))),
        plots_dir=abspath(get(vals, "plots-dir", joinpath(DEFAULT_OUTPUT_DIR, "cygnss_estimator_validation_plots"))),
        profile=sym_value("profile", :quick),
        max_points=parse(Int, get(vals, "max-points", "0")),
        position_noise_m=parse(Float64, get(vals, "position-noise-m", "100.0")),
        velocity_noise_mps=parse(Float64, get(vals, "velocity-noise-mps", "10.0")),
        process_accel_noise_mps2=parse(Float64, get(vals, "process-accel-noise-mps2", "1e-5")),
        initial_position_sigma_m=parse(Float64, get(vals, "initial-position-sigma-m", "1000.0")),
        initial_velocity_sigma_mps=parse(Float64, get(vals, "initial-velocity-sigma-mps", "10.0")),
        innovation_gate_sigma=parse(Float64, get(vals, "innovation-gate-sigma", "0.0")),
        position_innovation_gate_km=parse(Float64, get(vals, "position-innovation-gate-km", "25.0")),
        velocity_innovation_gate_mps=parse(Float64, get(vals, "velocity-innovation-gate-mps", "100.0")),
        truth_time_offset_s=parse(Float64, get(vals, "truth-time-offset-s", "0.0")),
        include_velocity_measurements=bool_value("include-velocity-measurements", false),
        run_posvel_case=bool_value("run-posvel-case", true),
        run_sensitivity=bool_value("sensitivity", false),
        generate_plots=bool_value("plots", true)
    )
end

function run_cygnss_estimator_validation_cli(args::Vector{String}=copy(ARGS))::CYGNSSEstimatorValidationResult
    if any(arg -> startswith(arg, "--offset-sweep"), args)
        vals = Dict{String, String}()
        flags = Set{String}()
        for arg in args
            if startswith(arg, "--") && occursin("=", arg)
                key, val = split(arg[3:end], "=", limit=2)
                vals[key] = val
            elseif startswith(arg, "--")
                push!(flags, arg[3:end])
            end
        end
        bool_value(key, default) = haskey(vals, key) ? _safe_parse_bool(vals[key], default) : (key in flags ? true : default)
        sweep_result = run_cygnss_timing_offset_sweep(CYGNSSTimingOffsetSweepRequest(
            truth_path=abspath(get(vals, "truth", joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS", "cyg04_nasa_pvt_96hr.feather"))),
            measurement_path=abspath(get(vals, "measurements", joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS", "cygnss_data_48hr.feather"))),
            out_summary=abspath(get(vals, "offset-summary", joinpath(DEFAULT_OUTPUT_DIR, "cygnss_timing_offset_sweep.csv"))),
            plots_dir=abspath(get(vals, "plots-dir", joinpath(DEFAULT_OUTPUT_DIR, "cygnss_estimator_validation_plots"))),
            min_offset_s=parse(Float64, get(vals, "offset-min-s", "-120.0")),
            max_offset_s=parse(Float64, get(vals, "offset-max-s", "120.0")),
            step_s=parse(Float64, get(vals, "offset-step-s", "1.0")),
            max_points=parse(Int, get(vals, "max-points", "0")),
            generate_plots=bool_value("plots", true)
        ))
        if bool_value("offset-sweep-only", false)
            empty = DataFrame()
            return CYGNSSEstimatorValidationResult(
                summary=empty,
                errors=empty,
                summary_path=sweep_result.summary_path,
                errors_path="",
                plots_dir=sweep_result.plots_dir
            )
        end
        if bool_value("apply-best-offset", false) && !haskey(vals, "truth-time-offset-s")
            push!(args, "--truth-time-offset-s=$(sweep_result.best_raw_gps_rmse_offset_s)")
        end
    end
    return run_cygnss_estimator_validation(_cygnss_estimator_parse_cli(args))
end

function _cygnss_timing_offset_parse_cli(args::Vector{String})
    vals = Dict{String, String}()
    flags = Set{String}()
    for arg in args
        if startswith(arg, "--") && occursin("=", arg)
            key, val = split(arg[3:end], "=", limit=2)
            vals[key] = val
        elseif startswith(arg, "--")
            push!(flags, arg[3:end])
        else
            throw(ArgumentError("Unsupported CYGNSS timing offset sweep argument '$arg'. Use --key=value flags."))
        end
    end
    bool_value(key, default) = haskey(vals, key) ? _safe_parse_bool(vals[key], default) : (key in flags ? true : default)
    return CYGNSSTimingOffsetSweepRequest(
        truth_path=abspath(get(vals, "truth", joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS", "cyg04_nasa_pvt_96hr.feather"))),
        measurement_path=abspath(get(vals, "measurements", joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS", "cygnss_data_48hr.feather"))),
        out_summary=abspath(get(vals, "out-summary", joinpath(DEFAULT_OUTPUT_DIR, "cygnss_timing_offset_sweep.csv"))),
        plots_dir=abspath(get(vals, "plots-dir", joinpath(DEFAULT_OUTPUT_DIR, "cygnss_estimator_validation_plots"))),
        min_offset_s=parse(Float64, get(vals, "offset-min-s", "-120.0")),
        max_offset_s=parse(Float64, get(vals, "offset-max-s", "120.0")),
        step_s=parse(Float64, get(vals, "offset-step-s", "1.0")),
        max_points=parse(Int, get(vals, "max-points", "0")),
        generate_plots=bool_value("plots", true)
    )
end

function run_cygnss_timing_offset_sweep_cli(args::Vector{String}=copy(ARGS))::CYGNSSTimingOffsetSweepResult
    return run_cygnss_timing_offset_sweep(_cygnss_timing_offset_parse_cli(args))
end
