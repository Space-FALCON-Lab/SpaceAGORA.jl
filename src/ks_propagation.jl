function ks_satellite_rhs_sundman_time!(
    du_sat::AbstractVector{<:Real},
    p::AbstractVector{<:Real},
    q::AbstractVector{<:Real},
    h::Real,
    accel::AbstractVector{<:Real},
)
    rho = ks_radius_scalar(p)
    rho > 0.0 || throw(ArgumentError("KS radius must be positive."))

    v = ks_velocity(p, q)
    s_mat = ks_s_matrix(p)
    accel4 = ks_augmented_acceleration(accel)

    du_sat[1:4] .= Float64.(q)
    du_sat[5:8] .= (-0.25 * float(h)) .* Float64.(p) .+ 0.5 .* rho .* (s_mat * accel4)
    du_sat[9] = -2.0 * rho * dot(v, Float64.(accel))
    du_sat[10] = rho
    return nothing
end

function ks_constellation_rhs_with_pair_data(
    u::AbstractVector{<:Real},
    pair_data::AbstractVector;
    mu::Real = MU_EARTH_KM3_S2,
    j2::Real = J2_EARTH,
    earth_radius_km::Real = R_EARTH_KM,
    cr::Real = 2.0,
    laser_power_w::Real = 1.0e4,
    satellite_mass_kg::Real = 300.0,
    control_amplitudes::AbstractVector{<:Real},
)
    nsat = ks_satellite_count(u)
    du = zeros(Float64, length(u))
    laser_accels = ks_total_laser_accelerations(
        pair_data,
        nsat;
        control_amplitudes = control_amplitudes,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
    )

    for sat in 1:nsat
        p_sat, q_sat, h_sat, _ = ks_state_components(u, sat)
        r_sat = ks_position(p_sat)
        total_accel = ks_j2_acceleration(
            r_sat;
            mu = mu,
            j2 = j2,
            earth_radius_km = earth_radius_km,
        ) + laser_accels[:, sat]

        ks_satellite_rhs_sundman_time!(
            view(du, ks_state_slice(sat)),
            p_sat,
            q_sat,
            h_sat,
            total_accel,
        )
    end

    return du
end

function ks_constellation_dynamics!(du, u, p, t)
    nsat = p.nsat
    pair_data = ks_pair_geometry_data(
        u;
        pair_list = p.pair_list,
        link_max_range_km = p.link_max_range_km,
        atmosphere_top_km = p.atmosphere_top_km,
    )

    du .= ks_constellation_rhs_with_pair_data(
        u,
        pair_data;
        mu = p.mu,
        j2 = p.j2,
        earth_radius_km = p.earth_radius_km,
        cr = p.cr,
        laser_power_w = p.laser_power_w,
        satellite_mass_kg = p.satellite_mass_kg,
        control_amplitudes = p.control_amplitudes,
    )
    return nothing
end

function propagate_ks_constellation(
    u0::AbstractVector{<:Real},
    physical_tspan::Tuple{<:Real, <:Real};
    mu::Real = MU_EARTH_KM3_S2,
    j2::Real = J2_EARTH,
    earth_radius_km::Real = R_EARTH_KM,
    link_max_range_km::Real = 200.0,
    atmosphere_top_km::Real = 100.0,
    cr::Real = 2.0,
    laser_power_w::Real = 1.0e4,
    satellite_mass_kg::Real = 300.0,
    control_amplitudes::Union{Nothing, AbstractVector{<:Real}} = nothing,
    saveat_sundman::Union{Nothing, Real} = nothing,
    reltol::Real = 1e-9,
    abstol::Real = 1e-9,
    solver = Tsit5(),
)
    nsat = ks_satellite_count(u0)
    pair_list = ks_pair_indices(nsat)
    default_controls = ones(Float64, nsat)
    controls = isnothing(control_amplitudes) ? default_controls : Float64.(control_amplitudes)
    length(controls) == nsat ||
        throw(ArgumentError("control_amplitudes must have one entry per satellite."))

    t0 = float(physical_tspan[1])
    tf = float(physical_tspan[2])
    tf >= t0 || throw(ArgumentError("physical_tspan must satisfy tf >= t0."))
    initial_min_radius = minimum(ks_radius_scalar(ks_state_components(u0, sat)[1]) for sat in 1:nsat)
    initial_min_radius > 0.0 || throw(ArgumentError("Initial KS radius must be positive."))
    s_upper = max((tf - t0) / initial_min_radius * 2.5, 1e-6)
    saveat_val = isnothing(saveat_sundman) ? max((tf - t0) / initial_min_radius / 100.0, 1e-4) : float(saveat_sundman)

    p = (
        nsat = nsat,
        pair_list = pair_list,
        mu = float(mu),
        j2 = float(j2),
        earth_radius_km = float(earth_radius_km),
        link_max_range_km = float(link_max_range_km),
        atmosphere_top_km = float(atmosphere_top_km),
        cr = float(cr),
        laser_power_w = float(laser_power_w),
        satellite_mass_kg = float(satellite_mass_kg),
        control_amplitudes = controls,
    )

    condition(u, s, integrator) = ks_min_physical_time(u) - tf
    affect!(integrator) = terminate!(integrator)
    callback = ContinuousCallback(condition, affect!; save_positions = (true, true))

    prob = ODEProblem(ks_constellation_dynamics!, Float64.(u0), (0.0, s_upper), p)
    return solve(
        prob,
        solver;
        saveat = saveat_val,
        reltol = reltol,
        abstol = abstol,
        callback = callback,
    )
end
