function ks_velocity_jacobians(
    p::AbstractVector{<:Real},
    q::AbstractVector{<:Real},
)
    rho = ks_radius_scalar(p)
    rho > 0.0 || throw(ArgumentError("KS radius must be positive."))

    lambda_p = ks_lambda_matrix(p)
    lambda_q = ks_lambda_matrix(q)
    vel_from_pq = lambda_p * Float64.(q)

    dv_dp = (2.0 / rho) .* lambda_q .- (4.0 / rho^2) .* (vel_from_pq * transpose(Float64.(p)))
    dv_dq = (2.0 / rho) .* lambda_p
    return dv_dp, dv_dq
end

function ks_force_partials_wrt_p(
    p::AbstractVector{<:Real},
    q::AbstractVector{<:Real},
    jr::AbstractMatrix{<:Real},
    jv::AbstractMatrix{<:Real},
)
    dr_dp = 2.0 .* ks_lambda_matrix(p)
    dv_dp, _ = ks_velocity_jacobians(p, q)
    return Matrix{Float64}(jr) * dr_dp + Matrix{Float64}(jv) * dv_dp
end

function ks_force_partials_wrt_q(
    p::AbstractVector{<:Real},
    q::AbstractVector{<:Real},
    jv::AbstractMatrix{<:Real},
)
    _, dv_dq = ks_velocity_jacobians(p, q)
    return Matrix{Float64}(jv) * dv_dq
end

function ks_sa_derivative_columns(
    p::AbstractVector{<:Real},
    q::AbstractVector{<:Real},
    accel::AbstractVector{<:Real},
    jr::AbstractMatrix{<:Real},
    jv::AbstractMatrix{<:Real},
)
    s_mat = ks_s_matrix(p)
    accel4 = ks_augmented_acceleration(accel)
    dap_dp = ks_force_partials_wrt_p(p, q, jr, jv)

    dcols = zeros(Float64, 4, 4)
    for j in 1:4
        dap4 = zeros(Float64, 4)
        dap4[1:3] .= dap_dp[:, j]
        dcols[:, j] .= KS_DSP[j] * accel4 + s_mat * dap4
    end

    return dcols
end

function ks_satellite_continuous_jacobian_block(
    p::AbstractVector{<:Real},
    q::AbstractVector{<:Real},
    h::Real,
    accel::AbstractVector{<:Real};
    jr::AbstractMatrix{<:Real},
    jv::AbstractMatrix{<:Real} = zeros(Float64, 3, 3),
)
    rho = ks_radius_scalar(p)
    rho > 0.0 || throw(ArgumentError("KS radius must be positive."))

    p_vec = Float64.(p)
    q_vec = Float64.(q)
    accel_vec = Float64.(accel)
    v = ks_velocity(p_vec, q_vec)
    s_mat = ks_s_matrix(p_vec)
    s_a = s_mat * ks_augmented_acceleration(accel_vec)

    dv_dp, dv_dq = ks_velocity_jacobians(p_vec, q_vec)
    dap_dp = ks_force_partials_wrt_p(p_vec, q_vec, jr, jv)
    dap_dq = ks_force_partials_wrt_q(p_vec, q_vec, jv)
    dsa_dp = ks_sa_derivative_columns(p_vec, q_vec, accel_vec, jr, jv)

    a = zeros(Float64, 10, 10)
    a[1:4, 5:8] .= KS_I4

    a[5:8, 1:4] .= (-0.25 * float(h)) .* KS_I4 .+ s_a * transpose(0.5 .* (2.0 .* p_vec)) .+ 0.5 .* rho .* dsa_dp

    dapdq4 = zeros(Float64, 4, 4)
    dapdq4[1:3, :] .= dap_dq
    a[5:8, 5:8] .= 0.5 .* rho .* (s_mat * dapdq4)
    a[5:8, 9] .= -0.25 .* p_vec

    a[9, 1:4] .= -2.0 .* ((dot(v, accel_vec)) .* (2.0 .* p_vec) .+ rho .* (transpose(dv_dp) * accel_vec + transpose(dap_dp) * v))
    a[9, 5:8] .= -2.0 .* rho .* (transpose(dv_dq) * accel_vec + transpose(dap_dq) * v)
    a[10, 1:4] .= 2.0 .* p_vec

    return a
end

function ks_satellite_force_coupling_block(
    target_p::AbstractVector{<:Real},
    target_q::AbstractVector{<:Real},
    source_p::AbstractVector{<:Real},
    source_q::AbstractVector{<:Real};
    jr::AbstractMatrix{<:Real},
    jv::AbstractMatrix{<:Real} = zeros(Float64, 3, 3),
)
    rho_target = ks_radius_scalar(target_p)
    rho_target > 0.0 || throw(ArgumentError("KS radius must be positive."))

    s_target = ks_s_matrix(target_p)
    v_target = ks_velocity(target_p, target_q)
    dap_dp = ks_force_partials_wrt_p(source_p, source_q, jr, jv)
    dap_dq = ks_force_partials_wrt_q(source_p, source_q, jv)

    a = zeros(Float64, 10, 10)

    dapdp4 = zeros(Float64, 4, 4)
    dapdp4[1:3, :] .= dap_dp
    a[5:8, 1:4] .= 0.5 .* rho_target .* (s_target * dapdp4)
    a[9, 1:4] .= -2.0 .* rho_target .* (transpose(dap_dp) * v_target)

    dapdq4 = zeros(Float64, 4, 4)
    dapdq4[1:3, :] .= dap_dq
    a[5:8, 5:8] .= 0.5 .* rho_target .* (s_target * dapdq4)
    a[9, 5:8] .= -2.0 .* rho_target .* (transpose(dap_dq) * v_target)

    return a
end

function ks_continuous_jacobians(
    u::AbstractVector{<:Real};
    mu::Real = MU_EARTH_KM3_S2,
    j2::Real = J2_EARTH,
    earth_radius_km::Real = R_EARTH_KM,
    link_max_range_km::Real = 200.0,
    atmosphere_top_km::Real = 100.0,
    cr::Real = 2.0,
    laser_power_w::Real = 1.0e4,
    satellite_mass_kg::Real = 300.0,
    control_amplitudes::Union{Nothing, AbstractVector{<:Real}} = nothing,
    pair_list::AbstractVector{<:Tuple{Int, Int}} = ks_pair_indices(ks_satellite_count(u)),
)
    nsat = ks_satellite_count(u)
    nstate = 10 * nsat
    default_controls = ones(Float64, nsat)
    controls = isnothing(control_amplitudes) ? default_controls : Float64.(control_amplitudes)
    length(controls) == nsat ||
        throw(ArgumentError("control_amplitudes must have one entry per satellite."))

    pair_data = ks_pair_geometry_data(
        u;
        pair_list = pair_list,
        link_max_range_km = link_max_range_km,
        atmosphere_top_km = atmosphere_top_km,
    )
    ju_cart = ks_laser_cartesian_input_jacobian(
        pair_data,
        nsat;
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
    )
    laser_accels = ks_total_laser_accelerations(
        pair_data,
        nsat;
        control_amplitudes = controls,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
    )
    jr_laser_cart = ks_laser_cartesian_state_jacobian(
        pair_data,
        nsat;
        control_amplitudes = controls,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
    )

    a_c = zeros(Float64, nstate, nstate)
    b_c = zeros(Float64, nstate, nsat)

    for sat in 1:nsat
        state_rows = ks_state_slice(sat)
        base = first(state_rows) - 1
        p_sat, q_sat, h_sat, _ = ks_state_components(u, sat)
        r_sat = ks_position(p_sat)
        rho = ks_radius_scalar(p_sat)
        s_mat = ks_s_matrix(p_sat)
        v_sat = ks_velocity(p_sat, q_sat)
        cart_rows = (3 * (sat - 1) + 1):(3 * sat)

        a_j2 = ks_j2_acceleration(r_sat; mu = mu, j2 = j2, earth_radius_km = earth_radius_km)
        a_total = a_j2 + laser_accels[:, sat]
        jr = ks_j2_acceleration_jacobian(r_sat; mu = mu, j2 = j2, earth_radius_km = earth_radius_km) .+
            jr_laser_cart[cart_rows, cart_rows]
        jv = zeros(Float64, 3, 3)

        a_c[state_rows, state_rows] .= ks_satellite_continuous_jacobian_block(
            p_sat,
            q_sat,
            h_sat,
            a_total;
            jr = jr,
            jv = jv,
        )

        for source_sat in 1:nsat
            source_sat == sat && continue

            source_rows = ks_state_slice(source_sat)
            source_p, source_q, _, _ = ks_state_components(u, source_sat)
            source_cart_rows = (3 * (source_sat - 1) + 1):(3 * source_sat)
            jr_cross = jr_laser_cart[cart_rows, source_cart_rows]
            iszero(norm(jr_cross)) && continue

            a_c[state_rows, source_rows] .+= ks_satellite_force_coupling_block(
                p_sat,
                q_sat,
                source_p,
                source_q;
                jr = jr_cross,
            )
        end

        ju_sat = ju_cart[(3 * (sat - 1) + 1):(3 * sat), :]
        ju4 = zeros(Float64, 4, size(ju_sat, 2))
        ju4[1:3, :] .= ju_sat
        b_c[(base + 5):(base + 8), :] .= 0.5 .* rho .* (s_mat * ju4)
        b_c[base + 9, :] .= vec(-2.0 .* rho .* (transpose(v_sat) * ju_sat))
    end

    return (
        A = a_c,
        B = b_c,
        pair_data = pair_data,
        cartesian_input_jacobian = ju_cart,
        control_amplitudes = controls,
    )
end

function ks_state_components_view(u::AbstractVector, sat_id::Int)
    base = 10 * (sat_id - 1)
    p = view(u, (base + 1):(base + 4))
    q = view(u, (base + 5):(base + 8))
    return p, q, u[base + 9], u[base + 10]
end

function ks_lambda_matrix_generic(p::AbstractVector)
    T = promote_type(eltype(p), Float64)
    return T[
        p[1] -p[2] -p[3] p[4]
        p[2] p[1] -p[4] -p[3]
        p[3] p[4] p[1] p[2]
    ]
end

function ks_l_matrix_generic(p::AbstractVector)
    T = promote_type(eltype(p), Float64)
    return T[
        p[1] -p[2] -p[3] p[4]
        p[2] p[1] -p[4] -p[3]
        p[3] p[4] p[1] p[2]
        p[4] -p[3] p[2] -p[1]
    ]
end

ks_s_matrix_generic(p::AbstractVector) = transpose(ks_l_matrix_generic(p))
ks_radius_scalar_generic(p::AbstractVector) = dot(p, p)

function ks_position_generic(p::AbstractVector)
    return ks_lambda_matrix_generic(p) * p
end

function ks_velocity_generic(p::AbstractVector, q::AbstractVector)
    rho = ks_radius_scalar_generic(p)
    rho > zero(rho) || throw(ArgumentError("KS radius must be positive."))
    return (2.0 / rho) .* (ks_lambda_matrix_generic(p) * q)
end

function ks_augmented_acceleration_generic(accel::AbstractVector)
    T = promote_type(eltype(accel), Float64)
    accel4 = zeros(T, 4)
    accel4[1:3] .= accel
    return accel4
end

function ks_j2_acceleration_generic(
    r::AbstractVector;
    mu::Real = MU_EARTH_KM3_S2,
    j2::Real = J2_EARTH,
    earth_radius_km::Real = R_EARTH_KM,
)
    T = promote_type(eltype(r), typeof(float(mu)), typeof(float(j2)), typeof(float(earth_radius_km)))
    x = r[1]
    y = r[2]
    z = r[3]
    r2 = x * x + y * y + z * z
    rmag = sqrt(r2)
    inv_r5 = inv(r2 * r2 * rmag)
    z2 = z * z
    z_ratio = T(5.0) * z2 / r2
    scale = T(1.5) * T(float(j2)) * T(float(mu)) * T(float(earth_radius_km))^2 * inv_r5
    return scale .* T[
        x * (z_ratio - one(T)),
        y * (z_ratio - one(T)),
        z * (z_ratio - T(3.0)),
    ]
end

function ks_satellite_rhs_sundman_time_generic!(
    du_sat::AbstractVector,
    p::AbstractVector,
    q::AbstractVector,
    h,
    accel::AbstractVector,
)
    T = promote_type(eltype(du_sat), eltype(p), eltype(q), eltype(accel), typeof(h))
    rho = ks_radius_scalar_generic(p)
    rho > zero(rho) || throw(ArgumentError("KS radius must be positive."))

    v = ks_velocity_generic(p, q)
    s_mat = ks_s_matrix_generic(p)
    accel4 = ks_augmented_acceleration_generic(accel)

    du_sat[1:4] .= q
    du_sat[5:8] .= (-T(0.25) * h) .* p .+ T(0.5) .* rho .* (s_mat * accel4)
    du_sat[9] = -T(2.0) * rho * dot(v, accel)
    du_sat[10] = rho
    return nothing
end

function ks_constellation_rhs_smooth(
    u::AbstractVector;
    pair_data::AbstractVector,
    mu::Real = MU_EARTH_KM3_S2,
    j2::Real = J2_EARTH,
    earth_radius_km::Real = R_EARTH_KM,
    cr::Real = 2.0,
    laser_power_w::Real = 1.0e4,
    satellite_mass_kg::Real = 300.0,
    control_amplitudes::AbstractVector,
)
    nsat = ks_satellite_count(u)
    T = promote_type(eltype(u), eltype(control_amplitudes), Float64)
    du = zeros(T, length(u))
    positions = zeros(T, 3, nsat)

    for sat in 1:nsat
        p_sat, _, _, _ = ks_state_components_view(u, sat)
        positions[:, sat] .= ks_position_generic(p_sat)
    end

    laser_accels = ks_total_laser_accelerations_from_positions(
        positions,
        pair_data;
        control_amplitudes = control_amplitudes,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
    )

    for sat in 1:nsat
        p_sat, q_sat, h_sat, _ = ks_state_components_view(u, sat)
        total_accel = ks_j2_acceleration_generic(
            view(positions, :, sat);
            mu = mu,
            j2 = j2,
            earth_radius_km = earth_radius_km,
        ) + laser_accels[:, sat]

        ks_satellite_rhs_sundman_time_generic!(
            view(du, ks_state_slice(sat)),
            p_sat,
            q_sat,
            h_sat,
            total_accel,
        )
    end

    return du
end

function ks_forwarddiff_continuous_jacobians(
    u::AbstractVector{<:Real};
    mu::Real = MU_EARTH_KM3_S2,
    j2::Real = J2_EARTH,
    earth_radius_km::Real = R_EARTH_KM,
    link_max_range_km::Real = 200.0,
    atmosphere_top_km::Real = 100.0,
    cr::Real = 2.0,
    laser_power_w::Real = 1.0e4,
    satellite_mass_kg::Real = 300.0,
    control_amplitudes::Union{Nothing, AbstractVector{<:Real}} = nothing,
    pair_list::AbstractVector{<:Tuple{Int, Int}} = ks_pair_indices(ks_satellite_count(u)),
)
    pair_data = ks_pair_geometry_data(
        u;
        pair_list = pair_list,
        link_max_range_km = link_max_range_km,
        atmosphere_top_km = atmosphere_top_km,
    )
    controls = isnothing(control_amplitudes) ? ones(Float64, ks_satellite_count(u)) : Float64.(control_amplitudes)

    rhs_state = x -> ks_constellation_rhs_smooth(
        x;
        pair_data = pair_data,
        mu = mu,
        j2 = j2,
        earth_radius_km = earth_radius_km,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
        control_amplitudes = controls,
    )
    rhs_control = uc -> ks_constellation_rhs_smooth(
        u;
        pair_data = pair_data,
        mu = mu,
        j2 = j2,
        earth_radius_km = earth_radius_km,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
        control_amplitudes = uc,
    )

    return (
        A = Matrix{Float64}(ForwardDiff.jacobian(rhs_state, Float64.(u))),
        B = Matrix{Float64}(ForwardDiff.jacobian(rhs_control, controls)),
        pair_data = pair_data,
        control_amplitudes = controls,
    )
end

function ks_radial_unit_and_jacobian(r::AbstractVector{<:Real})
    r_vec = Float64.(r)
    rmag = norm(r_vec)
    rmag > 0.0 || throw(ArgumentError("Position magnitude must be positive."))

    rhat = r_vec ./ rmag
    jr = (Matrix{Float64}(I, 3, 3) .- rhat * transpose(rhat)) ./ rmag
    return rhat, jr
end

function ks_radial_control_continuous_jacobians(
    u::AbstractVector{<:Real};
    mu::Real = MU_EARTH_KM3_S2,
    j2::Real = J2_EARTH,
    earth_radius_km::Real = R_EARTH_KM,
    radial_control_accelerations::Union{Nothing, AbstractVector{<:Real}} = nothing,
)
    nsat = ks_satellite_count(u)
    nstate = 10 * nsat
    controls = isnothing(radial_control_accelerations) ? zeros(Float64, nsat) : Float64.(radial_control_accelerations)
    length(controls) == nsat ||
        throw(ArgumentError("radial_control_accelerations must have one entry per satellite."))

    a_c = zeros(Float64, nstate, nstate)
    b_c = zeros(Float64, nstate, nsat)
    radial_directions = zeros(Float64, 3, nsat)

    for sat in 1:nsat
        state_rows = ks_state_slice(sat)
        base = first(state_rows) - 1
        p_sat, q_sat, h_sat, _ = ks_state_components(u, sat)
        r_sat = ks_position(p_sat)
        rho = ks_radius_scalar(p_sat)
        s_mat = ks_s_matrix(p_sat)
        v_sat = ks_velocity(p_sat, q_sat)
        rhat, jr_radial_direction = ks_radial_unit_and_jacobian(r_sat)
        radial_directions[:, sat] .= rhat

        a_j2 = ks_j2_acceleration(r_sat; mu = mu, j2 = j2, earth_radius_km = earth_radius_km)
        a_total = a_j2 .+ controls[sat] .* rhat
        jr = ks_j2_acceleration_jacobian(r_sat; mu = mu, j2 = j2, earth_radius_km = earth_radius_km) .+
            controls[sat] .* jr_radial_direction

        a_c[state_rows, state_rows] .= ks_satellite_continuous_jacobian_block(
            p_sat,
            q_sat,
            h_sat,
            a_total;
            jr = jr,
        )

        radial4 = zeros(Float64, 4)
        radial4[1:3] .= rhat
        b_c[(base + 5):(base + 8), sat] .= 0.5 .* rho .* (s_mat * radial4)
        b_c[base + 9, sat] = -2.0 * rho * dot(v_sat, rhat)
    end

    return (
        A = a_c,
        B = b_c,
        radial_directions = radial_directions,
        radial_control_accelerations = controls,
    )
end

function ks_radial_unit_generic(r::AbstractVector)
    rmag = norm(r)
    rmag > zero(rmag) || throw(ArgumentError("Position magnitude must be positive."))
    return r ./ rmag
end

function ks_constellation_rhs_radial_control_smooth(
    u::AbstractVector;
    mu::Real = MU_EARTH_KM3_S2,
    j2::Real = J2_EARTH,
    earth_radius_km::Real = R_EARTH_KM,
    radial_control_accelerations::AbstractVector,
)
    nsat = ks_satellite_count(u)
    length(radial_control_accelerations) == nsat ||
        throw(ArgumentError("radial_control_accelerations length must match the number of satellites."))

    T = promote_type(eltype(u), eltype(radial_control_accelerations), Float64)
    du = zeros(T, length(u))

    for sat in 1:nsat
        p_sat, q_sat, h_sat, _ = ks_state_components_view(u, sat)
        r_sat = ks_position_generic(p_sat)
        rhat = ks_radial_unit_generic(r_sat)
        total_accel = ks_j2_acceleration_generic(
            r_sat;
            mu = mu,
            j2 = j2,
            earth_radius_km = earth_radius_km,
        ) .+ radial_control_accelerations[sat] .* rhat

        ks_satellite_rhs_sundman_time_generic!(
            view(du, ks_state_slice(sat)),
            p_sat,
            q_sat,
            h_sat,
            total_accel,
        )
    end

    return du
end

function ks_forwarddiff_radial_control_continuous_jacobians(
    u::AbstractVector{<:Real};
    mu::Real = MU_EARTH_KM3_S2,
    j2::Real = J2_EARTH,
    earth_radius_km::Real = R_EARTH_KM,
    radial_control_accelerations::Union{Nothing, AbstractVector{<:Real}} = nothing,
)
    controls = isnothing(radial_control_accelerations) ? zeros(Float64, ks_satellite_count(u)) : Float64.(radial_control_accelerations)

    rhs_state = x -> ks_constellation_rhs_radial_control_smooth(
        x;
        mu = mu,
        j2 = j2,
        earth_radius_km = earth_radius_km,
        radial_control_accelerations = controls,
    )
    rhs_control = uc -> ks_constellation_rhs_radial_control_smooth(
        u;
        mu = mu,
        j2 = j2,
        earth_radius_km = earth_radius_km,
        radial_control_accelerations = uc,
    )

    return (
        A = Matrix{Float64}(ForwardDiff.jacobian(rhs_state, Float64.(u))),
        B = Matrix{Float64}(ForwardDiff.jacobian(rhs_control, controls)),
        radial_control_accelerations = controls,
    )
end

function ks_finite_difference_radial_control_continuous_jacobians(
    u::AbstractVector{<:Real};
    mu::Real = MU_EARTH_KM3_S2,
    j2::Real = J2_EARTH,
    earth_radius_km::Real = R_EARTH_KM,
    radial_control_accelerations::Union{Nothing, AbstractVector{<:Real}} = nothing,
    state_step::Real = 1e-6,
    control_step::Real = 1e-8,
)
    x0 = Float64.(u)
    controls = isnothing(radial_control_accelerations) ? zeros(Float64, ks_satellite_count(u)) : Float64.(radial_control_accelerations)
    nstate = length(x0)
    ncontrols = length(controls)

    rhs_state = x -> ks_constellation_rhs_radial_control_smooth(
        x;
        mu = mu,
        j2 = j2,
        earth_radius_km = earth_radius_km,
        radial_control_accelerations = controls,
    )
    rhs_control = uc -> ks_constellation_rhs_radial_control_smooth(
        x0;
        mu = mu,
        j2 = j2,
        earth_radius_km = earth_radius_km,
        radial_control_accelerations = uc,
    )

    a_c = zeros(Float64, nstate, nstate)
    b_c = zeros(Float64, nstate, ncontrols)
    dx = float(state_step)
    du = float(control_step)

    for col in 1:nstate
        xp = copy(x0)
        xm = copy(x0)
        xp[col] += dx
        xm[col] -= dx
        a_c[:, col] .= (rhs_state(xp) .- rhs_state(xm)) ./ (2.0 * dx)
    end

    for col in 1:ncontrols
        up = copy(controls)
        um = copy(controls)
        up[col] += du
        um[col] -= du
        b_c[:, col] .= (rhs_control(up) .- rhs_control(um)) ./ (2.0 * du)
    end

    return (
        A = a_c,
        B = b_c,
        radial_control_accelerations = controls,
        state_step = dx,
        control_step = du,
    )
end

function ks_cartesian_state_from_ks(u::AbstractVector{<:Real})
    nsat = ks_satellite_count(u)
    x = zeros(Float64, 6 * nsat)

    for sat in 1:nsat
        p_sat, q_sat, _, _ = ks_state_components(u, sat)
        base = 6 * (sat - 1)
        x[(base + 1):(base + 3)] .= ks_position(p_sat)
        x[(base + 4):(base + 6)] .= ks_velocity(p_sat, q_sat)
    end

    return x
end

function ks_cartesian_state_jacobian_from_ks(u::AbstractVector{<:Real})
    nsat = ks_satellite_count(u)
    h_mat = zeros(Float64, 6 * nsat, 10 * nsat)

    for sat in 1:nsat
        p_sat, q_sat, _, _ = ks_state_components(u, sat)
        cart_base = 6 * (sat - 1)
        ks_base = 10 * (sat - 1)

        dv_dp, dv_dq = ks_velocity_jacobians(p_sat, q_sat)
        h_mat[(cart_base + 1):(cart_base + 3), (ks_base + 1):(ks_base + 4)] .= 2.0 .* ks_lambda_matrix(p_sat)
        h_mat[(cart_base + 4):(cart_base + 6), (ks_base + 1):(ks_base + 4)] .= dv_dp
        h_mat[(cart_base + 4):(cart_base + 6), (ks_base + 5):(ks_base + 8)] .= dv_dq
    end

    return h_mat
end

function ks_pq_state_indices(nsat::Int)
    indices = Int[]
    for sat in 1:nsat
        base = 10 * (sat - 1)
        append!(indices, (base + 1):(base + 8))
    end
    return indices
end

function ks_cartesian_state_jacobian_from_ks_pq(u::AbstractVector{<:Real})
    return ks_cartesian_state_jacobian_from_ks(u)[:, ks_pq_state_indices(ks_satellite_count(u))]
end

function ks_pq_time_scaling(u::AbstractVector{<:Real})
    nsat = ks_satellite_count(u)
    scaling = zeros(Float64, 8 * nsat, 8 * nsat)
    for sat in 1:nsat
        p_sat, _, _, _ = ks_state_components(u, sat)
        rho = ks_radius_scalar(p_sat)
        rows = (8 * (sat - 1) + 1):(8 * sat)
        scaling[rows, rows] .= (1.0 / rho) .* Matrix{Float64}(I, 8, 8)
    end
    return scaling
end

function ks_physical_time_rhs_from_sundman(
    u::AbstractVector{<:Real},
    sundman_rhs::AbstractVector{<:Real},
)
    nsat = ks_satellite_count(u)
    time_rhs = zeros(Float64, length(u))
    for sat in 1:nsat
        p_sat, _, _, _ = ks_state_components(u, sat)
        rho = ks_radius_scalar(p_sat)
        rows = ks_state_slice(sat)
        time_rhs[rows] .= Float64.(sundman_rhs[rows]) ./ rho
    end
    return time_rhs
end

function ks_rtn_cartesian_transform_matrix_from_ks(u::AbstractVector{<:Real})
    nsat = ks_satellite_count(u)
    transform = zeros(Float64, 6 * nsat, 6 * nsat)

    for sat in 1:nsat
        r_sat, v_sat, _, _ = ks_cartesian_components(u, sat)
        m = ks_rtn_basis(r_sat, v_sat)
        rows = (6 * (sat - 1) + 1):(6 * sat)
        transform[rows, rows] .= [
            transpose(m) zeros(Float64, 3, 3)
            zeros(Float64, 3, 3) transpose(m)
        ]
    end

    return transform
end

function ks_pq_reduction_basis(
    h_pq::AbstractMatrix{<:Real};
    rtol::Real = 1e-9,
    atol::Real = 0.0,
)
    fact = svd(Matrix{Float64}(h_pq))
    isempty(fact.S) && return zeros(Float64, size(h_pq, 2), 0)
    max_sv = maximum(fact.S)
    threshold = max(float(atol), float(rtol) * max_sv)
    rank = count(>(threshold), fact.S)
    return fact.V[:, 1:rank]
end

function ks_rtn_t_matrix_from_ks_pq(
    u::AbstractVector{<:Real},
    n_basis::AbstractMatrix{<:Real},
)
    r_transform = ks_rtn_cartesian_transform_matrix_from_ks(u)
    h_pq = ks_cartesian_state_jacobian_from_ks_pq(u)
    return r_transform * h_pq * Matrix{Float64}(n_basis)
end

function ks_cartesian_pair_geometry_data(
    x::AbstractVector{<:Real};
    pair_list::AbstractVector{<:Tuple{Int, Int}},
    link_max_range_km::Real,
    atmosphere_top_km::Real,
)
    length(x) % 6 == 0 || throw(ArgumentError("Cartesian state length must be divisible by 6."))
    nsat = length(x) ÷ 6
    positions = Vector{Vector{Float64}}(undef, nsat)
    for sat in 1:nsat
        base = 6 * (sat - 1)
        positions[sat] = Float64.(x[(base + 1):(base + 3)])
    end

    occlusion_radius_km = R_EARTH_KM + float(atmosphere_top_km)
    data = Vector{NamedTuple}(undef, length(pair_list))
    for (idx, (i, j)) in enumerate(pair_list)
        dvec = positions[j] - positions[i]
        distance = norm(dvec)
        khat = distance > 0.0 ? dvec ./ distance : zeros(Float64, 3)
        active = distance > 0.0 &&
            distance <= float(link_max_range_km) &&
            ks_has_direct_los(positions[i], positions[j]; occlusion_radius_km = occlusion_radius_km)

        data[idx] = (
            pair = (i, j),
            distance_km = distance,
            khat = khat,
            zeta = active ? 1.0 : 0.0,
        )
    end

    return data
end

function ks_twobody_acceleration_generic(r::AbstractVector; mu::Real = MU_EARTH_KM3_S2)
    T = promote_type(eltype(r), typeof(float(mu)), Float64)
    r2 = dot(r, r)
    r2 > zero(r2) || throw(ArgumentError("Position magnitude must be positive."))
    return .-T(float(mu)) .* r ./ (r2 * sqrt(r2))
end

function ks_cartesian_constellation_rhs_sundman_smooth(
    x::AbstractVector;
    pair_data::AbstractVector,
    mu::Real = MU_EARTH_KM3_S2,
    j2::Real = J2_EARTH,
    earth_radius_km::Real = R_EARTH_KM,
    cr::Real = 2.0,
    laser_power_w::Real = 1.0e4,
    satellite_mass_kg::Real = 300.0,
    control_amplitudes::AbstractVector,
)
    length(x) % 6 == 0 || throw(ArgumentError("Cartesian state length must be divisible by 6."))
    nsat = length(x) ÷ 6
    length(control_amplitudes) == nsat ||
        throw(ArgumentError("control_amplitudes must have one entry per satellite."))

    T = promote_type(eltype(x), eltype(control_amplitudes), Float64)
    dx = zeros(T, length(x))
    positions = zeros(T, 3, nsat)

    for sat in 1:nsat
        base = 6 * (sat - 1)
        positions[:, sat] .= x[(base + 1):(base + 3)]
    end

    laser_accels = ks_total_laser_accelerations_from_positions(
        positions,
        pair_data;
        control_amplitudes = control_amplitudes,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
    )

    for sat in 1:nsat
        base = 6 * (sat - 1)
        r_sat = view(positions, :, sat)
        v_sat = view(x, (base + 4):(base + 6))
        rho = norm(r_sat)
        total_accel = ks_twobody_acceleration_generic(r_sat; mu = mu) .+
            ks_j2_acceleration_generic(
                r_sat;
                mu = mu,
                j2 = j2,
                earth_radius_km = earth_radius_km,
            ) .+
            laser_accels[:, sat]

        dx[(base + 1):(base + 3)] .= rho .* v_sat
        dx[(base + 4):(base + 6)] .= rho .* total_accel
    end

    return dx
end

function ks_forwarddiff_cartesian_continuous_jacobians(
    x::AbstractVector{<:Real};
    mu::Real = MU_EARTH_KM3_S2,
    j2::Real = J2_EARTH,
    earth_radius_km::Real = R_EARTH_KM,
    link_max_range_km::Real = 200.0,
    atmosphere_top_km::Real = 100.0,
    cr::Real = 2.0,
    laser_power_w::Real = 1.0e4,
    satellite_mass_kg::Real = 300.0,
    control_amplitudes::Union{Nothing, AbstractVector{<:Real}} = nothing,
    pair_list::AbstractVector{<:Tuple{Int, Int}} = ks_pair_indices(length(x) ÷ 6),
)
    length(x) % 6 == 0 || throw(ArgumentError("Cartesian state length must be divisible by 6."))
    nsat = length(x) ÷ 6
    controls = isnothing(control_amplitudes) ? ones(Float64, nsat) : Float64.(control_amplitudes)
    pair_data = ks_cartesian_pair_geometry_data(
        x;
        pair_list = pair_list,
        link_max_range_km = link_max_range_km,
        atmosphere_top_km = atmosphere_top_km,
    )

    rhs_state = y -> ks_cartesian_constellation_rhs_sundman_smooth(
        y;
        pair_data = pair_data,
        mu = mu,
        j2 = j2,
        earth_radius_km = earth_radius_km,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
        control_amplitudes = controls,
    )
    rhs_control = uc -> ks_cartesian_constellation_rhs_sundman_smooth(
        x;
        pair_data = pair_data,
        mu = mu,
        j2 = j2,
        earth_radius_km = earth_radius_km,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
        control_amplitudes = uc,
    )

    return (
        A = Matrix{Float64}(ForwardDiff.jacobian(rhs_state, Float64.(x))),
        B = Matrix{Float64}(ForwardDiff.jacobian(rhs_control, controls)),
        pair_data = pair_data,
        control_amplitudes = controls,
    )
end

function ks_controllability_matrix(
    a::AbstractMatrix{<:Real},
    b::AbstractMatrix{<:Real};
    horizon::Union{Nothing, Int} = nothing,
)
    n = size(a, 1)
    size(a, 2) == n || throw(ArgumentError("A must be square."))
    size(b, 1) == n || throw(ArgumentError("B row count must match A."))

    steps = isnothing(horizon) ? n : horizon
    cmat = zeros(Float64, n, steps * size(b, 2))
    block = Matrix{Float64}(b)
    a_mat = Matrix{Float64}(a)

    for k in 1:steps
        cols = ((k - 1) * size(b, 2) + 1):(k * size(b, 2))
        cmat[:, cols] .= block
        block = a_mat * block
    end

    return cmat
end

function ks_controllability_basis(
    a::AbstractMatrix{<:Real},
    b::AbstractMatrix{<:Real};
    horizon::Union{Nothing, Int} = nothing,
    rtol::Real = 1e-9,
    atol::Real = 0.0,
)
    n = size(a, 1)
    steps = isnothing(horizon) ? n : horizon
    a_mat = Matrix{Float64}(a)
    q = ks_orthonormal_column_basis(b; rtol = rtol, atol = atol)

    for _ in 2:steps
        size(q, 2) == 0 && return q
        q = ks_orthonormal_column_basis([q a_mat * q]; rtol = rtol, atol = atol)
    end

    return q
end

function ks_numerical_rank(m::AbstractMatrix{<:Real}; rtol::Real = 1e-9, atol::Real = 0.0)
    isempty(m) && return 0
    singular_values = svdvals(Matrix{Float64}(m))
    isempty(singular_values) && return 0
    max_sv = maximum(singular_values)
    max_sv == 0.0 && return 0
    threshold = max(float(atol), float(rtol) * max_sv)
    return count(>(threshold), singular_values)
end

function ks_orthonormal_column_basis(m::AbstractMatrix{<:Real}; rtol::Real = 1e-9, atol::Real = 0.0)
    m_float = Matrix{Float64}(m)
    isempty(m_float) && return zeros(Float64, size(m_float, 1), 0)
    fact = svd(m_float)
    isempty(fact.S) && return zeros(Float64, size(m_float, 1), 0)
    max_sv = maximum(fact.S)
    max_sv == 0.0 && return zeros(Float64, size(m_float, 1), 0)
    threshold = max(float(atol), float(rtol) * max_sv)
    rank = count(>(threshold), fact.S)
    return fact.U[:, 1:rank]
end

function ks_subspace_residual(
    source::AbstractMatrix{<:Real},
    target::AbstractMatrix{<:Real};
    rtol::Real = 1e-9,
    atol::Real = 0.0,
)
    q_source = ks_orthonormal_column_basis(source; rtol = rtol, atol = atol)
    q_target = ks_orthonormal_column_basis(target; rtol = rtol, atol = atol)
    size(q_source, 2) == 0 && return 0.0
    size(q_target, 2) == 0 && return opnorm(q_source)
    return opnorm(q_source .- q_target * (transpose(q_target) * q_source))
end

"""Full physical-time KS tangent and analytical Cartesian/RTN push-forward.

As, Bs, fs are Sundman-time quantities, including energy and clock rows.
The horizontal lift G preserves the velocity gauge and Kepler energy constraint.
"""
function ks_physical_linearization(u, As, Bs, fs; mu = MU_EARTH_KM3_S2)
    nsat = ks_satellite_count(u)
    At, Bt = Matrix{Float64}(As), Matrix{Float64}(Bs)
    ft = ks_physical_time_rhs_from_sundman(u, fs)
    H = ks_cartesian_state_jacobian_from_ks(u)
    Hd = zeros(6nsat, 10nsat)
    G = zeros(10nsat, 6nsat)
    R, Rd = zeros(6nsat, 6nsat), zeros(6nsat, 6nsat)
    cart_rhs = H * ft
    for sat in 1:nsat
        rows = ks_state_slice(sat)
        cols = (6(sat-1)+1):(6sat)
        p, q, _, _ = ks_state_components(u, sat)
        r, v, _, _ = ks_cartesian_components(u, sat)
        rho = dot(p, p)
        At[rows, :] ./= rho
        At[rows, rows[1:4]] .-= fs[rows] * (2p)' / rho^2
        Bt[rows, :] ./= rho
        pd, qd = ft[rows[1:4]], ft[rows[5:8]]
        rhod = 2dot(p, pd)
        acceleration = cart_rhs[cols[4:6]]
        hd, g = view(Hd, cols, rows), view(G, rows, cols)
        hd[1:3, 1:4] .= 2ks_lambda_matrix(pd)
        hd[4:6, 1:4] .= 2ks_lambda_matrix(qd)/rho -
            2rhod*ks_lambda_matrix(q)/rho^2 -
            2(acceleration*p' + v*pd')/rho + 2rhod*v*p'/rho^2
        hd[4:6, 5:8] .= 2ks_lambda_matrix(pd)/rho - 2rhod*ks_lambda_matrix(p)/rho^2
        g[1:4, 1:3] .= ks_lambda_matrix(p)'/(2rho)
        Qp = 0.5hcat([KS_DSP[j]*vcat(v, 0.0) for j in 1:4]...)
        g[5:8, :] .= Qp*g[1:4, :]
        g[5:8, 4:6] .+= 0.5ks_lambda_matrix(p)'
        g[9, :] .= vcat(-2mu*r/rho^3, -2v)
        frame = ks_rtn_basis(r, v)
        framed = ks_rtn_basis_derivative(r, v, acceleration)
        for block in (cols[1:3], cols[4:6])
            R[block, block] .= frame'
            Rd[block, block] .= framed'
        end
    end
    Ac, Bc = (Hd + H*At)*G, H*Bt
    return (A_time=At, B_time=Bt, rhs_time=ft, H=H, Hdot=Hd, G=G,
        R=R, Rdot=Rd, A_cartesian=Ac, B_cartesian=Bc,
        A_RTN=Rd*R' + R*Ac*R', B_RTN=R*Bc)
end

function ks_cartesian_controllability_comparison(
    u::AbstractVector{<:Real};
    mu::Real = MU_EARTH_KM3_S2,
    j2::Real = J2_EARTH,
    earth_radius_km::Real = R_EARTH_KM,
    link_max_range_km::Real = 200.0,
    atmosphere_top_km::Real = 100.0,
    cr::Real = 2.0,
    laser_power_w::Real = 1.0e4,
    satellite_mass_kg::Real = 300.0,
    control_amplitudes::Union{Nothing, AbstractVector{<:Real}} = nothing,
    pair_list::AbstractVector{<:Tuple{Int, Int}} = ks_pair_indices(ks_satellite_count(u)),
    rtol::Real = 1e-9,
    atol::Real = 0.0,
)
    x = ks_cartesian_state_from_ks(u)
    controls = isnothing(control_amplitudes) ? ones(Float64, ks_satellite_count(u)) : Float64.(control_amplitudes)

    ks_lin = ks_continuous_jacobians(
        u;
        mu = mu,
        j2 = j2,
        earth_radius_km = earth_radius_km,
        link_max_range_km = link_max_range_km,
        atmosphere_top_km = atmosphere_top_km,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
        control_amplitudes = controls,
        pair_list = pair_list,
    )
    cart_lin = ks_forwarddiff_cartesian_continuous_jacobians(
        x;
        mu = mu,
        j2 = j2,
        earth_radius_km = earth_radius_km,
        link_max_range_km = link_max_range_km,
        atmosphere_top_km = atmosphere_top_km,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
        control_amplitudes = controls,
        pair_list = pair_list,
    )

    fs = ks_constellation_rhs_smooth(u; pair_data=ks_lin.pair_data, mu=mu,
        j2=j2, earth_radius_km=earth_radius_km, cr=cr,
        laser_power_w=laser_power_w, satellite_mass_kg=satellite_mass_kg,
        control_amplitudes=controls)
    physical = ks_physical_linearization(u, ks_lin.A, ks_lin.B, fs; mu=mu)
    # Independent Cartesian reference: convert its Sundman Jacobian by the quotient rule.
    cart_fs = ks_cartesian_constellation_rhs_sundman_smooth(x;
        pair_data=cart_lin.pair_data, mu=mu, j2=j2, earth_radius_km=earth_radius_km,
        cr=cr, laser_power_w=laser_power_w, satellite_mass_kg=satellite_mass_kg,
        control_amplitudes=controls)
    for sat in 1:ks_satellite_count(u)
        rows = (6(sat-1)+1):(6sat)
        r = x[rows[1:3]]
        rho = norm(r)
        cart_lin.A[rows, :] ./= rho
        cart_lin.A[rows, rows[1:3]] .-= cart_fs[rows]*r'/rho^3
        cart_lin.B[rows, :] ./= rho
    end
    ks_sundman = ks_lin
    h_mat = physical.H
    # Freeze the physical Cartesian tangent, then embed it in KS coordinates.
    # H*K(A_time,B_time) is NOT a valid LTI comparison when H changes with time.
    ks_lin = merge(ks_lin, (A=physical.G*physical.A_cartesian*h_mat,
        B=physical.G*physical.B_cartesian,
        representation=:frozen_physical_tangent_lift))
    cart_cmat = ks_controllability_matrix(cart_lin.A, cart_lin.B)
    ks_cmat = ks_controllability_matrix(ks_lin.A, ks_lin.B)
    projected_ks_cmat = h_mat * ks_cmat
    ks_basis = ks_controllability_basis(ks_lin.A, ks_lin.B; rtol = rtol, atol = atol)
    cart_basis = ks_controllability_basis(cart_lin.A, cart_lin.B; rtol = rtol, atol = atol)
    projected_ks_basis = ks_orthonormal_column_basis(h_mat * ks_basis; rtol = rtol, atol = atol)

    cart_rank = size(cart_basis, 2)
    projected_ks_rank = size(projected_ks_basis, 2)
    ks_rank = size(ks_basis, 2)

    return (
        ks = ks_lin,
        ks_sundman = ks_sundman,
        physical = physical,
        cartesian = cart_lin,
        cartesian_state = x,
        cartesian_projection = h_mat,
        ks_controllability = ks_cmat,
        cartesian_controllability = cart_cmat,
        projected_ks_controllability = projected_ks_cmat,
        ks_controllability_basis = ks_basis,
        cartesian_controllability_basis = cart_basis,
        projected_ks_cartesian_basis = projected_ks_basis,
        ks_rank = ks_rank,
        cartesian_rank = cart_rank,
        projected_ks_cartesian_rank = projected_ks_rank,
        ks_to_cartesian_residual = ks_subspace_residual(projected_ks_basis, cart_basis; rtol = rtol, atol = atol),
        cartesian_to_ks_residual = ks_subspace_residual(cart_basis, projected_ks_basis; rtol = rtol, atol = atol),
        input_projection_error = maximum(abs.(h_mat * ks_lin.B .- cart_lin.B)),
        rank_match = projected_ks_rank == cart_rank,
        pair_data = ks_lin.pair_data,
    )
end

function ks_rtn_basis(
    r::AbstractVector{<:Real},
    v::AbstractVector{<:Real},
)
    r_vec = Float64.(r)
    v_vec = Float64.(v)
    r_norm = norm(r_vec)
    h_vec = cross(r_vec, v_vec)
    h_norm = norm(h_vec)
    r_norm > 0.0 || throw(ArgumentError("Position norm must be positive."))
    h_norm > 0.0 || throw(ArgumentError("Angular momentum norm must be positive."))

    r_hat = r_vec ./ r_norm
    n_hat = h_vec ./ h_norm
    t_hat = cross(n_hat, r_hat)
    return hcat(r_hat, t_hat, n_hat)
end

function ks_rtn_basis_derivative(
    r::AbstractVector{<:Real},
    v::AbstractVector{<:Real},
    accel::AbstractVector{<:Real},
)
    r_vec = Float64.(r)
    v_vec = Float64.(v)
    accel_vec = Float64.(accel)
    m = ks_rtn_basis(r_vec, v_vec)
    r_hat = m[:, 1]
    t_hat = m[:, 2]
    n_hat = m[:, 3]

    r_norm = norm(r_vec)
    h_vec = cross(r_vec, v_vec)
    h_norm = norm(h_vec)

    r_hat_dot = (Matrix{Float64}(I, 3, 3) .- r_hat * transpose(r_hat)) * v_vec ./ r_norm
    h_dot = cross(r_vec, accel_vec)
    n_hat_dot = (Matrix{Float64}(I, 3, 3) .- n_hat * transpose(n_hat)) * h_dot ./ h_norm
    t_hat_dot = cross(n_hat_dot, r_hat) + cross(n_hat, r_hat_dot)

    return hcat(r_hat_dot, t_hat_dot, n_hat_dot)
end

function ks_relative_position_selector_rtn1(
    m1::AbstractMatrix{<:Real},
    m2::AbstractMatrix{<:Real},
)
    size(m1) == (3, 3) || throw(ArgumentError("m1 must be 3x3."))
    size(m2) == (3, 3) || throw(ArgumentError("m2 must be 3x3."))
    i3 = Matrix{Float64}(I, 3, 3)
    z3 = zeros(Float64, 3, 3)
    q12 = transpose(Matrix{Float64}(m1)) * Matrix{Float64}(m2)
    return hcat(-i3, z3, q12, z3)
end

function ks_relative_state_selector_rtn1(
    m1::AbstractMatrix{<:Real},
    m2::AbstractMatrix{<:Real},
)
    p_rho = ks_relative_position_selector_rtn1(m1, m2)
    z3 = zeros(Float64, 3, 3)
    i3 = Matrix{Float64}(I, 3, 3)
    q12 = transpose(Matrix{Float64}(m1)) * Matrix{Float64}(m2)
    p_rho_dot = hcat(z3, -i3, z3, q12)
    return vcat(p_rho, p_rho_dot)
end

function ks_barycenter_relative_transform_rtn1(
    m1::AbstractMatrix{<:Real},
    m2::AbstractMatrix{<:Real},
)
    size(m1) == (3, 3) || throw(ArgumentError("m1 must be 3x3."))
    size(m2) == (3, 3) || throw(ArgumentError("m2 must be 3x3."))
    i3 = Matrix{Float64}(I, 3, 3)
    z3 = zeros(Float64, 3, 3)
    q12 = transpose(Matrix{Float64}(m1)) * Matrix{Float64}(m2)
    p_cm_r = hcat(0.5 .* i3, z3, 0.5 .* q12, z3)
    p_cm_v = hcat(z3, 0.5 .* i3, z3, 0.5 .* q12)
    p_rho = hcat(-i3, z3, q12, z3)
    p_rho_dot = hcat(z3, -i3, z3, q12)
    return vcat(p_cm_r, p_cm_v, p_rho, p_rho_dot)
end

function ks_relative_position_controllability_eigensystem(
    controllability_rtn::AbstractMatrix{<:Real},
    m1::AbstractMatrix{<:Real},
    m2::AbstractMatrix{<:Real},
)
    size(controllability_rtn, 1) == 12 ||
        throw(ArgumentError("Two-satellite RTN controllability matrix must have 12 rows."))
    selector = ks_relative_position_selector_rtn1(m1, m2)
    controllability_relative = selector * Matrix{Float64}(controllability_rtn)
    gramian_relative = controllability_relative * transpose(controllability_relative)
    eig = eigen(Symmetric(0.5 .* (gramian_relative .+ transpose(gramian_relative))))
    dominant_index = argmax(eig.values)
    dominant_vector = eig.vectors[:, dominant_index]
    dominant_vector ./= norm(dominant_vector)
    return (
        eigenvalues = eig.values,
        eigenvectors = eig.vectors,
        dominant_value = eig.values[dominant_index],
        dominant_vector = dominant_vector,
        gramian = gramian_relative,
        controllability = controllability_relative,
        selector = selector,
    )
end

function ks_cartesian_rtn_transform(
    x::AbstractVector{<:Real};
    pair_data::AbstractVector,
    mu::Real = MU_EARTH_KM3_S2,
    j2::Real = J2_EARTH,
    earth_radius_km::Real = R_EARTH_KM,
    cr::Real = 2.0,
    laser_power_w::Real = 1.0e4,
    satellite_mass_kg::Real = 300.0,
    control_amplitudes::AbstractVector{<:Real},
)
    length(x) % 6 == 0 || throw(ArgumentError("Cartesian state length must be divisible by 6."))
    nsat = length(x) ÷ 6
    length(control_amplitudes) == nsat ||
        throw(ArgumentError("control_amplitudes must have one entry per satellite."))

    positions = zeros(Float64, 3, nsat)
    velocities = zeros(Float64, 3, nsat)
    for sat in 1:nsat
        base = 6 * (sat - 1)
        positions[:, sat] .= Float64.(x[(base + 1):(base + 3)])
        velocities[:, sat] .= Float64.(x[(base + 4):(base + 6)])
    end

    laser_accels = ks_total_laser_accelerations_from_positions(
        positions,
        pair_data;
        control_amplitudes = control_amplitudes,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
    )

    transform = zeros(Float64, 6 * nsat, 6 * nsat)
    transform_dot = zeros(Float64, 6 * nsat, 6 * nsat)
    sundman_to_time = zeros(Float64, 6 * nsat, 6 * nsat)
    accelerations = zeros(Float64, 3, nsat)

    for sat in 1:nsat
        base = 6 * (sat - 1)
        r_sat = positions[:, sat]
        v_sat = velocities[:, sat]
        accel_sat = ks_twobody_acceleration_generic(r_sat; mu = mu) .+
            ks_j2_acceleration(
                r_sat;
                mu = mu,
                j2 = j2,
                earth_radius_km = earth_radius_km,
            ) .+
            laser_accels[:, sat]
        accelerations[:, sat] .= accel_sat

        m = ks_rtn_basis(r_sat, v_sat)
        m_dot = ks_rtn_basis_derivative(r_sat, v_sat, accel_sat)
        rows = (base + 1):(base + 6)
        transform[rows, rows] .= [
            transpose(m) zeros(Float64, 3, 3)
            zeros(Float64, 3, 3) transpose(m)
        ]
        transform_dot[rows, rows] .= [
            transpose(m_dot) zeros(Float64, 3, 3)
            zeros(Float64, 3, 3) transpose(m_dot)
        ]
        sundman_to_time[rows, rows] .= (1.0 / norm(r_sat)) .* Matrix{Float64}(I, 6, 6)
    end

    return (
        transform = transform,
        transform_dot = transform_dot,
        inverse_transform = transpose(transform),
        sundman_to_time = sundman_to_time,
        accelerations = accelerations,
    )
end

function ks_laser_rtn_linearization(
    u::AbstractVector{<:Real};
    mu::Real = MU_EARTH_KM3_S2,
    j2::Real = J2_EARTH,
    earth_radius_km::Real = R_EARTH_KM,
    link_max_range_km::Real = 200.0,
    atmosphere_top_km::Real = 100.0,
    cr::Real = 2.0,
    laser_power_w::Real = 1.0e12,
    satellite_mass_kg::Real = 300.0,
    control_amplitudes::Union{Nothing, AbstractVector{<:Real}} = nothing,
    rtol::Real = 1e-9,
)
    nsat = ks_satellite_count(u)
    controls = isnothing(control_amplitudes) ? ones(Float64, nsat) : Float64.(control_amplitudes)
    comparison = ks_cartesian_controllability_comparison(
        u;
        mu = mu,
        j2 = j2,
        earth_radius_km = earth_radius_km,
        link_max_range_km = link_max_range_km,
        atmosphere_top_km = atmosphere_top_km,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
        control_amplitudes = controls,
        rtol = rtol,
    )

    rtn = ks_cartesian_rtn_transform(
        comparison.cartesian_state;
        pair_data = comparison.pair_data,
        mu = mu,
        j2 = j2,
        earth_radius_km = earth_radius_km,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
        control_amplitudes = controls,
    )

    physical = comparison.physical
    a_rtn, b_rtn = physical.A_RTN, physical.B_RTN
    eig = eigen(a_rtn)
    controllability = ks_controllability_matrix(a_rtn, b_rtn)
    controllability_basis = ks_controllability_basis(a_rtn, b_rtn; rtol = rtol)

    return (
        A_RTN = a_rtn,
        B_RTN = b_rtn,
        eigenvalues = eig.values,
        eigenvectors = eig.vectors,
        controllability = controllability,
        controllability_basis = controllability_basis,
        controllability_rank = size(controllability_basis, 2),
        comparison = comparison,
        rtn_transform = rtn,
        physical = physical,
        H = physical.H,
        G = physical.G,
        Hdot = physical.Hdot,
    )
end

function ks_midpoint_discretization(
    a_c::AbstractMatrix{<:Real},
    b_c::AbstractMatrix{<:Real},
    dt::Real,
)
    n = size(a_c, 1)
    i_n = Matrix{Float64}(I, n, n)
    dtf = float(dt)
    lhs = i_n .- 0.5 .* dtf .* Matrix{Float64}(a_c)
    rhs = i_n .+ 0.5 .* dtf .* Matrix{Float64}(a_c)
    a_d = lhs \ rhs
    b_d = lhs \ (dtf .* Matrix{Float64}(b_c))
    return a_d, b_d
end

function compute_ks_linearization_sequences(
    sol;
    mu::Real = MU_EARTH_KM3_S2,
    j2::Real = J2_EARTH,
    earth_radius_km::Real = R_EARTH_KM,
    link_max_range_km::Real = 200.0,
    atmosphere_top_km::Real = 100.0,
    cr::Real = 2.0,
    laser_power_w::Real = 1.0e4,
    satellite_mass_kg::Real = 300.0,
    control_amplitudes::Union{Nothing, AbstractVector{<:Real}} = nothing,
)
    nsamples = length(sol.t)
    nsamples >= 2 || throw(ArgumentError("Solution must contain at least two time samples."))

    a_c_seq = Matrix{Float64}[]
    b_c_seq = Matrix{Float64}[]
    a_d_seq = Matrix{Float64}[]
    b_d_seq = Matrix{Float64}[]
    pair_data_seq = Vector{NamedTuple}[]

    for k in 1:(nsamples - 1)
        dt = float(sol.t[k + 1] - sol.t[k])
        x_mid = 0.5 .* (Float64.(sol.u[k]) .+ Float64.(sol.u[k + 1]))

        cont = ks_continuous_jacobians(
            x_mid;
            mu = mu,
            j2 = j2,
            earth_radius_km = earth_radius_km,
            link_max_range_km = link_max_range_km,
            atmosphere_top_km = atmosphere_top_km,
            cr = cr,
            laser_power_w = laser_power_w,
            satellite_mass_kg = satellite_mass_kg,
            control_amplitudes = control_amplitudes,
        )
        a_d, b_d = ks_midpoint_discretization(cont.A, cont.B, dt)

        push!(a_c_seq, cont.A)
        push!(b_c_seq, cont.B)
        push!(a_d_seq, a_d)
        push!(b_d_seq, b_d)
        push!(pair_data_seq, cont.pair_data)
    end

    return (
        A_c = a_c_seq,
        B_c = b_c_seq,
        A_d = a_d_seq,
        B_d = b_d_seq,
        pair_data = pair_data_seq,
    )
end

function rollout_ks_discrete_linearization(
    reference_sol,
    nonlinear_sol;
    relinearize::Bool,
    linearization = nothing,
    mu::Real = MU_EARTH_KM3_S2,
    j2::Real = J2_EARTH,
    earth_radius_km::Real = R_EARTH_KM,
    link_max_range_km::Real = 200.0,
    atmosphere_top_km::Real = 100.0,
    cr::Real = 2.0,
    laser_power_w::Real = 1.0e4,
    satellite_mass_kg::Real = 300.0,
    control_amplitudes::Union{Nothing, AbstractVector{<:Real}} = nothing,
)
    length(reference_sol.t) == length(nonlinear_sol.t) ||
        throw(ArgumentError("Reference and nonlinear solutions must share the same sample grid."))

    lin = isnothing(linearization) ? compute_ks_linearization_sequences(
        reference_sol;
        mu = mu,
        j2 = j2,
        earth_radius_km = earth_radius_km,
        link_max_range_km = link_max_range_km,
        atmosphere_top_km = atmosphere_top_km,
        cr = cr,
        laser_power_w = laser_power_w,
        satellite_mass_kg = satellite_mass_kg,
        control_amplitudes = control_amplitudes,
    ) : linearization

    nsamples = length(reference_sol.t)
    predictions = Vector{Vector{Float64}}(undef, nsamples)
    delta = Float64.(nonlinear_sol.u[1]) .- Float64.(reference_sol.u[1])
    predictions[1] = Float64.(reference_sol.u[1]) .+ delta

    for k in 1:(nsamples - 1)
        a_d = relinearize ? lin.A_d[k] : lin.A_d[1]
        delta = a_d * delta
        predictions[k + 1] = Float64.(reference_sol.u[k + 1]) .+ delta
    end

    return predictions
end

function ks_position_error_timeseries(
    truth_states::AbstractVector{<:AbstractVector{<:Real}},
    approx_states::AbstractVector{<:AbstractVector{<:Real}},
)
    length(truth_states) == length(approx_states) ||
        throw(ArgumentError("truth_states and approx_states must have the same length."))

    nsamples = length(truth_states)
    nsat = ks_satellite_count(truth_states[1])
    errors = zeros(Float64, nsat, nsamples)

    for k in 1:nsamples
        for sat in 1:nsat
            r_truth, _, _, _ = ks_cartesian_components(truth_states[k], sat)
            r_approx, _, _, _ = ks_cartesian_components(approx_states[k], sat)
            errors[sat, k] = norm(r_truth - r_approx)
        end
    end

    return errors
end

function ks_rms_position_error(errors::AbstractMatrix{<:Real})
    nsat, nsamples = size(errors)
    rms = zeros(Float64, nsamples)
    for k in 1:nsamples
        rms[k] = sqrt(sum(abs2, errors[:, k]) / nsat)
    end
    return rms
end
