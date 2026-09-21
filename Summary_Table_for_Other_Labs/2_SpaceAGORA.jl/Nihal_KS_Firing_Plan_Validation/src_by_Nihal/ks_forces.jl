function ks_has_direct_los(
    r1::AbstractVector{<:Real},
    r2::AbstractVector{<:Real};
    occlusion_radius_km::Real,
)
    sx = float(r2[1]) - float(r1[1])
    sy = float(r2[2]) - float(r1[2])
    sz = float(r2[3]) - float(r1[3])

    seg_len2 = sx * sx + sy * sy + sz * sz
    seg_len2 == 0.0 && return false

    t_star = -(float(r1[1]) * sx + float(r1[2]) * sy + float(r1[3]) * sz) / seg_len2
    t_clamped = clamp(t_star, 0.0, 1.0)

    cx = float(r1[1]) + t_clamped * sx
    cy = float(r1[2]) + t_clamped * sy
    cz = float(r1[3]) + t_clamped * sz
    return cx * cx + cy * cy + cz * cz > float(occlusion_radius_km)^2
end

function ks_laser_acceleration_gain(
    cr::Real,
    laser_power_w::Real,
    satellite_mass_kg::Real,
)
    satellite_mass_kg > 0 || throw(ArgumentError("satellite_mass_kg must be positive."))
    c_m_s = 299_792_458.0
    return float(cr) * float(laser_power_w) / (float(satellite_mass_kg) * c_m_s) / 1000.0
end

function ks_j2_acceleration(
    r::AbstractVector{<:Real};
    mu::Real = MU_EARTH_KM3_S2,
    j2::Real = J2_EARTH,
    earth_radius_km::Real = R_EARTH_KM,
)
    x = float(r[1])
    y = float(r[2])
    z = float(r[3])
    r2 = x * x + y * y + z * z
    rmag = sqrt(r2)
    inv_r5 = inv(r2 * r2 * rmag)
    z2 = z * z
    z_ratio = 5.0 * z2 / r2
    scale = 1.5 * float(j2) * float(mu) * float(earth_radius_km)^2 * inv_r5
    return scale .* [
        x * (z_ratio - 1.0),
        y * (z_ratio - 1.0),
        z * (z_ratio - 3.0),
    ]
end

function ks_j2_acceleration_jacobian(
    r::AbstractVector{<:Real};
    mu::Real = MU_EARTH_KM3_S2,
    j2::Real = J2_EARTH,
    earth_radius_km::Real = R_EARTH_KM,
)
    x = float(r[1])
    y = float(r[2])
    z = float(r[3])
    r2 = x * x + y * y + z * z
    r2 > 0.0 || throw(ArgumentError("Position magnitude must be non-zero for J2 Jacobian evaluation."))

    c = 1.5 * float(j2) * float(mu) * float(earth_radius_km)^2
    iszero(c) && return zeros(Float64, 3, 3)

    r4 = r2 * r2
    z2 = z * z
    factor = c / (r2^4 * sqrt(r2))

    m = Float64[
        -r4 + 5.0 * r2 * (x * x + z2) - 35.0 * x * x * z2 5.0 * x * y * (r2 - 7.0 * z2) 5.0 * x * z * (3.0 * r2 - 7.0 * z2)
        5.0 * x * y * (r2 - 7.0 * z2) -r4 + 5.0 * r2 * (y * y + z2) - 35.0 * y * y * z2 5.0 * y * z * (3.0 * r2 - 7.0 * z2)
        5.0 * x * z * (3.0 * r2 - 7.0 * z2) 5.0 * y * z * (3.0 * r2 - 7.0 * z2) -3.0 * r4 + 30.0 * r2 * z2 - 35.0 * z2 * z2
    ]
    return factor .* m
end

function ks_pair_geometry_data(
    u::AbstractVector{<:Real};
    pair_list::AbstractVector{<:Tuple{Int, Int}} = ks_pair_indices(ks_satellite_count(u)),
    link_max_range_km::Real = 200.0,
    atmosphere_top_km::Real = 100.0,
)
    nsat = ks_satellite_count(u)
    positions = Vector{Vector{Float64}}(undef, nsat)
    for sat in 1:nsat
        positions[sat], _, _, _ = ks_cartesian_components(u, sat)
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

function ks_laser_cartesian_input_jacobian(
    pair_data::AbstractVector,
    nsat::Int;
    cr::Real,
    laser_power_w::Real,
    satellite_mass_kg::Real,
)
    ncontrols = nsat
    ju = zeros(Float64, 3 * nsat, ncontrols)
    gain = ks_laser_acceleration_gain(cr, laser_power_w, satellite_mass_kg)

    for data in pair_data
        i, j = data.pair
        direction = (data.zeta * gain) .* data.khat
        ju[(3 * (i - 1) + 1):(3 * i), i] .-= direction
        ju[(3 * (j - 1) + 1):(3 * j), j] .+= direction
    end

    return ju
end

function ks_total_laser_accelerations(
    pair_data::AbstractVector,
    nsat::Int;
    control_amplitudes::AbstractVector{<:Real},
    cr::Real,
    laser_power_w::Real,
    satellite_mass_kg::Real,
)
    length(control_amplitudes) == nsat ||
        throw(ArgumentError("control_amplitudes length must match the number of satellites."))

    gain = ks_laser_acceleration_gain(cr, laser_power_w, satellite_mass_kg)
    accelerations = zeros(Float64, 3, nsat)

    for data in pair_data
        i, j = data.pair
        base_direction = (data.zeta * gain) .* data.khat
        accelerations[:, i] .-= float(control_amplitudes[i]) .* base_direction
        accelerations[:, j] .+= float(control_amplitudes[j]) .* base_direction
    end

    return accelerations
end

function ks_laser_cartesian_state_jacobian(
    pair_data::AbstractVector,
    nsat::Int;
    control_amplitudes::AbstractVector{<:Real},
    cr::Real,
    laser_power_w::Real,
    satellite_mass_kg::Real,
)
    length(control_amplitudes) == nsat ||
        throw(ArgumentError("control_amplitudes length must match the number of satellites."))

    jr = zeros(Float64, 3 * nsat, 3 * nsat)
    gain = ks_laser_acceleration_gain(cr, laser_power_w, satellite_mass_kg)
    i3 = Matrix{Float64}(I, 3, 3)

    for data in pair_data
        data.zeta == 0.0 && continue

        distance = float(data.distance_km)
        distance > 0.0 || continue

        i, j = data.pair
        khat = Float64.(data.khat)
        scale = data.zeta * gain / distance
        pair_jac = scale .* (i3 .- khat * transpose(khat))

        rows_i = (3 * (i - 1) + 1):(3 * i)
        rows_j = (3 * (j - 1) + 1):(3 * j)

        amplitude_i = float(control_amplitudes[i])
        pair_jac_i = amplitude_i .* pair_jac
        jr[rows_i, rows_i] .+= pair_jac_i
        jr[rows_i, rows_j] .-= pair_jac_i

        amplitude_j = float(control_amplitudes[j])
        pair_jac_j = amplitude_j .* pair_jac
        jr[rows_j, rows_i] .-= pair_jac_j
        jr[rows_j, rows_j] .+= pair_jac_j
    end

    return jr
end

function ks_total_laser_accelerations_from_positions(
    positions::AbstractMatrix,
    pair_data::AbstractVector;
    control_amplitudes::AbstractVector,
    cr::Real,
    laser_power_w::Real,
    satellite_mass_kg::Real,
)
    nsat = size(positions, 2)
    length(control_amplitudes) == nsat ||
        throw(ArgumentError("control_amplitudes length must match the number of satellites."))

    T = promote_type(eltype(positions), eltype(control_amplitudes), Float64)
    gain = convert(T, ks_laser_acceleration_gain(cr, laser_power_w, satellite_mass_kg))
    accelerations = zeros(T, 3, nsat)

    for data in pair_data
        zeta = convert(T, data.zeta)
        iszero(zeta) && continue

        i, j = data.pair
        dvec = positions[:, j] .- positions[:, i]
        distance = norm(dvec)
        iszero(distance) && continue

        base_direction = (zeta * gain / distance) .* dvec
        accelerations[:, i] .-= convert(T, control_amplitudes[i]) .* base_direction
        accelerations[:, j] .+= convert(T, control_amplitudes[j]) .* base_direction
    end

    return accelerations
end
