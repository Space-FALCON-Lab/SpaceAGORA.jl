function initialize_coplanar_ks_constellation(
    N::Int;
    altitude_km::Real = 550.0,
    inclination_deg::Real = 53.0,
    eccentricity::Real = 0.0,
    argp_deg::Real = 0.0,
    raan_deg::Real = 0.0,
    nu0_deg::Real = 0.0,
    anomaly_spacing_deg::Real = 0.5,
    anomaly_offsets_deg::Union{Nothing, AbstractVector{<:Real}} = nothing,
    altitude_offsets_km::Union{Nothing, AbstractVector{<:Real}} = nothing,
    mu::Real = MU_EARTH_KM3_S2,
)
    N > 0 || throw(ArgumentError("N must be positive."))
    altitude_km > 0 || throw(ArgumentError("altitude_km must be positive."))

    offsets = isnothing(anomaly_offsets_deg) ? zeros(Float64, N) : Float64.(anomaly_offsets_deg)
    length(offsets) == N || throw(ArgumentError("anomaly_offsets_deg must have length N."))
    altitude_offsets = isnothing(altitude_offsets_km) ? zeros(Float64, N) : Float64.(altitude_offsets_km)
    length(altitude_offsets) == N || throw(ArgumentError("altitude_offsets_km must have length N."))

    a_km = R_EARTH_KM + float(altitude_km)
    i_rad = deg2rad(float(inclination_deg))
    argp_rad = deg2rad(float(argp_deg))
    raan_rad = deg2rad(float(raan_deg))
    nu0_rad = deg2rad(float(nu0_deg))
    spacing_rad = deg2rad(float(anomaly_spacing_deg))

    u0 = zeros(Float64, 10 * N)
    anomalies_deg = zeros(Float64, N)

    for sat in 1:N
        sat_a_km = a_km + altitude_offsets[sat]
        sat_a_km > R_EARTH_KM || throw(ArgumentError("Satellite semi-major axis must exceed Earth radius."))
        nu_rad = nu0_rad + (sat - 1) * spacing_rad + deg2rad(offsets[sat])
        r_eci, v_eci = coe_to_rv(sat_a_km, eccentricity, i_rad, raan_rad, argp_rad, nu_rad; mu = mu)
        x_ks = cartesian_to_ks_state(r_eci, v_eci; mu = mu, t0 = 0.0)
        u0[ks_state_slice(sat)] .= x_ks
        anomalies_deg[sat] = rad2deg(nu_rad)
    end

    return (
        u0 = u0,
        N = N,
        a_km = a_km,
        altitude_km = float(altitude_km),
        inclination_deg = float(inclination_deg),
        eccentricity = float(eccentricity),
        argp_deg = float(argp_deg),
        raan_deg = float(raan_deg),
        nu_deg = anomalies_deg,
        anomaly_spacing_deg = float(anomaly_spacing_deg),
        altitude_offsets_km = altitude_offsets,
    )
end
