"""
    coe_to_rv(a_km, e, i_rad, raan_rad, argp_rad, nu_rad; mu=MU_EARTH_KM3_S2)

Convert classical orbital elements to an Earth-centered inertial Cartesian state.
All distances are in km and time is in seconds.
"""
function coe_to_rv(
    a_km::Real,
    e::Real,
    i_rad::Real,
    raan_rad::Real,
    argp_rad::Real,
    nu_rad::Real;
    mu::Real = MU_EARTH_KM3_S2,
)
    p = a_km * (1 - e^2)
    r_mag = p / (1 + e * cos(nu_rad))

    r_pf = [r_mag * cos(nu_rad), r_mag * sin(nu_rad), 0.0]
    v_pf = sqrt(mu / p) * [-sin(nu_rad), e + cos(nu_rad), 0.0]

    cO = cos(raan_rad)
    sO = sin(raan_rad)
    ci = cos(i_rad)
    si = sin(i_rad)
    cw = cos(argp_rad)
    sw = sin(argp_rad)

    q_p_to_i = [
        cO * cw - sO * sw * ci -cO * sw - sO * cw * ci sO * si
        sO * cw + cO * sw * ci -sO * sw + cO * cw * ci -cO * si
        sw * si cw * si ci
    ]

    r_eci = q_p_to_i * r_pf
    v_eci = q_p_to_i * v_pf
    return r_eci, v_eci
end

"""
    initialize_walker_delta(N, P, F; kwargs...)

Create initial states for a Walker-delta constellation `i : N / P / F` around Earth.
Returns a named tuple with a packed state vector `u0` (`[x y z vx vy vz]` for each satellite).
"""
function initialize_walker_delta(
    N::Int,
    P::Int,
    F::Int;
    altitude_km::Real = 550.0,
    inclination_deg::Real = 53.0,
    eccentricity::Real = 0.0,
    argp_deg::Real = 0.0,
    raan0_deg::Real = 0.0,
    nu0_deg::Real = 0.0,
    mu::Real = MU_EARTH_KM3_S2,
)
    N > 0 || throw(ArgumentError("N must be positive."))
    P > 0 || throw(ArgumentError("P (number of planes) must be positive."))
    N % P == 0 || throw(ArgumentError("N must be divisible by P for a uniform Walker delta constellation."))
    0 <= F < P || throw(ArgumentError("F must satisfy 0 <= F < P."))
    altitude_km > 0 || throw(ArgumentError("altitude_km must be positive."))
    0 <= eccentricity < 1 || throw(ArgumentError("eccentricity must satisfy 0 <= e < 1."))

    sats_per_plane = div(N, P)
    a_km = R_EARTH_KM + float(altitude_km)

    i_rad = deg2rad(float(inclination_deg))
    argp_rad = deg2rad(float(argp_deg))
    raan0_rad = deg2rad(float(raan0_deg))
    nu0_rad = deg2rad(float(nu0_deg))

    u0 = Vector{Float64}(undef, 6 * N)

    sat_index = 0
    for plane in 0:(P - 1)
        raan_rad = raan0_rad + 2pi * plane / P

        for slot in 0:(sats_per_plane - 1)
            # Walker-delta phasing between planes, expressed in true anomaly for circular orbits.
            nu_rad = nu0_rad + 2pi * slot / sats_per_plane + 2pi * F * plane / N
            r_eci, v_eci = coe_to_rv(a_km, eccentricity, i_rad, raan_rad, argp_rad, nu_rad; mu = mu)

            base = 6 * sat_index
            u0[base + 1] = r_eci[1]
            u0[base + 2] = r_eci[2]
            u0[base + 3] = r_eci[3]
            u0[base + 4] = v_eci[1]
            u0[base + 5] = v_eci[2]
            u0[base + 6] = v_eci[3]

            sat_index += 1
        end
    end

    return (u0 = u0, N = N, P = P, F = F, sats_per_plane = sats_per_plane, a_km = a_km)
end

function recommended_walker_plane_count(N::Int)
    N > 0 || throw(ArgumentError("N must be positive."))

    # Pick the largest divisor at or below sqrt(N) for a balanced plane layout.
    guess = floor(Int, sqrt(float(N)))
    for d in guess:-1:1
        if N % d == 0
            return d
        end
    end

    return 1
end

"""
    define_walker_delta_constellation(N; kwargs...)

High-level Walker-delta initializer where `N` is the only required argument.
If `P` is omitted, a balanced plane count is selected from divisors of `N`.
If `F` is omitted, a default phasing of 1 (or 0 when `P == 1`) is used.
"""
function define_walker_delta_constellation(
    N::Int;
    P::Union{Nothing, Int} = nothing,
    F::Union{Nothing, Int} = nothing,
    altitude_km::Real = 550.0,
    inclination_deg::Real = 53.0,
    eccentricity::Real = 0.0,
    argp_deg::Real = 0.0,
    raan0_deg::Real = 0.0,
    nu0_deg::Real = 0.0,
    mu::Real = MU_EARTH_KM3_S2,
)
    P_val = isnothing(P) ? recommended_walker_plane_count(N) : P
    F_default = P_val == 1 ? 0 : 1
    F_val = isnothing(F) ? F_default : mod(F, P_val)

    return initialize_walker_delta(
        N,
        P_val,
        F_val;
        altitude_km = altitude_km,
        inclination_deg = inclination_deg,
        eccentricity = eccentricity,
        argp_deg = argp_deg,
        raan0_deg = raan0_deg,
        nu0_deg = nu0_deg,
        mu = mu,
    )
end

"""
    define_flower_constellation(N; kwargs...)

High-level Flower-like initializer where `N` is the only required argument.
This builds a phased set of similar orbits whose angular elements are coupled
to create a rosette/flower-style spatial pattern.
"""
function define_flower_constellation(
    N::Int;
    petals::Union{Nothing, Int} = nothing,
    altitude_km::Real = 1200.0,
    inclination_deg::Real = 63.4,
    eccentricity::Real = 0.15,
    argp0_deg::Real = 270.0,
    raan0_deg::Real = 0.0,
    nu0_deg::Real = 0.0,
    mu::Real = MU_EARTH_KM3_S2,
)
    N > 0 || throw(ArgumentError("N must be positive."))
    altitude_km > 0 || throw(ArgumentError("altitude_km must be positive."))
    0 <= eccentricity < 1 || throw(ArgumentError("eccentricity must satisfy 0 <= e < 1."))

    petals_val = if isnothing(petals)
        min(N, max(3, round(Int, sqrt(float(N)))))
    else
        petals
    end
    petals_val > 0 || throw(ArgumentError("petals must be positive."))

    a_km = R_EARTH_KM + float(altitude_km)
    i_rad = deg2rad(float(inclination_deg))
    argp0_rad = deg2rad(float(argp0_deg))
    raan0_rad = deg2rad(float(raan0_deg))
    nu0_rad = deg2rad(float(nu0_deg))

    u0 = Vector{Float64}(undef, 6 * N)

    for sat in 0:(N - 1)
        phase = 2pi * sat / N

        # Coupled angular phasing produces a flower-like rosette structure.
        raan_rad = raan0_rad + phase
        argp_rad = argp0_rad + petals_val * phase
        nu_rad = nu0_rad - petals_val * phase

        r_eci, v_eci = coe_to_rv(a_km, eccentricity, i_rad, raan_rad, argp_rad, nu_rad; mu = mu)

        base = 6 * sat
        u0[base + 1] = r_eci[1]
        u0[base + 2] = r_eci[2]
        u0[base + 3] = r_eci[3]
        u0[base + 4] = v_eci[1]
        u0[base + 5] = v_eci[2]
        u0[base + 6] = v_eci[3]
    end

    return (u0 = u0, N = N, petals = petals_val, a_km = a_km)
end
