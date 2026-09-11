"""
This module provides functions to convert between position/velocity (r, v) and classical orbital elements (coe).
"""

"""
    position & velocity (r, v) -> classical orbital elements (coe) (robust to circular/equatorial)

    Inputs:
        r: position vector (3-vector)
        v: velocity vector (3-vector)
        μ: gravitational parameter
        tol: tolerance for numerical stability (default 1e-10)

    Returns:
        a: semi-major axis
        e: eccentricity
        i: inclination
        Ω: right ascension of ascending node (RAAN)
        ω: argument of perigee
        ν: true anomaly
        u: argument of latitude (for circular orbits)
        p: semi-latus rectum
"""
function rv2coe(r::SVector{3,Float64}, v::SVector{3,Float64}, μ::Float64; tol=1e-10)

    rnorm = norm(r); v2 = dot(v,v); h = cross(r,v); hnorm = norm(h) #
    ĥ = h / (hnorm+tol) # ω/ν/u sign must be measured about the orbit normal, not Ẑ
    i = acos(clamp(h[3]/(hnorm+tol), -1.0, 1.0)) # inclination; acos(-1, 1)
    n = cross(Ẑ, h); nnorm = norm(n)
    evec = (cross(v,h)/μ) - r/(rnorm+tol); e = norm(evec) # eccentricity vector
    ε = 0.5*v2 - μ/rnorm # specific orbital energy
    a = -μ/(2ε) # semi-major axis
    p = hnorm^2/μ # semi-latus rectum

    Ω = nnorm > tol ? atan(n[2], n[1]) : 0.0 # RAAN

    if e > 1e-8 && nnorm > tol # non-circular, non-equatorial
        ω = atan(dot(ĥ, cross(n, evec)) / (nnorm*e + tol),
                dot(n, evec) / (nnorm*e + tol)) # argument of perigee
        ν = atan(dot(ĥ, cross(evec, r)) / (e*rnorm + tol),
                dot(evec, r) / (e*rnorm + tol)) # true anomaly
        u = NaN # u is the argument of latitude, not defined for non-circular orbits
    else
        ω = 0.0
        u = nnorm > tol ? atan(dot(ĥ, cross(n, r))/(nnorm*rnorm + tol),
                                dot(n, r)/(nnorm*rnorm + tol)) :
                        atan(r[2], r[1]) # argument of latitude
        ν = u  # periapsis defined at ascending node (ω=0) for circular orbits
    end
    return (a=a, e=e, i=i, Ω=Ω, ω=ω, ν=ν, u=u, p=p)
end

"""
    rv2coe_2pi — identical to rv2coe but all angles are returned in [0, 2π) instead of (-π, π].
    Intended for difference / time-series analysis where the ±180° branch-cut jump in rv2coe
    produces artificial ±360° spikes in Δangle plots.
    The only remaining discontinuity is at 0/2π (once per orbit), which is easy to unwrap if needed.
"""
function rv2coe_2pi(r::SVector{3,Float64}, v::SVector{3,Float64}, μ::Float64; tol=1e-10)

    rnorm = norm(r); v2 = dot(v,v); h = cross(r,v); hnorm = norm(h)
    ĥ  = h / (hnorm+tol) # ω/ν/u sign must be measured about the orbit normal, not Ẑ
    i  = acos(clamp(h[3]/(hnorm+tol), -1.0, 1.0))
    n  = cross(Ẑ, h); nnorm = norm(n)
    evec = (cross(v,h)/μ) - r/(rnorm+tol); e = norm(evec)
    ε  = 0.5*v2 - μ/rnorm
    a  = -μ/(2ε)
    p  = hnorm^2/μ

    Ω = nnorm > tol ? mod(atan(n[2], n[1]), 2π) : 0.0

    if e > 1e-8 && nnorm > tol  # non-circular, non-equatorial
        ω = mod(atan(dot(ĥ, cross(n, evec)) / (nnorm*e + tol),
                     dot(n, evec)            / (nnorm*e + tol)), 2π)
        ν = mod(atan(dot(ĥ, cross(evec, r)) / (e*rnorm + tol),
                     dot(evec, r)            / (e*rnorm + tol)), 2π)
        u = NaN
    else                        # circular (and/or equatorial)
        ω = 0.0
        u = nnorm > tol ? mod(atan(dot(ĥ, cross(n, r)) / (nnorm*rnorm + tol),
                                   dot(n, r)            / (nnorm*rnorm + tol)), 2π) :
                          mod(atan(r[2], r[1]), 2π)
        ν = u
    end
    return (a=a, e=e, i=i, Ω=Ω, ω=ω, ν=ν, u=u, p=p)
end

"""
    classical orbital elements (coe) -> position & velocity (r, v)

    Inputs:
        a: semi-major axis
        e: eccentricity
        i: inclination
        Ω: right ascension of ascending node (RAAN)
        ω: argument of perigee
        ν: true anomaly
        μ: gravitational parameter

    Returns:
        r: position vector (3-vector)
        v: velocity vector (3-vector)
"""
function coe2rv(a::Float64, e::Float64, i::Float64, Ω::Float64, ω::Float64, ν::Float64, μ::Float64) # μ is the gravitational parameter

    p = a*(1 - e^2)
    r_pf = @SVector [p*cos(ν)/(1 + e*cos(ν)), p*sin(ν)/(1 + e*cos(ν)), 0.0]
    v_pf = @SVector [-sqrt(μ/p)*sin(ν), sqrt(μ/p)*(e + cos(ν)), 0.0]
    cΩ=cos(Ω); sΩ=sin(Ω); ci=cos(i); si=sin(i); cω=cos(ω); sω=sin(ω)
    # perifocal->ECI rotation (3-1-3)
    R = @SMatrix [ cΩ*cω - sΩ*sω*ci   -cΩ*sω - sΩ*cω*ci    sΩ*si;
                sΩ*cω + cΩ*sω*ci   -sΩ*sω + cΩ*cω*ci   -cΩ*si;
                sω*si               cω*si               ci    ]
    r = R * r_pf
    v = R * v_pf
    return r, v
end

"""
    Build state from OE (degrees or radians) # OE means orbital elements

        Inputs:
            a_m: semi-major axis (meters)
            e: eccentricity
            i_deg: inclination (degrees, default 0.0)
            Ω_deg: right ascension of ascending node (degrees, default 0.0)
            ω_deg: argument of perigee (degrees, default 0.0)
            ν_deg: true anomaly (degrees, default 0.0)
            μ: gravitational parameter (default MU)

        Returns:
            state tuple: (x, y, z, vx, vy, vz)
"""
function state_from_OE(a_m; e=0.0, i_deg=0.0, Ω_deg=0.0, ω_deg=0.0, ν_deg=0.0, μ=MU)

    i = deg2rad(i_deg); Ω = deg2rad(Ω_deg); ω = deg2rad(ω_deg); ν = deg2rad(ν_deg) # convert format
    r, v = coe2rv(a_m, e, i, Ω, ω, ν, μ) # core function is coe2rr
    return (r[1],r[2],r[3], v[1],v[2],v[3]) #convert to tuple format
end

"""
    This function builds a 6N state vector u by concatenating states generated from a list of
    orbital element named-tuples. Each element must have fields:
    a_m, e, i_deg, Ω_deg, ω_deg, ν_deg

        Inputs:
            oe_list: list of named tuples with orbital elements
            μ: gravitational parameter (default MU)

        Returns:
            u: state vector (6N)
"""
function build_u_from_oe(oe_list; μ=MU)
    u = Vector{Float64}()
    for oe in oe_list
        s = state_from_OE(oe.a_m; e=oe.e, i_deg=oe.i_deg, Ω_deg=oe.Ω_deg, ω_deg=oe.ω_deg, ν_deg=oe.ν_deg, μ=μ)
        append!(u, s)
    end
    return u
end