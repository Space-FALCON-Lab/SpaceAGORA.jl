# ------------------------------------------------------------------
# z-hat unit vector (inertial pole), used by rv2coe / rv2coe_2pi
# ------------------------------------------------------------------
if !@isdefined(Ẑ); const Ẑ = SVector(0.0, 0.0, 1.0); end

"""
    rv2coe(r, v, μ) → NamedTuple(a, e, i, Ω, ω, ν, u, p)

Convert ECI position/velocity to classical orbital elements.
Robust to circular and equatorial orbits (returns u for circular, ν=NaN).
"""
function rv2coe(r::SVector{3,Float64}, v::SVector{3,Float64}, μ::Float64; tol=1e-10)
    rnorm = norm(r); v2 = dot(v,v); h = cross(r,v); hnorm = norm(h)
    i  = acos(clamp(h[3]/(hnorm+tol), -1.0, 1.0))
    n  = cross(Ẑ, h); nnorm = norm(n)
    evec = (cross(v,h)/μ) - r/(rnorm+tol); e = norm(evec)
    ε  = 0.5*v2 - μ/rnorm
    a  = -μ/(2ε)
    p  = hnorm^2/μ
    Ω  = nnorm > tol ? atan(n[2], n[1]) : 0.0
    if e > 1e-8 && nnorm > tol
        ω = atan(dot(Ẑ, cross(n,evec))/(nnorm*e+tol), dot(n,evec)/(nnorm*e+tol))
        ν = atan(dot(Ẑ, cross(evec,r))/(e*rnorm+tol), dot(evec,r)/(e*rnorm+tol))
        u = NaN
    else
        ω = 0.0
        u = nnorm > tol ? atan(dot(Ẑ,cross(n,r))/(nnorm*rnorm+tol), dot(n,r)/(nnorm*rnorm+tol)) :
                          atan(r[2], r[1])
        ν = u
    end
    return (a=a, e=e, i=i, Ω=Ω, ω=ω, ν=ν, u=u, p=p)
end

"""
    rv2coe_2pi(r, v, μ) — identical to rv2coe but all angles returned in [0, 2π).
    Avoids ±180° branch-cut jumps in time-series difference plots.
"""
function rv2coe_2pi(r::SVector{3,Float64}, v::SVector{3,Float64}, μ::Float64; tol=1e-10)
    rnorm = norm(r); v2 = dot(v,v); h = cross(r,v); hnorm = norm(h)
    i  = acos(clamp(h[3]/(hnorm+tol), -1.0, 1.0))
    n  = cross(Ẑ, h); nnorm = norm(n)
    evec = (cross(v,h)/μ) - r/(rnorm+tol); e = norm(evec)
    ε  = 0.5*v2 - μ/rnorm
    a  = -μ/(2ε)
    p  = hnorm^2/μ
    Ω  = nnorm > tol ? mod(atan(n[2], n[1]), 2π) : 0.0
    if e > 1e-8 && nnorm > tol
        ω = mod(atan(dot(Ẑ,cross(n,evec))/(nnorm*e+tol), dot(n,evec)/(nnorm*e+tol)), 2π)
        ν = mod(atan(dot(Ẑ,cross(evec,r))/(e*rnorm+tol), dot(evec,r)/(e*rnorm+tol)), 2π)
        u = NaN
    else
        ω = 0.0
        u = nnorm > tol ? mod(atan(dot(Ẑ,cross(n,r))/(nnorm*rnorm+tol), dot(n,r)/(nnorm*rnorm+tol)), 2π) :
                          mod(atan(r[2], r[1]), 2π)
        ν = u
    end
    return (a=a, e=e, i=i, Ω=Ω, ω=ω, ν=ν, u=u, p=p)
end

# ------------------------------------------------------------------
# RTN frame unit vectors
# ------------------------------------------------------------------
@inline function _rtn_basis(r::SVector{3, Float64}, v::SVector{3, Float64})
    rhat = r / norm(r)
    nhat = cross(r, v)
    nhat = nhat / norm(nhat)
    that = cross(nhat, rhat)
    return rhat, that, nhat
end

# rv to orbit elements
function _rv_to_elements(r::SVector{3, Float64}, v::SVector{3, Float64}, mu::Float64)
    rmag = norm(r)
    v2 = dot(v, v)
    h = cross(r, v)
    hmag = norm(h)
    n = cross(SVector{3, Float64}(0.0, 0.0, 1.0), h)
    nmag = norm(n)
    evec = cross(v, h) / mu - r / rmag
    e = norm(evec)
    a = -mu / (2.0 * (0.5 * v2 - mu / rmag))
    inc = acos(clamp(h[3] / hmag, -1.0, 1.0))
    raan = nmag <= 1e-12 ? 0.0 : mod(atan(n[2], n[1]), 2pi)
    return (a=a, e=e, i=inc, raan=raan)
end
