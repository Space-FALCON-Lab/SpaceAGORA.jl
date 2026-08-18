# Named-tuple orbital element converters for the ORACLE analysis module.
# rv2coe / rv2coe_2pi take μ::Float64 and return named tuples — a different
# interface from Kuang-ver2's rvtoorbitalelement(r,v,planet) → SVector{6}.
# Both are needed because the analysis code is built on named-tuple field access.

using LinearAlgebra
using StaticArrays

const _Ẑ_OE = SVector(0.0, 0.0, 1.0)

export rv2coe, rv2coe_2pi

"""
    rv2coe(r, v, μ) → NamedTuple(a, e, i, Ω, ω, ν, u, p)

Convert ECI position/velocity to classical orbital elements.
Robust to circular and equatorial orbits (returns u for circular, ν=NaN).
"""
function rv2coe(r::SVector{3,Float64}, v::SVector{3,Float64}, μ::Float64; tol=1e-10)
    rnorm = norm(r); v2 = dot(v,v); h = cross(r,v); hnorm = norm(h)
    i  = acos(clamp(h[3]/(hnorm+tol), -1.0, 1.0))
    n  = cross(_Ẑ_OE, h); nnorm = norm(n)
    evec = (cross(v,h)/μ) - r/(rnorm+tol); e = norm(evec)
    ε  = 0.5*v2 - μ/rnorm
    a  = -μ/(2ε)
    p  = hnorm^2/μ
    Ω  = nnorm > tol ? atan(n[2], n[1]) : 0.0
    if e > 1e-8 && nnorm > tol
        ω = atan(dot(_Ẑ_OE, cross(n,evec))/(nnorm*e+tol), dot(n,evec)/(nnorm*e+tol))
        ν = atan(dot(_Ẑ_OE, cross(evec,r))/(e*rnorm+tol), dot(evec,r)/(e*rnorm+tol))
        u = NaN
    else
        ω = 0.0
        u = nnorm > tol ? atan(dot(_Ẑ_OE,cross(n,r))/(nnorm*rnorm+tol), dot(n,r)/(nnorm*rnorm+tol)) :
                          atan(r[2], r[1])
        ν = u
    end
    return (a=a, e=e, i=i, Ω=Ω, ω=ω, ν=ν, u=u, p=p)
end

"""
    rv2coe_2pi(r, v, μ) — identical to rv2coe but all angles in [0, 2π).
"""
function rv2coe_2pi(r::SVector{3,Float64}, v::SVector{3,Float64}, μ::Float64; tol=1e-10)
    rnorm = norm(r); v2 = dot(v,v); h = cross(r,v); hnorm = norm(h)
    i  = acos(clamp(h[3]/(hnorm+tol), -1.0, 1.0))
    n  = cross(_Ẑ_OE, h); nnorm = norm(n)
    evec = (cross(v,h)/μ) - r/(rnorm+tol); e = norm(evec)
    ε  = 0.5*v2 - μ/rnorm
    a  = -μ/(2ε)
    p  = hnorm^2/μ
    Ω  = nnorm > tol ? mod(atan(n[2], n[1]), 2π) : 0.0
    if e > 1e-8 && nnorm > tol
        ω = mod(atan(dot(_Ẑ_OE,cross(n,evec))/(nnorm*e+tol), dot(n,evec)/(nnorm*e+tol)), 2π)
        ν = mod(atan(dot(_Ẑ_OE,cross(evec,r))/(e*rnorm+tol), dot(evec,r)/(e*rnorm+tol)), 2π)
        u = NaN
    else
        ω = 0.0
        u = nnorm > tol ? mod(atan(dot(_Ẑ_OE,cross(n,r))/(nnorm*rnorm+tol), dot(n,r)/(nnorm*rnorm+tol)), 2π) :
                          mod(atan(r[2], r[1]), 2π)
        ν = u
    end
    return (a=a, e=e, i=i, Ω=Ω, ω=ω, ν=ν, u=u, p=p)
end
