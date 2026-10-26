"""
    j2_secular_rates(a, e, i, planet) -> (Ωdot, ωdot)

First-order secular drift rates for the J2-only zonal gravity model, using the
standard expressions given in Vallado and Montenbruck & Gill.

Returns the right ascension of ascending node drift `Ωdot` and argument of
periapsis drift `ωdot` in rad/s.
"""
@inline function j2_secular_rates(
    a::Float64,
    e::Float64,
    i::Float64,
    planet
)::Tuple{Float64, Float64}
    if !isfinite(planet.J2) || planet.J2 == 0.0 || !isfinite(a) || a <= 0.0 || !isfinite(e) || e < 0.0 || e >= 1.0
        return 0.0, 0.0
    end
    p = a * (1.0 - e^2)
    if !isfinite(p) || p <= 0.0
        return 0.0, 0.0
    end
    n = sqrt(planet.μ / a^3)
    if !isfinite(n) || n <= 0.0
        return 0.0, 0.0
    end
    scale = planet.J2 * (planet.Rp_e / p)^2
    Ωdot = -1.5 * n * scale * cos(i)
    ωdot = 0.75 * n * scale * (5.0 * cos(i)^2 - 1.0)
    return Ωdot, ωdot
end
