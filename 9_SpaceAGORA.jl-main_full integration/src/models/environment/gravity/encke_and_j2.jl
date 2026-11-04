"""
    encke_point_mass_difference(μ, r, ρ)

Point-mass gravity difference `g(r + ρ) - g(r)` with `g(r) = -μ r / |r|^3`, evaluated without forming
the two large terms: `-μ/|r|^3 (ρ - f(σ)(r + ρ))`, `σ = ρ·(2r + ρ)/|r|^2`,
`f(σ) = σ(3 + 3σ + σ^2) / ([1 + (1 + σ)^(3/2)] (1 + σ)^(3/2))` (Encke's method). The factor `σ` is
exact in the numerator, so the result keeps full relative precision as `ρ → 0` (zero for `ρ = 0`).
Generic over the element type and allocation-free. Internal to the gravity code.
"""
@inline function encke_point_mass_difference(μ, r, ρ)
    r2 = dot(r, r)
    σ = dot(ρ, 2r + ρ) / r2
    s3 = sqrt(1 + σ)^3
    f = σ * (3 + σ * (3 + σ)) / ((1 + s3) * s3)
    return -μ / (r2 * sqrt(r2)) * (ρ - f * (r + ρ))
end

# J2 term only (the second term of `_inverse_squared_j2_gravity_accel`), planet-fixed frame.
@inline function _j2_only_gravity_accel(pos_pp::SVector{3, Float64}, planet)::SVector{3, Float64}
    r = norm(pos_pp)
    μ = Float64(planet.μ)
    J2 = Float64(planet.J2)
    Rp_m = Float64(planet.Rp_e)
    x, y, z = pos_pp
    r_squared = r^2
    j2_term = SVector{3, Float64}(
        x / r * (5 * z^2 / r_squared - 1),
        y / r * (5 * z^2 / r_squared - 1),
        z / r * (5 * z^2 / r_squared - 3),
    )
    return 3 / 2 * J2 * μ * Rp_m^2 / r^4 * j2_term
end
