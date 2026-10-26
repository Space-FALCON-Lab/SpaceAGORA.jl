## Kustaanheimo-Stiefel (KS) <-> Cartesian maps
#
# With the 3×4 KS matrix Λ(u), the position is r = Λ(u) u (|r| = |u|²), the
# physical-time velocity is v = (2/|u|²) Λ(u) u′, where u′ = du/ds is the
# derivative in the Sundman fictitious time s (dt/ds = |r|), and conversely
# u′ = ½ Λ(u)ᵀ v. The 4×3 matrix Λ(u)ᵀ also lifts Cartesian accelerations into
# the KS equations of motion.

@inline _ks_coordinate(u) = SVector{4}(u[1], u[2], u[3], u[4])

"""Return the 3×4 KS matrix Λ(u), so that `r = Λ(u) u`."""
@inline function ks_lambda_matrix(u)
    return @SMatrix [
        u[1] -u[2] -u[3]  u[4]
        u[2]  u[1] -u[4] -u[3]
        u[3]  u[4]  u[1]  u[2]
    ]
end

"""Return the Cartesian position of a four-component KS coordinate."""
function ks_position(u)
    u4 = _ks_coordinate(u)
    return ks_lambda_matrix(u4) * u4
end

"""Return the physical-time Cartesian velocity of a KS coordinate and its derivative."""
function ks_velocity(u, u_prime)
    radius = dot(u, u)
    radius > eps(Float64) || throw(ArgumentError("KS velocity is undefined at the origin."))
    return ks_velocity_with_radius(u, u_prime, radius)
end

@inline ks_velocity_with_radius(u, u_prime, radius) =
    (2 / radius) * (ks_lambda_matrix(u) * _ks_coordinate(u_prime))

# One of the KS coordinates of a position (the KS map is one-to-many). The
# branch keeps the divisor away from zero: it sets u₄ = 0 when x₁ ≥ 0 and u₃ = 0
# otherwise.
function cartesian_position_to_ks_coordinate(x)
    x = SVector{3,Float64}(x)
    r = norm(x)
    r > eps(Float64) || throw(ArgumentError("KS coordinates are undefined at the origin."))
    if x[1] >= 0.0
        u1 = sqrt(max(0.0, 0.5 * (r + x[1])))
        return SVector(u1, x[2] / (2.0 * u1), x[3] / (2.0 * u1), 0.0)
    end
    u2 = sqrt(max(0.0, 0.5 * (r - x[1])))
    return SVector(x[2] / (2.0 * u2), u2, 0.0, x[3] / (2.0 * u2))
end

# u′ = ½ Λ(u)ᵀ v, which also satisfies the KS bilinear constraint.
cartesian_velocity_to_ks_derivative(velocity, u) =
    0.5 * (transpose(ks_lambda_matrix(u)) * SVector{3}(velocity[1], velocity[2], velocity[3]))
