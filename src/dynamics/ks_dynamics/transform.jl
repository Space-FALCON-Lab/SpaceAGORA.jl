"""Return Cartesian position from a four-component KS coordinate."""
function ks_position(u)
    return SVector(
        u[1]^2 - u[2]^2 - u[3]^2 + u[4]^2,
        2 * (u[1] * u[2] - u[3] * u[4]),
        2 * (u[1] * u[3] + u[2] * u[4]),
    )
end

"""Return physical-time Cartesian velocity from KS position and derivative."""
function ks_velocity(u, u_prime)
    radius = dot(u, u)
    radius > eps(Float64) || throw(ArgumentError("KS velocity is undefined at the origin."))
    return ks_velocity_with_radius(u, u_prime, radius)
end

@inline function ks_velocity_with_radius(u, u_prime, radius)
    u1, u2, u3, u4 = u
    up1, up2, up3, up4 = u_prime
    return SVector(
        (2 / radius) * (u1 * up1 - u2 * up2 - u3 * up3 + u4 * up4),
        (2 / radius) * (u2 * up1 + u1 * up2 - u4 * up3 - u3 * up4),
        (2 / radius) * (u3 * up1 + u4 * up2 + u1 * up3 + u2 * up4),
    )
end

function cartesian_position_to_ks_coordinate(x)
    x = SVector{3,Float64}(x)
    r = norm(x)
    r > eps(Float64) || throw(ArgumentError("KS coordinates are undefined at the origin."))
    if r + x[1] >= r - x[1]
        u1 = sqrt(max(0.0, 0.5 * (r + x[1])))
        denom = 2.0 * u1
        return SVector(u1, x[2] / denom, x[3] / denom, 0.0)
    end
    u2 = sqrt(max(0.0, 0.5 * (r - x[1])))
    denom = 2.0 * u2
    return SVector(x[2] / denom, u2, 0.0, x[3] / denom)
end

function cartesian_velocity_to_ks_derivative(velocity, u)
    u1, u2, u3, u4 = u
    xd1, xd2, xd3 = velocity
    return SVector(
        0.5 * (u1 * xd1 + u2 * xd2 + u3 * xd3),
        0.5 * (-u2 * xd1 + u1 * xd2 + u4 * xd3),
        0.5 * (-u3 * xd1 - u4 * xd2 + u1 * xd3),
        0.5 * (u4 * xd1 - u3 * xd2 + u2 * xd3),
    )
end

function ks_lambda_matrix(p)
    return @SMatrix [
        p[1] -p[2] -p[3] p[4]
        p[2] p[1] -p[4] -p[3]
        p[3] p[4] p[1] p[2]
    ]
end
