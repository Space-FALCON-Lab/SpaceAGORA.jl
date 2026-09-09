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
    u1, u2, u3, u4 = u
    up1, up2, up3, up4 = u_prime
    r = dot(u, u)
    r > eps(Float64) || throw(ArgumentError("KS velocity is undefined at the origin."))
    return SVector(
        (2 / r) * (u1 * up1 - u2 * up2 - u3 * up3 + u4 * up4),
        (2 / r) * (u2 * up1 + u1 * up2 - u4 * up3 - u3 * up4),
        (2 / r) * (u3 * up1 + u4 * up2 + u1 * up3 + u2 * up4),
    )
end

function _ks_coordinate_from_position(x)
    x = Float64.(x)
    r = norm(x)
    r > eps(Float64) || throw(ArgumentError("KS coordinates are undefined at the origin."))
    if r + x[1] >= r - x[1]
        u1 = sqrt(max(0.0, 0.5 * (r + x[1])))
        denom = 2.0 * u1
        return [u1, x[2] / denom, x[3] / denom, 0.0]
    end
    u2 = sqrt(max(0.0, 0.5 * (r - x[1])))
    denom = 2.0 * u2
    return [x[2] / denom, u2, 0.0, x[3] / denom]
end

function _ks_derivative_from_velocity(velocity, u)
    u1, u2, u3, u4 = u
    xd1, xd2, xd3 = velocity
    return [
        0.5 * (u1 * xd1 + u2 * xd2 + u3 * xd3),
        0.5 * (-u2 * xd1 + u1 * xd2 + u4 * xd3),
        0.5 * (-u3 * xd1 - u4 * xd2 + u1 * xd3),
        0.5 * (u4 * xd1 - u3 * xd2 + u2 * xd3),
    ]
end

function _ks_L(p)
    return @SMatrix [
        p[1] -p[2] -p[3] p[4]
        p[2] p[1] -p[4] -p[3]
        p[3] p[4] p[1] p[2]
        p[4] -p[3] p[2] -p[1]
    ]
end

function _ks_lambda(p)
    return @SMatrix [
        p[1] -p[2] -p[3] p[4]
        p[2] p[1] -p[4] -p[3]
        p[3] p[4] p[1] p[2]
    ]
end

function _ks_phi(omega, delta_s)
    I4 = Matrix(I, 4, 4)
    if abs(omega) <= eps(Float64)
        return [I4 delta_s * I4; zeros(4, 4) I4]
    end
    return [
        cos(omega * delta_s) * I4 (sin(omega * delta_s) / omega) * I4
        (-omega * sin(omega * delta_s) * I4) cos(omega * delta_s) * I4
    ]
end
