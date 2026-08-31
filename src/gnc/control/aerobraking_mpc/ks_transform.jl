#=
"""
    KS transformation equations.

    These equations are kept from the supplied source codes so the SpaceAGORA
    implementation stays traceable to the original KS formulation.
"""
=#

function x_from_u(u)
    return SVector( u[1]^2 - u[2]^2 - u[3]^2 + u[4]^2,
                     2*(u[1]*u[2] - u[3]*u[4]),
                     2*(u[1]*u[3] + u[2]*u[4]) )
end

function xdot_from_u(u,u_prime)
    up = u_prime
    u1,u2,u3,u4 = u
    up1,up2,up3,up4 = up
    r = dot(u,u)
    x1d = (2/r)*(u1*up1 - u2*up2 - u3*up3 + u4*up4)
    x2d = (2/r)*(u2*up1 + u1*up2 - u4*up3 - u3*up4)
    x3d = (2/r)*(u3*up1 + u4*up2 + u1*up3 + u2*up4)
    return SVector(x1d,x2d,x3d)
end

function u_from_x(x)
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

function uprime_from_xdot(xdot, u)
    u1,u2,u3,u4 = u
    xd1,xd2,xd3 = xdot
    up1 = 0.5*( u1*xd1 + u2*xd2 + u3*xd3 )
    up2 = 0.5*(-u2*xd1 + u1*xd2 + u4*xd3 )
    up3 = 0.5*(-u3*xd1 - u4*xd2 + u1*xd3 )
    up4 = 0.5*( u4*xd1 - u3*xd2 + u2*xd3 )
    return [up1; up2; up3; up4]
end

function ks_L(p)
    return @SMatrix [p[1] -p[2] -p[3]  p[4];
                     p[2]  p[1] -p[4] -p[3];
                     p[3]  p[4]  p[1]  p[2];
                     p[4] -p[3]  p[2] -p[1]]
end

function Λ(p)
    return @SMatrix [p[1] -p[2] -p[3]  p[4];
                     p[2]  p[1] -p[4] -p[3];
                     p[3]  p[4]  p[1]  p[2]]
end

function Φ(ω, Δs)
    I4 = Matrix(I, 4, 4)
    return [cos(ω*Δs)*I4 (sin(ω*Δs)/ω)*I4;
            (-ω*sin(ω*Δs)*I4) cos(ω*Δs)*I4]
end
