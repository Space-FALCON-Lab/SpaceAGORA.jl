module ClassicalElementSensitivity

using LinearAlgebra
export classical_element_jacobian, classical_element_values, cross_matrix

cross_matrix(x) = [0.0 -x[3] x[2]; x[3] 0.0 -x[1]; -x[2] x[1] 0.0]

function classical_element_values(r, v; mu)
    radius = norm(r)
    h = cross(r, v)
    node = cross([0.0, 0.0, 1.0], h)
    ev = cross(v, h) / mu - r / radius
    a = inv(2 / radius - dot(v, v) / mu)
    inclination = atan(norm(h[1:2]), h[3])
    raan = atan(node[2], node[1])
    omega = atan(dot(h / norm(h), cross(node, ev)), dot(node, ev))
    return [a, norm(ev), inclination, raan, omega]
end

"""Analytical d[a,e,i,Omega,omega]/d[r,v]; km, km/s, radians."""
function classical_element_jacobian(r, v; mu, include_true_anomaly=false)
    radius = norm(r)
    h = cross(r, v)
    hm = norm(h)
    node = cross([0.0, 0.0, 1.0], h)
    nm = norm(node)
    ev = cross(v, h) / mu - r / radius
    eccentricity = norm(ev)
    eccentricity > 1e-10 && nm / hm > 1e-10 || throw(ArgumentError(
        "Classical e/omega and i/Omega derivatives need a nonequatorial, noncircular orbit.",
    ))
    eye = Matrix{Float64}(I, 3, 3)
    dr = hcat(eye, zeros(3, 3))
    dv = hcat(zeros(3, 3), eye)
    dh = -cross_matrix(v) * dr + cross_matrix(r) * dv
    dn = cross_matrix([0.0, 0.0, 1.0]) * dh
    de = (-cross_matrix(h) * dv + cross_matrix(v) * dh) / mu -
        (eye / radius - r * r' / radius^3) * dr
    hh, nn, ee = h / hm, node / nm, ev / eccentricity
    dhh = (eye - hh * hh') * dh / hm
    dnn = (eye - nn * nn') * dn / nm
    dee = (eye - ee * ee') * de / eccentricity
    a = inv(2 / radius - dot(v, v) / mu)
    da = 2a^2 / mu * (mu * r' / radius^3 * dr + v' * dv)
    di = (h[3] * (h[1] * dh[1:1, :] + h[2] * dh[2:2, :]) / nm -
        nm * dh[3:3, :]) / hm^2
    dOmega = (node[1] * dn[2:2, :] - node[2] * dn[1:1, :]) / nm^2
    x = dot(nn, ee)
    y = dot(hh, cross(nn, ee))
    dx = ee' * dnn + nn' * dee
    dy = cross(nn, ee)' * dhh + hh' * (-cross_matrix(ee) * dnn + cross_matrix(nn) * dee)
    domega = (x * dy - y * dx) / (x^2 + y^2)
    result = vcat(da, ee' * de, di, dOmega, domega)
    if include_true_anomaly
        xnu, ynu = dot(ev, r), dot(hh, cross(ev, r))
        dxnu = r' * de + ev' * dr
        dynu = cross(ev, r)' * dhh + hh' * (-cross_matrix(r) * de + cross_matrix(ev) * dr)
        dnu = (xnu * dynu - ynu * dxnu) / (xnu^2 + ynu^2)
        return vcat(result, dnu)
    end
    return result
end

end
