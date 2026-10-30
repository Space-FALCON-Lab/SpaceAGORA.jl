module PeriapsisCorridorMeanElements

using LinearAlgebra

export mean_to_osculating, osculating_to_mean, elements_to_cartesian, angle_difference

"""
    PeriapsisCorridorMeanElements

Local, example-scoped Vallado–Kwok first-order J2 short-period element map.

Elements are named tuples `(a, e, i, raan, omega, nu)` in SI units and radians.
The map corrects only the J2 short-period terms; it is neither a general mean
element theory nor a replacement for Cartesian truth propagation. The sources are Vallado and Kwok's first-order J2 short-period
theory; see `docs/src/user/periapsis_corridor.md` for limitations.
"""

const _TWO_PI = 2π

_wrap_angle(x::Real) = mod(Float64(x), _TWO_PI)
angle_difference(a::Real, b::Real) = mod(Float64(a) - Float64(b) + π, _TWO_PI) - π

function _true_to_mean_anomaly(e, f)
    E = _wrap_angle(2atan(sqrt(1 - e) * sin(f / 2), sqrt(1 + e) * cos(f / 2)))
    return _wrap_angle(E - e * sin(E))
end

function _elements_tuple(elements)
    required = (:a, :e, :i, :raan, :omega, :nu)
    all(k -> hasproperty(elements, k), required) ||
        throw(ArgumentError("elements must provide a, e, i, raan, omega, and nu"))
    nt = NamedTuple{required}(Tuple(Float64(getproperty(elements, k)) for k in required))
    all(isfinite, values(nt)) || throw(ArgumentError("all orbital elements must be finite"))
    return nt
end

function _validate_constants(elements, mu, radius, j2)
    el = _elements_tuple(elements)
    isfinite(mu) && mu > 0 || throw(ArgumentError("mu must be finite and positive"))
    isfinite(radius) && radius > 0 || throw(ArgumentError("radius must be finite and positive"))
    isfinite(j2) || throw(ArgumentError("j2 must be finite"))
    el.a > radius || throw(ArgumentError("semimajor axis must exceed the reference radius"))
    0 < el.e < 1 || throw(ArgumentError("only noncircular elliptic orbits (0 < e < 1) are supported"))
    0 < el.i < π || throw(ArgumentError("equatorial and singular inclinations are not supported"))
    return el, Float64(mu), Float64(radius), Float64(j2)
end

"""Convert elliptic osculating elements to inertial Cartesian position and velocity."""
function elements_to_cartesian(elements; mu::Real)
    el = _elements_tuple(elements)
    isfinite(mu) && mu > 0 || throw(ArgumentError("mu must be finite and positive"))
    el.a > 0 && 0 < el.e < 1 || throw(ArgumentError("Cartesian conversion requires an elliptic orbit"))
    0 ≤ el.i ≤ π || throw(ArgumentError("inclination must lie in [0, π]"))
    p = el.a * (1 - el.e^2)
    nu = el.nu
    radius = p / (1 + el.e * cos(nu))
    speed_scale = sqrt(Float64(mu) / p)
    r_pqw = (radius * cos(nu), radius * sin(nu), 0.0)
    v_pqw = (-speed_scale * sin(nu), speed_scale * (el.e + cos(nu)), 0.0)

    cO, sO = cos(el.raan), sin(el.raan)
    cw, sw = cos(el.omega), sin(el.omega)
    ci, si = cos(el.i), sin(el.i)
    R11, R12 = cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci
    R21, R22 = sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci
    R31, R32 = sw * si, cw * si
    r = (R11 * r_pqw[1] + R12 * r_pqw[2],
         R21 * r_pqw[1] + R22 * r_pqw[2],
         R31 * r_pqw[1] + R32 * r_pqw[2])
    v = (R11 * v_pqw[1] + R12 * v_pqw[2],
         R21 * v_pqw[1] + R22 * v_pqw[2],
         R31 * v_pqw[1] + R32 * v_pqw[2])
    return r, v
end

"""
    mean_to_osculating(elements; mu, radius, j2)

Apply the first-order Vallado–Kwok J2 short-period map at the input mean true
anomaly. Supports elliptic, noncircular, non-equatorial classical elements
within the map's mathematical domain. This coordinate conversion does not
certify corridor safety; `corridor_safety_certified` is always false.
"""
function mean_to_osculating(elements; mu::Real, radius::Real, j2::Real)
    el, μ, R, J₂ = _validate_constants(elements, mu, radius, j2)
    a, e, inc, Ω, ω, f = el.a, el.e, el.i, _wrap_angle(el.raan), _wrap_angle(el.omega), _wrap_angle(el.nu)
    e2 = e * e
    p = a * (1 - e2)
    p2 = p * p
    M = _true_to_mean_anomaly(e, f)
    KJ₂ = J₂ * R * R

    si, ci = sin(inc), cos(inc)
    sf, cf = sin(f), cos(f)
    u = ω + f
    s2u, c2u = sincos(2u)
    swf, cwf = sincos(2ω + f)
    sw3f, cw3f = sincos(2ω + 3f)
    si2, ci2 = si * si, ci * ci
    ecosf, esinf = e * cf, e * sf
    aux1 = 3c2u + 3e * cwf + e * cw3f
    aux2 = sqrt(1 - e2)
    aux3 = 3ci2 - 1
    aux4 = -aux3 / (1 + aux2)

    # Adapted from SatelliteToolboxPropagators.jl 1.2.0, _osculating_elements (MIT).
    δi = KJ₂ * si * ci / (4p2) * aux1
    δp = KJ₂ * si2 / (2p) * aux1
    Δf_M = rem2pi(f - M, RoundNearest)
    δΩ = -KJ₂ * ci / (4p2) * (
        6 * (Δf_M + esinf) - 3s2u - 3e * swf - e * sw3f
    )
    δr = -KJ₂ / (4p) * (
        aux3 * (2aux2 / (1 + ecosf) + ecosf / (1 + aux2) + 1) - si2 * c2u
    )
    δrdot = KJ₂ * sqrt(μ) / (4 * sqrt(p^5)) * (
        aux3 * esinf * (aux2 + ((1 + ecosf)^2) / (1 + aux2)) -
        2si2 * (1 + ecosf)^2 * s2u
    )
    δu = KJ₂ / (8p2) * (
        (6 - 30ci2) * Δf_M +
        4esinf * (1 - 6ci2 + aux4) +
        aux4 * e2 * sin(2f) +
        (5ci2 - 2) * (2e) * swf +
        (7ci2 - 1) * s2u +
        2ci2 * e * sw3f
    )

    r = p / (1 + ecosf)
    rdot = sqrt(μ / p) * esinf
    rosc = r + δr
    rdot_osc = rdot + δrdot
    posc = p + δp
    A = posc / rosc - 1
    B = sqrt(posc / μ) * rdot_osc
    eosc = hypot(A, B)
    eosc < 1 || throw(DomainError(eosc, "J2 short-period map produced a non-elliptic state"))
    aosc = posc / (1 - eosc^2)
    iosc = inc + δi
    Ωosc = Ω + δΩ
    uosc = u + δu
    fosc = atan(B, A)
    ωosc = uosc - fosc

    result = (a=aosc, e=eosc, i=iosc, raan=_wrap_angle(Ωosc),
              omega=_wrap_angle(ωosc), nu=_wrap_angle(fosc))
    return (elements=result, converged=true, iterations=0, residual_norm=0.0,
            corridor_safety_certified=false)
end

_scaled_vector(elements, radius) = [elements.a / radius, elements.e, elements.i,
                                    _wrap_angle(elements.raan), _wrap_angle(elements.omega),
                                    _wrap_angle(elements.nu)]

function _vector_elements(x, radius)
    return (a=x[1] * radius, e=x[2], i=x[3], raan=_wrap_angle(x[4]),
            omega=_wrap_angle(x[5]), nu=_wrap_angle(x[6]))
end

function _scaled_residual(x, target, mu, radius, j2)
    mapped = mean_to_osculating(_vector_elements(x, radius); mu, radius, j2).elements
    y = _scaled_vector(mapped, radius)
    return [y[1] - target[1], y[2] - target[2], y[3] - target[3],
            angle_difference(y[4], target[4]), angle_difference(y[5], target[5]),
            angle_difference(y[6], target[6])]
end

function _inverse_result(x, radius, iteration, residual_norm)
    mean_elements = _vector_elements(x, radius)
    return (elements=mean_elements, converged=true, iterations=iteration,
            residual_norm, corridor_safety_certified=false)
end

"""
    osculating_to_mean(elements; mu, radius, j2, initial_guess=elements,
                       tolerance=1e-12, max_iterations=30)

Numerically invert [`mean_to_osculating`](@ref) with a damped finite-difference
Newton solve. The returned residual is dimensionless (`a` residual scaled by
`radius`; remaining residuals are in radians or unitless). Failure to converge
throws an error and never returns a partial estimate. This coordinate inversion
does not certify corridor safety; `corridor_safety_certified` is always false.
"""
function osculating_to_mean(elements; mu::Real, radius::Real, j2::Real,
                            initial_guess=elements, tolerance::Real=1e-12,
                            max_iterations::Integer=30)
    target_elements, μ, R, J₂ = _validate_constants(elements, mu, radius, j2)
    guess, _, _, _ = _validate_constants(initial_guess, mu, radius, j2)
    isfinite(tolerance) && tolerance > 0 || throw(ArgumentError("tolerance must be positive and finite"))
    max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    target = _scaled_vector(target_elements, R)
    x = _scaled_vector(guess, R)
    residual = _scaled_residual(x, target, μ, R, J₂)
    norm_residual = norm(residual, Inf)
    norm_residual ≤ tolerance && return _inverse_result(x, R, 0, norm_residual)

    for iteration in 1:max_iterations
        jacobian = Matrix{Float64}(undef, 6, 6)
        for column in 1:6
            step = 1e-7 * max(abs(x[column]), 1.0)
            xp, xm = copy(x), copy(x)
            xp[column] += step
            xm[column] -= step
            jacobian[:, column] = (
                _scaled_residual(xp, target, μ, R, J₂) -
                _scaled_residual(xm, target, μ, R, J₂)
            ) / (2step)
        end
        correction = try
            -(jacobian \ residual)
        catch err
            throw(ErrorException("osculating-to-mean Newton solve failed at iteration $iteration: $(sprint(showerror, err))"))
        end
        all(isfinite, correction) || throw(ErrorException("osculating-to-mean Newton step became non-finite"))

        accepted = false
        α = 1.0
        for _ in 1:20
            candidate = x .+ α .* correction
            try
                candidate_residual = _scaled_residual(candidate, target, μ, R, J₂)
                candidate_norm = norm(candidate_residual, Inf)
                if candidate_norm < norm_residual
                    x, residual, norm_residual = candidate, candidate_residual, candidate_norm
                    accepted = true
                    break
                end
            catch err
                err isa ArgumentError || err isa DomainError || rethrow()
            end
            α /= 2
        end
        accepted || throw(ErrorException("osculating-to-mean Newton solve stalled at iteration $iteration (residual=$norm_residual)"))
        norm_residual ≤ tolerance && return _inverse_result(x, R, iteration, norm_residual)
    end
    throw(ErrorException("osculating-to-mean conversion did not converge in $max_iterations iterations (residual=$norm_residual)"))
end

end # module PeriapsisCorridorMeanElements

# Adapted from SatelliteToolboxPropagators.jl, Copyright (c) 2017-2026 Ronan
# Arraes Jardim Chagas. Permission is hereby granted, free of charge, to any
# person obtaining a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so, subject to the following
# conditions: The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software. THE SOFTWARE
# IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR
# A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
# IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
