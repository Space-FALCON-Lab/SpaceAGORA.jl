module OrbitalGeometry
using StaticArrays

# Inertial axes to perifocal axes for radian orbital angles. Keep the original
# component expressions and their evaluation order; callers own vector scaling
# and any direction normalization. This is an internal, non-exported owner.
@inline function _perifocal_from_inertial(i::Float64, ω::Float64, Ω::Float64)::SMatrix{3, 3, Float64, 9}
    return @SMatrix [
        -sin(Ω)*cos(i)*sin(ω)+cos(Ω)*cos(ω)  cos(Ω)*cos(i)*sin(ω)+sin(Ω)*cos(ω)  sin(i)*sin(ω)
        -sin(Ω)*cos(i)*cos(ω)-cos(Ω)*sin(ω)  cos(Ω)*cos(i)*cos(ω)-sin(Ω)*sin(ω)  sin(i)*cos(ω)
         sin(Ω)*sin(i)                       -cos(Ω)*sin(i)                       cos(i)
    ]
end

# The true anomaly (rad, in [0, 2π)) at mean anomaly `M` (rad) on an ellipse of
# eccentricity 0 <= e < 1. The apsides are exact. Elsewhere Kepler's equation
# M = E - e sin E is solved by Newton's method, started at M (or at π for a very
# eccentric orbit, where M is a poor start). Its result stands where the step has
# converged and 1 - e cos E >= 0.05, which holds for every e <= 0.95: there the
# rounding of E - e sin E - M moves E by a few ulps at most. Near a parabolic orbit,
# with E (or 2π - E) small, Newton crawls toward a nearly triple root and that
# rounding swamps E: the loop can stop on a wrong E or not stop at all, and
# sqrt(1 - e) then amplifies the error. There `_eccentric_anomaly_safeguarded`
# solves the equation instead.
function _true_anomaly_from_mean(M::Float64, e::Float64)::Float64
    0.0 <= e < 1.0 || throw(ArgumentError("The mean anomaly needs an elliptic orbit, 0 <= e < 1; got e = $e."))
    isfinite(M) || throw(ArgumentError("The mean anomaly must be finite; got M = $M."))
    M = mod(M, 2pi)
    M == 0.0 && return 0.0
    M == Float64(pi) && return Float64(pi)
    E = e < 0.8 ? M : Float64(pi)
    for _ in 1:50
        step = (E - e * sin(E) - M) / (1.0 - e * cos(E))
        E -= step
        if abs(step) <= 1e-14   # Newton is quadratic: E is then exact to rounding
            1.0 - e * cos(E) >= 0.05 && return _true_anomaly_from_eccentric(E, e)
            break
        end
    end
    # Kepler's equation is symmetric about π: E(2π - M) = 2π - E(M).
    M < Float64(pi) && return _true_anomaly_from_eccentric(_eccentric_anomaly_safeguarded(M, e), e)
    ν = 2pi - _true_anomaly_from_eccentric(_eccentric_anomaly_safeguarded(2pi - M, e), e)
    return ν < 2pi ? ν : 0.0
end

@inline _true_anomaly_from_eccentric(E::Float64, e::Float64)::Float64 =
    mod(2.0 * atan(sqrt(1.0 + e) * sin(E / 2.0), sqrt(1.0 - e) * cos(E / 2.0)), 2pi)

# The eccentric anomaly in (0, π) at mean anomaly 0 < M < π. Newton's method stays
# inside a bracket of the root: f(E) = E - e sin E - M rises (f' = 1 - e cos E > 0)
# from f(0) < 0 to f(π) > 0, and a step that would leave the bracket bisects it
# instead. f is evaluated as (1 - e) E + e (E - sin E) - M, and f' as
# (1 - e) + 2 e sin²(E/2): near a parabolic orbit both terms are small and keep their
# precision, where E - e sin E would cancel. The start is the parabolic root
# cbrt(6M). It throws rather than return an unconverged E.
function _eccentric_anomaly_safeguarded(M::Float64, e::Float64)::Float64
    lo, hi = 0.0, Float64(pi)
    E = min(cbrt(6.0 * M), Float64(pi))
    for _ in 1:2000
        f = (1.0 - e) * E + e * _x_minus_sin(E) - M
        f == 0.0 && return E
        f < 0.0 ? (lo = E) : (hi = E)
        next = E - f / ((1.0 - e) + 2.0 * e * sin(E / 2.0)^2)
        lo < next < hi || (next = 0.5 * (lo + hi))
        (abs(next - E) <= 2.0 * eps(E) || hi <= nextfloat(lo)) && return next
        E = next
    end
    error("Kepler's equation did not converge for M = $M, e = $e.")
end

# x - sin x, by its series below 1 where the direct difference cancels.
@inline function _x_minus_sin(x::Float64)::Float64
    abs(x) >= 1.0 && return x - sin(x)
    x2 = x * x
    return x * x2 * evalpoly(x2, (1 / 6, -1 / 120, 1 / 5040, -1 / 362880, 1 / 39916800, -1 / 6227020800,
                                  1 / 1307674368000, -1 / 355687428096000))
end
end
