# oracle_closed_forms.jl
#
# Closed forms of Sections IV and V of the ORACLE SciTech 2027 paper:
# one target and a ring of N_h equally spaced helpers on a circular orbit.
# The file needs no package (only Base and the Statistics standard library).
#
# ------------------------------------------------------------------------------------------------
# HOW TO USE IT
# ------------------------------------------------------------------------------------------------
# 1. Put this file next to your plotting script and load it:
#
#        include("oracle_closed_forms.jl")
#
# 2. Set the parameters to the values of your code. The defaults are the values of Table 1 of the
#    paper. The one you must check is R_E, the Earth radius that your code uses to convert an
#    altitude into an orbit radius:
#
#        p = Params(R_E = 6378137.0)
#
# 3. Prediction for ONE case. The rule is :range, :amax or :imax.
#
#        r = predict(:amax, 300, 1050.0; inc_deg = 0.5, n_orbits = 100, p = p)
#        r.dvR, r.dvT, r.dvNp, r.dvNm      # velocity changes, m/s  (dvNp > 0, dvNm < 0)
#        r.da_km, r.di_deg, r.dOmega_deg   # changes of semi-major axis, inclination and node
#
# 4. Prediction for ONE CURVE of your figures (one rule, one target altitude, one inclination),
#    against the number of helpers. Use many values of N_h so that the line is smooth:
#
#        Nh    = [1; 25:25:800]
#        curve = [predict(:amax, N, 1050.0; inc_deg = 0.5, p = p) for N in Nh]
#        dvT   = [c.dvT for c in curve]        # plot dvT against Nh as a line
#
# 5. ALL the predictions for the grid of the paper in one CSV file (3 rules, 3 target altitudes,
#    3 inclinations, N_h = 1, 25, 50, ..., 800). Load it in your plotting script as a table:
#
#        prediction_table("oracle_predictions.csv"; p = p)
#
# 6. What to put on each panel of the two figures (simulation as markers, prediction as a line of
#    the same color):
#
#        panel of the figure          field of predict
#        Delta v_R                    dvR
#        Delta v_T                    dvT
#        Delta v_N+                   dvNp
#        Delta v_N-                   dvNm
#        Delta a                      da_km
#        Delta i                      di_deg
#        Delta Omega (laser part)     dOmega_deg   (compare with the run minus the links-off run)
#        Delta e                      no prediction, simulation only
#
# 7. Long runs. Shell separation d = a_target - r_ring against time, Delta a-max rule:
#
#        t = (1:1700) .* (2pi / mean_motion(r_h(p), p))     # one value per orbit, s
#        d = d_numeric(t, 0.0, 800, p)                      # any number of helpers, m
#        d_dense(t, 0.0, p), d_sparse(t, 50.0e3, 100, p)    # the two closed forms of the paper
#
# ------------------------------------------------------------------------------------------------
# WHEN THE PREDICTION DOES NOT APPLY
# ------------------------------------------------------------------------------------------------
# predict assumes that the helpers drift past the target, which needs a shell separation.
# For a target that starts ON the ring (target altitude = ring altitude):
#   - Delta a-max rule with fewer than 2 pi r_h / L_max helpers (232 at 1000 km): the result
#     depends on where the target starts between two helpers. Use
#         leapfrog(N_h, b0; n_orbits = 100, p = p)
#     with b0 the initial along-track distance (m) of the nearest helper BEHIND the target.
#     For all four RTN panels, also pass inc_deg and u0_deg (initial argument of latitude)
#     to use the phase-resolved, numerical extension instead of the scalar event formula.
#     If b0 > L_max no link ever forms and the velocity change is zero.
#   - Delta i-max rule: predict gives the average over all the starting positions. One simulation
#     can differ from it by up to a factor of two with 300 helpers. inc_change_fixed_offset gives
#     the lower bound for a given starting offset.
#
# ------------------------------------------------------------------------------------------------
# NOTATION (same as the paper), SI units (m, s, rad) unless the name ends in _km or _deg
# ------------------------------------------------------------------------------------------------
#   r_h      radius of the ring                 r_t = r_h + d   radius of the target orbit
#   d        shell separation, positive when the target is above the ring
#   s        arc between two adjacent helpers, 2 pi r_h / N_h
#   L_max    link range                         ell = sqrt(L_max^2 - d^2), half of the in-range arc
#   a_OCL    acceleration of one link, beta B P / (c m)
#   inc      relative inclination Delta i
#   u        argument of latitude of the target, from its ascending node on the plane of the ring
#   eR, eT, eN   average projections of the link direction (link off counts as zero)

using Statistics

const C_LIGHT = 299792458.0

Base.@kwdef struct Params
    B::Float64      = 100.0            # magnification
    P::Float64      = 1.0e4            # power, W
    beta::Float64   = 1.0              # transmission factor
    m::Float64      = 227.0            # spacecraft mass, kg
    L_max::Float64  = 200.0e3          # link range, m
    h_ring::Float64 = 1000.0e3         # altitude of the ring, m
    R_E::Float64    = 6378137.0        # Earth radius used to convert altitude to radius, m
    mu::Float64     = 398600.4418e9    # gravitational parameter, m^3/s^2
end

a_OCL(p::Params) = p.beta * p.B * p.P / (C_LIGHT * p.m)
r_h(p::Params) = p.R_E + p.h_ring
spacing(N_h, p::Params) = 2pi * r_h(p) / N_h
mean_motion(r, p::Params) = sqrt(p.mu / r^3)

# ------------------------------------------------------------------------------------------------
# Section IV.C: coplanar target (Delta i = 0)
# ------------------------------------------------------------------------------------------------

"Range-based rule (nearest helper within range). Returns (eR, eT)."
function range_based(d, s, p::Params)
    D = abs(d)
    (D >= p.L_max || D == 0.0) && return (0.0, 0.0)
    ell = sqrt(p.L_max^2 - D^2)
    return ((2d / s) * asinh(min(s / 2, ell) / D), 0.0)
end

"Delta a-max rule (farthest helper behind the target, within range). Returns (eR, eT)."
function amax(d, s, p::Params)
    D = abs(d)
    D >= p.L_max && return (0.0, 0.0)
    ell = sqrt(p.L_max^2 - D^2)
    eT = s >= ell ? (p.L_max - D) / s : (p.L_max - sqrt(D^2 + (ell - s)^2)) / s
    eR = D == 0.0 ? 0.0 : (d / s) * (asinh(ell / D) - asinh((ell - min(s, ell)) / D))
    return (eR, eT)
end

"Helper count at which a helper is always within range. The Delta a-max curves bend at 2 N_sat."
N_sat(d, p::Params) = pi * r_h(p) / sqrt(p.L_max^2 - d^2)

# ------------------------------------------------------------------------------------------------
# Sections IV.D and V.C: target with a relative inclination, averages over the argument of latitude
# ------------------------------------------------------------------------------------------------

"""
    averages(rule, d, s, inc, p; nu = 4000)

Average projections for a shell separation d (m) and a relative inclination inc (rad).
Returns (eR, eT, eNp, eNm, k_inc, k_node) with
    eNp, eNm   average of the positive and of the negative normal projection (eNm <= 0)
    k_inc      <e_N cos u>,      so that  Delta i(T) - Delta i(0) = a_OCL T k_inc / (n r_t)
    k_node     <sin^2 u / rho>,  so that  Delta Omega            = a_OCL T k_node / n
"""
function averages(rule::Symbol, d, s, inc, p::Params; nu::Int = 4000)
    r_t = r_h(p) + d
    eR = 0.0; eT = 0.0; eNp = 0.0; eNm = 0.0; k_inc = 0.0; k_node = 0.0
    for k in 1:nu
        u = (k - 0.5) * 2pi / nu
        z = r_t * sin(inc) * sin(u)                 # height of the target above the plane of the ring
        D = max(sqrt(d^2 + z^2), 1.0e-3)            # transverse separation D(u)
        D >= p.L_max && continue                    # ring out of range at this u
        ell = sqrt(p.L_max^2 - D^2)                 # half of the in-range arc at this u
        if rule == :range || rule == :imax
            (rule == :imax && sin(u) * cos(u) <= 0.0) && continue      # link off in these quarters
            W = (2 / s) * asinh(min(s / 2, ell) / D)                   # average of 1/rho, nearest helper
        elseif rule == :amax
            W = (1 / s) * (asinh(ell / D) - asinh((ell - min(s, ell)) / D))
            eT += s >= ell ? (p.L_max - D) / s : (p.L_max - sqrt(D^2 + (ell - s)^2)) / s
        else
            error("rule must be :range, :amax or :imax")
        end
        eN = z * W
        eR += d * W
        eNp += max(eN, 0.0)
        eNm += min(eN, 0.0)
        k_inc += eN * cos(u)
        k_node += W * sin(u)^2
    end
    return (eR = eR / nu, eT = eT / nu, eNp = eNp / nu, eNm = eNm / nu, k_inc = k_inc / nu, k_node = k_node / nu)
end

# ------------------------------------------------------------------------------------------------
# Prediction of one simulation case
# ------------------------------------------------------------------------------------------------

"""
    predict(rule, N_h, h_t_km; inc_deg = 0.0, n_orbits = 100, p = Params(), dt = 1000.0)

Velocity changes and orbit changes of the target after n_orbits (orbits of the initial target orbit).
rule is :range, :amax or :imax. h_t_km is the initial altitude of the target in km.

The shell separation is advanced in time with d_dot = (2 a_OCL / n) eT, so the result includes the rise
of the target under the Delta a-max rule.

Returns (dvR, dvT, dvNp, dvNm, da_km, di_deg, dOmega_deg, ceiling), velocities in m/s.
ceiling is a_OCL T, the largest velocity change that the links can give in the time T.

Not included: under the Delta a-max rule the relative inclination also decreases at second order, by
about (Delta i / 2) dvT / v, with v the orbital speed.
"""
function predict(rule::Symbol, N_h, h_t_km; inc_deg = 0.0, n_orbits = 100, p::Params = Params(), dt = 1000.0)
    rule in (:range, :amax, :imax) || error("rule must be :range, :amax or :imax")
    s = spacing(N_h, p)
    inc = deg2rad(inc_deg)
    d = p.R_E + h_t_km * 1.0e3 - r_h(p)
    d0 = d
    T = n_orbits * 2pi / mean_motion(r_h(p) + d0, p)
    # only the Delta a-max rule changes the shell separation, so only this rule needs time steps
    steps = rule == :amax ? max(round(Int, T / dt), 1) : 1
    h = T / steps
    dvR = 0.0; dvT = 0.0; dvNp = 0.0; dvNm = 0.0; di = 0.0; dOm = 0.0
    for _ in 1:steps
        r_t = r_h(p) + d
        n = mean_motion(r_t, p)
        if inc == 0.0
            eR, eT = rule == :range ? range_based(d, s, p) : rule == :amax ? amax(d, s, p) : (0.0, 0.0)
            av = (eR = eR, eT = eT, eNp = 0.0, eNm = 0.0, k_inc = 0.0, k_node = 0.0)
        else
            av = averages(rule, d, s, inc, p; nu = rule == :amax ? 2000 : 20000)
        end
        dvR += a_OCL(p) * av.eR * h
        dvT += a_OCL(p) * av.eT * h
        dvNp += a_OCL(p) * av.eNp * h
        dvNm += a_OCL(p) * av.eNm * h
        di += a_OCL(p) * av.k_inc / (n * r_t) * h
        dOm += a_OCL(p) * av.k_node / n * h
        d += 2 * a_OCL(p) * av.eT / n * h
    end
    return (dvR = dvR, dvT = dvT, dvNp = dvNp, dvNm = dvNm, da_km = (d - d0) / 1.0e3,
            di_deg = rad2deg(di), dOmega_deg = rad2deg(dOm), ceiling = a_OCL(p) * T)
end

"""
    prediction_table(path; p = Params(), n_orbits = 100, altitudes_km, inclinations_deg, helpers)

Writes the predictions for the whole grid of the paper to a CSV file, one row per case.
Columns: rule, h_t_km, inc_deg, N_h, dvR, dvT, dvNp, dvNm, da_km, di_deg, dOmega_deg.
"""
function prediction_table(path::AbstractString; p::Params = Params(), n_orbits = 100,
                          altitudes_km = (1000.0, 1050.0, 1100.0), inclinations_deg = (0.0, 0.5, 1.0),
                          helpers = [1; 25:25:800])
    open(path, "w") do io
        println(io, "rule,h_t_km,inc_deg,N_h,dvR,dvT,dvNp,dvNm,da_km,di_deg,dOmega_deg")
        for rule in (:range, :amax, :imax), h_t in altitudes_km, inc in inclinations_deg, N_h in helpers
            r = predict(rule, N_h, h_t; inc_deg = inc, n_orbits = n_orbits, p = p)
            println(io, join((rule, h_t, inc, N_h, r.dvR, r.dvT, r.dvNp, r.dvNm, r.da_km, r.di_deg, r.dOmega_deg), ","))
        end
    end
    return path
end

# ------------------------------------------------------------------------------------------------
# Section V.A: shell separation against time under the Delta a-max rule (for the long runs)
# ------------------------------------------------------------------------------------------------

"Sparse ring (s >= ell), target above the ring: d(t) = L_max - (L_max - d0) exp(-t / t_c)."
function d_sparse(t, d0, N_h, p::Params)
    t_c = pi * r_h(p) * mean_motion(r_h(p), p) / (a_OCL(p) * N_h)
    L = p.L_max
    return L .- (L - d0) .* exp.(-t ./ t_c)
end

"Dense ring: d(t) = L_max sin(phi0 + 2 a_OCL t / (n L_max)). Upper bound for any number of helpers."
function d_dense(t, d0, p::Params)
    L = p.L_max
    rate = 2 * a_OCL(p) / (mean_motion(r_h(p), p) * L)
    return L .* sin.(min.(asin(d0 / L) .+ rate .* t, pi / 2))
end

"Shell separation at the increasing times t for any number of helpers (numerical integration of d_dot)."
function d_numeric(t, d0, N_h, p::Params; dt = 50.0)
    s = spacing(N_h, p)
    out = Float64[]
    d = float(d0); now = 0.0
    for t_end in t
        while now < t_end
            d += 2 * a_OCL(p) * amax(d, s, p)[2] / mean_motion(r_h(p) + d, p) * dt
            now += dt
        end
        push!(out, d)
    end
    return out
end

# ------------------------------------------------------------------------------------------------
# Section V.B: phase
# ------------------------------------------------------------------------------------------------

"Change of the argument of latitude (rad) from a radial and from an along-track velocity change after the time T."
phase_change(dvR, dvT, T, r_t, p::Params) = (-2 * dvR / (mean_motion(r_t, p) * r_t), -3 * dvT * T / (2 * r_t))

# ------------------------------------------------------------------------------------------------
# Section V.D: target that starts on the ring
# ------------------------------------------------------------------------------------------------

"""
    leapfrog(N_h, b0; n_orbits = 100, p = Params(), inc_deg = nothing, u0_deg = 0.0, dt = 1.0)

Delta a-max rule, target on the ring, nearest helper initially at the along-track distance b0 (m)
BEHIND the target. Returns (dvT in m/s, altitude gained in km, number of completed links).
Use it when N_h < 2 pi r_h / L_max. With more helpers use predict.

Passing `inc_deg` (including zero) selects the phase-resolved RTN extension. It returns
`dvR, dvT, dvNp, dvNm` (m/s), `da_km`, and `links` (completed active intervals). This
integrates a reduced circular-orbit model with fixed inclination, not the paper's scalar
closed form or the J2 trajectory simulator. `u0_deg` sets the starting orbital phase;
`dt` is the maximum step in seconds. All impulses and the shell rise use the same link.
"""
function leapfrog(N_h, b0; n_orbits = 100, p::Params = Params(), inc_deg = nothing, u0_deg = 0.0, dt = 1.0)
    !isnothing(inc_deg) && return leapfrog_rtn(N_h, b0, inc_deg, u0_deg, n_orbits, dt, p)
    n = mean_motion(r_h(p), p)
    T = n_orbits * 2pi / n
    s = spacing(N_h, p)
    lam = 4 * a_OCL(p) / (3 * n^2)
    b0 > p.L_max && return (dvT = 0.0, da_km = 0.0, links = 0)
    t = 0.0; d = 0.0; rho = float(b0); k = 0
    while true
        gain = rho - d
        x = 0.5 * (sqrt(d^2 + 4 * lam * gain) - d)            # x (x + d) = lam (rho - d)
        b = sqrt(max(rho^2 - d^2, 1.0))
        t_on = n * x / (2 * a_OCL(p)) * b / gain              # duration of the link
        if t + t_on >= T
            d += x * (T - t) / t_on
            return (dvT = n * d / 2, da_km = d / 1.0e3, links = k)
        end
        t += t_on; d += x; k += 1
        t_gap = 2 * (s - p.L_max + x^2 / (2 * lam)) / (3 * n * d)   # coast until the next helper is within range
        t + t_gap >= T && return (dvT = n * d / 2, da_km = d / 1.0e3, links = k)
        t += t_gap; rho = p.L_max
    end
end

# Resolve the helper phase instead of averaging over one helper spacing as in averages().
function leapfrog_rtn(N_h, b0, inc_deg, u0_deg, n_orbits, dt, p::Params)
    isinteger(N_h) && N_h >= 1 || throw(ArgumentError("expected a positive helper count"))
    rh, s, acc = r_h(p), spacing(N_h, p), a_OCL(p)
    s > p.L_max || throw(ArgumentError("leapfrog requires sparse helpers; use predict for a dense ring"))
    0 < b0 <= s || throw(ArgumentError("b0 must be the positive arc to the nearest trailing helper"))
    isfinite(dt) && dt > 0 && isfinite(n_orbits) && n_orbits >= 0 &&
        isfinite(inc_deg) && isfinite(u0_deg) || throw(ArgumentError("invalid time or orbital angle"))
    n0 = mean_motion(rh, p)
    duration = n_orbits * 2pi / n0
    sini, u0 = sind(inc_deg), deg2rad(u0_deg)

    function rates(t, d, lag)
        n = mean_motion(rh + d, p)
        # lag is the unwrapped arc ahead of the *initial* helper. Wrapping picks
        # the current trailing helper. With s > L_max at most one can contribute.
        b = mod(lag, s)
        drift = rh * (n - n0)
        u = u0 + n0*t + (lag - b0)/rh
        z = (rh + d) * sini * sin(u)
        # Same local RTN sight-line approximation (d, b, z) as the averaged model.
        rho = sqrt(d*d + b*b + z*z)
        if b == 0 || rho > p.L_max
            return (; R=0.0, T=0.0, N=0.0, rise=0.0, drift, helper=nothing)
        end
        g = acc / rho
        (; R=g*d, T=g*b, N=g*z, rise=2*g*b/n, drift, helper=floor(Int, lag/s))
    end

    d = dvR = dvT = dvNp = dvNm = 0.0
    lag = float(b0)
    links = 0
    if b0 <= p.L_max && duration > 0
        steps = ceil(Int, duration/dt)
        h = duration/steps
        previous = nothing
        for j in 0:(steps-1)
            t = j*h
            start = rates(t, d, lag)
            mid = rates(t + h/2, d + h/2*start.rise, lag + h/2*start.drift)
            links += !isnothing(previous) && mid.helper != previous
            previous = mid.helper
            dvR += h*mid.R
            dvT += h*mid.T
            dvNp += h*max(mid.N, 0.0)
            dvNm += h*min(mid.N, 0.0)
            d += h*mid.rise
            lag += h*mid.drift
        end
        links += !isnothing(previous) && rates(duration, d, lag).helper != previous
    end
    (; dvR, dvT, dvNp, dvNm, da_km=d/1e3, links)
end

"""
    inc_change_fixed_offset(y, inc_deg; n_orbits = 100, p = Params())

Delta i-max rule, target on the ring, nearest helper at the fixed along-track offset y (m).
Returns the change of inclination in deg. For a target that starts at the offset y0 the simulation lies
between this value at y0 and the maximum a_OCL T / (pi n r_t).
"""
function inc_change_fixed_offset(y, inc_deg; n_orbits = 100, p::Params = Params())
    n = mean_motion(r_h(p), p)
    T = n_orbits * 2pi / n
    A = r_h(p) * sin(deg2rad(inc_deg))
    return rad2deg(a_OCL(p) * T / (pi * n * r_h(p)) * (hypot(A, y) - abs(y)) / A)
end

# ------------------------------------------------------------------------------------------------
# Run this file directly (julia oracle_closed_forms.jl) to print a few cases as a check
# ------------------------------------------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    p = Params()
    println("a_OCL = ", round(a_OCL(p), sigdigits = 4), " m/s^2")
    println("rule    N_h   h_t   inc |    dvR    dvT   dvN+   dvN- |  da_km  di_mdeg  dOmega_deg")
    for (rule, N_h, h_t, inc) in ((:amax, 100, 1050.0, 0.0), (:amax, 800, 1050.0, 0.0), (:amax, 300, 1100.0, 0.0),
                                  (:range, 300, 1050.0, 0.0), (:amax, 800, 1050.0, 1.0), (:imax, 800, 1050.0, 1.0),
                                  (:imax, 800, 1100.0, 0.5), (:range, 800, 1050.0, 0.5))
        r = predict(rule, N_h, h_t; inc_deg = inc, p = p)
        println(rpad(string(rule), 6), lpad(N_h, 5), lpad(h_t, 7), lpad(inc, 5), " |",
                lpad(round(r.dvR, digits = 2), 7), lpad(round(r.dvT, digits = 2), 7),
                lpad(round(r.dvNp, digits = 2), 7), lpad(round(r.dvNm, digits = 2), 7), " |",
                lpad(round(r.da_km, digits = 2), 7), lpad(round(1.0e3 * r.di_deg, digits = 2), 9),
                lpad(round(r.dOmega_deg, digits = 3), 12))
    end
    # expected output of the first two lines (default parameters):
    #   amax    100 1050.0  0.0 |   2.14   2.97    0.0    0.0 |   6.02      0.0         0.0
    #   amax    800 1050.0  0.0 |   3.22   8.78    0.0    0.0 |  17.84      0.0         0.0
end
