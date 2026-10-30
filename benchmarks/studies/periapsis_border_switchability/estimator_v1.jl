# Estimator v1 (task A2): per-pass mean-element propagator for the Earth heat-rate-band scenario.
#
# State: mean elements (a, k, h, i, Omega) with (k, h) = (e cos w, e sin w), advanced orbit by orbit from apoapsis to apoapsis.
# Physics, no fitted coefficient:
#   (a) single-averaged exact point-mass Sun and Moon tide: Gauss variational equations averaged over the Kepler ellipse
#       (uniform eccentric-anomaly nodes, weight dM = (1 - e cos E) dE), bodies held at their SPICE positions at the pass epoch;
#   (b) secular J2 rates of Omega, w and the J2 mean-motion rate (Brouwer/Kozai; Vallado, Fundamentals of Astrodynamics and
#       Applications, secular J2 section); J2 first-order short-period map (Vallado-Kwok, the corridor example's
#       PeriapsisCorridorMeanElements) between mean and osculating elements;
#   (c) drag: the Gauss equations integrated over the Kepler ellipse of the OSCULATING periapsis state with the real drag
#       acceleration -0.5 rho |v_rel| v_rel / beta (rho = the piecewise US76 fit at geodetic altitude, v_rel = v - w_E x r,
#       zero wind), Gauss-Legendre in eccentric anomaly over the part of the orbit within `drag_window_m` of periapsis;
#   (d) J3 and the Earth orientation (precession, polar motion) omitted; geodetic altitude from the repo's Bowring routine.
# Usage: included by run_selection_skill.jl. No main().
isdefined(@__MODULE__, :run_truth) || include(joinpath(@__DIR__, "common.jl"))   # guard: a driver that already loaded common.jl may include this file
using LinearAlgebra, StaticArrays, CSV, DataFrames, TOML

const SV5 = SVector{5, Float64}
const SV3 = SVector{3, Float64}

Base.@kwdef struct EstSettings
    n_tide::Int = 64
    n_drag_half::Int = 32
    drag_window_m::Float64 = 250.0e3
    body_epoch::Symbol = :periapsis       # :periapsis (window apo-apo centered on the pass) or :apoapsis (half an orbit later)
    tides::Bool = true
    j2::Bool = true
    drag::Bool = true
    horizon_s::Float64 = 60 * 86400.0
    margin_m::Float64 = 1000.0
    impact_altitude_m::Float64 = 50.0e3
end

with_settings(s::EstSettings; kw...) = EstSettings(; merge(NamedTuple{fieldnames(EstSettings)}(Tuple(getfield(s, f) for f in fieldnames(EstSettings))), NamedTuple(kw))...)

struct EstContext
    mu::Float64
    Re::Float64            # gravity reference radius
    J2::Float64
    Rp_e::Float64          # planet equatorial radius (impact test, initial state)
    planet::Any
    omega_e::SV3
    fit::Any               # CorridorAtmosphereFit
    beta::Float64
    q_low::Float64
    q_high::Float64
    q_hard::Float64
    w2m::Float64
    apo_alt_m::Float64
    body_names::Vector{String}
    body_mu::Vector{Float64}
    gl_x::Vector{Float64}
    gl_w::Vector{Float64}
    center::String         # SPICE name of the primary the third-body positions are taken relative to ("EARTH", "MARS")
    rot::Any               # J2000 -> body-equator inertial rotation (3x3), or `nothing` when the pole is the J2000 z axis (Earth)
end

# Legacy 17-argument form (Earth, J2000 z = pole): every pre-existing caller keeps working unchanged.
EstContext(mu, Re, J2, Rp_e, planet, omega_e, fit, beta, q_low, q_high, q_hard, w2m, apo_alt_m, body_names, body_mu, gl_x, gl_w) =
    EstContext(mu, Re, J2, Rp_e, planet, omega_e, fit, beta, q_low, q_high, q_hard, w2m, apo_alt_m, body_names, body_mu, gl_x, gl_w, "EARTH", nothing)

"""Gauss-Legendre nodes and weights on [0, 1] (Golub-Welsch)."""
function gauss_legendre01(n::Int)
    b = [i / sqrt(4i^2 - 1) for i in 1:n-1]
    ev = eigen(SymTridiagonal(zeros(n), b))
    return (ev.values .+ 1) ./ 2, ev.vectors[1, :] .^ 2
end

"""
Planet record (stage C, planet generality): every planet-dependent constant the estimator uses, as an input. `harm` is the truth gravity
harmonics model (GM, reference radius and J2 are read from it, as at Earth), `planet` the engine planet struct (equatorial radius, spin,
geodetic routine), `fit` the density fit (`CorridorAtmosphereFit`, single or piecewise layers), `rot` the J2000 -> body-equator inertial
rotation (`nothing` = the pole is the J2000 z axis, Earth as before).
"""
Base.@kwdef struct PlanetRecord
    name::String
    center::String
    planet::Any
    harm::Any
    fit::Any
    beta::Float64
    q_low::Float64
    q_high::Float64
    q_hard::Float64
    w2m::Float64
    apo_alt_m::Float64
    body_names::Vector{String}
    body_mu::Vector{Float64}
    rot::Any = nothing
end

"""The Earth record: exactly the constants `make_context` used before the parametrization."""
function earth_record(study)
    cfg = load_corridor_config(CORRIDOR_CONFIG_PATH)
    planet = make_no_gram_planet(:earth)
    harm = GravitationalHarmonicsModel(3, 0, joinpath(REPO_ROOT, String(cfg.values["gravity"]["coefficient_file"])), planet)
    b = study["band"]
    return PlanetRecord(name="Earth", center="EARTH", planet=planet, harm=harm, fit=fitted_atmosphere_parameters(cfg),
        beta=study_beta(study, cfg.values).beta, q_low=Float64(b["q_low_w_cm2"]), q_high=Float64(b["q_high_w_cm2"]), q_hard=Float64(b["q_hard_w_cm2"]),
        w2m=Float64(b["w_cm2_to_w_m2"]), apo_alt_m=Float64(cfg.values["orbit"]["apoapsis_altitude_m"]),
        body_names=String.(study["third_bodies"]["bodies"]), body_mu=[MU_BODY[String(b_)] for b_ in study["third_bodies"]["bodies"]], rot=nothing)
end

function make_context(rec::PlanetRecord, s::EstSettings)
    gx, gw = gauss_legendre01(s.n_drag_half)
    harm = rec.harm; planet = rec.planet
    return EstContext(harm.gm_m3s2, harm.reference_radius_m, -sqrt(5.0) * harm.C[3, 1], planet.Rp_e, planet,
        SV3(planet.ω), rec.fit, rec.beta, rec.q_low, rec.q_high, rec.q_hard, rec.w2m, rec.apo_alt_m,
        rec.body_names, rec.body_mu, gx, gw, rec.center, rec.rot)
end

make_context(study, s::EstSettings) = make_context(earth_record(study), s)

# ---- atmosphere / band ------------------------------------------------------------------------

"""Geodetic altitude (m) where the piecewise-exponential fit has density `rho` (closed-form inverse of the fit; the fit extends its last layer)."""
function altitude_of_density_fast(fit, rho::Float64)
    nl = length(fit.scale_heights_m); k = 1
    while k < nl && rho < fit.node_density_kg_m3[k + 1]; k += 1; end
    return fit.node_altitudes_m[k] - fit.scale_heights_m[k] * log(rho / fit.node_density_kg_m3[k])
end

"""`(h_at_q_high, h_at_q_low)` altitude interval of the band at relative speed `vrel` (same definition as `band_altitudes`)."""
band_edges(ctx::EstContext, vrel) = (altitude_of_density_fast(ctx.fit, 2 * ctx.q_high * ctx.w2m / vrel^3),
                                     altitude_of_density_fast(ctx.fit, 2 * ctx.q_low * ctx.w2m / vrel^3))

# ---- orbit geometry and the nonsingular Gauss equations ---------------------------------------

struct OrbitFrame
    a::Float64; e::Float64; i::Float64; w::Float64; cw::Float64; sw::Float64; sini::Float64; cosi::Float64
    P::SV3; Q::SV3; W::SV3; p::Float64; hm::Float64; eta::Float64; n::Float64
end

function OrbitFrame(a, k, h, i, Om, mu)
    e = hypot(k, h); w = atan(h, k); si, ci = sincos(i); sO, cO = sincos(Om)
    W = SV3(si * sO, -si * cO, ci); node = SV3(cO, sO, 0.0)
    cw, sw = e > 0 ? (k / e, h / e) : (1.0, 0.0)
    P = cw * node + sw * cross(W, node); Q = cross(W, P)
    p = a * (1 - e^2)
    return OrbitFrame(a, e, i, w, cw, sw, si, ci, P, Q, W, p, sqrt(mu * p), sqrt(1 - e^2), sqrt(mu / a^3))
end

"""Position, velocity, radius, cos/sin of true anomaly at eccentric anomaly `E` on the ellipse of `F`."""
@inline function node_state(F::OrbitFrame, E)
    sE, cE = sincos(E); d = 1 - F.e * cE
    r = F.a * (cE - F.e) * F.P + F.a * F.eta * sE * F.Q
    v = (F.n * F.a / d) * (-sE * F.P + F.eta * cE * F.Q)
    return r, v, F.a * d, (cE - F.e) / d, F.eta * sE / d, d
end

"""Gauss variational equations in (a, k, h, i, Omega) for a perturbing acceleration with radial/transverse/normal components (fR, fS, fW)."""
@inline function gauss_g(F::OrbitFrame, rm, cnu, snu, fR, fS, fW)
    p, hm, e = F.p, F.hm, F.e
    su = F.sw * cnu + F.cw * snu; cu = F.cw * cnu - F.sw * snu
    da = 2 * F.a^2 / hm * (e * snu * fR + p / rm * fS)
    de = (p * snu * fR + ((p + rm) * cnu + rm * e) * fS) / hm
    ew = (-p * cnu * fR + (p + rm) * snu * fS) / hm - e * rm * su * F.cosi / (hm * F.sini) * fW   # e * d(omega)/dt
    di = rm * cu / hm * fW
    dO = rm * su / (hm * F.sini) * fW
    return SV5(da, F.cw * de - F.sw * ew, F.sw * de + F.cw * ew, di, dO)
end

@inline tide_acc(r::SV3, rb::SV3, mub) = (d = rb - r; mub * (d / dot(d, d)^1.5 - rb / dot(rb, rb)^1.5))

"""Orbit-averaged (per second, uniform in mean anomaly) rates of (a, k, h, i, Omega) from point-mass bodies `rbs`, `mubs`."""
function tide_rates(F::OrbitFrame, rbs, mubs, N::Int)
    acc = zero(SV5)
    for j in 1:N
        E = 2π * (j - 0.5) / N
        r, _, rm, cnu, snu, d = node_state(F, E)
        f = zero(SV3)
        for b in eachindex(rbs); f += tide_acc(r, rbs[b], mubs[b]); end
        rh = r / rm; th = cross(F.W, rh)
        acc += d * gauss_g(F, rm, cnu, snu, dot(f, rh), dot(f, th), dot(f, F.W))
    end
    return acc / N
end

"""Secular J2 rates (dOmega, domega, dM) of the mean state (Vallado, secular J2 section)."""
function j2_rates(F::OrbitFrame, mu, J2, Re)
    f = F.n * J2 * (Re / F.p)^2
    return -1.5 * f * F.cosi, 0.75 * f * (5 * F.cosi^2 - 1), F.n * (1 + 0.75 * J2 * (Re / F.p)^2 * F.eta * (2 - 3 * F.sini^2))
end

"""Advance the mean state by `dt`: additive tide/drag increment `inc` (a, k, h, i, Omega per second), then the J2 apsidal rotation of (k, h)
by `dw * dt` (exact rotation: an additive step on a rotating vector would grow e by (dw dt)^2 / 2 per step) and `dOm * dt` on Omega."""
@inline function advance(m::SV5, inc::SV5, dOm, dw, dt)
    x = m + inc * dt; c, sn = cos(dw * dt), sin(dw * dt)
    return SV5(x[1], c * x[2] - sn * x[3], sn * x[2] + c * x[3], x[4], x[5] + dOm * dt)
end

mean_motion_rate(m::SV5, ctx::EstContext) = j2_rates(OrbitFrame(m[1], m[2], m[3], m[4], m[5], ctx.mu), ctx.mu, ctx.J2, ctx.Re)[3]

# ---- mean <-> osculating (first-order J2 short-period map, Vallado-Kwok) ----------------------

mean_state(el, ::EstContext) = SV5(el.a, el.e * cos(el.omega), el.e * sin(el.omega), el.i, el.raan)

"""Osculating elements at mean true anomaly 0 (the radial-event periapsis) of the mean state `m`."""
function osc_at_periapsis(m::SV5, ctx::EstContext)
    e = hypot(m[2], m[3]); w = mod(atan(m[3], m[2]), 2π)
    return mean_to_osculating((a=m[1], e=e, i=m[4], raan=mod(m[5], 2π), omega=w, nu=0.0); mu=ctx.mu, radius=ctx.Re, j2=ctx.J2).elements
end

"""Osculating periapsis state of `m` as the (a, k, h, i, Omega) vector (same coordinates as the mean state)."""
function osc_vec(m::SV5, ctx::EstContext)
    o = osc_at_periapsis(m, ctx)
    return SV5(o.a, o.e * cos(o.omega), o.e * sin(o.omega), o.i, o.raan)
end

wrap_pi(x) = mod(x + π, 2π) - π

"""
Increment `d` of (a, k, h, i, Omega) from the Gauss rates (computed on the ellipse `v`), re-expressed as increments of (a, r_p, omega, i, Omega):
dr_p = (1 - e) da - a de. Over a pass dr_p nearly cancels (drag at periapsis leaves r_p almost unchanged) while the product da de of the
two large increments does not (3 m against a 0.6 m first-order change), so the increments are applied to a and r_p and e is recovered
from them; the conversion uses the ellipse the rates were computed on, for the inbound and the outbound part alike.
"""
function rp_form(v::SV5, d::SV5)
    a, k, h = v[1], v[2], v[3]; e = hypot(k, h)
    de = (k * d[2] + h * d[3]) / e
    return SV5(d[1], (1 - e) * d[1] - a * de, (k * d[3] - h * d[2]) / e^2, d[4], d[5])
end

function apply_form(v::SV5, f::SV5)
    a, k, h = v[1], v[2], v[3]; e = hypot(k, h)
    a2 = a + f[1]; e2 = 1 - (a * (1 - e) + f[2]) / a2; w2 = atan(h, k) + f[3]
    return SV5(a2, e2 * cos(w2), e2 * sin(w2), v[4] + f[4], v[5] + f[5])
end

"""
Mean state after the increment `f` (an `rp_form`) of the OSCULATING periapsis elements (drag acts on the osculating orbit): solves
map(m') = apply_form(map(m), f) for the mean state m' by fixed-point iteration (the map's Jacobian differs from the identity by O(10^-2),
so three iterations leave an error of order 10^-6 of f). Adding the osculating increment to the mean state directly would err by the
Jacobian: tens of meters of periapsis radius per pass in the drag-only test.
"""
function shift_mean(m::SV5, f::SV5, ctx::EstContext)
    base = osc_vec(m, ctx); target = apply_form(base, f); x = m + (target - base)
    for _ in 1:3
        r = target - osc_vec(x, ctx)
        x += SV5(r[1], r[2], r[3], r[4], wrap_pi(r[5]))
    end
    return x
end

# ---- drag -------------------------------------------------------------------------------------

"""
Changes of (a, k, h, i, Omega) accumulated on the inbound (E < 0) and outbound (E > 0) parts of one pass through the
atmosphere, from the Gauss equations integrated over time (dt = (1 - e cos E) dE / n) on the Kepler ellipse of the osculating
periapsis state `osc`, over |E| <= Emax (altitude above periapsis radius <= window). `rho_fn(r_vec)` returns density at the node.
"""
function drag_deltas(osc, ctx::EstContext, s::EstSettings, rho_fn)
    F = OrbitFrame(osc.a, osc.e * cos(osc.omega), osc.e * sin(osc.omega), osc.i, osc.raan, ctx.mu)
    c = 1 - s.drag_window_m / (F.a * F.e)
    Emax = c <= -1 ? float(π) : acos(c)
    din = zero(SV5); dout = zero(SV5)
    for j in eachindex(ctx.gl_x), sgn in (-1.0, 1.0)
        E = sgn * Emax * ctx.gl_x[j]; wj = Emax * ctx.gl_w[j]
        r, v, rm, cnu, snu, d = node_state(F, E)
        rho = rho_fn(r)
        vrel = v - cross(ctx.omega_e, r); V = norm(vrel)
        acc = (-0.5 * rho * V / ctx.beta) * vrel
        rh = r / rm; th = cross(F.W, rh)
        g = (wj * d / F.n) * gauss_g(F, rm, cnu, snu, dot(acc, rh), dot(acc, th), dot(acc, F.W))
        if sgn < 0; din += g; else; dout += g; end
    end
    return din, dout
end

fit_rho(ctx::EstContext) = r -> corridor_density(ctx.fit, geodetic_altitude_latitude(r, ctx.planet)[1])

# ---- prediction -------------------------------------------------------------------------------

mutable struct Prediction
    t_s::Vector{Float64}; alt_m::Vector{Float64}; radius_m::Vector{Float64}; lat_rad::Vector{Float64}; vrel_mps::Vector{Float64}
    rho::Vector{Float64}; q_w_cm2::Vector{Float64}; a_m::Vector{Float64}; e::Vector{Float64}; i_rad::Vector{Float64}
    raan_rad::Vector{Float64}; omega_rad::Vector{Float64}; h_qhigh_m::Vector{Float64}; h_qlow_m::Vector{Float64}
    decayed::Bool; t_decay_s::Float64; ended::Symbol    # :horizon | :impact | :invalid
end
Prediction() = Prediction([Float64[] for _ in 1:14]..., false, NaN, :horizon)

body_positions(ctx::EstContext, et) = ctx.rot === nothing ?
    [SV3(spkpos(uppercase(b), et, "J2000", "NONE", ctx.center)[1]) .* 1e3 for b in ctx.body_names] :
    [SV3(ctx.rot' * (SV3(spkpos(uppercase(b), et, "J2000", "NONE", ctx.center)[1]) .* 1e3)) for b in ctx.body_names]

"""Period-averaged rates of the mean state at ephemeris time `et_body`: tide increment (a, k, h, i, Omega)/s and the J2 `(dOmega, domega)`."""
function mean_rates(m::SV5, ctx::EstContext, s::EstSettings, et_body)
    F = OrbitFrame(m[1], m[2], m[3], m[4], m[5], ctx.mu)
    dOm, dw, _ = s.j2 ? j2_rates(F, ctx.mu, ctx.J2, ctx.Re) : (0.0, 0.0, 0.0)
    inc = s.tides ? tide_rates(F, body_positions(ctx, et_body), ctx.body_mu, s.n_tide) : zero(SV5)
    return inc, dOm, dw
end

"""Initial osculating Keplerian state at apoapsis (true anomaly 180 deg), as the corridor example builds it."""
function initial_osculating(ctx::EstContext, seed_alt_m, inc_rad, raan_rad, omega_rad)
    rp = ctx.Rp_e + seed_alt_m; ra = ctx.Rp_e + ctx.apo_alt_m
    return (a=0.5 * (ra + rp), e=(ra - rp) / (ra + rp), i=inc_rad, raan=raan_rad, omega=omega_rad, nu=float(π))
end

true_to_mean_anomaly(e, f) = (E = 2atan(sqrt(1 - e) * sin(f / 2), sqrt(1 + e) * cos(f / 2)); mod(E - e * sin(E), 2π))

"""
Predict the periapsis history of one orbit from its seed (initial periapsis altitude) and geometry. `et0`: ephemeris time of the
run start. `max_passes` limits the number of passes (calibration uses 1). Rows are recorded at each radial-event periapsis.
"""
function predict_history(ctx::EstContext, s::EstSettings, seed_alt_m, inc_rad, raan_rad, omega_rad, et0; max_passes::Int=typemax(Int), rho_fn=fit_rho(ctx))
    out = Prediction()
    init = initial_osculating(ctx, seed_alt_m, inc_rad, raan_rad, omega_rad)
    mm = osculating_to_mean(init; mu=ctx.mu, radius=ctx.Re, j2=ctx.J2).elements
    m = mean_state(mm, ctx); M0 = true_to_mean_anomaly(mm.e, mm.nu)
    t_apo = 0.0; mdot = mean_motion_rate(m, ctx)
    for kpass in 1:max_passes
        mdot = mean_motion_rate(m, ctx)
        T = 2π / mdot
        dt_in = kpass == 1 ? (2π - M0) / mdot : 0.5 * T
        dt_out = 0.5 * T
        t_k = t_apo + dt_in
        t_k > s.horizon_s && break
        et_body = et0 + t_k + (s.body_epoch === :apoapsis ? dt_out : 0.0)
        i1, dO1, dw1 = mean_rates(m, ctx, s, et_body)
        i2, dO2, dw2 = mean_rates(advance(m, i1, dO1, dw1, dt_in), ctx, s, et_body)
        m_pre = advance(m, i2, dO2, dw2, dt_in)
        osc = try osc_at_periapsis(m_pre, ctx) catch; out.decayed = true; out.t_decay_s = t_k; out.ended = :invalid; break end
        din = zero(SV5); dout = zero(SV5)
        if s.drag; din, dout = drag_deltas(osc, ctx, s, rho_fn); end
        ov = osc_vec(m_pre, ctx); fin = rp_form(ov, din); fout = rp_form(ov, dout)
        m_peri = try shift_mean(m_pre, fin, ctx) catch; out.decayed = true; out.t_decay_s = t_k; out.ended = :invalid; break end
        osc = try osc_at_periapsis(m_peri, ctx) catch; out.decayed = true; out.t_decay_s = t_k; out.ended = :invalid; break end
        rp = osc.a * (1 - osc.e)
        if rp - ctx.Rp_e <= s.impact_altitude_m; out.decayed = true; out.t_decay_s = t_k; out.ended = :impact; break; end
        rv = elements_to_cartesian((a=osc.a, e=osc.e, i=osc.i, raan=osc.raan, omega=osc.omega, nu=0.0); mu=ctx.mu)
        r = SV3(rv[1]); v = SV3(rv[2])
        alt, lat = geodetic_altitude_latitude(r, ctx.planet)
        vrel = norm(v - cross(ctx.omega_e, r)); rho = corridor_density(ctx.fit, alt)
        hq, hl = band_edges(ctx, vrel)
        for (arr, val) in zip((out.t_s, out.alt_m, out.radius_m, out.lat_rad, out.vrel_mps, out.rho, out.q_w_cm2, out.a_m, out.e, out.i_rad, out.raan_rad, out.omega_rad, out.h_qhigh_m, out.h_qlow_m),
                              (t_k, alt, rp, lat, vrel, rho, 0.5 * rho * vrel^3 / ctx.w2m, osc.a, osc.e, osc.i, osc.raan, osc.omega, hq, hl))
            push!(arr, val)
        end
        m = try advance(shift_mean(m_peri, fout, ctx), i2, dO2, dw2, dt_out) catch; out.decayed = true; out.t_decay_s = t_k; out.ended = :invalid; break end
        t_apo = t_k + dt_out
    end
    return out
end

"""Seed (initial periapsis altitude) that makes the predicted first-pass peak Q equal the band's geometric mean (the truth runs' calibration rule), iterated as the truth calibration does."""
function calibrate_seed(ctx::EstContext, s::EstSettings, inc_rad, raan_rad, omega_rad, et0; seed0=140.0e3, tol_m=1.0, max_iter=8)
    q_target = sqrt(ctx.q_low * ctx.q_high); seed = seed0
    for _ in 1:max_iter
        p = predict_history(ctx, s, seed, inc_rad, raan_rad, omega_rad, et0; max_passes=1)
        isempty(p.t_s) && return NaN
        h_target = altitude_of_density_fast(ctx.fit, 2 * q_target * ctx.w2m / p.vrel_mps[1]^3)
        abs(h_target - p.alt_m[1]) < tol_m && return seed
        seed += h_target - p.alt_m[1]
    end
    return seed
end

# ---- double-averaged hand model (baseline B1), predictive ------------------------------------------

"""
Baseline B1: the earlier studies' double-averaged Kozai/Lidov formula, dr_p/dt = -a (15/8) (mu_b/(n r_b^3)) e sqrt(1-e^2) sin(2 w') sin^2(i'),
summed over Sun and Moon, with a, e, i held at the initial osculating values and w, Omega advanced by the secular J2 rates (no drag, no mean
elements); the perturber planes come from r_b x v_b at the pass epoch (`perturber_angles`). Altitude history = first-pass altitude (taken
from `alt1`, `vrel1`: the estimator's calibrated first pass) + the integrated r_p change; band and speed held at the first pass.
"""
function hand_history(ctx::EstContext, s::EstSettings, seed_alt_m, inc_rad, raan_rad, omega_rad, et0, alt1, vrel1)
    init = initial_osculating(ctx, seed_alt_m, inc_rad, raan_rad, omega_rad)
    a, e = init.a, init.e; n = sqrt(ctx.mu / a^3); T = 2π / n
    F = OrbitFrame(a, e * cos(omega_rad), e * sin(omega_rad), inc_rad, raan_rad, ctx.mu)
    f = n * ctx.J2 * (ctx.Re / F.p)^2
    dOm = -1.5 * f * F.cosi; dw = 0.75 * f * (5 * F.cosi^2 - 1)
    t = T / 2; drp = 0.0; tprev = 0.0
    ts = Float64[]; alts = Float64[]
    hq, hl = band_edges(ctx, vrel1)
    while t <= s.horizon_s
        rate = 0.0
        for (bi, b) in enumerate(ctx.body_names)
            pa = perturber_angles(b, et0 + t, inc_rad, raan_rad + dOm * t, omega_rad + dw * t; center=ctx.center, rot=ctx.rot)
            rate += -a * (15 / 8) * (pa.mu / pa.r_m^3 / n) * e * sqrt(1 - e^2) * sin(2pa.omega_p) * sin(pa.i_p)^2
        end
        drp += rate * (t - tprev)
        push!(ts, t); push!(alts, alt1 + drp)
        tprev = t; t += T
    end
    return (t_s=ts, alt_m=alts, q_w_cm2=[0.5 * corridor_density(ctx.fit, h) * vrel1^3 / ctx.w2m for h in alts], a_m=fill(a, length(ts)), e=fill(e, length(ts)),
            h_qhigh_m=fill(hq, length(ts)), h_qlow_m=fill(hl, length(ts)), vrel_mps=fill(vrel1, length(ts)))
end

# ---- frames of reference for the run -------------------------------------------------------------

epoch_et(study, idx::Int) = SimulationModel.ephemerides_time_seconds(parse_epoch(study["grid"]["epochs_utc"][idx]), SpiceEphemeridesModel())

function settings_from(toml::Dict)
    e = toml["estimator"]
    return EstSettings(n_tide=Int(e["n_tide_nodes"]), n_drag_half=Int(e["n_drag_nodes_per_half"]), drag_window_m=1000 * Float64(e["drag_window_km"]),
        body_epoch=Symbol(e["body_epoch"]), horizon_s=SECONDS_PER_DAY * Float64(e["horizon_days"]), margin_m=1000 * Float64(e["margin_km"]),
        impact_altitude_m=Float64(e["impact_altitude_m"]))
end
