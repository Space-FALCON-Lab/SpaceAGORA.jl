# Estimator v2: estimator v1 plus first-order SHORT-PERIOD (orbital-frequency) terms of the third-body tide, so that
#   (i)  the initial MEAN elements are the osculating elements minus the J2 and the Sun/Moon short-period terms, and
#   (ii) the osculating periapsis state of every pass is  mean + J2 short period + Sun/Moon short period.
# The long-period (apse-line, tens of days) behavior is the same single-averaged tide as v1. Derivation: see
# tasks/.../estimator_v2_results.md and shortperiod_derivation.py (closed forms, sympy).
#
# Third-body short period (derived, no fitted coefficient). With the perturber frozen at its position over one orbit (its motion
# enters at relative order nu = n_b/n, about 0.03 for the Moon), the first-order osculating-minus-mean difference of any element X is
#     dX(E) = (1/n) * int_{E*}^{E} (g(E') - <g>) (1 - e cos E') dE'  + const,     <dX>_M = 0,
# g = Gauss rate of the exact point-mass (direct + indirect) tide on the Kepler ellipse. The integrand is smooth and periodic in the
# eccentric anomaly E; with N uniform nodes (the same nodes as the averaged tide) the antiderivative is a Fourier sum and
#     dX(E=0)  = (1/n) [ - sum_k B_k / k  - e B_1 / 2 ],      dX(E=pi) = (1/n) [ - sum_k (-1)^k B_k / k - e B_1 / 2 ],
#     B_k = (2/N) sum_j (1 - e cos E_j) g_j sin(k E_j)                                         (the <g> term drops out of B_k).
# In the quadrupole limit this is Pi_3B times a quadratic form of the perturber direction cosines in the orbit frame (closed form
# `quad_sp_closed`, derived with sympy in shortperiod_derivation.py and verified against `tide_sp`).
# Usage: included by run_estimator_v2.jl. No main().
include(joinpath(@__DIR__, "estimator_v1.jl"))

# ---- options --------------------------------------------------------------------------------------------

"""
`sp_tide`: `:exact` (quadrature of the exact point-mass tide, all multipoles), `:quadrupole` (closed form Pi_3B form for (a, e) only), `:none` (= v1).
`sp_j2`: `:map` (the repo's Vallado-Kwok first-order map, as v1) or `:quad` (the same first-order quantity by quadrature of the J2 Gauss rates).
`n_sp`: number of eccentric-anomaly nodes of the short-period quadrature.
"""
Base.@kwdef struct V2Opts
    sp_tide::Symbol = :exact
    sp_j2::Symbol = :map
    n_sp::Int = 64
end

# ---- short-period quadrature --------------------------------------------------------------------------------

struct SPTable
    N::Int
    K::Int
    E::Vector{Float64}
    S::Matrix{Float64}      # S[j, k] = sin(k E_j)
    cosE::Vector{Float64}
    sinE::Vector{Float64}
end

function SPTable(N::Int)
    K = N ÷ 2 - 1
    E = [2π * (j - 0.5) / N for j in 1:N]
    return SPTable(N, K, E, [sin(k * Ej) for Ej in E, k in 1:K], cos.(E), sin.(E))
end

"""Node weights phi_j such that dX(E*) = (2/(N n)) sum_j (1 - e cos E_j) g_j phi_j  (`which` = :peri for E* = 0, :apo for E* = pi)."""
function sp_weights(tab::SPTable, e, which::Symbol)
    phi = zeros(tab.N)
    for k in 1:tab.K
        w = which === :peri ? -1.0 / k : -(isodd(k) ? -1.0 : 1.0) / k
        for j in 1:tab.N; phi[j] += w * tab.S[j, k]; end
    end
    for j in 1:tab.N; phi[j] -= 0.5 * e * tab.S[j, 1]; end
    return phi
end

"""Short-period difference (osculating minus mean) of (a, k, h, i, Omega) at periapsis (`:peri`) or apoapsis (`:apo`) of the mean orbit `F`
for the acceleration `accfun(r)` (frozen perturbers / J2), first order."""
function sp_quadrature(F::OrbitFrame, accfun, tab::SPTable, which::Symbol)
    phi = sp_weights(tab, F.e, which)
    acc = zero(SV5)
    for j in 1:tab.N
        E = tab.E[j]
        r, _, rm, cnu, snu, d = node_state(F, E)
        f = accfun(r)
        rh = r / rm; th = cross(F.W, rh)
        acc += (d * phi[j]) * gauss_g(F, rm, cnu, snu, dot(f, rh), dot(f, th), dot(f, F.W))
    end
    return acc * (2 / (tab.N * F.n))
end

"""Same for a general eccentric anomaly `Ev` (verification only: needs the cosine coefficients too). Returns dX(Ev)."""
function sp_quadrature_general(F::OrbitFrame, accfun, tab::SPTable, Ev)
    N, K = tab.N, tab.K
    A = [zero(SV5) for _ in 1:K]; B = [zero(SV5) for _ in 1:K]; gbar = zero(SV5)
    for j in 1:N
        E = tab.E[j]
        r, _, rm, cnu, snu, d = node_state(F, E)
        f = accfun(r); rh = r / rm; th = cross(F.W, rh)
        g = d * gauss_g(F, rm, cnu, snu, dot(f, rh), dot(f, th), dot(f, F.W))
        gbar += g / N                                       # <g>_M
        for k in 1:K
            A[k] += (2 / N) * cos(k * E) * g; B[k] += (2 / N) * sin(k * E) * g
        end
    end
    A[1] += F.e * gbar                  # the cos(E) coefficient of <g>(1 - e cos E), which the series of (g - <g>) d must not contain
    H = zero(SV5)
    for k in 1:K; H += (A[k] * sin(k * Ev) - B[k] * cos(k * Ev)) / k; end
    return (H - 0.5 * F.e * B[1]) / F.n
end

@inline j2_acc(r::SV3, J2, mu, Re) = (rn2 = dot(r, r); rn = sqrt(rn2); z2 = r[3]^2 / rn2; k = -1.5 * J2 * mu * Re^2 / rn2^2;
    k * SV3((1 - 5z2) * r[1] / rn, (1 - 5z2) * r[2] / rn, (3 - 5z2) * r[3] / rn))

tide_accfun(rbs, mubs) = r -> (f = zero(SV3); for b in eachindex(rbs); f += tide_acc(r, rbs[b], mubs[b]); end; f)

# ---- closed-form quadrupole short period of (a, e) (derived: shortperiod_derivation.py) ---------------------------------------

"""
Quadrupole-order short-period differences (osculating minus mean) `(da/a, de)` per unit Pi_3B = (mu_b/mu)(a/r_b)^3 at periapsis/apoapsis,
as functions of e and the perturber direction cosines (bx, by) along the apse line and the semilatus-rectum direction of the orbit.
"""
function quad_sp_coeffs(e, bx, by, which::Symbol)
    bx2, by2 = bx^2, by^2
    if which === :peri
        da = -(6bx2 * e^2 + 12bx2 * e - 3bx2 - 3by2 * e^2 + 3by2 - e^2 - 4e) / 2
        de = bx2 * (3e^3 / 8 - e^2 - 3e / 8 + 1) + by2 * (3e^3 / 8 + 4e^2 - 3e / 8 - 4) - e^3 / 4 - e^2 + e / 4 + 1
    else
        da = -(6bx2 * e^2 - 12bx2 * e - 3bx2 - 3by2 * e^2 + 3by2 - e^2 + 4e) / 2
        de = bx2 * (3e^3 / 8 + e^2 - 3e / 8 - 1) + by2 * (3e^3 / 8 - 4e^2 - 3e / 8 + 4) - e^3 / 4 + e^2 + e / 4 - 1
    end
    return da, de
end

"""Closed-form quadrupole short period as an SV5 (da, dk, dh, 0, 0) (dw, di, dOmega not included: they do not enter r_p)."""
function quad_sp_closed(F::OrbitFrame, rbs, mubs, mu, which::Symbol)
    da = 0.0; de = 0.0
    for b in eachindex(rbs)
        rb = norm(rbs[b]); u = rbs[b] / rb
        Pi = mubs[b] / mu * (F.a / rb)^3
        ca, ce = quad_sp_coeffs(F.e, dot(u, F.P), dot(u, F.Q), which)
        da += F.a * Pi * ca; de += Pi * ce
    end
    return SV5(da, F.cw * de, F.sw * de, 0.0, 0.0)
end

"""Delta r_p / a per unit Pi_3B at periapsis, closed form (derived)."""
function quad_drp_coeff(e, bx, by)
    bx2, by2 = bx^2, by^2
    return bx2 * (21e^3 / 8 + 4e^2 - 57e / 8 + 1 / 2) + by2 * (-15e^3 / 8 - 5e^2 / 2 + 15e / 8 + 5 / 2) - e^3 / 4 - e^2 / 2 + 7e / 4 - 1
end

# ---- closed-form slow (single-averaged) r_p rate, quadrupole + octupole ---------------------------------------------------

"""
Orbit-averaged rate of the periapsis radius from one perturber, divided by `n a Pi_3B`:  D2 + (a/r_b) D3, with the direction cosines
(bx, by) of the perturber along the apse line and the transverse in-plane direction (derived by Lagrange's equation from the averaged
quadrupole and octupole disturbing functions; D2 = (15/2) e eta bx by is the Kozai/Lidov form before averaging over the perturber).
"""
rp_rate_nd(e, bx, by, a_over_b) = (eta = sqrt(1 - e^2);
    7.5 * e * eta * bx * by + a_over_b * (-15 * by * eta * (30bx^2 * e^2 + 5bx^2 - 5by^2 * e^2 + 5by^2 - 3e^2 - 4) / 16))

# ---- prediction ---------------------------------------------------------------------------------------------------------

mean_rates_r(m::SV5, ctx::EstContext, s::EstSettings, rbs) = begin
    F = OrbitFrame(m[1], m[2], m[3], m[4], m[5], ctx.mu)
    dOm, dw, _ = s.j2 ? j2_rates(F, ctx.mu, ctx.J2, ctx.Re) : (0.0, 0.0, 0.0)
    inc = s.tides ? tide_rates(F, rbs, ctx.body_mu, s.n_tide) : zero(SV5)
    (inc, dOm, dw)
end

"""(a, e, i, Omega, omega) tuple of an (a, k, h, i, Omega) vector."""
vec_el(v::SV5) = (a=v[1], e=hypot(v[2], v[3]), i=v[4], raan=mod(v[5], 2π), omega=mod(atan(v[3], v[2]), 2π))

"""Third-body short-period vector at `which` of the mean orbit `m` (zero if disabled)."""
function tide_sp(m::SV5, ctx::EstContext, o::V2Opts, tab::SPTable, rbs, which::Symbol)
    (o.sp_tide === :none || isempty(rbs)) && return zero(SV5)
    F = OrbitFrame(m[1], m[2], m[3], m[4], m[5], ctx.mu)
    o.sp_tide === :quadrupole && return quad_sp_closed(F, rbs, ctx.body_mu, ctx.mu, which)
    return sp_quadrature(F, tide_accfun(rbs, ctx.body_mu), tab, which)
end

"""Osculating periapsis vector of the mean state: J2 short period (repo map or quadrature) + the supplied extra offset."""
function osc_vec2(m::SV5, ctx::EstContext, o::V2Opts, tab::SPTable, off::SV5)
    if o.sp_j2 === :map
        return osc_vec(m, ctx) + off
    end
    F = OrbitFrame(m[1], m[2], m[3], m[4], m[5], ctx.mu)
    return m + sp_quadrature(F, r -> j2_acc(r, ctx.J2, ctx.mu, ctx.Re), tab, :peri) + off
end

function shift_mean2(m::SV5, f::SV5, ctx::EstContext, o::V2Opts, tab::SPTable, off::SV5)
    base = osc_vec2(m, ctx, o, tab, off); target = apply_form(base, f); x = m + (target - base)
    for _ in 1:3
        r = target - osc_vec2(x, ctx, o, tab, off)
        x += SV5(r[1], r[2], r[3], r[4], wrap_pi(r[5]))
    end
    return x
end

"""
Predict the periapsis history (as `predict_history`) with the short-period terms of `o`. Rows: osculating periapsis state = mean + J2 SP + third-body SP.
`bodies_rb(et)` may be given to supply the perturber positions (default: SPICE through the context).
"""
function predict_history2(ctx::EstContext, s::EstSettings, o::V2Opts, tab::SPTable, seed_alt_m, inc_rad, raan_rad, omega_rad, et0;
                          max_passes::Int=typemax(Int), rho_fn=fit_rho(ctx), bodies_rb=(et -> body_positions(ctx, et)), init_override=nothing)
    out = Prediction()
    # init_override: osculating (a, e, i, raan, omega, nu) tuple at any anomaly (closed-loop predictor); default = the apoapsis seed state
    init = init_override === nothing ? initial_osculating(ctx, seed_alt_m, inc_rad, raan_rad, omega_rad) : init_override
    ivec = SV5(init.a, init.e * cos(init.omega), init.e * sin(init.omega), init.i, init.raan)
    rbs0 = s.tides ? bodies_rb(et0) : SV3[]
    d0 = zero(SV5); m = ivec; mm = init
    for _ in 1:(o.sp_tide === :none ? 1 : 3)
        el = vec_el(ivec - d0)
        mm = osculating_to_mean((a=el.a, e=el.e, i=el.i, raan=el.raan, omega=el.omega, nu=init.nu); mu=ctx.mu, radius=ctx.Re, j2=ctx.J2).elements
        m = mean_state(mm, ctx)
        d0 = tide_sp(m, ctx, o, tab, rbs0, :apo)
    end
    M0 = true_to_mean_anomaly(mm.e, mm.nu)
    t_apo = 0.0
    for kpass in 1:max_passes
        mdot = mean_motion_rate(m, ctx)
        T = 2π / mdot
        dt_in = kpass == 1 ? (2π - M0) / mdot : 0.5 * T
        dt_out = 0.5 * T
        t_k = t_apo + dt_in
        t_k > s.horizon_s && break
        et_body = et0 + t_k + (s.body_epoch === :apoapsis ? dt_out : 0.0)
        rbs = bodies_rb(et_body)
        i1, dO1, dw1 = mean_rates_r(m, ctx, s, rbs)
        i2, dO2, dw2 = mean_rates_r(advance(m, i1, dO1, dw1, dt_in), ctx, s, rbs)
        m_pre = advance(m, i2, dO2, dw2, dt_in)
        d3 = tide_sp(m_pre, ctx, o, tab, s.tides ? rbs : SV3[], :peri)
        ov = try osc_vec2(m_pre, ctx, o, tab, d3) catch; out.decayed = true; out.t_decay_s = t_k; out.ended = :invalid; break end
        osc = vec_el(ov)
        din = zero(SV5); dout = zero(SV5)
        if s.drag; din, dout = drag_deltas(osc, ctx, s, rho_fn); end
        fin = rp_form(ov, din); fout = rp_form(ov, dout)
        m_peri = try shift_mean2(m_pre, fin, ctx, o, tab, d3) catch; out.decayed = true; out.t_decay_s = t_k; out.ended = :invalid; break end
        osc = try vec_el(osc_vec2(m_peri, ctx, o, tab, d3)) catch; out.decayed = true; out.t_decay_s = t_k; out.ended = :invalid; break end
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
        m = try advance(shift_mean2(m_peri, fout, ctx, o, tab, d3), i2, dO2, dw2, dt_out) catch; out.decayed = true; out.t_decay_s = t_k; out.ended = :invalid; break end
        t_apo = t_k + dt_out
    end
    return out
end

"""Seed that makes the first-pass peak Q equal the band's geometric mean, with the v2 predictor (as `calibrate_seed`)."""
function calibrate_seed2(ctx::EstContext, s::EstSettings, o::V2Opts, tab::SPTable, inc_rad, raan_rad, omega_rad, et0; seed0=140.0e3, tol_m=1.0, max_iter=8)
    q_target = sqrt(ctx.q_low * ctx.q_high); seed = seed0
    for _ in 1:max_iter
        p = predict_history2(ctx, s, o, tab, seed, inc_rad, raan_rad, omega_rad, et0; max_passes=1)
        isempty(p.t_s) && return NaN
        h_target = altitude_of_density_fast(ctx.fit, 2 * q_target * ctx.w2m / p.vrel_mps[1]^3)
        abs(h_target - p.alt_m[1]) < tol_m && return seed
        seed += h_target - p.alt_m[1]
    end
    return seed
end

# ---- classical elements from a Cartesian state (verification / analysis) -------------------------------------------------

"""(a, k, h, i, Omega) vector of the osculating state (r, v)."""
function osc_from_state(r::SV3, v::SV3, mu)
    rn = norm(r); a = 1 / (2 / rn - dot(v, v) / mu)
    hv = cross(r, v); hn = norm(hv); W = hv / hn
    ev = cross(v, hv) / mu - r / rn
    i = acos(clamp(W[3], -1.0, 1.0)); Om = atan(W[1], -W[2])
    node = SV3(cos(Om), sin(Om), 0.0)
    k = dot(ev, node); h = dot(ev, cross(W, node))
    return SV5(a, k, h, i, Om)
end
