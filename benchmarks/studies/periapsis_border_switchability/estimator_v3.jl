# Estimator v3: estimator v2 plus CONDITIONAL higher-degree zonal harmonics (n >= 3) in the orbit-averaged dynamics.
#
# What is added (and nothing else): the zonal degrees n = 3..n_max that the Kaula-parameter criterion selects. D_n is the worst-case horizon effect of degree n
# on the periapsis radius (km, project_specs/kaula_harmonic_parameter.md section 3); the default rule `:tail` carries degree n iff the coherent sum of D_k over
# k >= n (up to the field degree) exceeds a budget (1 km), so that the omitted degrees together stay below the budget; `:per_degree` (D_n > budget) is the first
# form and under-carries when many degrees each stay below it. For every carried degree the exact zonal acceleration is averaged over the Kepler ellipse through
# the same Gauss variational equations and eccentric-anomaly nodes as the averaged third-body tide, and its rates of (a, k, h, i, Omega) are added to the
# single-averaged rates. This gives the secular and the long-period effect of every carried degree, odd (J3, J5: eccentricity and periapsis-radius drift) and even
# (J4, J6: apsidal and nodal drift) alike, in the same non-singular (k, h) state. Everything else is v2 unchanged (J2 secular rates and short-period map, third-body
# tide and short period, drag).
#
# No fitted coefficient: J_n = -sqrt(2n+1) Cbar_n0 of the SAME gravity model the truth uses (`PlanetRecord.harm`); a degree is carried only when
# the criterion says it matters; the criterion uses the geometry factor Gamma_n(e, i) of kaula_lib.py (re-derived here in Julia from the same
# integral; checked against kaula_lib in the selftest). With NO degree flagged, predict_history3 executes exactly the v2 arithmetic (the same
# mean_rates_r calls), so the output is bit-identical to predict_history2 (c0-style regression: c0_v3_regression.jl).
#
# Not included (stated limits): zonal short-period terms of degree >= 3 (about Pi_n a, 15 m for J3 at Mars), zonal contribution to the mean-anomaly
# rate (a timing effect of relative order Pi_n; at most about 1e-5 of the period), tesseral and sectorial terms (non-zero order), second order in J2.
# Usage: included by run_estimator_v3.jl, c0_v3_regression.jl and run_mars_closedloop.jl. No main().
isdefined(@__MODULE__, :predict_history2) || include(joinpath(@__DIR__, "estimator_v2.jl"))

# ---- zonal field of a planet record ------------------------------------------------------------------------------------

"""Zonal data of the truth gravity model: `J[n]` for n = 2..L (unnormalized, J_n = -sqrt(2n+1) Cbar_n0), reference radius, GM, and the Kaula
constant K of the model's own coefficients (geometric mean of A_n n^2 over 10 <= n <= L; NaN when L < 11)."""
struct ZonalField
    L::Int
    J::Vector{Float64}        # J[n], n = 1..L (J[1] = 0)
    Re::Float64
    mu::Float64
    K::Float64
end

function zonal_field(harm)
    L = Int(harm.L); J = zeros(L)
    for n in 2:L; J[n] = -sqrt(2n + 1.0) * harm.C[n + 1, 1]; end
    As = Float64[]
    for n in 10:L
        s = sum(harm.C[n + 1, m + 1]^2 + harm.S[n + 1, m + 1]^2 for m in 0:n)
        s > 0 && push!(As, sqrt(s / (2n + 1)) * n^2)
    end
    K = length(As) >= 2 ? exp(sum(log, As) / length(As)) : NaN
    return ZonalField(L, J, harm.reference_radius_m, harm.gm_m3s2, K)
end

"""J_n of the field; beyond the field's degree the Kaula rms value sqrt(2n+1) K / n^2 (NaN when the field has no K)."""
zonal_J(zf::ZonalField, n::Int) = n <= zf.L ? zf.J[n] : sqrt(2n + 1.0) * zf.K / n^2

# ---- exact zonal acceleration ---------------------------------------------------------------------------------------------

"""Acceleration of the zonal harmonics `ns` (with coefficients `Js`) at `r` (body-equator frame, z = pole): a = grad Phi,
Phi = -mu sum_n J_n Re^n P_n(z/r) / r^(n+1) (the J2/J3 forms of the repo are the n = 2, 3 cases; checked in the selftest)."""
@inline function zonal_acc(r::SV3, ns::Vector{Int}, Js::Vector{Float64}, mu, Re)
    x, y, z = r; rn = sqrt(x * x + y * y + z * z); s = z / rn; rho = Re / rn; nmax = ns[end]
    acc = zero(SV3); Pm = 1.0; P = s; dPm = 0.0; dP = 1.0; rhon = rho; q = 1; mur3 = mu / (rn * rn * rn)
    for k in 1:nmax
        if k > 1   # advance to degree k
            Pn = ((2k - 1) * s * P - (k - 1) * Pm) / k; dPn = k * P + s * dP
            Pm, P, dPm, dP = P, Pn, dP, dPn; rhon *= rho
        end
        if q <= length(ns) && ns[q] == k
            G = s * dP + (k + 1) * P
            acc += (-mur3 * Js[q] * rhon) * (-G * r + SV3(0.0, 0.0, rn * dP)); q += 1
        end
    end
    return acc
end

"""Orbit-averaged Gauss rates (per second, uniform in mean anomaly) of (a, k, h, i, Omega) from the zonal harmonics `ns`, `N` eccentric-anomaly nodes."""
function zonal_rates(F::OrbitFrame, ns::Vector{Int}, Js::Vector{Float64}, mu, Re, N::Int)
    acc = zero(SV5)
    for j in 1:N
        E = 2π * (j - 0.5) / N
        r, _, rm, cnu, snu, d = node_state(F, E)
        f = zonal_acc(r, ns, Js, mu, Re); rh = r / rm; th = cross(F.W, rh)
        acc += d * gauss_g(F, rm, cnu, snu, dot(f, rh), dot(f, th), dot(f, F.W))
    end
    return acc / N
end

# ---- the criterion (absolute Kaula-parameter form) -----------------------------------------------------------------------------

"""Options of the zonal terms. `budget_km`: carry degree n iff D_n > budget_km (ASSUMED, the oracle margin of 1 km, as the Kaula pre-registration).
`n_max`: highest degree scanned (0 = the field's degree). `gamma_nq`, `gamma_nw`: nodes of the Gamma_n quadrature in true anomaly and the omega grid
of its maximum (numerical resolutions). `nodes_min`, `nodes_per_degree`: eccentric-anomaly nodes of the averaged rate = max(nodes_min, nodes_per_degree * n_carried_max)
(numerical resolution, checked in the selftest). `force`: carry exactly these degrees (diagnostic; `Int[]` = none, i.e. v2). `rule`: `:tail` (default, cumulative: the Kaula Pi_K(N) argument, the sum of the omitted degrees must be below the budget) or `:per_degree` (the first form, D_n > budget; under-carries when many degrees each stay below the budget)."""
Base.@kwdef struct ZonalOpts
    budget_km::Float64 = 1.0
    n_max::Int = 0
    gamma_nq::Int = 256
    gamma_nw::Int = 256
    nodes_min::Int = 64
    nodes_per_degree::Int = 12
    force::Union{Nothing, Vector{Int}} = nothing
    rule::Symbol = :tail                  # :tail (carry n iff sum_{k>=n}^{L} D_k > budget: the omitted degrees add up) or :per_degree (carry n iff D_n > budget)
end

"""Legendre polynomial and derivative P_n(x), P_n'(x) by the three-term recurrences (no division by 1 - x^2)."""
@inline function legendre_pd(n::Int, x::Float64)
    n == 0 && return 1.0, 0.0
    Pm, P, dP = 1.0, x, 1.0
    for k in 1:n-1
        Pn = ((2k + 1) * x * P - k * Pm) / (k + 1); dP = (k + 1) * P + x * dP; Pm, P = P, Pn
    end
    return P, dP
end

"""max over omega of |Gamma_n(e, i, omega)|, Gamma_n = ((1-e^2)/e) <(1 + e cos nu)^(n-1) P_n'(sin i sin(omega+nu)) sin i cos(omega+nu)>_nu
(|de/dt| <= N |Pi_n| Gamma_max, Pi_n = J_n (Re/p)^n; kaula_lib.gamma_max). e is floored at 1e-4."""
function gamma_max(n::Int, e::Float64, inc::Float64, nq::Int, nw::Int)
    e = max(e, 1e-4); si = sin(inc); best = 0.0
    nu = [(j - 0.5) * 2π / nq for j in 1:nq]; w1 = [(1 + e * cos(v))^(n - 1) for v in nu]
    for iw in 0:nw-1
        w = 2π * iw / nw; acc = 0.0
        for j in 1:nq
            u = w + nu[j]; _, dP = legendre_pd(n, si * sin(u)); acc += w1[j] * dP * si * cos(u)
        end
        best = max(best, abs((1 - e^2) / e * acc / nq))
    end
    return best
end

const _GAMMA_CACHE = Dict{NTuple{4, Float64}, Float64}()   # (n, round(e, 3), round(i, 4), nq*1000+nw) -> Gamma_max; (e, i) are the only inputs
gamma_max_cached(n, e, inc, o::ZonalOpts) = get!(_GAMMA_CACHE, (Float64(n), round(e; digits=3), round(inc; digits=4), Float64(o.gamma_nq * 1000 + o.gamma_nw))) do
    gamma_max(n, round(e; digits=3), round(inc; digits=4), o.gamma_nq, o.gamma_nw)
end

"""Zonal set of an orbit: degrees n in 3..n_max selected by `o.rule` from D_n = a N |Pi_n| Gamma_max T_eff / 1000 (km), T_eff = min(T_h, 2/|omega_dot_J2|),
Pi_n = J_n (Re/p)^n, evaluated at the initial osculating elements (a, e, i). Returns (ns, Js, D) with D the full vector of D_n (n = 3..n_max, km)."""
function zonal_set(zf::ZonalField, o::ZonalOpts, a, e, inc, horizon_s)
    nmax = o.n_max > 0 ? o.n_max : zf.L; mu = zf.mu
    p = a * (1 - e^2); N = sqrt(mu / a^3)
    wdot2 = abs(0.75 * N * zf.J[2] * (zf.Re / p)^2 * (5 * cos(inc)^2 - 1))      # |J2 apsidal rate|, as j2_rates
    Teff = wdot2 * horizon_s < 2 ? horizon_s : 2 / wdot2
    D = zeros(max(nmax - 2, 0))
    for n in 3:nmax
        Jn = zonal_J(zf, n); isfinite(Jn) || continue
        D[n - 2] = a * N * abs(Jn) * (zf.Re / p)^n * gamma_max_cached(n, e, inc, o) * Teff / 1.0e3
    end
    tail = reverse(cumsum(reverse(D)))                          # tail[n - 2] = sum_{k >= n}^{nmax} D_k
    sel = o.rule === :tail ? tail : D
    ns = o.force === nothing ? [n for n in 3:nmax if sel[n - 2] > o.budget_km] : sort(unique(o.force))
    return ns, [zonal_J(zf, n) for n in ns], D
end

# ---- prediction ---------------------------------------------------------------------------------------------------------

"""Zonal part of the single-averaged rates in the form the stepper needs: `(inc, dw)`, `inc` an SV5 of (a, k, h, i, Omega) rates whose (k, h) part is purely
RADIAL (de/dt along the current apse direction) and `dw` the zonal apsidal rate. The raw Gauss (k, h) rate of a zonal degree is dominated by its apsidal part
(e * domega/dt is 1-100 times de/dt), and `advance` evaluates the second rate at the end-of-step state, so an additive (k, h) increment of the raw vector
would leak the large perpendicular part into e (13% of the J4 rate at Mars 25.2 deg, found by the Cowell self-test). The apsidal part is therefore returned as a
rate of omega and applied by the exact rotation of `advance`."""
function zonal_inc_dw(F::OrbitFrame, ns::Vector{Int}, Js::Vector{Float64}, mu, Re, nz::Int)
    r = zonal_rates(F, ns, Js, mu, Re, nz)
    F.e > 1e-8 || return r, 0.0
    ed = F.cw * r[2] + F.sw * r[3]; wz = (F.cw * r[3] - F.sw * r[2]) / F.e
    return SV5(r[1], F.cw * ed, F.sw * ed, r[4], r[5]), wz
end

"""Single-averaged rates with the carried zonal degrees added: `(inc, dOm, dw)` as `mean_rates_r`. Empty `ns`: exactly `mean_rates_r`."""
function mean_rates_r3(m::SV5, ctx::EstContext, s::EstSettings, rbs, ns::Vector{Int}, Js::Vector{Float64}, nz::Int)
    isempty(ns) && return mean_rates_r(m, ctx, s, rbs)
    F = OrbitFrame(m[1], m[2], m[3], m[4], m[5], ctx.mu)
    dOm, dw, _ = s.j2 ? j2_rates(F, ctx.mu, ctx.J2, ctx.Re) : (0.0, 0.0, 0.0)
    inc = s.tides ? tide_rates(F, rbs, ctx.body_mu, s.n_tide) : zero(SV5)
    zinc, wz = zonal_inc_dw(F, ns, Js, ctx.mu, ctx.Re, nz)
    return (inc + zinc, dOm, dw + wz)
end

"""
`predict_history2` with the conditional zonal terms. The extra return value `zonal` is `(ns, D_km)`: the carried degrees and the criterion values.
Same keywords as `predict_history2`; `zf` is the zonal field of the truth gravity model, `zo` the options. With no carried degree the arithmetic is v2's.
"""
function predict_history3(ctx::EstContext, s::EstSettings, o::V2Opts, tab::SPTable, zf::ZonalField, zo::ZonalOpts, seed_alt_m, inc_rad, raan_rad, omega_rad, et0;
                          max_passes::Int=typemax(Int), rho_fn=fit_rho(ctx), bodies_rb=(et -> body_positions(ctx, et)), init_override=nothing)
    out = Prediction()
    init = init_override === nothing ? initial_osculating(ctx, seed_alt_m, inc_rad, raan_rad, omega_rad) : init_override
    ns, Js, D = zonal_set(zf, zo, init.a, init.e, init.i, s.horizon_s)
    nz = isempty(ns) ? 0 : max(zo.nodes_min, zo.nodes_per_degree * ns[end])
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
        i1, dO1, dw1 = mean_rates_r3(m, ctx, s, rbs, ns, Js, nz)
        i2, dO2, dw2 = mean_rates_r3(advance(m, i1, dO1, dw1, dt_in), ctx, s, rbs, ns, Js, nz)
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
    return out, (ns, D)
end

"""Seed that makes the first-pass peak Q equal the band's geometric mean, with the v3 predictor (as `calibrate_seed2`)."""
function calibrate_seed3(ctx::EstContext, s::EstSettings, o::V2Opts, tab::SPTable, zf::ZonalField, zo::ZonalOpts, inc_rad, raan_rad, omega_rad, et0; seed0=140.0e3, tol_m=1.0, max_iter=8)
    q_target = sqrt(ctx.q_low * ctx.q_high); seed = seed0
    for _ in 1:max_iter
        p, _ = predict_history3(ctx, s, o, tab, zf, zo, seed, inc_rad, raan_rad, omega_rad, et0; max_passes=1)
        isempty(p.t_s) && return NaN
        h_target = altitude_of_density_fast(ctx.fit, 2 * q_target * ctx.w2m / p.vrel_mps[1]^3)
        abs(h_target - p.alt_m[1]) < tol_m && return seed
        seed += h_target - p.alt_m[1]
    end
    return seed
end
