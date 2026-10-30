# UNFROZEN DIAGNOSTIC (stage C2). Estimator v2 plus the orbit-averaged J3 zonal harmonic, to test the hypothesis that the omission of J3 explains the
# Mars failures. Nothing here is part of the frozen estimator and no result computed with it enters a pass/fail verdict.
#
# The J3 term is derived, not fitted: the exact J3 acceleration (zonal harmonic of degree 3, potential  mu J3 R^3 P3(z/r) / r^4  for the planet-fixed
# convention of the truth gravity file) is averaged over the Kepler ellipse through the same Gauss variational equations and the same eccentric-anomaly
# nodes as the averaged third-body tide (`tide_rates`); it adds secular and long-period rates of (a, k, h, i, Omega) to `mean_rates_r`.
# The J3 short-period terms and the mean-element convention for J3 are not included.
# Usage: included by run_mars_generality.jl. `predict_history2_j3` is a copy of `predict_history2` with the J3 rates added to the averaged rates.
isdefined(@__MODULE__, :predict_history2) || include(joinpath(@__DIR__, "estimator_v2.jl"))

"""Zonal J3 acceleration at the position `r` (body-equator frame, z = pole): -grad of Phi3 = mu J3 R^3 P3(z/r) / r^4."""
@inline function j3_acc(r::SV3, J3, mu, Re)
    x, y, z = r; r2 = x * x + y * y + z * z; rn = sqrt(r2); ri5 = 1 / (r2 * r2 * rn); ri7 = ri5 / r2; ri9 = ri7 / r2
    k = -mu * J3 * Re^3 / 2          # a = -mu J3 R^3 grad f,  f = (5 z^3 r^-7 - 3 z r^-5) / 2
    g = -35 * z^3 * ri9 + 15 * z * ri7
    return k * 2 * SV3(x * g / 2, y * g / 2, (30 * z^2 * ri7 - 35 * z^4 * ri9 - 3 * ri5) / 2)
end

"""Orbit-averaged Gauss rates (per second, uniform in mean anomaly) of (a, k, h, i, Omega) from the acceleration function `accfun`."""
function accel_rates(F::OrbitFrame, accfun, N::Int)
    acc = zero(SV5)
    for j in 1:N
        E = 2π * (j - 0.5) / N
        r, _, rm, cnu, snu, d = node_state(F, E)
        f = accfun(r); rh = r / rm; th = cross(F.W, rh)
        acc += d * gauss_g(F, rm, cnu, snu, dot(f, rh), dot(f, th), dot(f, F.W))
    end
    return acc / N
end

mean_rates_r_j3(m::SV5, ctx::EstContext, s::EstSettings, rbs, J3) = begin
    F = OrbitFrame(m[1], m[2], m[3], m[4], m[5], ctx.mu)
    dOm, dw, _ = s.j2 ? j2_rates(F, ctx.mu, ctx.J2, ctx.Re) : (0.0, 0.0, 0.0)
    inc = s.tides ? tide_rates(F, rbs, ctx.body_mu, s.n_tide) : zero(SV5)
    inc += accel_rates(F, r -> j3_acc(r, J3, ctx.mu, ctx.Re), s.n_tide)
    (inc, dOm, dw)
end

function predict_history2_j3(ctx::EstContext, s::EstSettings, o::V2Opts, tab::SPTable, J3, seed_alt_m, inc_rad, raan_rad, omega_rad, et0;
                             max_passes::Int=typemax(Int), rho_fn=fit_rho(ctx), bodies_rb=(et -> body_positions(ctx, et)))
    out = Prediction()
    init = initial_osculating(ctx, seed_alt_m, inc_rad, raan_rad, omega_rad)
    ivec = SV5(init.a, init.e * cos(init.omega), init.e * sin(init.omega), init.i, init.raan)
    rbs0 = s.tides ? bodies_rb(et0) : SV3[]
    d0 = zero(SV5); m = ivec; mm = init
    for _ in 1:(o.sp_tide === :none ? 1 : 3)
        el = vec_el(ivec - d0)
        mm = osculating_to_mean((a=el.a, e=el.e, i=el.i, raan=el.raan, omega=el.omega, nu=init.nu); mu=ctx.mu, radius=ctx.Re, j2=ctx.J2).elements
        m = mean_state(mm, ctx); d0 = tide_sp(m, ctx, o, tab, rbs0, :apo)
    end
    M0 = true_to_mean_anomaly(mm.e, mm.nu); t_apo = 0.0
    for kpass in 1:max_passes
        mdot = mean_motion_rate(m, ctx); T = 2π / mdot
        dt_in = kpass == 1 ? (2π - M0) / mdot : 0.5 * T; dt_out = 0.5 * T
        t_k = t_apo + dt_in
        t_k > s.horizon_s && break
        et_body = et0 + t_k + (s.body_epoch === :apoapsis ? dt_out : 0.0)
        rbs = bodies_rb(et_body)
        i1, dO1, dw1 = mean_rates_r_j3(m, ctx, s, rbs, J3)
        i2, dO2, dw2 = mean_rates_r_j3(advance(m, i1, dO1, dw1, dt_in), ctx, s, rbs, J3)
        m_pre = advance(m, i2, dO2, dw2, dt_in)
        d3 = tide_sp(m_pre, ctx, o, tab, s.tides ? rbs : SV3[], :peri)
        ov = try osc_vec2(m_pre, ctx, o, tab, d3) catch; out.decayed = true; out.t_decay_s = t_k; out.ended = :invalid; break end
        osc = vec_el(ov); din = zero(SV5); dout = zero(SV5)
        if s.drag; din, dout = drag_deltas(osc, ctx, s, rho_fn); end
        fin = rp_form(ov, din); fout = rp_form(ov, dout)
        m_peri = try shift_mean2(m_pre, fin, ctx, o, tab, d3) catch; out.decayed = true; out.t_decay_s = t_k; out.ended = :invalid; break end
        osc = try vec_el(osc_vec2(m_peri, ctx, o, tab, d3)) catch; out.decayed = true; out.t_decay_s = t_k; out.ended = :invalid; break end
        rp = osc.a * (1 - osc.e)
        if rp - ctx.Rp_e <= s.impact_altitude_m; out.decayed = true; out.t_decay_s = t_k; out.ended = :impact; break; end
        rv = elements_to_cartesian((a=osc.a, e=osc.e, i=osc.i, raan=osc.raan, omega=osc.omega, nu=0.0); mu=ctx.mu)
        r = SV3(rv[1]); v = SV3(rv[2]); alt, lat = geodetic_altitude_latitude(r, ctx.planet)
        vrel = norm(v - cross(ctx.omega_e, r)); rho = corridor_density(ctx.fit, alt); hq, hl = band_edges(ctx, vrel)
        for (arr, val) in zip((out.t_s, out.alt_m, out.radius_m, out.lat_rad, out.vrel_mps, out.rho, out.q_w_cm2, out.a_m, out.e, out.i_rad, out.raan_rad, out.omega_rad, out.h_qhigh_m, out.h_qlow_m),
                              (t_k, alt, rp, lat, vrel, rho, 0.5 * rho * vrel^3 / ctx.w2m, osc.a, osc.e, osc.i, osc.raan, osc.omega, hq, hl))
            push!(arr, val)
        end
        m = try advance(shift_mean2(m_peri, fout, ctx, o, tab, d3), i2, dO2, dw2, dt_out) catch; out.decayed = true; out.t_decay_s = t_k; out.ended = :invalid; break end
        t_apo = t_k + dt_out
    end
    return out
end
