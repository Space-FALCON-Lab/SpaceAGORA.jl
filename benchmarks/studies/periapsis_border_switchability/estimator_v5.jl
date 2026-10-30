# Estimator v5 = estimator v4 + CORRECTED PASS TIMING (study-local; v1 to v4 stay frozen).
#
# Why (estimator_v5_results.md section 2; measured by truth-state injection and 40-pass chains against a Cowell propagator of the same force subsets that reproduces the
# engine truth to 2 ms per orbit): the pass-time error of v1 to v4 is the accumulated error of the drag-driven change of the semi-major axis, not of the period formula (the
# two-body + J2 anomalistic period is exact to 0.01 s per orbit without drag). Three defects, in the order of their size:
#  T1 the drag increments of v1 to v4 are the first-order (Born) Gauss integral over the KEPLER ellipse of the undisturbed periapsis state. Along the real path the drag itself
#     bends the orbit (the speed falls by about 1 m/s during the inbound half, J2 changes the curvature by 0.35% of the radial acceleration), so the density seen and the work
#     done differ at second order: +0.29% (inbound) and +0.94% (outbound) too much change of a in the two-body limit at Mars, 0.84% (median) with J2. v5 integrates the real window
#     path numerically (point mass + J2 + drag, same fit and relative velocity), started from the undisturbed periapsis state propagated backwards with J2 only, and accumulates
#     the drag contributions to the invariants of the osculating orbit along it (exact for the actual path):
#         dE/dt = v . f,   dh/dt = r x f,   de/dt = (f x h + v x (r x f)) / mu     (f = drag acceleration)
#     split at the periapsis epoch into the inbound and outbound increments of (a, r_p, omega, i, Omega) that `shift_mean` consumes unchanged (`rp_form` of v1 to v4).
#  T2 the outbound half of the tesseral kick of v4 was evaluated on the PRE-drag ellipse (the drag removes 1e5 m of a per pass at Mars: apoapsis end 0.5% off, time against the
#     planet's rotation 0.7% slow); v5 evaluates it on the post-drag ellipse (mean elements without the inbound transient) with its own timing term.
#  T3 an uncarried pass (no kick: Earth) took both legs from the pre-pass period (v1 to v3); v5 uses half the post-drag period for the outbound leg, as the carried path of v4 does.
# Tested and REJECTED on the development cells (named option, default off): `zonal_sp`, the first-order short-period terms of the carried zonal degrees in the osculating state.
# Everything else (tide, J2, zonal rates, short-period map, tesseral kick, carry rule, Q band, calibration) is v4 unchanged. No fitted coefficient.
#
# Identities: `TimingOpts(window=:kepler, post_drag_period=false, kick_post_drag=false, zonal_sp=false)` returns exactly the v3 prediction where the carry rule is false and
# exactly the v4 prediction where it is true (selftest R1 of run_estimator_v5.jl).
# Usage: included by run_estimator_v5.jl, diag_v5_timing.jl, diag_v5_earth.jl and c0_v5_regression.jl. No main().
isdefined(@__MODULE__, :predict_history4) || include(joinpath(@__DIR__, "estimator_v4.jl"))

"""Options of the corrected drag arc. `window`: `:integrated` (this file) or `:kepler` (v4). `n_window`: RK4 steps of the drag integration over the window, nodes clustered
at periapsis by a sinh map with parameter `cluster` (numerical resolution; convergence checked in the diagnostics). `n_back`: RK4 steps of the J2-only back-propagation
from the undisturbed periapsis to the window entry. `window_j2`: include the J2 acceleration in the window (the interplay of J2 and drag was 0.5% of the pass effect)."""
Base.@kwdef struct TimingOpts
    window::Symbol = :integrated
    n_window::Int = 32
    n_back::Int = 8
    cluster::Float64 = 1.0
    window_j2::Bool = true
    kick_post_drag::Bool = true        # the OUTBOUND half of the tesseral kick is evaluated on the post-drag ellipse (`kick_arc` with the elements after the drag step)
    zonal_sp::Bool = false             # first-order SHORT-PERIOD terms of the carried zonal degrees (n >= 3) in the osculating periapsis/apoapsis state (the quadrature of v2, with the zonal acceleration)
    post_drag_period::Bool = true      # the outbound leg of an UNCARRIED pass lasts half the post-drag period (the carried path uses `TesseralOpts.post_drag_period`, as v4)
end

# ---- the window integration -------------------------------------------------------------------------------------------------

"""Right-hand side of the window integration: state y = (r, v, dE, dh, de) (3+3+1+3+3), `drag` switches the drag acceleration (and the accumulators) on."""
@inline function window_rhs(y::SVector{13, Float64}, ctx::EstContext, rho_fn, J2, drag::Bool)
    r = SV3(y[1], y[2], y[3]); v = SV3(y[4], y[5], y[6]); mu = ctx.mu
    rn = norm(r); a = -(mu / rn^3) * r
    J2 != 0 && (a += j2_acc(r, J2, mu, ctx.Re))
    if drag
        rho = rho_fn(r); vrel = v - cross(ctx.omega_e, r); V = norm(vrel)
        f = (-0.5 * rho * V / ctx.beta) * vrel; h = cross(r, v)
        dE = dot(v, f); dh = cross(r, f); de = (cross(f, h) + cross(v, dh)) / mu
        return SVector{13, Float64}(v[1], v[2], v[3], a[1] + f[1], a[2] + f[2], a[3] + f[3], dE, dh[1], dh[2], dh[3], de[1], de[2], de[3])
    end
    return SVector{13, Float64}(v[1], v[2], v[3], a[1], a[2], a[3], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

@inline function rk4_step(y::SVector{13, Float64}, h::Float64, ctx, rho_fn, J2, drag::Bool)
    k1 = window_rhs(y, ctx, rho_fn, J2, drag); k2 = window_rhs(y + (0.5 * h) * k1, ctx, rho_fn, J2, drag)
    k3 = window_rhs(y + (0.5 * h) * k2, ctx, rho_fn, J2, drag); k4 = window_rhs(y + h * k3, ctx, rho_fn, J2, drag)
    return y + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
end

"""Elements `(a, e, i, raan, omega)` from the invariants: specific energy `En`, angular momentum `h`, eccentricity vector `ev`."""
function invariants_to_elements(En, hv::SV3, ev::SV3, mu)
    a = -mu / (2 * En); e = norm(ev); hn = norm(hv); W = hv / hn
    i = acos(clamp(W[3], -1.0, 1.0)); Om = atan(W[1], -W[2]); node = SV3(cos(Om), sin(Om), 0.0)
    om = atan(dot(ev, cross(W, node)), dot(ev, node))
    return (a=a, e=e, i=i, raan=Om, omega=om)
end

"""
Drag increments of one pass `(fin, fout)` as `rp_form`-type SV5 `(da, drp, domega, di, dOmega)` of the osculating periapsis elements `osc` (the undisturbed periapsis
state of the pass, as v1 to v4's `drag_deltas` input), by integrating the real perturbed window path (see the file header). Inbound = before the periapsis epoch.
"""
function drag_window(osc, ctx::EstContext, s::EstSettings, rho_fn, to::TimingOpts)
    mu = ctx.mu; J2 = to.window_j2 ? ctx.J2 : 0.0
    F = OrbitFrame(osc.a, osc.e * cos(osc.omega), osc.e * sin(osc.omega), osc.i, osc.raan, mu)
    c = 1 - s.drag_window_m / (F.a * F.e); Emax = c <= -1 ? float(π) : acos(c)
    r0, v0, _, _, _, _ = node_state(F, 0.0)
    h0 = cross(r0, v0); ev0 = cross(v0, h0) / mu - r0 / norm(r0); En0 = 0.5 * dot(v0, v0) - mu / norm(r0)
    shp = sinh(to.cluster)
    tgrid(n, sgn) = [(E = sgn * Emax * sinh(to.cluster * j / n) / shp; (E - F.e * sin(E)) / F.n) for j in 0:n]       # |t| from 0 (periapsis) to the window edge
    # entry state: backward from the periapsis, J2 only (no drag), RK4 on the time grid
    y = SVector{13, Float64}(r0..., v0..., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    tb = tgrid(to.n_back, -1.0)
    for j in 1:to.n_back; y = rk4_step(y, tb[j + 1] - tb[j], ctx, rho_fn, J2, false); end
    # inbound: entry -> periapsis epoch with drag; outbound: periapsis epoch -> exit
    tin = -reverse(tgrid(to.n_window, 1.0))               # increasing times from -t_edge to 0
    for j in 1:to.n_window; y = rk4_step(y, tin[j + 1] - tin[j], ctx, rho_fn, J2, true); end
    yin = y
    tout = tgrid(to.n_window, 1.0)
    for j in 1:to.n_window; y = rk4_step(y, tout[j + 1] - tout[j], ctx, rho_fn, J2, true); end
    el(yy) = invariants_to_elements(En0 + yy[7], h0 + SV3(yy[8], yy[9], yy[10]), ev0 + SV3(yy[11], yy[12], yy[13]), mu)
    e0 = invariants_to_elements(En0, h0, ev0, mu); e1 = el(yin); e2 = el(y)
    form(x, z) = SV5(z.a - x.a, z.a * (1 - z.e) - x.a * (1 - x.e), wrap_pi(z.omega - x.omega), z.i - x.i, wrap_pi(z.raan - x.raan))
    return form(e0, e1), form(e1, e2)
end

# ---- the stepper -----------------------------------------------------------------------------------------------------------------

"""
`predict_history4` with the corrected pass timing. Same arguments plus `to::TimingOpts`; works with and without the tesseral kick (the carry rule of v4 decides, `tt.force`
overrides) so that it also corrects the zonal-only path used at Earth. Returns `(prediction, info)` as v4: `info = (ns, D, DT_km, carried)`.
Identities: `to.window == :kepler && !to.post_drag_period` returns the v3 object when the rule is false (v5 == v3 where no correction is on) and `to.window == :kepler`
reproduces v4 when the rule is true; `trace` (a Vector) receives `(k, t_apo_start, m_at_start)` per pass for the diagnostics.
"""
function predict_history5(ctx::EstContext, s::EstSettings, o::V2Opts, tab::SPTable, zf::ZonalField, zo::ZonalOpts, tf::TesseralField, tt::TesseralOpts, to::TimingOpts,
                          seed_alt_m, inc_rad, raan_rad, omega_rad, et0;
                          max_passes::Int=typemax(Int), rho_fn=fit_rho(ctx), bodies_rb=(et -> body_positions(ctx, et)), init_override=nothing, DT_override=nothing, theta_times=nothing, trace=nothing)
    base, info3 = predict_history3(ctx, s, o, tab, zf, zo, seed_alt_m, inc_rad, raan_rad, omega_rad, et0; max_passes=max_passes, rho_fn=rho_fn, bodies_rb=bodies_rb, init_override=init_override)
    legacy = to.window === :kepler && !to.post_drag_period && !to.zonal_sp
    tt.force === false && legacy && return base, (info3[1], info3[2], NaN, false)
    spin = norm(ctx.omega_e)
    DT = tt.force === false ? NaN : (DT_override !== nothing ? DT_override : (tt.force === true ? NaN : dt_of_prediction(tf, base, spin, tt.nodes_dt)))
    kick = tt.force === false ? false : (tt.force === true ? true : DT > tt.budget_km)
    (!kick && legacy) && return base, (info3[1], info3[2], DT, false)
    out = Prediction()
    Lt = tt.Lt > 0 ? min(tt.Lt, tf.L) : tf.L; Mt = tt.Mt > 0 ? min(tt.Mt, Lt) : Lt
    rotm = kick ? RotationModel(ctx, et0; kind=tt.rotation) : nothing
    init = init_override === nothing ? initial_osculating(ctx, seed_alt_m, inc_rad, raan_rad, omega_rad) : init_override
    ns, D = info3[1], info3[2]; Js = [zonal_J(zf, n) for n in ns]
    nz = isempty(ns) ? 0 : max(zo.nodes_min, zo.nodes_per_degree * ns[end])
    ivec = SV5(init.a, init.e * cos(init.omega), init.e * sin(init.omega), init.i, init.raan)
    rbs0 = s.tides ? bodies_rb(et0) : SV3[]
    ztab = (to.zonal_sp && !isempty(ns)) ? SPTable(nz) : nothing
    zsp(mean::SV5, which::Symbol) = ztab === nothing ? zero(SV5) :
        sp_quadrature(OrbitFrame(mean[1], mean[2], mean[3], mean[4], mean[5], ctx.mu), r -> zonal_acc(r, ns, Js, ctx.mu, ctx.Re), ztab, which)
    d0 = zero(SV5); m = ivec; mm = init
    for _ in 1:(o.sp_tide === :none && ztab === nothing ? 1 : 3)
        el = vec_el(ivec - d0)
        mm = osculating_to_mean((a=el.a, e=el.e, i=el.i, raan=el.raan, omega=el.omega, nu=init.nu); mu=ctx.mu, radius=ctx.Re, j2=ctx.J2).elements
        m = mean_state(mm, ctx)
        d0 = tide_sp(m, ctx, o, tab, rbs0, :apo) + zsp(m, :apo)
    end
    M0 = true_to_mean_anomaly(mm.e, mm.nu)
    f_start0 = mm.nu - 2π; f_start0 > 0 && (f_start0 -= 2π)
    post_period = kick ? tt.post_drag_period : to.post_drag_period
    t_apo = 0.0
    for kpass in 1:max_passes
        trace === nothing || push!(trace, (k=kpass, t_apo=t_apo, m=m))
        mdot = mean_motion_rate(m, ctx)
        T = 2π / mdot
        dt_in = kpass == 1 ? (2π - M0) / mdot : 0.5 * T
        dt_out = 0.5 * T
        t_k0 = t_apo + dt_in
        t_k0 > s.horizon_s && break
        kk = nothing; th_k = 0.0; wr_k = 0.0
        if kick
            t_rot = (theta_times !== nothing && kpass <= length(theta_times)) ? theta_times[kpass] : t_k0
            if to.kick_post_drag
                kk = try
                    el = mean_elements(m); th_k, wr_k = rotation_at(rotm, t_rot); kw = kick_work(tf.L); f0 = kpass == 1 ? f_start0 : -float(π)
                    nin = max(8, round(Int, tt.nodes * (-f0) / (2π)))
                    kin = kick_arc(tf, kw, Lt, Mt, el[1], el[2], el[3], el[4], el[5], th_k, wr_k, f0, 0.0, nin, tt.e_min)
                    (dah=kin.da, deh=kin.de, dih=kin.di, dOh=kin.dO, dwh=kin.dw, dt_in=1.5 * (-kin.datau) / el[1] - kin.dM / kin.Nm, drp=kin.da)
                catch; nothing end
            else
                kk = try tesseral_pass_kick(mean_elements(m), t_rot, tf, rotm; nodes=tt.nodes, Lt=Lt, Mt=Mt, f_start=kpass == 1 ? f_start0 : -float(π), e_min=tt.e_min) catch; nothing end
            end
            if kk === nothing || !isfinite(kk.drp)
                out.decayed = true; out.t_decay_s = t_k0; out.ended = :invalid; break
            end
        end
        dtk_in = kick && tt.dt_terms ? kk.dt_in : 0.0; dtk_arc = kick && tt.dt_terms && !to.kick_post_drag ? kk.dt_arc : 0.0
        et_body = et0 + t_k0 + (s.body_epoch === :apoapsis ? dt_out : 0.0)
        rbs = bodies_rb(et_body)
        i1, dO1, dw1 = mean_rates_r3(m, ctx, s, rbs, ns, Js, nz)
        i2, dO2, dw2 = mean_rates_r3(advance(m, i1, dO1, dw1, dt_in), ctx, s, rbs, ns, Js, nz)
        m_pre = advance(m, i2, dO2, dw2, dt_in)
        kick && (m_pre = add_kick(m_pre, kk.dah, kk.deh, kk.dih, kk.dOh, kk.dwh))
        t_k = t_k0 + dtk_in
        d3 = tide_sp(m_pre, ctx, o, tab, s.tides ? rbs : SV3[], :peri) + zsp(m_pre, :peri)
        ov = try osc_vec2(m_pre, ctx, o, tab, d3) catch; out.decayed = true; out.t_decay_s = t_k; out.ended = :invalid; break end
        osc = vec_el(ov)
        fin = zero(SV5); fout = zero(SV5)
        if s.drag
            if to.window === :integrated
                fin, fout = drag_window(osc, ctx, s, rho_fn, to)
            else
                din, dout = drag_deltas(osc, ctx, s, rho_fn); fin = rp_form(ov, din); fout = rp_form(ov, dout)
            end
        end
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
        kout = nothing; dt_ao = 0.0
        m = try
            m_post = shift_mean2(m_peri, fout, ctx, o, tab, d3)
            if post_period
                mb = kick && tt.base_period_without_kick ? add_kick(m_post, -kk.dah, -kk.deh, -kk.dih, -kk.dOh, -kk.dwh) : m_post
                dt_out = 0.5 * 2π / mean_motion_rate(mb, ctx)
            end
            if kick && to.kick_post_drag
                # outbound half on the POST-drag ellipse: the pre-drag arc is 0.5% too large at the apoapsis end (the drag removes 1e5 m of a per pass at Mars) and its time
                # runs 0.7% too slowly against the planet's rotation; the base orbit is the mean orbit without the transient inbound half kick (as the period above)
                el = mean_elements(add_kick(m_post, -kk.dah, -kk.deh, -kk.dih, -kk.dOh, -kk.dwh)); kw = kick_work(tf.L)
                kout = kick_arc(tf, kw, Lt, Mt, el[1], el[2], el[3], el[4], el[5], th_k, wr_k, 0.0, float(π), max(8, tt.nodes ÷ 2), tt.e_min)
                isfinite(kout.da) || error("non-finite outbound kick")
                Tend = π / kout.Nm
                tt.dt_terms && (dt_ao = 1.5 * (Tend * (kk.dah + kout.da) - kout.datau) / el[1] - kout.dM / kout.Nm)
            end
            advance(m_post, i2, dO2, dw2, dt_out)
        catch; out.decayed = true; out.t_decay_s = t_k; out.ended = :invalid; break end
        if kick
            if to.kick_post_drag
                m = add_kick(m, kout.da, kout.de, kout.di, kout.dO, kout.dw)
            else
                m = add_kick(m, kk.da - kk.dah, kk.de - kk.deh, kk.di - kk.dih, kk.dOmega - kk.dOh, kk.domega - kk.dwh)
            end
        end
        t_apo = t_k + dt_out + (to.kick_post_drag ? dt_ao : dtk_arc - dtk_in)
    end
    return out, (ns, D, DT, kick)
end

"""Self-calibrated seed with the v5 predictor (as `calibrate_seed4`)."""
function calibrate_seed5(ctx::EstContext, s::EstSettings, o::V2Opts, tab::SPTable, zf::ZonalField, zo::ZonalOpts, tf::TesseralField, tt::TesseralOpts, to::TimingOpts, inc_rad, raan_rad, omega_rad, et0;
                         seed0=140.0e3, tol_m=1.0, max_iter=8)
    s3 = calibrate_seed3(ctx, s, o, tab, zf, zo, inc_rad, raan_rad, omega_rad, et0; seed0=seed0, tol_m=tol_m, max_iter=max_iter)
    isfinite(s3) || return s3
    q_target = sqrt(ctx.q_low * ctx.q_high); seed = s3
    legacy = to.window === :kepler && !to.post_drag_period && !to.zonal_sp
    _, info = predict_history5(ctx, s, o, tab, zf, zo, tf, tt, to, s3, inc_rad, raan_rad, omega_rad, et0)      # the carry rule sees the full history of the v3 seed
    (!info[4] && legacy) && return s3
    ttf = with_force(tt, info[4])
    for _ in 1:max_iter
        p, _ = predict_history5(ctx, s, o, tab, zf, zo, tf, ttf, to, seed, inc_rad, raan_rad, omega_rad, et0; max_passes=1)
        isempty(p.t_s) && return NaN
        h_target = altitude_of_density_fast(ctx.fit, 2 * q_target * ctx.w2m / p.vrel_mps[1]^3)
        abs(h_target - p.alt_m[1]) < tol_m && return seed
        seed += h_target - p.alt_m[1]
    end
    return seed
end
