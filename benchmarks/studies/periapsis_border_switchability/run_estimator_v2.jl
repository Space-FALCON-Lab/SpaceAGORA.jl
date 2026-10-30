# Estimator v2 driver (Julia half): self-tests of the short-period terms against independent numerical integrations, the 336-cell predictions,
# attribution-cell series, the closed-form (Pi-form) oscillation predictor, the planet-generality evaluation, timing.
# Serial, single process (SPICE is not thread-safe). One Julia job at a time, memory capped:
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. run_estimator_v2.jl --phase=selftest|cells|attrib|osc|planets|timing [--variant=NAME] [--results-root=DIR]
# Reads <root>/border, <root>/border2, <root>/border3, <root>/oracle/oracle_cell_labels.csv; writes <root>/estimator_v2[_VARIANT]/.
# The Python half (analyze_selection_skill.py with --out, then analyze_estimator_v2.py) scores the written predictions.
using Statistics, OrdinaryDiffEq
include(joinpath(@__DIR__, "shortperiod_analytic.jl"))

function parse_cli()
    kv = Dict{String, String}()
    for a in ARGS
        m = match(r"^--([^=]+)=(.*)$", a); m === nothing || (kv[m.captures[1]] = m.captures[2])
    end
    return kv
end
const CLI = parse_cli()
const CFG2 = TOML.parsefile(get(CLI, "config", joinpath(STUDY_DIR, "estimator_v2_config.toml")))
const STUDY = load_study()
const ROOT = haskey(CLI, "results-root") ? abspath(CLI["results-root"]) : results_root(STUDY)
const VARIANT = get(CLI, "variant", "")
const OUT = haskey(CLI, "out") ? abspath(CLI["out"]) : joinpath(ROOT, isempty(VARIANT) ? "estimator_v2" : "estimator_v2_" * VARIANT)
const SETT = settings_from(CFG2)
const OPTS = begin
    v = CFG2["v2"]; n = Int(v["n_sp_nodes"])
    if VARIANT == "quadrupole"; V2Opts(sp_tide=:quadrupole, sp_j2=Symbol(v["sp_j2"]), n_sp=n)
    elseif VARIANT == "j2quad"; V2Opts(sp_tide=Symbol(v["sp_tide"]), sp_j2=:quad, n_sp=n)
    elseif VARIANT == "none"; V2Opts(sp_tide=:none, sp_j2=:map, n_sp=n)
    else V2Opts(sp_tide=Symbol(v["sp_tide"]), sp_j2=Symbol(v["sp_j2"]), n_sp=n) end
end
const TAB = SPTable(OPTS.n_sp)

cell_dir(tag, inc) = inc == 50.0 ? joinpath(ROOT, "border", "grid", tag[7:end]) : joinpath(ROOT, "border2", "grid", tag)

pred_df(p::Prediction) = DataFrame(k=1:length(p.t_s), t_s=p.t_s, alt_m=p.alt_m, radius_m=p.radius_m, lat_rad=p.lat_rad, vrel_mps=p.vrel_mps,
    rho=p.rho, q_w_cm2=p.q_w_cm2, a_m=p.a_m, e=p.e, i_rad=p.i_rad, raan_rad=p.raan_rad, omega_rad=p.omega_rad, h_qhigh_m=p.h_qhigh_m, h_qlow_m=p.h_qlow_m)

with_bodies(ctx::EstContext, names) = EstContext(ctx.mu, ctx.Re, ctx.J2, ctx.Rp_e, ctx.planet, ctx.omega_e, ctx.fit, ctx.beta, ctx.q_low, ctx.q_high, ctx.q_hard, ctx.w2m,
    ctx.apo_alt_m, String.(names), [MU_BODY[String(b)] for b in names], ctx.gl_x, ctx.gl_w, ctx.center, ctx.rot)

pv2(ctx, s, seed, inc, raan, om, et0; opts=OPTS, tab=TAB, kw...) = predict_history2(ctx, s, opts, tab, seed, inc, raan, om, et0; kw...)

"""Initial mean state of the v2 predictor (same code path as `predict_history2`): returns (m, mm) = ((a, k, h, i, Omega) mean vector, J2-mean element tuple)."""
function initial_mean(ctx, s, o, tab, seed, inc, raan, om, et0)
    init = initial_osculating(ctx, seed, inc, raan, om)
    ivec = SV5(init.a, init.e * cos(init.omega), init.e * sin(init.omega), init.i, init.raan)
    rbs0 = s.tides ? body_positions(ctx, et0) : SV3[]
    d0 = zero(SV5); m = ivec; mm = init
    for _ in 1:(o.sp_tide === :none ? 1 : 3)
        el = vec_el(ivec - d0)
        mm = osculating_to_mean((a=el.a, e=el.e, i=el.i, raan=el.raan, omega=el.omega, nu=init.nu); mu=ctx.mu, radius=ctx.Re, j2=ctx.J2).elements
        m = mean_state(mm, ctx); d0 = tide_sp(m, ctx, o, tab, rbs0, :apo)
    end
    return m, mm, d0
end

# ---------------------------------------------------------------------------------------------------------------
# self-tests: short-period terms against independent numerical integration (no repo engine, no truth data)

function integrate_osc(u0, accel, mu, ts)
    f(u, p, t) = (r = SV3(u[1], u[2], u[3]); v = SV3(u[4], u[5], u[6]); a = -mu * r / norm(r)^3 + accel(r); SVector(v[1], v[2], v[3], a[1], a[2], a[3]))
    sol = solve(ODEProblem(f, u0, (0.0, ts[end])), Vern9(); abstol=1e-7, reltol=1e-13, saveat=ts)
    return [osc_from_state(SV3(u[1], u[2], u[3]), SV3(u[4], u[5], u[6]), mu) for u in sol.u]
end

"""Model osculating (a, k, h, i, Omega) at the times `ts`: mean propagated with the secular rates (J2 + averaged tide) plus short periods (J2 map or quadrature, tide quadrature or none)."""
function model_osc(ctx, ts, osc0::SV5, E0, rbs, kind::Symbol; j2sp::Symbol, tidesp::Bool, j2_on::Bool, tides_on::Bool, tab=TAB)
    accj = r -> j2_acc(r, ctx.J2, ctx.mu, ctx.Re); acct = tide_accfun(rbs, ctx.body_mu)
    sp(F, E) = (j2sp === :quad ? sp_quadrature_general(F, accj, tab, E) : zero(SV5)) + (tidesp ? sp_quadrature_general(F, acct, tab, E) : zero(SV5))
    F0 = OrbitFrame(osc0[1], osc0[2], osc0[3], osc0[4], osc0[5], ctx.mu)
    if j2sp === :map
        e0 = F0.e; f0 = 2atan(sqrt(1 + e0) * sin(E0 / 2), sqrt(1 - e0) * cos(E0 / 2))
        mm = osculating_to_mean((a=F0.a, e=e0, i=F0.i, raan=mod(osc0[5], 2π), omega=mod(F0.w, 2π), nu=f0); mu=ctx.mu, radius=ctx.Re, j2=ctx.J2).elements
        m0 = mean_state(mm, ctx) - (tidesp ? sp_quadrature_general(F0, acct, tab, E0) : zero(SV5))
    else
        m0 = osc0 - sp(F0, E0)
    end
    Fm = OrbitFrame(m0[1], m0[2], m0[3], m0[4], m0[5], ctx.mu)
    M0 = E0 - Fm.e * sin(E0); ndot = j2_on ? mean_motion_rate(m0, ctx) : Fm.n
    inc, dOm, dw = mean_rates_r(m0, ctx, EstSettings(tides=tides_on, j2=j2_on), rbs)
    res = SV5[]
    for t in ts
        m = advance(m0, inc, dOm, dw, t)
        Mt = M0 + ndot * t; Et = Mt; Fe = OrbitFrame(m[1], m[2], m[3], m[4], m[5], ctx.mu)
        for _ in 1:60; Et -= (Et - Fe.e * sin(Et) - Mt) / (1 - Fe.e * cos(Et)); end
        if j2sp === :map
            f = 2atan(sqrt(1 + Fe.e) * sin(Et / 2), sqrt(1 - Fe.e) * cos(Et / 2))
            o = mean_to_osculating((a=Fe.a, e=Fe.e, i=Fe.i, raan=mod(m[5], 2π), omega=mod(Fe.w, 2π), nu=f); mu=ctx.mu, radius=ctx.Re, j2=ctx.J2).elements
            v = SV5(o.a, o.e * cos(o.omega), o.e * sin(o.omega), o.i, o.raan) + (tidesp ? sp_quadrature_general(Fe, acct, tab, Et) : zero(SV5))
        else
            v = m + sp(Fe, Et)
        end
        push!(res, v)
    end
    return res
end

function phase_selftest(ctx)
    println("selftest 1: closed-form quadrupole short period vs quadrature (quadrupole force, frozen perturber)")
    rng = (state = UInt64(20261016); () -> (state = state * 6364136223846793005 + 1442695040888963407; Float64(state >> 11) / 2.0^53))
    worst = 0.0
    for e in (0.3, 0.6, 0.823), trial in 1:3
        u = normalize(SV3(rng() - 0.5, rng() - 0.5, rng() - 0.5)); F = OrbitFrame(3.7e7, e * cos(1.0 + trial), e * sin(1.0 + trial), 0.4 + 0.5 * trial, 0.7, ctx.mu)
        Pi = 1e-5; K = Pi * ctx.mu / F.a^3
        accq = r -> K * (3 * u * dot(u, r) - r)
        for which in (:peri, :apo)
            d = sp_quadrature(F, accq, TAB, which); ca, ce = quad_sp_coeffs(F.e, dot(u, F.P), dot(u, F.Q), which)
            de = F.cw * d[2] + F.sw * d[3]
            worst = max(worst, abs(d[1] - F.a * Pi * ca) / abs(F.a * Pi * ca), abs(de - Pi * ce) / abs(Pi * ce))
            drp_q = (1 - F.e) * d[1] - F.a * de; drp_c = F.a * Pi * quad_drp_coeff(F.e, dot(u, F.P), dot(u, F.Q))
            which === :peri && (worst = max(worst, abs(drp_q - drp_c) / max(abs(drp_c), 1e-3 * F.a * Pi)))
        end
    end
    println("  max relative difference (a, e, r_p) = $worst"); @assert worst < 1e-8
    println("selftest 2: harmonic coefficients vs rp_rate_nd")
    wh = 0.0
    for (e, i, w, lam, aob) in ((0.82, 1.1, 0.4, 2.0, 0.095), (0.5, 0.3, 2.5, 5.0, 0.01), (0.9, 2.0, 4.0, 1.0, 0.05))
        hh = harmonic_coeffs(e, i, w); cw, sw = cos(w), sin(w); ci, si = cos(i), sin(i)
        P = SV3(cw, sw * ci, sw * si); Wv = SV3(0.0, -si, ci); Q = cross(Wv, P); b = SV3(cos(lam), sin(lam), 0.0)
        wh = max(wh, abs(rp_rate_nd(e, dot(b, P), dot(b, Q), aob) - rate_from_harmonics(e, i, w, lam, aob)))
    end
    println("  max |difference| = $wh"); @assert wh < 1e-12
    println("selftest 3: short-period terms vs numerical integration of the osculating elements (Vern9, rtol 1e-13); state a = 3.69e7 m, e = 0.823, i = 63.4 deg")
    rbs = [SV3(0.6, 0.7, 0.39) / norm(SV3(0.6, 0.7, 0.39)) * 3.844e8, SV3(-0.8, 0.5, 0.33) / norm(SV3(-0.8, 0.5, 0.33)) * 1.496e11]
    mub = [MU_BODY["Moon"], MU_BODY["Sun"]]
    ctx2 = with_bodies(ctx, ["Moon", "Sun"])
    for E0 in (0.0, 2.0, π), (inc, om) in ((deg2rad(63.4), 1.5), (deg2rad(23.4), 0.3))
        F0 = OrbitFrame(3.69e7, 0.823 * cos(om), 0.823 * sin(om), inc, 1.0, ctx.mu)
        r0, v0, _ = node_state(F0, E0); u0 = SVector(r0[1], r0[2], r0[3], v0[1], v0[2], v0[3])
        T = 2π / F0.n; ts = collect(range(0.0, 1.3T, length=40))
        osc0 = osc_from_state(r0, v0, ctx.mu)
        accj = r -> j2_acc(r, ctx.J2, ctx.mu, ctx.Re); acct = tide_accfun(rbs, mub)
        for (name, accel, jon, ton, kinds) in (("J2 only", accj, true, false, [(:none, false, "none"), (:map, false, "J2 map"), (:quad, false, "J2 quadrature")]),
                                     ("tide only (Moon+Sun frozen)", acct, false, true, [(:none, false, "none"), (:none, true, "tide quadrature")]),
                                     ("J2 + tide", r -> accj(r) + acct(r), true, true, [(:none, false, "none"), (:map, true, "J2 map + tide quad"), (:quad, true, "J2 quad + tide quad")])),
                (j2sp, tidesp, lab) in kinds
            truth = integrate_osc(u0, accel, ctx.mu, ts)
            mdl = model_osc(ctx2, ts, osc0, E0, rbs, :x; j2sp=j2sp, tidesp=tidesp, j2_on=jon, tides_on=ton)
            rp(v) = v[1] * (1 - hypot(v[2], v[3]))
            err_rp = maximum(abs(rp(a) - rp(b)) for (a, b) in zip(truth, mdl)); err_a = maximum(abs(a[1] - b[1]) for (a, b) in zip(truth, mdl))
            ptp_rp = maximum(rp.(truth)) - minimum(rp.(truth))
            println("  E0=$(round(E0, digits=2)) i=$(round(rad2deg(inc), digits=1)) $(name), model=$(lab): max |r_p error| = $(round(err_rp, digits=1)) m (r_p range $(round(ptp_rp, digits=1)) m), max |a error| = $(round(err_a, digits=1)) m")
        end
    end
    println("selftest 4: quadrature resolution (delta r_p at periapsis, m; tide = Moon at 3.844e8 m, J2 at the same orbit)")
    F = OrbitFrame(3.69e7, 0.823 * cos(1.5), 0.823 * sin(1.5), deg2rad(63.4), 1.0, ctx.mu)
    for N in (32, 64, 128, 256)
        tb = SPTable(N)
        dt = sp_quadrature(F, tide_accfun(rbs[1:1], mub[1:1]), tb, :peri); dj = sp_quadrature(F, r -> j2_acc(r, ctx.J2, ctx.mu, ctx.Re), tb, :peri)
        rpd(d) = (1 - F.e) * d[1] - F.a * (F.cw * d[2] + F.sw * d[3])
        println("  N=$N: tide $(rpd(dt)) m, J2 $(rpd(dj)) m")
    end
    println("selftest 5: J2 short period at periapsis, repo map vs quadrature (osculating minus mean, m; r_p and components)")
    for (e, inc, om) in ((0.823, 63.4, 1.5), (0.823, 23.4, 0.0), (0.823, 116.6, 3.0), (0.5, 50.0, 2.0), (0.2, 50.0, 2.0))
        a = 3.69e7; m = SV5(a, e * cos(om), e * sin(om), deg2rad(inc), 1.0); Fm = OrbitFrame(a, m[2], m[3], m[4], m[5], ctx.mu)
        dmap = osc_vec(m, ctx) - m; dq = sp_quadrature(Fm, r -> j2_acc(r, ctx.J2, ctx.mu, ctx.Re), TAB, :peri)
        rpd(d) = (1 - Fm.e) * d[1] - a * (Fm.cw * d[2] + Fm.sw * d[3])
        println("  e=$e i=$inc w=$om: r_p map $(round(rpd(dmap), digits=2)) quad $(round(rpd(dq), digits=2)); a map $(round(dmap[1], digits=2)) quad $(round(dq[1], digits=2)); di map $(dmap[4]) quad $(dq[4]); dOm map $(dmap[5]) quad $(dq[5])")
    end
    println("selftest ok")
end

# ---------------------------------------------------------------------------------------------------------------

function phase_cells(ctx)
    labels = CSV.read(joinpath(ROOT, "oracle", "oracle_cell_labels.csv"), DataFrame)
    for d in ("pred", "pred_selfcal", "hand"); mkpath(joinpath(OUT, d)); end
    rows = NamedTuple[]; first_pass = NamedTuple[]
    t0 = time()
    for (ci, c) in enumerate(eachrow(labels))
        inc, raan, om = deg2rad(c.inc_deg), deg2rad(c.raan_deg), deg2rad(c.omega_deg)
        et0 = epoch_et(STUDY, Int(c.epoch))
        seed_t = CSV.read(joinpath(cell_dir(c.tag, c.inc_deg), "calibration.csv"), DataFrame).seed_alt_m[end]
        w1 = @elapsed p1 = pv2(ctx, SETT, seed_t, inc, raan, om, et0)
        seed_s = calibrate_seed2(ctx, SETT, OPTS, TAB, inc, raan, om, et0; tol_m=Float64(CFG2["estimator"]["calibration_tol_m"]), max_iter=Int(CFG2["estimator"]["calibration_max_iter"]))
        p2 = pv2(ctx, SETT, seed_s, inc, raan, om, et0)
        CSV.write(joinpath(OUT, "pred", c.tag * ".csv"), pred_df(p1)); CSV.write(joinpath(OUT, "pred_selfcal", c.tag * ".csv"), pred_df(p2))
        h = hand_history(ctx, SETT, seed_t, inc, raan, om, et0, p1.alt_m[1], p1.vrel_mps[1])
        CSV.write(joinpath(OUT, "hand", c.tag * ".csv"), DataFrame(h))
        push!(rows, (tag=c.tag, seed_truth_m=seed_t, seed_selfcal_m=seed_s, n_pred=length(p1.t_s), decayed_pred=p1.decayed, t_decay_pred_days=p1.t_decay_s / SECONDS_PER_DAY,
            ended=String(p1.ended), n_pred_selfcal=length(p2.t_s), decayed_pred_selfcal=p2.decayed, t_decay_pred_selfcal_days=p2.t_decay_s / SECONDS_PER_DAY, wall_s_first_call=w1))
        ci % 56 == 0 && println("  $ci/$(nrow(labels)) cells, $(round(time() - t0, digits=1)) s")
    end
    CSV.write(joinpath(OUT, "pred_cells.csv"), DataFrame(rows))
    println("cells done in $(round(time() - t0, digits=1)) s")
end

# ---------------------------------------------------------------------------------------------------------------

function attribution_cells()
    cfg = TOML.parsefile(joinpath(STUDY_DIR, "attribution_config.toml"))["cells"]["list"]
    return [(inc=Float64(c[1]), raan=Float64(c[2]), omega=Float64(c[3]), epoch=Int(c[4]), src=String(c[5]),
             tag="I$(lpad(round(Int, 10 * c[1]), 4, '0'))_R$(lpad(round(Int, c[2]), 3, '0'))_W$(lpad(round(Int, c[3]), 3, '0'))_E$(Int(c[4]))") for c in cfg]
end

"""v1 and v2 series against the attribution-cell truth variants (v0 Sun+Moon drag on, v5 drag off, v3 neither, v1 Moon, v2 Sun, v8 Moon drag off, v9 Sun drag off)."""
function phase_attrib(ctx)
    ctx_m = with_bodies(ctx, ["Moon"]); ctx_s = with_bodies(ctx, ["Sun"])
    # (name, truth file tag, settings, context); the first four reproduce the variant names of run_selection_skill.jl (for the shared analysis script)
    variants = [("v0_drag_on", "v0", SETT, ctx), ("v5_drag_off", "v5", with_settings(SETT; drag=false), ctx),
                ("v5_drag_off_apoapsis_epoch", "v5", with_settings(SETT; body_epoch=:apoapsis, drag=false), ctx),
                ("v3_no_bodies_drag_on", "v3", with_settings(SETT; tides=false), ctx),
                ("v1_moon_drag_on", "v1", SETT, ctx_m), ("v2_sun_drag_on", "v2", SETT, ctx_s),
                ("v8_moon_drag_off", "v8", with_settings(SETT; drag=false), ctx_m), ("v9_sun_drag_off", "v9", with_settings(SETT; drag=false), ctx_s)]
    series = DataFrame()
    for c in attribution_cells()
        inc, raan, om = deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega); et0 = epoch_et(STUDY, c.epoch)
        seed = CSV.read(joinpath(cell_dir(c.tag, c.inc), "calibration.csv"), DataFrame).seed_alt_m[end]
        for (name, tv, s, cx) in variants, est in ("v1", "v2")
            tr = CSV.read(joinpath(ROOT, "border3", "runs", c.tag, "passes_$tv.csv"), DataFrame)
            p = est == "v1" ? predict_history(cx, s, seed, inc, raan, om, et0) : predict_history2(cx, s, OPTS, TAB, seed, inc, raan, om, et0)
            n = min(length(p.t_s), nrow(tr))
            append!(series, DataFrame(tag=c.tag, variant=name, est=est, k=1:n, t_truth_s=tr.time_s[1:n], t_pred_s=p.t_s[1:n], rp_truth_m=tr.radius_m[1:n], rp_pred_m=p.radius_m[1:n],
                alt_truth_m=tr.alt_m[1:n], alt_pred_m=p.alt_m[1:n], a_truth_m=tr.a_m[1:n], a_pred_m=p.a_m[1:n], e_truth=tr.e[1:n], e_pred=p.e[1:n],
                raan_truth=tr.raan_rad[1:n], raan_pred=p.raan_rad[1:n], om_truth=tr.omega_rad[1:n], om_pred=p.omega_rad[1:n],
                n_truth=nrow(tr), n_pred=length(p.t_s), decayed_pred=p.decayed))
        end
    end
    mkpath(OUT)
    CSV.write(joinpath(OUT, "attrib_series_all.csv"), series)
    CSV.write(joinpath(OUT, "attrib_series.csv"), select(series[(series.est .== "v2") .& (series.variant .∈ Ref(["v0_drag_on", "v5_drag_off", "v5_drag_off_apoapsis_epoch", "v3_no_bodies_drag_on"])), :], Not(:est)))
    println("attrib done")
end

# ---------------------------------------------------------------------------------------------------------------

"""Closed-form (Pi-form) predictor on the 12 attribution cells: series for Moon only, Sun only and both (drag-off truth V8, V9, V5) and the oscillation parameters."""
function phase_osc(ctx)
    ctx_m = with_bodies(ctx, ["Moon"]); ctx_s = with_bodies(ctx, ["Sun"])
    s_off = with_settings(SETT; drag=false)
    series = DataFrame(); harm = DataFrame()
    for c in attribution_cells()
        inc, raan, om = deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega); et0 = epoch_et(STUDY, c.epoch)
        seed = CSV.read(joinpath(cell_dir(c.tag, c.inc), "calibration.csv"), DataFrame).seed_alt_m[end]
        m0, mm, _ = initial_mean(ctx, s_off, OPTS, TAB, seed, inc, raan, om, et0)
        for (set, tv, bodies, cx) in (("moon", "v8", ["Moon"], ctx_m), ("sun", "v9", ["Sun"], ctx_s), ("both", "v5", ["Moon", "Sun"], ctx))
            tr = CSV.read(joinpath(ROOT, "border3", "runs", c.tag, "passes_$tv.csv"), DataFrame)
            p = predict_history2(cx, s_off, OPTS, TAB, seed, inc, raan, om, et0)
            n = min(length(p.t_s), nrow(tr)); tp = p.t_s[1:n]
            an = analytic_series(ctx, m0, et0, tp, bodies)
            append!(series, DataFrame(tag=c.tag, set=set, k=1:n, t_truth_s=tr.time_s[1:n], t_pred_s=tp, rp_truth_m=tr.radius_m[1:n], rp_v2_m=p.radius_m[1:n],
                rbar_m=an.rbar, sp_m=an.sp, a_m=tr.a_m[1:n], e_truth=tr.e[1:n]))
            if set != "both"
                op = oscillation_params(ctx, m0, et0, tp[1], tp[end], bodies[1])
                row = merge((tag=c.tag, body=bodies[1], i_deg=c.inc, raan_deg=c.raan, omega_deg=c.omega, epoch=c.epoch), NamedTuple(op))
                append!(harm, DataFrame([row]))
            end
        end
    end
    mkpath(OUT)
    CSV.write(joinpath(OUT, "osc_series.csv"), series); CSV.write(joinpath(OUT, "osc_harmonics.csv"), harm)
    println("osc done")
end

# ---------------------------------------------------------------------------------------------------------------

function phase_planets(ctx)
    P = CFG2["planets"]; au = Float64(P["au_m"]); rows = NamedTuple[]
    sun_mu = MU_BODY["Sun"]
    function evaluate(planet_name, planet, perts, orbits)
        for (a, alt_p) in orbits
            rp = planet.Rp_e + alt_p; e = 1 - rp / a; pm = a * (1 - e^2); n = sqrt(planet.μ / a^3)
            PiJ2 = planet.J2 * (planet.Rp_e / pm)^2
            # maximum of the quadrupole 2-lambda harmonic amplitude over the geometry (omega', i'), and of the periapsis short period over the perturber direction
            R2max = 0.0
            for i in range(0, π, length=91), w in range(0, 2π, length=181)
                h = harmonic_coeffs(e, i, w); R2max = max(R2max, hypot(h.c2, h.s2))
            end
            spmax = maximum(abs(quad_drp_coeff(e, r_ * cos(t), r_ * sin(t))) for t in range(0, 2π, length=361), r_ in 0.0:0.05:1.0)
            # J2 short period of r_p at periapsis by quadrature for i = 60 deg, omega = 90 deg (illustrative), in km
            F = OrbitFrame(a, e * cos(π / 2), e * sin(π / 2), deg2rad(60), 0.0, planet.μ)
            dj = sp_quadrature(F, r -> j2_acc(r, planet.J2, planet.μ, planet.Rp_e), TAB, :peri); drp_j2 = (1 - e) * dj[1] - a * (F.cw * dj[2] + F.sw * dj[3])
            for (pname, mu_b, ab, n_b, kind) in perts
                Pi = mu_b / planet.μ * (a / ab)^3
                nu = n_b / n
                push!(rows, (planet=planet_name, perturber=pname, a_km=a / 1e3, e=e, rp_alt_km=alt_p / 1e3, Pi_J2=PiJ2, Pi_3B=Pi, Gamma=Pi / PiJ2, a_over_ab=a / ab, nu=nu,
                    A2max_km=a * Pi * R2max / (2nu) / 1e3, R2max=R2max, sp3b_max_km=a * Pi * spmax / 1e3, sp_j2_km=drp_j2 / 1e3, octupole_A1_scale_km=a * Pi * (a / ab) / nu / 1e3,
                    regime=a / ab < 0.2 ? "outer perturber, expansion valid" : "a/ab >= 0.2: expansion not valid"))
            end
        end
    end
    earth = make_no_gram_planet(:earth); mars = make_no_gram_planet(:mars); venus = make_no_gram_planet(:venus)
    # n_b: Sun about the planet = sqrt(GM_sun / d^3) (planet's heliocentric mean motion, planet mass neglected); satellites sqrt((mu_cb + mu_s) / a^3)
    d_mars = Float64(P["mars_sun_distance_au"]) * au; d_venus = Float64(P["venus_sun_distance_au"]) * au
    d_earth_moon = 3.844e8; n_moon = 2π / (27.321661 * SECONDS_PER_DAY)
    evaluate("Earth", earth, [("Moon", MU_BODY["Moon"], d_earth_moon, n_moon, :outer), ("Sun", sun_mu, au, sqrt(sun_mu / au^3), :outer)],
             [(a, alt) for (a, alt) in ((36.7e6, 140.0e3), (20.0e6, 140.0e3), (10.0e6, 140.0e3))])
    evaluate("Mars", mars, [("Sun", sun_mu, d_mars, sqrt(sun_mu / d_mars^3), :outer)], [(Float64(o[1]), Float64(o[2])) for o in P["mars_orbits"]])
    evaluate("Venus", venus, [("Sun", sun_mu, d_venus, sqrt(sun_mu / d_venus^3), :outer)], [(Float64(o[1]), Float64(o[2])) for o in P["venus_orbits"]])
    # Mars satellites as perturbers (outer perturber form; valid only if a << a_moon)
    for (nm, gm, am) in (("Phobos", Float64(P["phobos_gm_m3s2"]), Float64(P["phobos_a_m"])), ("Deimos", Float64(P["deimos_gm_m3s2"]), Float64(P["deimos_a_m"])))
        evaluate("Mars", mars, [(nm, gm, am, sqrt((mars.μ + gm) / am^3), :outer)], [(Float64(o[1]), Float64(o[2])) for o in P["mars_orbits"]])
    end
    mkpath(OUT); CSV.write(joinpath(OUT, "planets_table.csv"), DataFrame(rows))
    println(DataFrame(rows))
end

# ---------------------------------------------------------------------------------------------------------------

function phase_timing(ctx)
    labels = CSV.read(joinpath(ROOT, "oracle", "oracle_cell_labels.csv"), DataFrame)
    sel = labels[round.(Int, range(1, nrow(labels), length=Int(CFG2["timing"]["n_cells"]))), :]
    runs = [(deg2rad(c.inc_deg), deg2rad(c.raan_deg), deg2rad(c.omega_deg), epoch_et(STUDY, Int(c.epoch)),
             CSV.read(joinpath(cell_dir(c.tag, c.inc_deg), "calibration.csv"), DataFrame).seed_alt_m[end]) for c in eachrow(sel)]
    o_exact = V2Opts(); o_quad = V2Opts(sp_tide=:quadrupole); o_jq = V2Opts(sp_j2=:quad)
    fns = Dict("v1" => r -> predict_history(ctx, SETT, r[5], r[1], r[2], r[3], r[4]),
               "v2_exact" => r -> predict_history2(ctx, SETT, o_exact, TAB, r[5], r[1], r[2], r[3], r[4]),
               "v2_quadrupole" => r -> predict_history2(ctx, SETT, o_quad, TAB, r[5], r[1], r[2], r[3], r[4]),
               "v2_j2quad" => r -> predict_history2(ctx, SETT, o_jq, TAB, r[5], r[1], r[2], r[3], r[4]))
    names = ["v1", "v2_exact", "v2_quadrupole", "v2_j2quad"]
    for nm in names; fns[nm](runs[1]); fns[nm](runs[1]); end   # warm-up (compilation)
    reps = Int(CFG2["timing"]["repeats"]); rows = NamedTuple[]
    for (ri, r) in enumerate(runs), rep in 1:reps, nm in names
        GC.gc(false)
        t = @elapsed p = fns[nm](r); al = @allocated fns[nm](r)
        push!(rows, (cell=sel.tag[ri], rep=rep, method=nm, wall_s=t, alloc_bytes=al, passes=length(p.t_s), decayed=p.decayed))
    end
    D = DataFrame(rows); mkpath(OUT); CSV.write(joinpath(OUT, "timing_candidates.csv"), D)
    open(joinpath(OUT, "timing_summary.toml"), "w") do io
        for nm in names
            g = D[D.method .== nm, :]
            println(io, "[$nm]\nn_runs = $(nrow(g))\nmedian_wall_ms = $(median(g.wall_s) * 1e3)\nmin_wall_ms = $(minimum(g.wall_s) * 1e3)\nmax_wall_ms = $(maximum(g.wall_s) * 1e3)\nmedian_alloc_mb = $(median(g.alloc_bytes) / 1e6)\nmedian_passes = $(median(g.passes))\n")
            println("$nm: median $(median(g.wall_s) * 1e3) ms (min $(minimum(g.wall_s) * 1e3), max $(maximum(g.wall_s) * 1e3)), alloc median $(median(g.alloc_bytes) / 1e6) MB")
        end
        println(io, "julia_threads = $(Threads.nthreads())")
    end
    # calibration cost with v2
    c = runs[1]; calibrate_seed2(ctx, SETT, OPTS, TAB, c[1], c[2], c[3], c[4]); tc = @elapsed calibrate_seed2(ctx, SETT, OPTS, TAB, c[1], c[2], c[3], c[4])
    println("warm v2 seed calibration: $(tc * 1e3) ms")
end

# ---------------------------------------------------------------------------------------------------------------

function phase_blind(ctx)
    B = CSV.read(joinpath(ROOT, "blind", "blind_cells.csv"), DataFrame); mkpath(joinpath(OUT, "blind_pred")); rows = NamedTuple[]
    for c in eachrow(B)
        inc, raan, om = deg2rad(c.inc_deg), deg2rad(c.raan_deg), deg2rad(c.omega_deg); et0 = epoch_et(STUDY, Int(c.epoch))
        seed = calibrate_seed2(ctx, SETT, OPTS, TAB, inc, raan, om, et0; tol_m=Float64(CFG2["estimator"]["calibration_tol_m"]), max_iter=Int(CFG2["estimator"]["calibration_max_iter"]))
        p = predict_history2(ctx, SETT, OPTS, TAB, seed, inc, raan, om, et0)
        CSV.write(joinpath(OUT, "blind_pred", c.tag * ".csv"), pred_df(p))
        push!(rows, (tag=c.tag, seed_m=seed, n_pred=length(p.t_s), decayed_pred=p.decayed, t_decay_pred_days=p.t_decay_s / SECONDS_PER_DAY))
    end
    CSV.write(joinpath(OUT, "blind_pred_cells.csv"), DataFrame(rows)); println("blind done: $(length(rows)) cells")
end

# ---------------------------------------------------------------------------------------------------------------

function main()
    phase = get(CLI, "phase", "selftest"); mkpath(OUT)
    nbody_planet()    # furnish the SPICE kernels
    ctx = make_context(STUDY, SETT)
    println("phase=$phase variant=$(isempty(VARIANT) ? "headline" : VARIANT) opts=$OPTS root=$ROOT out=$OUT")
    phase == "selftest" && phase_selftest(ctx)
    phase == "cells" && phase_cells(ctx)
    phase == "attrib" && phase_attrib(ctx)
    phase == "osc" && phase_osc(ctx)
    phase == "planets" && phase_planets(ctx)
    phase == "timing" && phase_timing(ctx)
    phase == "blind" && phase_blind(ctx)
end
main()
