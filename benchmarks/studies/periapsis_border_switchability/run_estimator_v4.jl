# Estimator v4 driver (Julia half): self-tests, Mars and Earth predictions, timing. Serial, single process (SPICE is not thread-safe). One Julia job at a time, memory capped:
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. run_estimator_v4.jl --phase=selftest|mars|earth|timing|passtime
#        [--results-root=DIR] [--config=MARS_TOML] [--name=NAME] [--lt=N] [--only=tag,tag] [--cells=FILE]
#   selftest   (a) Julia kick = Python reference (tesseral_lib.py), (b) = Lagrange-series kick, (c) single pass vs Cowell integration of the non-zonal acceleration along one arc,
#              (d) coupled stepper vs the truth (= the engine-based numerical predictor at 20x20) on 3 cells over 20 d, (f) no-tesseral limit = v3 bit for bit.
#              needs <root>/v4/reference_kicks.csv (python3 -I v4_reference_kicks.py) and, for (d), the Mars truth under <root>/mars.
#   mars       v4 predictions for the cells of <root>/mars (truth seeds from truth/<tag>/calibration.csv) -> pred_<name>/, pred_<name>_selfcal/, pred_cells_<name>.csv (default name v4)
#   earth      v4 predictions for the Earth cells of <root>/oracle/oracle_cell_labels.csv -> <root>/estimator_<name>/{pred,pred_selfcal,pred_cells.csv}
#   timing     warm cost per 60-day Mars prediction (v3, v4 by field size), single job alone
using Statistics, OrdinaryDiffEq
include(joinpath(@__DIR__, "mars_common.jl"))
include(joinpath(@__DIR__, "estimator_v4.jl"))

function parse_cli4()
    kv = Dict{String, String}()
    for a in ARGS
        m = match(r"^--([^=]+)=(.*)$", a); m === nothing || (kv[m.captures[1]] = m.captures[2])
    end
    return kv
end
const CLI4 = parse_cli4()
const STUDY4 = load_study(abspath(get(CLI4, "study", STUDY_CONFIG_PATH)))
const ROOT4 = haskey(CLI4, "results-root") ? abspath(CLI4["results-root"]) : results_root(STUDY4)
const CFG24 = TOML.parsefile(joinpath(STUDY_DIR, "estimator_v2_config.toml"))
const OPTS4 = (v = CFG24["v2"]; V2Opts(sp_tide=Symbol(v["sp_tide"]), sp_j2=Symbol(v["sp_j2"]), n_sp=Int(v["n_sp_nodes"])))
const TAB4 = SPTable(OPTS4.n_sp)
const SETT4 = settings_from(CFG24)
const ZO4 = ZonalOpts()
const NAME4 = get(CLI4, "name", "v4")
const TT4 = TesseralOpts(Lt=parse(Int, get(CLI4, "lt", "0")), Mt=parse(Int, get(CLI4, "mt", "0")), nodes=parse(Int, get(CLI4, "nodes", "64")),
    force=haskey(CLI4, "force") ? parse(Bool, CLI4["force"]) : nothing)

pred_df4(p::Prediction) = DataFrame(k=1:length(p.t_s), t_s=p.t_s, alt_m=p.alt_m, radius_m=p.radius_m, lat_rad=p.lat_rad, vrel_mps=p.vrel_mps,
    rho=p.rho, q_w_cm2=p.q_w_cm2, a_m=p.a_m, e=p.e, i_rad=p.i_rad, raan_rad=p.raan_rad, omega_rad=p.omega_rad, h_qhigh_m=p.h_qhigh_m, h_qlow_m=p.h_qlow_m)

function earth_record_deg(study, degree::Int, order::Int)
    rec = earth_record(study); cfg = load_corridor_config(CORRIDOR_CONFIG_PATH)
    harm = GravitationalHarmonicsModel(degree, order, joinpath(REPO_ROOT, String(cfg.values["gravity"]["coefficient_file"])), rec.planet)
    return PlanetRecord(; (f => getfield(rec, f) for f in fieldnames(PlanetRecord) if f != :harm)..., harm=harm)
end

function mars_numerics(out)
    nodes_f = joinpath(out, "estimator_numerics.toml")
    return isfile(nodes_f) ? with_settings(SETT4; n_drag_half=Int(TOML.parsefile(nodes_f)["n_drag_half"])) : SETT4
end

# ---------------------------------------------------------------------------------------------------------------
# self-tests

relerr(x, ref, floor) = abs(x - ref) / max(abs(ref), floor)

function selftest_ab(tf::TesseralField)
    R = CSV.read(joinpath(ROOT4, "v4", "reference_kicks.csv"), DataFrame)
    println("selftest (a): Julia kick (halves of 32 nodes) against the Python Gauss quadrature (512 nodes); (b) against the Lagrange series, subset field degree <= 6, order <= 4")
    quants = (:a, :e, :i, :raan, :omega, :rp); jn = (:da, :de, :di, :dOmega, :domega, :drp)
    rows = NamedTuple[]
    for fld in ("full20", "sub6x4")
        sub = R[R.field .== fld, :]; Lt, Mt = fld == "full20" ? (20, 20) : (6, 4)
        # floors: 1% of the largest state value of the same quantity (so a near-zero kick of one state does not dominate)
        floor = Dict(q => 0.01 * maximum(abs.(sub[!, Symbol("quad_", q)])) for q in quants)
        for r in eachrow(sub)
            rot = uniform_rotation(r.theta_p, 7.08823596e-5)
            kk = tesseral_pass_kick((r.a, r.e, deg2rad(r.inc_deg), deg2rad(r.raan_deg), deg2rad(r.omega_deg)), 0.0, tf, rot; nodes=64, Lt=Lt, Mt=Mt)
            e_quad = [relerr(getfield(kk, j), r[Symbol("quad_", q)], floor[q]) for (q, j) in zip(quants, jn)]
            e_half = relerr(kk.drp_half, r.quad_rp_half, 0.01 * maximum(abs.(sub.quad_rp_half)))
            e_ser = fld == "sub6x4" ? [relerr(getfield(kk, j), r[Symbol("series_", q)], floor[q]) for (q, j) in zip(quants, jn)] : fill(NaN, 6)
            push!(rows, (; field=fld, state=r.state, max_quad=maximum(e_quad), half_quad=e_half, max_series=fld == "sub6x4" ? maximum(e_ser) : NaN, rp_J=kk.drp, rp_py=r.quad_rp))
        end
    end
    D = DataFrame(rows)
    for fld in ("full20", "sub6x4")
        s = D[D.field .== fld, :]
        println("  $fld: max over $(nrow(s)) states of the relative difference to the Python quadrature (a, e, i, Omega, omega, rp): $(maximum(s.max_quad)); inbound-half rp: $(maximum(s.half_quad))" *
                (fld == "sub6x4" ? "; to the series (b): $(maximum(s.max_series))" : ""))
    end
    return D
end

"""Cowell integration (point mass + the non-zonal field, uniformly rotating planet) over one apoapsis-to-apoapsis arc starting at the apoapsis of the given osculating
elements; returns the osculating (a, k, h, i, Omega) and time at the periapsis and the second apoapsis events."""
function cowell_arc(tf::TesseralField, a, e, inc, raan, om, thp, w; Lt=tf.L, Mt=Lt)
    mu = tf.mu; T0 = 2π * sqrt(a^3 / mu); kw = kick_work(tf.L)
    r0, v0 = elements_to_cartesian((a=a, e=e, i=inc, raan=raan, omega=om, nu=float(π)); mu=mu)
    u0 = SVector(r0..., v0...)
    function f(u, p, t)
        th = thp + w * (t - T0 / 2); c, s = cos(th), sin(th)
        x, y, z = u[1], u[2], u[3]; rn = sqrt(x * x + y * y + z * z)
        abx, aby, abz = nz_acc(tf, kw, Lt, Mt, c * x + s * y, -s * x + c * y, z)
        return SVector(u[4], u[5], u[6], -mu * x / rn^3 + c * abx - s * aby, -mu * y / rn^3 + s * abx + c * aby, -mu * z / rn^3 + abz)
    end
    ev = NamedTuple[]
    cond(u, t, integ) = u[1] * u[4] + u[2] * u[5] + u[3] * u[6]
    up!(integ) = t_ok(integ.t) && push!(ev, (kind=:peri, t=integ.t, el=osc_from_state(SV3(integ.u[1], integ.u[2], integ.u[3]), SV3(integ.u[4], integ.u[5], integ.u[6]), mu), r=norm(SV3(integ.u[1], integ.u[2], integ.u[3]))))
    down!(integ) = t_ok(integ.t) && (push!(ev, (kind=:apo, t=integ.t, el=osc_from_state(SV3(integ.u[1], integ.u[2], integ.u[3]), SV3(integ.u[4], integ.u[5], integ.u[6]), mu), r=norm(SV3(integ.u[1], integ.u[2], integ.u[3])))); terminate!(integ))
    t_ok(t) = t > 0.05 * T0
    cb = ContinuousCallback(cond, up!, down!; save_positions=(false, false), abstol=1e-9, reltol=0.0)
    solve(ODEProblem(f, u0, (0.0, 1.5 * T0)), Vern9(); abstol=1e-9, reltol=1e-13, callback=cb, save_everystep=false)
    el0 = osc_from_state(SV3(r0...), SV3(v0...), mu)
    return el0, ev, T0
end

function selftest_c(tf::TesseralField)
    R = CSV.read(joinpath(ROOT4, "v4", "reference_kicks.csv"), DataFrame); sub = R[R.field .== "full20", :]
    println("selftest (c): one-arc kick (halves) against Vern9 Cowell integration of the non-zonal acceleration (rtol 1e-13), 6 states; element differences at the periapsis and the second apoapsis event")
    nm = ("a", "e", "i", "Omega", "omega", "rp")
    worst_full = 0.0; worst_half = 0.0; worst_dt = 0.0; rows = NamedTuple[]
    for r in eachrow(sub)[1:6]
        a, e, inc, raan, om, thp, w = r.a, r.e, deg2rad(r.inc_deg), deg2rad(r.raan_deg), deg2rad(r.omega_deg), r.theta_p, 7.08823596e-5
        el0, ev, T0 = cowell_arc(tf, a, e, inc, raan, om, thp, w)
        peri = only(filter(x -> x.kind === :peri, ev)); apo = only(filter(x -> x.kind === :apo, ev))
        kk = tesseral_pass_kick((a, e, inc, raan, om), 0.0, tf, uniform_rotation(thp, w); nodes=64)   # theta(t_peri = 0) = thp, as the Cowell run
        d(ev_) = (ev_.el[1] - el0[1], hypot(ev_.el[2], ev_.el[3]) - hypot(el0[2], el0[3]), ev_.el[4] - el0[4], wrap_pi(ev_.el[5] - el0[5]), wrap_pi(atan(ev_.el[3], ev_.el[2]) - atan(el0[3], el0[2])))
        dc_full = d(apo); dc_half = d(peri)
        jf = (kk.da, kk.de, kk.di, kk.dOmega, kk.domega); jh = (kk.dah, kk.deh, kk.dih, kk.dOh, kk.dwh)
        ef = [abs(jf[q] - dc_full[q]) / max(abs(dc_full[q]), 1e-12) for q in 1:5]
        eh = [abs(jh[q] - dc_half[q]) / max(abs(dc_half[q]), 1e-12) for q in 1:5]
        # r_p kicks (m): full arc from the apoapsis elements, and the actual periapsis radius against the osculating a(1-e) plus the inbound half kick
        rp0 = el0[1] * (1 - hypot(el0[2], el0[3])); rp_full_c = apo.el[1] * (1 - hypot(apo.el[2], apo.el[3])) - rp0
        rp_phys_c = peri.r - rp0
        e_rpf = abs(kk.drp - rp_full_c) / abs(rp_full_c); e_rph = abs(kk.drp_half - rp_phys_c) / abs(rp_phys_c)
        dt_in_c = peri.t - T0 / 2; dt_arc_c = apo.t - T0
        e_dtin = abs(kk.dt_in - dt_in_c); e_dtarc = abs(kk.dt_arc - dt_arc_c)
        println("  state $(r.state): rp full kick $(round(kk.drp; digits=2)) m (Cowell $(round(rp_full_c; digits=2))), physical periapsis radius change $(round(kk.drp_half; digits=2)) m (Cowell $(round(rp_phys_c; digits=2))); ",
                "dt periapsis $(round(kk.dt_in; digits=3)) s (Cowell $(round(dt_in_c; digits=3))), dt apoapsis $(round(kk.dt_arc; digits=3)) s (Cowell $(round(dt_arc_c; digits=3)))")
        push!(rows, (; state=r.state, rp_full_J=kk.drp, rp_full_C=rp_full_c, rp_phys_J=kk.drp_half, rp_phys_C=rp_phys_c, e_rp_full=e_rpf, e_rp_phys=e_rph, dt_in_J=kk.dt_in, dt_in_C=dt_in_c, dt_arc_J=kk.dt_arc, dt_arc_C=dt_arc_c,
            err_a=ef[1], err_e=ef[2], err_i=ef[3], err_O=ef[4], err_w=ef[5], errh_a=eh[1], errh_e=eh[2], errh_i=eh[3], errh_O=eh[4], errh_w=eh[5]))
    end
    D = DataFrame(rows); show(stdout, "text/plain", D); println()
    return D
end

# (d) coupled stepper against the truth (= the engine-based numerical predictor at 20x20: C1 of numerical_predictor_results.md, max difference 0 m)
function selftest_d(tf::TesseralField, ms, mc, S)
    out = mars_results_dir(ROOT4, mc); cells = CSV.read(joinpath(out, "cells.csv"), DataFrame)
    ep(i) = parse_epoch(mc["grid"]["epochs_utc"][i]); etm(i) = SimulationModel.ephemerides_time_seconds(ep(i), SpiceEphemeridesModel())
    chosen = [first(cells[cells.inc .== inc, :]) for inc in (25.2, 63.4, 116.6)]     # rule: the first cell of three inclinations in cells.csv order
    zf = zonal_field(ms.harm); rows = NamedTuple[]
    println("selftest (d): coupled stepper (v4, truth seed) against the truth passes (= engine-based numerical predictor at 20x20, 0 m by construction), 3 cells, first 20 d")
    for c in chosen
        ei = Int(c.epoch_idx); ctx = make_context(mars_record(ms, ep(ei)), S)
        seed_t = CSV.read(joinpath(out, "truth", c.tag, "calibration.csv"), DataFrame).seed_alt_m[end]
        tr = CSV.read(joinpath(out, "truth", c.tag, "passes_3b.csv"), DataFrame); tr = tr[tr.time_s .<= 20 * 86400.0, :]
        inc, raan, om = deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega); et0 = etm(ei)
        S20 = S      # the carry rule sees the full 60-d history (a 20-d horizon would truncate the commensurability sweep and lower D_T); the first 20 d are compared
        p3, _ = predict_history3(ctx, S20, OPTS4, TAB4, zf, ZO4, seed_t, inc, raan, om, et0)
        p4, info = predict_history4(ctx, S20, OPTS4, TAB4, zf, ZO4, tf, TesseralOpts(), seed_t, inc, raan, om, et0)
        for (nm, p) in (("v3", p3), ("v4", p4))
            n = min(sum(p.t_s .<= 20 * 86400.0), nrow(tr)); e = p.alt_m[1:n] .- tr.alt_m[1:n]; dt = p.t_s[1:n] .- tr.time_s[1:n]
            push!(rows, (; tag=c.tag, inc=c.inc, pred=nm, n_pred=length(p.alt_m), n_truth=nrow(tr), DT_km=info[3], carried=info[4], pass1_err_m=e[1],
                rms_km_0_5=sqrt(mean(e[tr.time_s[1:n] .< 5 * 86400] .^ 2)) / 1e3,
                rms_km_5_10=sqrt(mean(e[5 * 86400 .<= tr.time_s[1:n] .< 10 * 86400] .^ 2)) / 1e3,
                rms_km_10_20=sqrt(mean(e[10 * 86400 .<= tr.time_s[1:n]] .^ 2)) / 1e3, rms_dt_s=sqrt(mean(dt .^ 2)), max_dt_s=maximum(abs.(dt))))
        end
    end
    D = DataFrame(rows); show(stdout, "text/plain", D); println()
    return D
end

# (f) the no-tesseral limits: force=false returns the v3 object; a zero-coefficient field through the whole stepper reproduces v3 bit for bit
function selftest_f(tf::TesseralField, ms, mc, S)
    out = mars_results_dir(ROOT4, mc); cells = CSV.read(joinpath(out, "cells.csv"), DataFrame)
    ep(i) = parse_epoch(mc["grid"]["epochs_utc"][i]); etm(i) = SimulationModel.ephemerides_time_seconds(ep(i), SpiceEphemeridesModel())
    zf = zonal_field(ms.harm); ok1 = true; ok2 = true; n = 0; tfix = Float64[]
    for c in eachrow(cells[1:32:end, :])
        ei = Int(c.epoch_idx); ctx = make_context(mars_record(ms, ep(ei)), S)
        seed_t = CSV.read(joinpath(out, "truth", c.tag, "calibration.csv"), DataFrame).seed_alt_m[end]; inc, raan, om = deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega); et0 = etm(ei)
        p3, _ = predict_history3(ctx, S, OPTS4, TAB4, zf, ZO4, seed_t, inc, raan, om, et0)
        pa, _ = predict_history4(ctx, S, OPTS4, TAB4, zf, ZO4, tf, TesseralOpts(force=false), seed_t, inc, raan, om, et0)
        pz, _ = predict_history4(ctx, S, OPTS4, TAB4, zf, ZO4, zero_field(tf), TesseralOpts(force=true, post_drag_period=false), seed_t, inc, raan, om, et0)
        pt, _ = predict_history4(ctx, S, OPTS4, TAB4, zf, ZO4, zero_field(tf), TesseralOpts(force=true), seed_t, inc, raan, om, et0)
        n_ = min(length(p3.t_s), length(pt.t_s)); n_ > 0 && push!(tfix, maximum(abs.(pt.t_s[1:n_] .- p3.t_s[1:n_])))
        same(p, q) = all(isequal(getfield(p, f), getfield(q, f)) for f in fieldnames(Prediction))
        ok1 &= same(p3, pa); ok2 &= same(p3, pz); n += 1
    end
    println("selftest (f): $n Mars cells; force=false identical to v3: $ok1; zero-coefficient field through the full stepper identical to v3 (every field of Prediction, post_drag_period off): $ok2; with the outbound-leg timing fix the pass times differ from v3 by up to $(round.(tfix; digits=0)) s per cell (the fix only)")
    return ok1, ok2
end

function phase_selftest()
    mcp = abspath(get(CLI4, "config", MARS_CONFIG_PATH)); mc = load_mars(mcp); ms = mars_setup(mc); tf = tesseral_field(ms.harm); S = mars_numerics(mars_results_dir(ROOT4, mc))
    which = split(get(CLI4, "which", "a,c,d,f"), ",")
    out = mkpath(joinpath(ROOT4, "v4"))
    "a" in which && (CSV.write(joinpath(out, "selftest_ab.csv"), selftest_ab(tf)))
    "c" in which && (CSV.write(joinpath(out, "selftest_c.csv"), selftest_c(tf)))
    "d" in which && (CSV.write(joinpath(out, "selftest_d.csv"), selftest_d(tf, ms, mc, S)))
    "f" in which && selftest_f(tf, ms, mc, S)
    "n" in which && selftest_nodes(tf)
    println("selftest done")
end

"""Node-count convergence of the Julia kick against its own 512-node value (numerical resolution; used to fix the default `nodes`)."""
function selftest_nodes(tf::TesseralField)
    R = CSV.read(joinpath(ROOT4, "v4", "reference_kicks.csv"), DataFrame); sub = R[R.field .== "full20", :]
    println("node convergence of the kick (relative difference of the r_p kick, the inbound-half kick and dt_arc to the 512-node value; worst of the 10 states)")
    for nn in (16, 24, 32, 48, 64, 96, 128)
        w = 0.0; wh = 0.0; wt = 0.0
        for r in eachrow(sub)
            rot = uniform_rotation(r.theta_p, 7.08823596e-5); el = (r.a, r.e, deg2rad(r.inc_deg), deg2rad(r.raan_deg), deg2rad(r.omega_deg))
            k = tesseral_pass_kick(el, 0.0, tf, rot; nodes=nn); k0 = tesseral_pass_kick(el, 0.0, tf, rot; nodes=512)
            w = max(w, abs(k.drp - k0.drp) / abs(k0.drp)); wh = max(wh, abs(k.drp_half - k0.drp_half) / abs(k0.drp_half)); wt = max(wt, abs(k.dt_arc - k0.dt_arc) / max(abs(k0.dt_arc), 1e-3))
        end
        println("  nodes $nn: rp $w, rp_half $wh, dt_arc $wt")
    end
end

# ---------------------------------------------------------------------------------------------------------------

function cell_row_v4(tag, seed_t, seed_s, p1, p2, info1, info2)
    ns, D = info1
    d(n) = n - 2 <= length(D) ? D[n - 2] : NaN
    return (tag=tag, seed_truth_m=seed_t, seed_selfcal_m=seed_s, n_pred=length(p1.t_s), decayed_pred=p1.decayed, t_decay_pred_days=p1.t_decay_s / SECONDS_PER_DAY, ended=String(p1.ended),
        n_pred_selfcal=length(p2.t_s), decayed_pred_selfcal=p2.decayed, t_decay_pred_selfcal_days=p2.t_decay_s / SECONDS_PER_DAY,
        n_pred_v1=0, decayed_pred_v1=false, t_decay_pred_v1_days=NaN, wall_s_first_call=NaN,
        carried=join(ns, ";"), D3_km=d(3), D4_km=d(4), D5_km=d(5), D6_km=d(6), D8_km=d(8), D12_km=d(12), D20_km=d(20),
        DT_km=info2[3], tesseral=info2[4], DT_km_selfcal=NaN)
end

function phase_mars()
    mcp = abspath(get(CLI4, "config", MARS_CONFIG_PATH)); mc = load_mars(mcp); ms = mars_setup(mc); out = mars_results_dir(ROOT4, mc)
    ep(i) = parse_epoch(mc["grid"]["epochs_utc"][i]); etm(i) = SimulationModel.ephemerides_time_seconds(ep(i), SpiceEphemeridesModel())
    cf = get(CLI4, "cells", "cells.csv"); cells = CSV.read(joinpath(out, cf), DataFrame); g = Float64(mc["orbit"]["initial_periapsis_altitude_guess_m"])
    S = mars_numerics(out); zf = zonal_field(ms.harm); tf = tesseral_field(ms.harm); sfx = NAME4
    for d in ("pred_$sfx", "pred_$(sfx)_selfcal"); mkpath(joinpath(out, d)); end
    println("estimator v4 at Mars: zonal L=$(zf.L), tesseral field L=$(tf.L), kick field Lt=$(TT4.Lt), Mt=$(TT4.Mt), nodes=$(TT4.nodes), force=$(TT4.force), budget=$(TT4.budget_km) km")
    rows = NamedTuple[]; ctxs = Dict{Int, EstContext}(); t0 = time()
    for c in eachrow(cells)
        ei = Int(c.epoch_idx); ctx = get!(ctxs, ei) do; make_context(mars_record(ms, ep(ei)), S); end
        calf = joinpath(out, "truth", c.tag, "calibration.csv"); isfile(calf) || continue
        seed_t = CSV.read(calf, DataFrame).seed_alt_m[end]; inc, raan, om = deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega); et0 = etm(ei)
        tth = haskey(CLI4, "oracle-time") ? CSV.read(joinpath(out, "truth", c.tag, "passes_3b.csv"), DataFrame).time_s : nothing     # UNFROZEN diagnostic (name it v4x_*)
        p1, info1 = predict_history4(ctx, S, OPTS4, TAB4, zf, ZO4, tf, TT4, seed_t, inc, raan, om, et0; theta_times=tth)
        seed_s = calibrate_seed4(ctx, S, OPTS4, TAB4, zf, ZO4, tf, TT4, inc, raan, om, et0; seed0=g, tol_m=Float64(CFG24["estimator"]["calibration_tol_m"]), max_iter=Int(CFG24["estimator"]["calibration_max_iter"]))
        p2, info2 = predict_history4(ctx, S, OPTS4, TAB4, zf, ZO4, tf, TT4, seed_s, inc, raan, om, et0; theta_times=tth)
        CSV.write(joinpath(out, "pred_$sfx", c.tag * ".csv"), pred_df4(p1)); CSV.write(joinpath(out, "pred_$(sfx)_selfcal", c.tag * ".csv"), pred_df4(p2))
        push!(rows, cell_row_v4(c.tag, seed_t, seed_s, p1, p2, info1, info1))
    end
    CSV.write(joinpath(out, "pred_cells_$sfx.csv"), DataFrame(rows)); println("mars v4 done: $(length(rows)) cells in $(round(time() - t0; digits=1)) s (includes compilation and SPICE; not a timing result)")
end

function phase_earth()
    nbody_planet(); degree = parse(Int, get(CLI4, "degree", "20")); order = parse(Int, get(CLI4, "order", "20"))
    rec = earth_record_deg(STUDY4, degree, order); zf = zonal_field(rec.harm); tf = tesseral_field(rec.harm)
    out = haskey(CLI4, "out") ? abspath(CLI4["out"]) : joinpath(ROOT4, "estimator_" * NAME4)
    for d in ("pred", "pred_selfcal"); mkpath(joinpath(out, d)); end
    ctx = make_context(rec, SETT4)
    labels = CSV.read(joinpath(ROOT4, "oracle", "oracle_cell_labels.csv"), DataFrame)
    cell_dir4(tag, inc) = inc == 50.0 ? joinpath(ROOT4, "border", "grid", tag[7:end]) : joinpath(ROOT4, "border2", "grid", tag)
    println("estimator v4 at Earth: zonal L=$(zf.L), tesseral L=$(tf.L), force=$(TT4.force), budget=$(TT4.budget_km) km")
    rows = NamedTuple[]; t0 = time()
    for c in eachrow(labels)
        inc, raan, om = deg2rad(c.inc_deg), deg2rad(c.raan_deg), deg2rad(c.omega_deg); et0 = epoch_et(STUDY4, Int(c.epoch))
        seed_t = CSV.read(joinpath(cell_dir4(c.tag, c.inc_deg), "calibration.csv"), DataFrame).seed_alt_m[end]
        p1, info1 = predict_history4(ctx, SETT4, OPTS4, TAB4, zf, ZO4, tf, TT4, seed_t, inc, raan, om, et0)
        seed_s = calibrate_seed4(ctx, SETT4, OPTS4, TAB4, zf, ZO4, tf, TT4, inc, raan, om, et0; tol_m=Float64(CFG24["estimator"]["calibration_tol_m"]), max_iter=Int(CFG24["estimator"]["calibration_max_iter"]))
        p2, info2 = predict_history4(ctx, SETT4, OPTS4, TAB4, zf, ZO4, tf, TT4, seed_s, inc, raan, om, et0)
        CSV.write(joinpath(out, "pred", c.tag * ".csv"), pred_df4(p1)); CSV.write(joinpath(out, "pred_selfcal", c.tag * ".csv"), pred_df4(p2))
        push!(rows, cell_row_v4(c.tag, seed_t, seed_s, p1, p2, info1, info1))
    end
    CSV.write(joinpath(out, "pred_cells.csv"), DataFrame(rows)); println("earth v4 done: $(length(rows)) cells in $(round(time() - t0; digits=1)) s -> $out")
end

"""Warm per-prediction cost (ms) of v3 and of v4 by kick field, on every 5th Mars cell, 5 repeats, single thread, truth seed, 60-day horizon. Run ONLY with the machine
otherwise idle (one Julia job at a time). The v4 times include the carry decision (v3 history + D_T); `decision_only` is v3 + D_T without the stepper."""
function phase_timing()
    mcp = abspath(get(CLI4, "config", MARS_CONFIG_PATH)); mc = load_mars(mcp); ms = mars_setup(mc); out = mars_results_dir(ROOT4, mc); S = mars_numerics(out)
    ep(i) = parse_epoch(mc["grid"]["epochs_utc"][i]); etm(i) = SimulationModel.ephemerides_time_seconds(ep(i), SpiceEphemeridesModel())
    zf = zonal_field(ms.harm); tf = tesseral_field(ms.harm); cells = CSV.read(joinpath(out, "cells.csv"), DataFrame)
    ctxs = Dict(i => make_context(mars_record(ms, ep(i)), S) for i in 1:length(mc["grid"]["epochs_utc"]))
    jobs = [(ctxs[Int(c.epoch_idx)], deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega), etm(Int(c.epoch_idx)),
             CSV.read(joinpath(out, "truth", c.tag, "calibration.csv"), DataFrame).seed_alt_m[end]) for c in eachrow(cells[1:parse(Int, get(CLI4, "stride", "5")):end, :])]
    rows = NamedTuple[]
    variants = Any[("v3", j -> predict_history3(j[1], S, OPTS4, TAB4, zf, ZO4, j[6], j[2], j[3], j[4], j[5])),
        ("decision_only(v3+D_T)", j -> (p = predict_history3(j[1], S, OPTS4, TAB4, zf, ZO4, j[6], j[2], j[3], j[4], j[5])[1]; dt_of_prediction(tf, p, norm(j[1].omega_e), 96)))]
    for Lt in (6, 8, 12, 20)
        tt = TesseralOpts(Lt=Lt, force=true)
        push!(variants, ("v4_kick_only_L$Lt(force)", j -> predict_history4(j[1], S, OPTS4, TAB4, zf, ZO4, tf, tt, j[6], j[2], j[3], j[4], j[5])))
        tt2 = TesseralOpts(Lt=Lt)
        push!(variants, ("v4_total_L$Lt(rule)", j -> predict_history4(j[1], S, OPTS4, TAB4, zf, ZO4, tf, tt2, j[6], j[2], j[3], j[4], j[5])))
    end
    for (name, f) in variants
        for j in jobs; f(j); end        # warm-up
        ts = Float64[]
        for rep in 1:5, j in jobs; push!(ts, 1e3 * @elapsed f(j)); end
        push!(rows, (variant=name, n_predictions=length(ts), median_ms=median(ts), q1_ms=quantile(ts, 0.25), q3_ms=quantile(ts, 0.75), mean_ms=mean(ts)))
        println(rows[end])
    end
    CSV.write(joinpath(ROOT4, "v4", get(CLI4, "timing-file", "estimator_v4_timing.csv")), DataFrame(rows))
end

function main4()
    phase = get(CLI4, "phase", "selftest")
    phase == "selftest" && phase_selftest()
    phase == "mars" && phase_mars()
    phase == "earth" && phase_earth()
    phase == "timing" && phase_timing()
end
main4()
