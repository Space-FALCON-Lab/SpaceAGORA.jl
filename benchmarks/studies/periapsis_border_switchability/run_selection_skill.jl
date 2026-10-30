# Task A2 driver (Julia half): estimator v1 predictions for the 336 truth cells, attribution-cell checks, timing, candidate scan.
# Serial, single process (SPICE is not thread-safe). One Julia job at a time, memory capped:
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. run_selection_skill.jl --phase=selftest|cells|attrib|timing|scan [--results-root=DIR] [--out=DIR]
# Reads: <root>/border, <root>/border2 (truth), <root>/oracle/oracle_cell_labels.csv (cell list), <root>/border3 (attribution cells).
# Writes under --out (default <root>/selection). The Python half (analyze_selection_skill.py) applies the A1 policies and computes metrics.
include(joinpath(@__DIR__, "estimator_v1.jl"))
using Statistics

function parse_cli()
    kv = Dict{String, String}()
    for a in ARGS
        m = match(r"^--([^=]+)=(.*)$", a); m === nothing || (kv[m.captures[1]] = m.captures[2])
    end
    return kv
end
const CLI = parse_cli()

const SEL_CONFIG = TOML.parsefile(get(CLI, "config", joinpath(STUDY_DIR, "selection_skill_config.toml")))
const STUDY = load_study()
const ROOT = haskey(CLI, "results-root") ? abspath(CLI["results-root"]) : results_root(STUDY)
const OUT = haskey(CLI, "out") ? abspath(CLI["out"]) : joinpath(ROOT, "selection")
const SETT = settings_from(SEL_CONFIG)
old_dir(tag) = tag[7:end]    # I0500_R000_W000_E1 -> R000_W000_E1

cell_dir(tag, inc) = inc == 50.0 ? joinpath(ROOT, "border", "grid", old_dir(tag)) : joinpath(ROOT, "border2", "grid", tag)

pred_df(p::Prediction) = DataFrame(k=1:length(p.t_s), t_s=p.t_s, alt_m=p.alt_m, radius_m=p.radius_m, lat_rad=p.lat_rad, vrel_mps=p.vrel_mps,
    rho=p.rho, q_w_cm2=p.q_w_cm2, a_m=p.a_m, e=p.e, i_rad=p.i_rad, raan_rad=p.raan_rad, omega_rad=p.omega_rad, h_qhigh_m=p.h_qhigh_m, h_qlow_m=p.h_qlow_m)

# ---------------------------------------------------------------------------------------------------------------

function phase_selftest(ctx)
    println("selftest: closed-form inverse of the density fit vs bisection (`altitude_of_density`)")
    worst = 0.0
    for rho in exp.(range(log(1e-10), log(2e-6), length=400))
        worst = max(worst, abs(altitude_of_density_fast(ctx.fit, rho) - altitude_of_density(ctx.fit, rho)))
    end
    println("  max |difference| = $(worst) m"); @assert worst < 1e-3
    println("selftest: nonsingular Gauss equations vs the classical (a, e, i, Omega, omega) form")
    F = OrbitFrame(3.7e7, 0.8 * cos(1.1), 0.8 * sin(1.1), 0.5, 2.0, ctx.mu); worstc = 0.0
    for E in (0.3, 1.7, 3.0, 4.4)
        r, _, rm, cnu, snu, _ = node_state(F, E); fR, fS, fW = 1e-6, -2e-6, 3e-6; p, hm, e = F.p, F.hm, F.e; u = F.w + atan(snu, cnu)
        da = 2 * F.a^2 / hm * (e * snu * fR + p / rm * fS); de = (p * snu * fR + ((p + rm) * cnu + rm * e) * fS) / hm
        dw = (-p * cnu * fR + (p + rm) * snu * fS) / (hm * e) - rm * sin(u) * F.cosi / (hm * F.sini) * fW
        g = gauss_g(F, rm, cnu, snu, fR, fS, fW)
        dk = cos(F.w) * de - e * sin(F.w) * dw; dh = sin(F.w) * de + e * cos(F.w) * dw
        worstc = max(worstc, abs(g[1] - da) / abs(da), abs(g[2] - dk) / abs(dk), abs(g[3] - dh) / abs(dh))
    end
    println("  max relative difference = $worstc"); @assert worstc < 1e-12
    println("selftest: drag integral vs the textbook circular-orbit result (constant density, no rotation, whole orbit)")
    ctx0 = EstContext(ctx.mu, ctx.Re, ctx.J2, ctx.Rp_e, ctx.planet, zero(SV3), ctx.fit, ctx.beta, ctx.q_low, ctx.q_high, ctx.q_hard, ctx.w2m, ctx.apo_alt_m,
        ctx.body_names, ctx.body_mu, ctx.gl_x, ctx.gl_w)
    ec = 1e-4; a = 7.0e6; rho0 = 1e-11
    din, dout = drag_deltas((a=a, e=ec, i=0.5, raan=1.0, omega=0.7), ctx0, EstSettings(drag_window_m=1e12, n_drag_half=ctx.gl_x |> length), r -> rho0)
    da_num = (din + dout)[1]; da_ref = -2π * a^2 * rho0 / ctx.beta
    println("  da per orbit: numeric $da_num m, -2 pi a^2 rho / beta = $da_ref m, ratio $(da_num / da_ref)"); @assert abs(da_num / da_ref - 1) < 2e-3
    println("selftest: first pass of a truth cell vs the estimator (seed from calibration.csv)")
    tag = "I0234_R000_W090_E1"; cd_ = cell_dir(tag, 23.4)
    seed = CSV.read(joinpath(cd_, "calibration.csv"), DataFrame).seed_alt_m[end]; tr = CSV.read(joinpath(cd_, "passes_3b.csv"), DataFrame)
    p = predict_history(ctx, SETT, seed, deg2rad(23.4), 0.0, deg2rad(90.0), epoch_et(STUDY, 1))
    println("  pass 1: predicted alt $(p.alt_m[1]) m, truth $(tr.alt_m[1]) m; t $(p.t_s[1]) vs $(tr.time_s[1]) s; pass $(length(p.t_s)): alt $(p.alt_m[end]) vs $(tr.alt_m[min(end, length(p.t_s))])")
    println("selftest ok")
end

# ---------------------------------------------------------------------------------------------------------------

function phase_cells(ctx)
    labels = CSV.read(joinpath(ROOT, "oracle", "oracle_cell_labels.csv"), DataFrame)
    for d in ("pred", "pred_selfcal", "hand"); mkpath(joinpath(OUT, d)); end
    rows = NamedTuple[]
    t0 = time()
    for (ci, c) in enumerate(eachrow(labels))
        inc, raan, om = deg2rad(c.inc_deg), deg2rad(c.raan_deg), deg2rad(c.omega_deg)
        et0 = epoch_et(STUDY, Int(c.epoch))
        seed_t = CSV.read(joinpath(cell_dir(c.tag, c.inc_deg), "calibration.csv"), DataFrame).seed_alt_m[end]
        w1 = @elapsed p1 = predict_history(ctx, SETT, seed_t, inc, raan, om, et0)
        seed_s = calibrate_seed(ctx, SETT, inc, raan, om, et0; tol_m=Float64(SEL_CONFIG["estimator"]["calibration_tol_m"]), max_iter=Int(SEL_CONFIG["estimator"]["calibration_max_iter"]))
        p2 = predict_history(ctx, SETT, seed_s, inc, raan, om, et0)
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

function phase_attrib(ctx)
    series = DataFrame(); jtab = DataFrame()
    variants = [("v0_drag_on", SETT, "v0"), ("v5_drag_off", with_settings(SETT; drag=false), "v5"),
                ("v5_drag_off_apoapsis_epoch", with_settings(SETT; body_epoch=:apoapsis, drag=false), "v5"),
                ("v3_no_bodies_drag_on", with_settings(SETT; tides=false), "v3")]
    for c in attribution_cells()
        inc, raan, om = deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega); et0 = epoch_et(STUDY, c.epoch)
        seed = CSV.read(joinpath(cell_dir(c.tag, c.inc), "calibration.csv"), DataFrame).seed_alt_m[end]
        for (name, s, tv) in variants
            tr = CSV.read(joinpath(ROOT, "border3", "runs", c.tag, "passes_$tv.csv"), DataFrame)
            p = predict_history(ctx, s, seed, inc, raan, om, et0)
            n = min(length(p.t_s), nrow(tr))
            append!(series, DataFrame(tag=c.tag, variant=name, k=1:n, t_truth_s=tr.time_s[1:n], t_pred_s=p.t_s[1:n], rp_truth_m=tr.radius_m[1:n], rp_pred_m=p.radius_m[1:n],
                alt_truth_m=tr.alt_m[1:n], alt_pred_m=p.alt_m[1:n], a_truth_m=tr.a_m[1:n], a_pred_m=p.a_m[1:n], e_truth=tr.e[1:n], e_pred=p.e[1:n],
                raan_truth=tr.raan_rad[1:n], raan_pred=p.raan_rad[1:n], om_truth=tr.omega_rad[1:n], om_pred=p.omega_rad[1:n],
                n_truth=nrow(tr), n_pred=length(p.t_s), decayed_pred=p.decayed))
        end
        # J2 secular investigation on the V3 (no Sun, no Moon) truth: static osculating rates, static mean rates, estimator with drag on, tides off
        tr = CSV.read(joinpath(ROOT, "border3", "runs", c.tag, "passes_v3.csv"), DataFrame); T = tr.time_s[end] - tr.time_s[1]
        rate_of(a, e, i) = (F = OrbitFrame(a, e, 0.0, i, 0.0, ctx.mu); j2_rates(F, ctx.mu, ctx.J2, ctx.Re)[1:2])
        osc_r = rate_of(tr.a_m[1], tr.e[1], tr.i_rad[1])
        mm = osculating_to_mean((a=tr.a_m[1], e=tr.e[1], i=tr.i_rad[1], raan=tr.raan_rad[1], omega=tr.omega_rad[1], nu=0.0); mu=ctx.mu, radius=ctx.Re, j2=ctx.J2).elements
        mean_r = rate_of(mm.a, mm.e, mm.i)
        pj = predict_history(ctx, with_settings(SETT; tides=false), seed, inc, raan, om, et0); nj = min(length(pj.t_s), nrow(tr))
        push!(jtab, (tag=c.tag, days=T / SECONDS_PER_DAY, truth_dOmega_deg=rad2deg(tr.raan_rad[end] - tr.raan_rad[1]), truth_domega_deg=rad2deg(tr.omega_rad[end] - tr.omega_rad[1]),
            static_osc_dOmega_deg=rad2deg(osc_r[1] * T), static_osc_domega_deg=rad2deg(osc_r[2] * T),
            static_mean_dOmega_deg=rad2deg(mean_r[1] * T), static_mean_domega_deg=rad2deg(mean_r[2] * T),
            est_dOmega_deg=rad2deg(pj.raan_rad[nj] - pj.raan_rad[1]), est_domega_deg=rad2deg(pj.omega_rad[nj] - pj.omega_rad[1]),
            truth_a_end_over_start=tr.a_m[end] / tr.a_m[1], est_a_end_over_start=pj.a_m[nj] / pj.a_m[1]))
    end
    mkpath(OUT)
    CSV.write(joinpath(OUT, "attrib_series.csv"), series); CSV.write(joinpath(OUT, "j2_investigation_raw.csv"), jtab)
    println("attrib done")
end

# ---------------------------------------------------------------------------------------------------------------

"""Numerical convergence of the three resolution settings on the 12 attribution cells: max |periapsis altitude difference| to a refined reference."""
function phase_convergence(ctx)
    ref = with_settings(SETT; n_tide=256, n_drag_half=128, drag_window_m=400.0e3)
    ctx_ref = make_context(STUDY, ref)
    variants = vcat([("n_tide=$n", with_settings(SETT; n_tide=n), ctx) for n in (8, 16, 32, 64, 128)],
                    [("n_drag_half=$n", with_settings(SETT; n_drag_half=n), make_context(STUDY, with_settings(SETT; n_drag_half=n))) for n in (8, 16, 32, 64, 128)],
                    [("drag_window_km=$(w)", with_settings(SETT; drag_window_m=1000.0 * w), ctx) for w in (100, 150, 250, 400, 600)])
    rows = NamedTuple[]
    for c in attribution_cells()
        inc, raan, om = deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega); et0 = epoch_et(STUDY, c.epoch)
        seed = CSV.read(joinpath(cell_dir(c.tag, c.inc), "calibration.csv"), DataFrame).seed_alt_m[end]
        pr = predict_history(ctx_ref, ref, seed, inc, raan, om, et0)
        for (name, s_, cx) in variants
            p = predict_history(cx, s_, seed, inc, raan, om, et0); n = min(length(p.alt_m), length(pr.alt_m))
            push!(rows, (tag=c.tag, variant=name, max_alt_diff_m=maximum(abs.(p.alt_m[1:n] .- pr.alt_m[1:n])), n_passes=length(p.alt_m), n_ref=length(pr.alt_m)))
        end
    end
    D = DataFrame(rows); CSV.write(joinpath(OUT, "convergence_raw.csv"), D)
    println(combine(groupby(D, :variant), :max_alt_diff_m => maximum => :worst_over_12_cells_m))
end

# ---------------------------------------------------------------------------------------------------------------

function phase_timing(ctx)
    labels = CSV.read(joinpath(ROOT, "oracle", "oracle_cell_labels.csv"), DataFrame)
    sel = labels[round.(Int, range(1, nrow(labels), length=Int(SEL_CONFIG["timing"]["n_cells"]))), :]
    runs = [(deg2rad(c.inc_deg), deg2rad(c.raan_deg), deg2rad(c.omega_deg), epoch_et(STUDY, Int(c.epoch)),
             CSV.read(joinpath(cell_dir(c.tag, c.inc_deg), "calibration.csv"), DataFrame).seed_alt_m[end]) for c in eachrow(sel)]
    one(r) = predict_history(ctx, SETT, r[5], r[1], r[2], r[3], r[4])
    one(runs[1]); one(runs[1])   # warm-up (compilation)
    reps = Int(SEL_CONFIG["timing"]["repeats"]); rows = NamedTuple[]
    for (ri, r) in enumerate(runs), rep in 1:reps
        GC.gc(false)
        t = @elapsed p = one(r); al = @allocated one(r)
        push!(rows, (cell=sel.tag[ri], rep=rep, wall_s=t, alloc_bytes=al, passes=length(p.t_s), decayed=p.decayed))
    end
    D = DataFrame(rows); CSV.write(joinpath(OUT, "timing_candidates.csv"), D)
    ts = Float64[]
    println("warm wall per 60-day candidate: median $(median(D.wall_s) * 1e3) ms, min $(minimum(D.wall_s) * 1e3) ms, max $(maximum(D.wall_s) * 1e3) ms; alloc median $(median(D.alloc_bytes) / 1e6) MB; Threads=$(Threads.nthreads())")
    # self-calibration cost
    c = runs[1]; calibrate_seed(ctx, SETT, c[1], c[2], c[3], c[4]); tc = @elapsed calibrate_seed(ctx, SETT, c[1], c[2], c[3], c[4])
    println("warm seed calibration: $(tc * 1e3) ms")
    open(joinpath(OUT, "timing_summary.toml"), "w") do io
        println(io, "n_runs = $(nrow(D))\nmedian_wall_ms = $(median(D.wall_s) * 1e3)\nmin_wall_ms = $(minimum(D.wall_s) * 1e3)\nmax_wall_ms = $(maximum(D.wall_s) * 1e3)")
        println(io, "median_alloc_mb = $(median(D.alloc_bytes) / 1e6)\nmedian_passes = $(median(D.passes))\ncalibration_wall_ms = $(tc * 1e3)\njulia_threads = $(Threads.nthreads())")
    end
end

# ---------------------------------------------------------------------------------------------------------------

function phase_scan(ctx)
    sc = SEL_CONFIG["scan"]
    incs = Float64.(sc["inclination_deg"]); raans = collect(Float64(sc["raan_start_deg"]):Float64(sc["raan_step_deg"]):Float64(sc["raan_stop_deg"]))
    oms = collect(Float64(sc["omega_start_deg"]):Float64(sc["omega_step_deg"]):Float64(sc["omega_stop_deg"])); eps_ = Int.(sc["epochs"])
    mkpath(OUT)
    ets = Dict(e => epoch_et(STUDY, e) for e in eps_)
    summ = NamedTuple[]; hist = DataFrame()
    n = length(incs) * length(raans) * length(oms) * length(eps_); println("scan: $n candidates")
    t0 = time(); cnt = 0
    for e in eps_, inc in incs, raan in raans, om in oms
        i_, r_, o_ = deg2rad(inc), deg2rad(raan), deg2rad(om)
        seed = calibrate_seed(ctx, SETT, i_, r_, o_, ets[e]; tol_m=Float64(SEL_CONFIG["estimator"]["calibration_tol_m"]), max_iter=Int(SEL_CONFIG["estimator"]["calibration_max_iter"]))
        p = predict_history(ctx, SETT, seed, i_, r_, o_, ets[e])
        push!(summ, (inc_deg=inc, raan_deg=raan, omega_deg=om, epoch=e, seed_m=seed, n_pred=length(p.t_s), decayed_pred=p.decayed, t_decay_days=p.t_decay_s / SECONDS_PER_DAY))
        append!(hist, DataFrame(inc_deg=inc, raan_deg=raan, omega_deg=om, epoch=e, k=1:length(p.t_s), t_s=p.t_s, alt_m=p.alt_m, q_w_cm2=p.q_w_cm2, a_m=p.a_m, e=p.e, h_qhigh_m=p.h_qhigh_m, h_qlow_m=p.h_qlow_m))
        cnt += 1; cnt % 500 == 0 && println("  $cnt/$n, $(round(time() - t0, digits=1)) s")
    end
    wall = time() - t0
    CSV.write(joinpath(OUT, "scan_summary.csv"), DataFrame(summ)); CSV.write(joinpath(OUT, "scan_hist.csv.gz"), hist; compress=true)
    open(joinpath(OUT, "scan_timing.toml"), "w") do io
        println(io, "candidates = $n\nwall_s_estimator_and_calibration_and_collection = $wall\nwall_ms_per_candidate = $(1e3 * wall / n)")
    end
    println("scan done: $n candidates in $(round(wall, digits=1)) s")
end

# ---------------------------------------------------------------------------------------------------------------

function main()
    phase = get(CLI, "phase", "selftest"); mkpath(OUT)
    nbody_planet()    # furnish the SPICE kernels
    ctx = make_context(STUDY, SETT)
    println("phase=$phase root=$ROOT out=$OUT mu=$(ctx.mu) Re=$(ctx.Re) J2=$(ctx.J2) beta=$(ctx.beta)")
    phase == "selftest" && phase_selftest(ctx)
    phase == "cells" && phase_cells(ctx)
    phase == "attrib" && phase_attrib(ctx)
    phase == "timing" && phase_timing(ctx)
    phase == "convergence" && phase_convergence(ctx)
    phase == "scan" && phase_scan(ctx)
end
main()
