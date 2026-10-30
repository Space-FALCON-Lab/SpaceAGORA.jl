# Estimator v3 driver (Julia half): self-tests, Mars and Earth predictions. Serial, single process (SPICE is not thread-safe). One Julia job at a time, memory capped:
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. run_estimator_v3.jl --phase=selftest|mars|earth|mars_attrib|timing [--results-root=DIR] [--config=MARS_TOML]
#        [--force=3,4,5 --name=NAME] [--study=STUDY_TOML] [--budget-km=1.0] [--rule=tail|per_degree] [--out=DIR]
#   selftest   (a) averaged zonal rates vs the disturbing-function (Phi_n) rates, (b) vs Vern9 integration, (c) critical inclinations, acceleration checks. No truth data.
#   mars       v3 predictions for the cells of <root>/mars (truth seeds from truth/<tag>/calibration.csv) -> pred_v3/, pred_v3_selfcal/, pred_cells_v3.csv
#   mars_attrib  gravity-term attribution: predictors vs truth variants zonal/deg3/j2only (needs run_mars_generality.jl --phase=diag) -> mars/attrib_gravity.csv
#   earth      v3 predictions for the 336 Earth cells of <root>/oracle/oracle_cell_labels.csv -> <root>/estimator_v3/{pred,pred_selfcal,pred_cells.csv}
#   (the criterion values D_3, D_4, D_5, D_6, D_8, D_12, D_20 and the carried degrees are columns of pred_cells*.csv)
#   --force    carry exactly these degrees (diagnostic); output folder suffix --name. Default: the criterion decides (budget --budget-km, default 1 km).
using Statistics, OrdinaryDiffEq
include(joinpath(@__DIR__, "mars_common.jl"))
include(joinpath(@__DIR__, "estimator_v3.jl"))

function parse_cli3()
    kv = Dict{String, String}()
    for a in ARGS
        m = match(r"^--([^=]+)=(.*)$", a); m === nothing || (kv[m.captures[1]] = m.captures[2])
    end
    return kv
end
const CLI3 = parse_cli3()
const STUDY3 = load_study(abspath(get(CLI3, "study", STUDY_CONFIG_PATH)))
const ROOT3 = haskey(CLI3, "results-root") ? abspath(CLI3["results-root"]) : results_root(STUDY3)
const CFG23 = TOML.parsefile(joinpath(STUDY_DIR, "estimator_v2_config.toml"))
const OPTS3 = (v = CFG23["v2"]; V2Opts(sp_tide=Symbol(v["sp_tide"]), sp_j2=Symbol(v["sp_j2"]), n_sp=Int(v["n_sp_nodes"])))
const TAB3 = SPTable(OPTS3.n_sp)
const SETT3 = settings_from(CFG23)
const ZO3 = ZonalOpts(budget_km=parse(Float64, get(CLI3, "budget-km", "1.0")), rule=Symbol(get(CLI3, "rule", "tail")),
    force=haskey(CLI3, "force") ? (isempty(CLI3["force"]) ? Int[] : parse.(Int, split(CLI3["force"], ","))) : nothing)
const NAME3 = get(CLI3, "name", "")

pred_df3(p::Prediction) = DataFrame(k=1:length(p.t_s), t_s=p.t_s, alt_m=p.alt_m, radius_m=p.radius_m, lat_rad=p.lat_rad, vrel_mps=p.vrel_mps,
    rho=p.rho, q_w_cm2=p.q_w_cm2, a_m=p.a_m, e=p.e, i_rad=p.i_rad, raan_rad=p.raan_rad, omega_rad=p.omega_rad, h_qhigh_m=p.h_qhigh_m, h_qlow_m=p.h_qlow_m)

"""Earth record with the truth-gravity model of the requested degree (the estimator reads GM, Re, J2 and the zonal J_n from it)."""
function earth_record_deg(study, degree::Int, order::Int)
    rec = earth_record(study); cfg = load_corridor_config(CORRIDOR_CONFIG_PATH)
    harm = GravitationalHarmonicsModel(degree, order, joinpath(REPO_ROOT, String(cfg.values["gravity"]["coefficient_file"])), rec.planet)
    return PlanetRecord(; (f => getfield(rec, f) for f in fieldnames(PlanetRecord) if f != :harm)..., harm=harm)
end

# ---------------------------------------------------------------------------------------------------------------
# self-tests (no truth data)

"""Disturbing-function secular rates (kaula_lib.zonal_rates): (de/dt, domega/dt, dOmega/dt) of the single-averaged zonal degree n with Pi_n = J_n (Re/p)^n signed; quadrature in true anomaly."""
function phi_rates(n::Int, e, inc, w, Pi_n, N; nq=4096)
    nu = [(j - 0.5) * 2π / nq for j in 1:nq]; si, ci = sin(inc), cos(inc)
    Phi = 0.0; Phi_e = 0.0; Phi_i = 0.0; Phi_w = 0.0
    for v in nu
        c = cos(v); u = w + v; x = si * sin(u); P, dP = legendre_pd(n, x); S = (1 + e * c)^(n - 1)
        Phi += S * P; Phi_e += n >= 2 ? (n - 1) * (1 + e * c)^(n - 2) * c * P : 0.0
        Phi_i += S * dP * ci * sin(u); Phi_w += S * dP * si * cos(u)
    end
    Phi, Phi_e, Phi_i, Phi_w = Phi / nq, Phi_e / nq, Phi_i / nq, Phi_w / nq; q = 1 - e^2
    return N * Pi_n * q / e * Phi_w, -N * Pi_n * (q / e * ((2n - 1) * e / q * Phi + Phi_e) - ci / si * Phi_i), -N * Pi_n * Phi_i / si
end

function avg_rates_edwO(F::OrbitFrame, ns, Js, mu, Re, N)
    r = zonal_rates(F, ns, Js, mu, Re, N); e = F.e; k, h = F.e * F.cw, F.e * F.sw
    return (k * r[2] + h * r[3]) / e, (k * r[3] - h * r[2]) / e^2, r[5], r[1], r[4]       # de/dt, domega/dt, dOmega/dt, da/dt, di/dt
end

"""Cowell periapsis-radius series under J2 + the zonal degrees in `Jd` (Dict n -> J_n), sampled at the apoapsis crossings (r.v = 0, downward) up to the time `tend`:
returns `(t_apo, rp)`. Vern9, rtol 1e-13. Sampling at a fixed orbital phase keeps the short-period terms of the high degrees (large near periapsis) out of the comparison."""
function cowell_rp(Jd::Dict{Int, Float64}, mu, Re, a, e, inc, raan, om, tend)
    ns = sort(collect(keys(Jd))); Js = [Jd[n] for n in ns]
    r0, v0 = elements_to_cartesian((a=a, e=e, i=inc, raan=raan, omega=om, nu=π); mu=mu)
    u0 = SVector(r0..., v0...)
    f(u, p, t) = (r = SV3(u[1], u[2], u[3]); acc = -mu * r / norm(r)^3 + zonal_acc(r, ns, Js, mu, Re); SVector(u[4], u[5], u[6], acc[1], acc[2], acc[3]))
    cond(u, t, integ) = u[1] * u[4] + u[2] * u[5] + u[3] * u[6]
    cb = ContinuousCallback(cond, nothing, integ -> nothing; save_positions=(false, true))
    sol = solve(ODEProblem(f, u0, (0.0, tend)), Vern9(); abstol=1e-7, reltol=1e-13, callback=cb, save_everystep=false, save_start=false, save_end=false)
    rp = [(el = osc_from_state(SV3(u[1], u[2], u[3]), SV3(u[4], u[5], u[6]), mu); el[1] * (1 - hypot(el[2], el[3]))) for u in sol.u]
    return sol.t, rp
end

"""Mean-element periapsis radius under J2 secular + zonal rates (the v3 averaged dynamics, no tide, no drag), sampled at `ts`: rk2 steps of half an orbit as predict_history3."""
function mean_rp_series(zf::ZonalField, ns, Js, nz, a, e, inc, raan, om, ts)
    ctx_J2 = zf.J[2]; mu = zf.mu; Re = zf.Re
    # initial MEAN elements from the osculating apoapsis state through the J2 first-order map, exactly as predict_history3 (and v2) does
    mm = osculating_to_mean((a=a, e=e, i=inc, raan=raan, omega=om, nu=π); mu=mu, radius=Re, j2=ctx_J2).elements
    m = SV5(mm.a, mm.e * cos(mm.omega), mm.e * sin(mm.omega), mm.i, mm.raan); res = Float64[]; t = 0.0; T = 2π * sqrt(m[1]^3 / mu)
    function rates(m)
        F = OrbitFrame(m[1], m[2], m[3], m[4], m[5], mu); dOm, dw, _ = j2_rates(F, mu, ctx_J2, Re)
        isempty(ns) && return zero(SV5), dOm, dw
        zinc, wz = zonal_inc_dw(F, ns, Js, mu, Re, nz); return zinc, dOm, dw + wz
    end
    for tt in ts
        while t < tt - 1e-6
            dt = min(0.5 * T, tt - t)
            i1, dO1, dw1 = rates(m); i2, dO2, dw2 = rates(advance(m, i1, dO1, dw1, dt)); m = advance(m, i2, dO2, dw2, dt); t += dt
        end
        push!(res, m[1] * (1 - hypot(m[2], m[3])))
    end
    return res
end

function phase_selftest()
    println("selftest 0: zonal acceleration vs the repo's J2 and the diagnostic J3 forms")
    include(joinpath(@__DIR__, "estimator_v2_j3_diag.jl"))
    mu, Re = 4.2828370371e13, 3.3942e6; worst = 0.0
    for r in (SV3(1.1e7, 2.3e6, 6.0e6), SV3(-4.0e6, 3.1e6, -2.5e6), SV3(3.6e6, 0.4e6, 0.2e6))
        a3 = zonal_acc(r, [3], [3.1575e-5], mu, Re); b3 = Base.invokelatest(j3_acc, r, 3.1575e-5, mu, Re)
        a2 = zonal_acc(r, [2], [1.9586e-3], mu, Re); b2 = j2_acc(r, 1.9586e-3, mu, Re)
        worst = max(worst, norm(a3 - b3) / norm(b3), norm(a2 - b2) / norm(b2))
    end
    println("  max relative difference of the n = 2 and n = 3 accelerations: $worst"); @assert worst < 1e-12
    # numerical gradient of the potential for a high degree
    Phi(r) = (rn = norm(r); s = r[3] / rn; -mu * sum(J * Re^n * legendre_pd(n, s)[1] / rn^(n + 1) for (n, J) in ((4, -1.6e-5), (7, 3.0e-6), (12, -4.0e-7))))
    unit(k) = SV3(ntuple(j -> j == k ? 1.0 : 0.0, 3))
    r = SV3(9.0e6, -3.0e6, 5.0e6); h = 1.0; g = SV3([(Phi(r + h * unit(k)) - Phi(r - h * unit(k))) / (2h) for k in 1:3])
    ga = zonal_acc(r, [4, 7, 12], [-1.6e-5, 3.0e-6, -4.0e-7], mu, Re); println("  gradient of the potential (n = 4, 7, 12): relative difference $(norm(g - ga) / norm(ga))"); @assert norm(g - ga) / norm(ga) < 1e-7

    println("selftest 1 (a): single-averaged Gauss rates of the exact zonal acceleration vs the disturbing-function (Phi_n) rates; Mars-like orbit, unit Pi_n scaling checked per case")
    worst = 0.0; rows = 0
    for (e, inc, w, Om) in ((0.78, deg2rad(25.2), 0.7, 1.0), (0.82, deg2rad(93.5), 2.4, 0.3), (0.6, deg2rad(116.6), 4.5, 2.0), (0.82, deg2rad(40.0), 1.1, 5.0))
        a = 1.6e7; F = OrbitFrame(a, e * cos(w), e * sin(w), inc, Om, mu); p = a * (1 - e^2)
        for n in (2, 3, 4, 5, 6, 8, 12, 20)
            Jn = 1e-5; Pi_n = Jn * (Re / p)^n
            ed, wd, Od, da, di = avg_rates_edwO(F, [n], [Jn], mu, Re, 1024)
            ee, ww, OO = phi_rates(n, e, inc, w, Pi_n, F.n)
            sc = F.n * abs(Pi_n)      # natural rate scale
            d = max(abs(ed - ee), abs(wd - ww), abs(Od - OO)) / sc; worst = max(worst, d); rows += 1
            (n == 3 || n == 5 || n == 20) && println("  e=$e i=$(round(rad2deg(inc), digits=1)) n=$n: de/dt Gauss $(ed / sc) vs Phi $(ee / sc) (units N|Pi_n|); dw $(wd / sc) vs $(ww / sc); dOm $(Od / sc) vs $(OO / sc); da/dt $(da / (F.n * a * abs(Pi_n)))")
            @assert abs(da) < 1e-9 * F.n * a * abs(Pi_n) "secular a-rate should vanish"
        end
    end
    println("  $rows cases, max |difference| / (N |Pi_n|) = $worst"); @assert worst < 5e-6

    println("selftest 2: node-count convergence of the averaged rates (e = 0.82, i = 63.4 and 25.2 deg; difference to 4096 nodes in units of N|Pi_n|)")
    for inc in (deg2rad(25.2), deg2rad(63.4)), n in (3, 6, 12, 20)
        e = 0.82; w = 1.1; F = OrbitFrame(1.6e7, e * cos(w), e * sin(w), inc, 0.5, mu); p = F.a * (1 - e^2); Pi_n = 1e-5 * (Re / p)^n
        ref = avg_rates_edwO(F, [n], [1e-5], mu, Re, 4096)
        out = []
        for N in (64, 128, 256, 12n)
            r = avg_rates_edwO(F, [n], [1e-5], mu, Re, N); push!(out, "N=$N: $(round(max(abs(r[1] - ref[1]), abs(r[2] - ref[2]), abs(r[3] - ref[3])) / (F.n * abs(Pi_n)); sigdigits=2))")
        end
        println("  i=$(round(rad2deg(inc), digits=1)) n=$n: ", join(out, "  "))
    end

    println("selftest 3 (b): J2 + zonal degrees vs Vern9 Cowell integration of the osculating elements (Mars-like, i = 25.2 and 93.5 deg): effect of adding degree n on the periapsis radius (m), averaged theory vs Cowell")
    zfM = zonal_field_file(joinpath(REPO_ROOT, "data/planets/gravity/Mars50c.csv"), 20)
    for inc_deg in (25.2, 93.522), (nadd, T_d) in ((3, 20.0), (4, 60.0), (5, 60.0), (6, 60.0))
        a, e = 1.603e7, 0.78; inc = deg2rad(inc_deg); T = T_d * 86400.0
        base = Dict{Int, Float64}(2 => zfM.J[2]); with = merge(base, Dict(n => zfM.J[n] for n in 3:nadd))
        without = merge(base, Dict(n => zfM.J[n] for n in 3:nadd-1))
        tw, cw = cowell_rp(with, zfM.mu, zfM.Re, a, e, inc, 0.5, 0.7, T); to, cwo = cowell_rp(without, zfM.mu, zfM.Re, a, e, inc, 0.5, 0.7, T)
        ts = tw[1:min(length(tw), length(to))]; cw = cw[1:length(ts)]; cwo = cwo[1:length(ts)]
        # the Cowell series are osculating at the apoapsis crossings; the averaged series is the mean one at the same times: the DIFFERENCE with/without degree n is compared
        nsw = sort(filter(>(2), collect(keys(with)))); nso = sort(filter(>(2), collect(keys(without)))); nzf(ns) = isempty(ns) ? 64 : max(64, 12 * ns[end])
        aw = mean_rp_series(zfM, nsw, [with[n] for n in nsw], nzf(nsw), a, e, inc, 0.5, 0.7, ts)
        ao = mean_rp_series(zfM, nso, [without[n] for n in nso], nzf(nso), a, e, inc, 0.5, 0.7, ts)
        dc = cw .- cwo; da_ = aw .- ao
        println("  i=$inc_deg add degree $nadd, $T_d d ($(length(ts)) apoapsis samples): Cowell effect at the end $(round(dc[end]; digits=1)) m, averaged $(round(da_[end]; digits=1)) m; max over samples |difference| $(round(maximum(abs.(dc .- da_)); digits=1)) m of a max effect $(round(maximum(abs.(dc)); digits=1)) m ($(round(100 * maximum(abs.(dc .- da_)) / maximum(abs.(dc)); digits=1)) %)")
    end

    println("selftest 4 (c): critical inclinations. J3 eccentricity rate and J2 apsidal rate at i = asin(sqrt(4/5)) (63.435 deg), and the J3 rate at 90 deg")
    for inc in (asin(sqrt(0.8)), π - asin(sqrt(0.8)))
        for w in (0.0, 1.0, 2.5)
            a = 1.6e7; e = 0.78; F = OrbitFrame(a, e * cos(w), e * sin(w), inc, 0.4, mu)
            ed, wd, Od, _, _ = avg_rates_edwO(F, [3], [3.0e-5], mu, Re, 512); scl = F.n * 3.0e-5 * (Re / F.p)^3
            dOm, dw, _ = j2_rates(F, mu, 1.96e-3, Re)
            println("  i=$(round(rad2deg(inc), digits=3)) w=$w: |de/dt_J3|/(N|Pi_3|) = $(abs(ed) / scl) ; J2 domega/dt = $dw rad/s")
            @assert abs(ed) / scl < 1e-9 && abs(dw) < 1e-12 * F.n
        end
    end
    F = OrbitFrame(1.6e7, 0.78, 0.0, deg2rad(90.0), 0.4, mu); ed, = avg_rates_edwO(F, [3], [3.0e-5], mu, Re, 512)
    println("  i=90 deg, w=0: |de/dt_J3|/(N|Pi_3|) = $(abs(ed) / (F.n * 3.0e-5 * (Re / F.p)^3)) (Gamma_3 = (3/8)(1-e^2) sin i (5 sin^2 i - 4) = $(0.375 * (1 - 0.78^2)))")
    println("selftest 5: Gamma_max (Julia) vs kaula_lib values; D_n for a Mars-like and an Earth-like cell")
    for (n, e, inc_deg) in ((3, 0.78, 25.2), (4, 0.78, 25.2), (5, 0.82, 93.5), (12, 0.82, 63.4))
        println("  n=$n e=$e i=$inc_deg: Gamma_max = $(gamma_max(n, e, deg2rad(inc_deg), 256, 256)) (nq=4096, nw=256: $(gamma_max(n, e, deg2rad(inc_deg), 4096, 256)))")
    end
    println("selftest ok")
end

"""Zonal field straight from a gravity CSV (selftest only; planets use `zonal_field(record.harm)`)."""
function zonal_field_file(path, L)
    meta = Dict{String, String}()
    for line in eachline(path)
        startswith(line, "#") || break
        if occursin(":", line); kv = split(line[2:end], ":"; limit=2); meta[strip(kv[1])] = strip(kv[2]); end
    end
    df = CSV.read(path, DataFrame; comment="#"); J = zeros(L)
    for n in 2:L; J[n] = -sqrt(2n + 1.0) * df.C[findfirst((df.degree .== n) .& (df.order .== 0))]; end
    return ZonalField(L, J, parse(Float64, meta["reference_radius_m"]), parse(Float64, meta["gm_m3s2"]), NaN)
end

# ---------------------------------------------------------------------------------------------------------------

function cell_row_v3(tag, seed_t, seed_s, p1, p2, info)
    ns, D = info
    d(n) = n - 2 <= length(D) ? D[n - 2] : NaN
    return (tag=tag, seed_truth_m=seed_t, seed_selfcal_m=seed_s, n_pred=length(p1.t_s), decayed_pred=p1.decayed, t_decay_pred_days=p1.t_decay_s / SECONDS_PER_DAY, ended=String(p1.ended),
        n_pred_selfcal=length(p2.t_s), decayed_pred_selfcal=p2.decayed, t_decay_pred_selfcal_days=p2.t_decay_s / SECONDS_PER_DAY,
        n_pred_v1=0, decayed_pred_v1=false, t_decay_pred_v1_days=NaN, wall_s_first_call=NaN,
        carried=join(ns, ";"), D3_km=d(3), D4_km=d(4), D5_km=d(5), D6_km=d(6), D8_km=d(8), D12_km=d(12), D20_km=d(20))
end

function phase_mars()
    mcp = abspath(get(CLI3, "config", MARS_CONFIG_PATH)); mc = load_mars(mcp); ms = mars_setup(mc); out = mars_results_dir(ROOT3, mc)
    ep(i) = parse_epoch(mc["grid"]["epochs_utc"][i]); etm(i) = SimulationModel.ephemerides_time_seconds(ep(i), SpiceEphemeridesModel())
    cells = CSV.read(joinpath(out, "cells.csv"), DataFrame); g = Float64(mc["orbit"]["initial_periapsis_altitude_guess_m"])
    nodes_f = joinpath(out, "estimator_numerics.toml"); S = isfile(nodes_f) ? with_settings(SETT3; n_drag_half=Int(TOML.parsefile(nodes_f)["n_drag_half"])) : SETT3
    zf = zonal_field(ms.harm); sfx = isempty(NAME3) ? "v3" : NAME3
    for d in ("pred_$sfx", "pred_$(sfx)_selfcal"); mkpath(joinpath(out, d)); end
    println("estimator v3 at Mars: zonal field L=$(zf.L), K=$(zf.K), forced=$(ZO3.force), budget=$(ZO3.budget_km) km, drag nodes $(S.n_drag_half)")
    rows = NamedTuple[]; ctxs = Dict{Int, EstContext}(); t0 = time()
    for c in eachrow(cells)
        ei = Int(c.epoch_idx); ctx = get!(ctxs, ei) do; make_context(mars_record(ms, ep(ei)), S); end
        calf = joinpath(out, "truth", c.tag, "calibration.csv"); isfile(calf) || continue
        seed_t = CSV.read(calf, DataFrame).seed_alt_m[end]; inc, raan, om = deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega); et0 = etm(ei)
        p1, info = predict_history3(ctx, S, OPTS3, TAB3, zf, ZO3, seed_t, inc, raan, om, et0)
        seed_s = calibrate_seed3(ctx, S, OPTS3, TAB3, zf, ZO3, inc, raan, om, et0; seed0=g, tol_m=Float64(CFG23["estimator"]["calibration_tol_m"]), max_iter=Int(CFG23["estimator"]["calibration_max_iter"]))
        p2, _ = predict_history3(ctx, S, OPTS3, TAB3, zf, ZO3, seed_s, inc, raan, om, et0)
        CSV.write(joinpath(out, "pred_$sfx", c.tag * ".csv"), pred_df3(p1)); CSV.write(joinpath(out, "pred_$(sfx)_selfcal", c.tag * ".csv"), pred_df3(p2))
        push!(rows, cell_row_v3(c.tag, seed_t, seed_s, p1, p2, info))
    end
    CSV.write(joinpath(out, "pred_cells_$sfx.csv"), DataFrame(rows)); println("mars v3 done: $(length(rows)) cells in $(round(time() - t0; digits=1)) s")
end


"""Attribution of the Mars error to gravity terms: predictors (frozen v2, v3 with J3 forced, v3 by criterion, v3 with every zonal degree up to the variant's degree forced)
against truth variants of the same gravity (full 20x20, zonal 20x0, deg3 3x0, j2only 2x0). Needs truth/<tag>/passes_<variant>.csv (run_mars_generality.jl --phase=diag).
RMS and max altitude error over the common passes, per cell -> <root>/mars/attrib_gravity.csv."""
function phase_mars_attrib()
    mcp = abspath(get(CLI3, "config", MARS_CONFIG_PATH)); mc = load_mars(mcp); ms = mars_setup(mc); out = mars_results_dir(ROOT3, mc)
    ep(i) = parse_epoch(mc["grid"]["epochs_utc"][i]); etm(i) = SimulationModel.ephemerides_time_seconds(ep(i), SpiceEphemeridesModel())
    cells = CSV.read(joinpath(out, "cells.csv"), DataFrame)
    nodes_f = joinpath(out, "estimator_numerics.toml"); S = isfile(nodes_f) ? with_settings(SETT3; n_drag_half=Int(TOML.parsefile(nodes_f)["n_drag_half"])) : SETT3
    zf_full = zonal_field(ms.harm); Lfull = zf_full.L
    variants = (("full", "passes_3b.csv", Lfull), ("zonal", "passes_zonal.csv", Lfull), ("deg3", "passes_deg3.csv", 3), ("j2only", "passes_j2only.csv", 2))
    ctxs = Dict{Int, EstContext}(); rows = NamedTuple[]
    for c in eachrow(cells)
        ei = Int(c.epoch_idx); ctx = get!(ctxs, ei) do; make_context(mars_record(ms, ep(ei)), S); end
        calf = joinpath(out, "truth", c.tag, "calibration.csv"); isfile(calf) || continue
        seed_t = CSV.read(calf, DataFrame).seed_alt_m[end]; inc, raan, om = deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega); et0 = etm(ei)
        for (vname, tfile, L) in variants
            tf = joinpath(out, "truth", c.tag, tfile); isfile(tf) || continue
            tr = CSV.read(tf, DataFrame); zf = ZonalField(L, zf_full.J[1:L], zf_full.Re, zf_full.mu, zf_full.K)
            preds = (("v2", ZonalOpts(force=Int[])), ("v3_J3", ZonalOpts(force=L >= 3 ? [3] : Int[], n_max=L)), ("v3_crit", ZonalOpts(n_max=L)), ("v3a_crit", ZonalOpts(n_max=L, rule=:per_degree)),
                     ("v3_all", ZonalOpts(force=collect(3:L), n_max=L)))
            for (pname, zo) in preds
                p, info = predict_history3(ctx, S, OPTS3, TAB3, zf, zo, seed_t, inc, raan, om, et0)
                n = min(length(p.alt_m), nrow(tr)); e = n == 0 ? Float64[] : p.alt_m[1:n] .- tr.alt_m[1:n]
                push!(rows, (tag=c.tag, inc_deg=c.inc, variant=vname, predictor=pname, carried=join(info[1], ";"), n_common=n, n_truth=nrow(tr), n_pred=length(p.alt_m),
                    truth_t_end_days=tr.time_s[end] / SECONDS_PER_DAY, pred_decayed=p.decayed,
                    rms_km=n == 0 ? NaN : sqrt(mean(e .^ 2)) / 1e3, max_km=n == 0 ? NaN : maximum(abs.(e)) / 1e3, end_err_km=n == 0 ? NaN : e[end] / 1e3))
            end
        end
    end
    CSV.write(joinpath(out, "attrib_gravity.csv"), DataFrame(rows)); println("mars attribution done: $(length(rows)) rows")
end

function phase_earth()
    nbody_planet(); degree = parse(Int, get(CLI3, "degree", "20")); order = parse(Int, get(CLI3, "order", "20"))
    rec = earth_record_deg(STUDY3, degree, order); zf = zonal_field(rec.harm)
    out = haskey(CLI3, "out") ? abspath(CLI3["out"]) : joinpath(ROOT3, isempty(NAME3) ? "estimator_v3" : "estimator_" * NAME3)
    for d in ("pred", "pred_selfcal", "hand"); mkpath(joinpath(out, d)); end
    ctx = make_context(rec, SETT3)
    labels = CSV.read(joinpath(ROOT3, "oracle", "oracle_cell_labels.csv"), DataFrame)
    cell_dir3(tag, inc) = inc == 50.0 ? joinpath(ROOT3, "border", "grid", tag[7:end]) : joinpath(ROOT3, "border2", "grid", tag)
    println("estimator v3 at Earth: zonal field L=$(zf.L), K=$(zf.K), forced=$(ZO3.force), budget=$(ZO3.budget_km) km")
    rows = NamedTuple[]; t0 = time()
    for c in eachrow(labels)
        inc, raan, om = deg2rad(c.inc_deg), deg2rad(c.raan_deg), deg2rad(c.omega_deg); et0 = epoch_et(STUDY3, Int(c.epoch))
        seed_t = CSV.read(joinpath(cell_dir3(c.tag, c.inc_deg), "calibration.csv"), DataFrame).seed_alt_m[end]
        w1 = @elapsed (p1, info) = predict_history3(ctx, SETT3, OPTS3, TAB3, zf, ZO3, seed_t, inc, raan, om, et0)
        seed_s = calibrate_seed3(ctx, SETT3, OPTS3, TAB3, zf, ZO3, inc, raan, om, et0; tol_m=Float64(CFG23["estimator"]["calibration_tol_m"]), max_iter=Int(CFG23["estimator"]["calibration_max_iter"]))
        p2, _ = predict_history3(ctx, SETT3, OPTS3, TAB3, zf, ZO3, seed_s, inc, raan, om, et0)
        CSV.write(joinpath(out, "pred", c.tag * ".csv"), pred_df3(p1)); CSV.write(joinpath(out, "pred_selfcal", c.tag * ".csv"), pred_df3(p2))
        isempty(p1.alt_m) || CSV.write(joinpath(out, "hand", c.tag * ".csv"), DataFrame(hand_history(ctx, SETT3, seed_t, inc, raan, om, et0, p1.alt_m[1], p1.vrel_mps[1])))   # B1 baseline
        r = cell_row_v3(c.tag, seed_t, seed_s, p1, p2, info); push!(rows, merge(r, (wall_s_first_call=w1,)))
    end
    CSV.write(joinpath(out, "pred_cells.csv"), DataFrame(rows)); println("earth v3 done: $(length(rows)) cells in $(round(time() - t0; digits=1)) s -> $out")
end

"""Warm per-prediction cost of v2 against v3 (criterion) on every 14th Earth cell and every 5th Mars cell, 5 repeats, single thread. Run ONLY with the machine otherwise idle
(one Julia job at a time); the criterion's Gamma table is filled by the warm-up call, so its cold cost (first call per (e, i)) is timed separately."""
function phase_timing()
    out = ROOT3; rows = NamedTuple[]
    nbody_planet(); rec = earth_record_deg(STUDY3, 20, 20); zfE = zonal_field(rec.harm); ctxE = make_context(rec, SETT3)
    labels = CSV.read(joinpath(ROOT3, "oracle", "oracle_cell_labels.csv"), DataFrame)
    cdir(tag, inc) = inc == 50.0 ? joinpath(ROOT3, "border", "grid", tag[7:end]) : joinpath(ROOT3, "border2", "grid", tag)
    jobsE = [(ctxE, zfE, deg2rad(c.inc_deg), deg2rad(c.raan_deg), deg2rad(c.omega_deg), epoch_et(STUDY3, Int(c.epoch)), CSV.read(joinpath(cdir(c.tag, c.inc_deg), "calibration.csv"), DataFrame).seed_alt_m[end], "Earth") for c in eachrow(labels[1:14:end, :])]
    mcp = abspath(get(CLI3, "config", MARS_CONFIG_PATH)); mc = load_mars(mcp); ms = mars_setup(mc); zfM = zonal_field(ms.harm); cells = CSV.read(joinpath(mars_results_dir(ROOT3, mc), "cells.csv"), DataFrame)
    ep(i) = parse_epoch(mc["grid"]["epochs_utc"][i]); ctxM = Dict(i => make_context(mars_record(ms, ep(i)), SETT3) for i in 1:length(mc["grid"]["epochs_utc"]))
    jobsM = [(ctxM[Int(c.epoch_idx)], zfM, deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega), SimulationModel.ephemerides_time_seconds(ep(Int(c.epoch_idx)), SpiceEphemeridesModel()),
              CSV.read(joinpath(mars_results_dir(ROOT3, mc), "truth", c.tag, "calibration.csv"), DataFrame).seed_alt_m[end], "Mars") for c in eachrow(cells[1:5:end, :])]
    for (name, f) in (("v2", (j) -> predict_history2(j[1], SETT3, OPTS3, TAB3, j[7], j[3], j[4], j[5], j[6])),
                      ("v3_empty", (j) -> predict_history3(j[1], SETT3, OPTS3, TAB3, j[2], ZonalOpts(force=Int[]), j[7], j[3], j[4], j[5], j[6])),
                      ("v3_tail", (j) -> predict_history3(j[1], SETT3, OPTS3, TAB3, j[2], ZonalOpts(), j[7], j[3], j[4], j[5], j[6])))
        for j in vcat(jobsE, jobsM); f(j); end          # warm-up (compilation, Gamma table)
        for planet in ("Earth", "Mars")
            ts = Float64[]
            for rep in 1:5, j in (planet == "Earth" ? jobsE : jobsM); push!(ts, 1e3 * @elapsed f(j)); end
            push!(rows, (predictor=name, planet=planet, n_predictions=length(ts), median_ms=median(ts), q1_ms=quantile(ts, 0.25), q3_ms=quantile(ts, 0.75)))
        end
    end
    empty!(_GAMMA_CACHE); j = jobsM[1]; tcold = @elapsed predict_history3(j[1], SETT3, OPTS3, TAB3, j[2], ZonalOpts(), j[7], j[3], j[4], j[5], j[6])
    push!(rows, (predictor="v3_tail_cold_gamma_cache", planet="Mars", n_predictions=1, median_ms=1e3 * tcold, q1_ms=NaN, q3_ms=NaN))
    CSV.write(joinpath(out, "estimator_v3_timing.csv"), DataFrame(rows)); println(DataFrame(rows))
end

function main3()
    phase = get(CLI3, "phase", "selftest")
    phase == "selftest" && phase_selftest()
    phase == "mars" && phase_mars()
    phase == "earth" && phase_earth()
    phase == "mars_attrib" && phase_mars_attrib()
    phase == "timing" && phase_timing()
end
main3()
