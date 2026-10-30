# Stage I: the paper corridor (band 0.05 to 0.25 W/cm2, target 0.15). Reuses the Stage F, G and H case functions unchanged (they read the band from the config path in the seed, so a config with the
# paper band is all they need) and adds only (1) the calibration case, which calibrates the first-pass Q to the paper target instead of the geometric mean of the band, and (2) a dispatcher so one
# run_monte_carlo campaign can mix the fixed, Stage F adaptive, Stage G asymmetric and Stage H arms. No src/ or examples/ change; no Stage F, G or H file is edited.
include(joinpath(@__DIR__, "stagef_common.jl"))
include(joinpath(@__DIR__, "stageg_common.jl"))
include(joinpath(@__DIR__, "stageh_common.jl"))

"""View of the flat per-planet config for the calibration helpers: gram_calibrate and pi_rule_altitude take their target Q as sqrt(q_low * q_high) (Stage F's rule), so q_low is replaced by
q_target^2 / q_high, which makes that geometric mean exactly the paper target q_target. Used ONLY by the calibration case; the flights read the real band."""
function stagei_cal_view(sc)
    v = deepcopy(sc); b = v["band"]; qt = Float64(b["q_target_w_cm2"]); b["q_low_w_cm2"] = qt^2 / Float64(b["q_high_w_cm2"])
    return v
end

"""Calibration case (a copy of stagef_cal_task with the calibration helpers given stagei_cal_view); writes <seed.out>/<tag>/{guess.csv, nodes.csv, calibration.csv}. Plain-data row; never throws."""
function stagei_cal_task(seed)
    t0 = time()
    try
        meta = stagef_meta(seed.config_path); sc = load_stagef(seed.planet; path=seed.config_path); scc = stagei_cal_view(sc); ms = stagef_setup(sc)
        epoch = parse_epoch(sc["grid"]["epochs_utc"][seed.epoch_idx]); rot = mars_rot(ms.planet, epoch)
        et0 = SimulationModel.ephemerides_time_seconds(epoch, SpiceEphemeridesModel()); stagef_env!(false, meta)
        model = stagef_gram_model(sc, epoch, stagef_seed(seed.planet * seed.tag)); pargs = gram_pargs(ms.planet, Float64(sc["atmosphere_gram"]["entry_interface_km"]))
        g = pi_rule_altitude(scc, ms, rot, et0, model, pargs, seed.inc, seed.raan, seed.omega)
        f = gram_exp_fit(sc, model, pargs, g.h_g, g.site)
        temp = Float64(sc["atmosphere"]["temperature_k"]); top = Float64(sc["atmosphere"]["fit_top_altitude_m"])
        pwf = pw_fit(f.pw_h, f.pw_rho, temp, "piecewise"); msc = with_fit(ms, exp_fit(f.rho_g, g.h_g, f.H, top, temp, "cell fit"))
        dir = joinpath(seed.out, seed.tag); mkpath(dir)
        seed_alt, recs = gram_calibrate(msc, scc, model, rot, epoch, dir, g.h_g; inc_deg=seed.inc, raan_deg=seed.raan, omega_deg=seed.omega, pwfit=pwf)
        write_rows_m(joinpath(dir, "calibration.csv"), recs)
        write_rows_m(joinpath(dir, "nodes.csv"), [(h_m=a, rho_kg_m3=b) for (a, b) in zip(f.pw_h, f.pw_rho)])
        core_t = model.core.initial_time
        row = (; planet=seed.planet, tag=seed.tag, status="ok", h_guess_m=g.h_g, seed_alt_m=seed_alt, vrel_first_mps=recs[end].vrel_mps, q_first_w_cm2=recs[end].q_w_cm2, q_target_w_cm2=recs[end].q_target,
            site_lat_rad=g.site.lat, site_lon_rad=g.site.lon, fit_H_m=f.H, pw_max_abs_lnrho=f.pw_max_abs_lnrho, hloc_at_hg_m=f.hloc_at_hg, n_cal_iter=length(recs),
            gram_core_epoch=string(core_t.year, "-", core_t.month, "-", core_t.day, " ", core_t.hour, ":", core_t.minute, ":", core_t.second), wall_s=time() - t0, error="")
        write_rows_m(joinpath(dir, "guess.csv"), [row]); return row
    catch err
        return (; planet=seed.planet, tag=seed.tag, status="cal_error", h_guess_m=NaN, seed_alt_m=NaN, vrel_first_mps=NaN, q_first_w_cm2=NaN, q_target_w_cm2=NaN, site_lat_rad=NaN, site_lon_rad=NaN, fit_H_m=NaN,
            pw_max_abs_lnrho=NaN, hloc_at_hg_m=NaN, n_cal_iter=0, gram_core_epoch="", wall_s=time() - t0, error=replace(first(sprint(showerror, err), 400), "\n" => " "))
    end
end

"""One flight of any Stage I arm. seed.ctrl: "fixed" | "adaptive" (Stage F two-sided) -> stagef_run_task; "adaptive_asym" (Stage G) -> stageg_run_task; "adaptive_h025" (Stage H) -> stageh_run_task."""
function stagei_task(seed)
    seed.ctrl in ("fixed", "adaptive") && return stagef_run_task(seed)
    seed.ctrl == "adaptive_asym" && return stageg_run_task(seed)
    startswith(seed.ctrl, "adaptive_h") && return stageh_run_task(seed)
    throw(ArgumentError("unknown ctrl $(seed.ctrl)"))
end
