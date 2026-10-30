# Stage F: probabilistic density buffer with native GRAM density as truth at Venus, Mars and Earth. Reuses venus_gram_common.jl (planet-generic GRAM helpers: Pi-rule altitude,
# per-cell fits, GRAM truth configuration, pass table, calibration) and the frozen estimator v5 unchanged. New code only: per-planet setup from stagef_config.toml, the
# perturbed-density arm (opt-in SPACEAGORA_GRAM_DENSITY_PERTURBATION=pass, read, not edited), the in-run density measurement, the fixed / adaptive margin policy wrapper, and the
# self-contained case functions that one run_monte_carlo campaign dispatches over all planets together. No src/ or examples/ change.
include(joinpath(@__DIR__, "venus_gram_common.jl"))
include(joinpath(@__DIR__, "estimator_v3.jl"))
include(joinpath(@__DIR__, "estimator_v5.jl"))

const STAGEF_CONFIG_PATH = joinpath(@__DIR__, "stagef_config.toml")
stagef_meta(path::AbstractString=STAGEF_CONFIG_PATH) = TOML.parsefile(path)
"""Flat per-planet dict (the schema the venus_gram_common.jl helpers read): [common.*] merged with [planets.<name>.*]."""
function load_stagef(planet::AbstractString; path::AbstractString=STAGEF_CONFIG_PATH)
    t = TOML.parsefile(path); sc = _deep_merge(t["common"], t["planets"][planet])
    sc["atmosphere_gram"]["planet_name"] = planet
    return sc
end

function stagef_planet(name::AbstractString)
    n = lowercase(name)
    n == "venus" && return Venus("", SPICE_PATH)
    n == "mars" && return Mars("", SPICE_PATH)
    n == "earth" && return Earth("", SPICE_PATH)
    throw(ArgumentError("unknown planet $name"))
end

"""MarsSetup (the planet-agnostic bundle) for a Stage F planet; the density fit is a placeholder replaced per cell."""
function stagef_setup(sc)
    p0 = stagef_planet(String(sc["planet"]["name"])); g = sc["gravity"]
    harm = GravitationalHarmonicsModel(Int(g["degree"]), Int(g["order"]), joinpath(REPO_ROOT, String(g["coefficient_file"])), p0)
    v = sc["vehicle"]
    area = Float64(v["bus_dims_m"][1]) * Float64(v["bus_dims_m"][3]) + 2 * Float64(v["panel_dims_m"][2]) * Float64(v["panel_dims_m"][3])
    dry = Float64(v["bus_mass_kg"]) + 2 * Float64(v["panel_mass_each_kg"]); wet = dry + Float64(v["prop_mass_kg"])
    o = sc["orbit"]; apo = haskey(o, "apoapsis_radius_m") ? Float64(o["apoapsis_radius_m"]) - p0.Rp_e : Float64(o["apoapsis_altitude_m"])
    return MarsSetup(sc, p0, harm, exp_fit(1.0e-9, 140000.0, 4000.0, Float64(sc["atmosphere"]["fit_top_altitude_m"]), Float64(sc["atmosphere"]["temperature_k"]), "placeholder"),
        wet / (Float64(v["cd"]) * area), area, wet, dry, apo)
end

"""FNV-1a 32-bit hash of a string mapped into GRAM's documented seed range 1 .. 2^24-1 (one fixed seed per planet and cell)."""
function stagef_seed(s::AbstractString)
    h = UInt32(2166136261)
    for b in codeunits(s); h = (h ⊻ UInt32(b)) * UInt32(16777619); end
    return 1 + Int(h % UInt32(16777215))
end

function stagef_gram_model(sc, epoch, seed::Int; kw...)       # kw: optional per-level GRAM keywords (Revision 17a); empty = the unchanged call
    _corridor_gram_load!(); _corridor_gram_require_direct_sampling()
    return Base.invokelatest(SimulationModel.GRAMAtmosphereModel; planet_name=String(sc["atmosphere_gram"]["planet_name"]), initial_time=epoch, seed=seed, kw...)
end

"""Set the perturbation environment of the NEXT run_simulation in this process (the callback reads it when the run builds its callback set)."""
function stagef_env!(perturbed::Bool, meta)
    ENV["SPACEAGORA_GRAM_DENSITY_PERTURBATION"] = perturbed ? String(meta["stagef"]["perturbation_mode"]) : "off"
    ENV["SPACEAGORA_GRAM_DENSITY_PERTURBATION_PASS_RESEED"] = String(meta["stagef"]["perturbation_reseed"])
    perturbed ? (ENV["SPACEAGORA_GRAM_DENSITY_PERTURBATION_LOG"] = "results_directory") : delete!(ENV, "SPACEAGORA_GRAM_DENSITY_PERTURBATION_LOG")
    return nothing
end

# ---- phase cal: Pi-rule altitude, start-site fit, mean-density peak-Q calibration (one case per planet and cell) ---------------------------------------

"""Calibration case; writes <seed.out>/<tag>/{guess.csv, nodes.csv, calibration.csv}. Plain-data row; never throws."""
function stagef_cal_task(seed)
    t0 = time()
    try
        meta = stagef_meta(seed.config_path); sc = load_stagef(seed.planet; path=seed.config_path); ms = stagef_setup(sc)
        epoch = parse_epoch(sc["grid"]["epochs_utc"][seed.epoch_idx]); rot = mars_rot(ms.planet, epoch)
        et0 = SimulationModel.ephemerides_time_seconds(epoch, SpiceEphemeridesModel()); stagef_env!(false, meta)
        model = stagef_gram_model(sc, epoch, stagef_seed(seed.planet * seed.tag)); pargs = gram_pargs(ms.planet, Float64(sc["atmosphere_gram"]["entry_interface_km"]))
        g = pi_rule_altitude(sc, ms, rot, et0, model, pargs, seed.inc, seed.raan, seed.omega)
        f = gram_exp_fit(sc, model, pargs, g.h_g, g.site)
        temp = Float64(sc["atmosphere"]["temperature_k"]); top = Float64(sc["atmosphere"]["fit_top_altitude_m"])
        pwf = pw_fit(f.pw_h, f.pw_rho, temp, "piecewise"); msc = with_fit(ms, exp_fit(f.rho_g, g.h_g, f.H, top, temp, "cell fit"))
        dir = joinpath(seed.out, seed.tag); mkpath(dir)
        seed_alt, recs = gram_calibrate(msc, sc, model, rot, epoch, dir, g.h_g; inc_deg=seed.inc, raan_deg=seed.raan, omega_deg=seed.omega, pwfit=pwf)
        write_rows_m(joinpath(dir, "calibration.csv"), recs)
        write_rows_m(joinpath(dir, "nodes.csv"), [(h_m=a, rho_kg_m3=b) for (a, b) in zip(f.pw_h, f.pw_rho)])
        core_t = model.core.initial_time
        row = (; planet=seed.planet, tag=seed.tag, status="ok", h_guess_m=g.h_g, seed_alt_m=seed_alt, vrel_first_mps=recs[end].vrel_mps, q_first_w_cm2=recs[end].q_w_cm2,
            site_lat_rad=g.site.lat, site_lon_rad=g.site.lon, fit_H_m=f.H, pw_max_abs_lnrho=f.pw_max_abs_lnrho, hloc_at_hg_m=f.hloc_at_hg, n_cal_iter=length(recs),
            gram_core_epoch=string(core_t.year, "-", core_t.month, "-", core_t.day, " ", core_t.hour, ":", core_t.minute, ":", core_t.second), wall_s=time() - t0, error="")
        write_rows_m(joinpath(dir, "guess.csv"), [row]); return row
    catch err
        return (; planet=seed.planet, tag=seed.tag, status="cal_error", h_guess_m=NaN, seed_alt_m=NaN, vrel_first_mps=NaN, q_first_w_cm2=NaN, site_lat_rad=NaN, site_lon_rad=NaN, fit_H_m=NaN,
            pw_max_abs_lnrho=NaN, hloc_at_hg_m=NaN, n_cal_iter=0, gram_core_epoch="", wall_s=time() - t0, error=replace(first(sprint(showerror, err), 400), "\n" => " "))
    end
end

# ---- density measurement, predictors, margin policy ---------------------------------------------------------------------------------------------

"""Held perturbation factor r at time `t` from the pass interpolant of the engine's perturbation state (1.0 when no perturbation mode is installed)."""
function stagef_pass_r(p, t::Float64)
    hasproperty(p.shared_buffers, :gram_density_perturbation) || return 1.0
    st = p.shared_buffers.gram_density_perturbation[]; st === nothing && return 1.0
    knots = st.pass_r[1]; n = length(knots); n == 0 && return 1.0; n == 1 && return Float64(knots[1])
    x = (t - st.pass_t0[1]) / st.pass_dt_s; x <= 0 && return Float64(knots[1]); x >= n - 1 && return Float64(knots[n])
    i = floor(Int, x); f = x - i
    return (1 - f) * Float64(knots[i + 1]) + f * Float64(knots[i + 2])
end

"""One callback at each apsis crossing (a second continuous callback with the same condition retriggers forever, so the example's `_apsis_observation_callback` record is reproduced here):
it counts periapses for the controller and takes the controller's per-pass density measurement y = ln(rho_eff / rho_start_fit(h)), rho_eff = GRAM mean density at the point times r (assumed exact; plan section 4)."""
function stagef_apsis_callback(obs::Vector{NamedTuple}, meas::Vector{NamedTuple}, args, controller, model, pargs, pwfit)
    env = args.environment_model; planet = env.planet; et0 = SimulationModel.ephemerides_time_seconds(args.initial_time, env.ephemerides_model); geodetic_of = _geodetic_of(args)
    condition(u, t, integrator) = dot(u.sc[1].pos, u.sc[1].vel)
    function record!(integrator, kind::String)
        integrator.t > 1e-3 || return nothing
        t = Float64(integrator.t); sc1 = integrator.u.sc[1]
        push!(obs, (time_s=t, kind=kind, radius_m=norm(sc1.pos), physical_altitude_m=geodetic_of(t, sc1.pos)[1]))
        kind == "periapsis" || return nothing
        controller.periapsis_count += 1
        lpi = SimulationModel.planet_frame_lpi(planet, et0 + t, env.ephemerides_model); rp = lpi * SVector{3, Float64}(sc1.pos); alt, lat = geodetic_altitude_latitude(rp, planet)
        lon = SimulationModel.SimulationCallbacks.rtolatlong(SVector{3, Float64}(rp), planet)[3]
        rho = gram_rho(model, pargs, alt, lat, lon, t); r = stagef_pass_r(integrator.p, t)
        push!(meas, (; t_s=t, alt_m=alt, lat_rad=lat, lon_rad=lon, rho_mean=rho, r=r, y=log(rho * r / corridor_density(pwfit, alt))))
        return nothing
    end
    return ContinuousCallback(condition, i -> record!(i, "periapsis"), i -> record!(i, "apoapsis"); save_positions=(false, false))
end

const STAGEF_PREDICTORS = ("P1", "P2_3", "P2_5", "P2_10", "P3")
"""One-pass-ahead prediction of the next residual from the measured past `ys` and periapsis longitudes `lons`; the same definitions as density_buffer_analysis.predict_run."""
function stagef_predict_next(name::AbstractString, ys::Vector{Float64}, lons::Vector{Float64})
    n = length(ys); n == 0 && return 0.0
    name == "P1" && return ys[n]
    if startswith(name, "P2_")
        a = 2.0 / (parse(Int, name[4:end]) + 1); ew = ys[1]
        for k in 2:n; ew = a * ys[k] + (1 - a) * ew; end
        return ew
    end
    if name == "P3"
        # python (0-based, target j = n): pairs k = max(1, j-10) .. j-1; dy = y[k]-y[k-1]; dx = wrapped lon[k]-lon[k-1]; needs >= 4 pairs
        wrap(x) = atan(sin(x), cos(x)); ks = max(1, n - 10):(n - 1)
        length(ks) >= 4 || return ys[n]
        dy = [ys[k + 1] - ys[k] for k in ks]; dx = [wrap(lons[k + 1] - lons[k]) for k in ks]
        den = sum(dx .^ 2); m = den > 0 ? sum(dy .* dx) / den : 0.0
        return ys[n] + m * wrap(lons[n] - lons[n - 1])
    end
    throw(ArgumentError("unknown predictor $name"))
end

mutable struct StageFPolicy <: SimulationModel.AbstractGuidanceModel
    inner::CorridorGuidancePolicy
    make_predictor::Any            # (p, r, v, mass, t) -> (dv -> next periapsis geodetic altitude, m)
    mode::Symbol                   # :fixed | :adaptive
    fit::Any                       # start piecewise fit (altitude-only)
    vrel_first::Float64
    q_low_w_m2::Float64
    q_high_w_m2::Float64
    margin_m::Float64
    meas::Vector{NamedTuple}
    pred_name::String
    e_hi::Float64                  # upper one-sided quantile of the one-pass-ahead error (density higher than predicted -> Q too high)
    e_lo::Float64                  # lower one-sided quantile
    min_band_m::Float64
    edge_log::Vector{NamedTuple}
end

"""Corridor edges the controller flies at this decision: fixed = band at the first-pass speed with the fixed margin; adaptive = edges moved by the predicted density and its error quantiles."""
function stagef_edges(pol::StageFPolicy)
    V3 = pol.vrel_first^3; hq0 = altitude_of_density_fast(pol.fit, 2 * pol.q_high_w_m2 / V3); hl0 = altitude_of_density_fast(pol.fit, 2 * pol.q_low_w_m2 / V3)
    n = length(pol.meas)
    if pol.mode === :fixed
        return (lo=hq0, hi=hl0, margin=pol.margin_m, target=(hq0 + hl0) / 2, yhat=NaN, collapsed=false, n=n)
    end
    yhat = stagef_predict_next(pol.pred_name, Float64[m.y for m in pol.meas], Float64[m.lon_rad for m in pol.meas])
    lo = altitude_of_density_fast(pol.fit, 2 * pol.q_high_w_m2 / V3 * exp(-(yhat + pol.e_hi)))
    hi = altitude_of_density_fast(pol.fit, 2 * pol.q_low_w_m2 / V3 * exp(-(yhat + pol.e_lo)))
    collapsed = !(hi - lo >= pol.min_band_m)
    if collapsed
        c = altitude_of_density_fast(pol.fit, 2 * sqrt(pol.q_low_w_m2 * pol.q_high_w_m2) / V3 * exp(-yhat)); lo = c - pol.min_band_m / 2; hi = c + pol.min_band_m / 2
    end
    return (lo=lo, hi=hi, margin=0.0, target=(lo + hi) / 2, yhat=yhat, collapsed=collapsed, n=n)
end

function SimulationModel.GuidanceHooks.calcGuidanceEffect!(wrapper::StageFPolicy, u, p, t::Float64, sat_idx::Int64)
    policy = wrapper.inner; commands = p.shared_buffers.maneuver_commands
    (sat_idx < 1 || sat_idx > length(commands)) && return nothing
    if !policy.armed
        commands[sat_idx] = SimulationModel.CommandTypes.PropulsiveManeuverCommand(valid=true, source_orbit=policy.periapsis_count + 1)
        return nothing
    end
    r = SVector{3, Float64}(u.sc[sat_idx].pos); v = SVector{3, Float64}(u.sc[sat_idx].vel); mass = Float64(u.sc[sat_idx].mass)
    vehicle = policy.config.values["vehicle"]; exhaust_velocity = Float64(vehicle["specific_impulse_s"]) * 9.80665
    consumed_at_limit = Float64(vehicle["thrust_n"]) * policy.max_burn_duration_s / exhaust_velocity
    duration_limited_dv = consumed_at_limit < mass ? exhaust_velocity * log(mass / (mass - consumed_at_limit)) : Inf
    wet_mass = Float64(vehicle["dry_mass_kg"]) + Float64(vehicle["propellant_mass_kg"])
    delivered_so_far = mass < wet_mass ? exhaust_velocity * log(wet_mass / mass) : 0.0
    policy.remaining_budget_mps = max(0.0, Float64(policy.config.values["comparison"]["total_delta_v_budget_mps"]) - delivered_so_far)
    max_allowed_dv = min(policy.max_dv_mps, policy.remaining_budget_mps, duration_limited_dv)
    predictor = wrapper.make_predictor(p, r, v, mass, Float64(t))
    ed = stagef_edges(wrapper)
    choice = select_corridor_maneuver(predictor; lower_m=ed.lo, upper_m=ed.hi, margin_m=ed.margin, target_m=ed.target, policy=policy.policy,
        max_dv_mps=max_allowed_dv, deadband_m=policy.deadband_m, monotone_tolerance_m=Float64(policy.config.values["comparison"]["predictor_monotone_tolerance_m"]))
    policy.last_selection = choice; target_pass = policy.periapsis_count + 1
    push!(wrapper.edge_log, (; request_time_s=Float64(t), target_pass, n_meas=ed.n, yhat=ed.yhat, lower_m=ed.lo, upper_m=ed.hi, collapsed=ed.collapsed, delta_v_mps=choice.delta_v_mps, reason=String(choice.reason)))
    push!(policy.decisions, (request_time_s=Float64(t), source_orbit=target_pass, target_pass=target_pass, policy=String(policy.policy), delta_v_mps=choice.delta_v_mps,
        predicted_periapsis_altitude_m=choice.predicted_periapsis_altitude_m, target_periapsis_altitude_m=choice.target_periapsis_altitude_m, trigger_reason=String(choice.reason),
        feasible=choice.feasible, predictor_evaluations=choice.evaluations))
    commands[sat_idx] = SimulationModel.CommandTypes.PropulsiveManeuverCommand(valid=true, delta_v_mps=choice.feasible ? abs(choice.delta_v_mps) : 0.0,
        direction_rad=choice.delta_v_mps >= 0.0 ? 0.0 : π, source_orbit=target_pass)
    policy.armed = false
    return nothing
end

"""Frozen estimator v5 one-pass prediction context for one cell (planet record with the cell's start fit)."""
function stagef_est_state(ms::MarsSetup, epoch)
    c2 = TOML.parsefile(joinpath(STUDY_DIR, "estimator_v2_config.toml")); s = settings_from(c2); v = c2["v2"]
    o = V2Opts(sp_tide=Symbol(v["sp_tide"]), sp_j2=Symbol(v["sp_j2"]), n_sp=Int(v["n_sp_nodes"])); rec = mars_record(ms, epoch)
    # order-0 gravity (Mars, Earth) has no tesseral terms: v5 reduces to v3 there and tesseral_field fails on a 1-column coefficient matrix, so v3 is used (plan section 4)
    return (; ctx=make_context(rec, s), s, o, tab=SPTable(o.n_sp), rot=rec.rot, zf=zonal_field(ms.harm), zo=ZonalOpts(), tf=ms.harm.M == 0 ? nothing : tesseral_field(ms.harm))
end

"""Next geodetic periapsis altitude after a tangential `dv` at the next apoapsis, from the measured osculating state (as run_mars_closedloop.jl, predictor v5)."""
function stagef_predict_alt(st, r0, v0, dv, et_now)
    ctx = st.ctx; oe = _corridor_element_summary(SVector{3, Float64}(st.rot' * r0), SVector{3, Float64}(st.rot' * v0), ctx.mu)
    ra = oe.a * (1 + oe.e); va = sqrt(ctx.mu * (2 / ra - 1 / oe.a)); a2 = 1 / (2 / ra - (va + dv)^2 / ctx.mu); e2 = ra / a2 - 1
    init = (a=oe.a + (a2 - oe.a), e=oe.e + (e2 - oe.e), i=oe.i, raan=oe.raan, omega=oe.omega, nu=oe.nu)
    p = st.tf === nothing ? predict_history3(ctx, st.s, st.o, st.tab, st.zf, st.zo, 0.0, 0.0, 0.0, 0.0, et_now; max_passes=1, init_override=init)[1] :
        predict_history5(ctx, st.s, st.o, st.tab, st.zf, st.zo, st.tf, TesseralOpts(force=true), TimingOpts(), 0.0, 0.0, 0.0, 0.0, et_now; max_passes=1, init_override=init)[1]
    isempty(p.alt_m) && throw(ErrorException("estimator predicted no next periapsis (ended=$(p.ended))"))
    return p.alt_m[1]
end

"""CorridorConfig-compatible values for the example's controller classes (band as altitudes at the first-pass speed, as mars_corridor_config, with the cell's piecewise fit)."""
function stagef_corridor_config(ms::MarsSetup, vrel_first, margin_m)
    mc = ms.mc; b = mc["band"]; w2m = Float64(b["w_cm2_to_w_m2"]); v = mc["vehicle"]
    hq = altitude_of_density_fast(ms.fit, 2 * Float64(b["q_high_w_cm2"]) * w2m / vrel_first^3); hl = altitude_of_density_fast(ms.fit, 2 * Float64(b["q_low_w_cm2"]) * w2m / vrel_first^3)
    mid = (hq + hl) / 2
    values = Dict{String, Any}(
        "mission" => Dict{String, Any}("corridor_min_altitude_m" => hq, "corridor_max_altitude_m" => hl, "reference_altitude_m" => mid, "safety_margin_m" => margin_m),
        "vehicle" => Dict{String, Any}("dry_mass_kg" => ms.dry_mass_kg, "propellant_mass_kg" => Float64(v["prop_mass_kg"]), "thrust_n" => Float64(v["thrust_n"]), "specific_impulse_s" => Float64(v["specific_impulse_s"])),
        "comparison" => Dict{String, Any}("target_periapsis_altitude_m" => mid, "max_burn_delta_v_mps" => Float64(v["max_burn_delta_v_mps"]), "max_burn_duration_s" => Float64(v["max_burn_duration_s"]),
            "total_delta_v_budget_mps" => Float64(v["total_delta_v_budget_mps"]), "planning_lead_deg" => Float64(v["planning_lead_deg"]),
            "predictor_monotone_tolerance_m" => Float64(v["predictor_monotone_tolerance_m"]), "near_passive_deadband_m" => margin_m))
    return CorridorConfig("stagef_config.toml", values), hq, hl
end

# ---- pass table with the perturbation factor and local solar time ---------------------------------------------------------------------------------

"""(t_s, r) rows of the engine's pass interpolant (log kind 3) from the perturbation log of one run; empty when none."""
function stagef_read_r_log(path::AbstractString)
    isfile(path) || return (Float64[], Float64[])
    ts = Float64[]; rs = Float64[]
    for line in eachline(path)
        startswith(line, "3,") || continue
        f = split(line, ','); push!(ts, parse(Float64, f[4])); push!(rs, parse(Float64, f[6]))
    end
    return ts, rs
end
function stagef_interp(ts, rs, t)
    isempty(ts) && return 1.0
    t <= ts[1] && return rs[1]; t >= ts[end] && return rs[end]
    j = searchsortedlast(ts, t); w = (t - ts[j]) / (ts[j + 1] - ts[j]); return (1 - w) * rs[j] + w * rs[j + 1]
end

"""Local solar time (hours) at planet-fixed longitude `lon` and ephemeris time `et` (Sun direction through the same planet-fixed rotation)."""
function stagef_local_time_h(planet, center::AbstractString, et, lon)
    s = spkezr("SUN", et, "J2000", "NONE", center)[1]; lpi = SimulationModel.planet_frame_lpi(planet, et, SpiceEphemeridesModel())
    sp = lpi * SVector{3, Float64}(s[1:3]); lon_sun = atan(sp[2], sp[1])
    return mod(12.0 + (atan(sin(lon - lon_sun), cos(lon - lon_sun))) * 12 / π, 24.0)
end

function stagef_table(sol, args, ms::MarsSetup, rot, sc, pwfit, logpath::AbstractString)
    rows = gram_periapsis_table(sol, args, ms, rot, sc, pwfit); ts, rs = stagef_read_r_log(logpath)
    env = args.environment_model; et0 = SimulationModel.ephemerides_time_seconds(args.initial_time, env.ephemerides_model); center = String(sc["planet"]["spice_center"]); out = NamedTuple[]
    for r in rows
        rr = stagef_interp(ts, rs, r.time_s)
        push!(out, merge(r, (; r_pert=rr, rho_eff=r.rho * rr, q_eff_w_cm2=r.q_w_cm2 * rr, lst_h=stagef_local_time_h(env.planet, center, et0 + r.time_s, r.lon_rad))))
    end
    return out
end

# ---- phase run: one (planet, cell, arm, control) case ---------------------------------------------------------------------------------------------

"""One run case. seed fields: planet, tag, inc, raan, omega, epoch_idx, arm ("mean" | "pert"), ctrl ("open" | "fixed" | "adaptive"), run_id, out, cal_root, horizon_days,
config_path, adapt_pred, adapt_e_hi, adapt_e_lo. Plain-data row; never throws."""
function stagef_run_task(seed)
    t0 = time()
    try
        meta = stagef_meta(seed.config_path); sc = load_stagef(seed.planet; path=seed.config_path); ms0 = stagef_setup(sc)
        epoch = parse_epoch(sc["grid"]["epochs_utc"][seed.epoch_idx]); rot = mars_rot(ms0.planet, epoch)
        cdir = joinpath(seed.cal_root, seed.tag); G = CSV.read(joinpath(cdir, "guess.csv"), DataFrame)[1, :]; N = CSV.read(joinpath(cdir, "nodes.csv"), DataFrame)
        pwf = pw_fit(N.h_m, N.rho_kg_m3, Float64(sc["atmosphere"]["temperature_k"]), "piecewise"); ms = with_fit(ms0, pwf)
        pert = seed.arm == "pert"; stagef_env!(pert, meta)
        model = stagef_gram_model(sc, epoch, stagef_seed(seed.planet * seed.tag)); core_t = model.core.initial_time
        pargs = gram_pargs(ms.planet, Float64(sc["atmosphere_gram"]["entry_interface_km"]))
        dir = joinpath(seed.out, seed.run_id); mkpath(dir); horizon = Float64(seed.horizon_days) * SECONDS_PER_DAY
        state0 = mars_initial_state(ms, rot, Float64(G.seed_alt_m), seed.inc, seed.raan, seed.omega)
        base = gram_args(ms, sc, model, epoch, horizon, dir; state=state0)
        meas = NamedTuple[]; pol = nothing; extra = ()
        if seed.ctrl == "open"
            args = base
        else
            margin = Float64(meta["stagef"]["fixed_margin_m"]); cfg, hq, hl = stagef_corridor_config(ms, Float64(G.vrel_first_mps), margin); c = cfg.values["comparison"]
            inner = CorridorGuidancePolicy(cfg, :favorable, Symbol(meta["stagef"]["policy"]), Float64(c["max_burn_delta_v_mps"]), Float64(c["max_burn_duration_s"]),
                Float64(c["total_delta_v_budget_mps"]), Float64(c["near_passive_deadband_m"]), Float64(c["target_periapsis_altitude_m"]), deg2rad(Float64(c["planning_lead_deg"])),
                0, false, nothing, NamedTuple[], NamedTuple[], String[])
            st = stagef_est_state(ms, epoch); et0 = SimulationModel.ephemerides_time_seconds(epoch, SpiceEphemeridesModel())
            maker = (p, r, v, mass, t) -> (dv -> stagef_predict_alt(st, r, v, Float64(dv), et0 + t))
            w2m = Float64(sc["band"]["w_cm2_to_w_m2"])
            adaptive = seed.ctrl == "adaptive"
            pol = StageFPolicy(inner, maker, adaptive ? :adaptive : :fixed, pwf, Float64(G.vrel_first_mps), Float64(sc["band"]["q_low_w_cm2"]) * w2m, Float64(sc["band"]["q_high_w_cm2"]) * w2m,
                margin, meas, adaptive ? String(seed.adapt_pred) : "", adaptive ? Float64(seed.adapt_e_hi) : 0.0, adaptive ? Float64(seed.adapt_e_lo) : 0.0,
                Float64(meta["analysis"]["adapt_min_band_m"]), NamedTuple[])
            v = cfg.values["vehicle"]; thruster = BaseThrusterModel([Float64(v["thrust_n"])], [0.0], [0.0], [-1.0], [-1.0], [Float64(v["specific_impulse_s"])])
            args = SimulationConfiguration(file_paths=base.file_paths, simulation_settings=base.simulation_settings, mission_configuration=base.mission_configuration,
                environment_model=base.environment_model, dynamics_model=base.dynamics_model, guidance_model=GuidanceModel((pol,), [1.0e12]), navigation_model=base.navigation_model,
                control_model=ControlModel((thruster,), [1.0e12]), initial_time=base.initial_time, integration_tolerances=base.integration_tolerances, solver_config=base.solver_config)
            obs = NamedTuple[]
            extra = (stagef_apsis_callback(obs, meas, args, inner, model, pargs, pwf), _preapoapsis_guidance_callback(inner))
        end
        tr = time(); sol = run_simulation(args; return_solution=true, isolate_state=false, extra_callbacks=extra); run_wall = time() - tr
        logpath = joinpath(dir, "gram_density_perturbation_log.csv")
        rows = stagef_table(sol, args, ms, rot, sc, pwf, logpath); isfile(logpath) && rm(logpath; force=true); isfile(logpath * ".summary.toml") && rm(logpath * ".summary.toml"; force=true)
        write_rows_m(joinpath(dir, "passes.csv"), rows)
        dv = NaN; burns = 0; stopped = false; failures = ""
        if pol !== nothing
            inner = pol.inner; ve = Float64(sc["vehicle"]["specific_impulse_s"]) * 9.80665; dv = ve * log(ms.wet_mass_kg / Float64(sol.u[end].sc[1].mass)); burns = length(inner.burn_plans)
            stopped = !isempty(inner.failures); failures = join(inner.failures, " | ")
            write_rows_m(joinpath(dir, "decisions.csv"), pol.edge_log); write_rows_m(joinpath(dir, "measured.csv"), meas)
        end
        ret = string(sol.retcode)
        return (; planet=seed.planet, tag=seed.tag, run_id=seed.run_id, arm=seed.arm, ctrl=seed.ctrl, status="ok", retcode=ret, passes=length(rows), t_end_days=sol.t[end] / SECONDS_PER_DAY,
            horizon_days=seed.horizon_days, completed=ret == "Success" && !stopped, controller_stopped=stopped, controller_failures=failures, total_dv_mps=dv, burns=burns,
            collapsed_decisions=pol === nothing ? 0 : count(e -> e.collapsed, pol.edge_log), gram_core_epoch=string(core_t.year, "-", core_t.month, "-", core_t.day, " ", core_t.hour, ":", core_t.minute, ":", core_t.second),
            run_wall_s=run_wall, wall_s=time() - t0, error="")
    catch err
        return (; planet=seed.planet, tag=seed.tag, run_id=seed.run_id, arm=seed.arm, ctrl=seed.ctrl, status="run_error", retcode="", passes=0, t_end_days=NaN, horizon_days=seed.horizon_days,
            completed=false, controller_stopped=false, controller_failures="", total_dv_mps=NaN, burns=0, collapsed_decisions=0, gram_core_epoch="", run_wall_s=NaN, wall_s=time() - t0,
            error=replace(first(sprint(showerror, err), 500), "\n" => " "))
    end
end
