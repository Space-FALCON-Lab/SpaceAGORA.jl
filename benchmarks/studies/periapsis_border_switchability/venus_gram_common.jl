# Venus with the native Venus-GRAM atmosphere as truth. Reuses venus_common.jl / mars_common.jl unchanged for the planet, rotation, initial state and pass-table
# conventions; adds (1) the Pi-based initial-periapsis rule and the per-cell exponential fit to GRAM MEAN density (data available at the cell's start only),
# (2) a truth configuration that samples the engine's GRAM density (CorridorBallisticDrag(:sampled), as the corridor example's GRAM backend), (3) the GRAM pass
# table (adds longitude, local scale height, GRAM band edges), (4) the cell task. No src/ or examples/ change. Included by venus_gram_generality.jl and by workers.
include(joinpath(@__DIR__, "venus_common.jl"))
using Roots
_corridor_gram_load!()   # load GRAMSuite and its extension at include time, so later calls (engine, workers) see its methods (world age)

const VENUS_GRAM_CONFIG_PATH = joinpath(@__DIR__, "venus_gram_config.toml")
function _deep_merge(a::AbstractDict, b::AbstractDict)
    out = Dict{String, Any}(a)
    for (k, v) in b; out[k] = (haskey(out, k) && out[k] isa AbstractDict && v isa AbstractDict) ? _deep_merge(out[k], v) : v; end
    return out
end
"""Base Venus config (unchanged) overlaid with the GRAM override tables."""
load_venus_gram(path::AbstractString=VENUS_GRAM_CONFIG_PATH) = _deep_merge(load_venus(VENUS_CONFIG_PATH), TOML.parsefile(path))

gram_pargs(planet, ei_km) = (args=(environment_model=(EI=ei_km, planet=planet), mission_configuration=(keplerian=false,)),)
"""GRAM mean density (kg/m3) at one point through the engine's own `getDensity` (the call the drag sample makes)."""
gram_rho(model, pargs, h, lat, lon, t) = Float64(Base.invokelatest(SimulationModel.EnvironmentModels.getDensity, model, Float64(h), Float64(lat), Float64(lon), Float64(t), true, pargs)[1])

"""Venus-GRAM model at `epoch` (direct sampling required, as the corridor GRAM backend; no fallback)."""
function gram_model(vgc, epoch)
    _corridor_gram_load!(); _corridor_gram_require_direct_sampling()
    return Base.invokelatest(SimulationModel.GRAMAtmosphereModel; planet_name=String(vgc["atmosphere_gram"]["planet_name"]), initial_time=epoch)
end

"""Venus bundle for GRAM truth: SPICE Venus (struct atmosphere constants unused), MGNP180U harmonics, beta; `fit` is a placeholder replaced per cell."""
function venus_gram_setup(vgc=load_venus_gram())
    p0 = Venus("", SPICE_PATH); g = vgc["gravity"]
    harm = GravitationalHarmonicsModel(Int(g["degree"]), Int(g["order"]), joinpath(REPO_ROOT, String(g["coefficient_file"])), p0)
    v = vgc["vehicle"]
    area = Float64(v["bus_dims_m"][1]) * Float64(v["bus_dims_m"][3]) + 2 * Float64(v["panel_dims_m"][2]) * Float64(v["panel_dims_m"][3])
    dry = Float64(v["bus_mass_kg"]) + 2 * Float64(v["panel_mass_each_kg"]); wet = dry + Float64(v["prop_mass_kg"])
    beta = wet / (Float64(v["cd"]) * area)
    return MarsSetup(vgc, p0, harm, exp_fit(1.0e-9, 140000.0, 4000.0, Float64(vgc["atmosphere"]["fit_top_altitude_m"]), Float64(vgc["atmosphere"]["temperature_k"]), "placeholder"), beta, area, wet, dry,
        Float64(vgc["orbit"]["apoapsis_radius_m"]) - p0.Rp_e)
end
with_fit(ms::MarsSetup, fit) = MarsSetup(ms.mc, ms.planet, ms.harm, fit, ms.beta, ms.area_m2, ms.wet_mass_kg, ms.dry_mass_kg, ms.apo_alt_m)

"""Single exponential layer as a `CorridorAtmosphereFit`: rho(h) = rho0 exp(-(h - h0)/H), extended above and below."""
exp_fit(rho0, h0, H, top, temp, src) = CorridorAtmosphereFit([h0, top], [rho0, rho0 * exp(-(top - h0) / H)], [H], :extend, :extend, temp, src, Float64[], Float64[])

"""Piecewise-exponential `CorridorAtmosphereFit` through GRAM nodes (ln-linear in altitude, continuous): layer k has H_k = (h_{k+1} - h_k) / ln(rho_k / rho_{k+1}); extended above and below."""
function pw_fit(hs, rhos, temp, src)
    Hs = [(hs[k + 1] - hs[k]) / log(rhos[k] / rhos[k + 1]) for k in 1:(length(hs) - 1)]
    return CorridorAtmosphereFit(Float64.(hs), Float64.(rhos), Hs, :extend, :extend, temp, src, Float64[], Float64[])
end

# ---- Pi-based initial periapsis and the per-cell fit (no truth trajectory) ---------------------------------------------

"""First-periapsis site and relative speed from two-body geometry at periapsis altitude `h` (apoapsis start, true anomaly 180 deg, half a Kepler period later)."""
function first_periapsis_site(ms::MarsSetup, rot, et0, h, inc_deg, raan_deg, omega_deg)
    mu = ms.harm.gm_m3s2; rp = ms.planet.Rp_e + h; ra = ms.planet.Rp_e + ms.apo_alt_m; a = 0.5 * (ra + rp); e = (ra - rp) / (ra + rp)
    t_p = pi * sqrt(a^3 / mu)
    r_b, v_b = elements_to_cartesian((a=a, e=e, i=deg2rad(inc_deg), raan=deg2rad(raan_deg), omega=deg2rad(omega_deg), nu=0.0); mu=mu)
    rJ = rot * SVector{3, Float64}(r_b); vJ = rot * SVector{3, Float64}(v_b)
    lpi = SimulationModel.planet_frame_lpi(ms.planet, et0 + t_p, SpiceEphemeridesModel())
    rpp = lpi * rJ; llh = SimulationModel.SimulationCallbacks.rtolatlong(SVector{3, Float64}(rpp), ms.planet)
    return (t_p=t_p, lat=llh[2], lon=llh[3], vrel=norm(lpi * vJ - cross(ms.planet.ω, rpp)), a=a, e=e)
end

"""Pi-based initial periapsis altitude: the root of 0.5 rho_GRAM(h, site) V_rel^3 = Q_target at the first-periapsis site, fixed-point on the site/speed (3 passes)."""
function pi_rule_altitude(vgc, ms, rot, et0, model, pargs, inc_deg, raan_deg, omega_deg)
    b = vgc["band"]; w2m = Float64(b["w_cm2_to_w_m2"]); q_t = sqrt(Float64(b["q_low_w_cm2"]) * Float64(b["q_high_w_cm2"])) * w2m
    g = vgc["guess"]; h = Float64(g["initial_altitude_m"]); lo, hi = Float64.(g["root_bracket_m"]); site = nothing
    for _ in 1:Int(g["fixed_point_iterations"])
        site = first_periapsis_site(ms, rot, et0, h, inc_deg, raan_deg, omega_deg)
        rho_t = 2 * q_t / site.vrel^3
        h = find_zero(x -> log(gram_rho(model, pargs, x, site.lat, site.lon, site.t_p)) - log(rho_t), (lo, hi), Roots.Bisection(); xatol=Float64(g["root_tolerance_m"]))
    end
    site = first_periapsis_site(ms, rot, et0, h, inc_deg, raan_deg, omega_deg)
    return (h_g=h, site=site, q_target_w_m2=q_t)
end

"""Estimator density inputs from GRAM MEAN density at the first-periapsis site over [h_g - below, h_g + above] (data available at the cell's start only):
(i) PRIMARY piecewise-exponential: GRAM values at nodes every `node_spacing_m`, ln-linear between; (ii) SECONDARY single exponential anchored at GRAM's density at h_g (the target-Q altitude), slope by least squares on ln rho over the same window.
Residuals are reported on a 500 m grid offset from the nodes (piecewise) and on the 1 km sample grid (single); H_local = -1/(d ln rho/dh) by central difference."""
function gram_exp_fit(vgc, model, pargs, h_g, site)
    f = vgc["fit"]; temp = Float64(vgc["atmosphere"]["temperature_k"]); top = Float64(vgc["atmosphere"]["fit_top_altitude_m"])
    lo = h_g - Float64(f["window_below_m"]); hi = h_g + Float64(f["window_above_m"]); rho_at(h) = gram_rho(model, pargs, h, site.lat, site.lon, site.t_p)
    hs = collect(lo:Float64(f["step_m"]):hi); y = [log(rho_at(h)) for h in hs]; x = hs .- h_g
    i0 = argmin(abs.(x)); icpt = y[i0]; slope = sum(x .* (y .- icpt)) / sum(x .^ 2)   # anchored at the Pi-rule altitude (x = 0), slope by least squares
    res = y .- (icpt .+ slope .* x); d = Float64(f["local_scale_height_step_m"])
    hl = [2d / (log(rho_at(h - d)) - log(rho_at(h + d))) for h in hs]
    H = -1 / slope; rho_g = exp(icpt)
    fit = exp_fit(rho_g, h_g, H, top, temp, "single exponential anchored at GRAM mean density at the Pi-rule altitude, LS slope over the window, first-periapsis site")
    nh = collect(lo:Float64(f["node_spacing_m"]):(hi + 1.0)); nr = [rho_at(h) for h in nh]
    pw = pw_fit(nh, nr, temp, "piecewise exponential through GRAM mean density at the first-periapsis site")
    mid = collect((lo + Float64(f["step_m"]) / 2):Float64(f["step_m"]):hi); pres = [log(corridor_density(pw, h)) - log(rho_at(h)) for h in mid]
    return (fit=fit, pw=pw, pw_h=nh, pw_rho=nr, pw_max_abs_lnrho=maximum(abs.(pres)), pw_rms_lnrho=sqrt(mean(pres .^ 2)),
        H=H, rho_g=rho_g, rms_lnrho=sqrt(mean(res .^ 2)), max_abs_lnrho=maximum(abs.(res)), res_at_hg=res[argmin(abs.(x))],
        hloc_min=minimum(hl), hloc_max=maximum(hl), hloc_at_hg=hl[argmin(abs.(x))], n_window=length(hs))
end

# ---- truth configuration (GRAM-sampled drag) -----------------------------------------------------------------------------

"""Truth `SimulationConfiguration` of a GRAM cell: 20x20 harmonics + Sun + `CorridorBallisticDrag(:sampled)` on the engine's GRAM density (EI = 2000 km)."""
function gram_args(ms::MarsSetup, vgc, model, epoch, horizon_s::Float64, dir::String; state)
    mc = vgc; v = mc["vehicle"]; planet = ms.planet
    spacecraft = make_three_body_spacecraft(bus_dims=Tuple(Float64.(v["bus_dims_m"])), panel_dims=Tuple(Float64.(v["panel_dims_m"])),
        bus_mass=Float64(v["bus_mass_kg"]), panel_mass_each=Float64(v["panel_mass_each_kg"]), panel_offset_y=Float64(v["panel_offset_y_m"]),
        ic=CartesianInitialCondition(state[1], state[2]), prop_mass=Float64(v["prop_mass_kg"]), id=1)
    effectors = Any[ms.harm, CorridorBallisticDrag(ms.beta, ms.wet_mass_kg, ms.fit, :sampled),
        NBodyGravityModel(body_names=Tuple(String.(mc["planet"]["third_bodies"])), primary_body_name=String(mc["planet"]["name"]), planet=planet)]
    base = make_example_config(planet=planet, spacecraft=spacecraft, mission_time=horizon_s, initial_time=epoch,
        dynamic_effectors=Tuple(effectors), density_model=model, ephemerides_model=SpiceEphemeridesModel(),
        orientation_sim=false, keplerian=false, EI_km=Float64(mc["atmosphere_gram"]["entry_interface_km"]), verbose=false, results=false, results_directory=dir)
    ig = mc["integration"]
    tol = IntegrationTolerances(reltol=Float64(ig["reltol"]), abstol=Float64(ig["abstol"]), reltol_orbit=Float64(ig["reltol"]), abstol_orbit=Float64(ig["abstol"]),
        reltol_atmosphere=Float64(ig["reltol"]), abstol_atmosphere=Float64(ig["abstol"]), dt_max=Float64(ig["dt_max_orbit_s"]),
        dt_max_orbit=Float64(ig["dt_max_orbit_s"]), dt_max_atmosphere=Float64(ig["dt_max_atmosphere_s"]))
    return SimulationConfiguration(file_paths=base.file_paths, simulation_settings=base.simulation_settings, mission_configuration=base.mission_configuration,
        environment_model=base.environment_model, dynamics_model=base.dynamics_model, guidance_model=base.guidance_model, navigation_model=base.navigation_model,
        control_model=base.control_model, initial_time=base.initial_time, integration_tolerances=tol, solver_config=SolverConfig(solver_mode=Symbol(ig["solver_mode"])))
end

"""Altitude where GRAM density at the site equals `rho_t`, by damped Newton steps in ln rho with the local scale height (3 steps; converged to well below 1 m for a smooth profile)."""
function gram_altitude_of_density(model, pargs, h0, lat, lon, t, rho_t; d=500.0)
    h = h0
    for _ in 1:4
        r = gram_rho(model, pargs, h, lat, lon, t); Hl = 2d / (log(gram_rho(model, pargs, h - d, lat, lon, t)) - log(gram_rho(model, pargs, h + d, lat, lon, t)))
        h += Hl * log(r / rho_t)
    end
    return h
end

"""One row per radial-event periapsis, GRAM truth: as `mars_periapsis_table` plus longitude, local scale height, the cell fit's density at the same altitude, GRAM band edges."""
function gram_periapsis_table(sol, args, ms::MarsSetup, rot, vgc, pwfit=nothing)
    env = args.environment_model; planet = env.planet; model = env.density_model; pargs = (args=args,)
    events = physical_apsis_events(sol, _geodetic_of(args))
    et0 = SimulationModel.ephemerides_time_seconds(args.initial_time, env.ephemerides_model)
    w2m = Float64(vgc["band"]["w_cm2_to_w_m2"]); mu = ms.harm.gm_m3s2; d = Float64(vgc["fit"]["local_scale_height_step_m"])
    qh = Float64(vgc["band"]["q_high_w_cm2"]); ql = Float64(vgc["band"]["q_low_w_cm2"]); rows = NamedTuple[]
    for ev in events
        ev.kind == "periapsis" || continue
        lpi = SimulationModel.planet_frame_lpi(planet, et0 + ev.time_s, env.ephemerides_model)
        rpp = lpi * ev.position_ii_m; vrel = norm(lpi * ev.velocity_ii_mps - cross(planet.ω, rpp))
        lon = SimulationModel.SimulationCallbacks.rtolatlong(SVector{3, Float64}(rpp), planet)[3]; lat = ev.latitude_rad; h = ev.physical_altitude_m
        rho = gram_rho(model, pargs, h, lat, lon, ev.time_s)
        hloc = 2d / (log(gram_rho(model, pargs, h - d, lat, lon, ev.time_s)) - log(gram_rho(model, pargs, h + d, lat, lon, ev.time_s)))
        k = findfirst(c -> c.kind == "apoapsis" && c.time_s > ev.time_s, events)
        k === nothing && (k = findlast(c -> c.kind == "apoapsis" && c.time_s < ev.time_s, events))
        el = _corridor_element_summary(SVector{3, Float64}(rot' * ev.position_ii_m), SVector{3, Float64}(rot' * ev.velocity_ii_mps), mu)
        hq = gram_altitude_of_density(model, pargs, h, lat, lon, ev.time_s, 2 * qh * w2m / vrel^3; d)
        hl = gram_altitude_of_density(model, pargs, h, lat, lon, ev.time_s, 2 * ql * w2m / vrel^3; d)
        push!(rows, (pass=length(rows) + 1, time_s=ev.time_s, alt_m=h, lat_rad=lat, radius_m=ev.radius_m, vrel_mps=vrel, rho=rho, q_w_cm2=0.5 * rho * vrel^3 / w2m,
            apo_alt_m=k === nothing ? NaN : events[k].physical_altitude_m, a_m=el.a, e=el.e, i_rad=el.i, raan_rad=el.raan, omega_rad=el.omega, mass_kg=ev.mass_kg,
            h_qhigh_m=hq, h_qlow_m=hl, lon_rad=lon, h_local_m=hloc, rho_fit=corridor_density(ms.fit, h), rho_pw=pwfit === nothing ? NaN : corridor_density(pwfit, h)))
    end
    return rows
end

function gram_run(ms::MarsSetup, vgc, model, rot, epoch, horizon_s, dir; seed_alt_m, inc_deg, raan_deg, omega_deg, pwfit=nothing)
    state = mars_initial_state(ms, rot, seed_alt_m, inc_deg, raan_deg, omega_deg)
    args = gram_args(ms, vgc, model, epoch, Float64(horizon_s), dir; state=state)
    t0 = time(); sol = run_simulation(args; return_solution=true, isolate_state=false); wall = time() - t0
    return (rows=gram_periapsis_table(sol, args, ms, rot, vgc, pwfit), sol, args, wall_s=wall, retcode=string(sol.retcode))
end

"""Calibrate the seed periapsis altitude so the first-pass peak Q equals the geometric mean of the band (the Earth/Mars rule), starting from the Pi-rule altitude:
exponential step with the cell fit's H, then secant on ln Q after two points."""
function gram_calibrate(ms::MarsSetup, vgc, model, rot, epoch, dir, h_start; inc_deg, raan_deg, omega_deg, pwfit=nothing)
    cal = vgc["calibration"]; b = vgc["band"]; q_target = sqrt(Float64(b["q_low_w_cm2"]) * Float64(b["q_high_w_cm2"]))
    seed = h_start; records = NamedTuple[]; H = Float64(ms.fit.scale_heights_m[1])
    for it in 1:Int(cal["max_iterations"])
        rp = ms.planet.Rp_e + seed; ra = ms.planet.Rp_e + ms.apo_alt_m
        horizon = Float64(cal["horizon_fraction_of_period"]) * 2π * sqrt(((rp + ra) / 2)^3 / ms.harm.gm_m3s2)
        res = gram_run(ms, vgc, model, rot, epoch, horizon, dir; seed_alt_m=seed, inc_deg, raan_deg, omega_deg, pwfit)
        isempty(res.rows) && throw(ErrorException("calibration run reached no periapsis (seed $seed m)"))
        r1 = res.rows[1]
        push!(records, (iteration=it, seed_alt_m=seed, alt_m=r1.alt_m, vrel_mps=r1.vrel_mps, q_w_cm2=r1.q_w_cm2, q_target=q_target))
        abs(r1.q_w_cm2 / q_target - 1) <= Float64(cal["q_tolerance_relative"]) && return seed, records
        step = H * log(r1.q_w_cm2 / q_target)
        if it >= 2
            s0, s1 = records[end - 1].seed_alt_m, records[end].seed_alt_m; l0, l1 = log(records[end - 1].q_w_cm2), log(records[end].q_w_cm2)
            sl = (l1 - l0) / (s1 - s0)
            isfinite(sl) && sl < 0 && (step = -log(r1.q_w_cm2 / q_target) / sl)   # secant on ln Q against seed altitude
        end
        seed += step
    end
    throw(ErrorException("peak-Q calibration did not converge: $(records[end])"))
end

"""One GRAM cell: calibrate, then fly the horizon. `seed` carries the cell's Pi-rule altitude and fit (from the guess phase). Plain-data row; never throws."""
function venus_gram_cell_task(seed)
    t0 = time()
    try
        vgc = load_venus_gram(seed.config_path); ms0 = venus_gram_setup(vgc)
        epoch = parse_epoch(vgc["grid"]["epochs_utc"][seed.epoch_idx]); rot = mars_rot(ms0.planet, epoch)
        ms = with_fit(ms0, exp_fit(seed.fit_rho_g, seed.h_guess_m, seed.fit_H_m, Float64(vgc["atmosphere"]["fit_top_altitude_m"]), Float64(vgc["atmosphere"]["temperature_k"]), "cell fit"))
        pwf = pw_fit(seed.pw_h, seed.pw_rho, Float64(vgc["atmosphere"]["temperature_k"]), "piecewise")
        model = gram_model(vgc, epoch); core_t = model.core.initial_time
        dir = joinpath(seed.out, seed.tag); mkpath(dir)
        seed_alt, recs = gram_calibrate(ms, vgc, model, rot, epoch, dir, seed.h_guess_m; inc_deg=seed.inc, raan_deg=seed.raan, omega_deg=seed.omega, pwfit=pwf)
        write_rows_m(joinpath(dir, "calibration.csv"), recs)
        cal_wall = time() - t0
        r = gram_run(ms, vgc, model, rot, epoch, Float64(seed.horizon_days) * SECONDS_PER_DAY, dir; seed_alt_m=seed_alt, inc_deg=seed.inc, raan_deg=seed.raan, omega_deg=seed.omega, pwfit=pwf)
        write_rows_m(joinpath(dir, "passes_3b.csv"), r.rows)
        return (; tag=seed.tag, variant="full", inc_deg=seed.inc, raan_deg=seed.raan, omega_deg=seed.omega, epoch_idx=seed.epoch_idx, status="ok",
            seed_alt_m=seed_alt, cal_wall_s=cal_wall, run_wall_s=r.wall_s, retcode=r.retcode, passes=length(r.rows),
            t_end_days=isempty(r.rows) ? NaN : r.rows[end].time_s / SECONDS_PER_DAY, gram_core_epoch=string(core_t.year, "-", core_t.month, "-", core_t.day, " ", core_t.hour, ":", core_t.minute, ":", core_t.second), error="")
    catch err
        return (; tag=seed.tag, variant="full", inc_deg=seed.inc, raan_deg=seed.raan, omega_deg=seed.omega, epoch_idx=seed.epoch_idx, status="run_error",
            seed_alt_m=NaN, cal_wall_s=NaN, run_wall_s=NaN, retcode="", passes=0, t_end_days=NaN, gram_core_epoch="", error=replace(first(sprint(showerror, err), 400), "\n" => " "))
    end
end

"""Process-pool launcher: `run_monte_carlo`, `threads=:auto`, `campaign_route_features(density_family="gram", ...)`; the route is asserted to be :process (native GRAM is not thread-safe)."""
function run_gram_campaign(task, seeds, workers::Int; mission_time_s::Float64, worker_files::Vector{String})
    pool = campaign_process_pool(); ids = ensure_process_workers!(pool, workers)
    Distributed.remotecall_eval(Main, ids, :(ENV["PC_SWEEP_WORKER"] = "1"))
    for f in worker_files; Distributed.remotecall_eval(Main, ids, :(include($f))); end
    features = SpaceAGORA.SimulationCampaigns.campaign_route_features(samples=length(seeds), n_sats=1, density_family="gram", mission_time_s=mission_time_s)
    tuning = SpaceAGORA.SimulationCampaigns.OuterRouteTuning(adaptive_enabled=false, process_workers_resident=length(ids), process_max_workers=workers)
    res = run_monte_carlo(task, seeds; threads=:auto, route_features=features, route_state=SpaceAGORA.SimulationCampaigns.OuterRouteState(), route_tuning=tuning)
    res.route === :process || throw(ErrorException("GRAM campaign ran on route $(res.route), not :process"))
    return res
end
