# Stage C (planet generality): Mars truth runs and the Mars planet record for the frozen estimator.
# Physics of the engine run is the corridor example's pattern (constant-beta drag from a density fit, J2/J3 harmonics) at Mars:
# the Mars planet from SPICE (IAU pole), Mars50c harmonics degree 3 order 0, the Sun as third body (SPICE), single-layer exponential
# atmosphere (the Mars() struct fields). Elements are in the Mars body-equator inertial frame (z = pole): the truth initial state is
# rotated into J2000 with the engine's own pole and the truth pass table is rotated back, so truth and estimator share one element frame.
# No src/ or examples/ change. Usage: included by run_mars_generality.jl and run_mars_closedloop.jl (and by process workers).
include(joinpath(@__DIR__, "estimator_v2.jl"))   # common.jl (corridor example, process helpers), estimator v1 and v2, the planet record

const MARS_CONFIG_PATH = joinpath(@__DIR__, "mars_generality_config.toml")
load_mars(path::AbstractString=MARS_CONFIG_PATH) = TOML.parsefile(path)
mars_results_dir(root::AbstractString, mc) = joinpath(root, String(mc["meta"]["subdir"]))

# ---- planet, fit, rotation ---------------------------------------------------------------------------------------

mars_cell_tag(inc, raan, omega, e) = "I$(lpad(round(Int, 10inc), 4, '0'))_R$(lpad(round(Int, raan), 3, '0'))_W$(lpad(round(Int, omega), 3, '0'))_E$(e)"

"""Per-process Mars bundle: SPICE-backed planet, the truth gravity harmonics, the exponential density fit and beta."""
struct MarsSetup
    mc::Dict{String, Any}
    planet::Any
    harm::Any
    fit::Any
    beta::Float64
    area_m2::Float64
    wet_mass_kg::Float64
    dry_mass_kg::Float64
    apo_alt_m::Float64
end

"""Single-layer exponential as a `CorridorAtmosphereFit`: rho(h) = rho_ref exp(-(h - h_ref) / H) at geodetic altitude, extended above and below."""
function mars_density_fit(planet, mc)
    top = Float64(mc["atmosphere"]["fit_top_altitude_m"]); h0 = Float64(planet.h_ref); H = Float64(planet.H); r0 = Float64(planet.ρ_ref)
    return CorridorAtmosphereFit([h0, top], [r0, r0 * exp(-(top - h0) / H)], [H], :extend, :extend, Float64(mc["atmosphere"]["temperature_k"]),
        "Mars() struct: rho_ref, h_ref, H (single-layer exponential, ExponentialAtmosphereModel(planet))", Float64[], Float64[])
end

function mars_setup(mc=load_mars(); degree::Int=Int(mc["gravity"]["degree"]))
    planet = Mars("", SPICE_PATH)
    harm = GravitationalHarmonicsModel(degree, Int(mc["gravity"]["order"]), joinpath(REPO_ROOT, String(mc["gravity"]["coefficient_file"])), planet)
    v = mc["vehicle"]
    area = Float64(v["bus_dims_m"][1]) * Float64(v["bus_dims_m"][3]) + 2 * Float64(v["panel_dims_m"][2]) * Float64(v["panel_dims_m"][3])
    dry = Float64(v["bus_mass_kg"]) + 2 * Float64(v["panel_mass_each_kg"]); wet = dry + Float64(v["prop_mass_kg"])
    beta_der = wet / (Float64(v["cd"]) * area); bo = Float64(v["beta_override_kg_m2"])
    return MarsSetup(mc, planet, harm, mars_density_fit(planet, mc), bo > 0 ? bo : beta_der, area, wet, dry,
        Float64(mc["orbit"]["apoapsis_radius_m"]) - planet.Rp_e)
end

"""Rotation J2000 -> Mars body-equator inertial frame at `epoch` (columns x, y, z of the frame in J2000; z = the engine's pole at that epoch;
x along the node of the J2000 equator on the Mars equator, as `_body_equator_frame_rotation` of the telemetry builder). `v_J2000 = rot * v_B`."""
function mars_rot(planet, epoch)
    model = SpiceEphemeridesModel()
    et = SimulationModel.ephemerides_time_seconds(epoch, model)
    l = SimulationModel.planet_frame_lpi(planet, et, model)
    zh = SVector{3, Float64}(l[3, 1], l[3, 2], l[3, 3]); zh = zh / norm(zh)
    node = SVector{3, Float64}(-zh[2], zh[1], 0.0); xh = node / norm(node); yh = cross(zh, xh)
    return SMatrix{3, 3, Float64}(hcat(xh, yh, zh))
end

"""Planet record of the frozen estimator for Mars at one epoch."""
function mars_record(ms::MarsSetup, epoch)
    mc = ms.mc; b = mc["band"]
    return PlanetRecord(name="Mars", center=String(mc["planet"]["spice_center"]), planet=ms.planet, harm=ms.harm, fit=ms.fit, beta=ms.beta,
        q_low=Float64(b["q_low_w_cm2"]), q_high=Float64(b["q_high_w_cm2"]), q_hard=Float64(b["q_hard_w_cm2"]), w2m=Float64(b["w_cm2_to_w_m2"]),
        apo_alt_m=ms.apo_alt_m, body_names=String.(mc["planet"]["third_bodies"]), body_mu=[MU_BODY[String(x)] for x in mc["planet"]["third_bodies"]],
        rot=mars_rot(ms.planet, epoch))
end

# ---- truth run ---------------------------------------------------------------------------------------------------

"""Truth variants. `full` = the configured gravity (degree/order from the config: 3x0 in `mars_generality_config.toml`, 20x20 in `mars_generality_config_g20.toml`), Sun, drag.
Diagnostic variants: nosun, nodrag, j2only (degree 2), nosun_nodrag; `zonal` = configured degree, order 0 (tesseral terms off); `deg3` = degree 3 order 0 (the Stage C
gravity). `degree`/`order` of `nothing` mean the configured values; an explicit order is capped at the degree."""
function variant_flags(name::AbstractString)
    name == "full" && return (sun=true, drag=true, degree=nothing, order=nothing)
    name == "nosun" && return (sun=false, drag=true, degree=nothing, order=nothing)
    name == "nodrag" && return (sun=true, drag=false, degree=nothing, order=nothing)
    name == "j2only" && return (sun=true, drag=true, degree=2, order=0)
    name == "nosun_nodrag" && return (sun=false, drag=false, degree=nothing, order=nothing)
    name == "zonal" && return (sun=true, drag=true, degree=nothing, order=0)
    name == "deg3" && return (sun=true, drag=true, degree=3, order=0)
    throw(ArgumentError("unknown truth variant $name"))
end

"""Cartesian J2000 state at the apoapsis seed: osculating elements in the body-equator frame (the study's inclination / RAAN / omega), rotated by `rot`."""
function mars_initial_state(ms::MarsSetup, rot, seed_alt_m, inc_deg, raan_deg, omega_deg)
    rp = ms.planet.Rp_e + seed_alt_m; ra = ms.planet.Rp_e + ms.apo_alt_m
    el = (a=0.5 * (ra + rp), e=(ra - rp) / (ra + rp), i=deg2rad(inc_deg), raan=deg2rad(raan_deg), omega=deg2rad(omega_deg),
          nu=deg2rad(Float64(ms.mc["orbit"]["true_anomaly_deg"])))
    r_b, v_b = elements_to_cartesian(el; mu=ms.harm.gm_m3s2)
    return SVector{3, Float64}(rot * SVector{3, Float64}(r_b)), SVector{3, Float64}(rot * SVector{3, Float64}(v_b))
end

"""Truth `SimulationConfiguration` of a Mars cell (or of a restart from `state` = (r, v) in J2000 with mass `mass_kg`, used by closed-loop predictors)."""
function mars_args(ms::MarsSetup, epoch, horizon_s::Float64, dir::String; state, mass_kg=nothing, variant::AbstractString="full", controller=nothing, thruster=nothing)
    flags = variant_flags(variant); mc = ms.mc; v = mc["vehicle"]
    planet = ms.planet
    vdeg = flags.degree === nothing ? Int(mc["gravity"]["degree"]) : flags.degree
    vord = min(flags.order === nothing ? Int(mc["gravity"]["order"]) : flags.order, vdeg)
    harm = (vdeg == Int(mc["gravity"]["degree"]) && vord == Int(mc["gravity"]["order"])) ? ms.harm :
        GravitationalHarmonicsModel(vdeg, vord, joinpath(REPO_ROOT, String(mc["gravity"]["coefficient_file"])), planet)
    prop = mass_kg === nothing ? Float64(v["prop_mass_kg"]) : Float64(mass_kg) - ms.dry_mass_kg
    spacecraft = make_three_body_spacecraft(bus_dims=Tuple(Float64.(v["bus_dims_m"])), panel_dims=Tuple(Float64.(v["panel_dims_m"])),
        bus_mass=Float64(v["bus_mass_kg"]), panel_mass_each=Float64(v["panel_mass_each_kg"]), panel_offset_y=Float64(v["panel_offset_y_m"]),
        ic=CartesianInitialCondition(state[1], state[2]), prop_mass=prop, id=1)
    effectors = Any[harm]
    flags.drag && push!(effectors, CorridorBallisticDrag(ms.beta, ms.wet_mass_kg, ms.fit, :fit))
    flags.sun && push!(effectors, NBodyGravityModel(body_names=Tuple(String.(mc["planet"]["third_bodies"])),
        primary_body_name=String(mc["planet"]["name"]), planet=planet))
    base = make_example_config(planet=planet, spacecraft=spacecraft, mission_time=horizon_s, initial_time=epoch,
        dynamic_effectors=Tuple(effectors), density_model=ExponentialAtmosphereModel(planet), ephemerides_model=SpiceEphemeridesModel(),
        orientation_sim=false, keplerian=false, EI_km=Float64(mc["environment"]["EI_km"]), verbose=false, results=false, results_directory=dir)
    ig = mc["integration"]
    tol = IntegrationTolerances(reltol=Float64(ig["reltol"]), abstol=Float64(ig["abstol"]), reltol_orbit=Float64(ig["reltol"]), abstol_orbit=Float64(ig["abstol"]),
        reltol_atmosphere=Float64(ig["reltol"]), abstol_atmosphere=Float64(ig["abstol"]), dt_max=Float64(ig["dt_max_orbit_s"]),
        dt_max_orbit=Float64(ig["dt_max_orbit_s"]), dt_max_atmosphere=Float64(ig["dt_max_atmosphere_s"]))
    guidance, control = controller === nothing ? (base.guidance_model, base.control_model) : (GuidanceModel((controller,), [1.0e12]), ControlModel((thruster,), [1.0e12]))
    return SimulationConfiguration(file_paths=base.file_paths, simulation_settings=base.simulation_settings, mission_configuration=base.mission_configuration,
        environment_model=base.environment_model, dynamics_model=base.dynamics_model, guidance_model=guidance, navigation_model=base.navigation_model,
        control_model=control, initial_time=base.initial_time, integration_tolerances=tol, solver_config=SolverConfig(solver_mode=Symbol(ig["solver_mode"])))
end

"""One row per radial-event periapsis (same columns as `periapsis_table`); elements in the body-equator frame, mu = the harmonics GM."""
function mars_periapsis_table(sol, args, ms::MarsSetup, rot)
    env = args.environment_model; planet = env.planet
    events = physical_apsis_events(sol, _geodetic_of(args))
    et0 = SimulationModel.ephemerides_time_seconds(args.initial_time, env.ephemerides_model)
    w2m = Float64(ms.mc["band"]["w_cm2_to_w_m2"]); mu = ms.harm.gm_m3s2
    rows = NamedTuple[]
    for ev in events
        ev.kind == "periapsis" || continue
        lpi = SimulationModel.planet_frame_lpi(planet, et0 + ev.time_s, env.ephemerides_model)
        vrel = norm(lpi * ev.velocity_ii_mps - cross(planet.ω, lpi * ev.position_ii_m))
        rho = corridor_density(ms.fit, ev.physical_altitude_m)
        k = findfirst(c -> c.kind == "apoapsis" && c.time_s > ev.time_s, events)
        k === nothing && (k = findlast(c -> c.kind == "apoapsis" && c.time_s < ev.time_s, events))
        el = _corridor_element_summary(SVector{3, Float64}(rot' * ev.position_ii_m), SVector{3, Float64}(rot' * ev.velocity_ii_mps), mu)
        hq = altitude_of_density_fast(ms.fit, 2 * Float64(ms.mc["band"]["q_high_w_cm2"]) * w2m / vrel^3)
        hl = altitude_of_density_fast(ms.fit, 2 * Float64(ms.mc["band"]["q_low_w_cm2"]) * w2m / vrel^3)
        push!(rows, (pass=length(rows) + 1, time_s=ev.time_s, alt_m=ev.physical_altitude_m, lat_rad=ev.latitude_rad, radius_m=ev.radius_m, vrel_mps=vrel, rho=rho,
            q_w_cm2=0.5 * rho * vrel^3 / w2m, apo_alt_m=k === nothing ? NaN : events[k].physical_altitude_m, a_m=el.a, e=el.e, i_rad=el.i, raan_rad=el.raan,
            omega_rad=el.omega, mass_kg=ev.mass_kg, h_qhigh_m=hq, h_qlow_m=hl))
    end
    return rows
end

function mars_run(ms::MarsSetup, rot, epoch, horizon_s, dir; seed_alt_m, inc_deg, raan_deg, omega_deg, variant="full")
    state = mars_initial_state(ms, rot, seed_alt_m, inc_deg, raan_deg, omega_deg)
    args = mars_args(ms, epoch, Float64(horizon_s), dir; state=state, variant=variant)
    t0 = time(); sol = run_simulation(args; return_solution=true, isolate_state=false); wall = time() - t0
    return (rows=mars_periapsis_table(sol, args, ms, rot), sol, args, wall_s=wall, retcode=string(sol.retcode))
end

"""Calibrate the seed periapsis altitude so the first-pass peak Q equals the geometric mean of the band (the Earth rule, `calibrate_peak_q`)."""
function mars_calibrate(ms::MarsSetup, rot, epoch, dir; inc_deg, raan_deg, omega_deg, variant="full")
    mc = ms.mc; cal = mc["calibration"]; b = mc["band"]
    q_target = sqrt(Float64(b["q_low_w_cm2"]) * Float64(b["q_high_w_cm2"])); w2m = Float64(b["w_cm2_to_w_m2"])
    seed = Float64(mc["orbit"]["initial_periapsis_altitude_guess_m"]); records = NamedTuple[]
    for it in 1:Int(cal["max_iterations"])
        rp = ms.planet.Rp_e + seed; ra = ms.planet.Rp_e + ms.apo_alt_m
        horizon = Float64(cal["horizon_fraction_of_period"]) * 2π * sqrt(((rp + ra) / 2)^3 / ms.harm.gm_m3s2)
        res = mars_run(ms, rot, epoch, horizon, dir; seed_alt_m=seed, inc_deg, raan_deg, omega_deg, variant)
        isempty(res.rows) && throw(ErrorException("calibration run reached no periapsis (seed $seed m)"))
        r1 = res.rows[1]
        push!(records, (iteration=it, seed_alt_m=seed, alt_m=r1.alt_m, vrel_mps=r1.vrel_mps, q_w_cm2=r1.q_w_cm2, q_target=q_target))
        abs(r1.q_w_cm2 / q_target - 1) <= Float64(cal["q_tolerance_relative"]) && return seed, records
        seed += altitude_of_density_fast(ms.fit, 2 * q_target * w2m / r1.vrel_mps^3) - r1.alt_m
    end
    throw(ErrorException("peak-Q calibration did not converge: $(records[end])"))
end

write_rows_m(path, rows) = (mkpath(dirname(path)); CSV.write(path, DataFrame(rows)))

"""One Mars cell: calibrate, then fly the horizon; writes calibration.csv and passes_<variant>.csv. Plain-data row; never throws (failure -> status)."""
function mars_cell_task(seed)
    t0 = time()
    try
        mc = load_mars(seed.config_path); ms = mars_setup(mc)
        epoch = parse_epoch(mc["grid"]["epochs_utc"][seed.epoch_idx]); rot = mars_rot(ms.planet, epoch)
        dir = joinpath(seed.out, seed.tag); mkpath(dir)
        suffix = seed.variant == "full" ? "3b" : seed.variant
        if seed.variant == "full"
            seed_alt, recs = mars_calibrate(ms, rot, epoch, dir; inc_deg=seed.inc, raan_deg=seed.raan, omega_deg=seed.omega)
            write_rows_m(joinpath(dir, "calibration.csv"), recs)
        else   # diagnostic variants start from the accepted full-model seed (same initial state)
            seed_alt = CSV.read(joinpath(seed.out, seed.tag, "calibration.csv"), DataFrame).seed_alt_m[end]
        end
        cal_wall = time() - t0
        r = mars_run(ms, rot, epoch, Float64(seed.horizon_days) * SECONDS_PER_DAY, dir; seed_alt_m=seed_alt, inc_deg=seed.inc, raan_deg=seed.raan, omega_deg=seed.omega, variant=seed.variant)
        write_rows_m(joinpath(dir, "passes_$(suffix).csv"), r.rows)
        return (; tag=seed.tag, variant=seed.variant, inc_deg=seed.inc, raan_deg=seed.raan, omega_deg=seed.omega, epoch_idx=seed.epoch_idx, status="ok",
            seed_alt_m=seed_alt, cal_wall_s=cal_wall, run_wall_s=r.wall_s, retcode=r.retcode, passes=length(r.rows),
            t_end_days=isempty(r.rows) ? NaN : r.rows[end].time_s / SECONDS_PER_DAY, error="")
    catch err
        return (; tag=seed.tag, variant=seed.variant, inc_deg=seed.inc, raan_deg=seed.raan, omega_deg=seed.omega, epoch_idx=seed.epoch_idx, status="run_error",
            seed_alt_m=NaN, cal_wall_s=NaN, run_wall_s=NaN, retcode="", passes=0, t_end_days=NaN, error=replace(first(sprint(showerror, err), 300), "\n" => " "))
    end
end

mars_cells(mc) = [(tag=mars_cell_tag(i, r, w, e), inc=Float64(i), raan=Float64(r), omega=Float64(w), epoch_idx=e)
                  for e in eachindex(mc["grid"]["epochs_utc"]) for i in mc["grid"]["inclination_deg"] for r in mc["grid"]["raan_deg"] for w in mc["grid"]["omega_deg"]]
