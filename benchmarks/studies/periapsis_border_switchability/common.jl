# Shared truth-run helpers for the periapsis border switchability study.
# Physics is the corridor example's (included, not copied); this file only adds Sun+Moon
# third-body gravity (NBodyGravityModel, SPICE), the heat-rate diagnostic, calibration to a
# target peak Q, and the perturber-plane angles.
include(joinpath(@__DIR__, "..", "..", "..", "examples", "AGORA_Earth_Periapsis_Corridor.jl"))
using Dates
using Distributed
using Roots
using SPICE

const STUDY_DIR = @__DIR__
const STUDY_CONFIG_PATH = joinpath(STUDY_DIR, "study_config.toml")
const CORRIDOR_CONFIG_PATH = joinpath(REPO_ROOT, "examples", "periapsis_corridor_config.toml")
const SECONDS_PER_DAY = 86_400.0

load_study(path::AbstractString=STUDY_CONFIG_PATH) = TOML.parsefile(path)

"""Results root: `SPACEAGORA_BORDER_RESULTS` if set, else `results_root` from study_config.toml; relative paths resolve against the repo root."""
function results_root(study)
    env = get(ENV, "SPACEAGORA_BORDER_RESULTS", "")
    return abspath(joinpath(REPO_ROOT, isempty(env) ? String(study["results_root"]) : env))
end

"""Ballistic coefficient (kg/m2) and its derivation, from the study vehicle table and the corridor wet mass."""
function study_beta(study, corridor_values)
    v = study["vehicle"]
    area = Float64(v["bus_ref_area_m2"]) + Float64(v["n_wings"]) * Float64(v["wing_ref_area_m2"])
    wet = Float64(corridor_values["vehicle"]["dry_mass_kg"]) + Float64(corridor_values["vehicle"]["propellant_mass_kg"])
    derived = wet / (Float64(v["cd"]) * area)
    override = Float64(v["beta_override_kg_m2"])
    return (beta=override > 0 ? override : derived, derived=derived, area_m2=area, wet_mass_kg=wet, cd=Float64(v["cd"]))
end

function make_study_config(study; raan_deg::Real, omega_deg::Real)
    cfg = load_corridor_config(CORRIDOR_CONFIG_PATH)
    values = deepcopy(cfg.values)
    values["atmosphere"]["backend"] = "us76_piecewise"
    values["vehicle"]["ballistic_coefficient_kg_m2"] = study_beta(study, values).beta
    values["orbit"]["inclination_deg"] = Float64(study["grid"]["inclination_deg"])
    values["orbit"]["raan_deg"] = Float64(raan_deg)
    values["orbit"]["argument_of_periapsis_deg"] = Float64(omega_deg)
    values["orbit"]["initial_periapsis_altitude_m"] = Float64(study["calibration"]["initial_periapsis_altitude_guess_m"])
    return CorridorConfig(cfg.path, values)
end

function parse_epoch(s::AbstractString)
    d = DateTime(s)
    return InitialTime(year=Dates.year(d), month=Dates.month(d), day=Dates.day(d),
        hour=Dates.hour(d), minute=Dates.minute(d), second=Dates.second(d))
end

# Furnishes the SPICE kernels (de430.bsp first, pck00011.tpc, naif0012.tls) and returns the planet
# NBodyGravityModel carries. Call once per process before a third-body run.
nbody_planet() = Earth("", SPICE_PATH)

"""Truth `SimulationConfiguration`: the corridor's args (J2/J3 + US76 drag) plus, optionally, Sun and Moon."""
function build_args(cfg::CorridorConfig, study, epoch, horizon_s::Float64, dir::String;
                    third_bodies::Bool, state=nothing, mass_kg=nothing)
    args = state === nothing ?
        _corridor_args(cfg, :favorable, horizon_s, dir; epoch=epoch) :
        _corridor_args(cfg, :favorable, horizon_s, dir; initial_state=state, initial_mass_kg=mass_kg, epoch=epoch)
    # The engine refuses N-body gravity under SimpleEphemeridesModel, so every run of this study (with and
    # without third bodies, so the pair differs only by Sun+Moon) uses SpiceEphemeridesModel: ITRF93 Earth
    # orientation instead of the example's analytic GMST frame. Geodetic altitude is unaffected by longitude.
    env = args.environment_model
    spice_env = EnvironmentModel(planet=env.planet, EI=env.EI, density_model=env.density_model,
        thermal_model=env.thermal_model, ephemerides_model=SpiceEphemeridesModel(), topography=false, wind=false)
    effectors = args.dynamics_model.dynamic_effectors
    # Optional truth-gravity override (Earth >= 20x20 truth, estimator v3 study): a study config with a [gravity] table (degree, order,
    # coefficient_file) replaces the example's degree-3 order-0 harmonics model. Absent table: unchanged (the earlier truth runs).
    if haskey(study, "gravity")
        g = study["gravity"]; pl = args.environment_model.planet
        new_h = GravitationalHarmonicsModel(Int(g["degree"]), Int(g["order"]), joinpath(REPO_ROOT, String(g["coefficient_file"])), pl)
        effectors = Tuple(e isa GravitationalHarmonicsModel ? new_h : e for e in effectors)
    end
    if third_bodies
        nbody = NBodyGravityModel(body_names=Tuple(String.(study["third_bodies"]["bodies"])),
            primary_body_name="Earth", planet=nbody_planet())
        effectors = (effectors..., nbody)
    else
        nbody_planet()  # kernels for the Spice frame
    end
    return SimulationConfiguration(
        file_paths=args.file_paths, simulation_settings=args.simulation_settings,
        mission_configuration=args.mission_configuration, environment_model=spice_env,
        dynamics_model=DynamicsModel(args.dynamics_model.spacecraft, effectors),
        guidance_model=args.guidance_model, navigation_model=args.navigation_model,
        control_model=args.control_model, initial_time=args.initial_time,
        integration_tolerances=args.integration_tolerances, solver_config=args.solver_config)
end

# ---- heat rate --------------------------------------------------------------------------------

heat_rate_w_m2(rho, v_rel) = 0.5 * rho * v_rel^3

"""Geodetic altitude (m) where the US76 fit density equals `rho` (the fit extends its last layer above 150 km)."""
altitude_of_density(fit, rho) = find_zero(h -> log(corridor_density(fit, h)) - log(rho), (90_000.0, 250_000.0), Roots.Bisection())

"""`(h_at_q_high, h_at_q_low)`: altitude interval of the heat-rate band at relative speed `v_rel` (m/s)."""
function band_altitudes(fit, v_rel, study)
    b = study["band"]; f = Float64(b["w_cm2_to_w_m2"])
    rho_of(q) = 2 * q * f / v_rel^3
    return (altitude_of_density(fit, rho_of(Float64(b["q_high_w_cm2"]))),
            altitude_of_density(fit, rho_of(Float64(b["q_low_w_cm2"]))))
end

# ---- pass table -------------------------------------------------------------------------------

"""One row per radial-event periapsis: geodetic altitude, atmosphere-relative speed, density, peak Q, elements."""
function periapsis_table(sol, args, cfg::CorridorConfig, study)
    env = args.environment_model; planet = env.planet
    events = physical_apsis_events(sol, _geodetic_of(args))
    fit = fitted_atmosphere_parameters(cfg)
    et0 = SimulationModel.ephemerides_time_seconds(args.initial_time, env.ephemerides_model)
    f = Float64(study["band"]["w_cm2_to_w_m2"])
    rows = NamedTuple[]
    for ev in events
        ev.kind == "periapsis" || continue
        lpi = SimulationModel.planet_frame_lpi(planet, et0 + ev.time_s, env.ephemerides_model)
        vpp = lpi * ev.velocity_ii_mps - cross(planet.ω, lpi * ev.position_ii_m)
        vrel = norm(vpp)
        rho = corridor_density(fit, ev.physical_altitude_m)
        q = heat_rate_w_m2(rho, vrel) / f
        k = findfirst(c -> c.kind == "apoapsis" && c.time_s > ev.time_s, events)
        k === nothing && (k = findlast(c -> c.kind == "apoapsis" && c.time_s < ev.time_s, events))
        el = _corridor_element_summary(ev.position_ii_m, ev.velocity_ii_mps, planet.μ)
        push!(rows, (pass=length(rows) + 1, time_s=ev.time_s, alt_m=ev.physical_altitude_m, lat_rad=ev.latitude_rad,
            radius_m=ev.radius_m, vrel_mps=vrel, rho=rho, q_w_cm2=q,
            apo_alt_m=k === nothing ? NaN : events[k].physical_altitude_m,
            a_m=el.a, e=el.e, i_rad=el.i, raan_rad=el.raan, omega_rad=el.omega, mass_kg=ev.mass_kg))
    end
    return rows
end

# ---- perturber planes -------------------------------------------------------------------------

const MU_BODY = Dict("Moon" => 4.9028005821478e12, "Sun" => 1.32712440041279419e20)  # same table as NBodyGravityModel

"""
Sun/Moon plane angles of an orbit. The perturber plane normal is n = r_b x v_b (geocentric J2000, SPICE
spkezr at ephemeris time `et`). With h the orbit normal, k = n x h is the node of the orbit on the
perturber plane; omega' = atan2(h . (k x e_hat), k . e_hat) (apse angle from that node, positive along
the motion), i' = acos(h . n). Reduces to the usual omega, i for n = +z.
"""
function perturber_angles(body::String, et::Float64, i, raan, omega; center::String="EARTH", rot=nothing)
    st = spkezr(body, et, "J2000", "NONE", center)[1]
    rb = SVector{3, Float64}(st[1:3]) .* 1e3; vb = SVector{3, Float64}(st[4:6]) .* 1e3
    rot === nothing || ((rb, vb) = (SVector{3, Float64}(rot' * rb), SVector{3, Float64}(rot' * vb)))   # J2000 -> body-equator inertial frame
    n = normalize(cross(rb, vb))
    h = SVector(sin(i) * sin(raan), -sin(i) * cos(raan), cos(i))
    node = SVector(cos(raan), sin(raan), 0.0)
    ehat = cos(omega) * node + sin(omega) * cross(h, node)
    k = normalize(cross(n, h))
    omega_p = atan(dot(h, cross(k, ehat)), dot(k, ehat))
    return (omega_p=omega_p, i_p=acos(clamp(dot(h, n), -1.0, 1.0)), r_m=norm(rb), mu=MU_BODY[body])
end

"""Add `<body>_omega_p_rad, _i_p_rad, _kozai_drp_mps` columns (double-averaged r_p rate, -a*de/dt, hand formula)."""
function add_perturber_columns(rows, args, study)
    et0 = SimulationModel.ephemerides_time_seconds(args.initial_time, args.environment_model.ephemerides_model)
    out = NamedTuple[]
    for r in rows
        extra = (;)
        n_mean = sqrt(args.environment_model.planet.μ / r.a_m^3)
        for body in String.(study["third_bodies"]["bodies"])
            pa = perturber_angles(body, et0 + r.time_s, r.i_rad, r.raan_rad, r.omega_rad)
            rate = -r.a_m * (15 / 8) * (pa.mu / pa.r_m^3 / n_mean) * r.e * sqrt(1 - r.e^2) * sin(2 * pa.omega_p) * sin(pa.i_p)^2
            key = lowercase(body)
            extra = merge(extra, NamedTuple{(Symbol(key, "_omega_p_rad"), Symbol(key, "_i_p_rad"), Symbol(key, "_kozai_drp_mps"))}(
                (pa.omega_p, pa.i_p, rate)))
        end
        push!(out, merge(r, extra))
    end
    return out
end

# ---- running ----------------------------------------------------------------------------------

"""Run a truth case and return `(rows, sol, args, wall_s)`; `rows` is the periapsis table."""
function run_truth(cfg, study, epoch, horizon_s, dir; third_bodies, state=nothing, mass_kg=nothing, planes=false)
    args = build_args(cfg, study, epoch, Float64(horizon_s), dir; third_bodies, state, mass_kg)
    t0 = time()
    sol = run_simulation(args; return_solution=true, isolate_state=false)
    wall = time() - t0
    rows = periapsis_table(sol, args, cfg, study)
    planes && (rows = add_perturber_columns(rows, args, study))
    return (; rows, sol, args, wall_s=wall, retcode=string(sol.retcode))
end

"""Calibrate the seed periapsis altitude so the first-pass peak Q equals the band's geometric mean."""
function calibrate_peak_q(cfg::CorridorConfig, study, epoch, dir; third_bodies::Bool)
    cal = study["calibration"]; b = study["band"]
    q_target = sqrt(Float64(b["q_low_w_cm2"]) * Float64(b["q_high_w_cm2"]))
    working = CorridorConfig(cfg.path, deepcopy(cfg.values))
    fit = fitted_atmosphere_parameters(working)
    gm = _earth_radius_and_gm(working)[2]
    records = NamedTuple[]
    for it in 1:Int(cal["max_iterations"])
        o = working.values["orbit"]
        rp = Float64(o["initial_periapsis_altitude_m"]) + _earth_radius_and_gm(working)[1]
        ra = Float64(o["apoapsis_altitude_m"]) + _earth_radius_and_gm(working)[1]
        horizon = Float64(cal["horizon_fraction_of_period"]) * 2π * sqrt(((rp + ra) / 2)^3 / gm)
        res = run_truth(working, study, epoch, horizon, dir; third_bodies)
        isempty(res.rows) && throw(ErrorException("calibration run reached no periapsis"))
        r1 = res.rows[1]
        push!(records, (iteration=it, seed_alt_m=Float64(o["initial_periapsis_altitude_m"]), alt_m=r1.alt_m, vrel_mps=r1.vrel_mps, q_w_cm2=r1.q_w_cm2, q_target=q_target))
        abs(r1.q_w_cm2 / q_target - 1) <= Float64(cal["q_tolerance_relative"]) && return working, records
        h_target = altitude_of_density(fit, 2 * q_target * Float64(b["w_cm2_to_w_m2"]) / r1.vrel_mps^3)
        o["initial_periapsis_altitude_m"] = Float64(o["initial_periapsis_altitude_m"]) + (h_target - r1.alt_m)
    end
    throw(ErrorException("peak-Q calibration did not converge: $(records[end])"))
end

# ---- band statistics --------------------------------------------------------------------------

"""Band dwell statistics of a pass table (peak Q per pass against [q_low, q_high])."""
function band_stats(rows, study)
    qlo, qhi, qhard = Float64(study["band"]["q_low_w_cm2"]), Float64(study["band"]["q_high_w_cm2"]), Float64(study["band"]["q_hard_w_cm2"])
    isempty(rows) && return (; passes=0)
    inband = [qlo <= r.q_w_cm2 <= qhi for r in rows]
    dt = [j < length(rows) ? rows[j + 1].time_s - rows[j].time_s : 0.0 for j in eachindex(rows)]
    first_exit = findfirst(!, inband)
    edge = first_exit === nothing ? "" : (rows[first_exit].q_w_cm2 > qhi ? "high_Q(periapsis too low)" : "low_Q(periapsis too high)")
    return (; passes=length(rows), passes_in_band=count(inband), time_in_band_days=sum(dt[inband]; init=0.0) / SECONDS_PER_DAY,
        first_exit_time_days=first_exit === nothing ? NaN : rows[first_exit].time_s / SECONDS_PER_DAY, first_exit_edge=edge,
        q_max=maximum(r.q_w_cm2 for r in rows), q_min=minimum(r.q_w_cm2 for r in rows),
        exceeds_hard=any(r.q_w_cm2 > qhard for r in rows))
end

"""Linear interpolation of `col` of a pass table at time `t_s` (NaN outside the table)."""
function interp_col(rows, col::Symbol, t_s)
    (isempty(rows) || t_s < rows[1].time_s || t_s > rows[end].time_s) && return NaN
    j = searchsortedlast([r.time_s for r in rows], t_s)
    j >= length(rows) && return Float64(getfield(rows[end], col))
    w = (t_s - rows[j].time_s) / (rows[j + 1].time_s - rows[j].time_s)
    return (1 - w) * getfield(rows[j], col) + w * getfield(rows[j + 1], col)
end

"""Process-pool helper shared by the grid and switchability drivers: run `task` on `seeds` via run_monte_carlo, process route."""
function run_process_campaign(task, seeds, workers::Int; mission_time_s::Float64, worker_files::Vector{String})
    pool = campaign_process_pool()
    ids = ensure_process_workers!(pool, workers)
    Distributed.remotecall_eval(Main, ids, :(ENV["PC_SWEEP_WORKER"] = "1"))
    for f in worker_files
        Distributed.remotecall_eval(Main, ids, :(include($f)))
    end
    # density_family="gram_point" is the router feature whose candidate set is [:none, :process]
    # (the same device the corridor sweep uses): SPICE is not thread-safe, so no thread route.
    features = SpaceAGORA.SimulationCampaigns.campaign_route_features(samples=length(seeds), density_family="gram_point", mission_time_s=mission_time_s)
    tuning = SpaceAGORA.SimulationCampaigns.OuterRouteTuning(adaptive_enabled=false, process_workers_resident=length(ids), process_max_workers=workers)
    res = run_monte_carlo(task, seeds; threads=:auto, route_features=features, route_state=SpaceAGORA.SimulationCampaigns.OuterRouteState(), route_tuning=tuning)
    res.route === :process || throw(ErrorException("campaign ran on route $(res.route), not :process"))
    return res
end
