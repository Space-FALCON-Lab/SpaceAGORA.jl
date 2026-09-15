# Apollo 11 powered descent to Tranquility Base, 1969-07-20, from powered
# descent initiation (PDI, 20:05:05 UTC, 15.24 km above the site, 480 km
# uprange on the 15 x 105 km descent orbit, flying west at 1693 m/s) to
# touchdown, with the LM descent stage flown 6-DOF: the Apollo quadratic
# guidance (P63 braking, P64 approach, P66 rate of descent) commands the
# descent engine's thrust and the attitude that points it, an RCS attitude
# controller tracks that attitude, and the thrust acts along the vehicle's
# actual engine axis. The terrain is the LROC NAC DTM of the site inside the
# LOLA global grid (fetched by scripts/dev/terrain/fetch_moon_site.py): it
# is the radar altimeter of the guidance and the ground the touchdown event
# fires on, and the page drapes the LROC imagery over it, sharpening from
# the WAC mosaic to the 26 cm NAC mosaic as the lander approaches.
#
#   python3 scripts/dev/terrain/fetch_moon_site.py --site 0.67416 23.47314 --out data/terrain/moon/apollo11
#   julia --project=. scripts/dev/viewer_demos/apollo11_landing.jl
include(joinpath(@__DIR__, "common.jl"))
using Arrow, DataFrames

const OUTDIR = demo_outdir("apollo11_landing")
const SITE_JSON = joinpath(REPO_ROOT, "data", "terrain", "moon", "apollo11", "site.json")
const MODEL = joinpath(MODELS_DIR, "apollo_lunar_module_nasa_3d_resources.glb")
const LM_ROTATION_DEG = (-180, 90, 90)   # legs along body +z (down through the engine), hatch and windows along body +x
isfile(SITE_JSON) || error("Site terrain not found at $(SITE_JSON); run scripts/dev/terrain/fetch_moon_site.py first.")

planet = Moon("", SPICE_PATH)
terrain, site = load_site_terrain(SITE_JSON)
println("site ", site.name, " at ", site.lat_deg, "N ", site.lon_deg, "E, ground ", round(site.height_m; digits=1), " m relative to the ", site.reference_radius_m / 1e3, " km sphere")

# ---- initial state at PDI --------------------------------------------------------------
et0 = et_of("1969-07-20T20:05:05")
initial_time = initial_time_of(et0)
const PDI_ALTITUDE_M = 15_240.0      # 50,000 ft above the site
const PDI_UPRANGE_M = 480.0e3       # along the approach, east of the site
const PDI_SPEED_MPS = 1_693.0       # perilune speed of the 15.24 x 105 km descent orbit
const APPROACH_AZIMUTH_DEG = 270.0  # flying west
let
    global ic, q_pdi
    az = deg2rad(APPROACH_AZIMUTH_DEG)
    # site axes in planet-fixed coordinates
    φ = deg2rad(site.lat_deg); λ = deg2rad(site.lon_deg)
    up_s = SVector(cos(φ) * cos(λ), cos(φ) * sin(λ), sin(φ))
    east_s = SVector(-sin(λ), cos(λ), 0.0)
    north_s = SVector(-sin(φ) * cos(λ), -sin(φ) * sin(λ), cos(φ))
    flight = normalize(cos(az) * north_s + sin(az) * east_s)
    # PDI point: rotate the site's up vector back along the approach by the uprange arc
    θ = PDI_UPRANGE_M / planet.Rp_e
    up_pdi = normalize(cos(θ) * up_s - sin(θ) * flight)
    flight_pdi = normalize(flight - dot(flight, up_pdi) * up_pdi)   # horizontal at PDI, still toward the site
    r_p = (planet.Rp_e + site.height_m + PDI_ALTITUDE_M) * up_pdi
    v_p = PDI_SPEED_MPS * flight_pdi                                # relative to the rotating Moon
    r_i, v_i = SM.GuidanceHooks.r_pintor_i(r_p, v_p, planet, et0)
    # engine retrograde with the windows up: the attitude the trim phase starts from
    q_pdi = descent_attitude_command(-normalize(v_i), normalize(r_i), normalize(v_i))
    ic = SM.CartesianInitialCondition(r_i, v_i; q=q_pdi)
    println("PDI at ", utc_of(et0), ": r=", round(norm(r_i) / 1e3; digits=2), " km, v=", round(norm(v_i); digits=1), " m/s inertial")
end

# ---- the lunar module ------------------------------------------------------------------
# Descent stage + ascent stage, ~15.1 t at PDI with 8.2 t of DPS propellant. Body frame:
# +z down through the descent engine, +x out of the windows. Inertia is the order of the
# LM's at PDI; thrusters are drawn by the viewer (the effector applies the forces).
dps = SM.Thruster(max_thrust=45_040.0, location=MVector{3, Float64}(0.0, 0.0, 1.5), direction=MVector{3, Float64}(0.0, 0.0, 1.0), Isp=311.0)
rcs = SM.Thruster[]
for sx in (-1.0, 1.0), sy in (-1.0, 1.0)
    quad = MVector{3, Float64}(1.65 * sx, 1.65 * sy, -0.3)
    for d in ((sx, 0.0, 0.0), (0.0, sy, 0.0), (0.0, 0.0, 1.0), (0.0, 0.0, -1.0))
        push!(rcs, SM.Thruster(max_thrust=445.0, location=copy(quad), direction=MVector{3, Float64}(d...), Isp=290.0))
    end
end
root = SM.Link(root=true, m=6_900.0, dims=MVector{3, Float64}(4.2, 4.2, 7.0), ref_area=1.0, q=MVector{4, Float64}(q_pdi...), thrusters=[dps; rcs])
sc = SM.SpacecraftModel(links=[root], root=root, prop_mass=8_200.0,
    inertia_tensor=SMatrix{3, 3, Float64}(2.2e4, 0, 0, 0, 2.4e4, 0, 0, 0, 2.0e4), initial_condition=ic, id=1)

# ---- guidance and control ---------------------------------------------------------------
braking, approach = apollo11_descent_targets()
gcfg = ApolloDescentConfig(site_lat_deg=site.lat_deg, site_lon_deg=site.lon_deg, approach_azimuth_deg=APPROACH_AZIMUTH_DEG,
    braking=braking, approach=approach)
state = ApolloDescentState(1)
guidance = ApolloDescentGuidanceModel(gcfg, state, terrain)
# the state point is the model's center, 3.7 m above the footpads in NASA's model at this scale
control = ApolloDescentControlModel(ApolloDescentControlConfig(touchdown_height_m=3.7), gcfg, state, terrain)

effectors = (
    GravitationalHarmonicsModel(50, 50, joinpath(HARMONICS_DIR, "LP165P.csv"), planet),
    NBodyGravityModel(body_names=("Earth", "Sun"), primary_body_name="Moon", planet=planet),
)
mission_time = 1_000.0   # the touchdown event ends the run earlier
base = make_example_config(planet=planet, spacecraft=sc, mission_time=mission_time, initial_time=initial_time,
    dynamic_effectors=effectors, density_model=NoAtmosphereModel(), orientation_sim=true,
    keplerian=false, EI_km=1.0, verbose=false, results=true, results_directory=OUTDIR)
args = SM.SimulationConfiguration(file_paths=base.file_paths, simulation_settings=base.simulation_settings,
    mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, keplerian=false, number_of_orbits=1,
        mission_time=mission_time, orientation_sim=true, num_steps_to_save=4000, data_rate=0.25),
    environment_model=base.environment_model, dynamics_model=base.dynamics_model,
    guidance_model=SM.GuidanceModel(guidance_effectors=(guidance,), guidance_rates=[2.0]),
    navigation_model=base.navigation_model,
    control_model=SM.ControlModel(control_effectors=(control,), control_rates=[0.05]),
    initial_time=base.initial_time,
    integration_tolerances=SM.IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=0.25,
        reltol_atmosphere=1e-9, abstol_atmosphere=1e-9, dt_max_atmosphere=0.25))

# Guidance diagnostics saved beside the state: thrust, throttle, time-to-go, radar altitude, phase.
phase_code(s::Symbol) = s === :braking ? 1.0 : s === :approach ? 2.0 : s === :vertical ? 3.0 : 4.0
save_fields = [
    SM.SimulationCallbacks.default_save_fields(args)...,
    SM.SimulationCallbacks.SaveField(:dps_thrust_n, (u, t, integ) -> [control.actuators.thrust_n[1]]; per_satellite=true),
    SM.SimulationCallbacks.SaveField(:throttle, (u, t, integ) -> [state.throttle[1]]; per_satellite=true),
    SM.SimulationCallbacks.SaveField(:time_to_go_s, (u, t, integ) -> [state.t_go_s[1]]; per_satellite=true),
    SM.SimulationCallbacks.SaveField(:radar_altitude_m, (u, t, integ) -> [state.radar_altitude_m[1]]; per_satellite=true),
    SM.SimulationCallbacks.SaveField(:guidance_phase, (u, t, integ) -> [phase_code(state.phase[1])]; per_satellite=true),
    SM.SimulationCallbacks.SaveField(:attitude_error_deg, (u, t, integ) -> [rad2deg(control.actuators.attitude_error_rad[1])]; per_satellite=true),
]
prefix = run_or_reuse!(args, OUTDIR; save_fields=save_fields, isolate_state=false)

df = DataFrame(Arrow.Table(prefix * ".feather"))
last = df[end, :]
println("rows=", nrow(df), "  span=", round(last.time; digits=1), " s (", round(last.time / 60; digits=2), " min)")
if isfinite(state.touchdown_s[1])
    v = state.touchdown_v_mps[1]
    println("touchdown at t=", round(state.touchdown_s[1]; digits=1), " s: vertical ", round(v[3]; digits=2), " m/s, horizontal ", round(hypot(v[1], v[2]); digits=2),
        " m/s, ", round(state.touchdown_miss_m[1]; digits=1), " m from the target; propellant used ", round(15_100.0 - last.sc1_mass; digits=0), " kg")
    for (k, name) in enumerate(("braking", "approach", "vertical", "landed"))
        isfinite(state.phase_start_s[k, 1]) && println("  ", name, " from t=", round(state.phase_start_s[k, 1]; digits=1), " s")
    end
else
    println("final: radar altitude ", round(last.sc1_radar_altitude_m; digits=1), " m, phase ", last.sc1_guidance_phase, " (results reused or no touchdown)")
end

html = export_visualization(prefix; max_frames=4000, trail_orbits=1, texture_resolution="4k",
    title="AGORA Apollo 11 · powered descent to Tranquility Base",
    models=Dict(1 => MODEL), model_scale=1.45, model_rotation_deg=Dict(1 => LM_ROTATION_DEG), terrain=SITE_JSON)
println("html: ", html, " ", filesize(html))
cdn = build_cdn_page(html, joinpath(OUTDIR, "artifact.html"), "AGORA Apollo 11 Landing",
    "AGORA Apollo 11 · powered descent to Tranquility Base, 1969-07-20 20:05 UTC",
    "PDI 15.24 km above the site, 480 km uprange; LP165P 50x50, Earth + Sun; quadratic guidance P63/P64/P66, DPS 10-60% + full, RCS attitude control",
    "PDI to touchdown (≈$(round(last.time / 60; digits=1)) min)",
    "The LM flies 6-DOF; the terrain is the LROC NAC DTM inside the LOLA grid, draped with LROC WAC and NAC imagery down to 0.65 m/px. Click values for their history, click the LM for a face, F follows.")
println("cdn: ", cdn, " ", filesize(cdn))
