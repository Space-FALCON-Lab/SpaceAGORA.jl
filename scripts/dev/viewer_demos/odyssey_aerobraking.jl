# 2001 Mars Odyssey aerobraking, two orbits from the main phase of the
# campaign (2001-10-26 to 2002-01-11), with the reconstructed navigation
# trajectory from the NAIF ODY archive (m01_ab_v2.bsp, the kernel behind the
# telemetry validation study) drawn as a ghost.
#
# Initial state: the first apoapsis in the kernel after 2001-11-06 00:00 UTC
# (the epoch of the reconstruction record, orbit P19/P20, periapsis near
# 100 km and a 17 h period), straight from SPICE (J2000, Mars-centered).
# Force model: Mars-50c gravity to degree and order 50, Sun third body,
# solar radiation pressure, Mars-GRAM (TES mapping year 2) density through
# the free-molecular coefficient model. The spacecraft is the reconstruction
# record's bus-plus-two-wings composition (2.2 x 2.6 x 1.7 m bus, 1.95 x
# 1.7 m wings beside it, 461 kg), 11 m² broadside to the flow, the array on
# the +y side of the bus. NASA's Odyssey model is posed the same way: the
# array normal along the flow (model Z to body -x, cells away from the
# flow), the array offset along body +y, and the gamma-ray spectrometer boom
# trailing along body -x. The model shows the boom and the high-gain antenna
# deployed; both were stowed through aerobraking and deployed in 2002.
#
# Mars-GRAM's TES mapping-year-2 climatology gives about half the density
# Odyssey's accelerometers measured on these passes (the 2001 dust storm was
# still decaying), so the simulation loses about half the energy the kernel
# does and trails the ghost by a thousand kilometers after two orbits.
# SPACEAGORA_DEMO_ODYSSEY_ATMOSPHERE=accelerometer replaces the climatology
# with the accelerometer-derived per-pass density profiles the telemetry
# validation study uses (data/telemetry/Odyssey/), into
# `<case>_accelerometer/`, and the separation drops to the kilometer level.
#
#   julia --project=. scripts/dev/viewer_demos/odyssey_aerobraking.jl
include(joinpath(@__DIR__, "common.jl"))
setup_gram_example!()

const ACCELEROMETER = get(ENV, "SPACEAGORA_DEMO_ODYSSEY_ATMOSPHERE", "gram") == "accelerometer"
const OUTDIR = demo_outdir(demo_case_name("odyssey_aerobraking") * (ACCELEROMETER ? "_accelerometer" : ""))
const DENSITY_TABLE = joinpath(REPO_ROOT, "data", "telemetry", "Odyssey", "odyssey_accelerometer_density.feather")
const MODEL = joinpath(MODELS_DIR, "mars_odyssey_nasa_3d_resources.glb")
# array normal (model Z) along the flow with the boom trailing (model +Z to body -x); array offset (model +Y) to body +y
const ODYSSEY_ROTATION_DEG = (0, -90, 0)
kernel = ensure_mission_kernel("m01_ab_v2.bsp", "https://naif.jpl.nasa.gov/pub/naif/pds/data/ody-m-spice-6-v1.0/odsp_1000/data/spk/m01_ab_v2.bsp")
planet = Mars("", SPICE_PATH)
ody = furnish!(MissionSpice("Odyssey (SPICE)", "-53", "MARS", [kernel]))

et_apo = first_apoapsis_et_after(ody, et_of("2001-11-06T00:00:00"))
r0, v0 = spice_state_m(ody, et_apo)
a0 = 1.0 / (2.0 / norm(r0) - dot(v0, v0) / planet.μ)
period = 2pi * sqrt(a0^3 / planet.μ)
println("Odyssey apoapsis epoch ", utc_of(et_apo), "  r=", round(norm(r0) / 1e3; digits=1), " km  a=", round(a0 / 1e3; digits=1), " km  period=", round(period / 60; digits=1), " min")
initial_time = initial_time_of(et_apo)
mission_time = 2.0 * period
density_model = if ACCELEROMETER
    TV = SpaceAGORA.TelemetryVerification
    TV._make_tabulated_flight_density_model(initial_time, TV.AtmosphereTruthConfig(assumption_id="odyssey_accelerometer_demo",
        atmosphere_model="tabulated_flight", atmosphere_dataset="Odyssey accelerometer per-pass density profiles (PDS odya_0001, Tolson)",
        tabulated_flight_file=DENSITY_TABLE, tabulated_flight_sigma=0.0))
else
    GRAMAtmosphereModel(planet_name="mars", mars_map_year=2)
end

sc = make_three_body_spacecraft(
    bus_dims=(2.2, 2.6, 1.7), panel_dims=(0.01, 1.945, 1.7), bus_mass=391.0, panel_mass_each=10.0, panel_offset_y=2.2725,
    ic=cartesian_ic_at(ody, et_apo), reflection_coefficient=0.9, prop_mass=50.0, id=1)
effectors = (
    GravitationalHarmonicsModel(50, 50, joinpath(HARMONICS_DIR, "Mars50c.csv"), planet),
    NBodyGravityModel(body_names=("Sun",), primary_body_name="Mars", planet=planet),
    SolarRadiationPressureModel(1.3, 11.0),
    demo_aero_effector(MODEL, OUTDIR; scale=1.0, rotation_deg=ODYSSEY_ROTATION_DEG, wall_temperature_k=300.0),
)
base = make_example_config(planet=planet, spacecraft=sc, mission_time=mission_time, initial_time=initial_time,
    dynamic_effectors=effectors, density_model=density_model, orientation_sim=false,
    keplerian=true, EI_km=250.0, verbose=false, results=true, results_directory=OUTDIR)
args = SM.SimulationConfiguration(file_paths=base.file_paths, simulation_settings=base.simulation_settings,
    mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, keplerian=true, number_of_orbits=1,
        mission_time=mission_time, orientation_sim=false, num_steps_to_save=1000, data_rate=5.0),
    environment_model=base.environment_model, dynamics_model=base.dynamics_model, guidance_model=base.guidance_model,
    navigation_model=base.navigation_model, control_model=base.control_model, initial_time=base.initial_time,
    integration_tolerances=SM.IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=20.0,
        reltol_atmosphere=1e-9, abstol_atmosphere=1e-9, dt_max_atmosphere=0.5))

prefix = run_or_reuse!(args, OUTDIR)
summarize_run(prefix, planet; alt_m=250e3)
ghost = spice_reference(ody, prefix; name="Odyssey (SPICE)", color="#ff8c69")
reference_separation(prefix, ghost)

mesh = get(ENV, "SPACEAGORA_DEMO_MESH_AERO", "0") == "1"
atmosphere = ACCELEROMETER ? "accelerometer-derived density" : "Mars-GRAM"
html = export_visualization(prefix; max_frames=5000, trail_orbits=1, texture_resolution="4k",
    title="AGORA Odyssey · aerobraking at Mars, two orbits with the SPICE ghost" * (ACCELEROMETER ? ", accelerometer density" : "") * (mesh ? " (mesh aerodynamics)" : ""),
    models=Dict(1 => MODEL), model_scale=1.0, model_rotation_deg=Dict(1 => ODYSSEY_ROTATION_DEG), references=[ghost])
println("html: ", html, " ", filesize(html))
cdn = build_cdn_page(html, joinpath(OUTDIR, "artifact.html"), "AGORA Odyssey Aerobraking",
    "AGORA Odyssey · aerobraking at Mars, two orbits from 2001-11-06",
    "SPICE state at apoapsis; Mars-50c 50x50, Sun, SRP, " * atmosphere * " drag" * (mesh ? " on the CAD mesh" : " on the box model") * "; array broadside, boom trailing",
    "2 orbits (≈$(round(mission_time / 3600; digits=1)) h)",
    "The solid spacecraft is the simulation; the translucent copy follows the Odyssey navigation kernel (m01_ab_v2.bsp). Click the solid one to read the separation. Drag to orbit, wheel to zoom, Space to pause, F to follow.")
println("cdn: ", cdn, " ", filesize(cdn))
