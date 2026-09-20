# Magellan aerobraking at Venus, two orbits from the start of the campaign
# (1993-05-25 to 1993-08-16), with the real trajectory from the Magellan
# navigation SPK (NAIF MGN archive, AEROBRAK.BSP) drawn as a ghost.
#
# Initial state: the first apoapsis in the kernel after 1993-05-26 00:00 UTC,
# straight from SPICE (J2000, Venus-centerd). Force model: MGNP180U gravity
# to degree and order 70, Sun third body, solar radiation pressure, Venus-GRAM
# density through the free-molecular coefficient model. The spacecraft is the
# usual bus-plus-two-wings composition sized to Magellan (3.7 m HGA, two
# 2.5 m square wings, about 1000 kg at this point of the mission) with its
# long (bus + antenna) axis along the flow and the wings broadside (the NASA
# model's canted cruise wings are articulated to that pose for both the mesh
# aerodynamics and the picture). The bus
# and antenna end leads by default (model +Y, the antenna axis, to body +x);
# SPACEAGORA_DEMO_MAGELLAN_ANTENNA=aft flies it the other way round, antenna
# trailing, into `<case>_antenna_aft`. The box model has no fore/aft
# distinction, so only the mesh surrogate and the picture change.
#
#   julia --project=. scripts/dev/viewer_demos/magellan_aerobraking.jl
include(joinpath(@__DIR__, "common.jl"))
function _mission_main(argv)

antenna = get(ENV, "SPACEAGORA_DEMO_MAGELLAN_ANTENNA", "forward")
antenna in ("forward", "aft") || throw(ArgumentError("antenna option must be forward or aft"))
ANTENNA_AFT = antenna == "aft"
case_name = demo_case_name("magellan_aerobraking") * (ANTENNA_AFT ? "_antenna_aft" : "")
MODEL = demo_model_path("magellan_nasa_3d_resources.glb")
# antenna axis (model +Y) along the flow: to body +x (bus and antenna forward) or -x (antenna aft); wings along body y
MAGELLAN_ROTATION_DEG = ANTENNA_AFT ? (0, 0, 90) : (-180, 0, -90)
# The NASA model holds the wings in a cruise pose, canted 43.5 deg between the
# antenna axis (model Y) and model Z (area-weighted wing normal (0.01, 0.70, 0.66)).
# For aerobraking Magellan turned them broadside to the flow: rotate each wing
# about its own boom (model X) so its normal lies along model Y, the flow axis.
# Regions are in model units before the viewer transform; y_max keeps the
# boom struts (at x < -1.9, y >= 1) out of the wing region.
MAGELLAN_WING_TILT_DEG = -43.5
MAGELLAN_ARTICULATIONS = (
    (region=(x_min=1.9, y_max=1.0), axis=(1.0, 0.0, 0.0), angle_deg=MAGELLAN_WING_TILT_DEG),
    (region=(x_max=-1.9, y_max=1.0), axis=(1.0, 0.0, 0.0), angle_deg=MAGELLAN_WING_TILT_DEG),
)
kernel = ensure_mission_kernel("mgn_aerobrak.bsp", "https://naif.jpl.nasa.gov/pub/naif/MGN/kernels/spk/nav/AEROBRAK.BSP")
planet = Venus("", SPICE_PATH)
mgn = furnish!(MissionSpice("Magellan (SPICE)", "-18", "VENUS", [kernel]))

et_apo = first_apoapsis_et_after(mgn, et_of("1993-05-26T00:00:00"))
initial_time = initial_time_of(et_apo)
et0 = et_of(initial_time) # Match the engine clock and the saved SPICE reference.
r0, v0 = spice_state_m(mgn, et0)
a0 = 1.0 / (2.0 / norm(r0) - dot(v0, v0) / planet.μ)
period = 2pi * sqrt(a0^3 / planet.μ)
println("Magellan apoapsis epoch ", utc_of(et_apo), "  r=", round(norm(r0) / 1e3; digits=1), " km  a=", round(a0 / 1e3; digits=1), " km  period=", round(period / 60; digits=1), " min")
options = viewer_demo_options(case_name, 2.0 * period; argv=argv)
OUTDIR = options.output_dir
mission_time = options.duration_s

sc = make_three_body_spacecraft(
    bus_dims=(3.6, 3.0, 3.0), panel_dims=(0.01, 2.5, 2.5), bus_mass=960.0, panel_mass_each=35.0, panel_offset_y=3.2,
    ic=cartesian_ic_at(mgn, et0), reflection_coefficient=0.9, prop_mass=30.0, id=1)
effectors = (
    GravitationalHarmonicsModel(70, 70, joinpath(HARMONICS_DIR, "MGNP180U.csv"), planet),
    NBodyGravityModel(body_names=("Sun",), primary_body_name="Venus", planet=planet),
    SolarRadiationPressureModel(1.3, 23.0),
    demo_aero_effector(MODEL, OUTDIR; scale=1.0, rotation_deg=MAGELLAN_ROTATION_DEG, wall_temperature_k=300.0, articulations=MAGELLAN_ARTICULATIONS),
)
base = make_example_config(planet=planet, spacecraft=sc, mission_time=mission_time, initial_time=initial_time,
    dynamic_effectors=effectors, density_model=GRAMAtmosphereModel(planet_name="venus", initial_time=initial_time), orientation_sim=false,
    keplerian=true, EI_km=250.0, verbose=false, results=true, results_directory=OUTDIR)
args = SM.SimConfig._with_configuration(base;
    mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, keplerian=true, number_of_orbits=1,
        mission_time=mission_time, orientation_sim=false, num_steps_to_save=1000, data_rate=5.0),
    solver_config=SM.SolverConfig(solver_mode=:tsit5),
    integration_tolerances=SM.IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=20.0,
        reltol_atmosphere=1e-9, abstol_atmosphere=1e-9, dt_max_atmosphere=0.5))

prefix = run_or_reuse!(args, OUTDIR)
summarize_run(prefix, planet; alt_m=250e3)
ghost = spice_reference(mgn, prefix; name="Magellan (SPICE)", color="#ff8c69")
reference_separation(prefix, ghost)

model = MODEL
html = export_visualization(prefix; max_frames=5000, trail_orbits=1, texture_resolution="4k",
    title="AGORA Magellan · aerobraking at Venus, modeled trajectory with the SPICE reference" * (get(ENV, "SPACEAGORA_DEMO_MESH_AERO", "0") == "1" ? " (mesh aerodynamics)" : ""),
    models=Dict(1 => model), model_scale=1.0, model_rotation_deg=Dict(1 => MAGELLAN_ROTATION_DEG), model_articulations=Dict(1 => MAGELLAN_ARTICULATIONS), references=[ghost])
println("html: ", html, " ", filesize(html))
cdn = build_cdn_page(html, joinpath(OUTDIR, "artifact.html"), "AGORA Magellan Aerobraking",
    "AGORA Magellan · aerobraking at Venus, from 1993-05-26",
    "SPICE state at apoapsis; MGNP180U 70x70, Sun, SRP, Venus-GRAM drag" * (get(ENV, "SPACEAGORA_DEMO_MESH_AERO", "0") == "1" ? " on the CAD mesh" : " on the box model") * (ANTENNA_AFT ? "; antenna aft" : "; bus and antenna forward"), "$(round(mission_time / 3600; digits=3)) h modeled span",
    "The solid spacecraft is the simulation; the translucent copy follows the Magellan navigation SPK (AEROBRAK.BSP). Click the solid one to read the separation. Drag to orbit, wheel to zoom, Space to pause, F to follow.")
cdn === nothing || println("cdn: ", cdn, " ", filesize(cdn))

return (args=args, prefix=prefix, html=html, cdn=cdn)

end # main

function main(argv=ARGS)
    setup_gram_example!()
    return Base.invokelatest(_mission_main, argv)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
