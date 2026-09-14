# Cassini's Titan flybys TA (2004-10-26, first targeted flyby, Huygens still
# attached, closest approach ~1174 km) and T5 (2005-04-16, ~1027 km, the low
# pass with measurable drag), Titan-centred, with the reconstructed
# trajectory from the Cassini SCPSE kernels drawn as a ghost.
#
# The simulation starts 2.5 h before closest approach from the SPICE state
# (J2000, Titan-centred) and runs 5 h. Force model: Titan gravity to degree
# and order 5 (Goossens et al. 2024 field), Saturn and Sun third bodies,
# solar radiation pressure, Titan-GRAM density through the free-molecular
# coefficient model. Cassini is one 4 x 4 x 6.8 m bus (the wing links are
# vestigial); the page draws NASA's Cassini-Huygens model over it, with or
# without the probe depending on the flyby.
#
#   julia --project=. scripts/dev/viewer_demos/cassini_titan_flyby.jl [TA|T5|all]
include(joinpath(@__DIR__, "common.jl"))
setup_gram_example!()

const CASSINI_SPK = "https://naif.jpl.nasa.gov/pub/naif/CASSINI/kernels/spk/"
const FLYBYS = Dict(
    "TA" => (kernel="050105R_SCPSE_04247_04336.bsp", ca_guess="2004-10-26T15:30:00", model="cassini_huygens_nasa_3d_resources_a.glb",
             mass=4600.0, note="first targeted Titan flyby, Huygens attached"),
    "T5" => (kernel="050513R_SCPSE_05097_05114.bsp", ca_guess="2005-04-16T19:12:00", model="cassini_nasa_3d_resources_a_without_huygens.glb",
             mass=4200.0, note="low flyby with measurable drag, Huygens released"),
)

function run_flyby(key::String; half_span_s::Float64=2.5 * 3600.0)
    fb = FLYBYS[key]
    outdir = demo_outdir("cassini_$(lowercase(key))")
    kernel = ensure_mission_kernel(fb.kernel, CASSINI_SPK * fb.kernel)
    planet = Titan("", SPICE_PATH)
    cas = furnish!(MissionSpice("Cassini (SPICE)", "CASSINI", "TITAN", [kernel]))
    et_ca = closest_approach_et(cas, et_of(fb.ca_guess); window_s=2.0 * 3600.0)
    r_ca, v_ca = spice_state_m(cas, et_ca)
    println(key, ": closest approach ", utc_of(et_ca), "  altitude ", round((norm(r_ca) - planet.Rp_e) / 1e3; digits=1),
        " km  speed ", round(norm(v_ca) / 1e3; digits=3), " km/s")
    et0 = et_ca - half_span_s
    initial_time = initial_time_of(et0)
    mission_time = 2.0 * half_span_s

    sc = make_three_body_spacecraft(
        bus_dims=(4.0, 4.0, 6.8), panel_dims=(0.01, 0.02, 0.02), bus_mass=fb.mass - 200.0, panel_mass_each=0.5, panel_offset_y=2.1,
        ic=cartesian_ic_at(cas, et0), reflection_coefficient=0.8, prop_mass=200.0, id=1)
    effectors = (
        GravitationalHarmonicsModel(5, 5, joinpath(HARMONICS_DIR, "titan5.csv"), planet),
        NBodyGravityModel(body_names=("Saturn", "Sun"), primary_body_name="Titan", planet=planet),
        SolarRadiationPressureModel(1.2, 20.0),
        AerodynamicCoefficientfM(),
    )
    base = make_example_config(planet=planet, spacecraft=sc, mission_time=mission_time, initial_time=initial_time,
        dynamic_effectors=effectors, density_model=GRAMAtmosphereModel(planet_name="titan"), orientation_sim=false,
        keplerian=true, EI_km=1500.0, verbose=false, results=true, results_directory=outdir)
    args = SM.SimulationConfiguration(file_paths=base.file_paths, simulation_settings=base.simulation_settings,
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=mission_time, orientation_sim=false, num_steps_to_save=1000, data_rate=5.0),
        environment_model=base.environment_model, dynamics_model=base.dynamics_model, guidance_model=base.guidance_model,
        navigation_model=base.navigation_model, control_model=base.control_model, initial_time=base.initial_time,
        integration_tolerances=SM.IntegrationTolerances(reltol_orbit=1e-10, abstol_orbit=1e-9, dt_max_orbit=10.0,
            reltol_atmosphere=1e-10, abstol_atmosphere=1e-9, dt_max_atmosphere=0.5))

    prefix = run_or_reuse!(args, outdir)
    summarize_run(prefix, planet; alt_m=1500e3)
    ghost = spice_reference(cas, prefix; name="Cassini (SPICE)", color="#7fe0ff")
    reference_separation(prefix, ghost)

    model = joinpath(MODELS_DIR, fb.model)
    html = export_visualization(prefix; max_frames=4000, trail_orbits=1, texture_resolution="4k",
        title="AGORA Cassini · Titan flyby $(key) with the SPICE ghost",
        models=Dict(1 => model), model_scale=1.0, model_rotation_deg=Dict(1 => (-90, -180, 0)), references=[ghost])
    println("html: ", html, " ", filesize(html))
    cdn = build_cdn_page(html, joinpath(outdir, "artifact.html"), "AGORA Cassini $(key) Flyby",
        "AGORA Cassini · Titan flyby $(key), $(fb.note)",
        "SPICE state 2.5 h before closest approach; Titan 5x5 field, Saturn + Sun, SRP, Titan-GRAM drag",
        "5 h around closest approach ($(utc_of(et_ca)[1:16])Z)",
        "The solid spacecraft is the simulation; the translucent copy follows the reconstructed Cassini SPK. Click the solid one to read the separation. Drag to orbit, wheel to zoom, Space to pause, F to follow.")
    println("cdn: ", cdn, " ", filesize(cdn))
    return (prefix=prefix, html=html, cdn=cdn, et_ca=et_ca)
end

which = isempty(ARGS) ? "TA" : uppercase(ARGS[1])
for key in (which == "ALL" ? ("TA", "T5") : (which,))
    haskey(FLYBYS, key) || throw(ArgumentError("unknown flyby $(key); use TA, T5 or all"))
    run_flyby(key)
end
