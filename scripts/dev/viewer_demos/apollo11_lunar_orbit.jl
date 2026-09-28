# The Apollo 11 lunar module in the parking orbit at undocking
# (1969-07-20 17:44 UTC): about 101 x 122 km, retrograde, 1.25 deg off the
# lunar equator, two orbits. No Apollo navigation kernel is supplied here; the nominal orbit is
# an illustrative reconstruction: elements in the Moon's body-fixed frame at the
# epoch, rotated into J2000 with SPICE. Force model: LP165P gravity to
# degree and order 50, Earth and Sun third bodies, solar radiation pressure,
# no atmosphere. The page draws NASA's Apollo Lunar Module model.
#
#   julia --project=. scripts/dev/viewer_demos/apollo11_lunar_orbit.jl
include(joinpath(@__DIR__, "common.jl"))


# Orbit in a body-fixed lunar frame at `et`, returned as a J2000 Cartesian initial condition.
function lunar_orbit_ic(planet, et::Float64; rp_alt_m, ra_alt_m, inc_deg, raan_deg, argp_deg, ta_deg, frame="IAU_MOON")
    rp = planet.Rp_e + rp_alt_m; ra = planet.Rp_e + ra_alt_m
    a = (rp + ra) / 2; e = (ra - rp) / (ra + rp); p = a * (1 - e^2)
    ν = deg2rad(ta_deg)
    r_pf = p / (1 + e * cos(ν)) .* SVector(cos(ν), sin(ν), 0.0)
    v_pf = sqrt(planet.μ / p) .* SVector(-sin(ν), e + cos(ν), 0.0)
    Rz(θ) = SMatrix{3, 3}(cos(θ), sin(θ), 0.0, -sin(θ), cos(θ), 0.0, 0.0, 0.0, 1.0)
    Rx(θ) = SMatrix{3, 3}(1.0, 0.0, 0.0, 0.0, cos(θ), sin(θ), 0.0, -sin(θ), cos(θ))
    R = Rz(deg2rad(raan_deg)) * Rx(deg2rad(inc_deg)) * Rz(deg2rad(argp_deg))
    M = lock(RuntimeServices.SPICE_LOCK) do
        SMatrix{3, 3, Float64}(pxform(frame, "J2000", et))
    end
    # The elements are inertial (osculating about the lunar pole); the
    # body-fixed frame only supplies the axes at the epoch, so no rotation
    # term is added to the velocity.
    return SM.CartesianInitialCondition(M * (R * r_pf), M * (R * v_pf))
end

function main(argv=ARGS)
    model = demo_model_path("apollo_lunar_module_nasa_3d_resources.glb")
    planet = Moon("", SPICE_PATH)
initial_time = initial_time_of(et_of("1969-07-20T17:44:00"))
et0 = et_of(initial_time) # Use the same epoch for the frame transform and engine.
ic = lunar_orbit_ic(planet, et0; rp_alt_m=100.9e3, ra_alt_m=122.4e3, inc_deg=178.75, raan_deg=0.0, argp_deg=0.0, ta_deg=0.0)
a0 = 1.0 / (2.0 / norm(ic.pos) - dot(ic.vel, ic.vel) / planet.μ)
period = 2pi * sqrt(a0^3 / planet.μ)
options = viewer_demo_options("apollo11_lunar_orbit", 2.0 * period; argv=argv)
OUTDIR = options.output_dir
mission_time = options.duration_s
println("Apollo 11 LM epoch ", utc_of(et0), "  a=", round(a0 / 1e3; digits=1), " km  period=", round(period / 60; digits=1), " min")

sc = make_three_body_spacecraft(
    bus_dims=(4.2, 4.2, 7.0), panel_dims=(0.01, 0.02, 0.02), bus_mass=15_000.0, panel_mass_each=0.5, panel_offset_y=2.2,
    ic=ic, reflection_coefficient=0.6, prop_mass=100.0, id=1)
effectors = (
    GravitationalHarmonicsModel(50, 50, joinpath(HARMONICS_DIR, "LP165P.csv"), planet),
    NBodyGravityModel(body_names=("Earth", "Sun"), primary_body_name="Moon", planet=planet),
    SolarRadiationPressureModel(1.3, 15.0),
)
base = make_example_config(planet=planet, spacecraft=sc, mission_time=mission_time, initial_time=initial_time,
    dynamic_effectors=effectors, density_model=NoAtmosphereModel(), orientation_sim=false,
    keplerian=true, EI_km=1.0, verbose=false, results=true, results_directory=OUTDIR)
args = SM.SimConfig._with_configuration(base;
    mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, keplerian=true, number_of_orbits=1,
        mission_time=mission_time, orientation_sim=false, num_steps_to_save=1000, data_rate=5.0),
    solver_config=SM.SolverConfig(solver_mode=:tsit5),
    integration_tolerances=SM.IntegrationTolerances(reltol_orbit=1e-10, abstol_orbit=1e-9, dt_max_orbit=10.0,
        reltol_atmosphere=1e-10, abstol_atmosphere=1e-9, dt_max_atmosphere=1.0))

prefix = run_or_reuse!(args, OUTDIR)
summarize_run(prefix, planet; alt_m=110e3)

html = export_visualization(prefix; max_frames=4000, trail_orbits=1, texture_resolution="4k",
    title="AGORA Apollo 11 · lunar module in the parking orbit",
    models=Dict(1 => model), model_scale=1.45, model_rotation_deg=Dict(1 => (-180, 90, 90)))
println("html: ", html, " ", filesize(html))
cdn = build_cdn_page(html, joinpath(OUTDIR, "artifact.html"), "AGORA Apollo 11 Lunar Orbit",
    "AGORA Apollo 11 · lunar module in the parking orbit, 1969-07-20 17:44 UTC",
    "101 x 122 km, retrograde, 1.25° off the lunar equator; LP165P 50x50, Earth + Sun, SRP", "$(round(mission_time / 3600; digits=3)) h modeled span",
    "Nominal Apollo 11 parking-orbit illustration; lunar orientation uses SPICE, but no flown LM reference is supplied. Drag to orbit, wheel to zoom, Space to pause, F to follow the lander.")
cdn === nothing || println("cdn: ", cdn, " ", filesize(cdn))

return (args=args, prefix=prefix, html=html, cdn=cdn)

end # main

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
