# ISS in its orbit with NASA's ISS (B) glTF standing in for the link boxes.
using SpaceAGORA, StaticArrays
const SM = SpaceAGORA.SimulationModel
planet = make_no_gram_planet(:earth)
bus = SM.Link(root=true, m=420_000.0, dims=MVector{3, Float64}(73.0, 109.0, 20.0), ref_area=1500.0)
sc = SM.SpacecraftModel(root=bus, initial_condition=SM.InitialCondition(ra=planet.Rp_e + 420e3, rp=planet.Rp_e + 408e3, i=51.64, ω=30.0, Ω=120.0, ν=0.0), id=1)
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
outdir = joinpath(get(ENV, "SPACEAGORA_VIEWER_DEMO_OUT", joinpath(REPO_ROOT, "output", "viewer_demos")), "iss_orbit"); rm(outdir; force=true, recursive=true)
import SpaceAGORA.TelemetryVerification: make_example_config
args = make_example_config(planet=planet, spacecraft=sc, mission_time=2.0 * 5560.0,
    initial_time=SM.InitialTime(year=2024, month=6, day=21, hour=10, minute=0, second=0.0),
    dynamic_effectors=(SM.InverseSquaredJ2GravityModel(),), density_model=SM.NoAtmosphereModel(), ephemerides_model=SM.SimpleEphemeridesModel(),
    orientation_sim=false, keplerian=true, EI_km=120.0, verbose=false, results=true, results_directory=outdir)
@time run_simulation(args; visualization=true)
iss = joinpath(REPO_ROOT, "data", "models", "iss_nasa_3d_resources_b.glb")
html = export_visualization(joinpath(outdir, "simulation_results"); max_frames=2500, trail_orbits=1, title="AGORA ISS · NASA 3D model in a 408 km orbit",
    models=Dict(1 => iss), model_scale=2.4, model_rotation_deg=Dict(1 => (-90, 0, -90)))
println("html: ", html, " ", filesize(html))
