# AGORA_Earth.jl configuration, propagated for four days so the viewer holds several orbits.
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
using SpaceAGORA
using SpaceAGORA.SimulationModel
using StaticArrays
planet = make_no_gram_planet(:earth)
main_bus = Link(root=true, m=620.0, ref_area=2.05 * 2.8, dims=MVector{3, Float64}(2.05, 2.05, 2.8))
left_panel = Link(root=false, m=10.0, ref_area=5.7 * 1.0 / 2.0, dims=MVector{3, Float64}(0.02, 5.7 / 2.0, 1.0), r=MVector{3, Float64}(0.0, -2.05 / 2.0 - 5.7 / 4.0, 0.0))
right_panel = Link(root=false, m=10.0, ref_area=5.7 * 1.0 / 2.0, dims=MVector{3, Float64}(0.02, 5.7 / 2.0, 1.0), r=MVector{3, Float64}(0.0, 2.05 / 2.0 + 5.7 / 4.0, 0.0))
ic = InitialCondition(ra=56_378.7978559e3, rp=planet.Rp_e + 200_590.0, i=89.876, ω=75.505, Ω=104.115, ν=175.0)
spacecraft = SpacecraftModel(Joint[], [main_bus, left_panel, right_panel], main_bus, true,
    main_bus.m + left_panel.m + right_panel.m, 200.0, main_bus.inertia, 0, 0, ic, 1)
include(joinpath(@__DIR__, "demo_options.jl"))
options = viewer_demo_options("earth_4day", 4.0 * 86400.0)
outdir = options.output_dir
args = SimulationConfiguration(
    simulation_settings=SimulationSettings(results=true, verbose=false, generate_plots=false, results_directory=outdir, normalize=false),
    mission_configuration=MissionConfiguration(mission_type=MissionTime, keplerian=true, number_of_orbits=1,
        mission_time=options.duration_s, orientation_sim=false, num_steps_to_save=1000),
    environment_model=EnvironmentModel(planet=planet, EI=300.0, density_model=NoAtmosphereModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet), topography=false, wind=false, ephemerides_model=SimpleEphemeridesModel()),
    dynamics_model=DynamicsModel([spacecraft], (InverseSquaredJ2GravityModel(),)),
    guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
    navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
    control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
    initial_time=InitialTime(year=2014, month=5, day=27, hour=5, minute=0, second=0.0),
    integration_tolerances=IntegrationTolerances(reltol_orbit=1e-8, abstol_orbit=1e-8, dt_max_orbit=15.0)
)
@time run_simulation(args; visualization=true)
html = export_visualization(joinpath(outdir, "simulation_results"); max_frames=4000, trail_orbits=3, title="AGORA Earth · $(options.duration_s / 86400) days")
println("html: ", html, " ", filesize(html))
