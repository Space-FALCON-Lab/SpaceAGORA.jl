const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

using SPICE
using SpaceAGORA
using SpaceAGORA.SimulationModel

const SimulationModel = SpaceAGORA.SimulationModel
const quat_mult = SimulationModel.quat_mult
const run_simulation = SpaceAGORA.run_simulation

spice_path = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
planet = Earth("", spice_path)

root = Link{0}(root=true, m=120.0, ref_area=1.0)
ic = InitialCondition(planet.Rp_e + 500e3, 0.0, 0.0, 0.0, 0.0, 0.0)
spacecraft = SpacecraftModel(
    Joint[],
    [root],
    root,
    true,
    root.m,
    0.0,
    root.inertia,
    0,
    0,
    ic,
    1
)

args = SimulationConfiguration(
    simulation_settings=SimulationSettings(
        results=false,
        verbose=false,
        generate_plots=false,
        normalize=false
    ),
    mission_configuration=MissionConfiguration(
        mission_type=MissionTime,
        keplerian=true,
        number_of_orbits=1,
        mission_time=5.0,
        orientation_sim=false,
        num_steps_to_save=10
    ),
    environment_model=EnvironmentModel(
        planet=planet,
        EI=120.0,
        density_model=NoAtmosphereModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false,
        wind=false
    ),
    dynamics_model=DynamicsModel([spacecraft], (InverseSquaredGravityModel(),)),
    guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
    navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
    control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
    initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
    integration_tolerances=IntegrationTolerances(reltol_orbit=1e-8, abstol_orbit=1e-8, dt_max_orbit=5.0)
)

mktempdir() do tmp
    cd(tmp) do
        run_simulation(args)
        if isfile("simulation_results.csv")
            error("Smoke run should not write simulation_results.csv when results=false")
        end
    end
end

println("clean_depot_smoke_ok")
