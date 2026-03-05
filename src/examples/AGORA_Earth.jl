if !isdefined(@__MODULE__, :SimulationModel)
    include("../simulation_model/SimulationModel.jl")
end
using .SimulationModel
using CSV
using DataFrames
using SPICE
using StaticArrays

# run_simulation.jl expects quat_mult in including scope.
const quat_mult = SimulationModel.quat_mult
if !isdefined(@__MODULE__, :run_simulation)
    include("../simulation/execution/run_simulation.jl")
end
if !isdefined(@__MODULE__, :make_example_config)
    include("typed_example_utils.jl")
end

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SPICE_PATH = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")

planet = Earth("", SPICE_PATH)

main_bus = Link{0}(root=true, m=620.0, ref_area=2.05 * 2.8)
left_panel = Link{0}(root=false, m=10.0, ref_area=5.7 * 1.0 / 2.0, r=MVector{3, Float64}(0.0, -2.05 / 2.0 - 5.7 / 4.0, 0.0))
right_panel = Link{0}(root=false, m=10.0, ref_area=5.7 * 1.0 / 2.0, r=MVector{3, Float64}(0.0, 2.05 / 2.0 + 5.7 / 4.0, 0.0))

ic = InitialCondition(
    ra=56_378.7978559e3,
    rp=planet.Rp_e + 200_590.0,
    i=89.876,
    ω=75.505,
    Ω=104.115,
    ν=175.0
)

spacecraft = SpacecraftModel(
    Joint[],
    [main_bus, left_panel, right_panel],
    main_bus,
    true,
    main_bus.m + left_panel.m + right_panel.m,
    200.0,
    main_bus.inertia,
    0,
    0,
    ic,
    1
)

args = SimulationConfiguration(
    simulation_settings=SimulationSettings(
        results=true,
        verbose=true,
        generate_plots=false,
        results_directory=joinpath(REPO_ROOT, "output"),
        normalize=false
    ),
    mission_configuration=MissionConfiguration(
        mission_type=MissionTime,
        keplerian=true,
        number_of_orbits=1,
        mission_time=12.0 * 3600.0,
        orientation_sim=false,
        num_steps_to_save=1000
    ),
    environment_model=EnvironmentModel(
        planet=planet,
        EI=300.0,
        density_model=NoAtmosphereModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false,
        wind=false
    ),
    dynamics_model=DynamicsModel([spacecraft], (InverseSquaredJ2GravityModel(),)),
    guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
    navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
    control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
    initial_time=InitialTime(
        year=2014,
        month=5,
        day=27,
        hour=5,
        minute=0,
        second=0.0
    ),
    integration_tolerances=IntegrationTolerances(
        reltol_orbit=1e-8,
        abstol_orbit=1e-8,
        dt_max_orbit=15.0
    )
)

run_and_report(args)
