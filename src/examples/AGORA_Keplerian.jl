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
    include("../simulation/run_simulation.jl")
end
if !isdefined(@__MODULE__, :make_example_config)
    include("typed_example_utils.jl")
end

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SPICE_PATH = joinpath(REPO_ROOT, "GRAM_Data", "SPICE")

planet = Mars("", SPICE_PATH)

main_bus = Link{0}(root=true, m=411.0, ref_area=2.2 * 1.7)
left_panel = Link{0}(root=false, m=25.0, ref_area=3.76 * 1.93 / 2.0, r=MVector{3, Float64}(0.0, -2.6 / 2.0 - 3.76 / 4.0, 0.0))
right_panel = Link{0}(root=false, m=25.0, ref_area=3.76 * 1.93 / 2.0, r=MVector{3, Float64}(0.0, 2.6 / 2.0 + 3.76 / 4.0, 0.0))

ic = InitialCondition(
    ra=28_038e3,
    rp=24_600e3,
    i=93.6,
    ω=0.0,
    Ω=0.0,
    ν=175.0
)

spacecraft = SpacecraftModel(
    Joint[],
    [main_bus, left_panel, right_panel],
    main_bus,
    true,
    main_bus.m + left_panel.m + right_panel.m,
    50.0,
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
        mission_time=8.0 * 3600.0,
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
        year=1993,
        month=5,
        day=25,
        hour=14,
        minute=21,
        second=28.0
    ),
    integration_tolerances=IntegrationTolerances(
        reltol_orbit=1e-8,
        abstol_orbit=1e-8,
        dt_max_orbit=20.0
    )
)

run_and_report(args)
