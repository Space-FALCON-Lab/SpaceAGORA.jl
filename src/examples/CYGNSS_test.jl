if !isdefined(@__MODULE__, :SimulationModel)
    include("../simulation_model/SimulationModel.jl")
end
using .SimulationModel
using SPICE
using StaticArrays

# run_simulation.jl expects quat_mult in the including scope.
const quat_mult = SimulationModel.quat_mult
if !isdefined(@__MODULE__, :run_simulation)
    include("../simulation/run_simulation.jl")
end
if !isdefined(@__MODULE__, :make_example_config)
    include("typed_example_utils.jl")
end

planet = Earth("", SPICE_PATH)

function make_cygnss_spacecraft(id::Int64, raan_deg::Float64)
    ic = InitialCondition(
        ra=planet.Rp_e + 540e3,
        rp=planet.Rp_e + 500e3,
        i=35.0,
        ω=0.0,
        Ω=raan_deg,
        ν=170.0
    )

    return make_three_body_spacecraft(
        bus_dims=(0.6, 0.5, 0.5),
        panel_dims=(0.01, 0.4, 0.3),
        bus_mass=28.0,
        panel_mass_each=1.0,
        panel_offset_y=0.35,
        ic=ic,
        prop_mass=0.0,
        id=id
    )
end

spacecraft = [
    make_cygnss_spacecraft(201, 0.0),
    make_cygnss_spacecraft(202, 20.0),
    make_cygnss_spacecraft(203, 40.0),
    make_cygnss_spacecraft(204, 60.0)
]

args = SimulationConfiguration(
    simulation_settings=SimulationSettings(
        results=true,
        verbose=true,
        generate_plots=false,
        results_directory=joinpath(REPO_ROOT, "output")
    ),
    mission_configuration=MissionConfiguration(
        mission_type="Time",
        keplerian=true,
        number_of_orbits=1,
        mission_time=1_800.0,
        orientation_sim=false,
        num_steps_to_save=1000
    ),
    environment_model=EnvironmentModel(
        planet=planet,
        EI=120.0,
        density_model=NoAtmosphereModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false,
        wind=false
    ),
    dynamics_model=DynamicsModel(spacecraft, (InverseSquaredJ2GravityModel(),)),
    guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
    navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
    control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
    initial_time=InitialTime(year=2018, month=12, day=15, hour=0, minute=0, second=0.0),
    integration_tolerances=IntegrationTolerances(
        reltol_orbit=1e-8,
        abstol_orbit=1e-8,
        dt_max_orbit=10.0
    )
)

run_and_report(args)
