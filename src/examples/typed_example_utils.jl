using CSV
using DataFrames
using StaticArrays

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SPICE_PATH = joinpath(REPO_ROOT, "GRAM_Data", "SPICE")

function make_three_body_spacecraft(;
    bus_dims::NTuple{3, Float64},
    panel_dims::NTuple{3, Float64},
    bus_mass::Float64,
    panel_mass_each::Float64,
    panel_offset_y::Float64,
    ic::InitialCondition,
    prop_mass::Float64=0.0,
    id::Int64=1
)
    main_bus = Link{0}(root=true, m=bus_mass, dims=MVector{3, Float64}(bus_dims...), ref_area=bus_dims[1] * bus_dims[3])
    left_panel = Link{0}(root=false, m=panel_mass_each, dims=MVector{3, Float64}(panel_dims...), ref_area=panel_dims[2] * panel_dims[3], r=MVector{3, Float64}(0.0, -panel_offset_y, 0.0))
    right_panel = Link{0}(root=false, m=panel_mass_each, dims=MVector{3, Float64}(panel_dims...), ref_area=panel_dims[2] * panel_dims[3], r=MVector{3, Float64}(0.0, panel_offset_y, 0.0))

    return SpacecraftModel(
        Joint[],
        [main_bus, left_panel, right_panel],
        main_bus,
        true,
        main_bus.m + left_panel.m + right_panel.m,
        prop_mass,
        main_bus.inertia,
        0,
        0,
        ic,
        id
    )
end

function make_example_config(;
    planet::AbstractPlanet,
    spacecraft::SpacecraftModel,
    mission_time::Float64,
    initial_time::InitialTime,
    dynamic_effectors::Tuple=(InverseSquaredJ2GravityModel(),),
    density_model::AbstractDensityModel=NoAtmosphereModel(),
    orientation_sim::Bool=false,
    keplerian::Bool=true,
    EI_km::Float64=300.0,
    verbose::Bool=true
)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=true,
            verbose=verbose,
            generate_plots=false,
            results_directory=joinpath(REPO_ROOT, "output")
        ),
        mission_configuration=MissionConfiguration(
            mission_type="Time",
            keplerian=keplerian,
            number_of_orbits=1,
            mission_time=mission_time,
            orientation_sim=orientation_sim,
            num_steps_to_save=1000
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=EI_km,
            density_model=density_model,
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false
        ),
        dynamics_model=DynamicsModel([spacecraft], dynamic_effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=initial_time,
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            dt_max_orbit=20.0,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_atmosphere=0.2
        )
    )
end

function run_and_report(args::SimulationConfiguration)
    t = @elapsed run_simulation(args)
    if isfile("simulation_results.csv")
        df = CSV.read("simulation_results.csv", DataFrame)
        println("Saved $(nrow(df)) samples to $(abspath("simulation_results.csv"))")
    end
    println("COMPUTATIONAL TIME = $(t) s")
end
