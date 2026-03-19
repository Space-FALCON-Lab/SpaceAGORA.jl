using CSV
using DataFrames
using StaticArrays

if !isdefined(@__MODULE__, :REPO_ROOT)
    const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
end
if !isdefined(@__MODULE__, :SPICE_PATH)
    const SPICE_PATH = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
end
if !isdefined(@__MODULE__, :SimulationModel)
    include(joinpath(REPO_ROOT, "src", "simulation_model", "SimulationModel.jl"))
end
if !isdefined(@__MODULE__, :SM)
    const SM = SimulationModel
end

@inline _example_smoke_enabled() = get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1"
@inline _example_smoke_results_enabled() = get(ENV, "SPACEAGORA_EXAMPLE_SMOKE_RESULTS", "0") == "1"

@inline function _example_smoke_mission_time(default_time::Float64)
    raw = get(ENV, "SPACEAGORA_EXAMPLE_SMOKE_MISSION_TIME", "120.0")
    parsed = tryparse(Float64, raw)
    if parsed === nothing || !(parsed > 0.0)
        return min(default_time, 120.0)
    end
    return min(default_time, parsed)
end

function _example_smoke_args(args::SM.SimulationConfiguration)
    if !_example_smoke_enabled()
        return args
    end

    mc = args.mission_configuration
    ss = args.simulation_settings
    keep_results = _example_smoke_results_enabled()
    mc_smoke = SM.MissionConfiguration(
        mission_type=mc.mission_type,
        keplerian=mc.keplerian,
        number_of_orbits=max(1, min(mc.number_of_orbits, 1)),
        mission_time=_example_smoke_mission_time(mc.mission_time),
        orientation_sim=mc.orientation_sim,
        num_steps_to_save=max(50, min(mc.num_steps_to_save, 200))
    )
    ss_smoke = SM.SimulationSettings(
        results=keep_results,
        verbose=false,
        # Keep smoke outputs local to the current run directory to avoid cross-run collisions.
        results_directory=joinpath(pwd(), "output"),
        generate_plots=false,
        generate_filenames=ss.generate_filenames,
        normalize=false,
        save_csv=keep_results
    )

    return SM.SimulationConfiguration(
        file_paths=args.file_paths,
        simulation_settings=ss_smoke,
        mission_configuration=mc_smoke,
        environment_model=args.environment_model,
        dynamics_model=args.dynamics_model,
        guidance_model=args.guidance_model,
        navigation_model=args.navigation_model,
        control_model=args.control_model,
        initial_time=args.initial_time,
        integration_tolerances=args.integration_tolerances
    )
end

function make_three_body_spacecraft(;
    bus_dims::NTuple{3, Float64},
    panel_dims::NTuple{3, Float64},
    bus_mass::Float64,
    panel_mass_each::Float64,
    panel_offset_y::Float64,
    ic::SM.InitialCondition,
    prop_mass::Float64=0.0,
    id::Int64=1
)
    main_bus = SM.Link{0}(root=true, m=bus_mass, dims=MVector{3, Float64}(bus_dims...), ref_area=bus_dims[1] * bus_dims[3])
    left_panel = SM.Link{0}(root=false, m=panel_mass_each, dims=MVector{3, Float64}(panel_dims...), ref_area=panel_dims[2] * panel_dims[3], r=MVector{3, Float64}(0.0, -panel_offset_y, 0.0))
    right_panel = SM.Link{0}(root=false, m=panel_mass_each, dims=MVector{3, Float64}(panel_dims...), ref_area=panel_dims[2] * panel_dims[3], r=MVector{3, Float64}(0.0, panel_offset_y, 0.0))

    return SM.SpacecraftModel(
        SM.Joint[],
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
    planet::SM.AbstractPlanet,
    spacecraft::SM.SpacecraftModel,
    # mission_time::Float64=1.0e6,
    orbits::Int64=50,
    data_rate::Float64=10.0,
    initial_time::SM.InitialTime,
    dynamic_effectors::Tuple=(SM.InverseSquaredJ2GravityModel(),),
    density_model::SM.AbstractDensityModel=SM.NoAtmosphereModel(),
    orientation_sim::Bool=false,
    keplerian::Bool=true,
    EI_km::Float64=300.0,
    verbose::Bool=true
)
    return SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(
            results=true,
            verbose=verbose,
            generate_plots=false,
            results_directory=joinpath(REPO_ROOT, "output"),
            normalize=false
        ),
        mission_configuration=SM.MissionConfiguration(
            mission_type=SM.MissionOrbits,
            keplerian=keplerian,
            number_of_orbits=orbits,
            # mission_time=mission_time,
            orientation_sim=orientation_sim,
            num_steps_to_save=1000,
            data_rate=data_rate
        ),
        environment_model=SM.EnvironmentModel(
            planet=planet,
            EI=EI_km,
            density_model=density_model,
            thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false
        ),
        dynamics_model=SM.DynamicsModel([spacecraft], dynamic_effectors),
        guidance_model=SM.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=SM.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=SM.ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=initial_time,
        integration_tolerances=SM.IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            dt_max_orbit=20.0,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_atmosphere=0.2
        )
    )
end

function make_example_config(;
    planet::SM.AbstractPlanet,
    spacecraft::SM.SpacecraftModel,
    mission_time::Float64=1.0e6,
    data_rate::Float64=10.0,
    initial_time::SM.InitialTime,
    dynamic_effectors::Tuple=(SM.InverseSquaredJ2GravityModel(),),
    density_model::SM.AbstractDensityModel=SM.NoAtmosphereModel(),
    orientation_sim::Bool=false,
    keplerian::Bool=true,
    EI_km::Float64=300.0,
    verbose::Bool=true
)
    return SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(
            results=true,
            verbose=verbose,
            generate_plots=false,
            results_directory=joinpath(REPO_ROOT, "output"),
            normalize=false
        ),
        mission_configuration=SM.MissionConfiguration(
            mission_type=SM.MissionTime,
            keplerian=keplerian,
            # number_of_orbits=orbits,
            mission_time=mission_time,
            orientation_sim=orientation_sim,
            num_steps_to_save=1000,
            data_rate=data_rate
        ),
        environment_model=SM.EnvironmentModel(
            planet=planet,
            EI=EI_km,
            density_model=density_model,
            thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false
        ),
        dynamics_model=SM.DynamicsModel([spacecraft], dynamic_effectors),
        guidance_model=SM.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=SM.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=SM.ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=initial_time,
        integration_tolerances=SM.IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            dt_max_orbit=20.0,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_atmosphere=0.2
        )
    )
end

function run_and_report(args::SM.SimulationConfiguration; save_fields=nothing)
    args_eff = _example_smoke_args(args)
    t = @elapsed run_simulation(args_eff; save_fields=save_fields)
    csv_path = joinpath(args_eff.simulation_settings.results_directory, "simulation_results.csv")
    if args_eff.simulation_settings.results && isfile(csv_path)
        df = CSV.read(csv_path, DataFrame)
        println("Saved $(nrow(df)) samples to $(abspath(csv_path))")
    end
    println("COMPUTATIONAL TIME = $(t) s")
end
