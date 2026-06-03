include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "aerobraking_mission_plot_utils.jl"))

using CSV
using DataFrames
using SPICE
using StaticArrays

setup_gram_example!()


struct ConstantDensityModel <: AbstractDensityModel
    rho::Float64
    temp::Float64
end

function SimulationModel.EnvironmentModels.getDensity(
    model::ConstantDensityModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p
)
    return model.rho, model.temp, SVector{3, Float64}(0.0, 0.0, 0.0)
end

function _load_aerobraking_maneuver_schedule(csv_path::String; burn_orbit_offset::Int=0)
    isfile(csv_path) || throw(ArgumentError("Maneuver schedule CSV not found at $(abspath(csv_path))."))
    df = CSV.read(csv_path, DataFrame)
    required = (:passage_number, :delta_v)
    missing_cols = setdiff(required, Symbol.(names(df)))
    isempty(missing_cols) || throw(ArgumentError("Maneuver schedule CSV missing columns: $(missing_cols)."))

    passage_number = Int.(df.passage_number)
    delta_v = Float64.(df.delta_v)
    keep = findall(!iszero, delta_v)
    return (
        maneuver_orbit_number=passage_number[keep] .+ burn_orbit_offset,
        maneuver_Δv=delta_v[keep]
    )
end

_set_aerobraking_mission_spice_config!(VEX_SPICE_CONFIG)

planet = Venus("", SPICE_PATH)
smoke_mode = get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1"
vex_schedule = _load_aerobraking_maneuver_schedule(joinpath(REPO_ROOT, "data", "Maneuver_plans", "vex.csv"))
requested_initial_time = InitialTime(year=2014, month=5, day=19, hour=14, minute=7, second=32.0)
vex_initial = _nearest_apoapsis_initial_time_and_condition_from_spice(
    requested_initial_time,
    SPICE_PATH,
    VEX_SPICE_CONFIG;
    search_window_s=10.0 * 24.0 * 3600.0
)
initial_time = vex_initial.initial_time
ic = vex_initial.ic

spacecraft = make_three_body_spacecraft(
    bus_dims=(2.05, 3.7, 2.8),
    panel_dims=(0.01, 5.7 / 1.0, 1.0),
    bus_mass=620.0,
    panel_mass_each=10.0,
    panel_offset_y=2.05 / 2.0 + 5.7 / 4.0,
    ic=ic,
    prop_mass=10.0,
    id=102
)

sun_gravity = NBodyGravityModel(body_names=("Sun",), primary_body_name="Venus", planet=planet)
solar_radiation_pressure = SolarRadiationPressureModel(spacecraft.root.reflection_coefficient, spacecraft.root.ref_area)
venus_harmonics_file = joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "MGNP180U.csv")
dynamic_effectors = smoke_mode ? (
    InverseSquaredGravityModel(),
    sun_gravity,
    solar_radiation_pressure
) : (
    GravitationalHarmonicsModel(50, 50, venus_harmonics_file, planet),
    sun_gravity,
    solar_radiation_pressure,
    AerodynamicCoefficientfM()
)

base_args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=smoke_mode ? 100_000.0 : 50.0 * 24.0 * 3_600.0,
    initial_time=initial_time,
    dynamic_effectors=dynamic_effectors,
    density_model=GRAMAtmosphereModel(planet_name="venus"),
    orientation_sim=false,
    keplerian=smoke_mode,
    EI_km=250.0
)

thruster = BaseThrusterModel(
    thrust=[4.0],
    direction=[0.0],
    Δv=[0.0],
    start_burn_time=[-1.0],
    stop_burn_time=[-1.0],
    Isp=[300.0]
)
guidance_effector = AerobrakingCampaignPropulsiveManeuverGuidanceModel(
    maneuver_orbit_number=vex_schedule.maneuver_orbit_number,
    maneuver_Δv=vex_schedule.maneuver_Δv
)
args = SimulationConfiguration(
    file_paths=base_args.file_paths,
    simulation_settings=base_args.simulation_settings,
    mission_configuration=base_args.mission_configuration,
    environment_model=base_args.environment_model,
    dynamics_model=base_args.dynamics_model,
    guidance_model=GuidanceModel(guidance_effectors=(guidance_effector,), guidance_rates=[30.0]),
    navigation_model=base_args.navigation_model,
    control_model=ControlModel(control_effectors=(thruster,), control_rates=[10.0]),
    initial_time=base_args.initial_time,
    integration_tolerances=IntegrationTolerances(
        reltol_orbit=1e-8,
        abstol_orbit=1e-8,
        dt_max_orbit=30.0,
        reltol_atmosphere=1e-8,
        abstol_atmosphere=1e-8,
        dt_max_atmosphere=1.0
    )
)

args_eff = SpaceAGORA.TelemetryVerification._example_smoke_args(args)
sim_elapsed = @elapsed sol = run_simulation(args_eff; return_solution=true)
csv_path = joinpath(args_eff.simulation_settings.results_directory, "simulation_results.csv")
if args_eff.simulation_settings.results && isfile(csv_path)
    df = CSV.read(csv_path, DataFrame)
    println("Saved $(nrow(df)) samples to $(abspath(csv_path))")
end
println("COMPUTATIONAL TIME = $(sim_elapsed) s")

_save_periapsis_event_table(
    sol,
    args_eff,
    args_eff.environment_model.planet;
    allow_empty=smoke_mode
)
_save_spice_periapsis_event_table(
    args_eff,
    args_eff.environment_model.planet,
    args_eff.initial_time,
    SPICE_PATH;
    allow_empty=smoke_mode
)
_save_spice_apoapsis_event_table(
    args_eff,
    args_eff.environment_model.planet,
    args_eff.initial_time,
    SPICE_PATH;
    allow_empty=smoke_mode
)
if _has_orbit_extrema_for_plots(args_eff, args_eff.environment_model.planet)
    _save_apoapsis_periapsis_plot(args_eff, args_eff.environment_model.planet, vex_schedule)
    _save_apoapsis_radius_plot(args_eff, args_eff.environment_model.planet, vex_schedule)
    _save_periapsis_latlon_plot(args_eff, args_eff.environment_model.planet, vex_schedule)
    _save_drag_along_velocity_plot(args_eff)
    _save_aero_sideways_components_plot(args_eff)
end
trajectory_plot_path = _save_trajectory_comparison_plot(args_eff, args_eff.environment_model.planet, args_eff.initial_time, SPICE_PATH)
planet_fixed_trajectory_plot_path = _save_planet_fixed_trajectory_comparison_plot(args_eff, args_eff.environment_model.planet, args_eff.initial_time, SPICE_PATH)
_save_planet_fixed_axis_trajectory_plot(args_eff, args_eff.environment_model.planet, args_eff.initial_time, SPICE_PATH)
_save_orbital_elements_comparison_plot(args_eff, args_eff.environment_model.planet, args_eff.initial_time, SPICE_PATH)
_save_inertial_position_error_plot(args_eff, args_eff.initial_time, SPICE_PATH)
_save_inertial_velocity_error_plot(args_eff, args_eff.initial_time, SPICE_PATH)
_save_planet_fixed_position_error_plot(args_eff, args_eff.environment_model.planet, args_eff.initial_time, SPICE_PATH)
_save_planet_fixed_velocity_error_plot(args_eff, args_eff.environment_model.planet, args_eff.initial_time, SPICE_PATH)
_save_rtn_position_error_plot(args_eff, args_eff.initial_time, SPICE_PATH)
_save_rtn_velocity_error_plot(args_eff, args_eff.initial_time, SPICE_PATH)
_open_plot_in_browser(trajectory_plot_path)
