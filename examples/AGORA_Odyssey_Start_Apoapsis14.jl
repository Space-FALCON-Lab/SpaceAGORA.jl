include(joinpath(@__DIR__, "common.jl"))
include(joinpath(REPO_ROOT, "src", "mission", "operations", "maneuver_plans.jl"))
using CSV
using DataFrames
using Plots
using PlotlyJS: Plot, scatter3d, surface, Layout, attr
using Roots
using SPICE
using StaticArrays
using LinearAlgebra
using Logging

include(joinpath(@__DIR__, "aerobraking_mission_plot_utils.jl"))

const ORIGINAL_MISSION_TIME_S = 3_600.0 * 400.0
const START_APOAPSIS_ORBIT = 14
const START_APOAPSIS_TIME_S = 735459.713412683

planet = Mars("", SPICE_PATH)
smoke_mode = get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1"
initial_time = InitialTime(year=2001, month=11, day=14, hour=22, minute=22, second=52.413270)

ic = _mars_odyssey_initial_condition_from_spice(initial_time, SPICE_PATH)

spacecraft = make_three_body_spacecraft(
    bus_dims=(2.2, 2.6, 1.7),
    panel_dims=(0.01, 5.5 / 1.35, 2.6),
    bus_mass=391.0,
    panel_mass_each=10.0,
    panel_offset_y=2.6 / 2.0 + 5.5 / 4.0,
    ic=ic,
    reflection_coefficient=0.9,
    prop_mass=50.0,
    id=100
)

mars_harmonics_file = joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "GMM3.csv")
sun_gravity = NBodyGravityModel(body_names=("Sun",), primary_body_name="Mars", planet=planet)
solar_radiation_pressure = SolarRadiationPressureModel(spacecraft.root.reflection_coefficient, spacecraft.root.ref_area)

dynamic_effectors = smoke_mode ? (
    InverseSquaredGravityModel(),
    sun_gravity
) : (
    InverseSquaredGravityModel(),
    sun_gravity,
    GravitationalHarmonicsModel(50, 50, mars_harmonics_file, planet),
    solar_radiation_pressure,
    AerodynamicCoefficientfM()
)
density_model = GRAMAtmosphereModel(planet_name="mars")
base_args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=ORIGINAL_MISSION_TIME_S,
    initial_time=initial_time,
    dynamic_effectors=dynamic_effectors,
    density_model=density_model,
    orientation_sim=false,
    keplerian=smoke_mode,
    EI_km=160.0,
    results_directory=joinpath(REPO_ROOT, "output", "odyssey_start_apoapsis_14")
)

odyssey_schedule_base = odyssey_campaign_maneuvers(1:311; planet=planet)
original_guidance_orbits = odyssey_schedule_base.maneuver_orbit_number .+ 1
# Convert the full-campaign guidance apoapsis number to this restart run's
# local orbit counter. The restart begins on START_APOAPSIS_ORBIT itself.
local_guidance_orbits = original_guidance_orbits .- (START_APOAPSIS_ORBIT - 1)
future_maneuvers = (original_guidance_orbits .>= START_APOAPSIS_ORBIT) .& (local_guidance_orbits .>= 1)
odyssey_schedule = (
    maneuver_orbit_number=local_guidance_orbits[future_maneuvers],
    maneuver_Δv=odyssey_schedule_base.maneuver_Δv[future_maneuvers]
)
odyssey_plot_schedule = (
    maneuver_orbit_number=odyssey_schedule.maneuver_orbit_number,
    maneuver_Δv=odyssey_schedule.maneuver_Δv
)
println("Starting from original apoapsis orbit $(START_APOAPSIS_ORBIT); shifted maneuver orbits: $(odyssey_schedule.maneuver_orbit_number)")

thruster = BaseThrusterModel(
    thrust=[4.0],
    direction=[0.0],
    Δv=[0.0],
    start_burn_time=[-1.0],
    stop_burn_time=[-1.0],
    Isp=[300.0]
)
guidance_effector = AerobrakingCampaignPropulsiveManeuverGuidanceModel(
    maneuver_orbit_number=odyssey_schedule.maneuver_orbit_number,
    maneuver_Δv=odyssey_schedule.maneuver_Δv
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
periapsis_event_path = _save_periapsis_event_table(
    sol,
    args_eff,
    args_eff.environment_model.planet;
    allow_empty=smoke_mode
)
spice_periapsis_event_path = _save_spice_periapsis_event_table(
    args_eff,
    args_eff.environment_model.planet,
    args_eff.initial_time,
    SPICE_PATH;
    allow_empty=smoke_mode
)
spice_apoapsis_event_path = _save_spice_apoapsis_event_table(
    args_eff,
    args_eff.environment_model.planet,
    args_eff.initial_time,
    SPICE_PATH;
    allow_empty=smoke_mode
)
if _has_orbit_extrema_for_plots(args_eff, args_eff.environment_model.planet)
    _save_apoapsis_periapsis_plot(args_eff, args_eff.environment_model.planet, odyssey_plot_schedule)
    _save_apoapsis_radius_plot(args_eff, args_eff.environment_model.planet, odyssey_plot_schedule)
    _save_periapsis_latlon_plot(args_eff, args_eff.environment_model.planet, odyssey_plot_schedule)
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
