#=
"""
    Aerobraking MPC example.

    All case-specific values are set here: mode, constraints, limits, weights,
    solver tolerances, spacecraft geometry, and reference settings.
"""
=#
include(joinpath(pwd(), "examples", "common.jl"))

using CSV
using DataFrames
using StaticArrays

planet = Mars()
smoke_mode = get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1"
mpc_mode_name = lowercase(get(ENV, "SPACEAGORA_MPC_MODE", "med"))
mpc_mode = mpc_mode_name in ("et", "target", "targeting", "target_energy") ?
    TargetEnergyMode() :
    MaxEnergyDepletionMode()

limit_heat_rate_w_cm2 = 0.20
limit_heat_load_j_cm2 = 30.0
limit_drag_n = 1.0
limit_area_slew_m2_s = 0.20
drag_coefficient = 2.2
area_weight = mpc_mode isa MaxEnergyDepletionMode ? 1.0e-5 : 0.0
area_slew_weight = 0.0
slack_weight = mpc_mode isa MaxEnergyDepletionMode ? 1.0e-3 : 1.0e5
target_energy_weight = 5.0e-10
max_depletion_energy_weight = 1.0
osqp_eps_abs = mpc_mode isa MaxEnergyDepletionMode ? 1.0e-5 : 1.0e-7
osqp_eps_rel = mpc_mode isa MaxEnergyDepletionMode ? 1.0e-5 : 1.0e-7
osqp_max_iter = 10_000
spacecraft_index = 1
controlled_panel_links = (2, 3)
prediction_latitude_rad = 0.0
prediction_longitude_rad = 0.0
prediction_wind = false
min_panel_alpha_rad = 1.0e-4
max_panel_alpha_rad = pi / 2
mpc_control_dt_s = 0.5
mpc_solve_interval_s = smoke_mode ? 60.0 : 10.0
build_reference_on_tick = true
qp_max_nodes = smoke_mode ? 40 : 120
reference_cutoff_altitude_m = 300.0e3
reference_delta_s = 1.7e-7
reference_max_coast_steps = 2_000_000
reference_max_pass_steps = 20_000
ra = 28_559.615e3
rp = planet.Rp_e + 110e3
# Mars-specific target: the energy of an orbit with the configured periapsis
# and a modest, single-pass-reachable apoapsis reduction of 50 km.
target_apoapsis_radius_m = ra - 50.0e3
target_energy_mj_kg = -planet.μ / (target_apoapsis_radius_m + rp) / 1.0e6
orbit_period_s = 2pi * sqrt(((ra + rp) / 2)^3 / planet.μ)

ic = InitialCondition(
    ra=ra,
    rp=rp,
    i=93.6,
    ω=109.7454,
    Ω=28.1517,
    ν=smoke_mode ? 345.0 : 180.0,
)

spacecraft = make_three_body_spacecraft(
    bus_dims=(2.2, 2.6, 1.7),
    panel_dims=(0.01, 3.76 / 2.0, 1.93),
    bus_mass=391.0,
    panel_mass_each=10.0,
    panel_offset_y=2.6 / 2.0 + 3.76 / 4.0,
    ic=ic,
    reflection_coefficient=0.9,
    prop_mass=50.0,
    id=101,
)

density_model = SpaceAGORA.SimulationModel.EnvironmentModels.PolynomialFitAtmosphereModel(planet)

base_args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=smoke_mode ? 900.0 : orbit_period_s,
    initial_time=InitialTime(year=2001, month=11, day=6, hour=19, minute=0, second=0.0),
    dynamic_effectors=(InverseSquaredJ2GravityModel(), AerodynamicCoefficientfM()),
    density_model=density_model,
    ephemerides_model=SimpleEphemeridesModel(),
    orientation_sim=false,
    keplerian=false,
    EI_km=160.0,
    results_directory=joinpath(REPO_ROOT, "output", "mpc_aerobraking_external"),
)

constraints = mpc_constraints(:heat_rate, :drag, :slew)
mpc_config = mpc_config_from_spaceagora(
    base_args;
    mode=mpc_mode,
    spacecraft_index=spacecraft_index,
    controlled_panel_links=controlled_panel_links,
    constraints=constraints,
    qdot_max_w_cm2=limit_heat_rate_w_cm2,
    heat_load_max_j_cm2=limit_heat_load_j_cm2,
    drag_max_n=limit_drag_n,
    area_slew_max_m2_s=limit_area_slew_m2_s,
    drag_coefficient=drag_coefficient,
    area_weight=area_weight,
    area_slew_weight=area_slew_weight,
    slack_weight=slack_weight,
    target_energy_mj_kg=target_energy_mj_kg,
    target_energy_weight=target_energy_weight,
    max_depletion_energy_weight=max_depletion_energy_weight,
    osqp_eps_abs=osqp_eps_abs,
    osqp_eps_rel=osqp_eps_rel,
    osqp_max_iter=osqp_max_iter,
)
mpc_reference = AerobrakingMPCReferenceConfig(
    h_cut_m=reference_cutoff_altitude_m,
    delta_s=reference_delta_s,
    max_coast_steps=reference_max_coast_steps,
    max_pass_steps=reference_max_pass_steps,
)
mpc_state = AerobrakingMPCState()
mpc_control = AerobrakingMPCControlModel(
    config=mpc_config,
    state=mpc_state,
    spacecraft_index=spacecraft_index,
    controlled_panel_links=controlled_panel_links,
    control_dt_s=mpc_control_dt_s,
    min_alpha_rad=min_panel_alpha_rad,
    max_alpha_rad=max_panel_alpha_rad,
    solve_interval_s=mpc_solve_interval_s,
    build_reference_on_tick=build_reference_on_tick,
    qp_max_nodes=qp_max_nodes,
    reference=mpc_reference,
    fallback_area_m2=mpc_config.bus_reference_area_m2 + mpc_config.controllable_area_m2,
    prediction_latitude_rad=prediction_latitude_rad,
    prediction_longitude_rad=prediction_longitude_rad,
    prediction_wind=prediction_wind,
)

args = SimulationConfiguration(
    file_paths=base_args.file_paths,
    simulation_settings=base_args.simulation_settings,
    mission_configuration=MissionConfiguration(
        mission_type=base_args.mission_configuration.mission_type,
        keplerian=base_args.mission_configuration.keplerian,
        number_of_orbits=base_args.mission_configuration.number_of_orbits,
        mission_time=base_args.mission_configuration.mission_time,
        orientation_sim=base_args.mission_configuration.orientation_sim,
        num_steps_to_save=base_args.mission_configuration.num_steps_to_save,
        data_rate=1.0,
    ),
    environment_model=base_args.environment_model,
    dynamics_model=base_args.dynamics_model,
    guidance_model=base_args.guidance_model,
    navigation_model=base_args.navigation_model,
    control_model=ControlModel(control_effectors=(mpc_control,), control_rates=[mpc_control.control_dt_s]),
    initial_time=base_args.initial_time,
    integration_tolerances=base_args.integration_tolerances,
)

println("active_constraints = ", constraint_names(constraints))
println("mpc_mode = ", objective_label(mpc_config.mode))
println("mpc_config = ", (
    mode=objective_kind(mpc_config.mode),
    mass_kg=mpc_config.mass_kg,
    bus_area_m2=mpc_config.bus_reference_area_m2,
    controllable_area_m2=mpc_config.controllable_area_m2,
))
println("full_area_alpha_rad = ", alpha_from_commanded_area(
    mpc_config,
    mpc_config.bus_reference_area_m2 + mpc_config.controllable_area_m2,
    min_alpha_rad=min_panel_alpha_rad,
    max_alpha_rad=max_panel_alpha_rad,
))

save_fields = vcat(default_save_fields(args), collect(mpc_control_save_fields(mpc_control)))
t = @elapsed run_simulation(args; save_fields=save_fields, isolate_state=false)
csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
if args.simulation_settings.results && isfile(csv_path)
    df = CSV.read(csv_path, DataFrame)
    println("Saved $(nrow(df)) samples to $(abspath(csv_path))")
    area_history = df.sc1_mpc_commanded_area_m2
    drag_history = filter(isfinite, df.sc1_mpc_drag_n)
    heat_rate_history = filter(isfinite, df.sc1_mpc_heat_rate_w_cm2)
    area_slew_history = abs.(diff(area_history)) ./ diff(df.time)
    @assert maximum(drag_history) <= limit_drag_n + 1.0e-6
    @assert maximum(heat_rate_history) <= limit_heat_rate_w_cm2 + 1.0e-6
    @assert maximum(area_slew_history) <= limit_area_slew_m2_s + 1.0e-6
    @assert maximum(area_history) - minimum(area_history) > 1.0e-6
    println("mpc_constraint_audit = ", (
        commanded_area_range_m2=extrema(area_history),
        distinct_area_commands=length(unique(area_history)),
        maximum_drag_n=maximum(drag_history),
        maximum_heat_rate_w_cm2=maximum(heat_rate_history),
        maximum_area_slew_m2_s=maximum(area_slew_history),
    ))
end
println("COMPUTATIONAL TIME = $(t) s")
println("mpc_runtime_status = ", (
    solver_status=mpc_state.solver_status,
    last_qp_status=mpc_state.last_solution === nothing ? :none : mpc_state.last_solution.solver_status,
    predicted_terminal_energy=mpc_state.predicted_terminal_energy,
    predicted_terminal_energy_units=mpc_mode isa TargetEnergyMode ? :MJ_kg : :J_kg,
    active_limit=mpc_state.active_limit,
    last_error=mpc_state.last_error,
))
