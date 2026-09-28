# Planet-selectable aerobraking campaign with maximum-depletion and
# energy-targeting MPC.
# Set SPACEAGORA_MPC_PLANET to earth, mars, or venus.
# SPACEAGORA_EXAMPLE_SMOKE=1 selects a shortened integration-test scenario.
include(joinpath(@__DIR__, "common.jl"))

using CSV
using DataFrames

# ---------------------------------------------------------------------------
# Mission and environment configuration
# ---------------------------------------------------------------------------
smoke_mode = get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1"
planet_name = Symbol(lowercase(get(ENV, "SPACEAGORA_MPC_PLANET", "earth")))
planet_name in (:earth, :mars, :venus) || throw(ArgumentError(
    "SPACEAGORA_MPC_PLANET must be earth, mars, or venus; got $(repr(planet_name))."))

planet_settings = if planet_name === :earth
    (
        planet=Earth("", SPICE_PATH),
        initial_time=InitialTime(year=2020, month=1, day=1),
        periapsis_altitude_m=125.0e3,
        initial_apoapsis_altitude_m=60_500.0e3,
        target_apoapsis_altitude_km="10000.0",
        inclination_deg=45.0,
        raan_deg=0.0,
        argument_of_periapsis_deg=0.0,
        interface_altitude_km=300.0,
        bus_dims=(2.05, 2.05, 2.8),
        panel_area_m2=5.70,
        bus_mass_kg=620.0,
        panel_mass_each_kg=10.0,
        panel_offset_y_m=2.05 / 2.0 + 5.70 / 4.0,
        prop_mass_kg=200.0,
        bus_ram_face=:legacy,
        heat_rate_limit_w_cm2=20.0,
        heat_load_limit_j_cm2=30.0,
        drag_limit_n=500.0,
        gravity_file="EarthGGM05C.csv",
        gravity_degree=50,
        third_bodies=("Sun", "Moon"),
    )
elseif planet_name === :mars
    (
        planet=Mars("", SPICE_PATH),
        initial_time=InitialTime(year=2020, month=1, day=1),
        periapsis_altitude_m=100.0e3,
        initial_apoapsis_altitude_m=6_000.0e3,
        target_apoapsis_altitude_km="1500.0",
        inclination_deg=93.6,
        raan_deg=0.0,
        argument_of_periapsis_deg=0.0,
        interface_altitude_km=160.0,
        bus_dims=(1.87, 2.0, 2.0),
        panel_area_m2=7.26,
        bus_mass_kg=441.0,
        prop_mass_kg=200.0,
        bus_ram_face=:legacy,
        panel_mass_each_kg=10.0,
        panel_offset_y_m=1.0 + 7.26 / 4.0,
        heat_rate_limit_w_cm2=0.15,
        heat_load_limit_j_cm2=30.0,
        drag_limit_n=12.1,
        gravity_file="Mars50c.csv",
        gravity_degree=50,
        third_bodies=("Sun",),
    )
else
    (
        planet=Venus("", SPICE_PATH),
        initial_time=InitialTime(year=2014, month=5, day=19),
        periapsis_altitude_m=136.0e3,
        initial_apoapsis_altitude_m=72_649.0e3 - Venus().Rp_e,
        target_apoapsis_altitude_km="10000.0",
        inclination_deg=89.876,
        raan_deg=104.115,
        argument_of_periapsis_deg=75.505,
        interface_altitude_km=200.0,
        bus_dims=(1.0, 2.0, 2.87),
        prop_mass_kg=0.0,
        bus_ram_face=:frontal,
        panel_area_m2=5.70,
        bus_mass_kg=630.0,
        panel_mass_each_kg=10.0,
        panel_offset_y_m=2.425,
        heat_rate_limit_w_cm2=0.29,
        heat_load_limit_j_cm2=10.0,
        drag_limit_n=7.55,
        gravity_file="MGNP180U.csv",
        gravity_degree=50,
        third_bodies=("Sun",),
    )
end

planet = planet_settings.planet
initial_time = planet_settings.initial_time
ephemerides_model = SpiceEphemeridesModel()
density_model = PolynomialFitAtmosphereModel(planet)

periapsis_altitude_m = planet_settings.periapsis_altitude_m
initial_apoapsis_altitude_m = planet_settings.initial_apoapsis_altitude_m
target_apoapsis_altitude_m = 1.0e3 * parse(Float64, get(
    ENV, "SPACEAGORA_CAMPAIGN_TARGET_APOAPSIS_KM",
    planet_settings.target_apoapsis_altitude_km))
campaign_orbits = parse(Int, get(ENV, "SPACEAGORA_CAMPAIGN_ORBITS",
    smoke_mode ? "1" : "50"))
minimum_depletion_passes = smoke_mode ? 1 : 2
environment_interface_altitude_km = planet_settings.interface_altitude_km
mpc_inbound_solve_altitude_m = 1.0e3 * environment_interface_altitude_km

target_periapsis_radius_m = planet.Rp_e + periapsis_altitude_m
target_apoapsis_radius_m = planet.Rp_e + target_apoapsis_altitude_m
target_energy_mj_kg = -planet.μ /
    (target_apoapsis_radius_m + target_periapsis_radius_m) / 1.0e6

# ---------------------------------------------------------------------------
# Constraint, actuator, and QP configuration
# ---------------------------------------------------------------------------
active_constraints = mpc_constraints(:heat_rate, :drag)
limit_heat_rate_w_cm2 = planet_settings.heat_rate_limit_w_cm2
limit_heat_load_j_cm2 = planet_settings.heat_load_limit_j_cm2
limit_drag_n = planet_settings.drag_limit_n
limit_area_slew_m2_s = 0.20
drag_coefficient = 2.2

med_area_weight = 1.0e-5
med_area_slew_weight = 0.0
med_slack_weight = 1.0e3
med_energy_weight = 1.0
target_area_weight = 0.0
target_area_slew_weight = 0.0
target_slack_weight = 1.0e5
target_energy_weight = 5.0e-10
osqp_eps_abs = 1.0e-6
osqp_eps_rel = 1.0e-6
osqp_max_iter = 10_000

spacecraft_index = 1
controlled_panel_links = (2, 3)
min_panel_alpha_rad = 1.0e-4
max_panel_alpha_rad = pi / 2
control_dt_s = 1.0
solve_interval_s = 10.0            # minimum separation; interface trigger gives one solve/pass
bracket_tolerance_mj_kg = 1.0e-5
prediction_latitude_rad = 0.0
prediction_longitude_rad = 0.0
prediction_wind = false
reference_delta_s = 1.7e-7
reference_max_coast_steps = 2_000_000
reference_max_pass_steps = 20_000
qp_max_nodes = smoke_mode ? 40 : 120

ENV["SPACEAGORA_SOLVER_MODE"] = get(ENV, "SPACEAGORA_SOLVER_MODE", "split_imex")
ENV["SPACEAGORA_SPLIT_IMEX_SOLVER"] = get(ENV, "SPACEAGORA_SPLIT_IMEX_SOLVER", "kencarp4")
ENV["SPACEAGORA_VACUUM_GRAM_CACHE"] = get(ENV, "SPACEAGORA_VACUUM_GRAM_CACHE", "1")

spacecraft = make_three_body_spacecraft(
    bus_dims=planet_settings.bus_dims,
    panel_dims=(0.01, planet_settings.panel_area_m2 / 2.0, 1.0),
    bus_mass=planet_settings.bus_mass_kg,
    panel_mass_each=planet_settings.panel_mass_each_kg,
    panel_offset_y=planet_settings.panel_offset_y_m,
    bus_ram_face=planet_settings.bus_ram_face,
    ic=InitialCondition(
        planet;
        ra=initial_apoapsis_altitude_m,
        hp=periapsis_altitude_m,
        i=planet_settings.inclination_deg,
        ω=planet_settings.argument_of_periapsis_deg,
        Ω=planet_settings.raan_deg,
        ν=smoke_mode ? 345.0 : 180.0,
        initial_time=initial_time,
        ephemerides_model=ephemerides_model,
    ),
    prop_mass=planet_settings.prop_mass_kg,
    id=1,
)

harmonics_file = joinpath(
    REPO_ROOT, "data", "Gravity_harmonics_data", planet_settings.gravity_file)
dynamic_effectors = (
    NBodyGravityModel(body_names=planet_settings.third_bodies,
        primary_body_name=planet.name, planet=planet),
    GravitationalHarmonicsModel(planet_settings.gravity_degree,
        planet_settings.gravity_degree, harmonics_file, planet),
    SolarRadiationPressureModel(spacecraft.root.reflection_coefficient,
        spacecraft.root.ref_area),
    AerodynamicCoefficientfM(),
)

orbital_period_s = 2pi * sqrt(spacecraft.initial_condition.a^3 / planet.μ)
mission_time_cap_s = smoke_mode ? 1_800.0 : campaign_orbits * orbital_period_s
base_args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=mission_time_cap_s,
    initial_time=initial_time,
    dynamic_effectors=dynamic_effectors,
    density_model=density_model,
    ephemerides_model=ephemerides_model,
    orientation_sim=false,
    keplerian=false,
    EI_km=environment_interface_altitude_km,
    verbose=true,
    results_directory=joinpath(
        REPO_ROOT, "output", "mpc_$(planet_name)_campaign"),
)

function campaign_mpc_config(mode; area_weight, area_slew_weight, slack_weight)
    return mpc_config_from_spaceagora(
        base_args;
        mode=mode,
        spacecraft_index=spacecraft_index,
        controlled_panel_links=controlled_panel_links,
        constraints=active_constraints,
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
        max_depletion_energy_weight=med_energy_weight,
        osqp_eps_abs=osqp_eps_abs,
        osqp_eps_rel=osqp_eps_rel,
        osqp_max_iter=osqp_max_iter,
    )
end

med_config = campaign_mpc_config(
    MaxEnergyDepletionMode();
    area_weight=med_area_weight,
    area_slew_weight=med_area_slew_weight,
    slack_weight=med_slack_weight,
)
target_config = campaign_mpc_config(
    TargetEnergyMode();
    area_weight=target_area_weight,
    area_slew_weight=target_area_slew_weight,
    slack_weight=target_slack_weight,
)
reference_config = AerobrakingMPCReferenceConfig(
    h_cut_m=mpc_inbound_solve_altitude_m,
    delta_s=reference_delta_s,
    max_coast_steps=reference_max_coast_steps,
    max_pass_steps=reference_max_pass_steps,
)

function pass_controller(config)
    return AerobrakingMPCControlModel(
        config=config,
        state=AerobrakingMPCState(),
        spacecraft_index=spacecraft_index,
        controlled_panel_links=controlled_panel_links,
        control_dt_s=control_dt_s,
        min_alpha_rad=min_panel_alpha_rad,
        max_alpha_rad=max_panel_alpha_rad,
        solve_interval_s=solve_interval_s,
        build_reference_on_tick=true,
        qp_max_nodes=qp_max_nodes,
        reference=reference_config,
        fallback_area_m2=config.bus_reference_area_m2 + config.controllable_area_m2,
        prediction_latitude_rad=prediction_latitude_rad,
        prediction_longitude_rad=prediction_longitude_rad,
        prediction_wind=prediction_wind,
        solve_trigger_altitude_m=mpc_inbound_solve_altitude_m,
    )
end

campaign_control = AerobrakingMPCCampaignControlModel(
    maximum_depletion_control=pass_controller(med_config),
    targeting_control=pass_controller(target_config),
    state=AerobrakingMPCCampaignState(),
    minimum_depletion_passes=minimum_depletion_passes,
    bracket_tolerance_mj_kg=bracket_tolerance_mj_kg,
)

args = SimulationConfiguration(
    file_paths=base_args.file_paths,
    simulation_settings=base_args.simulation_settings,
    mission_configuration=MissionConfiguration(
        mission_type=smoke_mode ? MissionTime : MissionOrbits,
        keplerian=false,
        number_of_orbits=campaign_orbits,
        mission_time=mission_time_cap_s,
        orientation_sim=false,
        num_steps_to_save=smoke_mode ? 900 : 25_000,
        data_rate=control_dt_s,
    ),
    environment_model=base_args.environment_model,
    dynamics_model=base_args.dynamics_model,
    guidance_model=base_args.guidance_model,
    navigation_model=base_args.navigation_model,
    control_model=ControlModel(control_effectors=(campaign_control,), control_rates=[control_dt_s]),
    initial_time=base_args.initial_time,
    integration_tolerances=IntegrationTolerances(
        reltol_orbit=1e-8,
        abstol_orbit=1e-8,
        dt_max_orbit=30.0,
        reltol_atmosphere=1e-8,
        abstol_atmosphere=1e-8,
        dt_max_atmosphere=5.0,
    ),
    solver_config=SolverConfig(solver_mode=Symbol(ENV["SPACEAGORA_SOLVER_MODE"])),
)

println("campaign_settings = ", (
    planet=typeof(planet),
    density_model=typeof(density_model),
    initial_apoapsis_altitude_km=initial_apoapsis_altitude_m / 1e3,
    periapsis_altitude_km=periapsis_altitude_m / 1e3,
    target_apoapsis_altitude_km=target_apoapsis_altitude_m / 1e3,
    target_energy_mj_kg=target_energy_mj_kg,
    constraints=constraint_names(active_constraints),
    minimum_depletion_passes=minimum_depletion_passes,
    inbound_solve_altitude_km=mpc_inbound_solve_altitude_m / 1e3,
))

save_fields = vcat(default_save_fields(args), collect(mpc_campaign_save_fields(campaign_control)))
elapsed_s = @elapsed run_simulation(args; save_fields=save_fields, isolate_state=false)
csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
if args.simulation_settings.results && isfile(csv_path)
    df = CSV.read(csv_path, DataFrame)
    area = df.sc1_mpc_commanded_area_m2
    drag = filter(isfinite, df.sc1_mpc_drag_n)
    heat_rate = filter(isfinite, df.sc1_mpc_heat_rate_w_cm2)
    slew = abs.(diff(area)) ./ diff(df.time)
    constraint_active(active_constraints, :drag) &&
        @assert maximum(drag) <= limit_drag_n + 1.0e-6
    constraint_active(active_constraints, :heat_rate) &&
        @assert maximum(heat_rate) <= limit_heat_rate_w_cm2 + 1.0e-6
    constraint_active(active_constraints, :slew) &&
        @assert maximum(slew) <= limit_area_slew_m2_s + 1.0e-6
    println("campaign_constraint_audit = ", (
        samples=nrow(df),
        area_range_m2=extrema(area),
        maximum_drag_n=maximum(drag),
        maximum_heat_rate_w_cm2=maximum(heat_rate),
        maximum_area_slew_m2_s=maximum(slew),
    ))
end
println("COMPUTATIONAL TIME = $(elapsed_s) s")
println("campaign_final_state = ", (
    phase=campaign_control.state.phase,
    completed_passes=campaign_control.state.completed_passes,
    bracket=(campaign_control.state.bracket_min_energy_mj_kg,
        campaign_control.state.bracket_max_energy_mj_kg),
    switch_time_s=campaign_control.state.switch_time_s,
    med_solves=campaign_control.maximum_depletion_control.state.solve_count,
    targeting_solves=campaign_control.targeting_control.state.solve_count,
    last_error=campaign_control.state.last_error,
))
