include(joinpath(@__DIR__, "common.jl"))

using CSV
using DataFrames
using LinearAlgebra
using StaticArrays

setup_gram_example!()

# ABTS VEX control-test values.
planet = Venus("", SPICE_PATH)
initial_time = InitialTime(year=2014, month=5, day=19, hour=4, minute=0, second=0.0)
ephemerides_model = SpiceEphemeridesModel()
limit_qdot = 0.29
limit_q = 15.0
limit_dyn_pressure = 0.3
ra_initial_m = 72_649e3
rp_initial_m = planet.Rp_e + 136e3
legacy_ra_target_m = 72_480e3
ra_target_m = parse(Float64, get(ENV, "SPACEAGORA_VEX_TARGET_RA_M", string(legacy_ra_target_m)))
orbit_period_s = 2pi * sqrt(((ra_initial_m + rp_initial_m) / 2)^3 / planet.μ)

ic = InitialCondition(
    ra=ra_initial_m,
    rp=rp_initial_m,
    i=89.876,
    ω=75.505,
    Ω=104.115,
    ν=180.0001,
)

spacecraft = make_three_body_spacecraft(
    bus_dims=(2.05, 3.7, 2.8),
    panel_dims=(0.01, 5.7 / 2.0, 1.0),
    bus_mass=620.0,
    panel_mass_each=10.0,
    panel_offset_y=2.05 / 2.0 + 5.7 / 4.0,
    ic=ic,
    reflection_coefficient=0.9,
    prop_mass=10.0,
    id=102,
)
@assert spacecraft.root.ref_area == 2.05 * 2.8
@assert spacecraft.links[2].ref_area + spacecraft.links[3].ref_area == 5.7 * 1.0

energy_state = AerobrakingEnergyDepletionState(num_sats=1)
energy_config = AerobrakingEnergyDepletionConfig(
    guidance_modes=(:targeting, :max_energy_depletion),
    max_energy_submodes=(:heat_rate, :heat_load),
    heat_load_switch_solver=:closed_form,
    controlled_panel_links=(2, 3),
    target_apoapsis_radius_m=ra_target_m,
    max_alpha_rad=pi / 2,
    min_alpha_rad=0.0,
    heat_rate_limit_w_cm2=limit_qdot,
    heat_load_limit_j_cm2=limit_q,
    structural_load_limit_pa=limit_dyn_pressure,
)
guidance = AerobrakingEnergyDepletionGuidanceModel(energy_config, energy_state)
control = AerobrakingEnergyDepletionControlModel(energy_config, energy_state)

sun_gravity = NBodyGravityModel(
    body_names=("Sun",),
    primary_body_name="Venus",
    planet=planet,
)
base_args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=orbit_period_s,
    initial_time=initial_time,
    dynamic_effectors=(
        InverseSquaredJ2GravityModel(),
        sun_gravity,
        AerodynamicCoefficientfM(),
    ),
    density_model=GRAMAtmosphereModel(planet_name="venus"),
    ephemerides_model=ephemerides_model,
    orientation_sim=false,
    keplerian=false,
    EI_km=200.0,
    results_directory=joinpath(
        REPO_ROOT,
        "output",
        "vex_energy_depletion",
        "targeting_nbody_$(round(Int, ra_target_m))m",
    ),
)

environment = EnvironmentModel(
    planet=base_args.environment_model.planet,
    EI=base_args.environment_model.EI,
    density_model=base_args.environment_model.density_model,
    ephemerides_model=base_args.environment_model.ephemerides_model,
    thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
    topography=false,
    topo_degree=0,
    topo_order=0,
    wind=true,
)
args = SimulationConfiguration(
    file_paths=base_args.file_paths,
    simulation_settings=base_args.simulation_settings,
    mission_configuration=MissionConfiguration(
        mission_type=base_args.mission_configuration.mission_type,
        keplerian=false,
        number_of_orbits=1,
        mission_time=orbit_period_s,
        orientation_sim=false,
        num_steps_to_save=base_args.mission_configuration.num_steps_to_save,
        data_rate=1.0,
    ),
    environment_model=environment,
    dynamics_model=base_args.dynamics_model,
    guidance_model=GuidanceModel(guidance_effectors=(guidance,), guidance_rates=[1 / 3]),
    navigation_model=base_args.navigation_model,
    control_model=ControlModel(control_effectors=(control,), control_rates=[0.1]),
    initial_time=base_args.initial_time,
    integration_tolerances=base_args.integration_tolerances,
    solver_config=SolverConfig(maxiters=2_000_000),
)

target_energy_save_field = SaveField(
    :edg_target_energy_jkg,
    (u, t, integrator) -> [
        integrator.p.args.control_model.control_effectors[1].state.target_energy_jkg[i]
        for i in eachindex(integrator.p.args.dynamics_model.spacecraft)
    ];
    per_satellite=true,
    column_prefix="edg_target_energy_jkg",
)
targeting_active_save_field = SaveField(
    :edg_targeting_active,
    (u, t, integrator) -> [
        integrator.p.args.control_model.control_effectors[1].state.targeting_active[i] ? 1.0 : 0.0
        for i in eachindex(integrator.p.args.dynamics_model.spacecraft)
    ];
    per_satellite=true,
    column_prefix="edg_targeting_active",
)
targeting_switch_save_field = SaveField(
    :edg_targeting_switch_s,
    (u, t, integrator) -> [
        integrator.p.args.control_model.control_effectors[1].state.targeting_switch_s[i]
        for i in eachindex(integrator.p.args.dynamics_model.spacecraft)
    ];
    per_satellite=true,
    column_prefix="edg_targeting_switch_s",
)
perturbation_energy_save_field = SaveField(
    :edg_target_perturbation_energy_change_jkg,
    (u, t, integrator) -> [
        integrator.p.args.control_model.control_effectors[1].state.target_perturbation_energy_change_jkg[i]
        for i in eachindex(integrator.p.args.dynamics_model.spacecraft)
    ];
    per_satellite=true,
    column_prefix="edg_target_perturbation_energy_change_jkg",
)
bracket_min_save_field = SaveField(
    :edg_bracket_min_energy_jkg,
    (u, t, integrator) -> [
        integrator.p.args.control_model.control_effectors[1].state.bracket_min_energy_jkg[i]
        for i in eachindex(integrator.p.args.dynamics_model.spacecraft)
    ];
    per_satellite=true,
    column_prefix="edg_bracket_min_energy_jkg",
)
bracket_max_save_field = SaveField(
    :edg_bracket_max_energy_jkg,
    (u, t, integrator) -> [
        integrator.p.args.control_model.control_effectors[1].state.bracket_max_energy_jkg[i]
        for i in eachindex(integrator.p.args.dynamics_model.spacecraft)
    ];
    per_satellite=true,
    column_prefix="edg_bracket_max_energy_jkg",
)

save_fields = vcat(
    default_save_fields(args),
    [
        target_energy_save_field,
        targeting_active_save_field,
        targeting_switch_save_field,
        perturbation_energy_save_field,
        bracket_min_save_field,
        bracket_max_save_field,
    ],
)

elapsed = @elapsed run_simulation(args; save_fields=save_fields)
csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
df = CSV.read(csv_path, DataFrame)
final_position = SVector{3, Float64}(df.sc1_pos_1[end], df.sc1_pos_2[end], df.sc1_pos_3[end])
final_velocity = SVector{3, Float64}(df.sc1_vel_1[end], df.sc1_vel_2[end], df.sc1_vel_3[end])
final_elements = SimulationModel.ControlHooks.rvtoorbitalelement(final_position, final_velocity, planet)
final_apoapsis_m = final_elements[1] * (1.0 + final_elements[2])
final_energy_jkg = dot(final_velocity, final_velocity) / 2.0 - planet.μ / norm(final_position)
target_energy_jkg = df.sc1_edg_target_energy_jkg[end]
perturbation_energy_change_jkg = df.sc1_edg_target_perturbation_energy_change_jkg[end]
nominal_target_energy_jkg = target_energy_jkg + perturbation_energy_change_jkg

println("Saved $(nrow(df)) samples to $(abspath(csv_path))")
println("Targeting active = $(df.sc1_edg_targeting_active[end] > 0.5)")
println("Targeting switch = $(df.sc1_edg_targeting_switch_s[end]) s")
println("Uncorrected two-body target energy = $(nominal_target_energy_jkg) J/kg")
println("N-body/J2-corrected target energy = $(target_energy_jkg) J/kg")
println("Exit-to-apoapsis perturbation energy change = $(perturbation_energy_change_jkg) J/kg")
println("Targeting energy bracket = ($(df.sc1_edg_bracket_min_energy_jkg[end]), $(df.sc1_edg_bracket_max_energy_jkg[end])) J/kg")
println("Final specific orbital energy = $(final_energy_jkg) J/kg")
println("Target apoapsis radius = $(ra_target_m) m")
println("Final apoapsis radius = $(final_apoapsis_m) m")
println("Apoapsis error = $(final_apoapsis_m - ra_target_m) m")
println("COMPUTATIONAL TIME = $(elapsed) s")
