#=
Planet-selectable maximum-energy-depletion MPC constraint cases.
Set SPACEAGORA_MPC_PLANET to earth, mars, or venus.

Case I   : heat-rate + drag-force constraints
Case II  : Case I + heat-load constraint
Case III : Case II + actuator area-slew constraint

The selected SpaceAGORA planet supplies the gravity constants and polynomial-fit
atmosphere. The script saves the MPC prediction and nonlinear KS verification.
=#

include(joinpath(@__DIR__, "common.jl"))

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

using CSV
using DataFrames
using LinearAlgebra
using Plots
using Printf

# =============================================================================
# CASE SETTINGS
# =============================================================================
const MPC_PLANET_NAME = Symbol(lowercase(get(ENV, "SPACEAGORA_MPC_PLANET", "earth")))
MPC_PLANET_NAME in (:earth, :mars, :venus) || throw(ArgumentError(
    "SPACEAGORA_MPC_PLANET must be earth, mars, or venus; got $(repr(MPC_PLANET_NAME))."))

const PLANET_SETTINGS = if MPC_PLANET_NAME === :earth
    (
        planet=Earth(),
        initial_time=InitialTime(year=2020, month=1, day=1),
        periapsis_altitude_m=110.0e3,
        apoapsis_radius_m=56_378.0e3,
        inclination_deg=89.876,
        raan_deg=104.115,
        argument_of_periapsis_deg=75.505,
        spacecraft_mass_kg=650.0,
        bus_area_m2=5.74,
        solar_panel_area_m2=5.70,
        heat_rate_limit_w_cm2=0.29,
        heat_load_limit_j_cm2=40.0,
        drag_force_limit_n=7.55,
        interface_altitude_m=300.0e3,
    )
elseif MPC_PLANET_NAME === :mars
    (
        planet=Mars(),
        initial_time=InitialTime(year=2020, month=1, day=1),
        periapsis_altitude_m=95.0e3,
        apoapsis_radius_m=9_396.2e3,
        inclination_deg=93.6,
        raan_deg=0.0,
        argument_of_periapsis_deg=0.0,
        spacecraft_mass_kg=461.0,
        bus_area_m2=3.74,
        solar_panel_area_m2=7.26,
        heat_rate_limit_w_cm2=0.15,
        heat_load_limit_j_cm2=30.0,
        drag_force_limit_n=12.1,
        interface_altitude_m=160.0e3,
    )
else
    (
        planet=Venus(),
        initial_time=InitialTime(year=2014, month=5, day=19),
        periapsis_altitude_m=136.0e3,
        apoapsis_radius_m=72_649.0e3,
        inclination_deg=89.876,
        raan_deg=104.115,
        argument_of_periapsis_deg=75.505,
        spacecraft_mass_kg=650.0,
        bus_area_m2=5.74,
        solar_panel_area_m2=5.70,
        heat_rate_limit_w_cm2=0.29,
        heat_load_limit_j_cm2=10.0,
        drag_force_limit_n=7.55,
        interface_altitude_m=200.0e3,
    )
end

const PLANET = PLANET_SETTINGS.planet
const INITIAL_TIME = PLANET_SETTINGS.initial_time
const PERIAPSIS_ALTITUDE_M = PLANET_SETTINGS.periapsis_altitude_m
const APOAPSIS_RADIUS_M = PLANET_SETTINGS.apoapsis_radius_m
const INCLINATION_DEG = PLANET_SETTINGS.inclination_deg
const RAAN_DEG = PLANET_SETTINGS.raan_deg
const ARGUMENT_OF_PERIAPSIS_DEG = PLANET_SETTINGS.argument_of_periapsis_deg
const INITIAL_TRUE_ANOMALY_DEG = 180.0
const INITIAL_CONDITION = InitialCondition(
    ra=APOAPSIS_RADIUS_M,
    rp=PLANET.Rp_e + PERIAPSIS_ALTITUDE_M,
    i=INCLINATION_DEG,
    Ω=RAAN_DEG,
    ω=ARGUMENT_OF_PERIAPSIS_DEG,
    ν=INITIAL_TRUE_ANOMALY_DEG,
)

const SPACECRAFT_MASS_KG = PLANET_SETTINGS.spacecraft_mass_kg
const BUS_AREA_M2 = PLANET_SETTINGS.bus_area_m2
const SOLAR_PANEL_AREA_M2 = PLANET_SETTINGS.solar_panel_area_m2
const DRAG_COEFFICIENT = 2.2

const HEAT_RATE_LIMIT_W_CM2 = PLANET_SETTINGS.heat_rate_limit_w_cm2
const HEAT_LOAD_LIMIT_J_CM2 = PLANET_SETTINGS.heat_load_limit_j_cm2
const DRAG_FORCE_LIMIT_N = PLANET_SETTINGS.drag_force_limit_n
const AREA_SLEW_LIMIT_M2_S = 0.20

const HEAT_RATE_QP_BACKOFF_W_CM2 = 1.0e-3
const HEAT_LOAD_QP_BACKOFF_J_CM2 = 5.0e-2
const DRAG_FORCE_QP_BACKOFF_N = 2.0e-2
const AREA_SLEW_QP_BACKOFF_M2_S = 2.0e-3

const AREA_WEIGHT = 1.0e-5
const AREA_SLEW_WEIGHT = 0.0
const SLACK_WEIGHT = 1.0e-2
const MAX_DEPLETION_ENERGY_WEIGHT = 1.0
const OSQP_EPS_ABS = 1.0e-5
const OSQP_EPS_REL = 1.0e-5
const OSQP_MAX_ITER = 100_000

const ATMOSPHERIC_CUTOFF_ALTITUDE_M = PLANET_SETTINGS.interface_altitude_m
const KS_FICTITIOUS_TIME_STEP = 1.7e-7
const MAX_COAST_STEPS = 2_000_000
const MAX_PASS_STEPS = 20_000
const QP_MAX_NODES = nothing

const DENSITY_SELECTION = :planet_polyfit
const OUTPUT_DIRECTORY = joinpath(
    @__DIR__, "..", "output", "mpc_$(MPC_PLANET_NAME)_med_constraint_cases")
const WRITE_PDF = true
const RUN_SPACEAGORA_CALLBACK_CASES = lowercase(get(
    ENV, "SPACEAGORA_MPC_RUN_CALLBACK_CASES", "true")) in ("1", "true", "yes")
const CONTROL_SAMPLE_TIME_S = 0.5
const CALLBACK_END_MARGIN_S = 5.0
const CARTESIAN_RELATIVE_TOLERANCE = 1.0e-9
const CARTESIAN_ABSOLUTE_TOLERANCE = 1.0e-9
const CARTESIAN_MAX_STEP_S = 0.2

const CASES = (
    (name=:case_I, label="Case I", constraints=mpc_constraints(:heat_rate, :drag)),
    (name=:case_II, label="Case II", constraints=mpc_constraints(:heat_rate, :heat_load, :drag)),
    (name=:case_III, label="Case III", constraints=mpc_constraints(:heat_rate, :heat_load, :drag, :slew)),
)

# =============================================================================
# ATMOSPHERE SELECTION
# =============================================================================
const DENSITY_MODEL = PolynomialFitAtmosphereModel(PLANET)
const DENSITY_CONTEXT = (
    args=(environment_model=(planet=PLANET, density_model=DENSITY_MODEL),),)
const DENSITY_FOR_MPC = density_function_from_spaceagora(DENSITY_CONTEXT)

# =============================================================================
# CASE CONSTRUCTION AND NONLINEAR VERIFICATION
# =============================================================================
function base_config()
    return AerobrakingMPCConfig(
        mode=MaxEnergyDepletionMode(),
        bus_reference_area_m2=BUS_AREA_M2,
        controllable_area_m2=SOLAR_PANEL_AREA_M2,
        mass_kg=SPACECRAFT_MASS_KG,
        drag_coefficient=DRAG_COEFFICIENT,
        qdot_max_w_cm2=HEAT_RATE_LIMIT_W_CM2 - HEAT_RATE_QP_BACKOFF_W_CM2,
        heat_load_max_j_cm2=HEAT_LOAD_LIMIT_J_CM2 - HEAT_LOAD_QP_BACKOFF_J_CM2,
        drag_max_n=DRAG_FORCE_LIMIT_N - DRAG_FORCE_QP_BACKOFF_N,
        area_slew_max_m2_s=AREA_SLEW_LIMIT_M2_S - AREA_SLEW_QP_BACKOFF_M2_S,
        use_constraints=true,
        use_slew_constraint=false,
        use_qdot_constraint=true,
        use_heat_load_constraint=false,
        use_drag_constraint=true,
        target_energy_mj_kg=0.0,
        area_weight=AREA_WEIGHT,
        area_slew_weight=AREA_SLEW_WEIGHT,
        slack_weight=SLACK_WEIGHT,
        target_energy_weight=0.0,
        max_depletion_energy_weight=MAX_DEPLETION_ENERGY_WEIGHT,
        osqp_eps_abs=OSQP_EPS_ABS,
        osqp_eps_rel=OSQP_EPS_REL,
        osqp_max_iter=OSQP_MAX_ITER,
    )
end

function refinement_comparison(case_label, baseline, refined)
    sample(field) = [interpolate_mpc_history(
        refined.absolute_time_s, getproperty(refined, field), time)
        for time in baseline.absolute_time_s]
    return (
        case=case_label,
        refined_delta_s=0.5 * KS_FICTITIOUS_TIME_STEP,
        maximum_altitude_difference_m=maximum(abs.(
            baseline.altitude_m .- sample(:altitude_m))),
        maximum_drag_difference_n=maximum(abs.(
            baseline.drag_n .- sample(:drag_n))),
        maximum_heat_rate_difference_w_cm2=maximum(abs.(
            baseline.heat_rate_w_cm2 .- sample(:heat_rate_w_cm2))),
        final_heat_load_difference_j_cm2=abs(
            last(baseline.heat_load_j_cm2) - last(refined.heat_load_j_cm2)),
        final_energy_difference_mj_kg=abs(
            last(baseline.energy_mj_kg) - last(refined.energy_mj_kg)),
    )
end

function case_summary(case, solution, solve_time_s, prediction, rollout, problem)
    predicted_heat_load = cumulative_mpc_heat_load(
        prediction[:, 3] ./ 1.0e4, problem.t)
    predicted_slew = vcat(NaN, diff(solution.area_m2) ./ diff(problem.t))
    finite_predicted_slew = filter(isfinite, abs.(predicted_slew))
    finite_rollout_slew = filter(isfinite, abs.(rollout.slew_m2_s))
    return (
        case=case.label,
        constraints=join(String.(constraint_names(case.constraints)), "+"),
        solver_status=String(solution.solver_status),
        solve_time_s,
        minimum_commanded_area_m2=minimum(solution.area_m2),
        maximum_commanded_area_m2=maximum(solution.area_m2),
        predicted_maximum_heat_rate_w_cm2=maximum(prediction[:, 3]) / 1.0e4,
        nonlinear_maximum_heat_rate_w_cm2=maximum(rollout.heat_rate_w_cm2),
        predicted_final_heat_load_j_cm2=last(predicted_heat_load),
        nonlinear_final_heat_load_j_cm2=last(rollout.heat_load_j_cm2),
        predicted_maximum_drag_n=maximum(prediction[:, 2]),
        nonlinear_maximum_drag_n=maximum(rollout.drag_n),
        predicted_maximum_area_slew_m2_s=maximum(finite_predicted_slew),
        nonlinear_maximum_area_slew_m2_s=maximum(finite_rollout_slew),
        predicted_final_energy_mj_kg=prediction[end, 4] / 1.0e6,
        nonlinear_final_energy_mj_kg=last(rollout.energy_mj_kg),
        maximum_altitude_prediction_error_m=maximum(abs.(
            prediction[:, 1] .- [rollout.altitude_m[argmin(abs.(
                rollout.absolute_time_s .- time))] for time in problem.t])),
        heat_rate_constraint_satisfied=
            maximum(rollout.heat_rate_w_cm2) <= HEAT_RATE_LIMIT_W_CM2,
        heat_load_constraint_satisfied=!constraint_active(case.constraints, :heat_load) ||
            last(rollout.heat_load_j_cm2) <= HEAT_LOAD_LIMIT_J_CM2,
        drag_constraint_satisfied=maximum(rollout.drag_n) <= DRAG_FORCE_LIMIT_N,
        slew_constraint_satisfied=!constraint_active(case.constraints, :slew) ||
            maximum(finite_rollout_slew) <= AREA_SLEW_LIMIT_M2_S,
    )
end

function output_dataframe(case_label, source, time_s, altitude_m, area_m2,
    heat_rate_w_cm2, heat_load_j_cm2, drag_n, slew_m2_s, energy_mj_kg)
    return DataFrame(
        case=fill(case_label, length(time_s)),
        source=fill(source, length(time_s)),
        time_s=time_s,
        altitude_km=altitude_m ./ 1.0e3,
        commanded_area_m2=area_m2,
        heat_rate_w_cm2=heat_rate_w_cm2,
        heat_load_j_cm2=heat_load_j_cm2,
        drag_force_n=drag_n,
        area_slew_m2_s=slew_m2_s,
        specific_energy_mj_kg=energy_mj_kg,
    )
end

function run_spaceagora_callback_case(case, reference)
    reference_initial = ks_state_to_cartesian(view(reference.states, 1, :))
    spacecraft = make_three_body_spacecraft(
        bus_dims=(1.0, 2.0, BUS_AREA_M2 / 2.0),
        panel_dims=(0.01, SOLAR_PANEL_AREA_M2 / 2.0, 1.0),
        bus_mass=SPACECRAFT_MASS_KG - 20.0,
        panel_mass_each=10.0,
        panel_offset_y=1.0 + SOLAR_PANEL_AREA_M2 / 4.0,
        ic=SpaceAGORA.SimulationModel.CartesianInitialCondition(
            reference_initial.position_ii_m,
            reference_initial.velocity_ii_m,
        ),
        prop_mass=0.0,
        id=1,
        bus_ram_face=:frontal,
    )
    density_model = DENSITY_MODEL
    results_directory = joinpath(OUTPUT_DIRECTORY, "callback_$(case.name)")
    base_args = make_example_config(
        planet=PLANET,
        spacecraft=spacecraft,
        mission_time=last(reference.time_s) - first(reference.time_s) +
            CALLBACK_END_MARGIN_S,
        initial_time=INITIAL_TIME,
        dynamic_effectors=(
            InverseSquaredJ2GravityModel(),
            AerodynamicCommandedAreaDragModel(
                drag_coefficient=DRAG_COEFFICIENT,
                controlled_link_indices=(2, 3),
            ),
        ),
        density_model=density_model,
        ephemerides_model=SimpleEphemeridesModel(),
        orientation_sim=false,
        keplerian=false,
        EI_km=ATMOSPHERIC_CUTOFF_ALTITUDE_M / 1.0e3,
        verbose=false,
        results_directory=results_directory,
    )
    config = mpc_config_from_spaceagora(
        base_args;
        mode=MaxEnergyDepletionMode(),
        spacecraft_index=1,
        controlled_panel_links=(2, 3),
        constraints=case.constraints,
        qdot_max_w_cm2=HEAT_RATE_LIMIT_W_CM2 - HEAT_RATE_QP_BACKOFF_W_CM2,
        heat_load_max_j_cm2=HEAT_LOAD_LIMIT_J_CM2 - HEAT_LOAD_QP_BACKOFF_J_CM2,
        drag_max_n=DRAG_FORCE_LIMIT_N - DRAG_FORCE_QP_BACKOFF_N,
        area_slew_max_m2_s=AREA_SLEW_LIMIT_M2_S - AREA_SLEW_QP_BACKOFF_M2_S,
        drag_coefficient=DRAG_COEFFICIENT,
        area_weight=AREA_WEIGHT,
        area_slew_weight=AREA_SLEW_WEIGHT,
        slack_weight=SLACK_WEIGHT,
        target_energy_mj_kg=0.0,
        target_energy_weight=0.0,
        max_depletion_energy_weight=MAX_DEPLETION_ENERGY_WEIGHT,
        osqp_eps_abs=OSQP_EPS_ABS,
        osqp_eps_rel=OSQP_EPS_REL,
        osqp_max_iter=OSQP_MAX_ITER,
    )
    @assert config.mass_kg ≈ SPACECRAFT_MASS_KG
    @assert config.bus_reference_area_m2 ≈ BUS_AREA_M2
    @assert config.controllable_area_m2 ≈ SOLAR_PANEL_AREA_M2

    state = AerobrakingMPCState()
    controller = AerobrakingMPCControlModel(
        config=config,
        state=state,
        spacecraft_index=1,
        controlled_panel_links=(2, 3),
        control_dt_s=CONTROL_SAMPLE_TIME_S,
        min_alpha_rad=0.0,
        max_alpha_rad=pi / 2,
        solve_interval_s=1.0e9,
        build_reference_on_tick=true,
        qp_max_nodes=QP_MAX_NODES,
        reference=AerobrakingMPCReferenceConfig(
            h_cut_m=ATMOSPHERIC_CUTOFF_ALTITUDE_M,
            delta_s=KS_FICTITIOUS_TIME_STEP,
            max_coast_steps=MAX_COAST_STEPS,
            max_pass_steps=MAX_PASS_STEPS,
        ),
        fallback_area_m2=BUS_AREA_M2 + SOLAR_PANEL_AREA_M2,
        prediction_latitude_rad=0.0,
        prediction_longitude_rad=0.0,
        prediction_wind=false,
        solve_trigger_altitude_m=ATMOSPHERIC_CUTOFF_ALTITUDE_M,
        heat_rate_model=:kinetic_energy_flux,
    )
    args = SimulationConfiguration(
        file_paths=base_args.file_paths,
        simulation_settings=base_args.simulation_settings,
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=false,
            number_of_orbits=1,
            mission_time=base_args.mission_configuration.mission_time,
            orientation_sim=false,
            num_steps_to_save=2_000,
            data_rate=CONTROL_SAMPLE_TIME_S,
        ),
        environment_model=base_args.environment_model,
        dynamics_model=base_args.dynamics_model,
        guidance_model=base_args.guidance_model,
        navigation_model=base_args.navigation_model,
        control_model=ControlModel(
            control_effectors=(controller,),
            control_rates=[CONTROL_SAMPLE_TIME_S],
        ),
        initial_time=base_args.initial_time,
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=CARTESIAN_RELATIVE_TOLERANCE,
            abstol_orbit=CARTESIAN_ABSOLUTE_TOLERANCE,
            dt_max_orbit=CARTESIAN_MAX_STEP_S,
            reltol_atmosphere=CARTESIAN_RELATIVE_TOLERANCE,
            abstol_atmosphere=CARTESIAN_ABSOLUTE_TOLERANCE,
            dt_max_atmosphere=CARTESIAN_MAX_STEP_S,
        ),
    )
    save_fields = vcat(default_save_fields(args),
        collect(mpc_control_save_fields(controller)))
    runtime_s = @elapsed run_simulation(
        args; save_fields=save_fields, isolate_state=false)
    csv_path = joinpath(results_directory, "simulation_results.csv")
    isfile(csv_path) || error("Callback simulation did not write $(csv_path).")
    table = CSV.read(csv_path, DataFrame)
    state.solve_count == 1 || error(
        "$(case.label) expected one interface solve, observed $(state.solve_count).")
    state.last_solution !== nothing || error(
        "$(case.label) callback did not retain an MPC solution.")
    state.last_solution.ok || error(
        "$(case.label) callback QP failed: $(state.last_solution.solver_status)")

    time_s = Float64.(table.time)
    area = Float64.(table.sc1_mpc_commanded_area_m2)
    alpha = Float64.(table.sc1_mpc_alpha_rad)
    controller_heat_rate = Float64.(table.sc1_mpc_heat_rate_w_cm2)
    controller_heat_load = Float64.(table.sc1_mpc_heat_load_j_cm2)
    controller_drag = Float64.(table.sc1_mpc_drag_n)
    cached_dynamics_drag = sqrt.(Float64.(table.sc1_drag_1).^2 .+
        Float64.(table.sc1_drag_2).^2 .+ Float64.(table.sc1_drag_3).^2)
    geodetic_altitude = Float64.(table.sc1_altitude)
    position = hcat(table.sc1_pos_1, table.sc1_pos_2, table.sc1_pos_3)
    velocity = hcat(table.sc1_vel_1, table.sc1_vel_2, table.sc1_vel_3)
    params = AerobrakingMPCParams(
        Re=PLANET.Rp_e, μ=PLANET.μ, J2=PLANET.J2, Ω=norm(PLANET.ω))
    # Evaluate constraints at the accepted Cartesian output states. The
    # dynamics drag cache is updated inside adaptive solver stages and can
    # therefore represent a rejected/internal stage rather than this CSV row.
    # It remains diagnostic telemetry below, but is not a truth constraint
    # audit. Use the simulation's geodetic altitude to sample the same planet
    # polyfit as the actual Cartesian aerodynamic effector.
    accepted_outputs = [evaluate_cartesian_mpc_outputs(
        view(position, index, :), view(velocity, index, :), time_s[index],
        area[index], params, config;
        density=(_, elapsed_time_s) -> getDensity(
            DENSITY_MODEL, geodetic_altitude[index], 0.0, 0.0,
            Float64(elapsed_time_s), false, DENSITY_CONTEXT)[1])
        for index in axes(position, 1)]
    accepted_altitude = getproperty.(accepted_outputs, :altitude_m)
    heat_rate = getproperty.(accepted_outputs, :heat_rate_w_cm2)
    heat_load = cumulative_mpc_heat_load(heat_rate, time_s)
    drag = getproperty.(accepted_outputs, :drag_n)
    energy = getproperty.(accepted_outputs, :specific_energy_mj_kg)
    slew = vcat(NaN, diff(area) ./ diff(time_s))
    reconstructed_area = [commanded_area_from_alpha(
        config, value; min_alpha_rad=0.0, max_alpha_rad=pi / 2)
        for value in alpha]
    mapping_errors = filter(isfinite, abs.(area .- reconstructed_area))
    altitude_reference_differences = filter(isfinite,
        abs.(geodetic_altitude .- accepted_altitude))
    controller_drag_errors = filter(isfinite, abs.(controller_drag .- drag))
    cached_drag_errors = filter(isfinite, abs.(cached_dynamics_drag .- drag))
    finite_heat_rate = filter(isfinite, heat_rate)
    finite_heat_load = filter(isfinite, heat_load)
    finite_drag = filter(isfinite, drag)
    area_mapping_error = maximum(mapping_errors)
    maximum_altitude_reference_difference = maximum(altitude_reference_differences)
    controller_drag_sampling_error = maximum(controller_drag_errors)
    cached_drag_sampling_error = maximum(cached_drag_errors)

    @assert area_mapping_error <= 1.0e-10
    @assert maximum(finite_heat_rate) <= HEAT_RATE_LIMIT_W_CM2
    @assert maximum(finite_drag) <= DRAG_FORCE_LIMIT_N
    constraint_active(case.constraints, :heat_load) &&
        @assert last(finite_heat_load) <= HEAT_LOAD_LIMIT_J_CM2
    constraint_active(case.constraints, :slew) &&
        @assert maximum(filter(isfinite, abs.(slew))) <= AREA_SLEW_LIMIT_M2_S

    raw_plan = DataFrame(
        time_s=state.plan_time_s,
        commanded_area_m2=state.plan_area_m2,
        commanded_alpha_rad=[alpha_from_commanded_area(
            config, value; min_alpha_rad=0.0, max_alpha_rad=pi / 2)
            for value in state.plan_area_m2],
    )
    CSV.write(joinpath(results_directory, "raw_mpc_plan.csv"), raw_plan)
    return (
        case=case,
        config=config,
        controller=controller,
        runtime_s=runtime_s,
        csv_path=csv_path,
        time_s=time_s,
        absolute_time_s=time_s,
        altitude_m=accepted_altitude,
        geodetic_altitude_m=geodetic_altitude,
        area_m2=area,
        alpha_rad=alpha,
        heat_rate_w_cm2=heat_rate,
        heat_load_j_cm2=heat_load,
        drag_n=drag,
        controller_heat_rate_w_cm2=controller_heat_rate,
        controller_heat_load_j_cm2=controller_heat_load,
        controller_drag_n=controller_drag,
        cached_dynamics_drag_n=cached_dynamics_drag,
        slew_m2_s=slew,
        energy_mj_kg=energy,
        area_mapping_error_m2=area_mapping_error,
        maximum_altitude_reference_difference_m=
            maximum_altitude_reference_difference,
        controller_drag_sampling_error_n=controller_drag_sampling_error,
        cached_drag_sampling_error_n=cached_drag_sampling_error,
        solve_count=state.solve_count,
        solver_status=state.last_solution.solver_status,
    )
end

# =============================================================================
# PRESENTATION PLOTS
# =============================================================================
function plot_histories(results, field, ylabel_text; limit=nothing, title_text="")
    styles = [:solid, :dash, :dot]
    colors = [:black, :blue, :green]
    plot_object = plot(xlabel="Elapsed atmospheric-pass time (s)",
        ylabel=ylabel_text, title=title_text, grid=true, legend=:best)
    for (index, result) in enumerate(results)
        plot!(plot_object, result.rollout.time_s, getproperty(result.rollout, field),
            label=result.case.label, color=colors[index], linestyle=styles[index], linewidth=2)
    end
    limit !== nothing && hline!(plot_object, [limit], label="Physical limit",
        color=:red, linestyle=:dash, linewidth=1.5)
    return plot_object
end

function save_plots(results)
    area_plot = plot_histories(results, :area_m2, "Commanded area (m²)";
        title_text="Commanded exposed area")
    hline!(area_plot, [BUS_AREA_M2], label="Minimum area",
        color=:red, linestyle=:dash)
    hline!(area_plot, [BUS_AREA_M2 + SOLAR_PANEL_AREA_M2], label="Maximum area",
        color=:red, linestyle=:dashdot)
    heat_rate_plot = plot_histories(results, :heat_rate_w_cm2,
        "Heat rate (W/cm²)"; limit=HEAT_RATE_LIMIT_W_CM2, title_text="Heat rate")
    heat_load_plot = plot_histories(results, :heat_load_j_cm2,
        "Heat load (J/cm²)"; limit=HEAT_LOAD_LIMIT_J_CM2, title_text="Integrated heat load")
    drag_plot = plot_histories(results, :drag_n, "Drag force (N)";
        limit=DRAG_FORCE_LIMIT_N, title_text="Drag force")
    slew_plot = plot_histories(results, :slew_m2_s, "Area rate (m²/s)";
        title_text="Actuator slew rate")
    hline!(slew_plot, [-AREA_SLEW_LIMIT_M2_S], label="Negative limit",
        color=:red, linestyle=:dash)
    hline!(slew_plot, [AREA_SLEW_LIMIT_M2_S], label="Positive limit",
        color=:red, linestyle=:dashdot)
    energy_plot = plot_histories(results, :energy_mj_kg,
        "Specific energy (MJ/kg)"; title_text="Maximum energy depletion")
    constraint_figure = plot(area_plot, heat_rate_plot, heat_load_plot, drag_plot,
        slew_plot, energy_plot; layout=(3, 2), size=(1300, 1050), dpi=180)

    prediction_error = plot(layout=(2, 2), size=(1250, 780), dpi=180)
    fields = ((2, :drag_n, "Drag prediction error (N)"),
        (3, :heat_rate_w_cm2, "Heat-rate prediction error (W/cm²)"),
        (4, :heat_load_j_cm2, "Heat-load prediction error (J/cm²)"),
        (5, :energy_mj_kg, "Energy prediction error (MJ/kg)"))
    for (panel, field, label) in fields
        for (index, result) in enumerate(results)
            rollout = result.rollout
            predicted = result.prediction
            prediction_values = field == :drag_n ? predicted[:, 2] :
                field == :heat_rate_w_cm2 ? predicted[:, 3] ./ 1.0e4 :
                field == :heat_load_j_cm2 ? result.predicted_heat_load_j_cm2 :
                predicted[:, 4] ./ 1.0e6
            nonlinear_at_nodes = [interpolate_mpc_history(
                rollout.absolute_time_s, getproperty(rollout, field), time)
                for time in result.problem.t]
            interior = 2:(length(result.problem.t) - 1)
            plot!(prediction_error[panel - 1],
                (result.problem.t .- first(result.problem.t))[interior],
                (nonlinear_at_nodes .- prediction_values)[interior],
                label=result.case.label, linewidth=2)
        end
        plot!(prediction_error[panel - 1], xlabel="Elapsed time (s)", ylabel=label,
            title=replace(label, " (" => "\n("), grid=true)
    end

    prefix = String(MPC_PLANET_NAME)
    for (name, figure) in ((prefix * "_med_constraint_cases", constraint_figure),
            (prefix * "_med_prediction_errors", prediction_error))
        savefig(figure, joinpath(OUTPUT_DIRECTORY, name * ".png"))
        WRITE_PDF && savefig(figure, joinpath(OUTPUT_DIRECTORY, name * ".pdf"))
    end
    return nothing
end

function save_callback_plots(callback_results)
    styles = [:solid, :dash, :dot]
    colors = [:black, :blue, :green]
    specifications = (
        (:area_m2, "Commanded area (m²)", "Callback commanded area"),
        (:alpha_rad, "Panel AOA (deg)", "Area-to-AOA callback output"),
        (:heat_rate_w_cm2, "Heat rate (W/cm²)", "Callback heat rate"),
        (:heat_load_j_cm2, "Heat load (J/cm²)", "Callback heat load"),
        (:drag_n, "Drag force (N)", "Callback drag force"),
        (:slew_m2_s, "Area rate (m²/s)", "Callback actuator slew"),
        (:energy_mj_kg, "Specific energy (MJ/kg)", "Cartesian truth energy"),
        (:altitude_m, "Altitude (km)", "Cartesian truth altitude"),
    )
    figure = plot(layout=(4, 2), size=(1300, 1350), dpi=180)
    for (panel, (field, label, title_text)) in enumerate(specifications)
        for (index, result) in enumerate(callback_results)
            values = getproperty(result, field)
            field === :alpha_rad && (values = rad2deg.(values))
            field === :altitude_m && (values = values ./ 1.0e3)
            plot!(figure[panel], result.time_s, values,
                label=result.case.label, color=colors[index],
                linestyle=styles[index], linewidth=2)
        end
        plot!(figure[panel], xlabel="Elapsed atmospheric-pass time (s)",
            ylabel=label, title=title_text, grid=true)
    end
    hline!(figure[1], [BUS_AREA_M2, BUS_AREA_M2 + SOLAR_PANEL_AREA_M2],
        label=["Minimum area" "Maximum area"], color=:red, linestyle=:dash)
    hline!(figure[3], [HEAT_RATE_LIMIT_W_CM2], label="Physical limit",
        color=:red, linestyle=:dash)
    hline!(figure[4], [HEAT_LOAD_LIMIT_J_CM2], label="Physical limit",
        color=:red, linestyle=:dash)
    hline!(figure[5], [DRAG_FORCE_LIMIT_N], label="Physical limit",
        color=:red, linestyle=:dash)
    hline!(figure[6], [-AREA_SLEW_LIMIT_M2_S, AREA_SLEW_LIMIT_M2_S],
        label=["Negative limit" "Positive limit"], color=:red, linestyle=:dash)
    for extension in (WRITE_PDF ? ("png", "pdf") : ("png",))
        savefig(figure, joinpath(
            OUTPUT_DIRECTORY,
            "$(MPC_PLANET_NAME)_med_spaceagora_callback_cases.$extension"))
    end
    return nothing
end

function main()
    mkpath(OUTPUT_DIRECTORY)
    density = DENSITY_FOR_MPC
    params = AerobrakingMPCParams(
        Re=PLANET.Rp_e, μ=PLANET.μ, J2=PLANET.J2, Ω=norm(PLANET.ω))
    position_0, velocity_0 = orbital_elements_to_cartesian(INITIAL_CONDITION, PLANET)
    config = base_config()
    reference_config = AerobrakingMPCReferenceConfig(
        h_cut_m=ATMOSPHERIC_CUTOFF_ALTITUDE_M,
        delta_s=KS_FICTITIOUS_TIME_STEP,
        max_coast_steps=MAX_COAST_STEPS,
        max_pass_steps=MAX_PASS_STEPS,
    )

    println("Phase 1/5: building the corrected hKS=-epsilon nonlinear KS reference")
    reference = build_reference_drag_pass(
        params, position_0, velocity_0;
        config=config,
        reference=reference_config,
        nominal_area_m2=BUS_AREA_M2 + SOLAR_PANEL_AREA_M2,
        density=density,
    )
    @assert all(diff(reference.time_s) .> 0.0)
    @assert minimum(reference.altitude_m) > 0.0

    println("Phase 2/5: constructing the common linearized KS MPC problem")
    problem = build_mpc_problem(
        reference, params, config; density=density, max_nodes=QP_MAX_NODES)
    @assert problem.N == size(reference.states, 1)
    @assert all(isfinite, problem.H)
    @assert all(isfinite, problem.Ybar)

    println("Phase 3/5: solving Cases I--III and running nonlinear verification")
    results = NamedTuple[]
    output_tables = DataFrame[]
    summaries = NamedTuple[]
    refinement_rows = NamedTuple[]
    for case in CASES
        case_config = apply_constraints(config, case.constraints)
        solve_time = @elapsed raw_solution = solve_mpc_qp(problem, case_config)
        raw_solution.ok || error("$(case.label) QP failed: $(raw_solution.solver_status)")
        # Propagate the optimizer output without post-solve clamping so the
        # nonlinear audit measures the QP command directly.
        area_plan = raw_solution.commanded_area_m2
        prediction = raw_solution.predicted_outputs
        rollout = propagate_ks_mpc_plan(
            reference, problem, area_plan, case_config, params;
            density=density, max_steps=MAX_PASS_STEPS)
        refined_rollout = propagate_ks_mpc_plan(
            reference, problem, area_plan, case_config, params;
            density=density,
            delta_s=0.5 * KS_FICTITIOUS_TIME_STEP,
            max_steps=MAX_PASS_STEPS)
        refined_max_slew = maximum(filter(isfinite,
            abs.(refined_rollout.slew_m2_s)))
        println("  $(case.label): status=$(raw_solution.solver_status), " *
            "solve=$(round(solve_time; digits=3)) s, " *
            "qdot=$(round(maximum(refined_rollout.heat_rate_w_cm2); digits=6)) W/cm^2, " *
            "Q=$(round(last(refined_rollout.heat_load_j_cm2); digits=6)) J/cm^2, " *
            "drag=$(round(maximum(refined_rollout.drag_n); digits=6)) N, " *
            "slew=$(round(refined_max_slew; digits=6)) m^2/s")
        @assert maximum(refined_rollout.heat_rate_w_cm2) <= HEAT_RATE_LIMIT_W_CM2
        @assert maximum(refined_rollout.drag_n) <= DRAG_FORCE_LIMIT_N
        if constraint_active(case.constraints, :heat_load)
            @assert last(refined_rollout.heat_load_j_cm2) <= HEAT_LOAD_LIMIT_J_CM2
        end
        if constraint_active(case.constraints, :slew)
            @assert refined_max_slew <= AREA_SLEW_LIMIT_M2_S
        end
        push!(refinement_rows,
            refinement_comparison(case.label, rollout, refined_rollout))
        predicted_heat_load = cumulative_mpc_heat_load(
            prediction[:, 3] ./ 1.0e4, problem.t)
        solution = (
            solver_status=raw_solution.solver_status,
            area_m2=area_plan,
        )
        result = (
            case=case,
            config=case_config,
            problem=problem,
            solution=solution,
            prediction=prediction,
            predicted_heat_load_j_cm2=predicted_heat_load,
            rollout=rollout,
        )
        push!(results, result)
        push!(summaries, case_summary(
            case, solution, solve_time, prediction, rollout, problem))
        predicted_slew = vcat(NaN, diff(area_plan) ./ diff(problem.t))
        push!(output_tables, output_dataframe(
            case.label, "linear_MPC_prediction",
            problem.t .- first(problem.t), prediction[:, 1], area_plan,
            prediction[:, 3] ./ 1.0e4, predicted_heat_load,
            prediction[:, 2], predicted_slew, prediction[:, 4] ./ 1.0e6))
        push!(output_tables, output_dataframe(
            case.label, "nonlinear_KS_verification",
            rollout.time_s, rollout.altitude_m, rollout.area_m2,
            rollout.heat_rate_w_cm2, rollout.heat_load_j_cm2,
            rollout.drag_n, rollout.slew_m2_s, rollout.energy_mj_kg))
    end

    summary = DataFrame(summaries)
    for row in eachrow(summary)
        @assert row.heat_rate_constraint_satisfied
        @assert row.heat_load_constraint_satisfied
        @assert row.drag_constraint_satisfied
        @assert row.slew_constraint_satisfied
    end

    callback_results = NamedTuple[]
    callback_summary = DataFrame()
    if RUN_SPACEAGORA_CALLBACK_CASES
        println("Phase 4/5: running Cases I--III through the SpaceAGORA control callback")
        callback_rows = NamedTuple[]
        callback_tables = DataFrame[]
        for (offline, case) in zip(results, CASES)
            callback = run_spaceagora_callback_case(case, reference)
            push!(callback_results, callback)
            raw_plan_difference = length(callback.controller.state.plan_area_m2) ==
                    length(offline.solution.area_m2) ?
                maximum(abs.(callback.controller.state.plan_area_m2 .-
                    offline.solution.area_m2)) : Inf
            finite_heat_rate = filter(isfinite, callback.heat_rate_w_cm2)
            finite_heat_load = filter(isfinite, callback.heat_load_j_cm2)
            finite_drag = filter(isfinite, callback.drag_n)
            finite_slew = filter(isfinite, abs.(callback.slew_m2_s))
            push!(callback_rows, (
                case=case.label,
                constraints=join(String.(constraint_names(case.constraints)), "+"),
                solver_status=String(callback.solver_status),
                solve_count=callback.solve_count,
                simulation_runtime_s=callback.runtime_s,
                minimum_commanded_area_m2=minimum(callback.area_m2),
                maximum_commanded_area_m2=maximum(callback.area_m2),
                maximum_heat_rate_w_cm2=maximum(finite_heat_rate),
                final_heat_load_j_cm2=last(finite_heat_load),
                maximum_drag_n=maximum(finite_drag),
                maximum_area_slew_m2_s=maximum(finite_slew),
                final_specific_energy_mj_kg=last(callback.energy_mj_kg),
                maximum_area_to_aoa_roundtrip_error_m2=
                    callback.area_mapping_error_m2,
                maximum_controller_drag_sampling_error_n=
                    callback.controller_drag_sampling_error_n,
                maximum_cached_dynamics_drag_sampling_error_n=
                    callback.cached_drag_sampling_error_n,
                maximum_geodetic_to_spherical_altitude_difference_m=
                    callback.maximum_altitude_reference_difference_m,
                raw_plan_difference_from_offline_qp_m2=raw_plan_difference,
            ))
            push!(callback_tables, DataFrame(
                case=fill(case.label, length(callback.time_s)),
                time_s=callback.time_s,
                altitude_km=callback.altitude_m ./ 1.0e3,
                geodetic_altitude_km=callback.geodetic_altitude_m ./ 1.0e3,
                commanded_area_m2=callback.area_m2,
                commanded_alpha_rad=callback.alpha_rad,
                commanded_alpha_deg=rad2deg.(callback.alpha_rad),
                heat_rate_w_cm2=callback.heat_rate_w_cm2,
                heat_load_j_cm2=callback.heat_load_j_cm2,
                accepted_state_heat_rate_w_cm2=callback.heat_rate_w_cm2,
                controller_sampled_heat_rate_w_cm2=
                    callback.controller_heat_rate_w_cm2,
                accepted_state_heat_load_j_cm2=callback.heat_load_j_cm2,
                controller_accumulated_heat_load_j_cm2=
                    callback.controller_heat_load_j_cm2,
                accepted_state_drag_n=callback.drag_n,
                controller_sampled_drag_n=callback.controller_drag_n,
                cached_internal_stage_drag_n=callback.cached_dynamics_drag_n,
                area_slew_m2_s=callback.slew_m2_s,
                specific_energy_mj_kg=callback.energy_mj_kg,
            ))
            println("  $(case.label) callback: solve_count=$(callback.solve_count), " *
                "qdot=$(round(maximum(finite_heat_rate); digits=6)) W/cm^2, " *
                "Q=$(round(last(finite_heat_load); digits=6)) J/cm^2, " *
                "drag=$(round(maximum(finite_drag); digits=6)) N, " *
                "slew=$(round(maximum(finite_slew); digits=6)) m^2/s, " *
                "area/AOA error=$(callback.area_mapping_error_m2) m^2")
        end
        callback_summary = DataFrame(callback_rows)
        CSV.write(joinpath(OUTPUT_DIRECTORY, "callback_summary.csv"), callback_summary)
        CSV.write(joinpath(OUTPUT_DIRECTORY, "callback_histories.csv"),
            reduce(vcat, callback_tables))
        save_callback_plots(callback_results)
    end

    println("Phase 5/5: writing audited results and presentation plots")
    CSV.write(joinpath(OUTPUT_DIRECTORY, "summary.csv"), summary)
    CSV.write(joinpath(OUTPUT_DIRECTORY, "step_refinement.csv"),
        DataFrame(refinement_rows))
    CSV.write(joinpath(OUTPUT_DIRECTORY, "histories.csv"), reduce(vcat, output_tables))
    metadata = DataFrame(
        parameter=[
            "planet", "planet_equatorial_radius_m", "planet_polar_radius_m",
            "planet_gravitational_parameter_m3_s2", "planet_j2",
            "periapsis_altitude_m", "apoapsis_radius_m",
            "inclination_deg", "raan_deg", "argument_of_periapsis_deg",
            "initial_true_anomaly_deg", "initial_position_m", "initial_velocity_m_s",
            "spacecraft_mass_kg", "bus_area_m2", "solar_panel_area_m2",
            "drag_coefficient", "heat_rate_limit_w_cm2", "heat_load_limit_j_cm2",
            "drag_force_limit_n", "area_slew_limit_m2_s", "density_selection",
            "density_model_type",
            "heat_rate_qp_backoff_w_cm2", "heat_load_qp_backoff_j_cm2",
            "drag_force_qp_backoff_n", "area_slew_qp_backoff_m2_s",
            "osqp_eps_abs", "osqp_eps_rel", "osqp_max_iter",
            "run_spaceagora_callback_cases", "control_sample_time_s",
            "cartesian_relative_tolerance", "cartesian_absolute_tolerance",
            "cartesian_max_step_s", "ks_energy_convention",
            "ks_fictitious_time_step", "reference_nodes",
            "atmospheric_pass_duration_s", "minimum_reference_altitude_m",
        ],
        value=string.([
            nameof(typeof(PLANET)), PLANET.Rp_e, PLANET.Rp_p,
            PLANET.μ, PLANET.J2,
            PERIAPSIS_ALTITUDE_M, APOAPSIS_RADIUS_M,
            INCLINATION_DEG, RAAN_DEG, ARGUMENT_OF_PERIAPSIS_DEG,
            INITIAL_TRUE_ANOMALY_DEG, repr(position_0), repr(velocity_0),
            SPACECRAFT_MASS_KG, BUS_AREA_M2, SOLAR_PANEL_AREA_M2,
            DRAG_COEFFICIENT, HEAT_RATE_LIMIT_W_CM2, HEAT_LOAD_LIMIT_J_CM2,
            DRAG_FORCE_LIMIT_N, AREA_SLEW_LIMIT_M2_S, DENSITY_SELECTION,
            nameof(typeof(DENSITY_MODEL)),
            HEAT_RATE_QP_BACKOFF_W_CM2, HEAT_LOAD_QP_BACKOFF_J_CM2,
            DRAG_FORCE_QP_BACKOFF_N, AREA_SLEW_QP_BACKOFF_M2_S,
            OSQP_EPS_ABS, OSQP_EPS_REL, OSQP_MAX_ITER,
            RUN_SPACEAGORA_CALLBACK_CASES, CONTROL_SAMPLE_TIME_S,
            CARTESIAN_RELATIVE_TOLERANCE, CARTESIAN_ABSOLUTE_TOLERANCE,
            CARTESIAN_MAX_STEP_S, "h_KS=-specific_energy",
            KS_FICTITIOUS_TIME_STEP, problem.N,
            last(reference.time_s) - first(reference.time_s),
            minimum(reference.altitude_m),
        ]),
    )
    CSV.write(joinpath(OUTPUT_DIRECTORY, "metadata.csv"), metadata)
    save_plots(results)

    show(stdout, MIME("text/plain"), summary; allcols=true)
    println("\n\nResults: $(abspath(OUTPUT_DIRECTORY))")
    return results, summary, callback_results, callback_summary
end

main()
