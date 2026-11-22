#=
Single energy-targeting MPC pass for Earth or Mars.
Set SPACEAGORA_MPC_PLANET to earth or mars.

The MPC runs in energy-targeting mode (TargetEnergyMode). The QP is solved once,
at the inbound interface crossing, about the pass flown with the arrays fully
deployed. Its cost drives the predicted end-of-pass specific energy to the
target energy. The heat-rate, heat-load and drag-force limits are soft
(penalized slacks); the area bounds and the area-rate limit are hard.

Optional settings:
  SPACEAGORA_MPC_PERIAPSIS_KM         periapsis altitude of the initial orbit, km
                                      (default: Earth 110, Mars 100)
  SPACEAGORA_MPC_TARGET_ENERGY_MJ_KG  target end-of-pass specific energy, MJ/kg
                                      (default: Earth -6.41; Mars the midpoint of
                                      the reachable range)

The reachable range is bounded by the end-of-pass energies of the passes flown
at constant area with the arrays fully deployed and with the bus alone. The
script stops if the target lies outside it.

Outputs, in output/mpc_<planet>_targeting_pass:
  histories.csv  time histories with the source names
                 linear_MPC_prediction      QP prediction at the nodes
                 nonlinear_KS_verification  nonlinear KS propagation of the QP plan
                 fully_deployed             pass at maximum area
                 bus_only                   pass at minimum area
  summary.csv    terminal energies and errors, peaks and their fraction of the
                 limits, prediction errors, slacks, nodes and solve time
  metadata.csv   all settings
  <planet>_targeting_pass.png/.pdf
=#

include(joinpath(@__DIR__, "common.jl"))

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

using CSV
using DataFrames
using LinearAlgebra
using Plots
using Printf

# =============================================================================
# PASS SETTINGS
# =============================================================================
const MPC_PLANET_NAME = Symbol(lowercase(get(ENV, "SPACEAGORA_MPC_PLANET", "earth")))
MPC_PLANET_NAME in (:earth, :mars) || throw(ArgumentError(
    "SPACEAGORA_MPC_PLANET must be earth or mars; got $(repr(MPC_PLANET_NAME))."))

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
        target_energy_mj_kg=-6.41,
    )
else
    (
        planet=Mars(),
        initial_time=InitialTime(year=2020, month=1, day=1),
        periapsis_altitude_m=100.0e3,
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
        target_energy_mj_kg=nothing,  # midpoint of the reachable range
    )
end

_env_float(name) = haskey(ENV, name) && !isempty(ENV[name]) ? parse(Float64, ENV[name]) : nothing

# Initial orbit: apoapsis radius and periapsis altitude above the equatorial
# radius, starting at apoapsis.
const PLANET = PLANET_SETTINGS.planet
const INITIAL_TIME = PLANET_SETTINGS.initial_time
const PERIAPSIS_ALTITUDE_M = something(
    _env_float("SPACEAGORA_MPC_PERIAPSIS_KM"), PLANET_SETTINGS.periapsis_altitude_m / 1.0e3) * 1.0e3
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
const REQUESTED_TARGET_ENERGY_MJ_KG = something(
    _env_float("SPACEAGORA_MPC_TARGET_ENERGY_MJ_KG"),
    PLANET_SETTINGS.target_energy_mj_kg, Some(nothing))

# Spacecraft: the bus area is the minimum exposed area, and the bus plus both
# solar arrays the maximum.
const SPACECRAFT_MASS_KG = PLANET_SETTINGS.spacecraft_mass_kg
const BUS_AREA_M2 = PLANET_SETTINGS.bus_area_m2
const SOLAR_PANEL_AREA_M2 = PLANET_SETTINGS.solar_panel_area_m2
const MAXIMUM_AREA_M2 = BUS_AREA_M2 + SOLAR_PANEL_AREA_M2
const DRAG_COEFFICIENT = 2.2

# Physical limits.
const HEAT_RATE_LIMIT_W_CM2 = PLANET_SETTINGS.heat_rate_limit_w_cm2
const HEAT_LOAD_LIMIT_J_CM2 = PLANET_SETTINGS.heat_load_limit_j_cm2
const DRAG_FORCE_LIMIT_N = PLANET_SETTINGS.drag_force_limit_n
const AREA_SLEW_LIMIT_M2_S = 0.20

# Margins: the QP limit is the physical limit minus the margin.
const HEAT_RATE_QP_BACKOFF_W_CM2 = 1.0e-3
const HEAT_LOAD_QP_BACKOFF_J_CM2 = 5.0e-2
const DRAG_FORCE_QP_BACKOFF_N = 2.0e-2
const AREA_SLEW_QP_BACKOFF_M2_S = 2.0e-3

# Active limits. With the area-rate limit, the first area change of the plan
# is zero and each later change is bounded by the QP rate limit times the node
# spacing.
const CONSTRAINTS = mpc_constraints(:heat_rate, :heat_load, :drag, :slew)

# QP weights and solver settings. The cost is
#   TARGET_ENERGY_WEIGHT * (E_N - E*)^2      with E in J/kg
# + SLACK_WEIGHT * sum(sigma^2)              sigma: violation divided by the QP limit
# + AREA_WEIGHT * sum(dA^2) + AREA_SLEW_WEIGHT * sum((dA_{k+1} - dA_k)^2).
# The QP adds 1e-10 to the diagonal of its Hessian to keep it positive definite.
const AREA_WEIGHT = 0.0
const AREA_SLEW_WEIGHT = 0.0
const SLACK_WEIGHT = 1.0e5
const TARGET_ENERGY_WEIGHT = 5.0e-10
const OSQP_EPS_ABS = 1.0e-7
const OSQP_EPS_REL = 1.0e-7
const OSQP_MAX_ITER = 100_000

# Prediction: nominal pass at maximum area, KS fictitious-time step in s/m,
# all nodes of the reference pass in the QP.
const ATMOSPHERIC_CUTOFF_ALTITUDE_M = PLANET_SETTINGS.interface_altitude_m
const KS_FICTITIOUS_TIME_STEP = 1.7e-7
const MAX_COAST_STEPS = 2_000_000
const MAX_PASS_STEPS = 20_000
const QP_MAX_NODES = nothing

const DENSITY_SELECTION = :planet_polyfit
const OUTPUT_DIRECTORY = joinpath(
    @__DIR__, "..", "output", "mpc_$(MPC_PLANET_NAME)_targeting_pass")
const WRITE_PDF = true

# =============================================================================
# ATMOSPHERE
# =============================================================================
# One polynomial fit of the planet's density for the prediction and for the
# nonlinear propagation.
const DENSITY_MODEL = PolynomialFitAtmosphereModel(PLANET)
const DENSITY_CONTEXT = (
    args=(environment_model=(planet=PLANET, density_model=DENSITY_MODEL),),)
const DENSITY_FOR_MPC = density_function_from_spaceagora(DENSITY_CONTEXT)

# =============================================================================
# CONFIGURATION AND PROPAGATION
# =============================================================================
function targeting_config(target_energy_mj_kg)
    return AerobrakingMPCConfig(
        mode=TargetEnergyMode(),
        bus_reference_area_m2=BUS_AREA_M2,
        controllable_area_m2=SOLAR_PANEL_AREA_M2,
        mass_kg=SPACECRAFT_MASS_KG,
        drag_coefficient=DRAG_COEFFICIENT,
        qdot_max_w_cm2=HEAT_RATE_LIMIT_W_CM2 - HEAT_RATE_QP_BACKOFF_W_CM2,
        heat_load_max_j_cm2=HEAT_LOAD_LIMIT_J_CM2 - HEAT_LOAD_QP_BACKOFF_J_CM2,
        drag_max_n=DRAG_FORCE_LIMIT_N - DRAG_FORCE_QP_BACKOFF_N,
        area_slew_max_m2_s=AREA_SLEW_LIMIT_M2_S - AREA_SLEW_QP_BACKOFF_M2_S,
        use_constraints=true,
        use_slew_constraint=constraint_active(CONSTRAINTS, :slew),
        use_qdot_constraint=constraint_active(CONSTRAINTS, :heat_rate),
        use_heat_load_constraint=constraint_active(CONSTRAINTS, :heat_load),
        use_drag_constraint=constraint_active(CONSTRAINTS, :drag),
        target_energy_mj_kg=target_energy_mj_kg,
        area_weight=AREA_WEIGHT,
        area_slew_weight=AREA_SLEW_WEIGHT,
        slack_weight=SLACK_WEIGHT,
        target_energy_weight=TARGET_ENERGY_WEIGHT,
        max_depletion_energy_weight=0.0,
        osqp_eps_abs=OSQP_EPS_ABS,
        osqp_eps_rel=OSQP_EPS_REL,
        osqp_max_iter=OSQP_MAX_ITER,
    )
end

reference_config() = AerobrakingMPCReferenceConfig(
    h_cut_m=ATMOSPHERIC_CUTOFF_ALTITUDE_M,
    delta_s=KS_FICTITIOUS_TIME_STEP,
    max_coast_steps=MAX_COAST_STEPS,
    max_pass_steps=MAX_PASS_STEPS,
)

# Reference pass, QP problem and nonlinear propagation of a pass flown at a
# constant area.
function constant_area_pass(params, position_0, velocity_0, config, area_m2)
    reference = build_reference_drag_pass(
        params, position_0, velocity_0;
        config=config, reference=reference_config(),
        nominal_area_m2=area_m2, density=DENSITY_FOR_MPC)
    problem = build_mpc_problem(reference, params, config; density=DENSITY_FOR_MPC)
    rollout = propagate_ks_mpc_plan(
        reference, problem, fill(area_m2, problem.N), config, params;
        density=DENSITY_FOR_MPC, max_steps=MAX_PASS_STEPS)
    return reference, problem, rollout
end

function output_dataframe(source, time_s, altitude_m, area_m2,
    heat_rate_w_cm2, heat_load_j_cm2, drag_n, slew_m2_s, energy_mj_kg)
    return DataFrame(
        planet=fill(String(MPC_PLANET_NAME), length(time_s)),
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

rollout_dataframe(source, r) = output_dataframe(source, r.time_s, r.altitude_m,
    r.area_m2, r.heat_rate_w_cm2, r.heat_load_j_cm2, r.drag_n, r.slew_m2_s,
    r.energy_mj_kg)

function plot_pass(rollout, full, target_energy_mj_kg)
    gray = RGB(0.6, 0.6, 0.6)
    tf, t = full.time_s, rollout.time_s
    panels = [
        (full.area_m2, rollout.area_m2, "Exposed area, m^2", [BUS_AREA_M2, MAXIMUM_AREA_M2]),
        (full.slew_m2_s, rollout.slew_m2_s, "Area rate, m^2/s", nothing),
        (full.heat_rate_w_cm2, rollout.heat_rate_w_cm2, "Heat rate, W/cm^2", [HEAT_RATE_LIMIT_W_CM2]),
        (full.heat_load_j_cm2, rollout.heat_load_j_cm2, "Heat load, J/cm^2", [HEAT_LOAD_LIMIT_J_CM2]),
        (full.drag_n, rollout.drag_n, "Drag force, N", [DRAG_FORCE_LIMIT_N]),
        (full.energy_mj_kg, rollout.energy_mj_kg, "Specific energy, MJ/kg", [target_energy_mj_kg]),
    ]
    plots = map(enumerate(panels)) do (k, (yf, y, label, levels))
        p = plot(tf, yf; color=gray, lw=1.5, label=k == 1 ? "Arrays fully deployed" : "",
            xlabel=k >= 5 ? "Time, s" : "", ylabel=label, legend=k == 1 ? :bottomleft : false,
            guidefontsize=9, tickfontsize=8, legendfontsize=8)
        plot!(p, t, y; color=:steelblue, lw=1.8, label=k == 1 ? "MPC, energy targeting" : "")
        levels === nothing || hline!(p, levels; color=:black, ls=:dot,
            label=k == 1 ? "Area bounds, limits, target" : "")
        p
    end
    figure = plot(plots...; layout=(3, 2), size=(900, 900), margin=4Plots.mm)
    for extension in (WRITE_PDF ? ("png", "pdf") : ("png",))
        savefig(figure, joinpath(OUTPUT_DIRECTORY, "$(MPC_PLANET_NAME)_targeting_pass.$extension"))
    end
    return nothing
end

function main()
    mkpath(OUTPUT_DIRECTORY)
    density = DENSITY_FOR_MPC
    params = KSPropagationParams(
        Re=PLANET.Rp_e, μ=PLANET.μ, J2=PLANET.J2, Ω=norm(PLANET.ω))
    position_0, velocity_0 = orbital_elements_to_cartesian(INITIAL_CONDITION, PLANET)

    println("Phase 1/4: reachable end-of-pass energies (fully deployed and bus alone)")
    probe_config = targeting_config(0.0)
    full_reference, _, full_rollout = constant_area_pass(
        params, position_0, velocity_0, probe_config, MAXIMUM_AREA_M2)
    _, _, bus_rollout = constant_area_pass(
        params, position_0, velocity_0, probe_config, BUS_AREA_M2)
    energy_full = last(full_rollout.energy_mj_kg)
    energy_bus = last(bus_rollout.energy_mj_kg)
    target_energy = something(REQUESTED_TARGET_ENERGY_MJ_KG, 0.5 * (energy_full + energy_bus))
    @printf("  fully deployed %.9f MJ/kg, bus alone %.9f MJ/kg, target %.9f MJ/kg\n",
        energy_full, energy_bus, target_energy)
    energy_full < target_energy < energy_bus || error(
        "Target $(target_energy) MJ/kg is outside the reachable range " *
        "($(energy_full), $(energy_bus)) MJ/kg of this pass.")

    println("Phase 2/4: linearizing about the fully deployed pass")
    config = targeting_config(target_energy)
    reference = full_reference
    problem = build_mpc_problem(
        reference, params, config; density=density, max_nodes=QP_MAX_NODES)
    @assert problem.N == size(reference.states, 1)
    @assert all(isfinite, problem.H)
    @assert all(isfinite, problem.Ybar)

    println("Phase 3/4: solving the targeting QP and propagating its plan")
    # The QP is solved twice and the second solve is timed, so the solve time
    # excludes Julia compilation.
    solve_mpc_qp(problem, config)
    solve_time = @elapsed solution = solve_mpc_qp(problem, config)
    solution.ok || error("Targeting QP failed: $(solution.solver_status)")
    area_plan = solution.commanded_area_m2
    prediction = solution.predicted_outputs
    rollout = propagate_ks_mpc_plan(
        reference, problem, area_plan, config, params;
        density=density, max_steps=MAX_PASS_STEPS)
    refined_rollout = propagate_ks_mpc_plan(
        reference, problem, area_plan, config, params;
        density=density, delta_s=0.5 * KS_FICTITIOUS_TIME_STEP, max_steps=MAX_PASS_STEPS)

    # Prediction at the nodes: outputs are [altitude m, drag N, heat rate W/m^2,
    # specific energy MJ/kg].
    predicted_time = problem.t .- first(problem.t)
    predicted_heat_rate = prediction[:, 3] ./ 1.0e4
    predicted_heat_load = cumulative_mpc_heat_load(predicted_heat_rate, problem.t)
    predicted_slew = vcat(NaN, diff(area_plan) ./ diff(problem.t))
    sample(field, r) = [interpolate_mpc_plan(r.time_s, getproperty(r, field), time)
        for time in predicted_time]
    node_spacing = diff(problem.t)
    N = problem.N
    # Slacks: N heat-rate entries, one heat-load entry, N drag entries.
    slacks = solution.slacks
    maxslack(range) = isempty(range) ? NaN : maximum(slacks[range])
    terminal_error(energy) = energy - target_energy
    peak_slew = maximum(filter(isfinite, abs.(rollout.slew_m2_s)))

    predicted_terminal = solution.terminal_energy
    nonlinear_terminal = last(rollout.energy_mj_kg)
    summary = DataFrame(
        quantity=[
            "solver_status", "solve_time_s", "qp_solves", "nodes",
            "node_spacing_min_s", "node_spacing_max_s",
            "minimum_altitude_km", "target_energy_mj_kg",
            "fully_deployed_energy_mj_kg", "bus_only_energy_mj_kg",
            "predicted_terminal_energy_mj_kg", "nonlinear_terminal_energy_mj_kg",
            "predicted_terminal_error_j_kg", "nonlinear_terminal_error_j_kg",
            "predicted_terminal_error_percent", "nonlinear_terminal_error_percent",
            "prediction_vs_nonlinear_terminal_energy_j_kg",
            "nonlinear_peak_heat_rate_w_cm2", "nonlinear_peak_heat_rate_percent_of_limit",
            "nonlinear_heat_load_j_cm2", "nonlinear_heat_load_percent_of_limit",
            "nonlinear_peak_drag_n", "nonlinear_peak_drag_percent_of_limit",
            "nonlinear_peak_area_rate_m2_s", "nonlinear_peak_area_rate_percent_of_limit",
            "minimum_area_m2", "maximum_area_m2",
            "max_prediction_error_heat_rate_percent_of_limit",
            "max_prediction_error_drag_percent_of_limit",
            "max_heat_rate_slack", "heat_load_slack", "max_drag_slack",
            "half_step_terminal_energy_difference_j_kg",
        ],
        value=Any[
            String(solution.solver_status), solve_time, 1, N,
            minimum(node_spacing), maximum(node_spacing),
            minimum(rollout.altitude_m) / 1.0e3, target_energy,
            energy_full, energy_bus,
            predicted_terminal, nonlinear_terminal,
            terminal_error(predicted_terminal) * 1.0e6, terminal_error(nonlinear_terminal) * 1.0e6,
            100 * abs(terminal_error(predicted_terminal) / target_energy),
            100 * abs(terminal_error(nonlinear_terminal) / target_energy),
            (predicted_terminal - nonlinear_terminal) * 1.0e6,
            maximum(rollout.heat_rate_w_cm2), 100 * maximum(rollout.heat_rate_w_cm2) / HEAT_RATE_LIMIT_W_CM2,
            last(rollout.heat_load_j_cm2), 100 * last(rollout.heat_load_j_cm2) / HEAT_LOAD_LIMIT_J_CM2,
            maximum(rollout.drag_n), 100 * maximum(rollout.drag_n) / DRAG_FORCE_LIMIT_N,
            peak_slew, 100 * peak_slew / AREA_SLEW_LIMIT_M2_S,
            minimum(area_plan), maximum(area_plan),
            100 * maximum(abs.(predicted_heat_rate .- sample(:heat_rate_w_cm2, rollout))) / HEAT_RATE_LIMIT_W_CM2,
            100 * maximum(abs.(prediction[:, 2] .- sample(:drag_n, rollout))) / DRAG_FORCE_LIMIT_N,
            maxslack(1:N), N + 1 <= length(slacks) ? slacks[N + 1] : NaN, maxslack(N + 2:length(slacks)),
            (nonlinear_terminal - last(refined_rollout.energy_mj_kg)) * 1.0e6,
        ],
    )

    println("Phase 4/4: writing results")
    histories = vcat(
        output_dataframe("linear_MPC_prediction", predicted_time, prediction[:, 1], area_plan,
            predicted_heat_rate, predicted_heat_load, prediction[:, 2], predicted_slew,
            prediction[:, 4]),
        rollout_dataframe("nonlinear_KS_verification", rollout),
        rollout_dataframe("fully_deployed", full_rollout),
        rollout_dataframe("bus_only", bus_rollout),
    )
    CSV.write(joinpath(OUTPUT_DIRECTORY, "histories.csv"), histories)
    CSV.write(joinpath(OUTPUT_DIRECTORY, "summary.csv"), summary)
    metadata = DataFrame(
        parameter=[
            "planet", "density_selection", "density_model_type",
            "periapsis_altitude_m", "apoapsis_radius_m", "inclination_deg", "raan_deg",
            "argument_of_periapsis_deg", "initial_true_anomaly_deg",
            "initial_position_m", "initial_velocity_m_s",
            "spacecraft_mass_kg", "bus_area_m2", "solar_panel_area_m2", "drag_coefficient",
            "heat_rate_limit_w_cm2", "heat_load_limit_j_cm2", "drag_force_limit_n",
            "area_slew_limit_m2_s", "heat_rate_qp_backoff_w_cm2", "heat_load_qp_backoff_j_cm2",
            "drag_force_qp_backoff_n", "area_slew_qp_backoff_m2_s", "constraints",
            "area_weight", "area_slew_weight", "slack_weight", "target_energy_weight",
            "osqp_eps_abs", "osqp_eps_rel", "osqp_max_iter", "target_energy_mj_kg",
            "target_energy_source", "nominal_area_m2", "interface_altitude_m",
            "ks_fictitious_time_step_s_per_m", "qp_max_nodes",
        ],
        value=Any[
            String(MPC_PLANET_NAME), DENSITY_SELECTION, nameof(typeof(DENSITY_MODEL)),
            PERIAPSIS_ALTITUDE_M, APOAPSIS_RADIUS_M, INCLINATION_DEG, RAAN_DEG,
            ARGUMENT_OF_PERIAPSIS_DEG, INITIAL_TRUE_ANOMALY_DEG,
            repr(position_0), repr(velocity_0),
            SPACECRAFT_MASS_KG, BUS_AREA_M2, SOLAR_PANEL_AREA_M2, DRAG_COEFFICIENT,
            HEAT_RATE_LIMIT_W_CM2, HEAT_LOAD_LIMIT_J_CM2, DRAG_FORCE_LIMIT_N,
            AREA_SLEW_LIMIT_M2_S, HEAT_RATE_QP_BACKOFF_W_CM2, HEAT_LOAD_QP_BACKOFF_J_CM2,
            DRAG_FORCE_QP_BACKOFF_N, AREA_SLEW_QP_BACKOFF_M2_S,
            join(String.(constraint_names(CONSTRAINTS)), "+"),
            AREA_WEIGHT, AREA_SLEW_WEIGHT, SLACK_WEIGHT, TARGET_ENERGY_WEIGHT,
            OSQP_EPS_ABS, OSQP_EPS_REL, OSQP_MAX_ITER, target_energy,
            REQUESTED_TARGET_ENERGY_MJ_KG === nothing ? "midpoint_of_reachable_range" : "setting",
            MAXIMUM_AREA_M2, ATMOSPHERIC_CUTOFF_ALTITUDE_M,
            KS_FICTITIOUS_TIME_STEP, something(QP_MAX_NODES, "all"),
        ],
    )
    CSV.write(joinpath(OUTPUT_DIRECTORY, "metadata.csv"), metadata)
    plot_pass(rollout, full_rollout, target_energy)

    @printf("  status=%s, solve=%.3f s, nodes=%d\n", solution.solver_status, solve_time, N)
    @printf("  terminal energy: predicted %.9f, nonlinear %.9f, target %.9f MJ/kg\n",
        predicted_terminal, nonlinear_terminal, target_energy)
    @printf("  peaks (nonlinear): heat rate %.4f W/cm^2, heat load %.4f J/cm^2, drag %.4f N\n",
        maximum(rollout.heat_rate_w_cm2), last(rollout.heat_load_j_cm2), maximum(rollout.drag_n))
    println("\nResults: $(abspath(OUTPUT_DIRECTORY))")
    return nothing
end

main()
