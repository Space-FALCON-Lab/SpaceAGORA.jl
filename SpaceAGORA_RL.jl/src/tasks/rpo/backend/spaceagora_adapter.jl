function _spaceagora_rpo_modules()
    spaceagora = _load_spaceagora!(load_gramsuite=false)
    simulation_model = getproperty(spaceagora, :SimulationModel)
    return (
        guidance=getproperty(simulation_model, :GuidanceHooks),
        navigation=getproperty(simulation_model, :NavigationHooks),
        control=getproperty(simulation_model, :ControlHooks),
    )
end

"""Build a default Gateway/CubeSat scenario using SpaceAGORA's station asset."""
function build_rpo_hypr_rl_scenario(;
    start_rtn=[-8.0, -4.0, 2.0],
    goal_rtn=[5.0, 0.0, 0.0],
    station_asset::Symbol=:gateway,
    station_points::Int=2_000,
    station_seed::Int=740,
    station_keepout_radius_m::Real=0.25,
    initial_attitude_rtn_to_body=[0.0, 0.0, 0.0, 1.0],
    final_attitude_rtn_to_body=[0.0, 0.0, 0.0, 1.0],
    pso_config=nothing,
    tracking_settings=nothing,
    rrt_settings=nothing,
)
    spaceagora = _load_spaceagora!(load_gramsuite=false)
    modules = _spaceagora_rpo_modules()
    points = Base.invokelatest(
        getproperty(spaceagora, :load_rpo_station_cad_pointcloud), station_asset;
        n_points=station_points, rng=MersenneTwister(station_seed),
    )
    station = Base.invokelatest(
        getproperty(modules.navigation, :RPOStationGeometry), points;
        keepout_radius_m=Float64(station_keepout_radius_m),
        name=String(station_asset),
    )
    geometry = Base.invokelatest(
        getproperty(modules.navigation, :RPOReferenceGeometry), station,
    )
    return RPOHyPRRLScenario(
        start_rtn=Vector{Float64}(start_rtn),
        goal_rtn=Vector{Float64}(goal_rtn),
        geometry=geometry,
        pso_config=pso_config,
        tracking_settings=tracking_settings,
        rrt_settings=rrt_settings,
        initial_attitude_rtn_to_body=Vector{Float64}(initial_attitude_rtn_to_body),
        final_attitude_rtn_to_body=Vector{Float64}(final_attitude_rtn_to_body),
    )
end

function _rpo_spaceagora_settings(scenario::RPOHyPRRLScenario,
                                  config::RPOHyPRRLConfig)
    modules = _spaceagora_rpo_modules()
    pso_config = scenario.pso_config === nothing ?
        Base.invokelatest(
            getproperty(modules.guidance, :rpo_740_mpc_final_pso_config);
            safe_distance_m=config.safe_distance_m,
        ) : scenario.pso_config
    pso_config = Base.invokelatest(
        getproperty(modules.guidance, :rpo_pso_config), pso_config;
        curve_type=:bezier, safe_distance_m=config.safe_distance_m,
    )
    tracking = scenario.tracking_settings === nothing ?
        Base.invokelatest(getproperty(modules.guidance, :RPOLQMPCTrackingSettings)) :
        scenario.tracking_settings
    rrt = scenario.rrt_settings === nothing ?
        Base.invokelatest(getproperty(modules.guidance, :RPORRTConnectSettings)) :
        scenario.rrt_settings
    return modules, pso_config, tracking, rrt
end

function _limit_seed_waypoints(path::AbstractMatrix{<:Real}, maximum_internal::Int)
    n = size(path, 2)
    n <= maximum_internal + 2 && return Matrix{Float64}(path)
    selected = unique(round.(Int, range(1, n; length=maximum_internal + 2)))
    return Matrix{Float64}(path[:, selected])
end

"""Construct the feasible warm start edited by HyPR-RL."""
function rpo_hypr_rl_seed_path(scenario::RPOHyPRRLScenario,
                               config::RPOHyPRRLConfig;
                               rng::AbstractRNG=Random.default_rng())
    scenario.seed_path_rtn === nothing ||
        return _limit_seed_waypoints(scenario.seed_path_rtn,
                                     config.max_translation_waypoints)
    modules, pso_config, _, rrt = _rpo_spaceagora_settings(scenario, config)
    plan = Base.invokelatest(
        getproperty(modules.guidance, :rpo_rrt_connect_bezier_plan_path),
        scenario.start_rtn, scenario.goal_rtn, scenario.geometry, pso_config;
        safe_distance_m=config.safe_distance_m, settings=rrt, rng=rng,
    )
    return _limit_seed_waypoints(plan.path, config.max_translation_waypoints)
end

function _allocate_thrusters(desired_force_body::AbstractVector{<:Real},
                             config::RPOHyPRRLConfig)
    directions = config.thruster_directions_body
    gram = directions * transpose(directions) + 1.0e-12I
    thrust = max.(0.0, transpose(directions) * (gram \ Vector{Float64}(desired_force_body)))
    thrust = min.(thrust, config.thruster_max_thrust_n)
    step = 0.9 / max(opnorm(directions)^2, eps(Float64))
    for _ in 1:12
        residual = directions * thrust - desired_force_body
        thrust .= clamp.(thrust .- step .* (transpose(directions) * residual),
                         0.0, config.thruster_max_thrust_n)
    end
    return thrust
end

function _thruster_translation_step(command_rtn, attitude_rtn_to_body, dt_s,
                                    config::RPOHyPRRLConfig)
    dt = Float64(dt_s)
    rotation_rtn_to_body = _quat_rotation_matrix(attitude_rtn_to_body)
    desired_acceleration = Vector{Float64}(command_rtn)
    desired_force_body = config.mass_kg .* rotation_rtn_to_body * desired_acceleration
    commanded_thrust = _allocate_thrusters(desired_force_body, config)
    impulse = zeros(6)
    for thruster_index in 1:6
        maximum_thrust = config.thruster_max_thrust_n[thruster_index]
        if maximum_thrust <= 0.0 || commanded_thrust[thruster_index] <= 0.0
            continue
        end
        on_time = commanded_thrust[thruster_index] / maximum_thrust * dt
        on_time = clamp(max(on_time, config.thruster_min_firing_time_s), 0.0, dt)
        impulse[thruster_index] = maximum_thrust * on_time
    end
    average_thrust = impulse ./ dt
    realized_force_body = config.thruster_directions_body * average_thrust
    realized_acceleration = transpose(rotation_rtn_to_body) * realized_force_body /
                            config.mass_kg
    propellant = sum(
        impulse[index] / (config.thruster_isp_s[index] * config.g0_mps2)
        for index in 1:6
    )
    error = norm(realized_acceleration - desired_acceleration)
    return (
        realized_acceleration_rtn=realized_acceleration,
        propellant_used_kg=propellant,
        thruster_impulse_ns=impulse,
        allocation_error_mps2=error,
        saturated=error > 1.0e-5,
    )
end

function _account_thruster_and_attitude_history(u_hist::AbstractMatrix{<:Real},
                                                q_reference::AbstractMatrix{<:Real},
                                                dt_s::Real,
                                                scenario::RPOHyPRRLScenario,
                                                config::RPOHyPRRLConfig)
    n_steps = size(u_hist, 2)
    impulse = zeros(6)
    propellant = 0.0
    allocation_error = 0.0
    saturated_steps = 0
    wheel_energy = 0.0
    wheel_momentum = zeros(3)
    wheel_peak = 0.0
    attitude = _quat_normalize(scenario.initial_attitude_rtn_to_body)
    angular_rate = zeros(3)
    dt = Float64(dt_s)
    for step_index in 1:n_steps
        desired_attitude = view(q_reference, :, min(step_index, size(q_reference, 2)))
        attitude, angular_rate, torque = _integrate_attitude(
            attitude, angular_rate, desired_attitude, dt, config,
        )
        wheel_momentum .-= torque .* dt
        wheel_peak = max(wheel_peak, maximum(abs, wheel_momentum))
        wheel_energy += abs(dot(torque, angular_rate)) * dt

        desired_acceleration_rtn = Vector{Float64}(view(u_hist, :, step_index))
        actuation = _thruster_translation_step(
            desired_acceleration_rtn, attitude, dt, config,
        )
        impulse .+= actuation.thruster_impulse_ns
        propellant += actuation.propellant_used_kg
        allocation_error += actuation.allocation_error_mps2 * dt
        saturated_steps += actuation.saturated ? 1 : 0
    end
    return (
        propellant_used_kg=propellant,
        allocation_error_impulse_mps=allocation_error,
        saturation_fraction=n_steps > 0 ? saturated_steps / n_steps : 0.0,
        thruster_impulse_ns=impulse,
        wheel_energy_j=wheel_energy,
        wheel_peak_momentum_nms=wheel_peak,
        final_attitude_rtn_to_body=attitude,
        final_angular_rate_radps=angular_rate,
    )
end

function _coupled_rpo_tracking(modules, t_ref, r_ref, v_ref, goal_rtn,
                               attitude_progress, attitude_quaternions,
                               tracking, scenario::RPOHyPRRLScenario,
                               config::RPOHyPRRLConfig)
    Q = Diagonal(vcat(fill(tracking.q_pos, 3), fill(tracking.q_vel, 3)))
    R = Diagonal(fill(tracking.r_accel, 3))
    Qf = Diagonal(vcat(fill(tracking.qf_pos, 3), fill(tracking.qf_vel, 3)))
    u_max = Vector{Float64}(tracking.u_max_mps2)
    controller = Base.invokelatest(
        getproperty(modules.control, :init_rpo_lqmpc),
        tracking.mean_motion_radps, tracking.dt_s, Q, R, Qf, tracking.horizon;
        u_min=-u_max, u_max=u_max,
    )
    n_plan_steps = max(size(r_ref, 2) - 1, 0)
    settle_steps = max(0, Int(ceil(tracking.settle_time_s / tracking.dt_s)))
    total_steps = n_plan_steps + settle_steps
    control_progress = total_steps == 0 ? Float64[] :
        collect(range(0.0, 1.0; length=total_steps))
    desired_attitudes = _attitude_reference(
        attitude_progress, attitude_quaternions, control_progress,
    )

    state = vcat(Vector{Float64}(r_ref[:, 1]), Vector{Float64}(v_ref[:, 1]))
    state_history = zeros(6, total_steps + 1)
    command_history = zeros(3, total_steps)
    realized_history = zeros(3, total_steps)
    actual_attitude_history = zeros(4, total_steps + 1)
    state_history[:, 1] .= state
    attitude = _quat_normalize(scenario.initial_attitude_rtn_to_body)
    angular_rate = zeros(3)
    actual_attitude_history[:, 1] .= attitude
    wheel_momentum = zeros(3)
    wheel_peak = 0.0
    wheel_energy = 0.0
    impulse = zeros(6)
    propellant = 0.0
    allocation_error = 0.0
    saturated_steps = 0
    min_clearance = Inf
    keepout_violations = 0
    dt = Float64(tracking.dt_s)

    for step_index in 1:total_steps
        reference_index = min(step_index, size(r_ref, 2))
        preview = Base.invokelatest(
            getproperty(modules.guidance, :rpo_lqmpc_reference_preview),
            r_ref, v_ref, reference_index, tracking.horizon,
        )
        command = Vector{Float64}(Base.invokelatest(
            getproperty(modules.control, :rpo_lqmpc_control),
            controller, state, preview,
        ))
        attitude, angular_rate, torque = _integrate_attitude(
            attitude, angular_rate, view(desired_attitudes, :, step_index), dt, config,
        )
        wheel_momentum .-= torque .* dt
        wheel_peak = max(wheel_peak, maximum(abs, wheel_momentum))
        wheel_energy += abs(dot(torque, angular_rate)) * dt

        actuation = _thruster_translation_step(command, attitude, dt, config)
        realized_acceleration = actuation.realized_acceleration_rtn
        impulse .+= actuation.thruster_impulse_ns
        propellant += actuation.propellant_used_kg
        allocation_error += actuation.allocation_error_mps2 * dt
        saturated_steps += actuation.saturated ? 1 : 0
        state .= controller.Ad * state .+ controller.Bd * realized_acceleration
        clearance = Float64(Base.invokelatest(
            getproperty(modules.navigation, :rpo_clearance_distance_to_station),
            view(state, 1:3), scenario.geometry,
        ))
        min_clearance = min(min_clearance, clearance)
        keepout_violations += clearance + 1.0e-9 < config.safe_distance_m ? 1 : 0
        state_history[:, step_index + 1] .= state
        command_history[:, step_index] .= command
        realized_history[:, step_index] .= realized_acceleration
        actual_attitude_history[:, step_index + 1] .= attitude
    end
    final_error = norm(state[1:3] - scenario.goal_rtn)
    wheel_feasible = wheel_peak <= config.reaction_wheel_max_momentum_nms + 1.0e-12
    success = final_error <= tracking.final_position_tol_m &&
              keepout_violations == 0 && wheel_feasible
    return (
        success=success,
        final_position_error_m=final_error,
        min_clearance_m=min_clearance,
        keepout_violations=keepout_violations,
        propellant_used_kg=propellant,
        allocation_error_impulse_mps=allocation_error,
        saturation_fraction=total_steps > 0 ? saturated_steps / total_steps : 0.0,
        thruster_impulse_ns=impulse,
        wheel_energy_j=wheel_energy,
        wheel_peak_momentum_nms=wheel_peak,
        wheel_feasible=wheel_feasible,
        state_history=state_history,
        command_history=command_history,
        realized_acceleration_history=realized_history,
        actual_attitude_history=actual_attitude_history,
        final_attitude_rtn_to_body=attitude,
        final_angular_rate_radps=angular_rate,
    )
end

"""
Evaluate a translation/attitude candidate after path retiming and LQ-MPC
tracking. Translation is allocated through all six fixed body thrusters; the
attitude is tracked using reaction-wheel torque and momentum limits.
"""
function evaluate_rpo_candidate(scenario::RPOHyPRRLScenario,
                                config::RPOHyPRRLConfig,
                                control_points_rtn::AbstractMatrix{<:Real},
                                attitude_progress::AbstractVector{<:Real},
                                attitude_quaternions::AbstractMatrix{<:Real})
    try
        modules, pso_config, tracking, _ = _rpo_spaceagora_settings(scenario, config)
        components = Base.invokelatest(
            getproperty(modules.guidance, :rpo_normalized_path_cost_components),
            control_points_rtn, scenario.geometry, pso_config;
            safe_distance_m=config.safe_distance_m,
        )
        if components.violation_count > 0 || !isfinite(components.total)
            return _failed_rpo_evaluation(reason=:path_infeasible)
        end
        retime_config = Base.invokelatest(
            getproperty(modules.guidance, :rpo_pso_config), pso_config;
            retime_dt_s=tracking.dt_s, mass_kg=tracking.mass_kg,
            isp_s=tracking.isp_s, g0_mps2=tracking.g0_mps2,
        )
        t_ref, r_ref, v_ref = Base.invokelatest(
            getproperty(modules.guidance, :rpo_reference_from_path),
            control_points_rtn, scenario.geometry, retime_config;
            safe_distance_m=config.safe_distance_m,
        )
        coupled = _coupled_rpo_tracking(
            modules, t_ref, r_ref, v_ref, scenario.goal_rtn,
            attitude_progress, attitude_quaternions, tracking, scenario, config,
        )
        reference_progress = isempty(t_ref) ? Float64[] :
            collect(range(0.0, 1.0; length=length(t_ref)))
        q_reference = _attitude_reference(
            attitude_progress, attitude_quaternions, reference_progress,
        )
        propellant_fraction = coupled.propellant_used_kg /
                              max(tracking.propellant_mass_kg, 1.0e-9)
        duration_scale = max(pso_config.tf_s, 1.0)
        objective = components.total +
            config.fuel_weight * propellant_fraction^2 +
            config.duration_weight * ((isempty(t_ref) ? 0.0 : t_ref[end]) / duration_scale)^2 +
            config.allocation_error_weight * coupled.allocation_error_impulse_mps^2 +
            config.wheel_weight * coupled.wheel_peak_momentum_nms /
                max(config.reaction_wheel_max_momentum_nms, 1.0e-9)
        return RPOHyPRRLEvaluation(
            coupled.success, objective, components.total, coupled.propellant_used_kg,
            isempty(t_ref) ? 0.0 : t_ref[end],
            coupled.allocation_error_impulse_mps,
            coupled.saturation_fraction,
            coupled.wheel_energy_j,
            coupled.wheel_peak_momentum_nms,
            coupled.min_clearance_m,
            coupled.final_position_error_m,
            Vector{Float64}(t_ref),
            Matrix{Float64}(r_ref),
            Matrix{Float64}(v_ref),
            q_reference,
            coupled.thruster_impulse_ns,
            (
                path_components=components,
                coupled_tracking_success=coupled.success,
                keepout_violations=coupled.keepout_violations,
                wheel_feasible=coupled.wheel_feasible,
                state_history=coupled.state_history,
                command_history=coupled.command_history,
                realized_acceleration_history=coupled.realized_acceleration_history,
                actual_attitude_history=coupled.actual_attitude_history,
                final_attitude_rtn_to_body=coupled.final_attitude_rtn_to_body,
                final_angular_rate_radps=coupled.final_angular_rate_radps,
            ),
        )
    catch error
        return _failed_rpo_evaluation(
            reason=:evaluation_exception,
            error_message=sprint(showerror, error),
        )
    end
end
