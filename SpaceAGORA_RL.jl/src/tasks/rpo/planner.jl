struct GreedyHyPRRLPolicy <: AbstractPolicy
    network::QNetwork
end

function policy_action_index(policy::GreedyHyPRRLPolicy, mdp::RPOHyPRRLMDP,
                             state::RPOHyPRRLState, observation, rng)
    values = predict_q(policy.network, observation)
    return _masked_argmax(values, valid_action_mask(mdp, state))
end

function load_hypr_rl_policy(path::AbstractString)
    payload = load_checkpoint(path)
    get(payload, :task, nothing) == :rpo_hypr_rl ||
        throw(ArgumentError("checkpoint is not a HyPR-RL checkpoint"))
    return GreedyHyPRRLPolicy(payload[:online])
end

function _hypr_rl_policy_action(policy, mdp, state, observation, mask, rng, test)
    if policy isa DDQNLearner
        return select_action(policy, observation; rng=rng, test=test, action_mask=mask)
    elseif policy isa Function
        return Int(policy(observation, mask, state))
    else
        return Int(policy_action_index(policy, mdp, state, observation, rng))
    end
end

"""Run the sequential HyPR-RL editor and return its best feasible plan."""
function hypr_rl_plan(mdp::RPOHyPRRLMDP, policy;
                      rng::AbstractRNG=Random.default_rng(), test::Bool=true,
                      terminal_evaluator=mdp.evaluator)
    state = reset_scenario(mdp, rng)
    while !state.stopped && state.edit_count < mdp.config.max_edits
        observation = observe_state(mdp, state)
        mask = valid_action_mask(mdp, state)
        action = _hypr_rl_policy_action(
            policy, mdp, state, observation, mask, rng, test,
        )
        result = step_scenario(mdp, state, action, rng)
        state = result.state
        (result.terminated || result.truncated) && break
    end
    edit_evaluation = state.best_evaluation
    evaluation = terminal_evaluator === mdp.evaluator ? edit_evaluation :
        Base.invokelatest(
            terminal_evaluator,
            mdp.scenario, mdp.config, state.best_control_points_rtn,
            state.best_attitude_progress, state.best_attitude_quaternions,
        )
    return RPOHyPRRLPlan(
        evaluation.feasible,
        state.best_control_points_rtn,
        evaluation.t_ref_s,
        evaluation.r_ref_rtn,
        evaluation.v_ref_rtn,
        evaluation.q_ref_rtn_to_body,
        evaluation.propellant_used_kg,
        evaluation.objective,
        state.seed_evaluation.objective,
        state.edit_count,
        state.stopped,
        (
            path_cost=evaluation.path_cost,
            duration_s=evaluation.duration_s,
            allocation_error_impulse_mps=evaluation.allocation_error_impulse_mps,
            thruster_saturation_fraction=evaluation.thruster_saturation_fraction,
            thruster_impulse_ns=evaluation.thruster_impulse_ns,
            wheel_energy_j=evaluation.wheel_energy_j,
            wheel_peak_momentum_nms=evaluation.wheel_peak_momentum_nms,
            min_clearance_m=evaluation.min_clearance_m,
            final_position_error_m=evaluation.final_position_error_m,
            evaluator=evaluation.diagnostics,
            edit_feasible=edit_evaluation.feasible,
            edit_cost=edit_evaluation.objective,
            edit_min_clearance_m=edit_evaluation.min_clearance_m,
            edit_evaluator=edit_evaluation.diagnostics,
        ),
    )
end

"""Evaluate one frozen HyPR-RL policy case with fast edits and one full terminal solve."""
function evaluate_hypr_rl_policy_case(config::RPOHyPRRLConfig,
                                      scenario::RPOHyPRRLScenario,
                                      policy, seed::Integer)
    mdp = RPOHyPRRLMDP(
        config, scenario; evaluator=evaluate_rpo_training_candidate,
    )
    elapsed_s = @elapsed plan = hypr_rl_plan(
        mdp, policy;
        rng=MersenneTwister(seed),
        test=true,
        terminal_evaluator=evaluate_rpo_candidate,
    )
    return (plan=plan, runtime_s=elapsed_s)
end

function _evaluate_hypr_pso_case(config::RPOHyPRRLConfig,
                                 scenario::RPOHyPRRLScenario,
                                 seed::Integer,
                                 modules,
                                 pso_config;
                                 planner::Symbol,
                                 objective_evaluator=nothing)
    planner_keywords = objective_evaluator === nothing ?
        (
            safe_distance_m=config.safe_distance_m,
            rng=MersenneTwister(seed),
        ) :
        (
            safe_distance_m=config.safe_distance_m,
            rng=MersenneTwister(seed),
            objective_evaluator=objective_evaluator,
        )
    planner_runtime_s = @elapsed plan_result = Base.invokelatest(
        getproperty(modules.guidance, :rpo_pso_plan_path),
        scenario.start_rtn,
        scenario.goal_rtn,
        scenario.geometry,
        pso_config;
        planner_keywords...,
    )
    progress, quaternions = _initial_attitude_knots(scenario)
    terminal_scenario = RPOHyPRRLScenario(
        start_rtn=scenario.start_rtn,
        goal_rtn=scenario.goal_rtn,
        geometry=scenario.geometry,
        pso_config=plan_result.config,
        tracking_settings=scenario.tracking_settings,
        rrt_settings=scenario.rrt_settings,
        initial_attitude_rtn_to_body=scenario.initial_attitude_rtn_to_body,
        final_attitude_rtn_to_body=scenario.final_attitude_rtn_to_body,
    )
    terminal_runtime_s = @elapsed evaluation = evaluate_rpo_candidate(
        terminal_scenario,
        config,
        plan_result.path,
        progress,
        quaternions,
    )
    seed_cost = isempty(plan_result.cost_history) ? plan_result.cost :
        first(plan_result.cost_history)
    uses_retimed_objective = objective_evaluator !== nothing
    edit_feasible = uses_retimed_objective ? Bool(get(
        plan_result.components, :retimed_feasible, false,
    )) :
        plan_result.components.violation_count == 0 &&
        plan_result.components.min_clearance + 1.0e-9 >= config.safe_distance_m &&
        isfinite(plan_result.cost)
    edit_evaluator_mode = uses_retimed_objective ?
        :retimed_feedforward_pso : :original_hypr_proxy_pso
    plan = RPOHyPRRLPlan(
        evaluation.feasible,
        Matrix{Float64}(plan_result.path),
        evaluation.t_ref_s,
        evaluation.r_ref_rtn,
        evaluation.v_ref_rtn,
        evaluation.q_ref_rtn_to_body,
        evaluation.propellant_used_kg,
        evaluation.objective,
        seed_cost,
        0,
        true,
        (
            planner=planner,
            path_cost=evaluation.path_cost,
            duration_s=evaluation.duration_s,
            allocation_error_impulse_mps=evaluation.allocation_error_impulse_mps,
            thruster_saturation_fraction=evaluation.thruster_saturation_fraction,
            thruster_impulse_ns=evaluation.thruster_impulse_ns,
            wheel_energy_j=evaluation.wheel_energy_j,
            wheel_peak_momentum_nms=evaluation.wheel_peak_momentum_nms,
            min_clearance_m=evaluation.min_clearance_m,
            final_position_error_m=evaluation.final_position_error_m,
            evaluator=evaluation.diagnostics,
            edit_feasible=edit_feasible,
            edit_cost=plan_result.cost,
            edit_min_clearance_m=plan_result.components.min_clearance,
            edit_evaluator=(
                evaluator_mode=edit_evaluator_mode,
                path_components=plan_result.components,
            ),
            pso_iterations=length(plan_result.cost_history),
            pso_refinement_improved=plan_result.refinement_improved,
            pso_early_stopped=plan_result.early_stopped,
            planner_runtime_s=planner_runtime_s,
            terminal_runtime_s=terminal_runtime_s,
        ),
    )
    return (
        plan=plan,
        runtime_s=planner_runtime_s + terminal_runtime_s,
        planner_runtime_s=planner_runtime_s,
        terminal_runtime_s=terminal_runtime_s,
    )
end

"""Evaluate baseline HYPR with the original geometric and proxy-fuel PSO objective."""
function evaluate_hypr_original_baseline_case(config::RPOHyPRRLConfig,
                                               scenario::RPOHyPRRLScenario,
                                               seed::Integer)
    modules = _spaceagora_rpo_modules()
    pso_config = scenario.pso_config === nothing ?
        Base.invokelatest(
            getproperty(modules.guidance, :rpo_740_mpc_final_pso_config);
            safe_distance_m=config.safe_distance_m,
        ) : scenario.pso_config
    pso_config = Base.invokelatest(
        getproperty(modules.guidance, :rpo_pso_config),
        pso_config;
        curve_type=:bezier,
        safe_distance_m=config.safe_distance_m,
    )
    return _evaluate_hypr_pso_case(
        config,
        scenario,
        seed,
        modules,
        pso_config;
        planner=:hypr_pso_original_proxy,
    )
end

"""Evaluate baseline HYPR with the same retimed fuel objective and terminal model as HyPR-RL."""
function evaluate_hypr_pso_baseline_case(config::RPOHyPRRLConfig,
                                         scenario::RPOHyPRRLScenario,
                                         seed::Integer)
    modules, pso_config, _, _ = _rpo_spaceagora_settings(scenario, config)
    return _evaluate_hypr_pso_case(
        config,
        scenario,
        seed,
        modules,
        pso_config;
        planner=:hypr_pso_retimed_fuel,
        objective_evaluator=RPOHyPRRLPSOObjectiveEvaluator(config, scenario),
    )
end


"""Evaluate original and retimed-fuel HYPR on one paired scenario and PSO seed."""
function evaluate_hypr_pso_comparison_case(config::RPOHyPRRLConfig,
                                           scenario::RPOHyPRRLScenario,
                                           seed::Integer)
    failure(error) = (
        plan=nothing,
        runtime_s=NaN,
        planner_runtime_s=NaN,
        terminal_runtime_s=NaN,
        error=sprint(showerror, error),
    )
    original = try
        evaluate_hypr_original_baseline_case(config, scenario, seed)
    catch error
        failure(error)
    end
    retimed_fuel = try
        evaluate_hypr_pso_baseline_case(config, scenario, seed)
    catch error
        failure(error)
    end
    return (
        original=original,
        retimed_fuel=retimed_fuel,
    )
end

"""SpaceAGORA-style convenience entry point for the HyPR-RL path planner."""
function rpo_hypr_rl_plan_path(start_rtn, goal_rtn, geometry, policy;
                               config::RPOHyPRRLConfig=RPOHyPRRLConfig(),
                               pso_config=nothing,
                               tracking_settings=nothing,
                               rrt_settings=nothing,
                               seed_path_rtn=nothing,
                               initial_attitude_rtn_to_body=[0.0, 0.0, 0.0, 1.0],
                               final_attitude_rtn_to_body=[0.0, 0.0, 0.0, 1.0],
                               evaluator=evaluate_rpo_candidate,
                               rng::AbstractRNG=Random.default_rng(),
                               test::Bool=true)
    scenario = RPOHyPRRLScenario(
        start_rtn=Vector{Float64}(start_rtn),
        goal_rtn=Vector{Float64}(goal_rtn),
        geometry=geometry,
        pso_config=pso_config,
        tracking_settings=tracking_settings,
        rrt_settings=rrt_settings,
        seed_path_rtn=seed_path_rtn === nothing ? nothing : Matrix{Float64}(seed_path_rtn),
        initial_attitude_rtn_to_body=Vector{Float64}(initial_attitude_rtn_to_body),
        final_attitude_rtn_to_body=Vector{Float64}(final_attitude_rtn_to_body),
    )
    mdp = RPOHyPRRLMDP(config, scenario; evaluator=evaluator)
    return hypr_rl_plan(mdp, policy; rng=rng, test=test)
end
