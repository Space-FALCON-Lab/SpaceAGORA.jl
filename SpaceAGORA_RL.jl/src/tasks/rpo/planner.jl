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
                      rng::AbstractRNG=Random.default_rng(), test::Bool=true)
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
    evaluation = state.best_evaluation
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
            evaluator=evaluation.diagnostics,
        ),
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
