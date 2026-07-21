module SpaceAGORA_RL

using Dates
using LinearAlgebra
using Printf
using Random
using Serialization
using SHA
using Statistics
using TOML

include("core/mdp_interface.jl")
include("core/backend_interface.jl")
include("core/policy_interface.jl")
include("core/transition_schema.jl")

include("tasks/aerobraking/mdp/actions.jl")
include("tasks/aerobraking/mdp/observations.jl")
include("tasks/aerobraking/mdp/normalization.jl")
include("tasks/aerobraking/mdp/rewards.jl")
include("tasks/aerobraking/mdp/termination.jl")
include("tasks/aerobraking/mdp/randomization.jl")

include("tasks/aerobraking/backend/maneuver_model.jl")
include("tasks/aerobraking/backend/simulation_config_builder.jl")
include("tasks/aerobraking/backend/pass_metrics.jl")
include("tasks/aerobraking/backend/paper_heat_rate_metrics.jl")
include("tasks/aerobraking/backend/spaceagora_core_adapter.jl")
include("tasks/aerobraking/backend/spaceagora_backend.jl")

include("tasks/aerobraking/baselines/no_maneuver.jl")
include("tasks/aerobraking/baselines/random_policy.jl")
include("tasks/aerobraking/baselines/fixed_corridor.jl")
include("tasks/aerobraking/baselines/aads_heuristic.jl")

include("tasks/aerobraking/evaluation/episode_logs.jl")
include("tasks/aerobraking/evaluation/metrics.jl")
include("tasks/aerobraking/evaluation/odyssey_comparison.jl")
include("tasks/aerobraking/evaluation/generalization_suites.jl")
include("tasks/aerobraking/evaluation/reports.jl")

include("tasks/aerobraking/AerobrakingMDP.jl")

include("algorithms/training_device.jl")
include("algorithms/ddqn/replay_buffer.jl")
include("algorithms/ddqn/network.jl")
include("algorithms/ddqn/epsilon_schedule.jl")
include("algorithms/ddqn/target_network.jl")
include("algorithms/ddqn/learner.jl")
include("algorithms/a2c/rollout_buffer.jl")
include("algorithms/a2c/learner.jl")
include("algorithms/ddqn/checkpoints.jl")
include("algorithms/ddqn/policy.jl")
include("algorithms/a2c/policy.jl")

include("runtime/config_loading.jl")
include("runtime/run_manifest.jl")
include("runtime/training_session.jl")

include("parallel/rollout_protocol.jl")
include("parallel/worker.jl")
include("parallel/master.jl")
include("parallel/distributed_backend.jl")

export AbstractMDP, AbstractRLBackend, AbstractPolicy
export reset_scenario, step_scenario, observe_state, normalize_observation
export policy_action_index, Transition, EpisodeSummary, empty_episode_summary

export PAPER_ACTIONS_MPS, AerobrakingAction, action_from_index, action_from_delta_v, nearest_action_index
export zero_action_index, action_count
export PaperObservation, paper_observation_names, raw_observation_vector
export NormalizationBounds, paper_normalization_bounds, normalize_value
export RewardConfig, TerminationConfig, ThermalStatus, thermal_low, thermal_nominal
export thermal_high, thermal_medium, thermal_hard, thermal_status, paper_reward
export TerminationFlags, classify_termination
export AerobrakingRandomizationConfig

export AerobrakingScenarioConfig, AerobrakingDecisionState, AerobrakingPassMetrics
export AerobrakingStepResult, SpaceAGORAAerobrakingBackend
export reset_scenario, step_scenario
export NoManeuverPolicy, RandomActionPolicy, FixedCorridorPolicy, AADSHeuristicPolicy

export evaluate_policy, evaluate_baselines, episode_metrics, aggregate_metrics
export paper_pr_drl_evaluation_config, paper_pr_drl_marsgram_evaluation_config
export paper_pr_drl_physics_evaluation_config, paper_odyssey_flight_evaluation_config
export generalization_suite_configs, write_evaluation_artifacts
export default_aerobraking_config, mars_odyssey_phase_constants

export ReplayBuffer, ReplayBatch, sample_batch
export QNetwork, init_q_network, predict_q, copy_network!
export EpsilonSchedule, epsilon_value
export DDQNConfig, DDQNLearner, select_action, observe!, maybe_train!, compute_ddqn_targets
export A2CConfig, A2CLearner, A2CRolloutBatch, compute_discounted_returns
export save_checkpoint, load_checkpoint, GreedyDDQNPolicy, GreedyPRDRLPolicy
export load_trained_ddqn_policy, load_trained_pr_drl_policy
export GreedyA2CPolicy, load_trained_a2c_policy
export CPUTrainingDevice, CUDATrainingDevice, resolve_training_device, training_device_name

export ResolvedConfig, load_config, resolve_config, default_config_path
export RunManifest, write_manifest, TrainingSession, build_training_session
export train_parallel!, setup_distributed_workers

end # module SpaceAGORA_RL
