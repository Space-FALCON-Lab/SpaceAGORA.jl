include(joinpath(@__DIR__, "Earth_RPO_CubeSat_MPC_Batch.jl"))

const RPO_RL_PROJECT = joinpath(REPO_ROOT, "SpaceAGORA_RL.jl")
RPO_RL_PROJECT in LOAD_PATH || pushfirst!(LOAD_PATH, RPO_RL_PROJECT)
using SpaceAGORA_RL

function _rpo_new_hypr_config(cfg)
    return RPOHyPRRLConfig(
        safe_distance_m=cfg.safe_distance_m,
        wheel_weight=0.0,
    )
end

function _rpo_new_hypr_scenario(case, geometry, cfg; pso_config=cfg.pso_config)
    return RPOHyPRRLScenario(
        start_rtn=Vector{Float64}(case.start_rtn),
        goal_rtn=Vector{Float64}(case.goal_rtn),
        geometry=geometry,
        pso_config=pso_config,
        tracking_settings=cfg.tracking,
    )
end

function _rpo_new_hypr_objective_evaluator(case, geometry, cfg)
    config = _rpo_new_hypr_config(cfg)
    scenario = _rpo_new_hypr_scenario(case, geometry, cfg)
    return RPOTimingAwareProxyPSOObjectiveEvaluator(config, scenario)
end

const RPO_50_OF_50_EVALUATION_SEED = 10_000_740
const RPO_50_OF_50_GEOMETRY_SEED = 740

function _rpo_50_of_50_case_generator(n_cases::Int, geometry)
    base_scenario = RPOHyPRRLScenario(
        start_rtn=zeros(3),
        goal_rtn=zeros(3),
        geometry=geometry,
    )
    sampler = build_rpo_hypr_rl_endpoint_sampler(
        base_scenario;
        station_asset=:gateway,
        safe_distance_m=SM.GuidanceHooks.RPO_PLANNER_COMPARISON_SAFE_DISTANCE_M,
        endpoint_clearance_margin_m=0.05,
        endpoint_max_clearance_m=1.0,
        min_separation_m=1.5,
        surrounded_max_distance_m=2.0,
        max_sampling_tries=4_000,
    )
    return [begin
        scenario_seed = RPO_50_OF_50_EVALUATION_SEED + 2 * case_index
        scenario = sample_rpo_hypr_rl_scenario(
            sampler, MersenneTwister(scenario_seed),
        )
        (
            case_id=case_index,
            seed=scenario_seed,
            start_rtn=scenario.start_rtn,
            goal_rtn=scenario.goal_rtn,
        )
    end for case_index in 1:n_cases]
end

function _rpo_50_of_50_planner_rng(_, case_index, _case)
    planner_seed = RPO_50_OF_50_EVALUATION_SEED + 2 * case_index + 1
    return MersenneTwister(planner_seed)
end

function _rpo_actuator_aware_tracking_evaluator(plan, case, geometry, cfg)
    config = _rpo_new_hypr_config(cfg)
    scenario = _rpo_new_hypr_scenario(
        case, geometry, cfg; pso_config=plan.config,
    )
    identity_quaternion = [0.0, 0.0, 0.0, 1.0]
    evaluation = evaluate_rpo_candidate(
        scenario,
        config,
        plan.path,
        [0.0, 1.0],
        hcat(identity_quaternion, identity_quaternion);
        optimize_pointing=true,
    )
    diagnostics = evaluation.diagnostics
    x_hist = get(diagnostics, :state_history, zeros(6, 0))
    u_hist = get(diagnostics, :command_history, zeros(3, 0))
    control_effort = sum(
        norm(view(u_hist, :, step_index)) * cfg.tracking.dt_s
        for step_index in axes(u_hist, 2);
        init=0.0,
    )
    actual_steps = max(size(x_hist, 2) - 1, 0)
    fuel_used_pct = 100.0 * evaluation.propellant_used_kg /
        max(cfg.tracking.propellant_mass_kg, eps(Float64))
    return (
        success=evaluation.feasible,
        fuel_used=evaluation.propellant_used_kg,
        fuel_used_pct=fuel_used_pct,
        control_effort_total=control_effort,
        translational_control_effort=control_effort,
        thrust_saturation_fraction=evaluation.thruster_saturation_fraction,
        min_clearance=evaluation.min_clearance_m,
        keepout_violations=Int(get(diagnostics, :keepout_violations, -1)),
        final_pos_error=evaluation.final_position_error_m,
        planned_travel_duration=evaluation.duration_s,
        actual_travel_duration=cfg.tracking.dt_s * actual_steps,
        t_ref=evaluation.t_ref_s,
        r_ref_rtn=evaluation.r_ref_rtn,
        v_ref_rtn=evaluation.v_ref_rtn,
        x_hist=x_hist,
        u_hist=u_hist,
    )
end

function run_rpo_new_hypr_planner_comparison_batch_cli(args::Vector{String}=copy(ARGS))
    kwargs = _rpo_planner_comparison_cli_kwargs(args)
    kwargs === nothing && return nothing
    return run_rpo_cubesat_mpc_planner_comparison_batch(
        ;
        kwargs...,
        seed=RPO_50_OF_50_EVALUATION_SEED,
        geometry_seed=RPO_50_OF_50_GEOMETRY_SEED,
        pso_sample_ds_m=_env_float(
            "SPACEAGORA_RPO_COMPARISON_HYPR_SAMPLE_DS", 0.5,
        ),
        pso_clearance_feasibility_tol_m=_env_float(
            "SPACEAGORA_RPO_COMPARISON_CLEARANCE_TOL", 0.005,
        ),
        mpc_horizon=60,
        mpc_u_max_mps2=0.0125,
        mpc_q_pos=10.0,
        mpc_q_vel=1.0,
        mpc_r_accel=0.1,
        mpc_qf_pos=50.0,
        mpc_qf_vel=5.0,
        mpc_settle_time_s=20.0,
        mpc_final_position_tol_m=0.25,
        hypr_objective_evaluator_factory=_rpo_new_hypr_objective_evaluator,
        tracking_evaluator=_rpo_actuator_aware_tracking_evaluator,
        comparison_case_generator=_rpo_50_of_50_case_generator,
        planner_rng_factory=_rpo_50_of_50_planner_rng,
        share_rrt_initial_path=true,
        hypr_objective_name=
            "retimed_hcw_continuous_thrust_linear_proxy_fuel_only",
        tracking_evaluator_name=
            "pwm_six_thruster_command_aligned_hcw_lqmpc_nonlinear_two_body",
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_rpo_new_hypr_planner_comparison_batch_cli(copy(ARGS))
end
