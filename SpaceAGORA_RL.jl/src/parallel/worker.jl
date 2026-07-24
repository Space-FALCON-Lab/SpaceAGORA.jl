function learner_policy_action!(learner::DDQNLearner, observation::Vector{Float32}, rng::AbstractRNG; test::Bool=false)
    return select_action(learner, observation; rng=rng, test=test)
end

function snapshot_policy_action(network::QNetwork, schedule::EpsilonSchedule,
                                ddqn_config::DDQNConfig,
                                observation::Vector{Float32}, step::Integer,
                                rng::AbstractRNG; test::Bool=false)
    eps = test ? 0.0 : epsilon_value(schedule, step)
    if !test && rand(rng) < eps
        return rand(rng, 1:ddqn_config.action_dim)
    end
    return argmax(predict_q(network, observation))
end

Base.@kwdef struct SpaceAGORAPhysicsPassEvent
    worker_id::Int
    episode_index::Int
    seed::Int
    transition::Union{Transition,Nothing} = nothing
    summary::EpisodeSummary = EpisodeSummary()
    result::Union{AerobrakingStepResult,Nothing} = nothing
    done::Bool = false
    protected::Bool = false
    error::Union{Nothing,String} = nothing
end

struct SpaceAGORAPhysicsWorkerHandle
    task::Task
    action_channel::Channel{Any}
end

Base.@kwdef mutable struct SpaceAGORAPhysicsCampaignRollout
    config::AerobrakingScenarioConfig
    schedule::EpsilonSchedule
    ddqn_config::DDQNConfig
    policy_snapshot::Union{Nothing,QNetwork}
    rng::MersenneTwister
    episode_index::Int
    worker_id::Int
    seed::Int
    max_passes_per_campaign::Int
    global_step_start::Int
    train::Bool = true
    state::AerobrakingDecisionState
    norm_obs::Vector{Float32}
    action_index::Int
    action::AerobrakingAction
    summary::EpisodeSummary
    transitions::Vector{Transition}
    stats::SpaceAGORAPhysicsPassStats = SpaceAGORAPhysicsPassStats()
    pass_start_time_s::Float64 = 0.0
    terminated::Bool = false
    streaming::Bool = false
    event_channel::Union{Nothing,Channel{Any}} = nothing
    action_channel::Union{Nothing,Channel{Any}} = nothing
    protected_next_transition::Bool = false
    protected_suppress_thermal_terminal::Bool = true
    last_integrator::Any = nothing
end

function _spaceagora_physics_campaign_mission_pass_cap(config::AerobrakingScenarioConfig,
                                                       max_passes_per_campaign::Integer)
    requested = max(1, Int(max_passes_per_campaign))
    configured = max(1, config.termination_config.max_passes)
    return min(requested, configured)
end

function _spaceagora_physics_campaign_step_index(rollout::SpaceAGORAPhysicsCampaignRollout)
    return rollout.global_step_start + length(rollout.transitions) + 1
end

function _spaceagora_physics_campaign_set_action!(rollout::SpaceAGORAPhysicsCampaignRollout,
                                                  action_index::Integer)
    rollout.action_index = Int(action_index)
    rollout.action = action_from_index(rollout.action_index)
    return rollout.action
end

function _spaceagora_physics_campaign_decode_action_command(command)
    if command isa NamedTuple
        hasproperty(command, :action_index) ||
            throw(ArgumentError("SpaceAGORA physics action command is missing action_index."))
        protected = hasproperty(command, :protected) ? Bool(getproperty(command, :protected)) : false
        return Int(getproperty(command, :action_index)), protected
    end
    return Int(command), false
end

function _spaceagora_physics_campaign_select_action!(rollout::SpaceAGORAPhysicsCampaignRollout)
    step = _spaceagora_physics_campaign_step_index(rollout)
    policy_snapshot = rollout.policy_snapshot
    policy_snapshot === nothing &&
        throw(ArgumentError("non-streaming campaign rollout requires a policy snapshot"))
    action_index = snapshot_policy_action(
        policy_snapshot,
        rollout.schedule,
        rollout.ddqn_config,
        rollout.norm_obs,
        step,
        rollout.rng;
        test = !rollout.train,
    )
    return _spaceagora_physics_campaign_set_action!(rollout, action_index)
end

function _spaceagora_physics_campaign_emit!(rollout::SpaceAGORAPhysicsCampaignRollout,
                                            transition::Union{Transition,Nothing},
                                            result::Union{AerobrakingStepResult,Nothing},
                                            done::Bool;
                                            protected::Bool=false,
                                            error::Union{Nothing,String}=nothing)
    rollout.streaming || return nothing
    channel = rollout.event_channel
    channel === nothing && return nothing
    put!(channel, SpaceAGORAPhysicsPassEvent(
        worker_id = rollout.worker_id,
        episode_index = rollout.episode_index,
        seed = rollout.seed,
        transition = transition,
        summary = rollout.summary,
        result = result,
        done = done,
        protected = protected,
        error = error,
    ))
    return nothing
end

function _spaceagora_protected_first_pass_flags(flags::TerminationFlags)
    thermal_only_terminal = flags.thermal_violation &&
                            flags.terminated &&
                            !(flags.success ||
                              flags.target_undershoot ||
                              flags.impact ||
                              flags.out_of_drag_passage)
    thermal_only_terminal || return flags
    return TerminationFlags(
        flags.success,
        flags.target_undershoot,
        flags.impact,
        flags.out_of_drag_passage,
        flags.thermal_violation,
        false,
        flags.truncated,
    )
end

function _spaceagora_physics_campaign_mark_modified!(spaceagora, integrator)
    callbacks = getproperty(getproperty(spaceagora, :SimulationModel), :SimulationCallbacks)
    try
        if isdefined(callbacks, :u_modified!)
            return Base.invokelatest(getproperty(callbacks, :u_modified!), integrator, true)
        end
    catch
    end
    try
        setproperty!(integrator, :u_modified, true)
    catch
    end
    return nothing
end

function _spaceagora_physics_campaign_reset_pass_state!(integrator)
    try
        sc = _spaceagora_solution_satellite_state(getproperty(integrator, :u))
        if hasproperty(sc, :heat_loads)
            fill!(getproperty(sc, :heat_loads), 0.0)
        end
    catch
    end
    try
        heat_rates = getproperty(getproperty(getproperty(integrator, :p), :shared_buffers), :heat_rates)
        if !isempty(heat_rates)
            fill!(heat_rates[1], 0.0)
        end
    catch
    end
    return nothing
end

function _spaceagora_physics_campaign_apply_action!(spaceagora,
                                                    rollout::SpaceAGORAPhysicsCampaignRollout,
                                                    integrator)
    _spaceagora_physics_campaign_reset_pass_state!(integrator)
    abs(rollout.action.delta_v_mps) <= 1e-12 && return nothing
    sc = _spaceagora_solution_satellite_state(getproperty(integrator, :u))
    vel = getproperty(sc, :vel)
    speed = norm(vel)
    speed > 0.0 || throw(ErrorException("SpaceAGORA campaign action cannot be applied to a zero-velocity state."))
    new_speed = max(speed + rollout.action.delta_v_mps, eps(Float64))
    vel .*= new_speed / speed
    _spaceagora_physics_campaign_mark_modified!(spaceagora, integrator)
    return nothing
end

function _spaceagora_physics_campaign_terminate!(spaceagora, integrator)
    callbacks = getproperty(getproperty(spaceagora, :SimulationModel), :SimulationCallbacks)
    terminate! = Base.invokelatest(getproperty, callbacks, :terminate!)
    return Base.invokelatest(terminate!, integrator)
end

function _spaceagora_physics_campaign_push_transition!(spaceagora,
                                                       rollout::SpaceAGORAPhysicsCampaignRollout,
                                                       args,
                                                       final_u,
                                                       final_time_s::Real;
                                                       force_truncated::Bool=false)
    length(rollout.transitions) < rollout.max_passes_per_campaign || return nothing
    elapsed_pass_s = max(0.0, Float64(final_time_s) - rollout.pass_start_time_s)
    periapsis_after_maneuver = periapsis_after_action_m(rollout.config, rollout.state, rollout.action)
    next_state = _spaceagora_physics_next_state_from_u(
        spaceagora,
        rollout.config,
        rollout.state,
        rollout.action,
        args,
        rollout.stats,
        final_u,
        elapsed_pass_s,
        periapsis_after_maneuver,
    )
    obs = observe_state(rollout.config, next_state)
    flags = classify_termination(obs, rollout.config;
                                  training=rollout.config.training,
                                  pass_count=next_state.pass_index)
    if rollout.protected_next_transition && rollout.protected_suppress_thermal_terminal
        flags = _spaceagora_protected_first_pass_flags(flags)
    end
    if force_truncated && !(flags.terminated || flags.truncated)
        flags = TerminationFlags(
            flags.success,
            flags.target_undershoot,
            flags.impact,
            flags.out_of_drag_passage,
            flags.thermal_violation,
            flags.terminated,
            true,
        )
    end
    reward = paper_reward(obs, rollout.config, rollout.action, flags, rollout.config.reward_config)
    normalized = normalize_observation(obs, rollout.config.normalization_bounds)
    metrics = pass_metrics_from_state(next_state)
    result = AerobrakingStepResult(
        next_state,
        rollout.action,
        obs,
        normalized,
        reward,
        flags,
        metrics,
        nothing,
    )
    transition = transition_from_step(
        rollout.norm_obs,
        rollout.action_index,
        result,
        length(rollout.transitions) + 1,
    )
    push!(rollout.transitions, transition)
    rollout.summary = update_episode_summary(rollout.summary, result)
    rollout.state = next_state
    rollout.norm_obs = normalized
    return result
end

function _spaceagora_physics_campaign_record_apoapsis!(spaceagora,
                                                       rollout::SpaceAGORAPhysicsCampaignRollout,
                                                       integrator,
                                                       idx::Int64)
    idx == 1 || return nothing
    rollout.terminated && return nothing
    t_now = Float64(getproperty(integrator, :t))
    elapsed_since_pass_start = t_now - rollout.pass_start_time_s
    elapsed_since_pass_start > 1.0 || return nothing
    engine = getproperty(spaceagora, :SimulationEngine)
    state_position_ii = getproperty(engine, :_state_position_ii)
    pos = state_position_ii(getproperty(integrator, :u), 1)
    planet = getproperty(getproperty(getproperty(integrator, :p), :args), :environment_model).planet
    # The radial-velocity root fires at periapsis and apoapsis; the RL decision
    # point is the high-altitude apoapsis root only.
    norm(pos) > planet.Rp_e + 250e3 || return nothing
    args = getproperty(getproperty(integrator, :p), :args)
    protected_transition = rollout.protected_next_transition
    result = _spaceagora_physics_campaign_push_transition!(
        spaceagora,
        rollout,
        args,
        getproperty(integrator, :u),
        t_now,
    )
    if result === nothing
        rollout.terminated = true
        _spaceagora_physics_campaign_emit!(rollout, nothing, nothing, true)
        return _spaceagora_physics_campaign_terminate!(spaceagora, integrator)
    end
    rollout.protected_next_transition = false

    done = result.flags.terminated ||
           result.flags.truncated ||
           length(rollout.transitions) >= rollout.max_passes_per_campaign
    transition = last(rollout.transitions)
    _spaceagora_physics_campaign_emit!(rollout, transition, result, done; protected=protected_transition)

    if done
        rollout.terminated = true
        return _spaceagora_physics_campaign_terminate!(spaceagora, integrator)
    end

    rollout.stats = SpaceAGORAPhysicsPassStats()
    rollout.pass_start_time_s = Float64(getproperty(integrator, :t))
    if rollout.streaming
        command = try
            take!(rollout.action_channel::Channel{Any})
        catch
            nothing
        end
        if command === nothing
            rollout.terminated = true
            return _spaceagora_physics_campaign_terminate!(spaceagora, integrator)
        end
        action_index, protected_next = _spaceagora_physics_campaign_decode_action_command(command)
        rollout.protected_next_transition = protected_next
        _spaceagora_physics_campaign_set_action!(rollout, action_index)
    else
        _spaceagora_physics_campaign_select_action!(rollout)
    end
    _spaceagora_physics_campaign_apply_action!(spaceagora, rollout, integrator)
    return nothing
end

function _spaceagora_physics_campaign_apoapsis_callback(spaceagora,
                                                        rollout::SpaceAGORAPhysicsCampaignRollout)
    engine = getproperty(spaceagora, :SimulationEngine)
    callbacks = getproperty(getproperty(spaceagora, :SimulationModel), :SimulationCallbacks)
    vector_callback = getproperty(callbacks, :VectorContinuousCallback)
    state_position_ii = getproperty(engine, :_state_position_ii)
    state_velocity_ii = getproperty(engine, :_state_velocity_ii)
    function condition!(out, u, t, integrator)
        pos = state_position_ii(u, 1)
        vel = state_velocity_ii(u, 1)
        out[1] = -dot(pos, vel)
        return nothing
    end
    affect!(integrator, idx::Int64) =
        _spaceagora_physics_campaign_record_apoapsis!(spaceagora, rollout, integrator, idx)
    return Base.invokelatest(vector_callback, condition!, affect!, nothing, 1)
end

function _spaceagora_physics_campaign_stats_callback(spaceagora,
                                                     rollout::SpaceAGORAPhysicsCampaignRollout)
    callbacks = getproperty(getproperty(spaceagora, :SimulationModel), :SimulationCallbacks)
    discrete_callback = getproperty(callbacks, :DiscreteCallback)
    state_position_ii = getproperty(getproperty(spaceagora, :SimulationEngine), :_state_position_ii)
    condition(u, t, integrator) = true
    function affect!(integrator)
        rollout.last_integrator = integrator
        return _record_spaceagora_physics_sample!(
            state_position_ii,
            spaceagora,
            rollout.stats,
            integrator,
        )
    end
    initialize = (cb, u, t, integrator) -> affect!(integrator)
    return Base.invokelatest(discrete_callback, condition, affect!; initialize=initialize)
end

function run_spaceagora_physics_campaign_worker_episode(config::AerobrakingScenarioConfig,
                                                        schedule::EpsilonSchedule,
                                                        ddqn_config::DDQNConfig,
                                                        policy_snapshot::QNetwork,
                                                        episode_index::Int,
                                                        worker_id::Int,
                                                        seed::Int,
                                                        max_passes_per_campaign::Int,
                                                        global_step_start::Int;
                                                        train::Bool=true)
    rng = MersenneTwister(seed)
    state = reset_scenario(config, rng)
    obs = observe_state(config, state)
    norm_obs = normalize_observation(obs, config.normalization_bounds)
    summary = empty_episode_summary(episode_index=episode_index, worker_id=worker_id, seed=seed)
    pass_cap = _spaceagora_physics_campaign_mission_pass_cap(config, max_passes_per_campaign)
    action_index = snapshot_policy_action(
        policy_snapshot,
        schedule,
        ddqn_config,
        norm_obs,
        global_step_start + 1,
        rng;
        test = !train,
    )
    action = action_from_index(action_index)
    rollout = SpaceAGORAPhysicsCampaignRollout(
        config = config,
        schedule = schedule,
        ddqn_config = ddqn_config,
        policy_snapshot = policy_snapshot,
        rng = rng,
        episode_index = episode_index,
        worker_id = worker_id,
        seed = seed,
        max_passes_per_campaign = pass_cap,
        global_step_start = global_step_start,
        train = train,
        state = state,
        norm_obs = norm_obs,
        action_index = action_index,
        action = action,
        summary = summary,
        transitions = Transition[],
    )

    spaceagora = _load_spaceagora!(; load_gramsuite=_spaceagora_live_needs_gramsuite(config))
    args, _ = _spaceagora_physics_simulation_configuration(
        config,
        state,
        action;
        prediction = false,
        campaign_max_passes = pass_cap,
    )
    stats_callback = _spaceagora_physics_campaign_stats_callback(spaceagora, rollout)
    apoapsis_callback = _spaceagora_physics_campaign_apoapsis_callback(spaceagora, rollout)
    ephemeris_callback =
        _spaceagora_rl_shared_ephemeris_callback(spaceagora, config, pass_cap)
    extra_callbacks = ephemeris_callback === nothing ?
                      (stats_callback, apoapsis_callback) :
                      (ephemeris_callback, stats_callback, apoapsis_callback)
    run_simulation_fn = Base.invokelatest(getproperty, spaceagora, :run_simulation)
    try
        Base.invokelatest(
            run_simulation_fn,
            args;
            return_solution = false,
            isolate_state = false,
            extra_callbacks = extra_callbacks,
        )
    catch err
        bt = catch_backtrace()
        throw(ErrorException(
            "SpaceAGORA physics campaign rollout failed while propagating a continuous aerobraking campaign. " *
            "Original error:\n$(sprint(showerror, err, bt))"
        ))
    end

    if !rollout.terminated && length(rollout.transitions) < pass_cap
        integrator = rollout.last_integrator
        integrator === nothing &&
            throw(ErrorException("SpaceAGORA physics campaign did not expose a final integrator state"))
        args = getproperty(getproperty(integrator, :p), :args)
        _spaceagora_physics_campaign_push_transition!(
            spaceagora,
            rollout,
            args,
            getproperty(integrator, :u),
            Float64(getproperty(integrator, :t));
            force_truncated = true,
        )
    end
    return finalize_episode_summary(rollout.summary, config), rollout.transitions
end

function run_spaceagora_physics_campaign_streaming_worker_episode(event_channel::Channel{Any},
                                                                  action_channel::Channel{Any},
                                                                  config::AerobrakingScenarioConfig,
                                                                  schedule::EpsilonSchedule,
                                                                  ddqn_config::DDQNConfig,
                                                                  policy_snapshot::Union{Nothing,QNetwork},
                                                                  simulation_template::SpaceAGORAPhysicsSimulationTemplate,
                                                                  state::AerobrakingDecisionState,
                                                                  norm_obs::Vector{Float32},
                                                                  action_index::Int,
                                                                  summary::EpisodeSummary,
                                                                  episode_index::Int,
                                                                  worker_id::Int,
                                                                  seed::Int,
                                                                  max_passes_per_campaign::Int,
                                                                  global_step_start::Int;
                                                                  protected_first_pass::Bool=false,
                                                                  protected_suppress_thermal_terminal::Bool=true)
    rollout = nothing
    try
        pass_cap = _spaceagora_physics_campaign_mission_pass_cap(config, max_passes_per_campaign)
        rollout = SpaceAGORAPhysicsCampaignRollout(
            config = config,
            schedule = schedule,
            ddqn_config = ddqn_config,
            policy_snapshot = policy_snapshot,
            rng = MersenneTwister(seed),
            episode_index = episode_index,
            worker_id = worker_id,
            seed = seed,
            max_passes_per_campaign = pass_cap,
            global_step_start = global_step_start,
            train = true,
            state = state,
            norm_obs = copy(norm_obs),
            action_index = action_index,
            action = action_from_index(action_index),
            summary = summary,
            transitions = Transition[],
            streaming = true,
            event_channel = event_channel,
            action_channel = action_channel,
            protected_next_transition = protected_first_pass,
            protected_suppress_thermal_terminal = protected_suppress_thermal_terminal,
        )

        spaceagora = _load_spaceagora!(; load_gramsuite=_spaceagora_live_needs_gramsuite(config))
        args, _ = _spaceagora_physics_simulation_configuration(
            config,
            state,
            rollout.action;
            prediction = false,
            campaign_max_passes = pass_cap,
            simulation_template = simulation_template,
        )
        stats_callback = _spaceagora_physics_campaign_stats_callback(spaceagora, rollout)
        apoapsis_callback = _spaceagora_physics_campaign_apoapsis_callback(spaceagora, rollout)
        ephemeris_callback =
            _spaceagora_rl_shared_ephemeris_callback(spaceagora, config, pass_cap)
        extra_callbacks = ephemeris_callback === nothing ?
                          (stats_callback, apoapsis_callback) :
                          (ephemeris_callback, stats_callback, apoapsis_callback)
        run_simulation_fn = Base.invokelatest(getproperty, spaceagora, :run_simulation)
        Base.invokelatest(
            run_simulation_fn,
            args;
            return_solution = false,
            isolate_state = false,
            extra_callbacks = extra_callbacks,
        )

        if !rollout.terminated && length(rollout.transitions) < pass_cap
            integrator = rollout.last_integrator
            integrator === nothing &&
                throw(ErrorException("SpaceAGORA physics campaign did not expose a final integrator state"))
            args = getproperty(getproperty(integrator, :p), :args)
            result = _spaceagora_physics_campaign_push_transition!(
                spaceagora,
                rollout,
                args,
                getproperty(integrator, :u),
                Float64(getproperty(integrator, :t));
                force_truncated = true,
            )
            if result !== nothing
                protected_transition = rollout.protected_next_transition
                rollout.protected_next_transition = false
                _spaceagora_physics_campaign_emit!(
                    rollout,
                    last(rollout.transitions),
                    result,
                    true;
                    protected=protected_transition,
                )
            end
        end
        return finalize_episode_summary(rollout.summary, config), rollout.transitions
    catch err
        bt = catch_backtrace()
        message = "SpaceAGORA physics campaign rollout failed while propagating a continuous aerobraking campaign. " *
                  "Original error:\n$(sprint(showerror, err, bt))"
        if rollout === nothing
            put!(event_channel, SpaceAGORAPhysicsPassEvent(
                worker_id = worker_id,
                episode_index = episode_index,
                seed = seed,
                summary = summary,
                done = true,
                error = message,
            ))
            return finalize_episode_summary(summary, config), Transition[]
        end
        _spaceagora_physics_campaign_emit!(rollout, nothing, nothing, true; error=message)
        return finalize_episode_summary(rollout.summary, config), rollout.transitions
    end
end

function start_spaceagora_physics_campaign_worker!(event_channel::Channel{Any},
                                                   config::AerobrakingScenarioConfig,
                                                   schedule::EpsilonSchedule,
                                                   ddqn_config::DDQNConfig,
                                                   policy_snapshot::Union{Nothing,QNetwork},
                                                   simulation_template::SpaceAGORAPhysicsSimulationTemplate,
                                                   state::AerobrakingDecisionState,
                                                   norm_obs::Vector{Float32},
                                                   action_index::Int,
                                                   summary::EpisodeSummary,
                                                   episode_index::Int,
                                                   worker_id::Int,
                                                   seed::Int,
                                                   max_passes_per_campaign::Int,
                                                   global_step_start::Int;
                                                   protected_first_pass::Bool=false,
                                                   protected_suppress_thermal_terminal::Bool=true)
    action_channel = Channel{Any}(1)
    task = Threads.@spawn run_spaceagora_physics_campaign_streaming_worker_episode(
        $event_channel,
        $action_channel,
        $config,
        $schedule,
        $ddqn_config,
        $policy_snapshot,
        $simulation_template,
        $state,
        $norm_obs,
        $action_index,
        $summary,
        $episode_index,
        $worker_id,
        $seed,
        $max_passes_per_campaign,
        $global_step_start;
        protected_first_pass = $protected_first_pass,
        protected_suppress_thermal_terminal = $protected_suppress_thermal_terminal,
    )
    return SpaceAGORAPhysicsWorkerHandle(task, action_channel)
end

function run_worker_episode!(session::TrainingSession, episode_index::Int, worker_id::Int, seed::Int;
                             train::Bool=true)
    rng = MersenneTwister(seed)
    config = session.config.scenario
    state = reset_scenario(session.backend, rng)
    obs = observe_state(config, state)
    norm_obs = normalize_observation(obs, config.normalization_bounds)
    summary = empty_episode_summary(episode_index=episode_index, worker_id=worker_id, seed=seed)
    transitions = Transition[]

    while length(transitions) < session.config.training.max_passes_per_campaign
        action_index = learner_policy_action!(session.learner, norm_obs, rng; test=!train)
        result = step_scenario(session.backend, state, action_index, rng)
        transition = transition_from_step(norm_obs, action_index, result, length(transitions) + 1)
        push!(transitions, transition)
        if train
            observe!(session.learner, transition)
            maybe_train!(session.learner, rng)
        end
        summary = update_episode_summary(summary, result)
        state = result.state
        obs = result.raw_observation
        norm_obs = result.normalized_observation
        if result.flags.terminated || result.flags.truncated
            break
        end
    end

    return finalize_episode_summary(summary, config), transitions
end

function run_threaded_worker_episode(config::AerobrakingScenarioConfig,
                                     schedule::EpsilonSchedule,
                                     ddqn_config::DDQNConfig,
                                     policy_snapshot::QNetwork,
                                     episode_index::Int,
                                     worker_id::Int,
                                     seed::Int,
                                     max_passes_per_campaign::Int,
                                     global_step_start::Int;
                                     train::Bool=true)
    if _is_spaceagora_live_backend(config.backend_mode)
        return run_spaceagora_physics_campaign_worker_episode(
            config,
            schedule,
            ddqn_config,
            policy_snapshot,
            episode_index,
            worker_id,
            seed,
            max_passes_per_campaign,
            global_step_start;
            train=train,
        )
    end

    rng = MersenneTwister(seed)
    state = reset_scenario(config, rng)
    obs = observe_state(config, state)
    norm_obs = normalize_observation(obs, config.normalization_bounds)
    summary = empty_episode_summary(episode_index=episode_index, worker_id=worker_id, seed=seed)
    transitions = Transition[]

    while length(transitions) < max_passes_per_campaign
        step = global_step_start + length(transitions) + 1
        action_index = snapshot_policy_action(policy_snapshot, schedule, ddqn_config, norm_obs, step, rng; test=!train)
        result = step_scenario(config, state, action_index, rng)
        transition = transition_from_step(norm_obs, action_index, result, length(transitions) + 1)
        push!(transitions, transition)
        summary = update_episode_summary(summary, result)
        state = result.state
        obs = result.raw_observation
        norm_obs = result.normalized_observation
        if result.flags.terminated || result.flags.truncated
            break
        end
    end

    return finalize_episode_summary(summary, config), transitions
end
