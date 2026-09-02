struct RunManifest
    run_id::String
    created_utc::DateTime
    phase::String
    seed::Int
    config_path::Union{Nothing,String}
    config_sha256::Union{Nothing,String}
    backend_mode::Symbol
    algorithm::Symbol
    device::Symbol
    action_table::Vector{Float64}
    normalization_names::Vector{Symbol}
    reward_config::RewardConfig
    ddqn_config::DDQNConfig
    td3_config::TD3Config
    a2c_config::A2CConfig
    a3c_config::A3CConfig
    actor_critic_action_config::ActorCriticActionConfig
end

function config_sha256(path::Union{Nothing,String})
    path === nothing && return nothing
    isfile(path) || return nothing
    return bytes2hex(sha256(read(path)))
end

function _run_config_name(config_path::Union{Nothing,String})
    config_path === nothing && return "inline_config"
    stem = splitext(basename(config_path))[1]
    return isempty(stem) ? "inline_config" : stem
end

function run_title(config::ResolvedConfig)
    return string(_run_config_name(config.source_path), "-", config.training.algorithm)
end

function run_title(manifest::RunManifest)
    return string(_run_config_name(manifest.config_path), "-", manifest.algorithm)
end

function _run_title_slug(title::AbstractString)
    slug = replace(lowercase(String(title)), r"[^a-z0-9_-]+" => "-")
    slug = replace(slug, r"^-+|-+$" => "")
    return isempty(slug) ? "run" : slug
end

function RunManifest(config::ResolvedConfig; run_id::Union{Nothing,String}=nothing)
    title = run_title(config)
    rid = run_id === nothing ?
          string(Dates.format(now(UTC), dateformat"yyyymmddTHHMMSS"), "_", _run_title_slug(title)) :
          run_id
    return RunManifest(
        rid,
        now(UTC),
        config.scenario.phase,
        config.training.seed,
        config.source_path,
        config_sha256(config.source_path),
        config.scenario.backend_mode,
        config.training.algorithm,
        config.training.device,
        copy(PAPER_ACTIONS_MPS),
        config.scenario.normalization_bounds.names,
        config.scenario.reward_config,
        config.ddqn,
        config.td3,
        config.a2c,
        config.a3c,
        config.actor_critic_action,
    )
end

function manifest_dict(manifest::RunManifest; gram_wind_mode=nothing)
    result = Dict{String,Any}(
        "run_id" => manifest.run_id,
        "title" => run_title(manifest),
        "created_utc" => string(manifest.created_utc),
        "phase" => manifest.phase,
        "seed" => manifest.seed,
        "config_path" => manifest.config_path === nothing ? "" : manifest.config_path,
        "config_sha256" => manifest.config_sha256 === nothing ? "" : manifest.config_sha256,
        "backend_mode" => string(manifest.backend_mode),
        "algorithm" => string(manifest.algorithm),
        "device" => string(manifest.device),
        "action_table_mps" => manifest.action_table,
        "normalization_names" => string.(manifest.normalization_names),
        "reward" => Dict(
            "heat_low_w_cm2" => manifest.reward_config.heat_low_w_cm2,
            "heat_high_w_cm2" => manifest.reward_config.heat_high_w_cm2,
            "heat_medium_w_cm2" => manifest.reward_config.heat_medium_w_cm2,
            "heat_hard_w_cm2" => manifest.reward_config.heat_hard_w_cm2,
        ),
        "ddqn" => Dict(
            "learning_rate" => manifest.ddqn_config.learning_rate,
            "discount" => manifest.ddqn_config.discount,
            "batch_size" => manifest.ddqn_config.batch_size,
            "train_frequency" => manifest.ddqn_config.train_frequency,
            "train_start" => manifest.ddqn_config.train_start,
            "target_update" => manifest.ddqn_config.target_update,
            "hidden_dim" => manifest.ddqn_config.hidden_dim,
        ),
        "a2c" => Dict(
            "learning_rate" => manifest.a2c_config.learning_rate,
            "discount" => manifest.a2c_config.discount,
            "segment_length" => manifest.a2c_config.segment_length,
            "train_start" => manifest.a2c_config.train_start,
            "entropy_coef" => manifest.a2c_config.entropy_coef,
            "value_coef" => manifest.a2c_config.value_coef,
            "normalize_advantages" => manifest.a2c_config.normalize_advantages,
            "gradient_clip_norm" => manifest.a2c_config.gradient_clip_norm,
            "hidden_dim" => manifest.a2c_config.hidden_dim,
        ),
    )
    if manifest.algorithm == :a2c
        result["a2c"] = Dict{String,Any}(result["a2c"])
        result["a2c"]["action_space"] = string(manifest.actor_critic_action_config.mode)
    end
    if manifest.algorithm == :td3
        result["td3"] = Dict(
            "actor_learning_rate" => manifest.td3_config.actor_learning_rate,
            "critic_learning_rate" => manifest.td3_config.critic_learning_rate,
            "discount" => manifest.td3_config.discount,
            "batch_size" => manifest.td3_config.batch_size,
            "train_frequency" => manifest.td3_config.train_frequency,
            "updates_per_step" => manifest.td3_config.updates_per_step,
            "train_start" => manifest.td3_config.train_start,
            "random_steps" => manifest.td3_config.random_steps,
            "replay_size" => manifest.td3_config.replay_size,
            "exploration_noise" => manifest.td3_config.exploration_noise,
            "target_policy_noise" => manifest.td3_config.target_policy_noise,
            "target_noise_clip" => manifest.td3_config.target_noise_clip,
            "policy_delay" => manifest.td3_config.policy_delay,
            "tau" => manifest.td3_config.tau,
            "gradient_clip_norm" => manifest.td3_config.gradient_clip_norm,
            "hidden_dim" => manifest.td3_config.hidden_dim,
            "bootstrap_truncated" => manifest.td3_config.bootstrap_truncated,
        )
    end
    if manifest.algorithm == :a3c
        result["a3c"] = Dict(
            "learning_rate" => manifest.a3c_config.learning_rate,
            "discount" => manifest.a3c_config.discount,
            "t_max" => manifest.a3c_config.t_max,
            "entropy_coef" => manifest.a3c_config.entropy_coef,
            "value_coef" => manifest.a3c_config.value_coef,
            "normalize_advantages" => manifest.a3c_config.normalize_advantages,
            "gradient_clip_norm" => manifest.a3c_config.gradient_clip_norm,
            "adam_beta1" => manifest.a3c_config.adam_beta1,
            "adam_beta2" => manifest.a3c_config.adam_beta2,
            "hidden_dim" => manifest.a3c_config.hidden_dim,
            "max_policy_lag" => manifest.a3c_config.max_policy_lag,
            "action_space" => string(manifest.actor_critic_action_config.mode),
        )
    end
    if manifest.algorithm == :td3
        result["continuous_action"] = Dict(
            "low_mps" => CONTINUOUS_ACTION_LOW_MPS,
            "high_mps" => CONTINUOUS_ACTION_HIGH_MPS,
        )
    elseif manifest.algorithm in (:a2c, :a3c) &&
           uses_continuous_actions(manifest.actor_critic_action_config)
        result["continuous_action"] = Dict(
            "low_mps" => CONTINUOUS_ACTION_LOW_MPS,
            "high_mps" => CONTINUOUS_ACTION_HIGH_MPS,
            "initial_log_std" => manifest.actor_critic_action_config.initial_log_std,
            "log_std_min" => manifest.actor_critic_action_config.log_std_min,
            "log_std_max" => manifest.actor_critic_action_config.log_std_max,
        )
    end
    if gram_wind_mode !== nothing
        result["gram_wind_mode"] = String(canonical_gram_wind_mode(gram_wind_mode))
    end
    return result
end

function write_manifest(path::AbstractString, manifest::RunManifest; gram_wind_mode=nothing)
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, manifest_dict(manifest; gram_wind_mode=gram_wind_mode))
    end
    return path
end
