struct RunManifest
    run_id::String
    created_utc::DateTime
    phase::String
    seed::Int
    config_path::Union{Nothing,String}
    config_sha256::Union{Nothing,String}
    backend_mode::Symbol
    action_table::Vector{Float64}
    normalization_names::Vector{Symbol}
    reward_config::RewardConfig
    ddqn_config::DDQNConfig
end

function config_sha256(path::Union{Nothing,String})
    path === nothing && return nothing
    isfile(path) || return nothing
    return bytes2hex(sha256(read(path)))
end

function RunManifest(config::ResolvedConfig; run_id::Union{Nothing,String}=nothing)
    rid = run_id === nothing ? Dates.format(now(UTC), dateformat"yyyymmddTHHMMSS") : run_id
    return RunManifest(
        rid,
        now(UTC),
        config.scenario.phase,
        config.training.seed,
        config.source_path,
        config_sha256(config.source_path),
        config.scenario.backend_mode,
        copy(PAPER_ACTIONS_MPS),
        config.scenario.normalization_bounds.names,
        config.scenario.reward_config,
        config.ddqn,
    )
end

function manifest_dict(manifest::RunManifest)
    return Dict{String,Any}(
        "run_id" => manifest.run_id,
        "created_utc" => string(manifest.created_utc),
        "phase" => manifest.phase,
        "seed" => manifest.seed,
        "config_path" => manifest.config_path === nothing ? "" : manifest.config_path,
        "config_sha256" => manifest.config_sha256 === nothing ? "" : manifest.config_sha256,
        "backend_mode" => string(manifest.backend_mode),
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
    )
end

function write_manifest(path::AbstractString, manifest::RunManifest)
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, manifest_dict(manifest))
    end
    return path
end
