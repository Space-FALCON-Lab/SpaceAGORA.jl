Base.@kwdef struct TrainingConfig
    seed::Int = 42
    episodes::Int = 4
    max_steps::Int = 1_000
    n_workers::Int = 1
    checkpoint_frequency::Int = 500
    progress_frequency::Int = 50
    output_dir::String = joinpath(package_root(), "outputs", "runs")
end

Base.@kwdef struct ReportConfig
    output_dir::String = joinpath(package_root(), "outputs", "reports")
    write_csv::Bool = true
    write_plots::Bool = true
end

struct ResolvedConfig
    source_path::Union{Nothing,String}
    raw::Dict{String,Any}
    scenario::AerobrakingScenarioConfig
    ddqn::DDQNConfig
    epsilon::EpsilonSchedule
    training::TrainingConfig
    reports::ReportConfig
end

package_root() = dirname(dirname(@__DIR__))
default_config_path() = joinpath(package_root(), "configs", "aerobraking", "paper_replication.toml")

function _table(raw::Dict{String,Any}, name::String)
    value = get(raw, name, Dict{String,Any}())
    value isa Dict || throw(ArgumentError("config section [$name] must be a table"))
    return value
end

_get(table::Dict, key::String, default) = haskey(table, key) ? table[key] : default

function load_config(path::AbstractString=default_config_path())
    return TOML.parsefile(path)
end

function reward_config_from_table(table::Dict)
    return RewardConfig(
        heat_low_w_cm2 = Float64(_get(table, "heat_low_w_cm2", 0.05)),
        heat_high_w_cm2 = Float64(_get(table, "heat_high_w_cm2", 0.25)),
        heat_medium_w_cm2 = Float64(_get(table, "heat_medium_w_cm2", 0.30)),
        heat_hard_w_cm2 = Float64(_get(table, "heat_hard_w_cm2", 0.45)),
    )
end

function resolve_config(raw::Dict{String,Any}; source_path::Union{Nothing,String}=nothing)
    scenario_table = _table(raw, "scenario")
    reward_table = _table(raw, "reward")
    term_table = _table(raw, "termination")
    random_table = _table(raw, "randomization")
    ddqn_table = _table(raw, "ddqn")
    eps_table = _table(raw, "epsilon")
    train_table = _table(raw, "training")
    report_table = _table(raw, "reports")

    phase = String(_get(scenario_table, "phase", "Main"))
    reward = reward_config_from_table(reward_table)
    termination = TerminationConfig(
        impact_periapsis_altitude_m = Float64(_get(term_table, "impact_periapsis_altitude_m", 85e3)),
        out_of_passage_periapsis_altitude_m = Float64(_get(term_table, "out_of_passage_periapsis_altitude_m", 145e3)),
        max_passes = Int(_get(term_table, "max_passes", 80)),
        terminal_on_thermal_violation = Bool(_get(term_table, "terminal_on_thermal_violation", true)),
    )
    randomization = AerobrakingRandomizationConfig(
        nominal = Bool(_get(random_table, "nominal", true)),
        apoapsis_jitter_m = Float64(_get(random_table, "apoapsis_jitter_m", 2.5e3)),
        periapsis_jitter_m = Float64(_get(random_table, "periapsis_jitter_m", 1.0e3)),
        angle_jitter_deg = Float64(_get(random_table, "angle_jitter_deg", 0.25)),
        process_noise = Bool(_get(random_table, "process_noise", false)),
        process_noise_scale = Float64(_get(random_table, "process_noise_scale", 0.0)),
    )
    scenario = default_aerobraking_config(
        phase = phase,
        nominal = randomization.nominal,
        max_passes = termination.max_passes,
        backend_mode = Symbol(_get(scenario_table, "backend_mode", "paper_surrogate")),
        training = Bool(_get(scenario_table, "training", true)),
        reward_config = reward,
        termination_config = termination,
        randomization_config = randomization,
    )

    ddqn = DDQNConfig(
        learning_rate = Float64(_get(ddqn_table, "learning_rate", 1e-4)),
        discount = Float64(_get(ddqn_table, "discount", 0.95)),
        batch_size = Int(_get(ddqn_table, "batch_size", 256)),
        train_frequency = Int(_get(ddqn_table, "train_frequency", 5)),
        train_start = Int(_get(ddqn_table, "train_start", 10_000)),
        replay_size = Int(_get(ddqn_table, "replay_size", 1_000_000)),
        target_update = Int(_get(ddqn_table, "target_update", 10_000)),
        gradient_clip_norm = Float64(_get(ddqn_table, "gradient_clip_norm", 10.0)),
        adam_epsilon = Float64(_get(ddqn_table, "adam_epsilon", 1e-6)),
        hidden_dim = Int(_get(ddqn_table, "hidden_dim", 1024)),
    )
    epsilon = EpsilonSchedule(
        start = Float64(_get(eps_table, "start", 1.0)),
        stop = Float64(_get(eps_table, "stop", 0.01)),
        decay_steps = Int(_get(eps_table, "decay_steps", 500_000)),
        decay_start_step = Int(_get(eps_table, "decay_start_step", ddqn.train_start)),
    )
    training = TrainingConfig(
        seed = Int(_get(train_table, "seed", 42)),
        episodes = Int(_get(train_table, "episodes", 4)),
        max_steps = Int(_get(train_table, "max_steps", 1_000)),
        n_workers = Int(_get(train_table, "n_workers", 1)),
        checkpoint_frequency = Int(_get(train_table, "checkpoint_frequency", 500)),
        progress_frequency = Int(_get(train_table, "progress_frequency", 50)),
        output_dir = String(_get(train_table, "output_dir", joinpath(package_root(), "outputs", "runs"))),
    )
    reports = ReportConfig(
        output_dir = String(_get(report_table, "output_dir", joinpath(package_root(), "outputs", "reports"))),
        write_csv = Bool(_get(report_table, "write_csv", true)),
        write_plots = Bool(_get(report_table, "write_plots", true)),
    )
    return ResolvedConfig(source_path, raw, scenario, ddqn, epsilon, training, reports)
end

function resolve_config(path::AbstractString=default_config_path())
    return resolve_config(load_config(path); source_path=String(path))
end
