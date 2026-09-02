Base.@kwdef struct TrainingConfig
    seed::Int = 42
    algorithm::Symbol = :pr_drl
    device::Symbol = :cpu
    global_steps::Int = 0
    episodes::Int = 4
    max_passes_per_campaign::Int = 1_000
    n_workers::Int = 1
    worker_backend::Symbol = :threads
    checkpoint_frequency::Int = 500
    progress_frequency::Int = 50
    validate_checkpoints::Bool = false
    validation_episodes::Int = PAPER_IID_EVALUATION_EPISODES
    validation_seed::Int = 1
    validation_checkpoint_stride::Int = 1
    successful_case_repetitions::Int = 0
    output_dir::String = joinpath(package_root(), "outputs", "runs")
    protected_first_pass::Bool = false
    protected_initial_corridor_maneuver::Bool = false
    protected_first_pass_suppress_thermal_terminal::Bool = true
    protected_corridor_low_w_cm2::Float64 = 0.025
    protected_corridor_high_w_cm2::Float64 = 0.40
end

const DDQN_FAMILY_ALGORITHMS = (:ddqn, :pr_drl)
const ACTOR_CRITIC_ALGORITHMS = (:a2c, :a3c)
const CONTINUOUS_OFF_POLICY_ALGORITHMS = (:td3,)

canonical_algorithm(value) = Symbol(replace(lowercase(String(value)), "-" => "_"))
canonical_spaceagora_atmosphere_model(value) = Symbol(replace(lowercase(strip(String(value))), "-" => "_"))
canonical_worker_backend(value) = Symbol(replace(lowercase(strip(String(value))), "-" => "_"))
is_ddqn_family_algorithm(algorithm::Symbol) = algorithm in DDQN_FAMILY_ALGORITHMS
is_actor_critic_algorithm(algorithm::Symbol) = algorithm in ACTOR_CRITIC_ALGORITHMS
is_continuous_off_policy_algorithm(algorithm::Symbol) =
    algorithm in CONTINUOUS_OFF_POLICY_ALGORITHMS
algorithm_report_name(algorithm::Symbol) = algorithm == :pr_drl ? "pr_drl" : string(algorithm)
algorithm_display_name(algorithm::Symbol) = algorithm == :pr_drl ? "PR-DRL" : uppercase(string(algorithm))

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
    td3::TD3Config
    a2c::A2CConfig
    a3c::A3CConfig
    actor_critic_action::ActorCriticActionConfig
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

function _optional_date_from_table(table::Dict, key::String)
    haskey(table, key) || return nothing
    value = table[key]
    value === nothing && return nothing
    value isa Date && return value
    value isa DateTime && return Date(value)
    if value isa AbstractString
        parsed = tryparse(Date, strip(String(value)))
        parsed === nothing && throw(ArgumentError("config value $key must be an ISO date, got $(repr(value))"))
        return parsed
    end
    throw(ArgumentError("config value $key must be an ISO date string, got $(typeof(value))"))
end

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

function resolve_config(raw::Dict{String,Any};
                        source_path::Union{Nothing,String}=nothing,
                        gram_wind_mode=nothing)
    scenario_table = _table(raw, "scenario")
    reward_table = _table(raw, "reward")
    term_table = _table(raw, "termination")
    random_table = _table(raw, "randomization")
    ddqn_table = _table(raw, "ddqn")
    td3_table = _table(raw, "td3")
    a2c_table = _table(raw, "a2c")
    a3c_table = _table(raw, "a3c")
    eps_table = _table(raw, "epsilon")
    train_table = _table(raw, "training")
    report_table = _table(raw, "reports")
    spaceagora_table = _table(raw, "spaceagora_physics")
    resolved_gram_wind_mode = canonical_gram_wind_mode(
        gram_wind_mode === nothing ?
        _get(spaceagora_table, "gram_wind_mode", "perturbed") :
        gram_wind_mode,
    )

    phase = String(_get(scenario_table, "phase", "Main"))
    reward = reward_config_from_table(reward_table)
    termination = TerminationConfig(
        impact_periapsis_altitude_m = Float64(_get(term_table, "impact_periapsis_altitude_m", 85e3)),
        out_of_passage_periapsis_altitude_m = Float64(_get(term_table, "out_of_passage_periapsis_altitude_m", 135e3)),
        max_passes = Int(_get(term_table, "max_passes", 1000)),
        terminal_on_thermal_violation = Bool(_get(term_table, "terminal_on_thermal_violation", true)),
    )
    randomization = AerobrakingRandomizationConfig(
        nominal = Bool(_get(random_table, "nominal", true)),
        apoapsis_jitter_m = Float64(_get(random_table, "apoapsis_jitter_m", 2.5e3)),
        periapsis_jitter_m = Float64(_get(random_table, "periapsis_jitter_m", 2.5e3)),
        angle_jitter_deg = Float64(_get(random_table, "angle_jitter_deg", 0.25)),
        nonnominal_inclination_low_deg = Float64(_get(random_table, "nonnominal_inclination_low_deg", 88.6)),
        nonnominal_inclination_high_deg = Float64(_get(random_table, "nonnominal_inclination_high_deg", 98.6)),
        nonnominal_aop_low_deg = Float64(_get(random_table, "nonnominal_aop_low_deg", 60.0)),
        nonnominal_aop_high_deg = Float64(_get(random_table, "nonnominal_aop_high_deg", 90.0)),
        nonnominal_raan_low_deg = Float64(_get(random_table, "nonnominal_raan_low_deg",
                                                _get(random_table, "nonnominal_angle_low_deg", 110.0))),
        nonnominal_raan_high_deg = Float64(_get(random_table, "nonnominal_raan_high_deg",
                                                 _get(random_table, "nonnominal_angle_high_deg", 120.0))),
        initial_date_start = _optional_date_from_table(random_table, "initial_date_start"),
        initial_date_days = Int(_get(random_table, "initial_date_days", 0)),
        randomize_initial_time_of_day = Bool(_get(random_table, "randomize_initial_time_of_day", true)),
        initial_true_anomaly_jitter_deg = Float64(_get(random_table, "initial_true_anomaly_jitter_deg", 0.0)),
        process_noise = Bool(_get(random_table, "process_noise", false)),
        process_noise_scale = Float64(_get(random_table, "process_noise_scale", 0.0)),
        aerodynamic_coefficient_dispersion = Bool(_get(random_table, "aerodynamic_coefficient_dispersion", false)),
        aerodynamic_coefficient_span = Float64(_get(random_table, "aerodynamic_coefficient_span", 0.10)),
        aerodynamic_cd_span = Float64(_get(random_table, "aerodynamic_cd_span",
                                           _get(random_table, "aerodynamic_coefficient_span", 0.10))),
        aerodynamic_cl_span = Float64(_get(random_table, "aerodynamic_cl_span",
                                           _get(random_table, "aerodynamic_coefficient_span", 0.10))),
        marsgram_perturbation_scale = Float64(_get(random_table, "marsgram_perturbation_scale", 1.0)),
        marsgram_seed_base = Int(_get(random_table, "marsgram_seed_base", 1001)),
    )
    scenario = default_aerobraking_config(
        phase = phase,
        nominal = randomization.nominal,
        max_passes = termination.max_passes,
        backend_mode = Symbol(_get(scenario_table, "backend_mode", "paper_surrogate")),
        training = Bool(_get(scenario_table, "training", true)),
        spaceagora_atmosphere_model = canonical_spaceagora_atmosphere_model(
            _get(spaceagora_table, "atmosphere_model", "gram")
        ),
        spaceagora_gram_wind_mode = resolved_gram_wind_mode,
        spaceagora_gram_once_per_step = Bool(_get(spaceagora_table, "gram_once_per_step", false)),
        spaceagora_tabulated_flight_file = String(_get(
            spaceagora_table,
            "tabulated_flight_file",
            "data/telemetry/Odyssey/odyssey_accelerometer_density.feather",
        )),
        spaceagora_tabulated_flight_sigma = Float64(_get(spaceagora_table, "tabulated_flight_sigma", 0.0)),
        spaceagora_gravity_harmonics_degree = Int(_get(spaceagora_table, "gravity_harmonics_degree", 50)),
        spaceagora_gravity_harmonics_order = Int(_get(spaceagora_table, "gravity_harmonics_order", 50)),
        spaceagora_gravity_harmonics_file = String(_get(
            spaceagora_table,
            "gravity_harmonics_file",
            "data/Gravity_harmonics_data/Mars50c.csv",
        )),
        reward_config = reward,
        termination_config = termination,
        randomization_config = randomization,
    )
    scenario.spaceagora_atmosphere_model in
        (:gram, :marsgram, :marsgram_surrogate, :tabulated_flight, :exponential) ||
        throw(ArgumentError(
            "spaceagora_physics.atmosphere_model must be \"gram\", \"marsgram\", " *
            "\"marsgram_surrogate\", \"tabulated_flight\", or \"exponential\""
        ))
    if scenario.spaceagora_gram_once_per_step &&
       !(scenario.spaceagora_atmosphere_model in (:gram, :marsgram))
        throw(ArgumentError(
            "spaceagora_physics.gram_once_per_step is only valid with native \"gram\" or \"marsgram\""
        ))
    end

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
    td3 = TD3Config(
        actor_learning_rate = Float64(_get(td3_table, "actor_learning_rate", 3e-4)),
        critic_learning_rate = Float64(_get(td3_table, "critic_learning_rate", 3e-4)),
        discount = Float64(_get(td3_table, "discount", 0.95)),
        batch_size = Int(_get(td3_table, "batch_size", 256)),
        train_frequency = Int(_get(td3_table, "train_frequency", 1)),
        updates_per_step = Int(_get(td3_table, "updates_per_step", 1)),
        train_start = Int(_get(td3_table, "train_start", 25_000)),
        random_steps = Int(_get(td3_table, "random_steps", 25_000)),
        replay_size = Int(_get(td3_table, "replay_size", 1_000_000)),
        exploration_noise = Float64(_get(td3_table, "exploration_noise", 0.1)),
        target_policy_noise = Float64(_get(td3_table, "target_policy_noise", 0.2)),
        target_noise_clip = Float64(_get(td3_table, "target_noise_clip", 0.5)),
        policy_delay = Int(_get(td3_table, "policy_delay", 2)),
        tau = Float64(_get(td3_table, "tau", 0.005)),
        gradient_clip_norm = Float64(_get(td3_table, "gradient_clip_norm", 10.0)),
        adam_epsilon = Float64(_get(td3_table, "adam_epsilon", 1e-6)),
        adam_beta1 = Float64(_get(td3_table, "adam_beta1", 0.9)),
        adam_beta2 = Float64(_get(td3_table, "adam_beta2", 0.999)),
        hidden_dim = Int(_get(td3_table, "hidden_dim", 1024)),
        bootstrap_truncated = Bool(_get(td3_table, "bootstrap_truncated", true)),
    )
    a2c = A2CConfig(
        learning_rate = Float64(_get(a2c_table, "learning_rate", 1e-4)),
        discount = Float64(_get(a2c_table, "discount", 0.95)),
        segment_length = Int(_get(a2c_table, "segment_length", 10)),
        train_start = Int(_get(a2c_table, "train_start", 0)),
        entropy_coef = Float64(_get(a2c_table, "entropy_coef", 0.1)),
        value_coef = Float64(_get(a2c_table, "value_coef", 0.5)),
        normalize_advantages = Bool(_get(a2c_table, "normalize_advantages", true)),
        gradient_clip_norm = Float64(_get(a2c_table, "gradient_clip_norm", 0.5)),
        adam_epsilon = Float64(_get(a2c_table, "adam_epsilon", 1e-6)),
        hidden_dim = Int(_get(a2c_table, "hidden_dim", 1024)),
    )
    a3c = A3CConfig(
        learning_rate = Float64(_get(a3c_table, "learning_rate", 1e-4)),
        discount = Float64(_get(a3c_table, "discount", 0.95)),
        t_max = Int(_get(a3c_table, "t_max", 10)),
        entropy_coef = Float64(_get(a3c_table, "entropy_coef", 0.01)),
        value_coef = Float64(_get(a3c_table, "value_coef", 0.5)),
        normalize_advantages = Bool(_get(a3c_table, "normalize_advantages", false)),
        gradient_clip_norm = Float64(_get(a3c_table, "gradient_clip_norm", 0.5)),
        adam_epsilon = Float64(_get(a3c_table, "adam_epsilon", 1e-6)),
        adam_beta1 = Float64(_get(a3c_table, "adam_beta1", 0.9)),
        adam_beta2 = Float64(_get(a3c_table, "adam_beta2", 0.99)),
        hidden_dim = Int(_get(a3c_table, "hidden_dim", 1024)),
        max_policy_lag = Int(_get(a3c_table, "max_policy_lag", -1)),
    )
    resolved_algorithm = canonical_algorithm(_get(train_table, "algorithm", "pr_drl"))
    actor_critic_table = resolved_algorithm == :a3c ? a3c_table : a2c_table
    actor_critic_action = ActorCriticActionConfig(
        mode = canonical_actor_critic_action_mode(
            _get(actor_critic_table, "action_space", "discrete"),
        ),
        initial_log_std = Float64(_get(actor_critic_table, "initial_log_std", -1.0)),
        log_std_min = Float64(_get(actor_critic_table, "log_std_min", -5.0)),
        log_std_max = Float64(_get(actor_critic_table, "log_std_max", 1.0)),
    )
    epsilon = EpsilonSchedule(
        start = Float64(_get(eps_table, "start", 1.0)),
        stop = Float64(_get(eps_table, "stop", 0.0)),
        decay_steps = Int(_get(eps_table, "decay_steps", 500_000)),
        decay_start_step = Int(_get(eps_table, "decay_start_step", ddqn.train_start)),
    )
    training = TrainingConfig(
        seed = Int(_get(train_table, "seed", 42)),
        algorithm = resolved_algorithm,
        device = Symbol(lowercase(String(_get(train_table, "device", "cpu")))),
        global_steps = Int(_get(train_table, "global_steps", 0)),
        episodes = Int(_get(train_table, "episodes", 4)),
        max_passes_per_campaign = Int(_get(train_table, "max_passes_per_campaign",
                                           _get(train_table, "max_steps", 1_000))),
        n_workers = Int(_get(train_table, "n_workers", 1)),
        worker_backend = canonical_worker_backend(_get(train_table, "worker_backend", "threads")),
        checkpoint_frequency = Int(_get(train_table, "checkpoint_frequency", 500)),
        progress_frequency = Int(_get(train_table, "progress_frequency", 50)),
        validate_checkpoints = Bool(_get(train_table, "validate_checkpoints", false)),
        validation_episodes = Int(_get(
            train_table,
            "validation_episodes",
            PAPER_IID_EVALUATION_EPISODES,
        )),
        validation_seed = Int(_get(train_table, "validation_seed", 1)),
        validation_checkpoint_stride = Int(_get(
            train_table,
            "validation_checkpoint_stride",
            1,
        )),
        successful_case_repetitions = Int(_get(
            train_table,
            "successful_case_repetitions",
            0,
        )),
        output_dir = String(_get(train_table, "output_dir", joinpath(package_root(), "outputs", "runs"))),
        protected_first_pass = Bool(_get(train_table, "protected_first_pass", false)),
        protected_initial_corridor_maneuver = Bool(_get(train_table, "protected_initial_corridor_maneuver", false)),
        protected_first_pass_suppress_thermal_terminal = Bool(_get(train_table, "protected_first_pass_suppress_thermal_terminal", true)),
        protected_corridor_low_w_cm2 = Float64(_get(train_table, "protected_corridor_low_w_cm2", 0.025)),
        protected_corridor_high_w_cm2 = Float64(_get(train_table, "protected_corridor_high_w_cm2", 0.40)),
    )
    reports = ReportConfig(
        output_dir = String(_get(report_table, "output_dir", joinpath(package_root(), "outputs", "reports"))),
        write_csv = Bool(_get(report_table, "write_csv", true)),
        write_plots = Bool(_get(report_table, "write_plots", true)),
    )
    (is_ddqn_family_algorithm(training.algorithm) ||
     is_actor_critic_algorithm(training.algorithm) ||
     is_continuous_off_policy_algorithm(training.algorithm)) ||
        throw(ArgumentError(
            "training.algorithm must be \"pr_drl\", \"ddqn\", \"td3\", \"a2c\", or \"a3c\"",
        ))
    training.worker_backend in (:threads, :processes) ||
        throw(ArgumentError("training.worker_backend must be \"threads\" or \"processes\""))
    training.validation_checkpoint_stride > 0 ||
        throw(ArgumentError("training.validation_checkpoint_stride must be positive"))
    training.successful_case_repetitions >= 0 ||
        throw(ArgumentError("training.successful_case_repetitions must be nonnegative"))
    if is_actor_critic_algorithm(training.algorithm)
        actor_critic_action.mode in ACTOR_CRITIC_ACTION_MODES || throw(ArgumentError(
            "$(training.algorithm).action_space must be \"discrete\" or \"continuous\"",
        ))
        actor_critic_action.log_std_min <= actor_critic_action.initial_log_std <=
            actor_critic_action.log_std_max || throw(ArgumentError(
                "$(training.algorithm).initial_log_std must be between log_std_min and log_std_max",
            ))
    end
    if training.algorithm == :a2c
        a2c.learning_rate > 0 || throw(ArgumentError("a2c.learning_rate must be positive"))
        0 <= a2c.discount <= 1 ||
            throw(ArgumentError("a2c.discount must be between 0 and 1"))
        a2c.segment_length > 0 || throw(ArgumentError("a2c.segment_length must be positive"))
        a2c.train_start >= 0 || throw(ArgumentError("a2c.train_start must be nonnegative"))
        a2c.entropy_coef >= 0 || throw(ArgumentError("a2c.entropy_coef must be nonnegative"))
        a2c.value_coef >= 0 || throw(ArgumentError("a2c.value_coef must be nonnegative"))
        a2c.gradient_clip_norm > 0 ||
            throw(ArgumentError("a2c.gradient_clip_norm must be positive"))
    end
    if training.algorithm == :a3c
        a3c.learning_rate > 0 || throw(ArgumentError("a3c.learning_rate must be positive"))
        0 <= a3c.discount <= 1 ||
            throw(ArgumentError("a3c.discount must be between 0 and 1"))
        a3c.t_max > 0 || throw(ArgumentError("a3c.t_max must be positive"))
        a3c.entropy_coef >= 0 || throw(ArgumentError("a3c.entropy_coef must be nonnegative"))
        a3c.value_coef >= 0 || throw(ArgumentError("a3c.value_coef must be nonnegative"))
        a3c.gradient_clip_norm > 0 ||
            throw(ArgumentError("a3c.gradient_clip_norm must be positive"))
        0 <= a3c.adam_beta1 < 1 ||
            throw(ArgumentError("a3c.adam_beta1 must be in [0, 1)"))
        0 <= a3c.adam_beta2 < 1 ||
            throw(ArgumentError("a3c.adam_beta2 must be in [0, 1)"))
        a3c.max_policy_lag >= -1 ||
            throw(ArgumentError("a3c.max_policy_lag must be -1 or nonnegative"))
    end
    if training.algorithm == :td3
        td3.actor_learning_rate > 0 ||
            throw(ArgumentError("td3.actor_learning_rate must be positive"))
        td3.critic_learning_rate > 0 ||
            throw(ArgumentError("td3.critic_learning_rate must be positive"))
        0 <= td3.discount <= 1 ||
            throw(ArgumentError("td3.discount must be between 0 and 1"))
        td3.batch_size > 0 || throw(ArgumentError("td3.batch_size must be positive"))
        td3.train_frequency > 0 ||
            throw(ArgumentError("td3.train_frequency must be positive"))
        td3.updates_per_step > 0 ||
            throw(ArgumentError("td3.updates_per_step must be positive"))
        td3.train_start >= 0 || throw(ArgumentError("td3.train_start must be nonnegative"))
        td3.random_steps >= 0 || throw(ArgumentError("td3.random_steps must be nonnegative"))
        td3.replay_size >= td3.batch_size || throw(ArgumentError(
            "td3.replay_size must be at least td3.batch_size",
        ))
        td3.exploration_noise >= 0 ||
            throw(ArgumentError("td3.exploration_noise must be nonnegative"))
        td3.target_policy_noise >= 0 ||
            throw(ArgumentError("td3.target_policy_noise must be nonnegative"))
        td3.target_noise_clip >= 0 ||
            throw(ArgumentError("td3.target_noise_clip must be nonnegative"))
        td3.policy_delay > 0 || throw(ArgumentError("td3.policy_delay must be positive"))
        0 < td3.tau <= 1 || throw(ArgumentError("td3.tau must be in (0, 1]"))
        td3.gradient_clip_norm > 0 ||
            throw(ArgumentError("td3.gradient_clip_norm must be positive"))
        td3.adam_epsilon > 0 || throw(ArgumentError("td3.adam_epsilon must be positive"))
        0 <= td3.adam_beta1 < 1 ||
            throw(ArgumentError("td3.adam_beta1 must be in [0, 1)"))
        0 <= td3.adam_beta2 < 1 ||
            throw(ArgumentError("td3.adam_beta2 must be in [0, 1)"))
        td3.hidden_dim > 0 || throw(ArgumentError("td3.hidden_dim must be positive"))
        td3.action_dim == 1 || throw(ArgumentError(
            "aerobraking TD3 requires td3.action_dim=1",
        ))
    end
    return ResolvedConfig(
        source_path,
        raw,
        scenario,
        ddqn,
        td3,
        a2c,
        a3c,
        actor_critic_action,
        epsilon,
        training,
        reports,
    )
end

function resolve_config(path::AbstractString=default_config_path(); gram_wind_mode=nothing)
    return resolve_config(
        load_config(path);
        source_path=String(path),
        gram_wind_mode=gram_wind_mode,
    )
end
