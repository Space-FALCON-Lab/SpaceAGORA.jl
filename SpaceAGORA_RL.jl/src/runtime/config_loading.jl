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
    output_dir::String = joinpath(package_root(), "outputs", "runs")
    protected_first_pass::Bool = false
    protected_initial_corridor_maneuver::Bool = false
    protected_first_pass_suppress_thermal_terminal::Bool = true
    protected_corridor_low_w_cm2::Float64 = 0.025
    protected_corridor_high_w_cm2::Float64 = 0.40
end

const DDQN_FAMILY_ALGORITHMS = (:ddqn, :pr_drl)

canonical_algorithm(value) = Symbol(replace(lowercase(String(value)), "-" => "_"))
canonical_spaceagora_atmosphere_model(value) = Symbol(replace(lowercase(strip(String(value))), "-" => "_"))
canonical_worker_backend(value) = Symbol(replace(lowercase(strip(String(value))), "-" => "_"))
is_ddqn_family_algorithm(algorithm::Symbol) = algorithm in DDQN_FAMILY_ALGORITHMS
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
    a2c::A2CConfig
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

function resolve_config(raw::Dict{String,Any}; source_path::Union{Nothing,String}=nothing)
    scenario_table = _table(raw, "scenario")
    reward_table = _table(raw, "reward")
    term_table = _table(raw, "termination")
    random_table = _table(raw, "randomization")
    ddqn_table = _table(raw, "ddqn")
    a2c_table = _table(raw, "a2c")
    eps_table = _table(raw, "epsilon")
    train_table = _table(raw, "training")
    report_table = _table(raw, "reports")
    spaceagora_table = _table(raw, "spaceagora_physics")

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
    scenario.spaceagora_atmosphere_model in (:gram, :marsgram, :marsgram_surrogate, :tabulated_flight) ||
        throw(ArgumentError(
            "spaceagora_physics.atmosphere_model must be \"gram\", \"marsgram\", " *
            "\"marsgram_surrogate\", or \"tabulated_flight\""
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
    a2c = A2CConfig(
        learning_rate = Float64(_get(a2c_table, "learning_rate", 1e-4)),
        discount = Float64(_get(a2c_table, "discount", 0.95)),
        segment_length = Int(_get(a2c_table, "segment_length", 10)),
        train_start = Int(_get(a2c_table, "train_start", 0)),
        entropy_coef = Float64(_get(a2c_table, "entropy_coef", 0.1)),
        value_coef = Float64(_get(a2c_table, "value_coef", 0.5)),
        gradient_clip_norm = Float64(_get(a2c_table, "gradient_clip_norm", 0.5)),
        adam_epsilon = Float64(_get(a2c_table, "adam_epsilon", 1e-6)),
        hidden_dim = Int(_get(a2c_table, "hidden_dim", 1024)),
    )
    epsilon = EpsilonSchedule(
        start = Float64(_get(eps_table, "start", 1.0)),
        stop = Float64(_get(eps_table, "stop", 0.01)),
        decay_steps = Int(_get(eps_table, "decay_steps", 500_000)),
        decay_start_step = Int(_get(eps_table, "decay_start_step", ddqn.train_start)),
    )
    training = TrainingConfig(
        seed = Int(_get(train_table, "seed", 42)),
        algorithm = canonical_algorithm(_get(train_table, "algorithm", "pr_drl")),
        device = Symbol(lowercase(String(_get(train_table, "device", "cpu")))),
        global_steps = Int(_get(train_table, "global_steps", 0)),
        episodes = Int(_get(train_table, "episodes", 4)),
        max_passes_per_campaign = Int(_get(train_table, "max_passes_per_campaign",
                                           _get(train_table, "max_steps", 1_000))),
        n_workers = Int(_get(train_table, "n_workers", 1)),
        worker_backend = canonical_worker_backend(_get(train_table, "worker_backend", "threads")),
        checkpoint_frequency = Int(_get(train_table, "checkpoint_frequency", 500)),
        progress_frequency = Int(_get(train_table, "progress_frequency", 50)),
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
    (is_ddqn_family_algorithm(training.algorithm) || training.algorithm == :a2c) ||
        throw(ArgumentError("training.algorithm must be \"pr_drl\", \"ddqn\", or \"a2c\""))
    training.worker_backend in (:threads, :processes) ||
        throw(ArgumentError("training.worker_backend must be \"threads\" or \"processes\""))
    return ResolvedConfig(source_path, raw, scenario, ddqn, a2c, epsilon, training, reports)
end

function resolve_config(path::AbstractString=default_config_path())
    return resolve_config(load_config(path); source_path=String(path))
end
