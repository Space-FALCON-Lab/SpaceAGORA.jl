using Test
using TOML
using DataFrames
using LinearAlgebra
using StaticArrays

const _COV_REPO_ROOT = isdefined(Main, :REPO_ROOT) ? Main.REPO_ROOT : normpath(joinpath(@__DIR__, "..", ".."))

if !isdefined(@__MODULE__, :SimulationModel)
    include(joinpath(_COV_REPO_ROOT, "src", "core", "simulation_model.jl"))
    using .SimulationModel
end
if !isdefined(@__MODULE__, :SimulationEngine)
    include(joinpath(_COV_REPO_ROOT, "src", "simulation", "engine", "simulation_engine.jl"))
end
if !isdefined(@__MODULE__, :TelemetryVerification)
    include(joinpath(_COV_REPO_ROOT, "src", "analysis", "verification", "telemetry_verification.jl"))
end

include(joinpath(_COV_REPO_ROOT, "src", "parallel", "routing", "parallel_profiles.jl"))

const PP = ParallelProfiles
const TV = TelemetryVerification

@testset "ParallelProfiles Coverage Probes" begin
    names_expected = Dict(
        PP.R0 => "R0",
        PP.R1_a => "R1_a",
        PP.R1_b => "R1_b",
        PP.R2 => "R2",
        PP.R3 => "R3",
        PP.R4 => "R4",
        PP.R5 => "R5"
    )
    for (profile, name_expected) in names_expected
        @test PP.parallel_profile_name(profile) == name_expected
        @test PP.parse_parallel_profile(name_expected) == profile
    end

    @test PP.parse_parallel_profile("outer_only_threads") == PP.R1_a
    @test PP.parse_parallel_profile("process_outer") == PP.R1_b
    @test PP.parse_parallel_profile("inner_only") == PP.R2
    @test PP.parse_parallel_profile("outer_inner_static") == PP.R3
    @test PP.parse_parallel_profile("auto_adaptive") == PP.R4
    @test PP.parse_parallel_profile("full_smart") == PP.R5
    @test PP.parse_parallel_profile("R4_full_auto") == PP.R5
    @test PP.parse_parallel_profile(:R0) == PP.R0
    @test PP.parse_parallel_profile(PP.R3) == PP.R3
    @test_throws ArgumentError PP.parse_parallel_profile("not_a_profile")

    cfg_r0 = PP.profile_config("R0")
    cfg_r3 = PP.profile_config("R3")
    cfg_r4 = PP.profile_config("R4")
    cfg_full = PP.profile_config("R5")
    @test cfg_r0.outer_backend == :none
    @test cfg_r3.outer_backend == :auto
    @test cfg_r4.inner_adaptive == true
    @test cfg_r4.thermal_mode == "auto"
    @test cfg_r4.inner_scheduler == "static"
    @test cfg_r4.adaptive_control_tail_guard == false
    @test cfg_r4.adaptive_measured_reward == false
    @test cfg_full.outer_route_adaptive == true
    @test cfg_full.label == "r5"
    # Both of these track measured reversals in profile_definitions.jl.
    # inner_scheduler went dynamic -> static in 622ae2a0 (the sweep now measures
    # the scheduler, so a profile constant must not pre-empt it) and this probe
    # was not updated with it. thermal_mode went on -> auto for the same reason:
    # R5 was the only profile forcing it, and "on" was never measured to beat
    # "auto" on any workload.
    @test cfg_full.thermal_mode == "auto"
    @test cfg_full.inner_scheduler == "static"
    @test cfg_full.adaptive_control_tail_guard == true
    @test cfg_full.adaptive_measured_reward == true
    @test cfg_full.adaptive_window == 4
    @test cfg_full.persistent_hints == true
    @test cfg_full.persistent_state_persist == true

    withenv("SPACEAGORA_PARALLEL_PROFILE" => "existing_profile") do
        env_pairs = PP.profile_env_pairs("R1_b"; preserve_existing=true, outer_parallel_active=true)
        env_map = Dict(env_pairs)
        @test env_map["SPACEAGORA_PARALLEL_PROFILE"] == "existing_profile"
        @test env_map["SPACEAGORA_OUTER_PARALLEL_ACTIVE"] == "1"
        @test env_map["SPACEAGORA_PERF_PARALLEL_BACKEND"] == "process"
    end
    env_pairs_override = PP.profile_env_pairs("R2"; preserve_existing=false, outer_parallel_active=false)
    env_map_override = Dict(env_pairs_override)
    @test env_map_override["SPACEAGORA_PARALLEL_PROFILE"] == "R2"
    @test env_map_override["SPACEAGORA_OUTER_PARALLEL_ACTIVE"] == "0"
    @test env_map_override["SPACEAGORA_PERF_PARALLEL_BACKEND"] == "none"
    @test env_map_override["SPACEAGORA_RHS_BATCH_PARALLEL"] == "auto"
    @test env_map_override["SPACEAGORA_PARALLEL_POLICY_WINDOW"] == "8"
    @test env_map_override["SPACEAGORA_PARALLEL_POLICY_CONTROL_TAIL_GUARD"] == "0"
    @test parse(Float64, env_map_override["SPACEAGORA_PARALLEL_POLICY_HINT_EXPLORATION"]) > 0.0
    @test parse(Int, env_map_override["SPACEAGORA_PARALLEL_POLICY_HINT_MIN_SAMPLES"]) >= 1
    env_pairs_auto = PP.profile_env_pairs("R5"; preserve_existing=false, outer_parallel_active=false)
    env_map_auto = Dict(env_pairs_auto)
    @test env_map_auto["SPACEAGORA_PERF_PARALLEL_BACKEND"] == "auto"
    @test env_map_auto["SPACEAGORA_THERMAL_CALLBACK_PARALLEL"] == "auto"
    @test env_map_auto["SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER"] == "static"
    @test env_map_auto["SPACEAGORA_PARALLEL_POLICY_WINDOW"] == "4"
    @test env_map_auto["SPACEAGORA_PARALLEL_POLICY_CONTROL_TAIL_GUARD"] == "1"
    @test env_map_auto["SPACEAGORA_PARALLEL_POLICY_MEASURED_REWARD"] == "1"
    @test env_map_auto["SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS"] == "1"
    @test env_map_auto["SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST"] == "1"
    @test env_map_auto["SPACEAGORA_RHS_BATCH_PARALLEL"] == "auto"
    @test parse(Float64, env_map_auto["SPACEAGORA_PARALLEL_POLICY_HINT_EXPLORATION"]) > 0.0
    @test parse(Int, env_map_auto["SPACEAGORA_PARALLEL_POLICY_HINT_MIN_SAMPLES"]) >= 2

    withenv("SPACEAGORA_RHS_BATCH_PARALLEL" => "on") do
        env_r0 = Dict(PP.profile_env_pairs("R0"; preserve_existing=true, outer_parallel_active=false))
        # R0 is forced true-serial and must not preserve a pre-existing RHS batch override.
        @test env_r0["SPACEAGORA_RHS_BATCH_PARALLEL"] == "off"
    end

    withenv("SPACEAGORA_PERF_HARDWARE_CLASS" => "small") do
        env_small = Dict(PP.profile_env_pairs("R5"; preserve_existing=false, outer_parallel_active=false))
        @test parse(Float64, env_small["SPACEAGORA_PARALLEL_POLICY_HINT_EXPLORATION"]) == 1.3
        @test parse(Int, env_small["SPACEAGORA_PARALLEL_POLICY_HINT_MIN_SAMPLES"]) == 2
    end
    withenv("SPACEAGORA_PERF_HARDWARE_CLASS" => "large") do
        env_large = Dict(PP.profile_env_pairs("R5"; preserve_existing=false, outer_parallel_active=false))
        @test parse(Float64, env_large["SPACEAGORA_PARALLEL_POLICY_HINT_EXPLORATION"]) == 1.8
        @test parse(Int, env_large["SPACEAGORA_PARALLEL_POLICY_HINT_MIN_SAMPLES"]) == 3
    end

    value_a = PP.with_parallel_profile("R1_a"; preserve_existing=false, outer_parallel_active=true) do
        (ENV["SPACEAGORA_PARALLEL_PROFILE"], ENV["SPACEAGORA_OUTER_PARALLEL_ACTIVE"])
    end
    @test value_a == ("R1_a", "1")
    value_b = PP.with_parallel_profile("R1_b", () -> ENV["SPACEAGORA_PERF_PARALLEL_BACKEND"]; preserve_existing=false)
    @test value_b == "process"

    f_light = PP.OuterRouteFeatures(category="deterministic", n_sats=1, n_links=1, mission_time_s=500.0)
    f_const = PP.OuterRouteFeatures(category="satellite_scaling", n_sats=6, n_links=8, mission_time_s=3_600.0, has_srp=true)
    f_heavy = PP.OuterRouteFeatures(category="deterministic", n_sats=3, n_links=6, mission_time_s=20_000.0, has_nbody=true, harmonics_degree=30)
    f_mc = PP.OuterRouteFeatures(category="montecarlo", n_sats=2, n_links=2, mission_time_s=6_000.0, montecarlo_samples=32)
    f_mc_small = PP.OuterRouteFeatures(category="montecarlo", n_sats=2, n_links=2, mission_time_s=500.0, montecarlo_samples=2)
    f_gram_point = PP.OuterRouteFeatures(
        category="deterministic",
        n_sats=1,
        n_links=1,
        mission_time_s=600.0,
        density_family="gram_point",
        gram_surrogate_enabled=false,
        gram_static_grid_enabled=false
    )
    f_gram_surrogate = PP.OuterRouteFeatures(
        category="deterministic",
        n_sats=1,
        n_links=1,
        mission_time_s=600.0,
        density_family="gram_surrogate",
        gram_surrogate_enabled=true,
        gram_static_grid_enabled=false
    )

    sig = PP.outer_route_signature(f_heavy)
    @test occursin("cat=deterministic", sig)
    @test occursin("harm=21p", sig)
    @test occursin("dens=", sig)
    @test occursin("solver=", sig)
    @test occursin("thermal=", sig)
    @test occursin("eff_cost=", sig)
    @test occursin("sat=1", PP.outer_route_signature(PP.OuterRouteFeatures(n_sats=1, n_links=1, mission_time_s=100.0)))
    @test occursin("sat=2", PP.outer_route_signature(PP.OuterRouteFeatures(n_sats=2, n_links=2, mission_time_s=100.0)))
    @test occursin("sat=5p", PP.outer_route_signature(PP.OuterRouteFeatures(n_sats=7, n_links=10, mission_time_s=10_000.0)))

    tune = PP.OuterRouteTuning(
        adaptive_enabled=true,
        adaptive_min_samples=2,
        spice_constellation_min_sats=4,
        mc_process_min_samples=16,
        mc_process_min_mission_s=3_600.0
    )
    @test PP.default_outer_route(f_light; tuning=tune, machine_class=:small, threads_available=true, parallel_enabled=true) == :none
    @test PP.default_outer_route(f_const; tuning=tune, machine_class=:small, threads_available=true, parallel_enabled=true) == :threads
    @test PP.default_outer_route(f_heavy; tuning=tune, machine_class=:large, threads_available=true, parallel_enabled=true) == :process
    @test PP.default_outer_route(f_heavy; tuning=tune, machine_class=:small, threads_available=false, parallel_enabled=true) == :none
    @test PP.default_outer_route(f_mc; tuning=tune, machine_class=:large, threads_available=true, parallel_enabled=true) == :process
    @test PP.default_outer_route(f_mc_small; tuning=tune, machine_class=:small, threads_available=true, parallel_enabled=true) == :threads
    @test PP.default_outer_route(f_mc; tuning=tune, machine_class=:large, threads_available=true, parallel_enabled=false) == :none
    @test PP.default_outer_route(f_gram_point; tuning=tune, machine_class=:large, threads_available=true, parallel_enabled=true) == :process
    @test PP.default_outer_route(f_gram_point; tuning=tune, machine_class=:small, threads_available=true, parallel_enabled=true) == :process
    @test PP.default_outer_route(f_gram_surrogate; tuning=tune, machine_class=:large, threads_available=true, parallel_enabled=true) == :none

    cand_heavy = PP.outer_route_candidates(f_heavy; tuning=tune, machine_class=:large, threads_available=true, parallel_enabled=true)
    cand_light = PP.outer_route_candidates(f_light; tuning=tune, machine_class=:small, threads_available=false, parallel_enabled=true)
    cand_mc = PP.outer_route_candidates(f_mc; tuning=tune, machine_class=:large, threads_available=true, parallel_enabled=true)
    cand_gram = PP.outer_route_candidates(f_gram_point; tuning=tune, machine_class=:small, threads_available=true, parallel_enabled=true)
    cand_gram_sur = PP.outer_route_candidates(f_gram_surrogate; tuning=tune, machine_class=:small, threads_available=true, parallel_enabled=true)
    @test :process in cand_heavy
    @test cand_light == [:none]
    @test :process in cand_mc
    @test :process in cand_gram
    @test !(:threads in cand_gram)
    @test !(:process in cand_gram_sur)

    state = PP.OuterRouteState()
    @test isempty(PP.outer_route_stats_snapshot(state, sig))
    PP.record_outer_route_feedback!(state, f_heavy; route=:bad_route, successes=1, failures=0, tuning=tune)
    @test isempty(PP.outer_route_stats_snapshot(state, sig))

    PP.record_outer_route_feedback!(state, f_heavy; route=:threads, successes=2, failures=0, elapsed_success_s=100.0, tuning=tune)
    PP.record_outer_route_feedback!(state, f_heavy; route=:process, successes=2, failures=0, elapsed_success_s=10.0, tuning=tune)
    PP.record_outer_route_feedback!(state, f_heavy; route=:none, successes=2, failures=0, elapsed_success_s=200.0, tuning=tune)
    snap = PP.outer_route_stats_snapshot(state, sig)
    @test haskey(snap, :threads)
    @test haskey(snap, :process)
    @test haskey(snap, :none)
    @test snap[:threads].samples == 2
    @test snap[:process].success_rate == 1.0

    state_var = PP.OuterRouteState()
    PP.record_outer_route_feedback!(
        state_var,
        f_heavy;
        route=:threads,
        successes=2,
        failures=0,
        elapsed_success_s=10.0,
        elapsed_success_sq_sum_s=58.0,
        tuning=tune
    )
    snap_var = PP._outer_route_stats_snapshot_internal(state_var, sig)
    @test isapprox(snap_var[:threads].mean_s, 5.0; atol=1e-12, rtol=0.0)
    @test snap_var[:threads].std_s > 0.0

    chosen_default = PP.select_outer_route!(
        state,
        f_heavy;
        tuning=PP.OuterRouteTuning(adaptive_enabled=false),
        machine_class=:large,
        threads_available=true,
        parallel_enabled=true
    )
    @test chosen_default == :process

    chosen_adaptive = PP.select_outer_route!(
        state,
        f_heavy;
        tuning=tune,
        machine_class=:large,
        threads_available=true,
        parallel_enabled=true
    )
    @test chosen_adaptive == :process

    f_variant = PP.OuterRouteFeatures(
        category="deterministic",
        n_sats=3,
        n_links=6,
        max_links_per_sat=3,
        mission_time_s=20_000.0,
        has_nbody=true,
        harmonics_degree=30,
        has_control=false,
        orientation_on=false,
        density_family="gram_point",
        solver_mode="split_imex:kencarp4",
        dt_max_orbit_s=0.2,
        thermal_enabled=true,
        dynamic_effector_count=5,
        effector_cost_class="heavy"
    )
    sig_variant = PP.outer_route_signature(f_variant)
    @test isempty(PP.outer_route_stats_snapshot(state, sig_variant))
    @test PP.default_outer_route(f_variant; tuning=tune, machine_class=:small, threads_available=true, parallel_enabled=true) == :process
    chosen_fallback = PP.select_outer_route!(
        state,
        f_variant;
        tuning=tune,
        machine_class=:small,
        threads_available=true,
        parallel_enabled=true
    )
    @test chosen_fallback == :process

    # Drive unhit route-bucket and adaptive-selection branches.
    @test_throws ArgumentError PP._profile_outer_backend_token(:unsupported_backend)
    @test PP._route_max_link_bucket(2) == "2"
    @test PP._route_max_link_bucket(6) == "5_8"
    @test PP._route_harmonics_bucket(5) == "1_10"
    @test PP._route_harmonics_bucket(15) == "11_20"
    @test PP._route_density_bucket("vacuum") == "none"
    @test PP._route_density_bucket("exp") == "exp"
    @test PP._route_density_bucket("polyfit") == "poly"
    @test PP._route_density_bucket("nrlmsise00") == "nrl"
    @test PP._route_density_bucket("") == "unknown"
    @test PP._route_solver_bucket("rodas5p") == "rodas"
    @test PP._route_solver_bucket("tsit5") == "tsit5"
    @test PP._route_solver_bucket("multirate") == "mrate"
    @test PP._route_solver_bucket("") == "auto"
    @test PP._route_solver_bucket("custom|solver") == "custom_solver"
    @test PP._route_interval_bucket(0.5) == "fast"
    @test PP._route_interval_bucket(5.0) == "med"
    @test PP._route_interval_bucket(20.0) == "slow"
    @test PP._route_count_bucket(1) == "1"
    @test PP._route_count_bucket(3) == "2_3"
    @test PP._route_count_bucket(9) == "7p"
    @test PP._route_effector_cost_bucket("") == "unknown"

    sparse_state = PP.OuterRouteState()
    lock(sparse_state.lock) do
        sparse_state.history["empty_sig"] = Dict(:threads => PP.OuterRouteStats(samples=0, successes=0, failures=0, elapsed_sum_s=0.0))
    end
    @test isempty(PP.outer_route_stats_snapshot(sparse_state, "empty_sig"))

    @test PP._priority_outer_route_montecarlo(
        f_mc,
        tune;
        machine_class=:large,
        threads_available=true,
        parallel_enabled=false
    ) == :none
    @test PP._priority_outer_route_montecarlo(
        PP.OuterRouteFeatures(category="montecarlo", n_sats=1, n_links=1, mission_time_s=100.0, montecarlo_samples=1),
        tune;
        machine_class=:large,
        threads_available=true,
        parallel_enabled=true
    ) == :none

    f_inner_gate = PP.OuterRouteFeatures(category="deterministic", n_sats=3, n_links=1, mission_time_s=2_000.0)
    tune_inner_gate = PP.OuterRouteTuning(inner_sat_threshold=2, inner_link_threshold=12)
    @test PP.default_outer_route(f_inner_gate; tuning=tune_inner_gate, machine_class=:large, threads_available=true, parallel_enabled=true) == :none

    f_sat_scaling = PP.OuterRouteFeatures(category="satellite_scaling", n_sats=5, n_links=2, mission_time_s=3_000.0)
    tune_sat_scaling = PP.OuterRouteTuning(inner_sat_threshold=10, inner_link_threshold=12)
    @test PP.default_outer_route(f_sat_scaling; tuning=tune_sat_scaling, machine_class=:small, threads_available=true, parallel_enabled=true) == :threads

    f_machine_choice = PP.OuterRouteFeatures(category="deterministic", n_sats=1, n_links=2, mission_time_s=3_000.0, has_control=true)
    @test PP.default_outer_route(f_machine_choice; tuning=tune, machine_class=:large, threads_available=true, parallel_enabled=true) == :process
    @test PP.default_outer_route(f_machine_choice; tuning=tune, machine_class=:small, threads_available=true, parallel_enabled=true) == :threads

    @test PP._route_ranked_candidates([:threads, :custom], :none) == [:threads, :custom]
    under_sampled = PP._under_sampled_candidate(
        [:threads, :none],
        Dict(
            :threads => (samples=3, mean_s=10.0, success_rate=1.0),
            :none => (samples=0, mean_s=Inf, success_rate=0.0)
        ),
        :threads,
        2
    )
    @test under_sampled == :none
    best_valid = PP._best_candidate(
        [:threads, :none],
        Dict(
            :threads => (samples=2, mean_s=Inf, success_rate=1.0),
            :none => (samples=2, mean_s=3.0, success_rate=0.5)
        )
    )
    @test best_valid == :none
    best_tie = PP._best_candidate(
        [:threads, :process],
        Dict(
            :threads => (samples=2, mean_s=2.0, success_rate=0.1),
            :process => (samples=2, mean_s=2.0, success_rate=0.9)
        )
    )
    @test best_tie == :process
    best_conf = PP._best_candidate_confidence(
        [:threads, :process],
        Dict(
            :threads => (samples=8, mean_s=5.0, success_rate=1.0, std_s=0.05),
            :process => (samples=8, mean_s=5.0, success_rate=1.0, std_s=1.5)
        ),
        :threads,
        1.25
    )
    @test best_conf.route == :process
    @test best_conf.confidence_s >= 0.0
    @test best_conf.regret_s >= 0.0

    state_explore = PP.OuterRouteState()
    PP.record_outer_route_feedback!(state_explore, f_machine_choice; route=:threads, successes=2, failures=0, elapsed_success_s=5.0, tuning=tune)
    chosen_default_not_candidate = PP.select_outer_route!(
        state_explore,
        f_machine_choice;
        tuning=PP.OuterRouteTuning(adaptive_enabled=true, adaptive_min_samples=2, trace=true),
        machine_class=:large,
        threads_available=true,
        parallel_enabled=true
    )
    @test chosen_default_not_candidate == :none

    mktempdir() do tmp
        cache_path = joinpath(tmp, "outer_route_state.toml")
        saved = PP.save_outer_route_state(state, cache_path; metadata=Dict("profile" => "quick"))
        @test isfile(cache_path)
        @test saved.rows >= 3
        @test saved.signatures >= 1
        parsed_cache = TOML.parsefile(cache_path)
        rows_cache = get(parsed_cache, "history", Any[])
        @test any(row -> haskey(get(row, "stats", Dict{String, Any}()), "elapsed_sq_sum_s"), rows_cache)

        loaded_state = PP.OuterRouteState()
        loaded = PP.load_outer_route_state!(loaded_state, cache_path; replace=true)
        @test loaded.rows >= 3
        snap_loaded = PP.outer_route_stats_snapshot(loaded_state, sig)
        @test snap_loaded[:process].samples == snap[:process].samples
        @test snap_loaded[:threads].samples == snap[:threads].samples

        PP.record_outer_route_feedback!(loaded_state, f_heavy; route=:process, successes=1, failures=0, elapsed_success_s=1.0, tuning=tune)
        merged = PP.load_outer_route_state!(loaded_state, cache_path; replace=false)
        @test merged.rows >= 3
        snap_merged = PP.outer_route_stats_snapshot(loaded_state, sig)
        @test snap_merged[:process].samples == 2 * snap[:process].samples + 1
    end

    PP.reset_outer_route_state!(state)
    @test isempty(PP.outer_route_stats_snapshot(state, sig))
end

@testset "TelemetryVerification Coverage Probes" begin
    @test TV._safe_parse_bool("on", false) == true
    @test TV._safe_parse_bool("off", true) == false
    @test TV._safe_parse_bool("unknown", true) == true
    @test TV._safe_parse_bool("unknown", false) == false

    withenv("TV_INT_ENV" => "") do
        @test TV._parse_positive_int_env("TV_INT_ENV", 7) == 7
    end
    withenv("TV_INT_ENV" => "11") do
        @test TV._parse_positive_int_env("TV_INT_ENV", 7) == 11
    end
    withenv("TV_INT_ENV" => "0") do
        @test_throws ArgumentError TV._parse_positive_int_env("TV_INT_ENV", 7)
    end
    withenv("TV_INT_ENV" => "bad") do
        @test_throws ArgumentError TV._parse_positive_int_env("TV_INT_ENV", 7)
    end

    withenv("TV_FLOAT_ENV" => "") do
        @test TV._parse_positive_float_env("TV_FLOAT_ENV") === nothing
    end
    withenv("TV_FLOAT_ENV" => "1.25") do
        @test TV._parse_positive_float_env("TV_FLOAT_ENV") == 1.25
    end
    withenv("TV_FLOAT_ENV" => "0.0") do
        @test_throws ArgumentError TV._parse_positive_float_env("TV_FLOAT_ENV")
    end
    withenv("TV_FLOAT_ENV" => "bad") do
        @test_throws ArgumentError TV._parse_positive_float_env("TV_FLOAT_ENV")
    end

    withenv("SPACEAGORA_TELEMETRY_SOLVER_MAXITERS" => "6000000") do
        @test TV._telemetry_solver_maxiters(:quick) == 6_000_000
    end
    withenv("SPACEAGORA_TELEMETRY_SOLVER_MAXITERS" => "") do
        @test TV._telemetry_solver_maxiters(:full) == TV.TELEMETRY_SOLVER_MAXITERS_FULL_DEFAULT
    end
    withenv("SPACEAGORA_TELEMETRY_SOLVER_MAXITERS_RETRY" => "7000000") do
        @test TV._telemetry_solver_retry_maxiters(6_000_000) == 7_000_000
    end
    withenv("SPACEAGORA_TELEMETRY_SOLVER_MAXITERS_RETRY" => "1000") do
        @test TV._telemetry_solver_retry_maxiters(6_000_000) == 6_000_001
    end

    withenv("SPACEAGORA_TELEMETRY_SOLVER_MODE" => "multirate", "SPACEAGORA_SOLVER_MODE" => "tsit5") do
        @test TV._telemetry_solver_mode() == "multirate"
    end
    withenv("SPACEAGORA_TELEMETRY_SOLVER_MODE" => "", "SPACEAGORA_SOLVER_MODE" => "rodas5p") do
        @test TV._telemetry_solver_mode() == "rodas5p"
    end
    withenv("SPACEAGORA_TELEMETRY_SOLVER_MODE" => "", "SPACEAGORA_SOLVER_MODE" => "") do
        @test TV._telemetry_solver_mode() == "auto_stiff"
    end

    @test TV._is_maxiters_error(ErrorException("MaxIters reached in solve"))
    @test !TV._is_maxiters_error(ErrorException("other"))
    @test TV._is_gram_library_missing_error(ErrorException("GRAM shared library not found"))
    @test !TV._is_gram_library_missing_error(ErrorException("random"))

    @test TV._parse_orbit_altitude_mode("vacuum", "ctx") == :vacuum
    @test TV._parse_orbit_altitude_mode("oblate", "ctx") == :oblate
    @test_throws ArgumentError TV._parse_orbit_altitude_mode("bad", "ctx")

    @test TV._parse_time_aligned_comparison_mode("time_aligned_state", "ctx") == :time_aligned_state
    @test TV._parse_time_aligned_comparison_mode("orbit_events", "ctx") == :orbit_events
    @test_throws ArgumentError TV._parse_time_aligned_comparison_mode("bad", "ctx")

    @test TV._parse_element_frame("j2000", "ctx") == :j2000
    @test TV._parse_element_frame("EME2000", "ctx") == :j2000
    @test TV._parse_element_frame("body_equator_inertial", "ctx") == :body_equator_inertial
    @test TV._parse_element_frame(" Body_Equator ", "ctx") == :body_equator_inertial
    @test_throws ArgumentError TV._parse_element_frame("venus", "ctx")

    empty_maneuvers = TV._parse_maneuver_config(Dict{String, Any}(), "ctx")
    @test isempty(empty_maneuvers.orbit_numbers)
    valid_maneuvers = TV._parse_maneuver_config(
        Dict("maneuvers" => Dict("orbit_numbers" => [1, 2], "delta_v_mps" => [0.1, 0.2])),
        "ctx"
    )
    @test valid_maneuvers.thrust_n == 4.0
    @test_throws ArgumentError TV._parse_maneuver_config(
        Dict("maneuvers" => Dict("orbit_numbers" => [1], "delta_v_mps" => [0.1, 0.2])),
        "ctx"
    )

    atmo_default = TV._parse_atmosphere_truth_config(Dict{String, Any}(), "ctx")
    @test atmo_default.atmosphere_model == "GRAM"
    atmo_custom = TV._parse_atmosphere_truth_config(
        Dict("atmosphere_truth" => Dict(
            "assumption_id" => "custom_truth",
            "atmosphere_model" => "GRAM",
            "atmosphere_dataset" => "dataset",
            "space_weather_model" => "weather",
            "solar_flux_model" => "solar",
            "gram_seed" => 123,
            "gram_perturbation_scales" => [0.01, 0.02, 0.03, 0.04],
            "gram_offline_surrogate" => "on",
            "gram_static_grid" => true,
            "gram_track_cache" => true,
            "gram_global_lock" => "off",
            "mars_map_year" => 2020,
            "mars_mgcm_dust_levels" => [1.0, 2.0, 3.0],
            "mars_dust_storm" => [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            "mars_f107" => 80.0,
            "mars_wind_scales" => [0.5, 1.5],
            "mars_mola_heights" => true,
            "mars_min_max" => 1
        )),
        "ctx"
    )
    @test atmo_custom.gram_seed == 123
    @test atmo_custom.gram_global_lock == "off"

    cal_default = TV._parse_calibration_config(Dict{String, Any}(), "ctx")
    @test cal_default.enabled == false
    cal_custom = TV._parse_calibration_config(
        Dict("calibration" => Dict(
            "enabled" => true,
            "profiles" => ["quick", "full"],
            "search_on_quick_subset" => false,
            "fit_cd_scale" => false,
            "cd_scale_min" => 1.0,
            "cd_scale_max" => 1.0,
            "cd_scale_steps" => 1,
            "fit_cr" => true,
            "cr_min" => 1.1,
            "cr_max" => 1.2,
            "cr_steps" => 2,
            "fit_bias" => true,
            "bias_abs_max_km" => 10.0,
            "objective" => "mean_rmse_km"
        )),
        "ctx"
    )
    @test cal_custom.enabled == true
    @test :quick in cal_custom.profiles
    @test_throws ArgumentError TV._parse_calibration_config(
        Dict("calibration" => Dict("objective" => "bad")),
        "ctx"
    )

    mktempdir() do tmp
        manifest_path = joinpath(tmp, "manifest.toml")
        peri_path = "data/telemetry/fake_peri.feather"
        apo_path = "data/telemetry/fake_apo.feather"
        tal_path = "data/telemetry/fake_time.feather"

        manifest = Dict(
            "version" => 1,
            "scenarios" => Any[
                Dict(
                    "name" => "oe_case",
                    "kind" => "orbit_events",
                    "planet" => "earth",
                    "events" => ["peri", "apo"],
                    "telemetry_peri" => peri_path,
                    "telemetry_apo" => apo_path,
                    "target_orbits_quick" => 2,
                    "target_orbits_full" => 3,
                    "compare_points_quick" => 2,
                    "compare_points_full" => 3,
                    "min_eval_points" => 1,
                    "ra_m" => 7.1e6,
                    "rp_altitude_m" => 120000.0,
                    "i_deg" => 30.0,
                    "aop_deg" => 20.0,
                    "raan_deg" => 10.0,
                    "ta_deg" => 170.0,
                    "gravity_model" => "inverse_squared",
                    "gravity_harmonics_degree" => 10,
                    "gravity_harmonics_order" => 10,
                    "nbody_bodies" => ["sun"],
                    "srp_enabled" => false,
                    "include_wind" => false,
                    "orbit_altitude_mode" => "vacuum",
                    "EI_km" => 120.0,
                    "maneuvers" => Dict(
                        "orbit_numbers" => [1],
                        "delta_v_mps" => [0.1]
                    ),
                    "initial_time" => Dict(
                        "year" => 2020,
                        "month" => 1,
                        "day" => 1,
                        "hour" => 0,
                        "minute" => 0,
                        "second" => 0.0
                    ),
                    "spacecraft" => Dict(
                        "bus_dims_m" => [1.0, 1.0, 1.0],
                        "panel_dims_m" => [0.1, 0.2, 0.3],
                        "bus_mass_kg" => 100.0,
                        "panel_mass_each_kg" => 5.0,
                        "panel_offset_y_m" => 0.5,
                        "prop_mass_kg" => 10.0,
                        "id" => 1
                    ),
                    "units" => Dict("x" => "orbit", "peri" => "km", "apo" => "km"),
                    "tolerances_quick" => Dict(
                        "peri" => Dict("max_abs_km" => 100.0, "max_nmae" => 1.0),
                        "apo" => Dict("max_abs_km" => 100.0, "max_nmae" => 1.0)
                    ),
                    "tolerances_full" => Dict(
                        "peri" => Dict("max_abs_km" => 80.0, "max_nmae" => 0.9, "max_rmse_km" => 50.0),
                        "apo" => Dict("max_abs_km" => 80.0, "max_nmae" => 0.9, "max_rmse_km" => 50.0)
                    )
                ),
                Dict(
                    "name" => "ta_case",
                    "kind" => "time_aligned_state",
                    "planet" => "earth",
                    "events" => ["peri", "apo"],
                    "comparison_mode" => "orbit_events",
                    "orbit_altitude_mode" => "oblate",
                    "telemetry" => tal_path,
                    "max_points_quick" => 16,
                    "max_points_full" => 24,
                    "min_eval_points" => 2,
                    "gravity_model" => "inverse_squared_j2",
                    "EI_km" => 120.0,
                    "initial_time" => Dict(
                        "year" => 2020,
                        "month" => 1,
                        "day" => 1,
                        "hour" => 0,
                        "minute" => 0,
                        "second" => 0.0
                    ),
                    "telemetry_columns" => Dict(
                        "time" => "t",
                        "altitude" => "alt",
                        "x" => "x",
                        "y" => "y",
                        "z" => "z",
                        "sma" => "sma",
                        "ecc" => "ecc",
                        "inc" => "inc",
                        "aop" => "aop",
                        "raan" => "raan",
                        "ta" => "ta"
                    ),
                    "spacecraft" => Dict(
                        "bus_dims_m" => [1.0, 1.0, 1.0],
                        "panel_dims_m" => [0.1, 0.2, 0.3],
                        "bus_mass_kg" => 100.0,
                        "panel_mass_each_kg" => 5.0,
                        "panel_offset_y_m" => 0.5,
                        "prop_mass_kg" => 10.0,
                        "id" => 2
                    ),
                    "units" => Dict("x" => "time", "peri" => "km", "apo" => "km"),
                    "tolerances_quick" => Dict(
                        "peri" => Dict("max_abs_km" => 100.0, "max_nmae" => 1.0),
                        "apo" => Dict("max_abs_km" => 100.0, "max_nmae" => 1.0)
                    ),
                    "tolerances_full" => Dict(
                        "peri" => Dict("max_abs_km" => 80.0, "max_nmae" => 0.9, "max_rmse_km" => 50.0),
                        "apo" => Dict("max_abs_km" => 80.0, "max_nmae" => 0.9, "max_rmse_km" => 50.0)
                    )
                )
            ]
        )
        open(manifest_path, "w") do io
            TOML.print(io, manifest)
        end

        scenarios = TV._load_scenarios_from_manifest(manifest_path)
        @test length(scenarios) == 2
        @test scenarios[1] isa TV.OrbitEventsScenarioConfig
        @test scenarios[2] isa TV.TimeAlignedScenarioConfig
        sc_oe = scenarios[1]
        sc_ta = scenarios[2]
        @test TV._axis_units(sc_oe) == "orbit"
        @test TV._axis_units(sc_ta) == "time"
        @test TV._value_units(sc_ta, "peri_speed") == "km/s"
        @test TV._source_file(sc_oe, "peri") == normpath(joinpath(TV.REPO_ROOT, peri_path))
        @test TV._source_file(sc_ta, "peri") == normpath(joinpath(TV.REPO_ROOT, tal_path))
        @test TV._orbit_altitude_mode(sc_ta) == "oblate"
        @test TV._maneuver_count(sc_oe) == 1
        @test TV._maneuver_count(sc_ta) == 0
        @test TV._min_eval_points(sc_oe) == 1
        @test TV._min_eval_points(sc_ta) == 2
        @test occursin("maneuvers=", TV._scenario_status_extra(sc_oe))
        @test occursin("comparison_mode=orbit_events", TV._scenario_status_extra(sc_ta))

        parsed_cli = TV.parse_cli([
            "--profile=quick",
            "--manifest=$(manifest_path)",
            "--out-summary=$(joinpath(tmp, "summary.csv"))",
            "--out-errors=$(joinpath(tmp, "errors.csv"))",
            "--enforce=1",
            "--plots=0"
        ])
        @test parsed_cli.profile == :quick
        @test parsed_cli.enforce == true
        @test parsed_cli.generate_plots == false
        @test_throws ArgumentError TV.parse_cli(["--unknown=1"])
    end

    req = TV.VerificationRequest(
        profile=:full,
        out_summary="output/s.csv",
        out_errors="output/e.csv",
        manifest_path="test/telemetry_benchmark_manifest.toml",
        enforce=true,
        generate_plots=false
    )
    cfg_from_req = TV._study_config(req)
    req_roundtrip = TV._request_from_study_config(cfg_from_req)
    @test cfg_from_req.profile == :full
    @test req_roundtrip.enforce == true

    earth = TV._planet_from_name("earth")
    @test TV._period_seconds(earth, earth.Rp_e + 500e3, earth.Rp_e + 400e3) > 0.0
    @test_throws ArgumentError TV._planet_from_name("pluto")
    @test TV._base_gravity_effector(:inverse_squared) isa TV.InverseSquaredGravityModel
    @test TV._base_gravity_effector(:inverse_squared_j2) isa TV.InverseSquaredJ2GravityModel
    @test_throws ArgumentError TV._base_gravity_effector(:bad)
    @test TV._harmonics_order(0, 0) == 0
    @test TV._harmonics_order(10, 0) == 0
    @test TV._harmonics_order(20, 10) == 10
    @test TV._nbody_primary_name("earth") == "Earth"
    @test_throws ArgumentError TV._nbody_primary_name("pluto")

    orb_summary, orb_errors = TV._compare_orbit_curve(
        "scenario",
        "peri",
        [1.0, 2.0, 3.0],
        [100.0, 101.0, 102.0],
        [100.0, 101.5, 102.5]
    )
    @test orb_summary.n_sim == 3
    @test nrow(orb_errors) == 3
    orb_empty_summary, orb_empty_errors = TV._compare_orbit_curve(
        "scenario",
        "peri",
        [1.0, 2.0],
        [100.0, 101.0],
        Float64[]
    )
    @test orb_empty_summary.n_sim == 0
    @test nrow(orb_empty_errors) == 0
    @test_throws ArgumentError TV._compare_orbit_curve(
        "scenario",
        "peri",
        [1.0, 2.0],
        [100.0, 101.0],
        [99.0, 100.0];
        sim_axis=[1.0]
    )

    time_summary, time_errors = TV._compare_time_series(
        "scenario",
        "state_x_time",
        [0.0, 1.0, 2.0],
        [10.0, 11.0, 12.0],
        [0.0, 1.0, 2.0],
        [10.2, 10.9, 12.1]
    )
    @test time_summary.n_sim == 3
    @test nrow(time_errors) == 3
    time_empty_summary, time_empty_errors = TV._compare_time_series(
        "scenario",
        "state_x_time",
        [0.0, 1.0],
        [10.0, 11.0],
        Float64[],
        Float64[]
    )
    @test time_empty_summary.n_sim == 0
    @test nrow(time_empty_errors) == 0

    bias_input = DataFrame(
        scenario=["s", "s"],
        event=["peri", "peri"],
        error_km=[5.0, -3.0]
    )
    bias = TV._estimate_event_biases([bias_input], 1.0)
    @test haskey(bias, "peri")
    @test bias["peri"] <= 1.0

    rows = [
        (nmae=0.2, rmse_km=1.0),
        (nmae=0.4, rmse_km=3.0)
    ]
    @test isapprox(TV._calibration_score(rows, "mean_nmae"), 0.3; atol=1e-12, rtol=0.0)
    @test TV._calibration_score(rows, "mean_rmse_km") == 2.0
    @test TV._calibration_score(rows, "max_nmae") == 0.4
    @test_throws ArgumentError TV._calibration_score(rows, "bad")

    @test TV._display_value_units("km/s") == "m/s"
    @test TV._display_value_units("km") == "km"
    @test TV._display_value_scale("km/s") == 1e3
    @test TV._display_value_scale("km") == 1.0

    summary_df = DataFrame(
        scenario=["s", "s"],
        event=["peri", "peri_speed"],
        value_units=["km", "km/s"],
        mae_km=[1.0, 0.1],
        rmse_km=[1.5, 0.2],
        max_abs_km=[2.0, 0.3],
        p95_abs_km=[1.8, 0.25],
        bias_km=[0.3, 0.02],
        limit_max_abs_km=[3.0, 0.4],
        limit_max_rmse_km=[2.0, 0.3],
        limit_nmae=[1.0, 1.0],
        nmae=[0.1, 0.2],
        n_sim=[5, 5],
        pass=[true, true]
    )
    TV._append_display_metric_columns!(summary_df)
    @test summary_df.value_units_display == ["km", "m/s"]
    @test summary_df.mae_display[2] == 100.0

    errors_df = DataFrame(
        scenario=["s", "s"],
        event=["peri", "peri_speed"],
        telemetry_value_km=[120.0, 7.5],
        sim_interp_value_km=[121.0, 7.6],
        error_km=[1.0, 0.1]
    )
    TV._append_display_error_columns!(errors_df, summary_df)
    @test errors_df.value_units_display == ["km", "m/s"]
    @test errors_df.error_display[2] == 100.0

    cfg_eval = TV.TimeAlignedScenarioConfig(
        name="eval",
        planet_name="earth",
        telemetry_path="unused",
        telemetry_time_col="t",
        telemetry_altitude_col="alt",
        telemetry_x_col="x",
        telemetry_y_col="y",
        telemetry_z_col="z",
        telemetry_sma_col="sma",
        telemetry_ecc_col="ecc",
        telemetry_inc_col="inc",
        telemetry_aop_col="aop",
        telemetry_raan_col="raan",
        telemetry_ta_col="ta",
        max_points_quick=10,
        max_points_full=20,
        min_eval_points=2,
        units_x="orbit",
        units_y=Dict("peri" => "km", "apo" => "km"),
        tolerances_quick=Dict(
            "peri" => (max_abs_km=10.0, max_nmae=1.0, max_rmse_km=5.0),
            "apo" => (max_abs_km=10.0, max_nmae=1.0, max_rmse_km=5.0)
        ),
        tolerances_full=Dict(
            "peri" => (max_abs_km=8.0, max_nmae=0.9, max_rmse_km=4.0),
            "apo" => (max_abs_km=8.0, max_nmae=0.9, max_rmse_km=4.0)
        ),
        initial_time=TV.InitialTime(),
        spacecraft=TV.SpacecraftConfig(
            bus_dims=(1.0, 1.0, 1.0),
            panel_dims=(0.1, 0.2, 0.3),
            bus_mass_kg=10.0,
            panel_mass_each_kg=1.0,
            panel_offset_y_m=0.2,
            prop_mass_kg=2.0,
            id=1
        ),
        gravity_model=:inverse_squared,
        EI_km=120.0,
        comparison_mode=:orbit_events
    )
    ok_row = (
        event="peri_speed",
        scenario="eval",
        n_sim=5,
        max_abs_km=1.0,
        nmae=0.1,
        rmse_km=0.5
    )
    gates = TV._evaluate_thresholds(ok_row, cfg_eval, :quick)
    @test gates.pass == true
    @test gates.limit_max_abs_km == 10.0
    @test_throws ArgumentError TV._evaluate_thresholds(
        (event="missing", scenario="eval", n_sim=1, max_abs_km=1.0, nmae=0.1, rmse_km=0.1),
        cfg_eval,
        :quick
    )

    @test TV._default_plots_outdir("output/a.csv", :quick) == normpath("output/telemetry_plots_quick")
end

@testset "Element-frame transform (kernel-anchored)" begin
    # Expected values below come from the mission SPICE kernels shipped in the
    # GRAMSuite SPICE directory: VEX osculating elements at 2014-05-19T04:00 UTC
    # from ORVV_T19_140501000000_00546.BSP and Odyssey's at 2001-11-06T19:00:32
    # from m01_ab_v2.bsp, expressed in J2000 axes. The manifest stores the same
    # states as body-mean-equator elements; the transform must map one to the other.

    # Identity when the body pole is the J2000 pole (Earth-like scenarios).
    R_id = TV._body_equator_frame_rotation(SVector(0.0, 0.0, 1.0))
    @test isapprox(R_id, SMatrix{3, 3, Float64}(1.0I); atol=1e-14)

    venus_pole = SVector(0.018691, -0.387709, 0.921592)  # IAU_VENUS +z in J2000 (pck00011)
    R_v = TV._body_equator_frame_rotation(venus_pole)
    @test isapprox(R_v' * R_v, SMatrix{3, 3, Float64}(1.0I); atol=1e-10)
    @test isapprox(R_v[:, 3], venus_pole / norm(venus_pole); atol=1e-6)

    j2000_elements = function (r_j::SVector{3, Float64}, v_j::SVector{3, Float64}, μ::Float64)
        n̂ = normalize(cross(r_j, v_j))
        inc = acosd(clamp(n̂[3], -1.0, 1.0))
        node = normalize(SVector(-n̂[2], n̂[1], 0.0))
        raan = mod(atand(node[2], node[1]), 360.0)
        e_vec = ((dot(v_j, v_j) - μ / norm(r_j)) .* r_j .- dot(r_j, v_j) .* v_j) ./ μ
        ê = normalize(e_vec)
        aop = acosd(clamp(dot(node, ê), -1.0, 1.0))
        e_vec[3] < 0.0 && (aop = 360.0 - aop)
        return (i=inc, raan=raan, aop=aop, e=norm(e_vec))
    end

    μ_venus = 3.24858599e14
    ic_vex = SimulationModel.InitialCondition(
        ra=7.2651770e7, rp=6.0518e6 + 186600.0,
        i=89.876, ω=75.505, Ω=104.115, ν=178.0
    )
    r_b, v_b = SimulationEngine.orbitalelemtorv(ic_vex, (μ=μ_venus,))
    oe_vex = j2000_elements(R_v * SVector{3, Float64}(r_b), R_v * SVector{3, Float64}(v_b), μ_venus)
    @test isapprox(oe_vex.i, 84.455; atol=0.05)
    @test isapprox(oe_vex.raan, 105.768; atol=0.05)
    @test isapprox(oe_vex.aop, 97.732; atol=0.06)
    @test isapprox(oe_vex.e, 0.84175; atol=5e-4)

    # Mars pole of date at the Odyssey epoch (IAU 2000 model).
    T_cent = (2452220.2920 - 2451545.0) / 36525.0
    ra_pole = deg2rad(317.68143 - 0.1061 * T_cent)
    dec_pole = deg2rad(52.88650 - 0.0609 * T_cent)
    mars_pole = SVector(cos(dec_pole) * cos(ra_pole), cos(dec_pole) * sin(ra_pole), sin(dec_pole))
    R_m = TV._body_equator_frame_rotation(mars_pole)
    μ_mars = 4.282837285418775e13
    ic_ody = SimulationModel.InitialCondition(
        ra=2.8559615e7, rp=3.396190e6 + 95000.0,
        i=93.522, ω=109.7454, Ω=28.1517, ν=175.0
    )
    r_bm, v_bm = SimulationEngine.orbitalelemtorv(ic_ody, (μ=μ_mars,))
    oe_ody = j2000_elements(R_m * SVector{3, Float64}(r_bm), R_m * SVector{3, Float64}(v_bm), μ_mars)
    @test isapprox(oe_ody.i, 125.45; atol=0.2)
    @test isapprox(oe_ody.raan, 83.0; atol=0.3)
    @test isapprox(oe_ody.aop, 130.1; atol=0.4)

    # Full builder path with the SPICE-derived pole, when kernels are available.
    if isdir(joinpath(TV.SPICE_PATH, "lsk"))
        planet = TV._planet_from_name("venus")
        epoch = TV.InitialTime(year=2014, month=5, day=19, hour=4, minute=0, second=0.0f0)
        @test TV._initial_condition_in_j2000(ic_vex, planet, epoch, :j2000) === ic_vex
        cic = TV._initial_condition_in_j2000(ic_vex, planet, epoch, :body_equator_inertial)
        @test cic isa SimulationModel.CartesianInitialCondition
        r_expected = R_v * SVector{3, Float64}(SimulationEngine.orbitalelemtorv(ic_vex, planet)[1])
        @test norm(cic.pos - r_expected) / norm(r_expected) < 1e-3
        @test_throws ArgumentError TV._initial_condition_in_j2000(ic_vex, planet, epoch, :bogus)
    else
        @info "SPICE kernels not present; skipping SPICE-backed element-frame test"
    end
end

@testset "Kernel Cartesian initial-state override" begin
    # Optional orbit-events key initial_state_j2000_m = [x, y, z, vx, vy, vz] (SI):
    # when present the builder must use the exact Cartesian state and bypass the
    # element-based initial condition; when absent the element path is preserved.
    oe_scenario = Dict{String, Any}(
        "name" => "ic_case",
        "kind" => "orbit_events",
        "planet" => "earth",
        "events" => ["peri", "apo"],
        "telemetry_peri" => "data/telemetry/fake_peri.feather",
        "telemetry_apo" => "data/telemetry/fake_apo.feather",
        "target_orbits_quick" => 2,
        "target_orbits_full" => 3,
        "compare_points_quick" => 2,
        "compare_points_full" => 3,
        "min_eval_points" => 1,
        "ra_m" => 7.1e6,
        "rp_altitude_m" => 120000.0,
        "i_deg" => 30.0,
        "aop_deg" => 20.0,
        "raan_deg" => 10.0,
        "ta_deg" => 170.0,
        "gravity_model" => "inverse_squared",
        "EI_km" => 120.0,
        "initial_time" => Dict(
            "year" => 2020, "month" => 1, "day" => 1,
            "hour" => 0, "minute" => 0, "second" => 0.0
        ),
        "spacecraft" => Dict(
            "bus_dims_m" => [1.0, 1.0, 1.0],
            "panel_dims_m" => [0.1, 0.2, 0.3],
            "bus_mass_kg" => 100.0,
            "panel_mass_each_kg" => 5.0,
            "panel_offset_y_m" => 0.5,
            "prop_mass_kg" => 10.0,
            "id" => 1
        ),
        "units" => Dict("x" => "orbit", "peri" => "km", "apo" => "km"),
        "tolerances_quick" => Dict(
            "peri" => Dict("max_abs_km" => 100.0, "max_nmae" => 1.0),
            "apo" => Dict("max_abs_km" => 100.0, "max_nmae" => 1.0)
        ),
        "tolerances_full" => Dict(
            "peri" => Dict("max_abs_km" => 80.0, "max_nmae" => 0.9),
            "apo" => Dict("max_abs_km" => 80.0, "max_nmae" => 0.9)
        )
    )
    state = (-1.0e6, -4.3e6, -7.3e5, -2371.2, -764.6, -3175.1)

    mktempdir() do tmp
        manifest_path = joinpath(tmp, "manifest.toml")
        write_manifest = scenario -> open(manifest_path, "w") do io
            TOML.print(io, Dict("version" => 1, "scenarios" => Any[scenario]))
        end

        # Absent key: parser leaves the override empty and the builder follows
        # the element path (j2000 elements pass through unchanged).
        write_manifest(oe_scenario)
        cfg_elements = only(TV._load_scenarios_from_manifest(manifest_path))
        @test cfg_elements isa TV.OrbitEventsScenarioConfig
        @test cfg_elements.initial_state_j2000_m === nothing
        planet = TV._planet_from_name("earth")
        ic_el = TV._scenario_initial_condition(cfg_elements, planet)
        @test ic_el isa SimulationModel.InitialCondition
        rp_expected = planet.Rp_e + 120000.0
        @test isapprox(ic_el.a, (7.1e6 + rp_expected) / 2.0; rtol=1e-12)
        @test isapprox(ic_el.e, (7.1e6 - rp_expected) / (7.1e6 + rp_expected); rtol=1e-12)
        @test isapprox(ic_el.ν, deg2rad(170.0); rtol=1e-12)

        # Key present: parser returns the 6-tuple and the builder produces a
        # CartesianInitialCondition carrying exactly the manifest values.
        oe_scenario["initial_state_j2000_m"] = collect(state)
        write_manifest(oe_scenario)
        cfg_state = only(TV._load_scenarios_from_manifest(manifest_path))
        @test cfg_state.initial_state_j2000_m == state
        ic_cart = TV._scenario_initial_condition(cfg_state, planet)
        @test ic_cart isa SimulationModel.CartesianInitialCondition
        @test ic_cart.pos == SVector{3, Float64}(state[1], state[2], state[3])
        @test ic_cart.vel == SVector{3, Float64}(state[4], state[5], state[6])

        # Malformed key: wrong arity is rejected at parse time.
        oe_scenario["initial_state_j2000_m"] = [1.0, 2.0, 3.0]
        write_manifest(oe_scenario)
        @test_throws ArgumentError TV._load_scenarios_from_manifest(manifest_path)
    end

    # The shipped manifest carries the NAV-kernel state for the odyssey campaign.
    shipped = TV._load_scenarios_from_manifest(TV.DEFAULT_MANIFEST_PATH)
    ody = only(filter(s -> s.name == "odyssey", shipped))
    @test ody.initial_state_j2000_m !== nothing
    @test length(ody.initial_state_j2000_m) == 6
end

@testset "Planet-frame transport term uses the body pole" begin
    # planet.ω is the spin vector in planet-fixed axes; the co-rotation term must
    # be applied after rotating into that frame. A point at rest in the rotating
    # frame of a TILTED-pole planet must have zero planet-frame velocity.
    θ = deg2rad(30.0)
    L_PI = SMatrix{3, 3, Float64}(
        1.0, 0.0, 0.0,
        0.0, cos(θ), sin(θ),
        0.0, -sin(θ), cos(θ)
    )
    ω_pf = SVector(0.0, 0.0, 7.0e-5)
    planet_tilted = (L_PI=L_PI, ω=ω_pf)

    r_i = SVector(7.0e6, 1.0e6, 2.0e6)
    Ω_j2000 = L_PI' * ω_pf            # spin vector expressed in J2000
    v_i = cross(Ω_j2000, r_i)         # inertial velocity of a frame-fixed point

    r_p, v_p = SimulationEngine.r_intor_p!(r_i, v_i, planet_tilted)
    @test isapprox(r_p, L_PI * r_i; atol=1e-9)
    @test norm(v_p) < 1e-9 * norm(v_i)

    # Round trip restores the inertial state.
    r_back, v_back = SimulationEngine.r_pintor_i(r_p, v_p, planet_tilted)
    @test isapprox(r_back, r_i; atol=1e-6)
    @test isapprox(v_back, v_i; atol=1e-9)

    # z-aligned pole reproduces the legacy behavior exactly.
    planet_z = (L_PI=SMatrix{3, 3, Float64}(1.0I), ω=ω_pf)
    r_pz, v_pz = SimulationEngine.r_intor_p!(r_i, v_i, planet_z)
    @test isapprox(v_pz, v_i - cross(ω_pf, r_i); atol=1e-9)
end

@testset "Apoapsis decay-rate diagnostic" begin
    # Sim decaying exactly twice as fast as telemetry: ratios must be 2.
    orbits = collect(0.0:1.0:10.0)
    tele = 1000.0 .- 10.0 .* orbits
    sim = 1000.0 .- 20.0 .* orbits
    d = TV._apo_decay_diagnostic(orbits, tele, orbits, sim, Float64[])
    @test isapprox(d.drag_decay_ratio_median, 2.0; atol=1e-12)
    @test isapprox(d.drag_decay_ratio_total, 2.0; atol=1e-12)
    @test d.drag_decay_n == 10

    # A maneuver inside an interval drops exactly that interval.
    d_man = TV._apo_decay_diagnostic(orbits, tele, orbits, sim, [4.5])
    @test d_man.drag_decay_n == 9
    @test isapprox(d_man.drag_decay_ratio_median, 2.0; atol=1e-12)

    # Truncated sim axis restricts the comparison to the overlap.
    d_trunc = TV._apo_decay_diagnostic(orbits, tele, orbits[1:6], sim[1:6], Float64[])
    @test d_trunc.drag_decay_n == 5

    # Flat telemetry (below the rate floor) produces no ratios but a NaN-safe result.
    flat = fill(1000.0, length(orbits))
    d_flat = TV._apo_decay_diagnostic(orbits, flat, orbits, sim, Float64[])
    @test isnan(d_flat.drag_decay_ratio_median)
    @test d_flat.drag_decay_n == 0

    # Too-short series return the schema placeholder.
    d_short = TV._apo_decay_diagnostic([0.0, 1.0], [1.0, 2.0], orbits, sim, Float64[])
    @test d_short.drag_decay_n == 0

    # Uneven telemetry spacing: the total weights each interval by its span
    # (sums altitude deltas), so a 3-orbit gap counts three times a 1-orbit one.
    gap_orbit = [0.0, 1.0, 4.0]
    gap_tele = [1000.0, 990.0, 960.0]        # rates -10, -10
    gap_sim = [1000.0, 970.0, 940.0]         # rates -30, -10
    d_gap = TV._apo_decay_diagnostic(gap_orbit, gap_tele, gap_orbit, gap_sim, Float64[])
    @test isapprox(d_gap.drag_decay_ratio_total, 1.5; atol=1e-12)   # (-30-30)/(-10-30)
    @test isapprox(d_gap.drag_decay_ratio_median, 2.0; atol=1e-12)  # median(3, 1)
end

@testset "Maneuver orbit-number offset" begin
    base = Dict("maneuvers" => Dict(
        "orbit_numbers" => [7, 10, 25, 146],
        "delta_v_mps" => [-0.5, -0.5, 0.1, 1.0],
        "orbit_number_offset" => 19
    ))
    m = TV._parse_maneuver_config(base, "ctx")
    @test m.orbit_numbers == Int64[6, 127]
    @test m.delta_v_mps == [0.1, 1.0]
    # The full campaign-numbered list survives for diagnostics on campaign axes.
    @test m.orbit_numbers_campaign == Int64[7, 10, 25, 146]

    no_offset = Dict("maneuvers" => Dict(
        "orbit_numbers" => [7, 10],
        "delta_v_mps" => [-0.5, -0.5]
    ))
    m0 = TV._parse_maneuver_config(no_offset, "ctx")
    @test m0.orbit_numbers == Int64[7, 10]

    @test_throws ArgumentError TV._parse_maneuver_config(Dict("maneuvers" => Dict(
        "orbit_numbers" => [7], "delta_v_mps" => [1.0], "orbit_number_offset" => -1
    )), "ctx")
end

@testset "Speed-event tolerance tables" begin
    ttbl = Dict(
        "tolerances_full" => Dict(
            "peri" => Dict("max_abs_km" => 2.0, "max_nmae" => 0.04, "max_rmse_km" => 1.0),
            "apo" => Dict("max_abs_km" => 184.0, "max_nmae" => 0.30),
            "apo_speed" => Dict("max_abs_km" => 0.012, "max_nmae" => 0.80, "max_rmse_km" => 0.008)
        )
    )
    tmap = TV._parse_tolerances(ttbl, "tolerances_full", ["peri", "apo"], "ctx")
    # Explicit speed table is parsed alongside its base event...
    @test tmap["apo_speed"].max_abs_km == 0.012
    @test tmap["apo_speed"].max_nmae == 0.80
    @test tmap["apo_speed"].max_rmse_km == 0.008
    # ...while an absent one is simply not present (evaluation falls back to base).
    @test !haskey(tmap, "peri_speed")
    @test tmap["peri"].max_abs_km == 2.0
    @test tmap["apo"].max_rmse_km == Inf

    # A malformed explicit speed table still errors loudly.
    bad = Dict("tolerances_full" => Dict(
        "peri" => Dict("max_abs_km" => 2.0, "max_nmae" => 0.04),
        "apo" => Dict("max_abs_km" => 184.0, "max_nmae" => 0.30),
        "peri_speed" => Dict("max_abs_km" => 0.010)
    ))
    @test_throws ArgumentError TV._parse_tolerances(bad, "tolerances_full", ["peri", "apo"], "ctx")
end

@testset "Libraryless surrogate honors off" begin
    # The library-missing fallback must respect the scenario's opt-out even
    # when the env gate would otherwise allow it.
    env_key = "SPACEAGORA_TELEMETRY_ALLOW_GRAM_OFFLINE_NO_LIB"
    old_env = get(ENV, env_key, nothing)
    ENV[env_key] = "1"
    try
        truth_off = TV.AtmosphereTruthConfig(gram_offline_surrogate="off")
        @test TV._try_libraryless_gram_surrogate("earth", truth_off) === nothing

        # Non-"off" proceeds past the opt-out; in the test environment the
        # surrogate payload is absent, so the file checks return nothing too —
        # this asserts the gate ordering, not payload loading.
        truth_auto = TV.AtmosphereTruthConfig(gram_offline_surrogate="auto")
        @test TV._try_libraryless_gram_surrogate("nonexistent_planet", truth_auto) === nothing

        # And the env gate still applies for non-"off" configs.
        ENV[env_key] = "0"
        @test TV._try_libraryless_gram_surrogate("earth", truth_auto) === nothing
    finally
        old_env === nothing ? delete!(ENV, env_key) : (ENV[env_key] = old_env)
    end
end

@testset "Tabulated flight atmosphere (certification mode)" begin
    # parse: model name validation and key coupling
    @test_throws ArgumentError TV._parse_atmosphere_truth_config(Dict("atmosphere_truth" => Dict(
        "atmosphere_model" => "bogus", "atmosphere_dataset" => "d",
        "space_weather_model" => "s", "solar_flux_model" => "f")), "ctx")
    @test_throws ArgumentError TV._parse_atmosphere_truth_config(Dict("atmosphere_truth" => Dict(
        "atmosphere_model" => "tabulated_flight", "atmosphere_dataset" => "d",
        "space_weather_model" => "s", "solar_flux_model" => "f")), "ctx")   # missing file
    @test_throws ArgumentError TV._parse_atmosphere_truth_config(Dict("atmosphere_truth" => Dict(
        "atmosphere_model" => "GRAM", "atmosphere_dataset" => "d",
        "space_weather_model" => "s", "solar_flux_model" => "f",
        "tabulated_flight_file" => "x.feather")), "ctx")                    # file without mode
    cfgt = TV._parse_atmosphere_truth_config(Dict("atmosphere_truth" => Dict(
        "atmosphere_model" => "tabulated_flight", "atmosphere_dataset" => "d",
        "space_weather_model" => "s", "solar_flux_model" => "f",
        "tabulated_flight_file" => "x.feather", "tabulated_flight_sigma" => -1.0)), "ctx")
    @test cfgt.atmosphere_model == "tabulated_flight"
    @test cfgt.tabulated_flight_sigma == -1.0

    # interpolation math: log-linear interior, exponential tails, sigma shift
    alts = [100.0e3, 110.0e3, 120.0e3]
    logs = [log(1e-7), log(1e-8), log(1e-9)]     # H = 10 km / ln(10)
    sigs = [0.1, 0.2, 0.4]
    m = SimulationModel.TabulatedFlightAtmosphereModel(
        [0.0], [(alts, alts)], [(logs, logs)], [(sigs, sigs)], 0.0, 3.4, 188.92)
    rho_mid, T_mid, w = SimulationModel.getDensity(m, 105.0e3, 0.0, 0.0, 10.0, false)
    @test isapprox(rho_mid, sqrt(1e-7 * 1e-8); rtol=1e-10)                 # geometric mean
    @test w == SimulationModel.SVector{3, Float64}(0.0, 0.0, 0.0)
    rho_above, _, _ = SimulationModel.getDensity(m, 130.0e3, 0.0, 0.0, 10.0, false)
    @test isapprox(rho_above, 1e-10; rtol=1e-6)                            # tail continues H
    rho_below, _, _ = SimulationModel.getDensity(m, 90.0e3, 0.0, 0.0, 10.0, false)
    @test isapprox(rho_below, 1e-6; rtol=1e-6)
    mp = SimulationModel.TabulatedFlightAtmosphereModel(
        [0.0], [(alts, alts)], [(logs, logs)], [(sigs, sigs)], 1.0, 3.4, 188.92)
    rho_p, _, _ = SimulationModel.getDensity(mp, 105.0e3, 0.0, 0.0, 10.0, false)
    @test isapprox(rho_p / rho_mid, exp(0.15); rtol=1e-10)                 # +1 sigma of interp sigma_log
    # leg selection: before the pass periapsis uses inbound profile
    logs_out = [log(2e-7), log(2e-8), log(2e-9)]
    m2 = SimulationModel.TabulatedFlightAtmosphereModel(
        [100.0], [(alts, alts)], [(logs, logs_out)], [(sigs, sigs)], 0.0, 3.4, 188.92)
    rin, _, _ = SimulationModel.getDensity(m2, 105.0e3, 0.0, 0.0, 50.0, false)
    rout, _, _ = SimulationModel.getDensity(m2, 105.0e3, 0.0, 0.0, 150.0, false)
    @test isapprox(rout / rin, 2.0; rtol=1e-10)

    # noise-inverted top bins: the tail scale height clamps to the physical
    # range and density still vanishes far above the profile (no constant
    # blanket along the orbit).
    logs_inv = [log(1e-7), log(1e-8), log(1.5e-8)]
    m3 = SimulationModel.TabulatedFlightAtmosphereModel(
        [0.0], [(alts, alts)], [(logs_inv, logs_inv)], [(sigs, sigs)], 0.0, 3.4, 188.92)
    rho_top, _, _ = SimulationModel.getDensity(m3, 120.0e3, 0.0, 0.0, 10.0, false)
    rho_far, _, _ = SimulationModel.getDensity(m3, 400.0e3, 0.0, 0.0, 10.0, false)
    @test rho_far < rho_top * exp(-(400.0e3 - 120.0e3) / 12000.0) * 1.0001
    @test rho_far < 1e-17
end

@testset "Diagnostic replay scaling (flight apoapsis ratio)" begin
    # Default: mode is delta_v, no flight-apoapsis data.
    m0 = TV._parse_maneuver_config(Dict("maneuvers" => Dict(
        "orbit_numbers" => [25], "delta_v_mps" => [0.1]
    )), "ctx")
    @test m0.replay_scale_mode == "delta_v"
    @test isempty(m0.flight_apoapsis_alt_m)

    # Diagnostic mode: altitudes parsed (km -> m) and filtered in sync with
    # the pre-epoch burn drop.
    m = TV._parse_maneuver_config(Dict("maneuvers" => Dict(
        "orbit_numbers" => [7, 25, 146],
        "delta_v_mps" => [-0.5, 0.1, 1.0],
        "orbit_number_offset" => 19,
        "replay_scale_mode" => "flight_apoapsis_ratio",
        "flight_apoapsis_alt_km" => [26000.0, 23000.0, 9000.0]
    )), "ctx")
    @test m.replay_scale_mode == "flight_apoapsis_ratio"
    @test m.orbit_numbers == Int64[6, 127]
    @test m.flight_apoapsis_alt_m == [23000.0e3, 9000.0e3]

    # Length mismatch, bad values, and data-without-mode all reject.
    @test_throws ArgumentError TV._parse_maneuver_config(Dict("maneuvers" => Dict(
        "orbit_numbers" => [25, 146], "delta_v_mps" => [0.1, 1.0],
        "replay_scale_mode" => "flight_apoapsis_ratio",
        "flight_apoapsis_alt_km" => [23000.0]
    )), "ctx")
    @test_throws ArgumentError TV._parse_maneuver_config(Dict("maneuvers" => Dict(
        "orbit_numbers" => [25], "delta_v_mps" => [0.1],
        "replay_scale_mode" => "flight_apoapsis_ratio",
        "flight_apoapsis_alt_km" => [-1.0]
    )), "ctx")
    @test_throws ArgumentError TV._parse_maneuver_config(Dict("maneuvers" => Dict(
        "orbit_numbers" => [25], "delta_v_mps" => [0.1],
        "flight_apoapsis_alt_km" => [23000.0]
    )), "ctx")
    @test_throws ArgumentError TV._parse_maneuver_config(Dict("maneuvers" => Dict(
        "orbit_numbers" => [25], "delta_v_mps" => [0.1],
        "replay_scale_mode" => "bogus"
    )), "ctx")

    # Scale math: flight/sim apoapsis-radius ratio, clamped, safe fallbacks.
    GM = SimulationModel.GuidanceHooks
    @test GM._flight_ratio_scale(2.0e7, 1.0e7) == 2.0
    @test GM._flight_ratio_scale(1.0e7, 2.0e7) == 0.5
    @test GM._flight_ratio_scale(1.0e7, NaN) == 1.0
    @test GM._flight_ratio_scale(-1.0, 1.0e7) == 1.0
    @test GM._flight_ratio_scale(1.0e9, 1.0) == 10.0
end

@testset "Robust calibration bias" begin
    # Early-orbit datum offset is recovered; late-mission drift is ignored.
    early = fill(-3.0, 10)
    drift = collect(range(-3.0, -4000.0, length=90))
    df = DataFrame(event=fill("apo", 100), error_km=vcat(early, drift))
    biases = TV._estimate_event_biases([df], 500.0)
    @test isapprox(biases["apo"], 3.0; atol=1e-12)

    # A bias at/beyond the cap signals model mismatch and is NOT applied.
    df_sat = DataFrame(event=fill("apo", 20), error_km=fill(-800.0, 20))
    biases_sat = TV._estimate_event_biases([df_sat], 500.0)
    @test biases_sat["apo"] == 0.0
end

@testset "Comparison masking and raw values" begin
    tele_axis = collect(1.0:10.0)
    tele = collect(100.0:-1.0:91.0)
    sim = collect(100.5:-1.0:96.5)   # five events, sim axis 1..5

    masked, mdf = TV._compare_orbit_curve(
        "s", "apo", tele_axis, tele, sim;
        sim_axis=collect(1.0:5.0), bias=2.0, mask_to_sim_span=true
    )
    @test masked.n_telemetry == 10
    @test masked.n_sim == 5
    @test isapprox(masked.coverage, 0.5; atol=1e-12)
    @test nrow(mdf) == 5

    # The mask bounds BOTH sides: telemetry preceding the simulated span is not
    # scored against the clamped first sim value.
    pre, pdf = TV._compare_orbit_curve(
        "s", "apo", tele_axis, tele, sim;
        sim_axis=collect(3.0:7.0), mask_to_sim_span=true
    )
    @test pdf.telemetry_axis[1] >= 2.5
    @test pdf.telemetry_axis[end] <= 7.5
    @test nrow(pdf) == 5
    @test isapprox(pre.coverage, 0.5; atol=1e-12)
    @test all(mdf.sim_interp_value_km .== mdf.sim_raw_interp_value_km .+ 2.0)
    @test isapprox(masked.max_abs_km, 2.5; atol=1e-9)  # constant 0.5 raw offset + 2.0 bias

    # Legacy behavior (no masking) still clamp-scores the full series.
    legacy, ldf = TV._compare_orbit_curve(
        "s", "apo", tele_axis, tele, sim;
        sim_axis=collect(1.0:5.0)
    )
    @test nrow(ldf) == 10
    @test legacy.max_abs_km > masked.max_abs_km  # clamped tail dominates
end

@testset "Airspeed sign contract" begin
    # Drag/heat models must use v_spacecraft MINUS v_atmosphere. The flipped
    # form (`vel_pp + wind_pp`) was copy-pasted across six files; pin all of them.
    files = [
        joinpath(_COV_REPO_ROOT, "src", "dynamics", "coupled", "aerodynamic_wrench_models.jl"),
        joinpath(_COV_REPO_ROOT, "src", "simulation", "callbacks", "thermal_callbacks.jl"),
        joinpath(_COV_REPO_ROOT, "src", "gnc", "guidance", "aerobraking", "t_edg", "trajectory_predictor.jl"),
        joinpath(_COV_REPO_ROOT, "src", "gnc", "guidance", "aerobraking", "t_edg", "eom_predictor.jl"),
        joinpath(_COV_REPO_ROOT, "src", "gnc", "control", "aerobraking", "constraint_tracking.jl"),
        joinpath(_COV_REPO_ROOT, "src", "gnc", "control", "aerobraking", "control_commands.jl"),
    ]
    for file in files
        live = [line for line in eachline(file) if !startswith(strip(line), "#")]
        @test !any(occursin(r"vel_pp_rw\s*=\s*(planet_frame\.)?vel_pp\s*\+\s*wind_pp", line) for line in live)
    end
    # Dynamic pressure in the coupled aero model must use the wind-relative speed,
    # in both the spacecraft-level and the per-link (rho_body/q_body) forms.
    aero_src = read(files[1], String)
    @test !occursin(r"q(_body)?\s*=\s*0\.5\s*\*\s*rho(_body)?\s*\*\s*vel_pp_mag\^2", aero_src)
    # The per-link Mach/molecular-speed-ratio in _aero_pure_wrench must also be
    # wind-relative. (The legacy calcForceTorque paths' spacecraft-level `mach =
    # vel_pp_mag / ...` is long-standing pinned behavior, deliberately not
    # covered here.)
    @test !occursin(r"mach_body\s*=\s*vel_pp_mag\s*/", aero_src)
end

@testset "Spacecraft builder ram-face convention" begin
    # The Hart free-molecular coefficients are normalized by the flow-normal
    # face; with the translational fixed attitude the flow runs along body +x,
    # so the consistent bus reference area is dims[2]*dims[3] (:frontal).
    # :legacy pins the historical dims[1]*dims[3] value so previously
    # calibrated scenarios are unchanged, and it must stay the default.
    ic = SimulationModel.InitialCondition(
        7.0e6, 1.0e-3, 35.0, 0.0, 0.0, 0.0,
        SVector{4, Float64}(0.0, 0.0, 0.0, 1.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    )
    dims = (0.2, 0.5, 0.6)
    kw = (panel_dims=(0.4, 0.001, 0.5), bus_mass=29.0, panel_mass_each=0.0,
          panel_offset_y=0.5, ic=ic)
    legacy = TelemetryVerification.make_three_body_spacecraft(; bus_dims=dims, kw...)
    @test legacy.links[1].ref_area ≈ dims[1] * dims[3]
    frontal = TelemetryVerification.make_three_body_spacecraft(; bus_dims=dims, bus_ram_face=:frontal, kw...)
    @test frontal.links[1].ref_area ≈ dims[2] * dims[3]
    # Panels already use the flow-normal dims[2]*dims[3] face in both modes.
    @test legacy.links[2].ref_area ≈ 0.001 * 0.5
    @test frontal.links[2].ref_area ≈ 0.001 * 0.5
    @test_throws ArgumentError TelemetryVerification.make_three_body_spacecraft(;
        bus_dims=dims, bus_ram_face=:sideways, kw...)
    # Manifest key default must remain :legacy.
    @test TelemetryVerification.SpacecraftConfig(
        bus_dims=dims, panel_dims=(0.4, 0.001, 0.5), bus_mass_kg=29.0,
        panel_mass_each_kg=0.0, panel_offset_y_m=0.5, prop_mass_kg=0.0, id=1
    ).bus_ram_face === :legacy
end

@testset "Replication guards: density source and tolerance env vars" begin
    # 1) Manifest whitelist accepts the NRLMSISE-00 source and still rejects
    #    unknown values.
    atm_tbl(model) = Dict{String, Any}(
        "atmosphere_truth" => Dict{String, Any}(
            "atmosphere_model" => model,
            "atmosphere_dataset" => "probe",
            "space_weather_model" => "probe",
            "solar_flux_model" => "probe"
        )
    )
    @test TV._parse_atmosphere_truth_config(atm_tbl("nrlmsise00"), "probe").atmosphere_model == "nrlmsise00"
    @test_throws ArgumentError TV._parse_atmosphere_truth_config(atm_tbl("msise"), "probe")

    # 2) The scenario builder wires that manifest value to the native
    #    NRLMSISE-00 model, and the source stays Earth-only.
    mk_cfg(planet) = TV.TimeAlignedScenarioConfig(
        name="probe", planet_name=planet, telemetry_path="",
        telemetry_time_col="t", telemetry_altitude_col="alt",
        telemetry_x_col="x", telemetry_y_col="y", telemetry_z_col="z",
        max_points_quick=10, max_points_full=10, min_eval_points=1,
        units_x="s", units_y=Dict{String, String}(),
        tolerances_quick=Dict{String, TV.EventTolerance}(),
        tolerances_full=Dict{String, TV.EventTolerance}(),
        initial_time=SimulationModel.InitialTime(year=2025, month=6, day=6),
        spacecraft=TV.SpacecraftConfig(
            bus_dims=(0.2, 0.5, 0.6), panel_dims=(0.4, 0.001, 0.5),
            bus_mass_kg=29.0, panel_mass_each_kg=0.0, panel_offset_y_m=0.5,
            prop_mass_kg=0.0, id=1
        ),
        gravity_model=:inverse_squared_j2,
        atmosphere_truth=TV.AtmosphereTruthConfig(atmosphere_model="nrlmsise00"),
        EI_km=600.0
    )
    @test TV._scenario_density_model(mk_cfg("earth")) isa SimulationModel.NRLMSISE00AtmosphereModel
    @test_throws ArgumentError TV._scenario_density_model(mk_cfg("mars"))

    # 3) The SPACEAGORA_TELEMETRY_{RELTOL,ABSTOL}_{ORBIT,ATM} env vars tighten
    #    the study tolerances but can never loosen them (they were historically
    #    set by callers and silently ignored).
    ic = SimulationModel.InitialCondition(
        7.0e6, 1.0e-3, 35.0, 0.0, 0.0, 0.0,
        SVector{4, Float64}(0.0, 0.0, 0.0, 1.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    )
    sc = TV.make_three_body_spacecraft(
        bus_dims=(0.2, 0.5, 0.6), panel_dims=(0.4, 0.001, 0.5), bus_mass=29.0,
        panel_mass_each=0.0, panel_offset_y=0.5, ic=ic
    )
    args = TV.make_example_config(
        planet=TV._planet_from_name("earth"), spacecraft=sc, mission_time=60.0,
        initial_time=SimulationModel.InitialTime(year=2025, month=6, day=6),
        verbose=false
    )
    base = withenv(
        "SPACEAGORA_TELEMETRY_RELTOL_ORBIT" => nothing,
        "SPACEAGORA_TELEMETRY_ABSTOL_ORBIT" => nothing,
        "SPACEAGORA_TELEMETRY_RELTOL_ATM" => nothing,
        "SPACEAGORA_TELEMETRY_ABSTOL_ATM" => nothing
    ) do
        TV._with_study_settings(args)
    end
    @test base.integration_tolerances.reltol_orbit == 1.0e-7
    @test base.integration_tolerances.abstol_atmosphere == 1.0e-9
    tightened = withenv(
        "SPACEAGORA_TELEMETRY_RELTOL_ORBIT" => "1e-9",
        "SPACEAGORA_TELEMETRY_ABSTOL_ATM" => "1e-12"
    ) do
        TV._with_study_settings(args)
    end
    @test tightened.integration_tolerances.reltol_orbit == 1.0e-9
    @test tightened.integration_tolerances.abstol_atmosphere == 1.0e-12
    loosened = withenv("SPACEAGORA_TELEMETRY_RELTOL_ORBIT" => "1e-3") do
        TV._with_study_settings(args)
    end
    @test loosened.integration_tolerances.reltol_orbit == 1.0e-7
end

@testset "Scenario builder branch probes (coverage backstop)" begin
    # These branches are otherwise exercised only by the scenario suites, whose
    # worker-process coverage can fall outside the CI gate's active window; pin
    # them from this in-process probe file so the gate is deterministic.

    # Planet and N-body primary name mappings.
    @test TV._planet_from_name("mars") isa SimulationModel.Mars
    @test TV._planet_from_name("venus") isa SimulationModel.Venus
    @test TV._planet_from_name("moon") isa SimulationModel.Moon
    @test_throws ArgumentError TV._planet_from_name("pluto")
    @test TV._nbody_primary_name("earth") == "Earth"
    @test TV._nbody_primary_name("mars") == "Mars"
    @test TV._nbody_primary_name("venus") == "Venus"
    @test TV._nbody_primary_name("moon") == "Moon"
    @test TV._nbody_primary_name("titan") == "Titan"
    @test_throws ArgumentError TV._nbody_primary_name("pluto")

    # Tabulated-flight density builder on a synthetic two-pass table with an
    # archive gap (P=1,3), both legs, and a nonpositive-density row that must
    # be skipped.
    mktempdir() do dir
        path = joinpath(dir, "flight_table.arrow")
        tbl = DataFrame(
            P=[1, 1, 1, 3, 3, 3],
            leg=["in", "out", "in", "in", "out", "out"],
            alt_km=[120.0, 130.0, 110.0, 121.0, 131.0, 141.0],
            rho_kgm3=[1.0e-9, 5.0e-10, 0.0, 1.1e-9, 6.0e-10, 3.0e-10],
            sigma_kgm3=[1.0e-10, 5.0e-11, 0.0, 1.1e-10, 6.0e-11, 3.0e-11],
            t_peri_utc=vcat(fill("2025-06-06T01:00:00", 3), fill("2025-06-06T03:00:00", 3))
        )
        TV.Arrow.write(path, tbl)
        truth = TV.AtmosphereTruthConfig(
            atmosphere_model="tabulated_flight",
            tabulated_flight_file=path,
            tabulated_flight_sigma=0.0
        )
        it = SimulationModel.InitialTime(year=2025, month=6, day=6)
        model = TV._make_tabulated_flight_density_model(it, truth)
        @test model isa SimulationModel.TabulatedFlightAtmosphereModel

        @test_throws ArgumentError TV._make_tabulated_flight_density_model(
            it,
            TV.AtmosphereTruthConfig(
                atmosphere_model="tabulated_flight",
                tabulated_flight_file=joinpath(dir, "missing.arrow")
            )
        )
        bad = joinpath(dir, "bad_columns.arrow")
        TV.Arrow.write(bad, DataFrame(P=[1], leg=["in"]))
        @test_throws ArgumentError TV._make_tabulated_flight_density_model(
            it,
            TV.AtmosphereTruthConfig(
                atmosphere_model="tabulated_flight",
                tabulated_flight_file=bad
            )
        )
    end
end

@testset "Time-tabulated density source (tabulated_time)" begin
    # Direct model semantics: log-linear interpolation in elapsed time, held
    # ends, multiplicative scale, constant temperature, zero wind.
    m = SimulationModel.TimeTabulatedAtmosphereModel(
        [0.0, 100.0], [1.0e-12, 4.0e-12]; scale=2.0, temperature_k=800.0
    )
    rho0, T0, w0 = SimulationModel.getDensity(m, 4.0e5, 0.0, 0.0, 0.0, false)
    @test rho0 ≈ 2.0e-12
    @test T0 == 800.0
    @test w0 == SVector{3, Float64}(0.0, 0.0, 0.0)
    rho_mid, _, _ = SimulationModel.getDensity(m, 4.0e5, 0.0, 0.0, 50.0, false)
    @test rho_mid ≈ 2.0 * exp(0.5 * (log(1.0e-12) + log(4.0e-12)))  # log-linear midpoint
    rho_end, _, _ = SimulationModel.getDensity(m, 4.0e5, 0.0, 0.0, 1.0e6, false)
    @test rho_end ≈ 8.0e-12  # held past the last node
    @test_throws ArgumentError SimulationModel.TimeTabulatedAtmosphereModel([0.0], [1.0e-12])
    @test_throws ArgumentError SimulationModel.TimeTabulatedAtmosphereModel([100.0, 0.0], [1.0e-12, 2.0e-12])
    @test_throws ArgumentError SimulationModel.TimeTabulatedAtmosphereModel([0.0, 1.0], [1.0e-12, -1.0e-12])
    @test_throws ArgumentError SimulationModel.TimeTabulatedAtmosphereModel([0.0, 1.0], [1.0e-12, 2.0e-12]; scale=0.0)

    # Manifest contract: whitelisted, file required with the mode, file key
    # rejected without the mode.
    atm_tbl(extra) = Dict{String, Any}(
        "atmosphere_truth" => merge(Dict{String, Any}(
            "atmosphere_model" => "tabulated_time",
            "atmosphere_dataset" => "probe",
            "space_weather_model" => "probe",
            "solar_flux_model" => "probe"
        ), extra)
    )
    @test_throws ArgumentError TV._parse_atmosphere_truth_config(atm_tbl(Dict{String, Any}()), "probe")
    parsed = TV._parse_atmosphere_truth_config(
        atm_tbl(Dict{String, Any}("tabulated_time_file" => "some.csv", "tabulated_time_scale" => 1.5)), "probe"
    )
    @test parsed.atmosphere_model == "tabulated_time"
    @test parsed.tabulated_time_scale == 1.5
    @test_throws ArgumentError TV._parse_atmosphere_truth_config(
        Dict{String, Any}("atmosphere_truth" => Dict{String, Any}(
            "atmosphere_model" => "GRAM", "atmosphere_dataset" => "probe",
            "space_weather_model" => "probe", "solar_flux_model" => "probe",
            "tabulated_time_file" => "some.csv"
        )), "probe"
    )

    # Builder loads a CSV table through the scenario dispatch.
    mktempdir() do dir
        path = joinpath(dir, "rho_t.csv")
        write(path, "time_s,rho_kgm3\n0.0,1.0e-12\n3600.0,2.0e-12\n")
        truth = TV.AtmosphereTruthConfig(
            atmosphere_model="tabulated_time",
            tabulated_time_file=path,
            tabulated_time_scale=1.0
        )
        model = TV._make_time_tabulated_density_model(truth)
        @test model isa SimulationModel.TimeTabulatedAtmosphereModel
        rho, _, _ = SimulationModel.getDensity(model, 4.0e5, 0.0, 0.0, 0.0, false)
        @test rho ≈ 1.0e-12
        bad = joinpath(dir, "bad.csv")
        write(bad, "t,rho\n0.0,1.0e-12\n")
        @test_throws ArgumentError TV._make_time_tabulated_density_model(
            TV.AtmosphereTruthConfig(atmosphere_model="tabulated_time", tabulated_time_file=bad)
        )
        @test_throws ArgumentError TV._make_time_tabulated_density_model(
            TV.AtmosphereTruthConfig(atmosphere_model="tabulated_time", tabulated_time_file=joinpath(dir, "missing.csv"))
        )
    end
end

@testset "Decay diagnostics: synthetic injection" begin
    # Synthetic 48 h LEO arc with gaps, a drifting-amplitude 2/rev harmonic
    # (the J2-precession leakage mechanism), and a known injected decay.
    mu = 3.98600436233e14
    a0 = 6.82e6
    P = 2.0 * pi * sqrt(a0^3 / mu)
    t_full = collect(0.0:10.0:172800.0)
    keep = [(mod(ti, 9000.0) > 1800.0) for ti in t_full]   # periodic gaps
    t = t_full[keep]
    harm(ti) = (1800.0 + 0.004 * ti) * sin(2.0 * pi * 2.0 * ti / P + 0.7) +
               400.0 * sin(2.0 * pi * 1.0 * ti / P - 0.3) +
               3.0 * sin(31.7 * ti)                         # deterministic "noise"
    slope_true_mpd = -55.0
    sma_meas = [a0 + (slope_true_mpd / 86400.0) * ti + harm(ti) for ti in t]
    sma_ref = [a0 + harm(ti) for ti in t]

    # The hazard: the drag-free reference alone shows a nonzero apparent slope
    # (drifting-amplitude leakage through the fixed-amplitude fit).
    ref_raw = TV.secular_sma_slope(t, sma_ref; period_s=P)
    @test abs(ref_raw) > 0.5

    # Zero-reference subtraction recovers the injected decay.
    z = TV.zero_referenced_decay(t, sma_meas, t, sma_ref; period_s=P)
    @test abs(z.decay_m_per_day - slope_true_mpd) < 0.02
    @test z.reference_m_per_day ≈ ref_raw

    # Estimator invariance of the corrected number (fixed vs drifting fit).
    zd = TV.zero_referenced_decay(t, sma_meas, t, sma_ref; period_s=P, drifting_amplitudes=true)
    @test abs(zd.decay_m_per_day - z.decay_m_per_day) < 0.1

    # vis-viva helper round-trip.
    v0 = sqrt(mu / a0)
    @test TV.visviva_sma(a0, v0, mu) ≈ a0
    @test TV.visviva_sma([a0, a0], [v0, v0], mu) ≈ [a0, a0]

    # Windowed density extraction: constant-density decay in, density out —
    # and the emitted table drives the tabulated_time source directly.
    rho_true = 2.0e-12
    cd_area = 0.167
    mass = 28.94
    adot = -rho_true * cd_area / mass * sqrt(mu * a0)      # m/s
    sma_decay = [a0 + adot * ti + harm(ti) for ti in t]
    tbl = TV.flight_density_table(
        t, sma_decay, t, sma_ref;
        mu=mu, cd_area_m2=cd_area, mass_kg=mass, period_s=P
    )
    @test names(tbl) == ["time_s", "rho_kgm3"]
    @test nrow(tbl) > 20
    @test all(abs.(tbl.rho_kgm3 .- rho_true) ./ rho_true .< 0.05)
    m = SimulationModel.TimeTabulatedAtmosphereModel(tbl.time_s, tbl.rho_kgm3)
    rho_mid, _, _ = SimulationModel.getDensity(m, 4.0e5, 0.0, 0.0, 86400.0, false)
    @test abs(rho_mid - rho_true) / rho_true < 0.05

    # Error paths.
    @test_throws ArgumentError TV.secular_sma_slope(t[1:5], sma_ref[1:5]; period_s=P)
    @test_throws ArgumentError TV.secular_sma_slope(t, sma_ref[1:10]; period_s=P)
    @test_throws ArgumentError TV.flight_density_table(
        t, sma_decay, t, sma_ref; mu=mu, cd_area_m2=0.0, mass_kg=mass, period_s=P
    )
end

# Under --code-coverage this testset is skipped: it executes the full
# verification pipeline in-process, which would pull telemetry_loading.jl,
# runner.jl, and error_tables.jl into the coverage gate's active window at
# the partial coverage of this file's calls (the per-file floors are meant
# for the scenario suites' coverage of those files, not this integration
# probe). The testset runs in the plain tests and tests-matrix jobs.
if Base.JLOptions().code_coverage != 0
    @info "Skipping IC-offset/mask/fit integration probes under coverage instrumentation (see comment)"
else
@testset "IC offsets, illumination mask, differential-correction fit" begin
    # Solar direction helper vs the CYGNSS campaign's screening values
    # (2025-06-06: RA 75.6 deg, dec 22.8 deg from ephemerides).
    it = SimulationModel.InitialTime(year=2025, month=6, day=6)
    s_hat = TV._sun_unit_vector_j2000(it, 0.0)
    ra = rad2deg(atan(s_hat[2], s_hat[1]))
    dec = rad2deg(asin(s_hat[3]))
    @test abs(ra - 75.6) < 1.5
    @test abs(dec - 22.8) < 0.5

    # IC offsets apply on the Cartesian path and throw on the Keplerian path.
    mk_cfg(; kw...) = TV.TimeAlignedScenarioConfig(;
        name="probe", planet_name="earth", telemetry_path="",
        telemetry_time_col="t", telemetry_altitude_col="alt",
        telemetry_x_col="x", telemetry_y_col="y", telemetry_z_col="z",
        max_points_quick=10, max_points_full=10, min_eval_points=1,
        units_x="s", units_y=Dict{String, String}(),
        tolerances_quick=Dict{String, TV.EventTolerance}(),
        tolerances_full=Dict{String, TV.EventTolerance}(),
        initial_time=it,
        spacecraft=TV.SpacecraftConfig(
            bus_dims=(0.2, 0.5, 0.6), panel_dims=(0.4, 0.001, 0.5),
            bus_mass_kg=29.0, panel_mass_each_kg=0.0, panel_offset_y_m=0.5,
            prop_mass_kg=0.0, id=1
        ),
        gravity_model=:inverse_squared_j2, EI_km=600.0, kw...
    )
    tel_cart = (x_ic_km=7000.0, y_ic_km=0.0, z_ic_km=0.0,
                vx_ic_kmps=0.0, vy_ic_kmps=7.5, vz_ic_kmps=0.0,
                sma_km=NaN, ecc=NaN, inc_deg=NaN, aop_deg=NaN, raan_deg=NaN, ta_deg=NaN)
    cfg_off = mk_cfg(ic_offset_m=(150.0, 0.0, -20.0), ic_offset_mps=(0.0, -0.003, 0.0))
    ic = TV._initial_condition_from_time_aligned_telemetry(cfg_off, tel_cart)
    @test ic isa SimulationModel.CartesianInitialCondition
    @test ic.pos[1] ≈ 7000.0e3 + 150.0
    @test ic.pos[3] ≈ -20.0
    @test ic.vel[2] ≈ 7500.0 - 0.003
    tel_kep = (x_ic_km=nothing, y_ic_km=nothing, z_ic_km=nothing,
               vx_ic_kmps=nothing, vy_ic_kmps=nothing, vz_ic_kmps=nothing,
               sma_km=7000.0, ecc=0.001, inc_deg=30.0, aop_deg=0.0, raan_deg=0.0, ta_deg=0.0)
    @test_throws ArgumentError TV._initial_condition_from_time_aligned_telemetry(cfg_off, tel_kep)
    @test_throws ArgumentError TV._parse_truth_mask("twilight", "probe")
    @test TV._parse_truth_mask("nightside", "probe") === :nightside

    # Synthetic circular two-body scenario shared by the mask and fit checks.
    planet = TV._planet_from_name("earth")
    mu = planet.μ
    a = 7.0e6
    n_mm = sqrt(mu / a^3)
    inc = deg2rad(30.0)
    t_samp = collect(0.0:60.0:2.0 * 2.0 * pi / n_mm)
    ci, si = cos(inc), sin(inc)
    xs = similar(t_samp); ys = similar(t_samp); zs = similar(t_samp)
    vxs = similar(t_samp); vys = similar(t_samp); vzs = similar(t_samp)
    for (i, ti) in enumerate(t_samp)
        th = n_mm * ti
        xs[i] = a * cos(th); ys[i] = a * sin(th) * ci; zs[i] = a * sin(th) * si
        v = a * n_mm
        vxs[i] = -v * sin(th); vys[i] = v * cos(th) * ci; vzs[i] = v * cos(th) * si
    end
    mktempdir() do dir
        tele = DataFrame(
            time_s=t_samp,
            altitude_km=(sqrt.(xs .^ 2 .+ ys .^ 2 .+ zs .^ 2) .- planet.Rp_e) .* 1e-3,
            x_km=xs .* 1e-3, y_km=ys .* 1e-3, z_km=zs .* 1e-3,
            x_ic_km=fill(xs[1] * 1e-3, length(t_samp)),
            y_ic_km=fill(ys[1] * 1e-3, length(t_samp)),
            z_ic_km=fill(zs[1] * 1e-3, length(t_samp)),
            vx_ic_kmps=fill(vxs[1] * 1e-3, length(t_samp)),
            vy_ic_kmps=fill(vys[1] * 1e-3, length(t_samp)),
            vz_ic_kmps=fill(vzs[1] * 1e-3, length(t_samp))
        )
        tele_path = joinpath(dir, "synthetic_circular.arrow")
        TV.Arrow.write(tele_path, tele)

        # Illumination mask splits the samples into complementary sets whose
        # kept members satisfy the screen's sign convention.
        cols = Dict(
            :telemetry_path => tele_path, :telemetry_time_col => "time_s",
            :telemetry_altitude_col => "altitude_km",
            :telemetry_x_col => "x_km", :telemetry_y_col => "y_km", :telemetry_z_col => "z_km",
            :telemetry_x_ic_col => "x_ic_km", :telemetry_y_ic_col => "y_ic_km",
            :telemetry_z_ic_col => "z_ic_km", :telemetry_vx_ic_col => "vx_ic_kmps",
            :telemetry_vy_ic_col => "vy_ic_kmps", :telemetry_vz_ic_col => "vz_ic_kmps"
        )
        cfg_none = mk_cfg(; cols...)
        cfg_night = mk_cfg(; truth_mask=:nightside, cols...)
        cfg_day = mk_cfg(; truth_mask=:dayside, cols...)
        full = TV._load_time_aligned_telemetry(cfg_none, 0)
        night = TV._load_time_aligned_telemetry(cfg_night, 0)
        day = TV._load_time_aligned_telemetry(cfg_day, 0)
        @test length(night.time_s) + length(day.time_s) == length(full.time_s)
        @test 0.3 < length(night.time_s) / length(full.time_s) < 0.7
        for i in eachindex(night.time_s)
            sh = TV._sun_unit_vector_j2000(it, night.time_s[i])
            r = SVector{3, Float64}(night.x_km[i], night.y_km[i], night.z_km[i])
            @test dot(r, sh) <= 0.0
        end

        # Differential-correction fit recovers a deliberately wrong IC.
        tol6 = Dict{String, Any}(k => Dict("max_abs_km" => 1.0e6, "max_nmae" => 1.0e6, "max_rmse_km" => 1.0e6)
            for k in ("altitude_time", "state_x_time", "state_y_time", "state_z_time"))
        scenario = Dict{String, Any}(
            "name" => "icfit_probe", "kind" => "time_aligned_state",
            "events" => Any["altitude_time", "state_x_time", "state_y_time", "state_z_time"],
            "telemetry" => tele_path,
            "telemetry_columns" => Dict{String, Any}(
                "time" => "time_s", "altitude" => "altitude_km",
                "x" => "x_km", "y" => "y_km", "z" => "z_km",
                "x_ic" => "x_ic_km", "y_ic" => "y_ic_km", "z_ic" => "z_ic_km",
                "vx_ic" => "vx_ic_kmps", "vy_ic" => "vy_ic_kmps", "vz_ic" => "vz_ic_kmps"),
            "max_points_quick" => 1000, "max_points_full" => 1000, "min_eval_points" => 2,
            "units" => Dict{String, Any}("x" => "s", "altitude_time" => "km",
                "state_x_time" => "km", "state_y_time" => "km", "state_z_time" => "km"),
            "tolerances_quick" => tol6, "tolerances_full" => tol6,
            "planet" => "earth", "gravity_model" => "inverse_squared",
            "nbody_bodies" => Any[], "drag_enabled" => false, "include_wind" => false,
            "EI_km" => 600.0,
            "spacecraft" => Dict{String, Any}(
                "bus_dims_m" => Any[0.2, 0.5, 0.6], "panel_dims_m" => Any[0.4, 0.001, 0.5],
                "bus_mass_kg" => 29.0, "panel_mass_each_kg" => 0.0,
                "panel_offset_y_m" => 0.5, "prop_mass_kg" => 0.0, "id" => 1),
            "initial_time" => Dict{String, Any}("year" => 2025, "month" => 6, "day" => 6,
                "hour" => 0, "minute" => 0, "second" => 0.0),
            "ic_offset_m" => Any[150.0, -80.0, 40.0],
            "ic_offset_mps" => Any[0.0, -0.003, 0.001]
        )
        manifest_path = joinpath(dir, "icfit_manifest.toml")
        open(manifest_path, "w") do io
            TOML.print(io, Dict{String, Any}("scenarios" => Any[scenario]))
        end
        # P1 regressions (Codex review): a mask that screens out the EARLY
        # samples must not shorten the propagation (the mission runs through
        # the last retained timestamp with the IC still anchored at t=0), and
        # a Keplerian-IC scenario must keep its original-epoch elements.
        masked = deepcopy(scenario)
        masked["name"] = "icfit_probe_masked"
        masked["truth_mask"] = "nightside"   # first samples are dayside -> dropped
        masked["ic_offset_m"] = Any[0.0, 0.0, 0.0]
        masked["ic_offset_mps"] = Any[0.0, 0.0, 0.0]
        masked_manifest = joinpath(dir, "masked_manifest.toml")
        open(masked_manifest, "w") do io
            TOML.print(io, Dict{String, Any}("scenarios" => Any[masked]))
        end
        req = TV.VerificationRequest(
            profile=:quick,
            out_summary=joinpath(dir, "masked_summary.csv"),
            out_errors=joinpath(dir, "masked_errors.csv"),
            manifest_path=masked_manifest, enforce=false, generate_plots=false)
        rmask = TV.run_verification(req)
        rows = rmask.summary[in.(rmask.summary.event, Ref(["state_x_time", "state_y_time", "state_z_time"])), :]
        # With the exact IC and pure two-body dynamics the masked comparison must
        # be tight EVERYWHERE — a shortened propagation clamps late samples to the
        # final state and produces km-scale garbage.
        @test maximum(Float64.(rows.max_abs_km)) < 0.01
        ed = TV.CSV.read(joinpath(dir, "masked_errors.csv"), DataFrame)
        ed = ed[ed.event .== "state_x_time", :]
        t_last = maximum(Float64.(ed.telemetry_axis))
        @test t_last > 0.75 * t_samp[end]     # late retained samples were compared
        @test abs(Float64(ed.error_km[argmax(Float64.(ed.telemetry_axis))])) < 0.01

        # Keplerian elements survive a first-sample-dropping mask at original
        # values (loader-level check with per-row varying element columns).
        kep = DataFrame(
            time_s=t_samp,
            altitude_km=(sqrt.(xs .^ 2 .+ ys .^ 2 .+ zs .^ 2) .- planet.Rp_e) .* 1e-3,
            x_km=xs .* 1e-3, y_km=ys .* 1e-3, z_km=zs .* 1e-3,
            sma_km=fill(7000.0, length(t_samp)) .+ collect(0:length(t_samp)-1),
            ecc=fill(0.001, length(t_samp)),
            inc_deg=fill(30.0, length(t_samp)),
            aop_deg=fill(0.0, length(t_samp)),
            raan_deg=fill(0.0, length(t_samp)),
            ta_deg=fill(0.0, length(t_samp))
        )
        kep_path = joinpath(dir, "synthetic_kep.arrow")
        TV.Arrow.write(kep_path, kep)
        cfg_kep = mk_cfg(; truth_mask=:nightside,
            telemetry_path=kep_path, telemetry_time_col="time_s",
            telemetry_altitude_col="altitude_km",
            telemetry_x_col="x_km", telemetry_y_col="y_km", telemetry_z_col="z_km",
            telemetry_sma_col="sma_km", telemetry_ecc_col="ecc",
            telemetry_inc_col="inc_deg", telemetry_aop_col="aop_deg",
            telemetry_raan_col="raan_deg", telemetry_ta_col="ta_deg")
        loaded_kep = TV._load_time_aligned_telemetry(cfg_kep, 0)
        @test loaded_kep.time_s[1] > 0.0       # the mask did drop the first samples
        @test loaded_kep.sma_km == 7000.0      # row 1 of the ORIGINAL data

        fitwork = joinpath(dir, "fitwork")
        mkpath(fitwork)
        fit = TV.fit_initial_state(manifest_path, "icfit_probe"; workdir=fitwork)
        @test fit.rmse_before_km > 0.05                       # the injected error is visible
        @test all(abs.(fit.offsets_total_m) .< 10.0)          # fit cancels the position offset
        # Position/velocity trade along the orbit (a ~0.15 m residual position
        # offset exchanges against ~0.5 mm/s over this 2-orbit two-body arc),
        # so the velocity gate is the degeneracy scale, not measurement zero.
        @test all(abs.(fit.offsets_total_mps) .< 1.0e-3)
        @test fit.rmse_after_km < 0.05 * fit.rmse_before_km   # validated propagation confirms
    end
end
end

@testset "Per-link atmosphere scoping (PR 49 review follow-up)" begin
    AE = SimulationModel.DynamicEffectors.AerodynamicEffectors
    # Defaults unchanged: field false, global false -> disabled.
    m = AerodynamicCoefficientfM()
    @test m.per_link_atmosphere === false
    @test AE._per_link_enabled(m) === false
    # Model-scoped enable works without touching the process-wide switch.
    m_on = AerodynamicCoefficientfM(per_link_atmosphere=true)
    @test AE._per_link_enabled(m_on) === true
    @test AE.PER_LINK_ATMOSPHERE_ENABLED[] === false
    # The process-wide switch still works (compatibility) and resets.
    AE.set_per_link_atmosphere!(true)
    @test AE._per_link_enabled(m) === true
    AE.set_per_link_atmosphere!(false)
    @test AE._per_link_enabled(m) === false
    # Field-less models fall back to the global without erroring.
    @test AE._per_link_enabled(AerodynamicCoefficientConstant()) === false
end

@testset "Magnetic field: dipole sign fix + IGRF source option" begin
    PE = SimulationModel.DynamicEffectors.PerturbationEffectors
    ES = parentmodule(PE.StateSample)
    # Constructor contract: default stays the tilted dipole; :igrf requires a
    # finite decimal year; unknown sources are rejected.
    @test PE.MagneticTorqueRodModel().field_model === :dipole
    m_igrf = PE.MagneticTorqueRodModel(field_model=:igrf, igrf_year=2025.4)
    @test m_igrf.igrf_year == 2025.4
    @test_throws ArgumentError PE.MagneticTorqueRodModel(field_model=:igrf)
    @test_throws ArgumentError PE.MagneticTorqueRodModel(field_model=:wmm)
    # The IGRF library hard-rejects epochs outside [1900, 2035); the
    # constructor must catch that at configuration time (Codex P2, PR 64).
    @test_throws ArgumentError PE.MagneticTorqueRodModel(field_model=:igrf, igrf_year=2050.0)
    @test_throws ArgumentError PE.MagneticTorqueRodModel(field_model=:igrf, igrf_year=1899.0)

    # Tilted-dipole SIGN pins. The pre-fix implementation used the north-pole
    # axis as the dipole moment and returned the antiparallel field (~170 deg
    # from IGRF at LEO). Physical pins: at the magnetic equator B = +B0 along
    # the pole axis (north); above the north magnetic pole B = -2 B0 (down).
    l_pi = SMatrix{3, 3, Float64, 9}(I)
    n_hat = SVector{3, Float64}(PE.M_HAT_ECEF)
    r_eq = 6.3712e6 * normalize(cross(n_hat, SVector(0.0, 0.0, 1.0)))
    @test isapprox(PE.get_magnetic_field_dipole(r_eq, MMatrix{3, 3, Float64}(l_pi)),
                   3.12e-5 * n_hat; atol=1e-9)
    @test isapprox(PE.get_magnetic_field_dipole(6.3712e6 * n_hat, MMatrix{3, 3, Float64}(l_pi)),
                   -2 * 3.12e-5 * n_hat; atol=1e-9)

    # IGRF-vs-dipole cross-validation at LEO sample points: LEO-band magnitude
    # and the two models roughly aligned (< 35 deg) — this catches any future
    # sign/frame regression in either path.
    earth = PE.Earth()
    m_dip = PE.MagneticTorqueRodModel()
    for (latd, lond) in ((0.0, 10.0), (30.0, 45.0), (-30.0, 135.0), (55.0, -100.0))
        r = 6898e3 * SVector(cosd(latd) * cosd(lond), cosd(latd) * sind(lond), sind(latd))
        alt, lat, lon = PE.rtolatlong(r, earth)
        B_i = PE._magnetic_field_inertial(m_igrf, l_pi, r, lat, lon, alt)
        B_d = PE._magnetic_field_inertial(m_dip, l_pi, r, lat, lon, alt)
        @test 1.5e-5 < norm(B_i) < 7e-5
        @test B_d == PE.get_magnetic_field_dipole(r, MMatrix{3, 3, Float64}(l_pi))
        @test acosd(clamp(dot(B_i, B_d) / (norm(B_i) * norm(B_d)), -1, 1)) < 35
    end

    # wrench plumbing: one magnet, identity attitude -> torque = m x B_ii for
    # BOTH field sources, zero force, and tau ⊥ m.
    ic = SimulationModel.InitialCondition(
        7.0e6, 1.0e-3, 35.0, 0.0, 0.0, 0.0,
        SVector{4, Float64}(0.0, 0.0, 0.0, 1.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    )
    sc = TV.make_three_body_spacecraft(
        bus_dims=(0.2, 0.5, 0.6), panel_dims=(0.4, 0.001, 0.5), bus_mass=29.0,
        panel_mass_each=0.0, panel_offset_y=0.5, ic=ic
    )
    MT = eltype(sc.links[1].magnets)
    m_vec = SVector(0.0, 0.0, 1.5)
    push!(sc.links[1].magnets, MT(m=MVector{3, Float64}(m_vec)))
    r = 6898e3 * SVector(cosd(20.0) * cosd(30.0), cosd(20.0) * sind(30.0), sind(20.0))
    alt, lat, lon = PE.rtolatlong(r, earth)
    pf = ES.PlanetFrameSample(l_pi, r, SVector(0.0, 0.0, 0.0), alt, lat, lon)
    x = PE.StateSample(r, SVector(7.6e3, 0.0, 0.0), 29.0;
                       q_ib=SVector(0.0, 0.0, 0.0, 1.0), spacecraft=sc)
    env = PE.EnvironmentSample(earth; planet_frame=pf)
    for model in (m_dip, m_igrf)
        B_ii = PE._magnetic_field_inertial(model, l_pi, r, lat, lon, alt)
        f, tq = PE.wrench(model, x, env, 0.0)
        @test f == SVector(0.0, 0.0, 0.0)
        @test isapprox(tq, cross(m_vec, B_ii); atol=1e-12)
        @test abs(dot(tq, m_vec)) < 1e-15
    end
end

@testset "LVLH cascade attitude controller" begin
    PE = SimulationModel.DynamicEffectors.PerturbationEffectors
    ES = parentmodule(PE.StateSample)
    # Generic round-number gains: probe values are NOT mission calibrations.
    mk(; kw...) = PE.LVLHCascadeAttitudeControlModel(;
        k_out=[0.01, 0.01, 0.01], w_max=1.5e-3, k_rate=[0.05, 0.05, 0.05], tau_max=1.0e-3, kw...)
    m = mk()
    @test m.q_cmd_lb == SVector(0.0, 0.0, 0.0, 1.0)
    @test_throws ArgumentError mk(w_max=0.0)
    @test_throws ArgumentError mk(tau_max=-1.0)
    @test_throws ArgumentError mk(k_out=[-0.01, 0.0, 0.0])
    @test_throws ArgumentError mk(q_cmd_lb=[0.0, 0.0, 0.0, 0.0])

    # Quaternion q with PE.rot(q) == R: Shepperd on the TRANSPOSE (PE.rot is
    # the transpose of the standard scalar-last DCM; pinned below).
    function dcm_to_quat(R_target)
        R = R_target'
        tr = R[1, 1] + R[2, 2] + R[3, 3]
        if tr > 0
            s = sqrt(tr + 1.0) * 2
            q = SVector((R[3, 2] - R[2, 3]) / s, (R[1, 3] - R[3, 1]) / s, (R[2, 1] - R[1, 2]) / s, s / 4)
        elseif R[1, 1] > R[2, 2] && R[1, 1] > R[3, 3]
            s = sqrt(1.0 + R[1, 1] - R[2, 2] - R[3, 3]) * 2
            q = SVector(s / 4, (R[1, 2] + R[2, 1]) / s, (R[1, 3] + R[3, 1]) / s, (R[3, 2] - R[2, 3]) / s)
        elseif R[2, 2] > R[3, 3]
            s = sqrt(1.0 + R[2, 2] - R[1, 1] - R[3, 3]) * 2
            q = SVector((R[1, 2] + R[2, 1]) / s, s / 4, (R[2, 3] + R[3, 2]) / s, (R[1, 3] - R[3, 1]) / s)
        else
            s = sqrt(1.0 + R[3, 3] - R[1, 1] - R[2, 2]) * 2
            q = SVector((R[1, 3] + R[3, 1]) / s, (R[2, 3] + R[3, 2]) / s, s / 4, (R[2, 1] - R[1, 2]) / s)
        end
        return q / norm(q)
    end
    lvlh_dcm(r, v) = begin
        z_l = -r / norm(r); y_l = normalize(cross(z_l, v)); x_l = cross(y_l, z_l)
        SMatrix{3, 3, Float64, 9}(x_l[1], y_l[1], z_l[1], x_l[2], y_l[2], z_l[2], x_l[3], y_l[3], z_l[3])
    end
    r0 = SVector(6.9e6, 0.0, 0.0)
    v0 = SVector(0.0, sqrt(3.98600436233e14 / 6.9e6), 0.0)
    R_li = lvlh_dcm(r0, v0)
    q_aligned = dcm_to_quat(R_li)
    @test maximum(abs.(PE.rot(q_aligned) - R_li)) < 1e-12    # dcm_to_quat consistent with rot
    w_lvlh_body = PE.rot(q_aligned) * (cross(r0, v0) / dot(r0, r0))

    # equilibrium: LVLH-aligned, LVLH-matched rates -> zero torque, zero force
    x_eq = PE.StateSample(r0, v0, 29.0; q_ib=q_aligned, ω_body=w_lvlh_body)
    env = PE.EnvironmentSample(nothing)
    f, tq = PE.wrench(m, x_eq, env, 0.0)
    @test f == SVector(0.0, 0.0, 0.0)
    @test norm(tq) < 1e-15

    # restoring: small roll offset (body vs LVLH), LVLH-matched rates ->
    # torque opposes the error, and only that axis responds
    q_off = SVector(sind(1.0), 0.0, 0.0, cosd(1.0))          # 2 deg about body x
    q_ib = dcm_to_quat(PE.rot(q_off) * R_li)
    @test maximum(abs.(PE.rot(q_ib) - PE.rot(q_off) * R_li)) < 1e-12  # composition pin
    # rate consistent with the PERTURBED attitude, so w_rel = 0 and the
    # response is purely the outer loop acting on the roll error
    x_roll = PE.StateSample(r0, v0, 29.0; q_ib=q_ib,
                            ω_body=PE.rot(q_ib) * (cross(r0, v0) / dot(r0, r0)))
    _, tq_r = PE.wrench(m, x_roll, env, 0.0)
    th_err = 2.0 * sind(1.0)  # small-angle error magnitude, rad-ish
    @test tq_r[1] * th_err < 0                                # restoring on the error axis
    @test abs(tq_r[2]) < 0.05 * abs(tq_r[1]) && abs(tq_r[3]) < 0.05 * abs(tq_r[1])

    # rate damping: aligned attitude, extra body rate -> tau = -k_rate .* dw
    dw = SVector(2.0e-4, -1.0e-4, 5.0e-5)
    x_rate = PE.StateSample(r0, v0, 29.0; q_ib=q_aligned, ω_body=w_lvlh_body + dw)
    _, tq_w = PE.wrench(m, x_rate, env, 0.0)
    @test isapprox(tq_w, -m.k_rate .* dw; rtol=1e-10)

    # saturation: large error saturates the rate command; torque respects tau_max
    q_big = dcm_to_quat(PE.rot(SVector(sind(45.0), 0.0, 0.0, cosd(45.0))) * R_li)
    x_big = PE.StateSample(r0, v0, 29.0; q_ib=q_big, ω_body=w_lvlh_body)
    _, tq_b = PE.wrench(m, x_big, env, 0.0)
    @test all(abs.(tq_b) .<= m.tau_max + 1e-18)
    @test abs(tq_b[1]) <= m.k_rate[1] * m.w_max * (1 + 1e-9) + m.k_rate[1] * norm(w_lvlh_body)

    # antipodal robustness (Codex review): at 180 deg error the vee-map
    # vanishes, but the rotation-log extraction must keep full authority —
    # nonzero, axis-dominant, restoring, rate command saturated.
    for ang in (180.0, 170.0)
        q_flip = dcm_to_quat(PE.rot(SVector(sind(ang / 2), 0.0, 0.0, cosd(ang / 2))) * R_li)
        x_flip = PE.StateSample(r0, v0, 29.0; q_ib=q_flip,
                                ω_body=PE.rot(q_flip) * (cross(r0, v0) / dot(r0, r0)))
        _, tq_f = PE.wrench(m, x_flip, env, 0.0)
        @test norm(tq_f) > 0.5 * m.k_rate[1] * m.w_max
        @test abs(tq_f[1]) > 5 * max(abs(tq_f[2]), abs(tq_f[3]))
    end

    # no attitude state -> zero wrench
    x_noq = PE.StateSample(r0, v0, 29.0)
    @test PE.wrench(m, x_noq, env, 0.0) == (SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0))

    # rigid-body invariant pin for the frame-consistency fix: torque-free
    # asymmetric tumble under the ENGINE's own kinematics + Euler equations
    # must conserve inertial angular momentum (the pre-fix inertial-rate
    # quaternion composition drifted it by ~68% on this exact test).
    DR = SimulationModel.DynamicsRotational
    I_ten = SMatrix{3, 3, Float64, 9}(1.5, 0, 0, 0, 1.0, 0, 0, 0, 2.0)
    q_t = SVector(0.0, 0.0, 0.0, 1.0); w_t = SVector(0.02, 0.03, 0.01)
    L0 = PE.rot(q_t)' * (I_ten * w_t)
    for _ in 1:40000
        w_t = w_t + DR.angular_acceleration(w_t, I_ten, SVector(0.0, 0.0, 0.0)) * 0.05
        q_t = q_t + DR.quaternion_derivative(w_t, q_t) * 0.05
        q_t = q_t / norm(q_t)
    end
    @test norm(PE.rot(q_t)' * (I_ten * w_t) - L0) / norm(L0) < 0.01

    # end-to-end convergence pin: rigid body + the engine's quaternion
    # kinematics + this controller, closed loop from a 5 deg offset.
    qdot(w, q) = DR.quaternion_derivative(SVector{3, Float64}(w), q)
    I_diag = SVector(1.5, 1.0, 2.0)
    n_orb = norm(cross(r0, v0)) / dot(r0, r0)
    q = dcm_to_quat(PE.rot(SVector(sind(2.5), 0.0, 0.0, cosd(2.5))) * R_li)  # 5 deg roll offset
    w_state = cross(r0, v0) / dot(r0, r0)
    dt = 0.5
    # circular-orbit propagation in the equatorial plane
    orbit_r(tK) = SVector(6.9e6 * cos(n_orb * tK), 6.9e6 * sin(n_orb * tK), 0.0)
    orbit_v(tK) = SVector(-v0[2] * sin(n_orb * tK), v0[2] * cos(n_orb * tK), 0.0)
    err_angle(q, r, v) = begin
        R_e = PE.rot(q) * lvlh_dcm(r, v)'
        acosd(clamp((R_e[1, 1] + R_e[2, 2] + R_e[3, 3] - 1) / 2, -1, 1))
    end
    e_start = err_angle(q, r0, v0)
    for k in 0:Int(2000 / dt)
        tK = k * dt
        r, v = orbit_r(tK), orbit_v(tK)
        tau = PE._lvlh_cascade_torque(m, r, v, q, w_state)
        w_state = w_state + (tau ./ I_diag) * dt
        q = q + qdot(w_state, q) * dt
        q = q / norm(q)
    end
    e_end = err_angle(q, orbit_r(2000.0), orbit_v(2000.0))
    @test e_start > 4.5
    @test e_end < 0.35                                        # >10x collapse, no divergence

    # tau_ff: default is zero and bit-identical; a constant disturbance holds
    # a steady PD offset of |tau_d|/(k_rate*k_out), and the matching
    # feedforward nulls it (closed-loop mini-sim below reuses the engine
    # kinematics; disturbance 2e-6 N m about x -> predicted offset ~0.23 deg).
    m_ff = mk(tau_ff=[2.0e-6, 0.0, 0.0])
    @test mk().tau_ff == SVector(0.0, 0.0, 0.0)
    @test PE.wrench(mk(), x_eq, env, 0.0)[2] == PE.wrench(m, x_eq, env, 0.0)[2]
    @test_throws ArgumentError mk(tau_ff=[Inf, 0.0, 0.0])
    function settle(ctrl, tau_d)
        qs = dcm_to_quat(R_li); ws = cross(r0, v0) / dot(r0, r0)
        for k in 0:Int(3000 / 0.5)
            tK = k * 0.5
            r, v = orbit_r(tK), orbit_v(tK)
            tau = PE._lvlh_cascade_torque(ctrl, r, v, qs, ws) + tau_d
            ws = ws + (tau ./ I_diag) * 0.5
            qs = qs + qdot(ws, qs) * 0.5
            qs = qs / norm(qs)
        end
        return err_angle(qs, orbit_r(3000.0), orbit_v(3000.0))
    end
    tau_d = SVector(-2.0e-6, 0.0, 0.0)
    off_pd = settle(mk(), tau_d)
    off_ff = settle(m_ff, tau_d)
    @test off_pd > 0.15                       # PD alone holds a steady offset
    @test off_ff < 0.2 * off_pd               # matching feedforward nulls it

end

@testset "Magnetic momentum manager (discrete ZOH control effector)" begin
    CH = SimulationModel.ControlHooks
    # mock integrator state: duck-typed sc accessor (pos/vel/q/ω), one sat
    r0 = SVector(6.9e6, 0.0, 0.0)
    v0 = SVector(0.0, 7.6e3, 0.0)
    q_id = SVector(0.0, 0.0, 0.0, 1.0)          # identity: body == inertial
    w0 = SVector(0.0, 0.0, 0.0)
    mock_u(q) = (sc = [(pos = r0, vel = v0, q = q, ω = w0)],)
    B0 = SVector(0.0, 0.0, 3.0e-5)              # +z inertial field, 30 uT
    bfun = (t, r) -> B0
    tau_const = SVector(1.0e-6, 0.0, 0.0)       # constant commanded torque
    taufun = (t, r, v, q, w) -> tau_const

    mk(; kw...) = CH.MagneticMomentumManagerModel(;
        mu_gain=0.01, h_wheels_0=SVector(0.02, 0.0, 0.0),
        commanded_torque=taufun, b_field_ii=bfun, kw...)

    # 1) first tick initializes the accumulator (no advance) and holds a command
    m = mk()
    CH.calcControlEffect!(m, mock_u(q_id), nothing, 0.0, 1)
    @test m.h_wheels == SVector(0.02, 0.0, 0.0)
    @test m.initialized
    f, tq = CH.calcControlForceTorque(m, nothing, nothing, 1, 0.0)
    @test f == SVector(0.0, 0.0, 0.0)
    # unloading: tau_rod = -mu * H_perp; H = (0.02,0,0) fully perp to +z field
    @test isapprox(tq, SVector(-0.01 * 0.02, 0.0, 0.0); rtol=1e-12)
    @test dot(tq, m.h_wheels) < 0.0

    # 2) ZOH: held command unchanged until the next tick
    tq_held = m.held_torque_body
    @test CH.calcControlForceTorque(m, nothing, nothing, 1, 0.37)[2] == tq_held

    # 3) accumulator advance: dH = -tau_cmd * dt at the tick
    CH.calcControlEffect!(m, mock_u(q_id), nothing, 1.0, 1)
    @test isapprox(m.h_wheels, SVector(0.02, 0.0, 0.0) - tau_const * 1.0; rtol=1e-12)

    # 4) field-parallel momentum is untouched by the law (tau ⟂ H_parallel)
    mpar = mk(h_wheels_0=SVector(0.0, 0.0, 0.05))     # H along +z == B direction
    CH.calcControlEffect!(mpar, mock_u(q_id), nothing, 0.0, 1)
    @test norm(mpar.held_torque_body) < 1e-15

    # 5) dipole cap: |m| clamped, torque scales down accordingly
    mcap = mk(m_max_am2=1.0)
    CH.calcControlEffect!(mcap, mock_u(q_id), nothing, 0.0, 1)
    @test isapprox(norm(mcap.held_dipole_am2), 1.0; rtol=1e-12)
    muncap = mk()
    CH.calcControlEffect!(muncap, mock_u(q_id), nothing, 0.0, 1)
    @test norm(muncap.held_dipole_am2) > 1.0   # confirms the cap actually bound

    # 6) body-frame consistency: rotated attitude gives the same PHYSICAL torque
    #    (rotate body 90 deg about +z; field/momentum vectors fixed in inertial
    #    space => body components rotate, torque back-rotated must match)
    ang = pi / 2
    q_z90 = SVector(0.0, 0.0, sin(ang / 2), cos(ang / 2))
    PE = SimulationModel.DynamicEffectors.PerturbationEffectors
    R = PE.rot(q_z90)
    mrot = mk(h_wheels_0=SVector{3, Float64}(R * SVector(0.02, 0.0, 0.0)))
    CH.calcControlEffect!(mrot, mock_u(q_z90), nothing, 0.0, 1)
    m0 = mk()
    CH.calcControlEffect!(m0, mock_u(q_id), nothing, 0.0, 1)
    @test isapprox(R' * mrot.held_torque_body, m0.held_torque_body; atol=1e-18)

    # 7) other-sat and degenerate-field guards; no propellant
    moff = mk(sat_idx=2)
    CH.calcControlEffect!(moff, mock_u(q_id), nothing, 0.0, 1)
    @test !moff.initialized
    @test CH.calcControlForceTorque(moff, nothing, nothing, 1, 0.0)[2] == SVector(0.0, 0.0, 0.0)
    mzero = mk(b_field_ii=(t, r) -> SVector(0.0, 0.0, 0.0))
    CH.calcControlEffect!(mzero, mock_u(q_id), nothing, 0.0, 1)
    @test mzero.held_torque_body == SVector(0.0, 0.0, 0.0)

    # 8) discrete boundedness of the law itself: constant secular disturbance,
    #    ticked at 1 Hz — with the manager the accumulated momentum plateaus at
    #    tau/mu; with mu=0 it grows linearly (drift = the unmanaged twin).
    #    (commanded torque = disturbance counter, i.e. wheels absorb +tau_d)
    tau_d = SVector(2.0e-6, 0.0, 0.0)
    # emulate the closed loop: at each tick the commanded torque the wheels
    # absorb is the disturbance counter MINUS the rod torque counter
    hold = SVector(0.0, 0.0, 0.0)
    mgr2 = mk(h_wheels_0=SVector(0.0, 0.0, 0.0),
              commanded_torque=(t, r, v, q, w) -> -tau_d - hold)
    drift = SVector(0.0, 0.0, 0.0)
    for k in 0:2000
        CH.calcControlEffect!(mgr2, mock_u(q_id), nothing, Float64(k), 1)
        hold = mgr2.held_torque_body
        drift = drift + tau_d * 1.0          # unmanaged accumulator, same ticks
    end
    @test norm(mgr2.h_wheels) < 0.5 * norm(drift)
    @test norm(mgr2.h_wheels) < 1.5 * norm(tau_d) / 0.01   # plateau ~ tau/mu
end

println("coverage_parallel_telemetry_probes_ok")
