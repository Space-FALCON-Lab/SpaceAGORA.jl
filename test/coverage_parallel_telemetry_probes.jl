using Test
using TOML
using DataFrames
using LinearAlgebra
using StaticArrays

const _COV_REPO_ROOT = isdefined(Main, :REPO_ROOT) ? Main.REPO_ROOT : normpath(joinpath(@__DIR__, ".."))

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
    @test cfg_full.thermal_mode == "on"
    @test cfg_full.inner_scheduler == "dynamic"
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
    @test env_map_auto["SPACEAGORA_THERMAL_CALLBACK_PARALLEL"] == "on"
    @test env_map_auto["SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER"] == "dynamic"
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
    # Dynamic pressure in the coupled aero model must use the wind-relative speed.
    aero_src = read(files[1], String)
    @test !occursin("q = 0.5 * rho * vel_pp_mag^2", aero_src)
end

println("coverage_parallel_telemetry_probes_ok")
