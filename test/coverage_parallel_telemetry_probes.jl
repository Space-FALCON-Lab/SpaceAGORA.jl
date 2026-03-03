using Test
using TOML
using DataFrames

const _COV_REPO_ROOT = isdefined(Main, :REPO_ROOT) ? Main.REPO_ROOT : normpath(joinpath(@__DIR__, ".."))

include(joinpath(_COV_REPO_ROOT, "src", "parallel", "ParallelProfiles.jl"))
include(joinpath(_COV_REPO_ROOT, "src", "analysis", "verification", "TelemetryVerification.jl"))

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
        PP.R4_full_auto => "R4_full_auto"
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
    @test PP.parse_parallel_profile("full_smart") == PP.R4_full_auto
    @test PP.parse_parallel_profile(:R0) == PP.R0
    @test PP.parse_parallel_profile(PP.R3) == PP.R3
    @test_throws ArgumentError PP.parse_parallel_profile("not_a_profile")

    cfg_r0 = PP.profile_config("R0")
    cfg_r4 = PP.profile_config("R4")
    cfg_full = PP.profile_config("R4_full_auto")
    @test cfg_r0.outer_backend == :none
    @test cfg_r4.inner_adaptive == true
    @test cfg_full.outer_route_adaptive == true
    @test cfg_full.label == "r4_full_auto"

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

    sig = PP.outer_route_signature(f_heavy)
    @test occursin("cat=deterministic", sig)
    @test occursin("harm=21p", sig)

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

    cand_heavy = PP.outer_route_candidates(f_heavy; tuning=tune, machine_class=:large, threads_available=true, parallel_enabled=true)
    cand_light = PP.outer_route_candidates(f_light; tuning=tune, machine_class=:small, threads_available=false, parallel_enabled=true)
    cand_mc = PP.outer_route_candidates(f_mc; tuning=tune, machine_class=:large, threads_available=true, parallel_enabled=true)
    @test :process in cand_heavy
    @test cand_light == [:none]
    @test :process in cand_mc

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
        peri_path = "telemetry/fake_peri.feather"
        apo_path = "telemetry/fake_apo.feather"
        tal_path = "telemetry/fake_time.feather"

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
    @test TV._harmonics_order(10, 0) == 10
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

println("coverage_parallel_telemetry_probes_ok")
