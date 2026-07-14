@testset "Monolith Split Runtime Contracts" begin
    TV = TelemetryVerification
    PP = SimulationModel.ParallelPolicy

    @testset "Telemetry manifest and request parity" begin
        manifest = joinpath(REPO_ROOT, "test", "telemetry_benchmark_manifest.toml")
        scenarios = TV._load_scenarios_from_manifest(manifest)
        @test !isempty(scenarios)
        @test any(sc -> sc isa TV.OrbitEventsScenarioConfig, scenarios)
        @test any(sc -> sc isa TV.TimeAlignedScenarioConfig, scenarios)

        parsed = TV.parse_cli([
            "quick",
            "--manifest=$(manifest)",
            "--out-summary=$(joinpath(REPO_ROOT, "output", "telemetry_tmp_summary.csv"))",
            "--out-errors=$(joinpath(REPO_ROOT, "output", "telemetry_tmp_errors.csv"))",
            "--enforce=false",
            "--plots=false",
        ])
        request = TV._request_from_study_config(parsed)
        reparsed = TV._study_config(request)
        @test reparsed.profile == parsed.profile
        @test reparsed.manifest_path == parsed.manifest_path
        @test reparsed.generate_plots == parsed.generate_plots
    end

    @testset "Telemetry comparison and calibration helpers" begin
        telemetry_axis = [1.0, 2.0, 3.0]
        telemetry_values = [100.0, 101.0, 102.0]
        sim_values = [100.5, 101.5, 102.5]
        summary, errors = TV._compare_orbit_curve("demo", "peri", telemetry_axis, telemetry_values, sim_values)
        @test summary.n_telemetry == 3
        @test summary.n_sim == 3
        @test nrow(errors) == 3
        @test isapprox(summary.bias_km, 0.5; atol=1e-12)

        bias_map = TV._estimate_event_biases([errors], 0.25)
        @test bias_map["peri"] == -0.25

        rows = [(nmae=0.3, rmse_km=2.0), (nmae=0.5, rmse_km=4.0)]
        @test TV._calibration_score(rows, "mean_nmae") == 0.4
        @test TV._calibration_score(rows, "mean_rmse_km") == 3.0
        @test TV._calibration_score(rows, "max_nmae") == 0.5
    end

    @testset "Parallel policy decision and persistent hint parity" begin
        mktempdir() do tmp
            state_path = joinpath(tmp, "hints.toml")
            withenv(
                "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "1",
                "SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => "1",
                "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => state_path,
                "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
                "SPACEAGORA_PARALLEL_POLICY_MEASURED_REWARD" => "1",
                "SPACEAGORA_PARALLEL_POLICY_HINT_MIN_SAMPLES" => "1",
            ) do
                PP.reset_persistent_hint_state!()
                PP.reset_policy_telemetry!()
                decision = PP.thread_policy_decision(16; mode=:auto, threshold=2, source=:probe_split)
                @test decision.num_items == 16
                @test decision.allotment >= 1
                PP.record_policy_observation!(:probe_split; mode=:auto, num_items=16, use_threads=decision.use_threads, elapsed_ns=10_000)
                lock(PP._persistent_hint_lock) do
                    PP._save_persistent_hint_state_locked!()
                end
                @test isfile(state_path)

                PP.reset_persistent_hint_state!()
                PP._ensure_persistent_hint_state_loaded!()
                rows = PP.hint_layer_stats_snapshot(profile=get(ENV, "SPACEAGORA_PARALLEL_PROFILE", "default"), machine=get(ENV, "SPACEAGORA_PERF_MACHINE_LABEL", "default"))
                @test rows !== nothing
                @test !isempty(rows)
                @test any(row -> row.layer == "probe_split", rows)
            end
        end
    end
end
