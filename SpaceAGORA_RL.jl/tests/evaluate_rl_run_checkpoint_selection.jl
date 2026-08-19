using Test

include(joinpath(@__DIR__, "..", "scripts", "evaluate_rl_run.jl"))

@testset "evaluation checkpoint selection" begin
    mktempdir() do run_dir
        checkpoint_10 = joinpath(run_dir, "checkpoint_10.jls")
        checkpoint_20 = joinpath(run_dir, "checkpoint_20.jls")
        checkpoint_final = joinpath(run_dir, "checkpoint_final.jls")
        touch(checkpoint_10)
        touch(checkpoint_20)
        touch(checkpoint_final)

        validation_dir = joinpath(run_dir, "checkpoint_validation")
        mkpath(validation_dir)
        open(joinpath(validation_dir, "best_validation_checkpoint.txt"), "w") do io
            println(io, "checkpoint=checkpoint_10.jls")
            println(io, "checkpoint_path=/stale/location/checkpoint_10.jls")
        end

        @test _checkpoint_path(run_dir, nothing) == checkpoint_10
        @test _checkpoint_path(run_dir, "best") == checkpoint_10
        @test _checkpoint_path(run_dir, "final") == checkpoint_final
        @test _checkpoint_path(run_dir, checkpoint_20) == checkpoint_20
    end

    mktempdir() do run_dir
        checkpoint_10 = joinpath(run_dir, "checkpoint_10.jls")
        checkpoint_20 = joinpath(run_dir, "checkpoint_20.jls")
        touch(checkpoint_10)
        touch(checkpoint_20)

        @test _checkpoint_path(run_dir, nothing) == checkpoint_20
        @test_throws ArgumentError _checkpoint_path(run_dir, "best")
    end
end

@testset "multi-run evaluation CLI and policy labels" begin
    options = _parse_cli([
        "runs/pr-drl",
        "--compare-run",
        "runs/a2c",
        "--compare-run",
        "runs/pr-drl-seed-2",
        "--output",
        "outputs/comparisons/test",
        "--episodes",
        "3",
    ])
    @test options.run_dir == "runs/pr-drl"
    @test options.comparison_runs == ["runs/a2c", "runs/pr-drl-seed-2"]
    @test options.output == "outputs/comparisons/test"
    @test options.episodes == 3
    @test_throws ArgumentError _parse_cli(["runs/pr-drl", "--compare-run"])
    generalization = _parse_cli(["runs/pr-drl", "--generalization-only"])
    @test generalization.generalization_only
    @test_throws ArgumentError _parse_cli([
        "runs/pr-drl",
        "--generalization-only",
        "--skip-generalization",
    ])
    @test_throws ArgumentError _parse_cli([
        "runs/pr-drl",
        "--generalization-only",
        "--compare-run",
        "runs/a2c",
    ])

    sources = [
        (algorithm=:pr_drl, run_id="pr-seed-1"),
        (algorithm=:a2c, run_id="a2c-seed-1"),
        (algorithm=:pr_drl, run_id="pr-seed-2"),
    ]
    specs = _comparison_policy_specs(sources)
    @test getfield.(specs, :key) == ["trained_pr_drl_1", "trained_a2c", "trained_pr_drl_2"]
    @test getfield.(specs, :label) == ["PR-DRL (pr-seed-1)", "A2C", "PR-DRL (pr-seed-2)"]
end

@testset "multi-policy flight comparison includes Odyssey" begin
    metrics = DataFrame([
        (
            policy=policy,
            maneuver_count=10,
            mission_duration_days=20.0,
            total_mission_delta_v_mps=25.0,
            thermal_violations=2,
            target_error_km=3.0,
        )
        for policy in ("trained_pr_drl", "trained_a2c", "aads_heuristic")
    ])
    specs = [
        (key="trained_pr_drl", label="PR-DRL"),
        (key="trained_a2c", label="A2C"),
        (key="aads_heuristic", label="AADS"),
    ]
    table = _flight_performance_table(metrics, specs)
    @test table.policy == ["PR-DRL", "A2C", "AADS", "Mars Odyssey"]
    @test table.target_reached_10km_count[1:3] == [1, 1, 1]
    mktempdir() do output_dir
        figure_path = joinpath(output_dir, "comparison.png")
        @test _paper_figure_11(metrics, figure_path, specs) == figure_path
        @test isfile(figure_path)
    end
end

@testset "generalization evaluation loss and Table VI record" begin
    ddqn_network = init_q_network(
        MersenneTwister(41);
        input_dim=3,
        hidden_dim=4,
        output_dim=2,
    )
    for field in (:W1, :b1, :W2, :b2, :W3, :b3)
        fill!(getfield(ddqn_network, field), 0f0)
    end
    transition = Transition(
        Float32[1, 2, 3],
        1,
        2f0,
        Float32[2, 3, 4],
        true,
        false,
        1,
    )
    ddqn_payload = Dict(
        :online => ddqn_network,
        :target => copy(ddqn_network),
        :config => DDQNConfig(obs_dim=3, action_dim=2),
    )
    @test _generalization_evaluation_loss(ddqn_payload, :pr_drl, [transition]) == 4.0
    @test _generalization_loss_definition(:pr_drl) == "ddqn_action_value_td_mse"

    critic = init_q_network(
        MersenneTwister(42);
        input_dim=3,
        hidden_dim=4,
        output_dim=1,
    )
    for field in (:W1, :b1, :W2, :b2, :W3, :b3)
        fill!(getfield(critic, field), 0f0)
    end
    a2c_payload = Dict(
        :critic => critic,
        :config => A2CConfig(obs_dim=3, action_dim=2),
    )
    @test _generalization_evaluation_loss(a2c_payload, :a2c, [transition]) == 4.0
    @test _generalization_loss_definition(:a2c) == "a2c_critic_value_td_mse"

    summary = EpisodeSummary(
        episode_reward=-3.0,
        success=true,
        thermal_violations=2,
        target_error_m=-2_000.0,
        mission_duration_days=5.0,
        total_delta_v_mps=4.0,
        pass_count=6,
        maneuver_count=3,
    )
    record = _generalization_record("nominal", (summaries=[summary],), 1.5, 1.0)
    @test record.generalization_gap == 0.5
    @test record.reached_goal_fraction == 1.0
    @test record.reached_goal_percent == 100.0
    @test record.mean_goal_distance_km == 2.0
end

@testset "generalization progress artifact" begin
    mktempdir() do output_dir
        progress_path = joinpath(output_dir, "progress.toml")
        started_at = time() - 2
        @test _write_generalization_progress(
            progress_path;
            status="running",
            current_case="nominal",
            case_index=2,
            case_count=7,
            case_episode=25,
            episodes_per_case=100,
            completed_episodes=125,
            total_episodes=700,
            started_at=started_at,
        ) == progress_path

        progress = TOML.parsefile(progress_path)
        @test progress["status"] == "running"
        @test progress["current_case"] == "nominal"
        @test progress["case_index"] == 2
        @test progress["case_episode"] == 25
        @test progress["completed_episodes"] == 125
        @test progress["percent_complete"] ≈ 100 * 125 / 700
        @test progress["elapsed_seconds"] >= 2
        @test progress["process_id"] == getpid()
        @test !isfile(progress_path * ".tmp")
    end
end
