using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

module EngineSandbox
include(joinpath(Main.REPO_ROOT, "src", "simulation", "runtime_services.jl"))
include(joinpath(Main.REPO_ROOT, "src", "core", "simulation_model.jl"))
include(joinpath(Main.REPO_ROOT, "src", "simulation", "engine", "simulation_engine.jl"))
end

const SE = EngineSandbox.SimulationEngine

@testset "typed_config_equivalence" begin
    cfg = SE.SimulationEngineConfig(
        parallel=SE.ParallelConfig(profile="r5", effector_parallel_mode="on", rhs_batch_parallel_mode="off"),
        solver=SE.SolverConfig(mode="rodas5p", maxiters=12345, gravity_backbone_dt_s=4.0),
        runtime_policy=SE.RuntimePolicyConfig(
            warn_normalize=false,
            allow_typed_normalize=true,
            gram_per_sat_instances=true,
            srp_ephemeris_cache=false,
            nbody_ephemeris_cache=false,
            planet_frame_cache=true,
            spice_rhs_memo=false
        ),
        artifacts=SE.ArtifactConfig(save_bundle=false, warn_deprecated_config=false)
    )

    SE._with_engine_env_overrides(cfg, () -> begin
        @test SE._solver_policy_mode() == :rodas5p
        @test SE._solver_maxiters() == 12345
        @test SE._gravity_backbone_fixed_dt_s((; integration_tolerances=(; dt_max_orbit=9.0))) == 4.0
        @test SE._typed_save_bundle_enabled() == false
        @test SE._typed_normalize_warning_enabled() == false
        @test SE._typed_allow_transition_normalize() == true
        @test SE._gram_per_sat_instances_enabled() == true
    end)

    withenv(
        "SPACEAGORA_SOLVER_MODE" => "split_imex",
        "SPACEAGORA_SOLVER_MAXITERS" => "777",
        "SPACEAGORA_GRAVITY_BACKBONE_DT_S" => "6.5",
        "SPACEAGORA_SAVE_BUNDLE" => "0",
        "SPACEAGORA_WARN_NORMALIZE" => "0"
    ) do
        from_env = SE.simulation_engine_config_from_env(ENV)
        @test from_env.solver.mode == "split_imex"
        @test from_env.solver.maxiters == 777
        @test from_env.solver.gravity_backbone_dt_s == 6.5
        @test from_env.artifacts.save_bundle == false
        @test from_env.runtime_policy.warn_normalize == false

        SE._with_engine_env_overrides(from_env, () -> begin
            @test SE._solver_policy_mode() == :split_imex
            @test SE._solver_maxiters() == 777
            @test SE._typed_save_bundle_enabled() == false
            @test SE._typed_normalize_warning_enabled() == false
        end)
    end
end

println("typed_config_equivalence_gate_ok")
