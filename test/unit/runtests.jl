include(joinpath(@__DIR__, "..", "helpers", "bootstrap.jl"))

@testset "Unit Test Placeholder" begin
    @test true
end

include("rpo_port_tests.jl")
include(joinpath(@__DIR__, "robotics", "runtests.jl"))

# Migrated out of the legacy test/suites/NN_*.jl files one domain at a time;
# see test/TEST_RESTRUCTURE_PLAN.md (Phase 2) for order/status.
@testset verbose=true "Unit Tests (migrated from legacy suites)" begin
    @testset "Core" begin
        include(joinpath(@__DIR__, "core", "api_convenience_constructor_tests.jl"))
    end
    @testset "Simulation Engine" begin
        include(joinpath(@__DIR__, "simulation_engine", "solver_env_helpers_tests.jl"))
        include(joinpath(@__DIR__, "simulation_engine", "run_metadata_and_rhs_tests.jl"))
        include(joinpath(@__DIR__, "simulation_engine", "callback_tests.jl"))
        include(joinpath(@__DIR__, "simulation_engine", "state_isolation_and_calibration_tests.jl"))
    end
    @testset "Dynamics" begin
        include(joinpath(@__DIR__, "dynamics", "rigid_body_and_rhs_tests.jl"))
        include(joinpath(@__DIR__, "dynamics", "orbital_energy_and_drift_tests.jl"))
        include(joinpath(@__DIR__, "dynamics", "rotational_tests.jl"))
        include(joinpath(@__DIR__, "dynamics", "aerodynamic_helper_tests.jl"))
    end
    @testset "Vehicle" begin
        include(joinpath(@__DIR__, "vehicle", "thruster_tests.jl"))
    end
    @testset "IO" begin
        include(joinpath(@__DIR__, "io", "persistence_and_checkpoint_tests.jl"))
    end
    @testset "Environment" begin
        include(joinpath(@__DIR__, "environment", "planet_constructor_tests.jl"))
        include(joinpath(@__DIR__, "environment", "no_gram_onboarding_tests.jl"))
    end
    @testset "Simulation" begin
        include(joinpath(@__DIR__, "simulation", "campaign_and_gram_lock_tests.jl"))
    end
    @testset "Analysis" begin
        include(joinpath(@__DIR__, "analysis", "telemetry_and_policy_tests.jl"))
    end
    @testset "CLI" begin
        include(joinpath(@__DIR__, "cli", "cli_and_asset_tests.jl"))
    end
    @testset "Parallel" begin
        include(joinpath(@__DIR__, "parallel", "policy_tests.jl"))
        include(joinpath(@__DIR__, "parallel", "adaptive_routing_tests.jl"))
    end
    @testset "GNC" begin
        include(joinpath(@__DIR__, "gnc", "control_effector_tests.jl"))
        include(joinpath(REPO_ROOT, "test", "gnc", "aerobraking", "e_edg_strategy_parity_tests.jl"))
        include(joinpath(REPO_ROOT, "test", "gnc", "aerobraking", "t_edg_strategy_parity_tests.jl"))
        include(joinpath(REPO_ROOT, "test", "mission", "aerobraking_policy_selector_stub_tests.jl"))
    end
    @testset "Mission" begin
        include(joinpath(@__DIR__, "mission", "maneuver_and_campaign_tests.jl"))
    end
    @testset "Coverage Probes (pending Phase 4 audit)" begin
        include(joinpath(REPO_ROOT, "test", "coverage_parallel_telemetry_probes.jl"))
        include(joinpath(REPO_ROOT, "test", "coverage_runtime_boundary_probes.jl"))
        include(joinpath(REPO_ROOT, "test", "coverage_targeted_90_probes.jl"))
    end
end

# Architecture/export/include-order contracts also live at
# test/contracts/architecture_and_export_contracts.jl (run there via
# test/contracts/runtests.jl); included here too so the default
# `test/runtests.jl` doesn't lose this coverage now that suite 01 is gutted.
@testset "Architecture Contracts" begin
    include(joinpath(REPO_ROOT, "test", "contracts", "architecture_and_export_contracts.jl"))
end
