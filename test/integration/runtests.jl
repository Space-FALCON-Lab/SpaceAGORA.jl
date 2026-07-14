include(joinpath(@__DIR__, "..", "helpers", "bootstrap.jl"))

# Each suite file is included inside its own nested @testset so a failure
# anywhere in one suite is recorded and aggregated rather than thrown
# immediately: Test.jl only raises once the outermost testset below finishes,
# so every suite always runs to completion and CI surfaces every failure from
# a single run instead of stopping at the first one.
@testset verbose=true "SpaceAGORA Suite Files" begin
    @testset "01 Contract and API Tests" begin
        include(joinpath(REPO_ROOT, "test", "suites", "01_contract_and_api_tests.jl"))
    end
    @testset "02 Callbacks Parallel and Smoke Tests" begin
        include(joinpath(REPO_ROOT, "test", "suites", "02_callbacks_parallel_and_smoke_tests.jl"))
    end
    @testset "03 Persistence Units and Rotational Tests" begin
        include(joinpath(REPO_ROOT, "test", "suites", "03_persistence_units_and_rotational_tests.jl"))
    end
    @testset "04 Solver Env and Regression Tests" begin
        include(joinpath(REPO_ROOT, "test", "suites", "04_solver_env_and_regression_tests.jl"))
    end
    @testset "05 Thruster Control and Quality Tests" begin
        include(joinpath(REPO_ROOT, "test", "suites", "05_thruster_control_and_quality_tests.jl"))
    end
    @testset "06 Monolith Split Runtime Tests" begin
        include(joinpath(REPO_ROOT, "test", "suites", "06_monolith_split_runtime_tests.jl"))
    end
    @testset "07 No GRAM Onboarding Tests" begin
        include(joinpath(REPO_ROOT, "test", "suites", "07_no_gram_onboarding_tests.jl"))
    end
    @testset "08 CLI and Assets Tests" begin
        include(joinpath(REPO_ROOT, "test", "suites", "08_cli_and_assets_tests.jl"))
    end
end
