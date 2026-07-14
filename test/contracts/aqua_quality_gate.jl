# Migrated out of test/suites/05_thruster_control_and_quality_tests.jl's
# "Aqua Package Quality" testset (see test/TEST_RESTRUCTURE_PLAN.md, Phase 2)
# -- a package-quality policy check belongs with the other test/contracts/
# gates, not mixed into thruster/control unit tests. Self-contained like the
# other ci_*_gate.jl files: doesn't depend on test/helpers/bootstrap.jl.
using Test
using SpaceAGORA

const HAS_AQUA_QUALITY_GATE = let
    try
        @eval using Aqua
        true
    catch err
        @info "Skipping Aqua-backed test checks" exception=(err, catch_backtrace())
        false
    end
end

@testset "Aqua Package Quality" begin
    if HAS_AQUA_QUALITY_GATE
        Aqua.test_all(
            SpaceAGORA.SimulationModel;
            ambiguities=false,
            stale_deps=false,
            deps_compat=false,
            project_extras=false,
            piracies=false,
            persistent_tasks=false,
            undocumented_names=false
        )
    else
        @test_skip "Aqua is not available in this test environment"
    end
end
