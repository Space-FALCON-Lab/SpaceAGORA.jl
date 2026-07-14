include(joinpath(@__DIR__, "..", "helpers", "bootstrap.jl"))

# See test/README.md for how to add a new scenario here, and
# test/TEST_RESTRUCTURE_PLAN.md for the harness's status/rationale.
@testset verbose=true "Golden Propagation Regression Tests" begin
    @testset "AGORA Earth Regression" begin
        include(joinpath(@__DIR__, "agora_earth", "config.jl"))
        run_and_compare_golden("agora_earth", build_agora_earth_config; spacecraft_count=1)
    end

    @testset "Constellation Regression" begin
        include(joinpath(@__DIR__, "constellation", "config.jl"))
        run_and_compare_golden("constellation", build_constellation_config; spacecraft_count=3)
    end

    @testset "GRAM Mars Regression" begin
        if !HAS_GRAMSUITE
            @test_skip "GRAMSuite is not available in this checkout"
        elseif !golden_scenario_enabled("gram_mars")
            @test_skip "gram_mars is a nightly-tier scenario; set SPACEAGORA_GOLDEN_TIER=all to include it"
        else
            include(joinpath(@__DIR__, "gram_mars", "config.jl"))
            run_and_compare_golden("gram_mars", build_gram_mars_config; spacecraft_count=1)
        end
    end
end
