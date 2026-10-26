using Test, SpaceAGORA, StarterName
using SpaceAGORA.Conformance: test_force_torque_model

@testset "StarterName" begin
    # Interface conformance of the study's own models with SpaceAGORA's contracts:
    # one call per model the study defines. A pass does not validate the science.
    test_force_torque_model(LowThrustPush(1.0e-3))

    @testset "short run of the case" begin
        solution = run_case(mission_time_s=10.0)
        @test last(solution.t) ≈ 10.0
        @test all(isfinite, solution.u[end])
    end

    # The study's own checks of its models and results go here.
end
