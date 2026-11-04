include(joinpath(@__DIR__, "..", "Simulation", "Code", "VerifyOutput.jl"))
using .StudyVerification
using Test

const check_signed_totals = StudyVerification.check_signed_totals

@testset "Signed along-track verification tolerances" begin
    @test check_signed_totals(0.0, 0.0, 0.0, 0.0, 0.0) == 0.0
    @test check_signed_totals(2.0, -1.0, 1.0, 1.0, -0.5) == 0.0

    previous_negative = -1.107049466416928e-5
    negative = -1.1052884253105596e-5
    previous_positive = 2.8701347220491766
    positive = 2.8726363467434877
    residual = check_signed_totals(positive, negative, positive + negative,
        previous_positive, previous_negative)
    @test residual ≈ 1.7610411063683556e-8
    @test 1e-8 < residual < 1e-7

    tolerance = StudyVerification.SIGNED_T_MONOTONICITY_TOLERANCE_MPS
    @test tolerance == 1e-7
    @test check_signed_totals(0.0, 0.0, 0.0, tolerance, 0.0) == tolerance
    @test check_signed_totals(0.0, 0.0, 0.0, 0.0, -tolerance) == tolerance
    @test_throws ErrorException check_signed_totals(0.0, 0.0, 0.0,
        nextfloat(tolerance), 0.0)
    @test_throws ErrorException check_signed_totals(0.0, 0.0, 0.0,
        0.0, -nextfloat(tolerance))

    sum_tolerance = StudyVerification.SIGNED_T_SUM_TOLERANCE_MPS
    @test sum_tolerance == 1e-8
    @test check_signed_totals(0.0, 0.0, prevfloat(sum_tolerance), 0.0, 0.0) == 0.0
    @test_throws ErrorException check_signed_totals(0.0, 0.0, sum_tolerance, 0.0, 0.0)
    @test_throws ErrorException check_signed_totals(positive, negative,
        positive + negative + 5e-8, previous_positive, previous_negative)
    @test_throws ErrorException check_signed_totals(NaN, 0.0, 0.0, 0.0, 0.0)
end
